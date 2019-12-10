/*
imgalt - Image Alignment Tool
Author: GreatAttractor

version 0.5
2014/05/22

This code can be freely distributed and used for any purpose.

File description:
    Utility functions implementation.

*/
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdexcept>
#include <fstream>
#include <cctype>
#include <algorithm>
#include <boost/format.hpp>
#include "comp.h"
using namespace boost;


int GetBytesPerPixel(PixelFormat_t pixFmt)
{
    switch (pixFmt)
    {
    case PIX_PAL8:
    case PIX_MONO8:
        return 1;
    case PIX_MONO16: return 2;
    case PIX_RGB24: return 3;
    case PIX_RGB48: return 3*2;
    case PIX_MONO32F: return sizeof(float);
    default: return -1;
    }
}


/// Returns the smallest power of 2 which is > n
unsigned GetClosestGPowerOf2(unsigned n)
{
    int msb = 0;
    while (n != 0)
    {
        n >>= 1;
        msb++;
    }

    return ((unsigned)1 << msb);
}

/// Returns the greatest power of 2 which is <= n
unsigned GetClosestLEPowerOf2(unsigned n)
{
    int msb = 0;
    while (n != 0)
    {
        n >>= 1;
        msb++;
    }

    return ((unsigned)1 << (msb-1));
}




/// Returns 0 for x=0, 1 for x=1
inline float CosineWindow(float x)
{
    return 0.5f*(1.0f - cosf(3.1415926535f*x));
}

/// Returns 0 for x=0, 1 for x=1
inline float BlackmanWindow(float x)
{
    const float A0 = 7938.0f/18608,
                A1 = 9240.0f/18608,
                A2 = 1430.0f/18608;

    return A0 - A1*cosf(3.1415926535f*x) + A2*cosf(2*3.1415926535f*x);
}

#define SQR(x) ((x)*(x))

/// Calculates window function and writes its values into 'buf' (a square array)
void CalcWindowFunction(
        int wndSize, ///< Window size
        float buf[] ///< Destination buffer, wndSize*wndSize elements
)
{
    // The window function is rotationally symmetrical, so calculate it only in a quarter of 'buf'
    #pragma omp parallel for
    for (int y = 0; y < wndSize/2; y++)
        for (int x = 0; x < wndSize/2; x++)
        {
            float value = 0.0f;
            float dist = sqrtf(SQR(wndSize/2 - x) + SQR(wndSize/2 - y));
            if (dist < wndSize/2)
                value = BlackmanWindow(1.0f-dist/(wndSize/2));

            // upper left
            buf[x + y*wndSize] = value;
            // upper right
            buf[wndSize-1-x + y*wndSize] = value;
            // lower right
            buf[wndSize-1-x + (wndSize-1-y)*wndSize] = value;
            // lower left
            buf[x + (wndSize-1-y)*wndSize] = value;
        }
}

/// Multiplies input image by window function; input and output images may be the same
void ApplyWindowFunction(float img[], ///< Input image, size*size elements; may be the same pointer as 'dest'
        float windowFunc[], ///< Window function values, size*size elements
        float dest[], ///< Output buffer, size*size elements; may be the same pointer as 'img'
        int size ///< Number of rows and columns in each array
)
{
    #pragma omp parallel for
    for (int y = 0; y < size; y++)
        for (int x = 0; x < size; x++)
            dest[x + y*size] = img[x + y*size] * windowFunc[x + y*size];
}

/// Resizes and translates image (or its fragment) by cropping and/or padding (with zeros) to the destination size and offset (there is no scaling)
void ResizeAndTranslate(
        void *input,          ///< Input buffer
        PixelFormat_t pixFmt, ///< Format of data in 'input' and 'output'
        int srcWidth,         ///< Width of input image
        int srcHeight,        ///< Height of input image
        int srcXmin,          ///< X min of input data in input image
        int srcYmin,          ///< Y min of input data in input image
        int srcXmax,          ///< X max of input data in input image
        int srcYmax,          ///< Y max of input data in input image
        void *output,         ///< Output buffer
        int destWidth,        ///< Width of output image
        int destHeight,       ///< Height of output image
        int xOfs,             ///< X offset of input data in output buffer
        int yOfs              ///< Y offset of input data in output buffer
)
{
    int bytesPP = GetBytesPerPixel(pixFmt);
    memset(output, 0, destWidth*destHeight*bytesPP); // Works also if 'dest' points to an array of floats; 32 zero bits represent a floating-point 0.0f

    // start and end (inclusive) coordinates to fill in the output buffer
    unsigned Xstart = (xOfs < 0) ? 0 : xOfs;
    unsigned Ystart = (yOfs < 0) ? 0 : yOfs;

    unsigned Xend = std::max(0, std::min((int)xOfs + srcXmax, destWidth-1));
    unsigned Yend = std::max(0, std::min((int)yOfs + srcYmax, destHeight-1));

    for (int y = Ystart; y <= Yend; y++)
    {
        memcpy((uint8_t *)output + (Xstart + y*destWidth)*bytesPP,
               (uint8_t *)input + (Xstart-xOfs+srcXmin + (y-yOfs + srcYmin)*srcWidth)*bytesPP,
               (Xend - Xstart + 1)*bytesPP);
    }
}

/// Blurs 'src' image using a 5x5 Gaussian kernel and writes the result to 'dest'; 2 pixel-wide borders are not blurred
void BlurImage(
    uint8_t *src, ///< Source image (8-bit luminance)
    int width,    ///< Image width
    int height,   ///< Image height
    uint8_t *dest ///< Destination image (same dimensions as 'src')
)
{
    const int KERNEL_SIZE = 5;
    const int kernel[KERNEL_SIZE][KERNEL_SIZE] =
    {
        { 2, 4,  5,  4,  2 },
        { 4, 9,  12, 9,  4 },
        { 5, 12, 15, 12, 5 },
        { 4, 9,  12, 9,  4 },
        { 2, 4,  5,  4,  2 }
    };

    for (int y = KERNEL_SIZE/2; y <= height - KERNEL_SIZE/2 - 1; y++)
        for (int x = KERNEL_SIZE/2; x <= width - KERNEL_SIZE/2 - 1; x++)
        {
            int sum = 0;
            for (int yofs = -KERNEL_SIZE/2; yofs <= KERNEL_SIZE/2; yofs++)
                for (int xofs = -KERNEL_SIZE/2; xofs <= KERNEL_SIZE/2; xofs++)
                    sum += (int)src[(x + xofs) + (y + yofs) * width] * kernel[xofs + KERNEL_SIZE/2][yofs + KERNEL_SIZE/2];

            dest[x + y*width] = sum/159;
        }
}

/// Calculates squared gradient lengths in the source image; 3 pixel-wide borders are skipped
void CalcGradients(
    uint8_t *src,  ///< Source image
    int width,     ///< Image width
    int height,    ///< Image height
    uint32_t *dest ///< Destination buffer (same dimensions as 'src')
)
{
    const int KERNEL_SIZE = 3;
    const int kernelX[KERNEL_SIZE][KERNEL_SIZE] =
    {
        { -1, 0, 1 },
        { -2, 0, 2 },
        { -1, 0, 1 }
    };

    const int kernelY[KERNEL_SIZE][KERNEL_SIZE] =
    {
        { 1,   2,  1 },
        { 0,   0,  0 },
        { -1, -2, -1 }
    };

    /// By skipping 3-pixel borders we skip the pixels that are also skipped by BlurImage()
    for (int y = 3; y <= height - 4; y++)
        for (int x = 3; x <= width - 4; x++)
        {
            int gradX = 0, gradY = 0;
            for (int yofs = -KERNEL_SIZE/2; yofs <= KERNEL_SIZE/2; yofs++)
                for (int xofs = -KERNEL_SIZE/2; xofs <= KERNEL_SIZE/2; xofs++)
                {
                    uint8_t srcVal = src[(x + xofs) + (y + yofs) * width];
                    gradX += srcVal * kernelX[xofs + KERNEL_SIZE/2][yofs + KERNEL_SIZE/2];
                    gradY += srcVal * kernelY[xofs + KERNEL_SIZE/2][yofs + KERNEL_SIZE/2];
                }

            dest[x + y*width] = gradX*gradX + gradY*gradY; // this is at most 2*(4*255)^2, so fits easily in uint32_t
        }
}

namespace BMP
{

#pragma pack(push)

#pragma pack(1)
typedef struct
{
    uint16_t bfType;
    uint32_t bfSize;
    uint16_t bfReserved1;
    uint16_t bfReserved2;
    uint32_t bfOffBits;
} BITMAPFILEHEADER_t; ///< BMP file header

#pragma pack(1)
typedef struct
{
   uint32_t biSize;
   int32_t biWidth;
   int32_t biHeight;
   uint16_t biPlanes;
   uint16_t biBitCount;
   uint32_t biCompression;
   uint32_t biSizeImage;
   int32_t biXPelsPerMeter;
   int32_t biYPelsPerMeter;
   uint32_t biClrUsed;
   uint32_t biClrImportant;
} BITMAPINFOHEADER_t; ///< BMP info header

#pragma pack(pop)

const uint32_t BMP_NO_COMPRESSION = 0;
const int BMP_PALETTE_SIZE = 256*4;

// returns the least multiple of 4 which is >= x
#define UP4MULT(x) (((x)+3)/4*4)

/// Reads a BMP image and returns pointer to the newly allocated buffer with pixel contents or 0 on error
void *ReadBmp(const char *fileName,
        int &imgWidth,  ///< Receives image width
        int &imgHeight, ///< Receives image height
        PixelFormat_t &pixFmt, ///< Receives the pixel format
        uint8_t palette[]    ///< If not NULL and reading from an 8-bit file, receives the palette (1024 bytes)
)
{
    std::ifstream file(fileName, std::ios_base::in | std::ios_base::binary);
    if (file.fail())
        return 0;

    BITMAPFILEHEADER_t bmpFileHdr;
    BITMAPINFOHEADER_t bmpInfoHdr;

    file.read((char *)&bmpFileHdr, sizeof(bmpFileHdr));
    file.read((char *)&bmpInfoHdr, sizeof(bmpInfoHdr));
    if (file.eof())
        return 0;

    imgWidth = bmpInfoHdr.biWidth;
    imgHeight = bmpInfoHdr.biHeight;
    if (imgWidth == 0 || imgHeight == 0 ||
        bmpFileHdr.bfType != 'B'+((int)'M'<<8) ||
        bmpInfoHdr.biPlanes != 1 ||
        bmpInfoHdr.biBitCount != 8 && bmpInfoHdr.biBitCount != 24 ||
        bmpInfoHdr.biCompression != BMP_NO_COMPRESSION)
    {
        return 0;
    }

    if (bmpInfoHdr.biBitCount == 8)
        pixFmt = PIX_PAL8;
    else if (bmpInfoHdr.biBitCount == 24)
        pixFmt = PIX_RGB24;

    int bytesPP = GetBytesPerPixel(pixFmt);

    void *pixels = malloc(imgWidth * imgHeight * bytesPP);

    if (bmpInfoHdr.biBitCount == 8)
    {
        unsigned bmpStride = UP4MULT(imgWidth); // line length in bytes in the BMP file's pixel data
        unsigned skip = bmpStride - imgWidth; // number of padding bytes at the end of a line

        int actualPalSize = bmpInfoHdr.biClrUsed == 0 ? BMP_PALETTE_SIZE : bmpInfoHdr.biClrUsed*4;

        // seek to the beginning of palette
        file.seekg(sizeof(bmpFileHdr) + bmpInfoHdr.biSize, std::ios_base::beg);

        if (palette != 0)
            file.read((char *)palette, actualPalSize);

        // Seek to the beginning of pixel values
        file.seekg(bmpFileHdr.bfOffBits, std::ios_base::beg);

        for (int y = imgHeight - 1; y >= 0; y--) // lines in BMP are stored bottom to top
        {
            file.read((char *)((uint8_t *)pixels + y*imgWidth), imgWidth);
            if (skip > 0)
                file.seekg(skip, std::ios_base::cur);
        }
    }
    else if (bmpInfoHdr.biBitCount == 24)
    {
        unsigned bmpStride = UP4MULT(imgWidth*3); // line length in bytes in the BMP file's pixel data
        unsigned skip = bmpStride - imgWidth*3; // number of padding bytes at the end of a row

        // Seek to the beginning of pixel values
        file.seekg(bmpFileHdr.bfOffBits, std::ios_base::beg);

        // read the lines directly into the buffer
        for (int y = imgHeight - 1; y >= 0; y--) // lines in BMP are stored bottom to top
        {
            file.read((char *)((uint8_t *)pixels + y*imgWidth*3), imgWidth*3);
            if (skip > 0)
                file.seekg(skip, std::ios_base::cur);
        }
    }

    return pixels;
}

/// Saves image in BMP format; returns 'false' on error
bool SaveBmp(const char *fileName, ///< Output file name
        int imgWidth,              ///< Image width
        int imgHeight,             ///< Image height
        PixelFormat_t pixFmt,      ///< Pixel format; has to be PIX_PAL8 or PIX_RGB24
        void *pixels,              ///< Pixel contents in 'pixFmt' format
        uint8_t palette[]          ///< Points to the palette (1024 bytes) to be saved if pixFmt is PIX_PAL8
)
{
    BITMAPFILEHEADER_t bmfh;
    BITMAPINFOHEADER_t bmih;
    int i;

    int bytesPP = GetBytesPerPixel(pixFmt);
    unsigned bmpLineWidth = UP4MULT(imgWidth * bytesPP);

    bmfh.bfType = 'B'+((int)'M'<<8);
    bmfh.bfSize = sizeof(bmfh) + sizeof(bmih) + imgHeight*bmpLineWidth;
    if (pixFmt == PIX_PAL8)
        bmfh.bfSize += BMP_PALETTE_SIZE;
    bmfh.bfReserved1 = 0;
    bmfh.bfReserved2 = 0;
    bmfh.bfOffBits = sizeof(bmih) + sizeof(bmfh);
    if (pixFmt == PIX_PAL8)
        bmfh.bfOffBits += BMP_PALETTE_SIZE;

    bmih.biSize = sizeof(bmih);
    bmih.biWidth = imgWidth;
    bmih.biHeight = imgHeight;
    bmih.biPlanes = 1;
    bmih.biBitCount = bytesPP * 8;
    bmih.biCompression = BMP_NO_COMPRESSION;
    bmih.biSizeImage = 0;
    bmih.biXPelsPerMeter = 1000;
    bmih.biYPelsPerMeter = 1000;
    bmih.biClrUsed = 0;
    bmih.biClrImportant = 0;

    std::ofstream file(fileName, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    if (file.fail())
        return false;

    file.write((const char *)&bmfh, sizeof(bmfh));
    file.write((const char *)&bmih, sizeof(bmih));
    if (pixFmt == PIX_PAL8)
        file.write((const char *)palette, BMP_PALETTE_SIZE);

    int skip = bmpLineWidth - imgWidth*bytesPP;

    for (i = imgHeight - 1; i >= 0; i--) // lines in BMP are stored bottom to top
    {
        file.write((const char*)pixels + i * imgWidth * bytesPP, imgWidth*bytesPP);
        if (skip > 0)
            file.write((const char *)pixels, skip); //this is just padding, so write anything
    }

    file.close();

    return true;
}

bool GetBmpDimensions(const char *fileName, unsigned &imgWidth, unsigned &imgHeight)
{
    std::ifstream file(fileName, std::ios_base::in | std::ios_base::binary);
    if (file.fail())
        return false;

    BITMAPFILEHEADER_t bmpFileHdr;
    BITMAPINFOHEADER_t bmpInfoHdr;

    file.read((char *)&bmpFileHdr, sizeof(bmpFileHdr));
    file.read((char *)&bmpInfoHdr, sizeof(bmpInfoHdr));
    if (file.eof())
        return false;

    imgWidth = bmpInfoHdr.biWidth;
    imgHeight = bmpInfoHdr.biHeight;

    return true;
}

} // namespace BMP

namespace TIFF
{

#pragma pack(push, 1)
typedef struct
{
    uint16_t tag;
    uint16_t type;
    uint32_t count;
    uint32_t value;
} TiffField_t;

typedef struct
{
    uint16_t id;
    uint16_t version;
    uint32_t dirOffset; // = offset of 'numDirEntries'
} TiffHeader_t;

#pragma pack(pop)

typedef enum { ttByte = 1, ttAscii, ttWord, ttDWord, ttRational } TagType_t;

const int TIFF_VERSION = 42;
const int TAG_IMAGE_WIDTH =                0x100;
const int TAG_IMAGE_HEIGHT =               0x101;
const int TAG_BITS_PER_SAMPLE =            0x102;
const int TAG_COMPRESSION =                0x103;
const int TAG_PHOTOMETRIC_INTERPRETATION = 0x106;
const int TAG_STRIP_OFFSETS =              0x111;
const int TAG_SAMPLES_PER_PIXEL =          0x115;
const int TAG_ROWS_PER_STRIP =             0x116;
const int TAG_STRIP_BYTE_COUNTS =          0x117;
const int TAG_PLANAR_CONFIGURATION =       0x11C;

const uint16_t NO_COMPRESSION = 1;
const uint16_t PLANAR_CONFIGURATION_CHUNKY = 1;
const uint16_t INTEL_BYTE_ORDER = ((uint16_t)'I' << 8) + 'I'; // little-endian
const uint16_t MOTOROLA_BYTE_ORDER = ((uint16_t)'M' << 8) + 'M'; // big-endian
const int PHMET_WHITE_IS_ZERO = 0;
const int PHMET_BLACK_IS_ZERO = 1;
const int PHMET_RGB = 2;

inline unsigned GetFieldTypeLength(TagType_t ttt)
{
    switch (ttt)
    {
    case ttByte: return 1; break;
    case ttAscii: return 1; break;
    case ttWord: return 2; break;
    case ttDWord: return 4; break;
    case ttRational: return 8; break;
    }
}

/// Conditionally swaps a 32-bit value
uint32_t inline SWAP32cnd(uint32_t x, bool swap)
{
    if (swap) return (x << 24) | ((x & 0x00FF0000) >> 8) | ((x & 0x0000FF00) << 8) | (x >> 24);
    else return x;
}

/// Conditionally swaps two lower bytes of a 32-bit value
uint32_t inline SWAP16in32cnd(uint32_t x, bool swap)
{
    if (swap) return ((x & 0xFF) << 8) | (x >> 8);
    else return x;
}

uint16_t inline SWAP16cnd(uint16_t x, bool swap)
{
    if (swap) return (x << 8) | (x >> 8);
    else return x;
}

/// Changes endianess of 16-bit words in the specified buffer
void SwapBufferWords(uint16_t *buf, int numWords)
{
    for (unsigned i = 0; i < numWords; i++)
        buf[i] = (buf[i] << 8) | (buf[i] >> 8);
}

/// Reverses values of an 8-bit grayscale buffer
void NegateGrayscale8(uint8_t *buf, int length)
{
    for (int i = 0; i < length; i++)
        buf[i] = 0xFF - buf[i];
}

/// Reverses values of a 16-bit grayscale buffer
void NegateGrayscale16(uint16_t *buf, int length)
{
    for (int i = 0; i < length; i++)
        buf[i] = 0xFFFF - buf[i];
}

bool SaveTiff(const char *fileName, ///< Output file name
              void *pixels,   ///< Pointer to the buffer with pixel data (left to right, top to bottom, no padding)
              PixelFormat_t pixFmt, ///< Pixel format of 'pixels'
              int imgWidth,
              int imgHeight
)
{
    if (pixFmt != PIX_MONO8 && pixFmt != PIX_MONO16 && pixFmt != PIX_RGB24 && pixFmt != PIX_RGB48)
        throw std::runtime_error("SaveTiff(): only grayscale and RGB, 8- and 16-bit formats are supported.");

    std::ofstream file(fileName, std::ios_base::trunc | std::ios_base::binary);

    if (file.fail())
        return false;

    TiffHeader_t tiffHeader;
    tiffHeader.id = INTEL_BYTE_ORDER;
    tiffHeader.version = TIFF_VERSION;
    tiffHeader.dirOffset = sizeof(tiffHeader);
    file.write((const char *)&tiffHeader, sizeof(tiffHeader));

    uint16_t numDirEntries = 10;
    file.write((const char *)&numDirEntries, sizeof(numDirEntries));

    uint32_t nextDirOffset = 0;

    TiffField_t field;

    field.tag = TAG_IMAGE_WIDTH;
    field.type = ttWord;
    field.count = 1;
    field.value = imgWidth;
    file.write((const char *)&field, sizeof(field));

    field.tag = TAG_IMAGE_HEIGHT;
    field.type = ttWord;
    field.count = 1;
    field.value = imgHeight;
    file.write((const char *)&field, sizeof(field));

    field.tag = TAG_BITS_PER_SAMPLE;
    field.type = ttWord;
    field.count = 1;
    switch (pixFmt)
    {
    case PIX_MONO8:
    case PIX_RGB24:
        field.value = 8; break;
    case PIX_MONO16:
    case PIX_RGB48:
        field.value = 16; break;
    }
    file.write((const char *)&field, sizeof(field));

    field.tag = TAG_COMPRESSION;
    field.type = ttWord;
    field.count = 1;
    field.value = NO_COMPRESSION;
    file.write((const char *)&field, sizeof(field));

    field.tag = TAG_PHOTOMETRIC_INTERPRETATION;
    field.type = ttWord;
    field.count = 1;
    switch (pixFmt)
    {
    case PIX_MONO8:
    case PIX_MONO16:
        field.value = PHMET_BLACK_IS_ZERO; break;
    case PIX_RGB24:
    case PIX_RGB48:
        field.value = PHMET_RGB; break;
    }
    file.write((const char *)&field, sizeof(field));

    field.tag = TAG_STRIP_OFFSETS;
    field.type = ttDWord;
    field.count = 1;
    // we write the header, num. of directory entries, 10 fields and a next directory offset (==0); pixel data starts next
    field.value = sizeof(tiffHeader) + sizeof(numDirEntries) + 10*sizeof(field) + sizeof(nextDirOffset);
    file.write((const char *)&field, sizeof(field));

    field.tag = TAG_SAMPLES_PER_PIXEL;
    field.type = ttWord;
    field.count = 1;
    switch (pixFmt)
    {
    case PIX_MONO8:
    case PIX_MONO16:
        field.value = 1; break;
    case PIX_RGB24:
    case PIX_RGB48:
        field.value = 3; break;
    }
    file.write((const char *)&field, sizeof(field));

    field.tag = TAG_ROWS_PER_STRIP;
    field.type = ttWord;
    field.count = 1;
    field.value = imgHeight; // there is only one strip for the whole image
    file.write((const char *)&field, sizeof(field));

    field.tag = TAG_STRIP_BYTE_COUNTS;
    field.type = ttDWord;
    field.count = 1;
    field.value = imgWidth * imgHeight * GetBytesPerPixel(pixFmt); // there is only one strip for the whole image
    file.write((const char *)&field, sizeof(field));

    field.tag = TAG_PLANAR_CONFIGURATION;
    field.type = ttWord;
    field.count = 1;
    field.value = PLANAR_CONFIGURATION_CHUNKY; // there is only one strip for the whole image
    file.write((const char *)&field, sizeof(field));

    // write the next directory offset (0 = no other directories)
    file.write((const char *)&nextDirOffset, sizeof(nextDirOffset));

    file.write((const char *)pixels, imgWidth * imgHeight * GetBytesPerPixel(pixFmt));

    file.close();
    return true;
}

/// Returns newly allocated buffer with contents of the specified TIFF file (returns 0 on error)
void *ReadTiff(const char *fileName, ///< Input file name
              PixelFormat_t &pixFmt, ///< Receives the pixel format
              int &imgWidth, ///< Receives image width
              int &imgHeight, ///< Receives image height
              std::string *errorMsg ///< If not null, receives error message (if any)
)
{
    std::ifstream file(fileName, std::ios_base::binary);

    if (file.fail())
        return 0;

    TiffHeader_t tiffHeader;
    file.read((char *)&tiffHeader, sizeof(tiffHeader));
    if (file.gcount() != sizeof(tiffHeader))
    {
        if (errorMsg) *errorMsg = "File header is incomplete.";
        return 0;
    }

    bool isBE = tiffHeader.id == MOTOROLA_BYTE_ORDER; // true if the file has big endian data

    if (SWAP16cnd(tiffHeader.version, isBE) != TIFF_VERSION)
    {
        if (errorMsg) *errorMsg = "Unknown TIFF version.";
        return 0;
    }

    // Seek to the first TIFF directory
    file.seekg(SWAP32cnd(tiffHeader.dirOffset, isBE), std::ios_base::beg);

    uint16_t numDirEntries;
    file.read((char *)&numDirEntries, sizeof(numDirEntries));
    numDirEntries = SWAP16cnd(numDirEntries, isBE);
    if (file.gcount() != sizeof(numDirEntries))
    {
        if (errorMsg) *errorMsg = "The number of TIFF directory entries tag is incomplete.";
        return 0;
    }

    unsigned numStrips = 0;
    unsigned bitsPerSample = 0;
    unsigned *stripOffsets = 0;
    unsigned *stripByteCounts = 0;
    unsigned rowsPerStrip = 0;
    int photometricInterpretation = -1;
    int samplesPerPixel = 0;

    std::fstream::pos_type nextFieldPos = file.tellg();
    for (unsigned i = 0; i < numDirEntries; i++)
    {
        TiffField_t tiffField;

        file.seekg(nextFieldPos, std::ios_base::beg);
        file.read((char *)&tiffField, sizeof(tiffField));
        if (file.gcount() != sizeof(tiffField))
        {
            if (errorMsg) *errorMsg = "TIFF field is incomplete.";
            return 0;
        }
        nextFieldPos = file.tellg();

        tiffField.tag = SWAP16cnd(tiffField.tag, isBE);
        tiffField.type = SWAP16cnd(tiffField.type, isBE);
        tiffField.count = SWAP32cnd(tiffField.count, isBE);
        if (tiffField.count > 1 || tiffField.type == ttDWord)
            tiffField.value = SWAP32cnd(tiffField.value, isBE);
        else if (tiffField.count == 1 && tiffField.type == ttWord)
            tiffField.value = SWAP16in32cnd(tiffField.value, isBE);

        switch (tiffField.tag)
        {
        case TAG_IMAGE_WIDTH: imgWidth = tiffField.value; break;

        case TAG_IMAGE_HEIGHT: imgHeight = tiffField.value; break;

        case TAG_BITS_PER_SAMPLE:
            if (tiffField.count == 1)
                bitsPerSample = tiffField.value;
            else
            {
                // Some files may have as many "bits per sample" values specified
                // as there are channels. Make sure they are all the same.

                file.seekg(tiffField.value, std::ios_base::beg);

                uint16_t *fieldBuf = new uint16_t[tiffField.count];
                file.read((char *)fieldBuf, tiffField.count * sizeof(uint16_t));

                bool allEqual = true;
                uint16_t first = fieldBuf[0];
                for (unsigned j = 1; j < tiffField.count; j++)
                    if (fieldBuf[j] != first)
                    {
                        allEqual = false;
                        break;
                    }

                 if (!allEqual)
                 {
                    if (errorMsg) *errorMsg = "Files with differing bit depts per channel are not supported.";
                    return 0;
                 }

                 bitsPerSample = SWAP16cnd(first, isBE);
            }

            if (bitsPerSample != 8 && bitsPerSample != 16)
            {
                if (errorMsg) *errorMsg = "Only 8 and 16 bits per channel files are supported.";
                return 0;
            }
            break;

        case TAG_COMPRESSION:
            if (tiffField.value != NO_COMPRESSION)
            {
                if (errorMsg) *errorMsg = "Compression is not supported.";
                return 0;
            }
            break;

        case TAG_PHOTOMETRIC_INTERPRETATION: photometricInterpretation = tiffField.value; break;

        case TAG_STRIP_OFFSETS:
            numStrips = tiffField.count;
            stripOffsets = new unsigned[numStrips];
            if (numStrips == 1)
                stripOffsets[0] = tiffField.value;
            else
            {
                file.seekg(tiffField.value, std::ios_base::beg);
                for (unsigned i = 0; i < numStrips; i++)
                {
                    file.read((char *)&stripOffsets[i], sizeof(stripOffsets[i]));
                    stripOffsets[i] = SWAP32cnd(stripOffsets[i], isBE);
                }
            }
            break;

        case TAG_SAMPLES_PER_PIXEL: samplesPerPixel = tiffField.value; break;

        case TAG_ROWS_PER_STRIP: rowsPerStrip = tiffField.value; break;

        case TAG_STRIP_BYTE_COUNTS:
            stripByteCounts = new unsigned[tiffField.count];
            if (tiffField.count == 1)
                stripByteCounts[0] = tiffField.value;
            else
            {
                file.seekg(tiffField.value, std::ios_base::beg);
                for (unsigned i = 0; i < tiffField.count; i++)
                {
                    file.read((char *)&stripByteCounts[i], sizeof(stripByteCounts[i]));
                    stripByteCounts[i] = SWAP32cnd(stripByteCounts[i], isBE);
                }
            }
            break;

        case TAG_PLANAR_CONFIGURATION:
            if (tiffField.value != PLANAR_CONFIGURATION_CHUNKY)
            {
                if (errorMsg) *errorMsg = "Files with planar configuration other than packed (chunky) are not supported.";
                return 0;
            }
            break;
        }
    }

    if (rowsPerStrip == 0 && numStrips == 1)
        // If there is only 1 strip, it contains all the rows
        rowsPerStrip = imgHeight;

    // Validate the values

    if (samplesPerPixel == 1 && photometricInterpretation != PHMET_BLACK_IS_ZERO && photometricInterpretation != PHMET_WHITE_IS_ZERO ||
        samplesPerPixel == 3 && photometricInterpretation != PHMET_RGB ||
        samplesPerPixel != 1 && samplesPerPixel != 3)
    {
        if (errorMsg) *errorMsg = "Only RGB and grayscale images are supported.";
        return 0;
    }

    if (samplesPerPixel == 1)
    {
        if (bitsPerSample == 8)
            pixFmt = PIX_MONO8;
        else if (bitsPerSample == 16)
            pixFmt = PIX_MONO16;
    }
    else if (samplesPerPixel == 3)
    {
        if (bitsPerSample == 8)
            pixFmt = PIX_RGB24;
        else if (bitsPerSample == 16)
            pixFmt = PIX_RGB48;
    }

    // Buffer with all image pixel values, left to right, top to bottom, without any padding
    void *pixels = malloc(imgWidth * imgHeight * GetBytesPerPixel(pixFmt));

    int bufOfs = 0;
    for (unsigned i = 0; i < numStrips; i++)
    {
        file.seekg(stripOffsets[i], std::ios_base::beg);
        file.read((char *)pixels + bufOfs, stripByteCounts[i]);
        bufOfs += stripByteCounts[i];
        if (file.gcount() != stripByteCounts[i])
        {
            if (errorMsg) *errorMsg = boost::str(boost::format("The file is incomplete: pixel data in strip %d is too short. Expected %d bytes, but read only %d.") % i % stripByteCounts[i] % file.gcount());
            free(pixels);
            return 0;
        }
    }

    if ((pixFmt == PIX_MONO16 || pixFmt == PIX_RGB48) && isBE)
        SwapBufferWords((uint16_t *)pixels, imgWidth*imgHeight*GetBytesPerPixel(pixFmt)/2);

    if (photometricInterpretation == PHMET_WHITE_IS_ZERO)
    {
        // Reverse the values so that "black" is zero, "white" is 255 or 65535.
        if (pixFmt == PIX_MONO8)
            NegateGrayscale8((uint8_t *)pixels, imgWidth*imgHeight);
        else if (pixFmt == PIX_MONO16)
            NegateGrayscale16((uint16_t *)pixels, imgWidth*imgHeight);
    }

    file.close();

    return pixels;
}

bool GetTiffDimensions(const char *fileName, unsigned &imgWidth, unsigned &imgHeight)
{
    std::ifstream file(fileName, std::ios_base::binary);

    if (file.fail())
        return false;

    TiffHeader_t tiffHeader;
    file.read((char *)&tiffHeader, sizeof(tiffHeader));
    if (file.gcount() != sizeof(tiffHeader))
        return false;

    bool isBE = tiffHeader.id == MOTOROLA_BYTE_ORDER; // true if the file has big endian data

    if (SWAP16cnd(tiffHeader.version, isBE) != TIFF_VERSION)
        return false;

    // Seek to the first TIFF directory
    file.seekg(SWAP32cnd(tiffHeader.dirOffset, isBE), std::ios_base::beg);

    uint16_t numDirEntries;
    file.read((char *)&numDirEntries, sizeof(numDirEntries));
    numDirEntries = SWAP16cnd(numDirEntries, isBE);
    if (file.gcount() != sizeof(numDirEntries))
        return false;

    imgWidth = imgHeight = -1;

    std::fstream::pos_type nextFieldPos = file.tellg();
    for (unsigned i = 0; i < numDirEntries; i++)
    {
        TiffField_t tiffField;

        file.seekg(nextFieldPos, std::ios_base::beg);
        file.read((char *)&tiffField, sizeof(tiffField));
        if (file.gcount() != sizeof(tiffField))
            return false;
        nextFieldPos = file.tellg();

        tiffField.tag = SWAP16cnd(tiffField.tag, isBE);
        tiffField.type = SWAP16cnd(tiffField.type, isBE);
        tiffField.count = SWAP32cnd(tiffField.count, isBE);
        if (tiffField.count > 1 || tiffField.type == ttDWord)
            tiffField.value = SWAP32cnd(tiffField.value, isBE);
        else if (tiffField.count == 1 && tiffField.type == ttWord)
            tiffField.value = SWAP16in32cnd(tiffField.value, isBE);

        switch (tiffField.tag)
        {
        case TAG_IMAGE_WIDTH: imgWidth = tiffField.value; break;
        case TAG_IMAGE_HEIGHT: imgHeight = tiffField.value; break;
        }

        if (imgWidth != -1 && imgHeight != -1)
            break;
    }

    return true;
}

} // namespace TIFF

/// Converts data in input buffer to the specified pixel format and writes it to the destination buffer; if formats are the same, does nothing
void ConvertPixelFormat(
    void *srcBuf,             ///< Source (input) buffer (pixels stored left to right, top to bottom, no padding)
    void *destBuf,            ///< Destination (output) buffer
    int width,                ///< Image width (number of columns in the buffers)
    int height,               ///< Image height (number of rows in the buffers)
    PixelFormat_t srcPixFmt,  ///< Pixel format in 'srcBuf'
    PixelFormat_t destPixFmt, ///< Desired pixel format in 'destBuf'; PIX_PAL8 is not supported
    uint8_t palette[]         ///< Pointer to 256-element RGB(+1) palette (256*4 bytes); required if 'srcPixFmt' or 'destPixFmt' is PIX_PAL8
)
{
    if (srcPixFmt == PIX_UNCHANGED || destPixFmt == PIX_UNCHANGED)
        throw std::runtime_error("ConvertPixelFormat(): specifying PIX_UNCHANGED is not allowed.");
    if (destPixFmt == PIX_PAL8 && srcPixFmt != PIX_PAL8)
        throw std::runtime_error("ConvertPixelFormat(): cannot convert to PIX_PAL8");
    if (srcPixFmt == PIX_PAL8 && !palette)
        throw std::runtime_error("ConvertPixelFormat(): palette required when converting from PIX_PAL8");

    if (srcPixFmt == destPixFmt)
        return;

    uint8_t *inpPtr = (uint8_t*)srcBuf,
            *outPtr = (uint8_t*)destBuf;

    int inpPtrStep = GetBytesPerPixel(srcPixFmt),
        outPtrStep = GetBytesPerPixel(destPixFmt);

    for (int i = 0; i < width*height; i++)
    {
        if (srcPixFmt == PIX_MONO8)
        {
            uint8_t src = *inpPtr;
            switch (destPixFmt)
            {
            case PIX_MONO16: *(uint16_t *)outPtr = (uint16_t)src << 8; break;
            case PIX_MONO32F: *(float *)outPtr = src * 1.0f/0xFF; break;
            case PIX_RGB24: outPtr[0] = outPtr[1] = outPtr[2] = src; break;
            case PIX_RGB48:
                ((uint16_t *)outPtr)[0] =
                    ((uint16_t *)outPtr)[1] =
                    ((uint16_t *)outPtr)[2] =  (uint16_t)src << 8;
                break;
            }
        }
        else if (srcPixFmt == PIX_MONO16)
        {
            uint16_t src = *(uint16_t *)inpPtr;
            switch (destPixFmt)
            {
            case PIX_MONO8: *outPtr = (uint8_t)(src >> 8); break;
            case PIX_MONO32F: *(float *)outPtr = src * 1.0f/0xFFFF; break;
            case PIX_RGB24: outPtr[0] = outPtr[1] = outPtr[2] = (uint8_t)(src >> 8); break;
            case PIX_RGB48:
                ((uint16_t *)outPtr)[0] =
                    ((uint16_t *)outPtr)[1] =
                    ((uint16_t *)outPtr)[2] = src;
                break;
            }
        }
        else if (srcPixFmt == PIX_MONO32F)
        {
            float src = *(float *)inpPtr;
            switch (destPixFmt)
            {
            case PIX_MONO8: *outPtr = (uint8_t)(src * 0xFF); break;
            case PIX_MONO16: *(uint16_t *)outPtr = (uint16_t)(src * 0xFFFF); break;
            case PIX_RGB24: outPtr[0] = outPtr[1] = outPtr[2] = (uint8_t)(src * 0xFF); break;
            case PIX_RGB48:
                ((uint16_t *)outPtr)[0] =
                    ((uint16_t *)outPtr)[1] =
                    ((uint16_t *)outPtr)[2] = (uint16_t)(src * 0xFFFF); break;
            }
        }
        // When converting from a color format to mono, use sum (scaled) of all channels as the pixel brightness.
        else if (srcPixFmt == PIX_PAL8)
        {
            uint8_t src = *inpPtr;
            switch (destPixFmt)
            {
            case PIX_MONO8: *outPtr = (uint8_t)(((int)palette[4*src] + palette[4*src+1] + palette[4*src+2])/3); break;
            case PIX_MONO16: *(uint16_t *)outPtr = ((uint16_t)palette[4*src] + palette[4*src+1] + palette[4*src+2])/3; break;
            case PIX_MONO32F: *(float *)outPtr = ((int)palette[4*src] + palette[4*src+1] + palette[4*src+2]) * 1.0f/(3*0xFF); break;
            case PIX_RGB24:
                outPtr[0] = palette[4*src];
                outPtr[1] = palette[4*src+1];
                outPtr[2] = palette[4*src+2];
                break;
            case PIX_RGB48:
                ((uint16_t *)outPtr)[0] = (uint16_t)palette[4*src] << 8;
                ((uint16_t *)outPtr)[1] = (uint16_t)palette[4*src+1] << 8;
                ((uint16_t *)outPtr)[2] = (uint16_t)palette[4*src+2] << 8;
                break;
            }
        }
        else if (srcPixFmt == PIX_RGB24)
        {
            switch (destPixFmt)
            {
            case PIX_MONO8: *outPtr = (uint8_t)(((int)inpPtr[0] + inpPtr[1] + inpPtr[2])/3); break;
            case PIX_MONO16: *(uint16_t *)outPtr = ((uint16_t)inpPtr[0] + inpPtr[1] + inpPtr[2])/3; break;
            case PIX_MONO32F: *(float *)outPtr = ((int)inpPtr[0] + inpPtr[1] + inpPtr[2]) * 1.0f/(3*0xFF); break;
            case PIX_RGB48:
                ((uint16_t *)outPtr)[0] = (uint16_t)inpPtr[0] << 8;
                ((uint16_t *)outPtr)[1] = (uint16_t)inpPtr[1] << 8;
                ((uint16_t *)outPtr)[2] = (uint16_t)inpPtr[2] << 8;
                break;
            }
        }
        else if (srcPixFmt == PIX_RGB48)
        {
            uint16_t *inpPtr16 = (uint16_t *)inpPtr;
            switch (destPixFmt)
            {
            case PIX_MONO8: *outPtr = (uint8_t)(((int)inpPtr16[0] + inpPtr16[1] + inpPtr16[2])/3); break;
            case PIX_MONO16: *(uint16_t *)outPtr = (uint16_t)(((int)inpPtr16[0] + inpPtr16[1] + inpPtr16[2])/3); break;
            case PIX_MONO32F: *(float *)outPtr = ((int)inpPtr16[0] + inpPtr16[1] + inpPtr16[2]) * 1.0f/(3*0xFFFF); break;
            case PIX_RGB24:
                outPtr[0] = (uint8_t)(inpPtr16[0] >> 8);
                outPtr[1] = (uint8_t)(inpPtr16[1] >> 8);
                outPtr[2] = (uint8_t)(inpPtr16[2] >> 8);
                break;
            }
        }

        inpPtr += inpPtrStep;
        outPtr += outPtrStep;
    }
}

/// Reads an image and returns pointer to the newly allocated buffer with pixel contents or 0 on error
void *ReadImageFile(std::string fileName,
        PixelFormat_t destPixFmt, ///< Desired pixel format of data in the returned buffer
        int &imgWidth,  ///< Receives image width
        int &imgHeight, ///< Receives image height
        PixelFormat_t *receivedPixFmt, ///< If not null and destPixFmt==PIX_UNCHANGED, receives the pixel format of the source image
        uint8_t palette[],     ///< If not null and reading from an 8-bit file with palette, receives the palette (1024 bytes)
        std::string *errorMsg  ///< If not null, receives error message (if any)
)
{
    std::string ext = fileName.substr(fileName.find_last_of('.'));
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    void *pixels = 0;

    PixelFormat_t dummy; // Used when the caller passed null as 'receivedPixFmt'
    if (receivedPixFmt == 0) receivedPixFmt = &dummy;

    uint8_t localPalette[256*4];
    uint8_t *palPtr = (palette != 0 ? palette : localPalette);

    if (ext == ".bmp")
        pixels = BMP::ReadBmp(fileName.c_str(), imgWidth, imgHeight, *receivedPixFmt, palPtr);
    else if (ext == ".tif" || ext == ".tiff")
        pixels = TIFF::ReadTiff(fileName.c_str(), *receivedPixFmt, imgWidth, imgHeight, errorMsg);
    else
        return 0;

    if (pixels && destPixFmt != PIX_UNCHANGED && destPixFmt != *receivedPixFmt)
    {
        void *convertedPixels = malloc(imgWidth * imgHeight * GetBytesPerPixel(destPixFmt));
        ConvertPixelFormat(pixels, convertedPixels, imgWidth, imgHeight, *receivedPixFmt, destPixFmt, palPtr);
        free(pixels);
        return convertedPixels;
    }
    else
        return pixels;
}

bool GetImageDimensions(std::string fileName, unsigned &imgWidth, unsigned &imgHeight)
{
    std::string ext = fileName.substr(fileName.find_last_of('.'));
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == ".bmp")
        return BMP::GetBmpDimensions(fileName.c_str(), imgWidth, imgHeight);
    else if (ext == ".tif" || ext == ".tiff")
        return TIFF::GetTiffDimensions(fileName.c_str(), imgWidth, imgHeight);
    else
        return false;
}

/// Saves image; returns 'false' on error
bool SaveImageFile(std::string fileName, ///< Output file name
        int imgWidth,              ///< Image width
        int imgHeight,             ///< Image height
        PixelFormat_t pixFmt,      ///< Pixel format
        void *pixels,              ///< Pixel contents in 'pixFmt' format
        uint8_t palette[]          ///< Points to the palette (1024 bytes) to be saved if pixFmt is PIX_PAL8;
)
{
    std::string ext = fileName.substr(fileName.find_last_of('.'));
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext == ".bmp")
        return BMP::SaveBmp(fileName.c_str(), imgWidth, imgHeight, pixFmt, pixels, palette);
    else if (ext == ".tif" || ext == ".tiff")
        return TIFF::SaveTiff(fileName.c_str(), pixels, pixFmt, imgWidth, imgHeight);
    else
        return false;
}


inline float ClampLuminance(float val, float maxVal)
{
    if (val < 0.0f)
        return 0.0f;
    else if (val > maxVal)
        return maxVal;
    else
        return val;
}

/// Cubic (Hermite) interpolation of 4 subsequent values fm1, f0, f1, f2 at location 0<=t<=1 between the middle elements (f0 and f1)
template<typename T>
inline float InterpolateCubic(float t, T fm1, T f0, T f1, T f2)
{
    float delta_k = (float)f1 - f0;
    float dk = ((float)f1 - fm1)*0.5f, dk1 = ((float)f2 - f0)*0.5f;

    float a0 = f0, a1 = dk, a2 = 3.0f*delta_k - 2.0f*dk - dk1,
        a3 = (float)dk + dk1 - 2.0f*delta_k;

    return t*(t*(a3*t + a2)+a1)+a0;
}

/// Performs sub-pixel image translation using cubic interpolation
template<typename Lum_t>
void SubpixelTranslationImpl(
        Lum_t *src,  ///< Source image
        Lum_t *dest, ///< Destination image
        int width,   ///< Image width
        int height,  ///< Image height
        int numChannels, ///< Number of channels (number of 'Lum_t' samples per pixel)
        float maxLum, ///< Max value to clamp the output values to
        float dx,   ///< Translation in X, |dx| < 1
        float dy    ///< Translation in Y, |dy| < 1
)
{
    if (fabs(dx) >= 1.0 || fabs(dy) >= 1.0)
        throw std::runtime_error("SubpixelTranslation(): the translation must be by less than 1 pixel in each direction.");

    int idx = dx < 0.0f ? 1 : -1;
    int idy = dy < 0.0f ? 1 : -1;

    dx = fabs(dx);
    dy = fabs(dy);

    // Skip 2-pixels borders on each side of the image
    #pragma omp parallel for
    for (int row = 2; row < height-2; row++)
    {
        for (int col = 2; col < width-2; col++)
        {
            for (int ch = 0; ch < numChannels; ch++)
            {
                float yvals[4];

                // Perform 4 interpolations at 4 adjacent rows, using X offsets -1, 0, 1, 2 (*idx)
                int y = row - idy;
                for (int relY = -1; relY <= 2; relY++)
                {
                    yvals[relY+1] = InterpolateCubic(dx,
                                     src[(col-idx     + y*width)*numChannels + ch],
                                     src[(col         + y*width)*numChannels + ch],
                                     src[(col+idx     + y*width)*numChannels + ch],
                                     src[(col+idx+idx + y*width)*numChannels + ch]);
                    y += idy;
                }

                // Perform the final vertical (column) interpolation of the 4 horizontal (row) values interpolated previously
                dest[(col + row*width)*numChannels + ch] = (Lum_t)ClampLuminance(InterpolateCubic(dy, yvals[0], yvals[1], yvals[2], yvals[3]), maxLum);
            }
        }
    }

    // Copy the 2-pixel borders without changes

    // 2 top rows
    memcpy(dest + (0 + 0*width) * numChannels, src + (0 + 0*width) * numChannels, width * numChannels * sizeof(Lum_t));
    memcpy(dest + (0 + 1*width) * numChannels, src + (0 + 1*width) * numChannels, width * numChannels * sizeof(Lum_t));
    // 2 bottom rows
    memcpy(dest + (0 + (height-1)*width) * numChannels, src + (0 + (height-1)*width) * numChannels, width * numChannels * sizeof(Lum_t));
    memcpy(dest + (0 + (height-2)*width) * numChannels, src + (0 + (height-2)*width) * numChannels, width * numChannels * sizeof(Lum_t));
    // 2 leftmost and 2 rightmost columns
    for (int row = 0; row < height; row++)
    {
        // 2 leftmost columns
        memcpy(dest + (0 + row*width) * numChannels, src + (0 + row*width) * numChannels,
                2  *numChannels * sizeof(Lum_t)); // copying 2 pixels, each is 'numChannels' elements
        // 2 rightmost columns
        memcpy(dest + (width-2 + row*width) * numChannels, src + (width-2 + row*width) * numChannels,
                2 * numChannels * sizeof(Lum_t)); // copying 2 pixels, each is 'numChannels' elements
    }
}

/// Performs a sub-pixel translation of an image using bicubic interpolation; PIX_PAL8 pixel format is not supported
void SubpixelTranslation(
    void *src,  ///< Source image
    void *dest, ///< Destination image
    int width,  ///< Image width
    int height, ///< Image height
    PixelFormat_t pixFmt, ///< Pixel format; PIX_PAL8 is not supported
    float dx,   ///< Translation in X, |dx| < 1
    float dy    ///< Translation in Y, |dy| < 1
)
{
    if (pixFmt == PIX_PAL8)
        throw std::runtime_error("SubpixelTranslation(): PIX_PAL8 is not supported.");

    switch (pixFmt)
    {
    case PIX_MONO8:   SubpixelTranslationImpl<uint8_t> ( (uint8_t *)src,  (uint8_t *)dest, width, height, 1,   (float)0xFF, dx, dy); break;
    case PIX_MONO16:  SubpixelTranslationImpl<uint16_t>((uint16_t *)src, (uint16_t *)dest, width, height, 1, (float)0xFFFF, dx, dy); break;
    case PIX_MONO32F: SubpixelTranslationImpl<float>   (   (float *)src,    (float *)dest, width, height, 1,          1.0f, dx, dy); break;
    case PIX_RGB24:   SubpixelTranslationImpl<uint8_t> ( (uint8_t *)src,  (uint8_t *)dest, width, height, 3,   (float)0xFF, dx, dy); break;
    case PIX_RGB48:   SubpixelTranslationImpl<uint16_t>((uint16_t *)src, (uint16_t *)dest, width, height, 3, (float)0xFFFF, dx, dy); break;
    }
}

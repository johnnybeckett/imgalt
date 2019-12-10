/*
imgalt - Image Alignment Tool
Author: GreatAttractor

version 0.5
2014/05/22

This code can be freely distributed and used for any purpose.

File description:
    Utility functions header.

*/
#ifndef COMP_H_
#define COMP_H_

#include <boost/cstdint.hpp>

/// Pixel format for image read/write operations
typedef enum {
    PIX_UNCHANGED, ///< Use the source pixel format
    PIX_PAL8, ///< 8-bit with palette (can be a grayscale palette)
    PIX_MONO8, ///< 8-bit grayscale
    PIX_RGB24, ///< 24-bit RGB (8 bits per channel)
    PIX_MONO16, ///< 16-bit grayscale
    PIX_RGB48, ///< 48-bit RGB (16 bits per channel)
    PIX_MONO32F ///< 32-bit floating point grayscale
} PixelFormat_t;

int GetBytesPerPixel(PixelFormat_t pixFmt);

/// Returns the smallest power of 2 which is > n
unsigned GetClosestGPowerOf2(unsigned n);

/// Returns the greatest power of 2 which is <= n
unsigned GetClosestLEPowerOf2(unsigned n);

/// Reads an image and returns pointer to the newly allocated buffer with pixel contents or 0 on error
void *ReadImageFile(std::string fileName,
        PixelFormat_t destPixFmt, ///< Desired pixel format of data in the returned buffer
        int &imgWidth,  ///< Receives image width
        int &imgHeight, ///< Receives image height
        PixelFormat_t *receivedPixFmt = 0, ///< If not null and destPixFmt==PIX_UNCHANGED, receives the pixel format of the source image
        uint8_t palette[] = 0,     ///< If not null and reading from an 8-bit file with palette, receives the palette (1024 bytes)
        std::string *errorMsg = 0  ///< If not null, receives error message (if any)
);

bool GetImageDimensions(std::string fileName, unsigned &imgWidth, unsigned &imgHeight);

/// Saves image; returns 'false' on error
bool SaveImageFile(std::string fileName, ///< Output file name
        int imgWidth,              ///< Image width
        int imgHeight,             ///< Image height
        PixelFormat_t pixFmt,      ///< Pixel format
        void *pixels,              ///< Pixel contents in 'pixFmt' format
        uint8_t palette[]          ///< Points to the palette (1024 bytes) to be saved if pixFmt is PIX_PAL8;
);


/// Multiplies input image by window function
void ApplyWindowFunction(float img[], ///< Input image, size*size elements
        float windowFunc[], ///< Window function values, size*size elements
        float dest[], ///< Output buffer, size*size elements
        int size ///< Number of rows and columns in each array
);

/// Calculates window function and writes its values into 'buf' (a square array)
void CalcWindowFunction(
        int wndSize, ///< Window size
        float buf[] ///< Destination buffer, wndSize*wndSize elements
);

/// Resizes and translates image by cropping and/or padding (with zeros) to the destination size and offset (there is no scaling)
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
);

/// Blurs 'src' image using a 5x5 Gaussian kernel and writes the result to 'dest'; 2-pixel borders are not blurred
void BlurImage(
    uint8_t *src, ///< Source image (8-bit luminance)
    int width,    ///< Image width
    int height,   ///< Image height
    uint8_t *dest ///< Destination image (same dimensions as 'src')
);

/// Calculates squared gradient lengths in the source image; 3 pixel-wide borders are skipped
void CalcGradients(
    uint8_t *src,  ///< Source image
    int width,     ///< Image width
    int height,    ///< Image height
    uint32_t *dest ///< Destination buffer (same dimensions as 'src')
);

/// Converts data in input buffer to the specified pixel format and writes it to the destination buffer; if formats are the same, does nothing
void ConvertPixelFormat(
    void *srcBuf,             ///< Source (input) buffer (pixels stored left to right, top to bottom, no padding)
    void *destBuf,            ///< Destination (output) buffer
    int width,                ///< Image width (number of columns in the buffers)
    int height,               ///< Image height (number of rows in the buffers)
    PixelFormat_t srcPixFmt,  ///< Pixel format in 'srcBuf'
    PixelFormat_t destPixFmt, ///< Desired pixel format in 'destBuf'; PIX_PAL8 is not supported
    uint8_t palette[]         ///< Pointer to 256-element RGB(+1) palette (256*4 bytes); required if 'srcPixFmt' or 'destPixFmt' is PIX_PAL8
);

/// Performs a sub-pixel translation of an image using bicubic interpolation; PIX_PAL8 pixel format is not supported
void SubpixelTranslation(
    void *src,  ///< Source image
    void *dest, ///< Destination image
    int width,  ///< Image width
    int height, ///< Image height
    PixelFormat_t pixFmt, ///< Pixel format; PIX_PAL8 is not supported
    float dx,   ///< Translation in X, |dx| < 1
    float dy    ///< Translation in Y, |dy| < 1
);

#endif /* COMP_H_ */

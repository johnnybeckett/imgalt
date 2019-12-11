/*
imgalt - Image Alignment Tool
Author: GreatAttractor

version 0.5
2014/05/22

This code can be freely distributed and used for any purpose.

File description:
    Main program file.

*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#if defined(_OPENMP)
  #include <omp.h>
#else
  #include "omp_stubs.h"
#endif
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/chrono.hpp>
#include "fft.h"
#include "comp.h"
#include "types.h"
namespace bfs = boost::filesystem;
namespace bch = boost::chrono;

const char *VERSION_STRING = "0.5";
const char *DATE_STRING = "2014/05/22";

namespace Options
{
    const char *INPUT_DIR = "--input-dir";
    const char *OUTPUT_DIR = "--output-dir";
    const char *HELP = "--help";
    const char *VERBOSE = "--verbose";
    const char *THREADS = "--threads";
    const char *NO_CROP = "--no-crop";
    const char *NO_SUBPIXEL = "--no-subpixel";
}

namespace Log
{
typedef enum { QUIET = 0, NORMAL, VERBOSE} LogLevel_t;

LogLevel_t loggingLevel = NORMAL;

/// Prints a message. Newline is NOT added by default.
void Print(LogLevel_t level, std::string msg)
{
    if (level <= loggingLevel)
    {
        std::cout << msg;
        std::cout.flush();
    }
}
}

namespace Vars
{
/// Input directory (default: current)
std::string inputDir = ".";
/// Input directory (default: current)
std::string outputDir = ".";
/// Number of threads in OpenMP parallel regions (default: number of logical CPUs detected by OpenMP)
int numThreads = omp_get_num_procs();
/// If 'true', aligned images are cropped to their intersection; otherwise they are are padded to the common bounding box size
bool cropImages = true;
/// If 'true', sub-pixel alignment is enabled
bool subpixelAlignment = true;
}

/// Lists all BMP images in 'inputDir' sorted by name
void ListImageFiles(std::string inputDir, std::vector<std::string> &list)
{
    list.clear();
    try
    {
        bfs::path p(inputDir);
        if (!bfs::exists(p))
        {
            std::cout << inputDir << " does not exist.\n";
            return;
        }

        if (!bfs::is_directory(p))
        {
            std::cout << inputDir << " is not a directory.\n";
            return;
        }

        bfs::directory_iterator dirIter(p);
        while (dirIter != bfs::directory_iterator())
        {
            if (bfs::is_regular_file(dirIter->path()))
            {
                std::string ext = dirIter->path().extension().string();
                std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
                if (ext == ".bmp" || ext == ".tif" || ext == ".tiff")
                    list.push_back(dirIter->path().string());
            }

            dirIter++;
        }
        std::sort(list.begin(), list.end());
    }
    catch (const bfs::filesystem_error &ex)
    {
        std::cout << ex.what() << std::endl;
    }
}

/// Returns maximum of width and height of all BMP and TIFF images in 'fileNames' (returns 0 on error)
unsigned GetMaxImageDimension(std::vector<std::string> &fileNames,
        std::vector<Point_t> &dimensions ///< Receives dimensions of files from 'fileNames'
)
{
    unsigned result = 0;
    dimensions.clear();
    for (unsigned i = 0; i < fileNames.size(); i++)
    {
        unsigned width, height;
        if (! GetImageDimensions(fileNames[i], width, height))
        {
            std::cout << "Cannot read " << fileNames[i] << ".\n";
            return 0;
        }
        if (width > result)
            result = width;
        if (height > result)
            result = height;

        dimensions.push_back(Point_t(width, height));
    }

    return result;
}

void PrintHelp()
{
    std::cout <<
    "Aligns all BMP (8- or 24-bit) and TIFF (8- or 16-bit, mono or RGB, not compressed) images (sorted by name) in the input directory (default: current) and saves them to output directory (default: current) under new names.\n"
    "Usage:\n\n"

    "imgalt [options] [[" << Options::INPUT_DIR << "] <input directory>]\n\n" <<

    " options:\n"
    "   " << Options::OUTPUT_DIR << " <directory>: output directory\n"
    "   " << Options::NO_CROP << ": do not crop to the intersection; instead, pad to bounding box\n"
    "   " << Options::THREADS << " <count>: number of threads to use\n"
    "   " << Options::NO_SUBPIXEL << ": disable sub-pixel alignment\n"
    "   " << Options::VERBOSE << ": print additional information during processing\n\n";
}

bool ParseCommandLineOptions(int argc, char *argv[])
{
    for (int i = 1; i < argc; i++)
    {
        if (0 == strcmp(argv[i], Options::INPUT_DIR))
        {
            Vars::inputDir = argv[++i];
        }
        else if (0 == strcmp(argv[i], Options::OUTPUT_DIR))
        {
            Vars::outputDir = argv[++i];
        }
        else if (0 == strcmp(argv[i], Options::HELP))
        {
            PrintHelp();
            return false;
        }
        else if (0 == strcmp(argv[i], Options::VERBOSE))
        {
            Log::loggingLevel = Log::VERBOSE;
        }
        else if (0 == strcmp(argv[i], Options::THREADS))
        {
            char *endp;
            i++;
            Vars::numThreads = strtol(argv[i], &endp, 10);
            if (Vars::numThreads == 0 || endp == argv[i])
            {
                std::cout << Options::THREADS << ": expected a positive integer value.\n";
                return false;
            }
        }
        else if (0 == strcmp(argv[i], Options::NO_CROP))
        {
            Vars::cropImages = false;
        }
        else if (0 == strcmp(argv[i], Options::NO_SUBPIXEL))
        {
            Vars::subpixelAlignment = false;
        }
        else // default parameter: input directory
        {
            Vars::inputDir = argv[i];
        }
    }

    return true;
}

void PrintInfo()
{
    std::cout << "Image Alignment Tool\n"
                 "(c) 2014 GreatAttractor\n"
                 "v. " << VERSION_STRING << " (" << DATE_STRING << ")\n"
                 "Free for all uses.\n\n";
}

/// Prints the time interval from 'tstart' to now and increases 'microseconds' (if not null) accordingly
void LogTimeInterval(bch::high_resolution_clock::time_point &tstart, bool newline = false, unsigned *microseconds = 0)
{
    std::string formatStr = "%.1f ms";
    if (newline)
        formatStr += "\n";

    unsigned dur = bch::duration_cast<bch::microseconds>(bch::high_resolution_clock::now() - tstart).count();
    Log::Print(Log::VERBOSE, boost::str(boost::format(formatStr) % (dur*0.001f)));

    if (microseconds)
        *microseconds += dur;
}


/// Determines (using phase correlation) the vector by which image 2 (given by its discrete Fourier transform 'img2FFT') is translated w.r.t. image 1 ('img1FFT')
RealPoint_t DetermineImageTranslation(
    unsigned N,  ///< FFT size (number of rows and columns in input images)
    std::complex<float> *img1FFT, ///< FFT of the first image (N*N elements)
    std::complex<float> *img2FFT, ///< FFT of the second image (N*N elements)
    bool subpixelAccuracy, ///< If 'true', the translation is determined down to sub-pixel accuracy
    unsigned *totalInvFFTmicroseconds ///< If not null, gets increased by the inv. FFT time
)
{
/*
    NOTE: If we do not want sub-pixel accuracy, we want this function to return an integer translation vector (i.e. no fractional part in X or Y).
          Otherwise the caller, DetermineTranslationVectors(), would accumulate the fractional parts in the translation vectors array, but since
          we wouldn't use sub-pixel positioning of output images, there would be a ragged 1-pixel jittering in the output sequence.

          In other words, we must either: return fractional vectors from here and use sub-pixel positioning of output images,
          or: return integer vectors from here and do only integer translations of output images.

          Using (accumulating) fractional vectors and then performing integer-only image translations gives poor results.
 */


    bch::high_resolution_clock::time_point tstart;


    // Cross-power spectrum
    std::complex<float> *cps = new std::complex<float>[N*N];
    // Cross-correlation
    std::complex<float> *cc = new std::complex<float>[N*N];

    Log::Print(Log::VERBOSE, std::string(" Cross-power spectrum: "));
    tstart = bch::high_resolution_clock::now();
    CalcCrossPowerSpectrum2D(img1FFT, img2FFT, cps, N);
    LogTimeInterval(tstart, true);

    Log::Print(Log::VERBOSE, std::string(" Cross-correlation: "));
    tstart = bch::high_resolution_clock::now();
    CalcFFTinv2D(cps, N, cc);
    LogTimeInterval(tstart, true, totalInvFFTmicroseconds);

    // Find the highest-Re element in cross-correlation array
    Log::Print(Log::VERBOSE, std::string(" Identifying the CC peak: "));
    tstart = bch::high_resolution_clock::now();
    int maxx = 0, maxy = 0;
    float maxval = 0.0f;
    for (unsigned y = 0; y < N; y++)
        for (unsigned x = 0; x < N; x++)
        {
            float currval = cc[x + y*N].real();
            if (currval > maxval)
            {
                maxval = currval;
                maxx = x;
                maxy = y;
            }
        }
    LogTimeInterval(tstart, true);

    // Infer the translation vector from the 'cc's largest element's indices

    int Tx, Ty; // prev->curr translation vector
    if (maxx < N/2)
        Tx = maxx;
    else
        Tx = maxx-N;

    if (maxy < N/2)
        Ty = maxy;
    else
        Ty = maxy-N;

    float subdx = 0.0f, subdy = 0.0f; // subpixel translation vector
    if (subpixelAccuracy)
    {
        //   Subpixel translation detection based on:
        //
        //     "Extension of Phase Correlation to Subpixel Registration"
        //     Hassan Foroosh, Josiane B. Zerubia, Marc Berthod
        //
        #define CLAMP(k) (((k)+N)%N)

        float ccXhi = cc[CLAMP(maxx+1) + maxy*N].real(),
              ccXlo = cc[CLAMP(maxx-1) + maxy*N].real(),
              ccYhi = cc[maxx + CLAMP(maxy+1)*N].real(),
              ccYlo = cc[maxx + CLAMP(maxy-1)*N].real(),
              ccPeak = cc[maxx + maxy*N].real();

        if (ccXhi > ccXlo)
        {
            float dx1 = ccXhi/(ccXhi + ccPeak),
                  dx2 = ccXhi/(ccXhi - ccPeak);

            if (dx1 > 0 && dx1 < 1.0f)
                subdx = dx1;
            else if (dx2 > 0 && dx2 < 1.0f)
                subdx = dx2;
        }
        else
        {
            float dx1 = ccXlo/(ccXlo + ccPeak),
                  dx2 = ccXlo/(ccXlo - ccPeak);

            if (dx1 > 0 && dx1 < 1.0f)
                subdx = -dx1;
            else if (dx2 > 0 && dx2 < 1.0f)
                subdx = -dx2;
        }

        if (ccYhi > ccYlo)
        {
            float dy1 = ccYhi/(ccYhi + ccPeak),
                  dy2 = ccYhi/(ccYhi - ccPeak);

            if (dy1 > 0 && dy1 < 1.0f)
                subdy = dy1;
            else if (dy2 > 0 && dy2 < 1.0f)
                subdy = dy2;
        }
        else
        {
            float dy1 = ccYlo/(ccYlo + ccPeak),
                  dy2 = ccYlo/(ccYlo - ccPeak);

            if (dy1 > 0 && dy1 < 1.0f)
                subdy = -dy1;
            else if (dy2 > 0 && dy2 < 1.0f)
                subdy = -dy2;
        }
    }

    free(cps);
    free(cc);

    return RealPoint_t(Tx + subdx, Ty + subdy);
}

/// Determines translation vectors of an image sequence
bool DetermineTranslationVectors(
        unsigned N, ///< FFT size
        std::vector<std::string> &inputFiles, ///< List of input file names
        /// Receives list of translation vectors between files in 'inputFiles'; each vector is a translation relative to the first image
        std::vector<RealPoint_t> &translation,
        /// Receives the bounding box (within the NxN working area) of all images after alignment
        Rectangle_t &bBox
        // (Note: an untranslated image starts in the working area at (N-imgWidth)/2, (N-imgHeight)/2)
)
{
    bch::high_resolution_clock::time_point tstart; // Starting point used for timing various actions
    unsigned totalFFTMicroseconds = 0, totalFFTinvMicroseconds = 0;

    int imgWidth, imgHeight; // Dimensions of the recently read image

    float *windowFunc = new float[N*N]; // values of window function
    CalcWindowFunction(N, windowFunc);
    // Window function smoothly varies from 0 at the array boundaries to 1 at the center and is used to
    // "blunt" the image, starting from the edges. Without it they would produce prominent
    // false peaks in the cross-correlation (as any sudden change in brightness generates
    // lots if high frequencies after FFT), making it very hard or impossible to detect
    // the true peak which corresponds to the actual image translation.

    float *prevImg = new float[N*N]; // previous image in the sequence (padded to N*N pixels and with window func. applied)
    float *currImg = new float[N*N]; // current image in the sequence (padded to N*N pixels and with window func. applied)

    std::complex<float> *prevFFT = new std::complex<float>[N*N];
    std::complex<float> *currFFT = new std::complex<float>[N*N];

    Log::Print(Log::NORMAL, boost::str(boost::format("Processing %d images (Ctrl-C to break)...\n\n") % inputFiles.size()));

    // Read the first image and calculate its FFT

    Log::Print(Log::VERBOSE, std::string("Reading ") + inputFiles[0] + ", padding and applying window func.: ");
    tstart = bch::high_resolution_clock::now();
    std::string errorMsg;
    float *pixels = (float *)ReadImageFile(inputFiles[0], PIX_MONO32F, imgWidth, imgHeight, 0, 0, &errorMsg);
    if (!pixels)
    {
        std::cout << "Could not read " << inputFiles[0] << ". " << errorMsg << std::endl;
        return false;
    }

    ResizeAndTranslate(pixels, PIX_MONO32F, imgWidth, imgHeight, 0, 0, imgWidth-1, imgHeight-1,
            prevImg, N, N, (N-imgWidth)/2, (N-imgHeight)/2);
    free(pixels);

    ApplyWindowFunction(prevImg, windowFunc, prevImg, N);

    LogTimeInterval(tstart);

    tstart = bch::high_resolution_clock::now();
    Log::Print(Log::VERBOSE, ". FFT: ");
    CalcFFT2D(prevImg, N, prevFFT);
    LogTimeInterval(tstart, true, &totalFFTMicroseconds);


    // Iterate over the remaining images and detect their translation

    translation.clear();
    translation.push_back(RealPoint_t(0.0f, 0.0f)); // first element corresponds with the first image - no translation

    // Bounding box of all the images after alignment (in coordinates
    // of the NxN working buffer, where an untranslated image
    // starts at (N-imgWidth)/2, (N-imgHeight)/2).
    //
    // Initially corresponds with dimensions and position of the first image.
    bBox.x = (N-imgWidth)/2;
    bBox.y = (N-imgHeight)/2;

    int xmax = bBox.x + imgWidth - 1;
    int ymax = bBox.y + imgHeight - 1;

    for (unsigned i = 1; i < inputFiles.size(); i++)
    {
        Log::Print(Log::NORMAL, boost::str(boost::format("(%d of %d) %s: ") % (i+1) % inputFiles.size() % inputFiles[i].c_str()));
        std::cout.flush();

        // Total processing time for i-th image
        bch::high_resolution_clock::time_point ttotal = bch::high_resolution_clock::now();

        // Read the current file
        Log::Print(Log::VERBOSE, std::string("Reading, padding and applying window func.: "));
        tstart = bch::high_resolution_clock::now();
        float *pixels = (float *)ReadImageFile(inputFiles[i], PIX_MONO32F, imgWidth, imgHeight, 0, 0, &errorMsg);
        if (!pixels)
        {
            std::cout << "Could not read " << inputFiles[i] << ". " << errorMsg << std::endl;
            return false;
        }

        ResizeAndTranslate(pixels, PIX_MONO32F, imgWidth, imgHeight, 0, 0, imgWidth-1, imgHeight-1,
                currImg, N, N, (N-imgWidth)/2, (N-imgHeight)/2);
        free(pixels);

        ApplyWindowFunction(currImg, windowFunc, currImg, N);

        LogTimeInterval(tstart, true);

        // Calculate the current image's FFT
        Log::Print(Log::VERBOSE, std::string(" FFT: "));
        tstart = bch::high_resolution_clock::now();
        CalcFFT2D(currImg, N, currFFT);
        LogTimeInterval(tstart, true, &totalFFTMicroseconds);

        RealPoint_t T = DetermineImageTranslation(N, prevFFT, currFFT, Vars::subpixelAlignment, &totalFFTinvMicroseconds);

        Log::Print(Log::NORMAL, boost::str(boost::format(" translated by %.2f, %.2f\n") % T.x % T.y));

        Log::Print(Log::VERBOSE, " Total: "); LogTimeInterval(ttotal, true); Log::Print(Log::VERBOSE, "\n");

        RealPoint_t Tprev = translation.back();
        translation.push_back(RealPoint_t(Tprev.x + T.x, Tprev.y + T.y));

        float intTx, intTy;
        modff(translation.back().x, &intTx);
        modff(translation.back().y, &intTy);

        bBox.x = std::min(bBox.x, (int)(N-imgWidth)/2 - (int)intTx);
        bBox.y = std::min(bBox.y, (int)(N-imgHeight)/2 - (int)intTy);

        int newXmax = (N-imgWidth)/2 - (int)intTx + imgWidth - 1;
        int newYmax = (N-imgHeight)/2 - (int)intTy + imgHeight - 1;

        if (newXmax > xmax) xmax = newXmax;
        if (newYmax > ymax) ymax = newYmax;

        // Swap pointers for the next iteration
        std::swap(prevImg, currImg);
        std::swap(prevFFT, currFFT);
    }

    bBox.width = xmax - bBox.x + 1;
    bBox.height = ymax - bBox.y + 1;

    delete[] windowFunc;
    delete[] prevImg;
    delete[] currImg;
    free(prevFFT);
    free(currFFT);

    Log::Print(Log::VERBOSE,
            boost::str(boost::format("\nAverage FFT time:      %.1f ms\n"
                                       "Average inv. FFT time: %.1f ms\n")
                % (totalFFTMicroseconds*0.001f / inputFiles.size())
                % (totalFFTinvMicroseconds*0.001f / (inputFiles.size()-1))));

    return true;
}

/// Returns the set-theoretic intersection, i.e. the largest shared area, of specified images
Rectangle_t DetermineImageIntersection(
        unsigned N,    ///< Dimensions of the (square) working buffer (i.e. FFT arrays)
        const Rectangle_t &bBox, ///< Bounding box of all aligned images
        std::vector<RealPoint_t> &translation, ///< Translation vectors relative to the first image
        std::vector<Point_t> imgSize       ///< Image sizes
)
{
    // Image intersection to be returned. Coordinates are relative to the NxN working buffer,
    // where an untranslated image starts at (N-imgWidth)/2, (N-imgHeight)/2.
    Rectangle_t result;

    // Set the intersection to cover the first image
    result.x = (N - imgSize[0].x)/2;
    result.y = (N - imgSize[0].y)/2;

    int xmax = result.x + imgSize[0].x - 1;
    int ymax = result.y + imgSize[0].y - 1;

    for (unsigned i = 1; i < translation.size(); i++)
    {
        float intTx, intTy;
        modff(translation[i].x, &intTx);
        modff(translation[i].y, &intTy);

        result.x = std::max(result.x, (int)(N-imgSize[i].x)/2 - (int)intTx);
        result.y = std::max(result.y, (int)(N-imgSize[i].y)/2 - (int)intTy);

        int newXmax = (N-imgSize[i].x)/2 - (int)intTx + imgSize[i].x - 1;
        int newYmax = (N-imgSize[i].y)/2 - (int)intTy + imgSize[i].y - 1;

        if (newXmax < xmax) xmax = newXmax;
        if (newYmax < ymax) ymax = newYmax;
    }

    result.width = xmax - result.x + 1;
    result.height = ymax - result.y + 1;

    return result;
}

int main(int argc, char *argv[])
{
    PrintInfo();
    if (!ParseCommandLineOptions(argc, argv))
    {
        std::cout << "\nRun \"imgalt --help\" for options.\n";
        return 3;
    }

    std::vector<std::string> inputFiles;
    ListImageFiles(Vars::inputDir, inputFiles);
    if (inputFiles.empty())
    {
        std::cout << "No BMP or TIFF files found.\nRun \"imgalt --help\" for options.\n";
        return 0;
    }
    else if (inputFiles.size() == 1)
    {
        std::cout << "Only 1 file found, nothing to do.\n";
        return 0;
    }

    int numProcs = omp_get_num_procs();
    Log::Print(Log::NORMAL, boost::str(boost::format("Found %d logical processor(s).\n") % numProcs));
    omp_set_num_threads(Vars::numThreads);
    if (Vars::numThreads != numProcs)
#if defined(_OPENMP)
        Log::Print(Log::NORMAL, boost::str(boost::format("Using %d thread(s).\n") % Vars::numThreads));
#else
        // in absence of OpenMP we ignore the --threads parameter
        Log::Print(Log::NORMAL, "No multithreading support in this build.\n");
#endif
    Log::Print(Log::NORMAL, "\n");

    std::vector<Point_t> imageSizes;
    unsigned maxDim = GetMaxImageDimension(inputFiles, imageSizes);
    if (maxDim == 0)
        return 1;

    /// Size of arrays used for FFT and cross-correlation
    unsigned N = GetClosestGPowerOf2(maxDim);

    Log::Print(Log::VERBOSE, boost::str(boost::format("Fourier Transform size: %dx%d.\n") % N % N));

    std::vector<RealPoint_t> translation; // Accumulated translation vectors between subsequent images (i.e. translations relative to the first image)
    Rectangle_t bBox; // Bounding box of all images after alignment

    if (!DetermineTranslationVectors(N, inputFiles, translation, bBox))
        return 2;

    // Intersection of all aligned images (the largest shared area)
    Rectangle_t imgIntersection = DetermineImageIntersection(N, bBox, translation, imageSizes);

    // Iterate again over all images, load, pad to the bounding box size or crop to intersection, translate and save
    Log::Print(Log::NORMAL, "\nSaving aligned images...\n");

    int outputWidth = Vars::cropImages ? imgIntersection.width : bBox.width;
    int outputHeight = Vars::cropImages ? imgIntersection.height : bBox.height;

    for (unsigned i = 0; i < inputFiles.size(); i++)
    {
        std::cout << i+1 << "/" << inputFiles.size() << " "; std::cout.flush();

        int imgWidth, imgHeight;

        PixelFormat_t pixFmt;
        uint8_t palette[256*4];
        std::string errorMsg;
        void *pixels = ReadImageFile(inputFiles[i], PIX_UNCHANGED, imgWidth, imgHeight, &pixFmt, palette, &errorMsg);
        if (!pixels)
        {
            std::cout << "Could not read " << inputFiles[i] << ". " << errorMsg << std::endl;
            return 2;
        }

        void *outputImg = malloc(outputWidth * outputHeight * GetBytesPerPixel(pixFmt));

        Point_t translationOrigin;
        if (Vars::cropImages)
        {
            translationOrigin.x = imgIntersection.x;
            translationOrigin.y = imgIntersection.y;
        }
        else
        {
            translationOrigin.x = bBox.x;
            translationOrigin.y = bBox.y;
        }

        float dxInt, dyInt;
        float dxFrac = modff(translation[i].x, &dxInt),
              dyFrac = modff(translation[i].y, &dyInt);

        ResizeAndTranslate(pixels, pixFmt, imgWidth, imgHeight, 0, 0, imgWidth-1, imgHeight-1,
                outputImg, outputWidth, outputHeight,
                (N-imgWidth)/2 - (int)dxInt - translationOrigin.x,
                (N-imgHeight)/2 - (int)dyInt - translationOrigin.y);
        free(pixels);

        if (Vars::subpixelAlignment)
        {
            // Interpolating of indexed-color (i.e. with palette) images is not supported; convert to RGB first
            if (pixFmt == PIX_PAL8)
            {
                void *rgbPixels = malloc(outputWidth * outputHeight * GetBytesPerPixel(PIX_RGB24));
                ConvertPixelFormat(outputImg, rgbPixels, outputWidth, outputHeight, PIX_PAL8, PIX_RGB24, palette);
                free(outputImg);
                outputImg = rgbPixels;
                pixFmt = PIX_RGB24;
            }

            void *xlatedOutput = malloc(outputWidth * outputHeight * GetBytesPerPixel(pixFmt));
            SubpixelTranslation(outputImg, xlatedOutput, outputWidth, outputHeight, pixFmt, -dxFrac, -dyFrac);

            free(outputImg);
            outputImg = xlatedOutput;
        }

        bfs::path p(inputFiles[i]);
        std::string newFileName = p.stem().string() + "_aligned" + p.extension().string();
        p = bfs::path(Vars::outputDir);
        p /= newFileName;

        if (!SaveImageFile(p.string(), outputWidth, outputHeight, pixFmt, outputImg, palette))
            std::cout << "Could not save output file " << p.string() << std::endl;
        free(outputImg);
    }

    std::cout << std::endl;

    return 0;
}

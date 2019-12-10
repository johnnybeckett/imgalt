/*
imgalt - Image Alignment Tool
Author: GreatAttractor

version 0.5
2014/05/22

This code can be freely distributed and used for any purpose.

File description:
    Fast Fourier Transform functions header.

*/
#ifndef IMGALT_FFT_HEADER
#define IMGALT_FFT_HEADER

#include <complex>

/// Calculates 1-dimensional discrete Fourier transform or its inverse
void CalcFFT1D(
    unsigned char input[], ///< Input vector
    unsigned N,      ///< Number of elements, has to be a power of two
    std::complex<float> output[] ///< Output vector
);

/// Calculates 1-dimensional inverse discrete Fourier transform
void CalcFFTinv1D(
    std::complex<float> input[], ///< Input vector
    unsigned N,                  ///< Number of elements, has to be a power of two
    std::complex<float> output[]); ///< Output vector

/// Calculates 2-dimensional discrete Fourier transform
/** Uses the row-column algorithm. */
void CalcFFT2D(
    float input[], ///< Input array containing N*N elements
    unsigned N, ///< Number of rows and columns; has to be a power of two
    std::complex<float> output[] ///< Output array containing N*N elements
        );

/// Calculates 2-dimensional inverse discrete Fourier transform
/** Uses the row-column algorithm. */
void CalcFFTinv2D(
    std::complex<float> input[], ///< Input array containing N*N elements
    unsigned N, ///< Number of rows and columns; has to be a power of two
    std::complex<float> output[] ///< Output array containing N*N elements
        );

/// Calculates cross-power spectrum of two 2D discrete Fourier transforms
void CalcCrossPowerSpectrum2D(
    std::complex<float> F1[], ///< First discrete Fourier transform (N*N elements)
    std::complex<float> F2[], ///< Second discrete Fourier transform (N*N elements)
    std::complex<float> output[], ///< Cross-power spectrum of F1 and F2 (N*N elements)
    unsigned N ///< Number of rows and columns, has to be a power of 2
    );

#endif

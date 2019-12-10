/*
imgalt - Image Alignment Tool
Author: GreatAttractor

version 0.5
2014/05/22

This code can be freely distributed and used for any purpose.

File description:
    Fast Fourier Transform functions implementation.

*/
#include "fft.h"
#include <cmath>
#include <stdlib.h>
#include <boost/cstdint.hpp>
using namespace std;
using namespace boost;

const float PI = 3.1415926536f;
const complex<float> I = complex<float>(0, 1);

/// Returns floor(log2(N))
inline int quickLog2(unsigned N)
{
    int result = 0;
    while (N > 0)
    {
        N >>= 1;
        result++;
    }
    return result - 1;
}

/*
/// Quick raising to integer power
template<typename T>
inline T binpow(T x, int n)
{
    T pow = 1;
    while (n != 0)
    {
        if (n & 1)
            pow *= x;
        n /= 2;
        x *= x;
    }
    return pow;
}

uint32_t ReverseBits(uint32_t n, int width)
{
    uint32_t result = 0;
    result |= (n & 1);
    for (int i = 1; i <= width-1; i++)
    {
        result <<= 1;
        n >>= 1;
        result |= (n & 1);
    }

    return result;
}*/

// Not using fft1d_Pease at the moment, because even with precalculated twiddle factors
// and without bit reversal it's not faster than the recursive approach.

/// Calculates 1-dimensional discrete Fourier transform
/** Uses in-place Pease algorithm (iterative). */
/*template<typename InputT>
void fft1d_Pease(
    InputT input[],  ///< Input vector of length 'N'
    unsigned N,      ///< Number of elements in 'input', has to be a power of two
    std::complex<float> output[] ///< Output vector of length 'N'
)
{
    int t = quickLog2(N);

    for (unsigned i = 0; i < N; i++)
        output[ReverseBits(i, t)] = input[i];

    complex<float> w = exp(-2.0f * PI * I * (1.0f/N));

    for (int c = t-1; c >= 0; c--)
        for (int r = 0; r < 1<<(t-1); r++)
        {
            unsigned r0 = r & ((1<<c) - 1);
            unsigned r1 = r >> c;
            unsigned a0 = (r0 << (t-c)) + r1;
            unsigned a1 = a0 + (1 << (t-c-1));
            complex<float> y0 = output[a0 + 1];
            complex<float> y1 = binpow<complex<float> >(w, r1 << c) * output[a1 + 1];
            output[a0 + 1] = y0 + y1;
            output[a1 + 1] = y0 - y1;
        }
}*/


/// Calculates 1-dimensional discrete Fourier transform or its inverse (not normalized by N, caller must do this)
template<typename InputT>
void fft1d(
    InputT input[],  ///< Input vector of length 'N'
    unsigned N,      ///< Number of elements in 'input', has to be a power of two
    std::complex<float> output[], ///< Output vector of length 'N'
    int stride,       ///< Stride of the input vector
    int outputStride, ///< Stride of the output vector

    /// Pointer to the initial twiddle factor for N, i.e. exp(-2*pi*i/N) (or exp(2*pi*i/N) for inverse transform).
    /// NOTE: (twiddlePtr-1) has to point to the lower twiddle factor, i.e. exp(+-2*pi*i/(N/2))
    complex<float> *twiddlePtr
)
{
    if (N == 1)
        output[0] = input[0];
    else
    {
        fft1d(input, N/2, output, 2*stride, outputStride, twiddlePtr - 1);
        fft1d(input + stride, N/2, output + N/2*outputStride, 2*stride, outputStride, twiddlePtr - 1);

        // Initial twiddle factor
        complex<float> tfactor0 = *twiddlePtr;

        complex<float> tfactor = 1.0f;
        for (unsigned k = 0; k <= N/2 - 1; k++)
        {
            complex<float> t = output[k*outputStride];
            complex<float> h = tfactor * output[(k + (N>>1))*outputStride];
            output[k*outputStride] = t + h;
            output[(k + (N>>1))*outputStride] = t - h;
            tfactor *= tfactor0; // in effect, tfactor = exp(-2*PI*I * k/N)
        }
    }
}

void CalcTwiddleFactors(unsigned N, complex<float> table[], bool inverse)
{
    for (int n = quickLog2(N); n >= 0; n--)
    {
        table[n] = (inverse ?
            exp(2.0f * PI * I * (1.0f/N)) :
            exp(-2.0f * PI * I * (1.0f/N)));

        N >>= 1;
    }
}

/// Calculates 1-dimensional discrete Fourier transform
void CalcFFT1D(
    uint8_t input[],
    unsigned N,
    std::complex<float> output[])
{
    complex<float> *twiddleFactors = (complex<float> *)malloc((quickLog2(N) + 1) * sizeof(complex<float>));
    CalcTwiddleFactors(N, twiddleFactors, false);

    fft1d<uint8_t>(input, N, output, 1, 1, twiddleFactors + quickLog2(N));

    free(twiddleFactors);
}

/// Calculates 1-dimensional inverse discrete Fourier transform
void CalcFFTinv1D(
    complex<float> input[],
    unsigned N,
    std::complex<float> output[])
{
    complex<float> *twiddleFactors = (complex<float> *)malloc((quickLog2(N) + 1) * sizeof(complex<float>));
    CalcTwiddleFactors(N, twiddleFactors, true);
    
    fft1d(input, N, output, 1, 1, twiddleFactors + quickLog2(N));
    
    // Normalize to obtain the inverse transform
    float Ninv = 1.0f/N;
    for (unsigned k = 0; k < N; k++)
        output[k] *= Ninv;

    free(twiddleFactors);
}


/// Calculates 2-dimensional discrete Fourier transform using row-column algorithm
void CalcFFT2D(
    float input[], ///< Input array containing N*N elements
    unsigned N, ///< Number of rows and columns; has to be a power of two
    std::complex<float> output[] ///< Output array containing N*N elements
        )
{
    int k;

    complex<float> *twiddleFactors = (complex<float> *)malloc((quickLog2(N) + 1) * sizeof(complex<float>));
    CalcTwiddleFactors(N, twiddleFactors, false);

    // Calculate 1-dimensional transforms of all the rows
    std::complex<float> *fftrows = (std::complex<float> *)malloc(N*N*sizeof(std::complex<float>));
    #pragma omp parallel for
    for (k = 0; k < N; k++)
        fft1d<float>(input + k*N, N, fftrows + k*N, 1, 1, twiddleFactors + quickLog2(N));

    // Calculate 1-dimensional transforms of all columns in 'fftrows' to get the final result
    #pragma omp parallel for
    for (k = 0; k < N; k++)
        fft1d<std::complex<float> >(fftrows + k, N, output + k, N, N, twiddleFactors + quickLog2(N));

    free(twiddleFactors);
    free(fftrows);
}

/// Calculates 2-dimensional inverse discrete Fourier transform using row-column algorithm
void CalcFFTinv2D(
    std::complex<float> input[], ///< Input array containing N*N elements
    unsigned N, ///< Number of rows and columns; has to be a power of two
    std::complex<float> output[] ///< Output array containing N*N elements
        )
{
    int k;
    float Ninv = 1.0f/N;

    complex<float> *twiddleFactors = (complex<float> *)malloc((quickLog2(N) + 1) * sizeof(complex<float>));
    CalcTwiddleFactors(N, twiddleFactors, true);

    // Calculate 1-dimensional inverse transforms of all the rows
    std::complex<float> *fftrows = (std::complex<float> *)malloc(N*N*sizeof(std::complex<float>));
    #pragma omp parallel for
    for (k = 0; k < N; k++)
        fft1d<std::complex<float> >(input + k*N, N, fftrows + k*N, 1, 1, twiddleFactors + quickLog2(N));

    for (k = 0; k < N*N; k++)
        fftrows[k] *= Ninv;

    // Calculate 1-dimensional inverse transforms of all columns in 'fftrows' to get the final result
    #pragma omp parallel for
    for (k = 0; k < N; k++)
        fft1d<std::complex<float> >(fftrows + k, N, output + k, N, N, twiddleFactors + quickLog2(N));

    free(twiddleFactors);
    free(fftrows);

    for (k = 0; k < N*N; k++)
        output[k] *= Ninv;
}

/// Calculates cross-power spectrum of two 2D discrete Fourier transforms
void CalcCrossPowerSpectrum2D(
    std::complex<float> F1[], ///< First discrete Fourier transform (N*N elements)
    std::complex<float> F2[], ///< Second discrete Fourier transform (N*N elements)
    std::complex<float> output[], ///< Cross-correlation of F1 and F2 (N*N elements)
    unsigned N ///< Number of rows and columns, has to be a power of 2
    )
{
    #pragma omp parallel for
    for (int i = 0; i < N*N; i++)
    {
        output[i] = std::conj(F1[i]) * F2[i];
        float magn = std::abs(output[i]);
        if (magn > 1.0e-8f)
            output[i] /= magn;
    }
}

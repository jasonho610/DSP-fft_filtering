//
//  DIT-FFT.h
//  DSP HW6
//
//  Created by 何冠勳 on 2021/1/23.
//

#ifndef DIT_FFT_h
#define DIT_FFT_h
#include <iostream>
#include <complex>
#include <vector>
#define pi 3.14159265358979323846

using namespace std;

const complex<double> J(0,1);
void split(complex<double> *inputs, double N);
void fft(complex<double>* x, double N, bool inverse);
void ifft(complex<double>* x, double N);
/* Moves all even indices to 1st half
   odd indices to 2nd half of inputs.
   inputs : an array of complex numbers
   N : length of inputs array  */
void split(complex<double> *inputs, int N)
{
    complex<double>* even = new complex<double>[N/2];
    complex<double>* odd = new complex<double>[N/2];
    int ei = 0;
    int oi = 0;

    for (int i=0; i < N; i++)
    {
        if ((i % 2) == 0)
            even[ei++] = inputs[i];
        else
            odd[oi++] = inputs[i];
    }

    int size = N/2;
    memcpy(inputs, even, sizeof(complex<double>)*size);
    memcpy(inputs+size, odd, sizeof(complex<double>)*size);

    delete[] even;
    delete[] odd;
}

/* This function computes the fast fourier
   transform of a list of complex numbers
   of length N.
   x : input array of complex number that represent
       a sampled function amplitude
   N : the length of x must be a power of 2
*/
void fft(complex<double>* x, double N, bool inverse = false)
{
    /* base of recursion */
    if (N == 1)
        return;
    
    /* no rounding needed if N is base 2 */
    int n = N/2;

    /* set primitive root of unity */
    complex<double> wn;
    if(inverse)
        wn = exp((2*pi*J)/N);
    else
        wn = exp(((-2)*pi*J)/N);
    complex<double> w = 1;

    /* move odd and evened indexed to each half
       of array x */
    split(x, 2*n);
    
    /* even and odd */
    fft(x, n, inverse);
    /* pass pointer starting
    at the n/2th element */
    fft(x+n, n, inverse);

    complex<double> even(0,0);
    complex<double> odd(0,0);

    for (int k = 0; k < n; k++)
    {
        /* code */
        even = x[k];
        odd = x[k+n]; /* k + N/2 */

        x[k] = even + w*odd;
        x[k+n] = even - w*odd;

        /* multiply same as k+1
           in exponent */
        w = w*wn;
    }
}

void ifft(complex<double>* x, double N)
{
    fft(x, N, true);
    for(int i=0;i<N;i++) x[i] = x[i]/N;
}

#endif /* DIT_FFT_h */

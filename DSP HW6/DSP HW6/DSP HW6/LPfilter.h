//
//  LPfilter.h
//  DSP HW6
//
//  Created by 何冠勳 on 2021/1/30.
//
#ifndef LPfilter_h
#define LPfilter_h

#include <string>
#include <iostream>
#include <vector>
#define pi 3.14159265358979323846
#define order 1023
using namespace std;

double sinc(double x)
{
    if(x == 0) return 1;
    else return sin(x)/x;
}

double impulse_response(int n, int M, double w_c, string window_opt)
{
    double w, centr = (double)(n - M);
    if(window_opt == "hanning") w = 0.5 + 0.5*cos(pi*centr/M);
    else if(window_opt == "hamming") w = 0.54 + 0.46*cos(pi*centr/M);
    else w = 1;
    
    return sinc(centr*w_c)*(w_c/pi)*w;
}

vector<double> lpfilter(int gain, double cutoff, int m = order)
{
    vector<double> h_n;
    for(int j=0;j<order+1;j++)
        h_n.push_back((double)gain*impulse_response(j, (order-1)/2, cutoff, "hamming"));
    return h_n;
}

#endif /* LPfilter_h */

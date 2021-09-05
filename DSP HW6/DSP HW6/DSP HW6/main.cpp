//
//  main.cpp
//  DSP HW6
//
//  Created by 何冠勳 on 2021/1/22.
//
#include <string>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <time.h>
#include "Wav.h"
#include "DIT-FFT.h"
#include "LPfilter.h"
#define pi 3.14159265358979323846
#define p 960
#define length 2048
using namespace std;

vector<short> resample(vector<short> b, int times, string mode)
{
    vector<short> a;
    if(mode == "up")
    {
        for(unsigned i=0;i<b.size();i++)
        {
            a.push_back(b[i]);
            for(int j=0; j<times-1; j++)
                a.push_back(0);
        }
    }
    else if (mode == "down")
    {
        for(unsigned i=0;i<b.size();i+=times)
            a.push_back(b[i]);
    }
    return a;
}

double magn(complex<double> z)
{
    return sqrt(pow(z.real(),2)+pow(z.imag(),2));
}

double approx_zero(double d)
{
    if (abs(d)<0.0000000000001)
        return 0;
    else
        return d;
}

void printlist(complex<double>* l, int N, bool extras)
{
    for (int i=0;i<N;i++)
        cout << " (" << approx_zero(l[i].real()) << "," << approx_zero(l[i].imag()) << ")";

    double mag = 0;
    double phase = 0;

    if(extras)
    {
        for(int i=0;i<N;i++)
        {
            mag = magn(l[i]);
            phase = atan(l[i].imag()/l[i].real());

            cout << endl;
            cout << "Frequency bin ["<< i <<"]"<<endl;
            cout << "Magnitude: " << mag << endl;
            cout << "Phase: " << phase << endl;
        }
    }
    cout << endl;
}

void initial(complex<double>* l, int N)
{
    for (int i=0;i<N;i++)
    {
        l[i].real(0.0);
        l[i].imag(0.0);
    }
}

/* ------ Main --------- */
int main(int argc, const char* argv[])
{
    string infname = "input.wav";
    string outfname = "output.wav";
    Stereo_Wav wavein = WaveRead(infname);
    Stereo_Wav waveout;
    vector<short> x_n;
    vector<short> y_n;
    vector<double> h_n = lpfilter(80, pi/441);

    complex<double>* x = new complex<double>[length];     //input
    complex<double>* h = new complex<double>[length];     //filter
    initial(x, length);
    initial(h, length);
    
    for(int i=0;i<h_n.size();i++) h[i].real(h_n[i]);
    //cout << "Low-pass filter (time-domain):" << endl;
    //printlist(h, length, false);
    //cout << "Low-pass filter (freq-domain):" << endl;
    fft(h, length);
    //printlist(h, length, true);
    
    time_t start = time(NULL);
    /*-----------------------------------*/
    /*------------Left Data--------------*/
    /*-----------------------------------*/
    cout << "---------------Left Data---------------" << endl;
    cout << "Upsampling 80" << endl;
    x_n = resample(wavein.left_data, 80, "up");
    int seg = (x_n.size()/p)+1;
    int size_of_y = seg*p+order-1;
    complex<double>* y = new complex<double>[size_of_y];     //output
    initial(y, size_of_y);
    for(int i=0;i<seg;i++)
    {
        initial(x, length);
        if(i%200==0) cout << ".";
        for(int j=0;j<p;j++)
            x[j].real((double)x_n[j+i*p]);
        
        //cout << "Input array (time-domain):" << endl;
        //printlist(x, length, false);
        //cout << "Input array (freq-domain):" << endl;
        fft(x, length);
        //printlist(x, length, false);
        
        complex<double>* temp = new complex<double>[length];
        for(int j=0;j<length;j++) temp[j] = h[j]*x[j];
        
        ifft(temp, length);
        
        for(int j=0;j<length;j++)
            y[j+i*p] = y[j+i*p] + temp[j];

        delete[] temp;
    }
    cout << endl;
    
    for(int i=0;i<size_of_y;i++)
    {
        y_n.push_back((short)y[i].real());
    }
    cout << "Downsampling 441" << endl;
    y_n = resample(y_n, 441, "down");
    waveout.left_data = y_n;
    
    x_n.clear(); y_n.clear();
    
    /*-----------------------------------*/
    /*------------Right Data-------------*/
    /*-----------------------------------*/
    cout << "---------------Right Data---------------" << endl;
    cout << "Upsampling 80" << endl;
    x_n = resample(wavein.right_data, 80, "up");
    initial(y, size_of_y);
    for(int i=0;i<seg;i++)
    {
        initial(x, length);
        if(i%200==0) cout << ".";
        for(int j=0;j<p;j++)
            x[j].real((double)x_n[j+i*p]);
        
        //cout << "Input array (time-domain):" << endl;
        //printlist(x, length, false);
        //cout << "Input array (freq-domain):" << endl;
        fft(x, length);
        //printlist(x, length, false);
        
        complex<double>* temp = new complex<double>[length];
        for(int j=0;j<length;j++) temp[j] = h[j]*x[j];
        
        ifft(temp, length);
        
        for(int j=0;j<length;j++)
            y[j+i*p] = y[j+i*p] + temp[j];

        delete[] temp;
    }
    cout << endl;
    
    for(int i=0;i<size_of_y;i++)
    {
        y_n.push_back((short)y[i].real());
    }
    cout << "Downsampling 441" << endl;
    y_n = resample(y_n, 441, "down");
    waveout.right_data = y_n;
    
    delete[] y; delete[] x; delete[] h;
    x_n.clear(); y_n.clear(); h_n.clear();
    time_t end = time(NULL);
    
    double diff = difftime(end, start);
    cout << endl << "Conversion Done, ";
    cout << "took " << diff << " sec." << endl;
    
    /*-----------------------------------*/
    /*-------Header & Writefile----------*/
    /*-----------------------------------*/
    waveout.header = make_WaveHeader(2, 8000, 16, waveout.left_data.size());
    WaveWrite(outfname, waveout);
    
    return 0;
}

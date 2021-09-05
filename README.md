# DSP HW6
<center> <font size=2> Jason [2021/01] </font> </center>

original hackmd : https://hackmd.io/@jasonho610/SkusPSBku

## Objective
Similar to HW5, we are going to convert a .wav from sampling rate 44.1k to 8k. Nevertheless, we use FFT method along with Overlap-add method instead, in order to raise the speed of filtering.

> A fast Fourier transform (FFT) is an algorithm that computes the discrete Fourier transform (DFT) of a sequence, or its inverse (IDFT). Fourier analysis converts a signal from its original domain (often time or space) to a representation in the frequency domain and vice versa. The DFT is obtained by decomposing a sequence of values into components of different frequencies.

As for overlap-add method, the algorithm is shown as following description:
![](https://i.imgur.com/FZcyglb.png)
A sequence of 5 plots depicts one cycle of the Overlap-add convolution algorithm. The first plot is a long sequence of data to be processed with a lowpass FIR filter. The 2nd plot is one segment of the data to be processed in piecewise fashion. The 3rd plot is the filtered segment, including the filter rise and fall transients. The 4th plot indicates where the new data will be added with the result of previous segments. The 5th plot is the updated output stream. The FIR filter is a boxcar lowpass with M=16 samples, the length of the segments is L=100 samples and the overlap is L+M-1=15 samples.

## Methodology
Therefore, we follow the algorithm of Overlap-add, and substitute convolution for FFT so as to achieve upgrade of speed.

For segmenting input signal $x[n]$:
![](https://i.imgur.com/Cippdz0.png)

For low-pass filter $h[n]$:
![](https://i.imgur.com/5IFGUNZ.png)

The definition of parameters *P, Q, N*:
![](https://i.imgur.com/bV8b5OU.png)


## Code Demonstration
The final code includes three parts:
1. main.cpp
2. Wav.h
3. DIT-FFT.h
4. LPfilter.h

### Stereo Wav Structure
Defined in Wav.h
```cpp
typedef struct WaveHeader
{
    // riff wave header
    char ChunkID[4] = {'R','I','F','F'};
    unsigned ChunkSize;        // 0 ~ FFFF,FFFF
    char Format[4] = {'W','A','V','E'};
    
    // fmt subchunk
    char SubChunk1ID[4] = {'f','m','t',' '};
    unsigned SubChunk1Size;    // 0 ~ FFFF,FFFF
    unsigned short AudioFormat;    // 0 ~ FFFF
    unsigned short NumChannels;    // 0 ~ FFFF
    unsigned SampleRate;       // 0 ~ FFFF,FFFF
    unsigned ByteRate;         // 0 ~ FFFF,FFFF
    unsigned short BlockAlign;     // 0 ~ FFFF
    unsigned short BitsPerSample;  // 0 ~ FFFF
    
    // data subchunk
    char SubChunk2ID[4] = {'d','a','t','a'};
    unsigned SubChunk2Size;    // 0 ~ FFFF,FFFF
} WaveHeader;

typedef struct Stereo_Wav
{
    WaveHeader header;
    vector<short> left_data;
    vector<short> right_data;
} Stereo_Wav;

WaveHeader make_WaveHeader(unsigned short const NumChannels, unsigned const SampleRate, unsigned short const BitsPerSample, long NoS)
{
    WaveHeader myWH;

    myWH.AudioFormat = 1;                  // 1 for PCM...
    myWH.SampleRate = SampleRate;
    myWH.NumChannels = NumChannels;        // 1 for Mono, 2 for Stereo
    myWH.BitsPerSample = BitsPerSample;
    myWH.ByteRate = (myWH.SampleRate * myWH.NumChannels * myWH.BitsPerSample)/8;
    myWH.BlockAlign = (myWH.NumChannels * myWH.BitsPerSample)/8;
    myWH.SubChunk2Size = (NoS * myWH.NumChannels * myWH.BitsPerSample)/8;
    myWH.SubChunk1Size = 16;               // 16 for PCM
    myWH.ChunkSize = 4+(8+myWH.SubChunk1Size)+(8+myWH.SubChunk2Size);

    return myWH;
}
```
### Wavfile Operation
Defined in Wav.h
```cpp
Stereo_Wav WaveRead(string filename)
{
    ifstream infile;
    WaveHeader header;
    vector<short> data;
    Stereo_Wav wavein;
    
    infile.open(filename, ofstream::binary|ios::in);
    if (!infile.is_open())
    {
        cerr << "Could not open the file" << endl;
        system("pause");
    }
    
    infile.read((char *) &header, sizeof(header));
    //cout << header.SampleRate << endl;
    wavein.header = header;

    while(!infile.eof())
    {
        short temp;        // data can't be greater than FFFF(65535).
        infile.read((char *) &temp, sizeof(temp));
        data.push_back(temp);
    }
    infile.close();

    /*------------------------------------*/
    /* Change data length here for testing*/
    for(unsigned i=0;i<data.size();i=i+2)
    {
        wavein.left_data.push_back(data[i]);
        wavein.right_data.push_back(data[i+1]);
    }
    
    return wavein;
}

void WaveWrite(string filename, Stereo_Wav waveout)
{
    ofstream outfile;
    outfile.open(filename, ofstream::binary|ofstream::out);
    if (!outfile.is_open())
    {
        cerr << "Could not open the file" << endl;
        system("pause");
    }
    
    outfile.write((char*)&waveout.header, sizeof(waveout.header));
    
    for(unsigned i=0;i<waveout.left_data.size();i++)
    {
        outfile.write((char*)&waveout.left_data[i], sizeof(waveout.left_data[i]));
        outfile.write((char*)&waveout.right_data[i], sizeof(waveout.right_data[i]));
    }
}
```
### Resample
Defined in main.cpp
```cpp
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
```
### Low-pass filter
Defined in LPfilter.h
```cpp
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
        h_n.push_back((double)gain*impulse_response(j, (order-1)/2, cutoff, "hanning"));
    return h_n;
}
```
### DIT-FFT
DIT-FFT stands for "Decimate in Time - Fast Fourier Transform".
```cpp
const complex<double> J(0,1);
void split(complex<double> *inputs, double N);
void fft(complex<double>* x, double N, bool inverse);
void ifft(complex<double>* x, double N);
```

```cpp
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
```
### Main Procedure
```cpp
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
```
## Results
### Output Wav Info
![](https://i.imgur.com/fVJwrLr.png)

![](https://i.imgur.com/z0fCzAv.png)

|  | Hex | Dec | 
| --- | :---: | --- |
| ChunkSize | 00 37 62 B8 | <font color = #F28C06>3629752 bytes</font> |
| Subchunk1Size | 00 00 00 10 | 16 bytes |
| AudioFormat | 00 01 | 1 -> PCM (i.e. Linear quantization) |
| NumChannels | 00 02 | 2 -> Stereo |
| SampleRate | 00 00 1F 40 | <font color = #03CF1F>8000 Hz</font> |
| ByteRate | 00 00 7D 00 | 32000 |
| BlockAlign | 00 04 | 4 |
| BitsPerSample | 00 10 | 16 bits |
| Subchunk2Size | 00 37 62 94 | <font color = #9803CF>3629716 bytes</font> |

Total file size: <font color = "red">3629760 bytes</font>
= <font color = #F28C06>3629752 bytes</font> + 8 (bytes of "RIFF") 
= <font color = #9803CF>3629716 bytes</font> + 44 (bytes of Wavheader)

It took 4131(s) = 1.14(hr) to complete the task.
This consequence is barely satisfactory; there are some improvement we can do, for example: Replace our recursive DIT-FFT.

According to Master Thereom, 
![](https://i.imgur.com/QF5KC9r.png)
the time complexity of our FFT is 
$$\begin{eqnarray}T(n)&=&2T(\frac{n}{2})+1 \\ &=&\theta(nlogn)\end{eqnarray}$$


//
//  Wav.h
//  DSP HW6
//
//  Created by 何冠勳 on 2021/1/23.
//

#ifndef Wav_h
#define Wav_h
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

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

WaveHeader make_WaveHeader(unsigned short const NumChannels, unsigned const SampleRate, unsigned short const BitsPerSample, unsigned NoS)
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

#endif /* Wav_h */

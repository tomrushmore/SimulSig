//
//  SetupMacros.h
//  FFT Chain
//
//  Created by Thomas Rushmore on 26/03/2013.
//  Copyright (c) 2013 Thomas Rushmore. All rights reserved.
//

#ifndef FFT_Chain_SetupMacros_h
#define FFT_Chain_SetupMacros_h

#include <limits>

static int calcMs(int tempo);

#define NSAMP 22050
#define FNMELBANDS 14
#define NMELBANDS FNMELBANDS * 2
#define MFCCGROUPING 6
#define MAXCLUSTITERATIONS 5000
#define NUMCLUSTERS 8000

#define MINFREQ 0
#define MAXFREQ 22050
#define FFTSIZE 4096
#define SR 44100
#define HOPMS 20
#define NUMMKVGEN 2
#define HOPOVRLAPMS 4096
#define TEMPO 120
//"/Users/TomRushmore/Final Project/Products/"
#define FOLDERPATH "/Users/thomasrushmore/Desktop/Products/"
#define SAMPMACHINEPATH FOLDERPATH
#define MACHINEPATH FOLDERPATH

static int calcMs(int tempo)
{
    int bps = tempo/60;
    int hop = SR/bps - FFTSIZE;
    bool r = true;
    int i = 2;
    while(r){
        hop % i == 0 && hop / i < FFTSIZE ? r = false : i++;
    }
    return hop / i;
}

#endif

SimulSig
========

A new approach to audio-mosaicing using probabilistic modelling

NOTE: The original repo containing commits from the last 4 months can be found here : https://github.com/tomrushmore/historic-simulsig-repo
Numerous clones were unremovable during git cleanup resulting in an unweidly size of 1.5gb


README
========

Once downloaded, open Simulsig.xcodeproj found in Markov Chains/FFT Chain/
1. In SetupMacros.h, found in the Main / Setup.h subfolder in Navigator window, change
the Macro FOLDERPATH to the included Products directory.
2. Populate the Analysis Corpus folder, found in Products/, with a selection of mono .wav files. 
3. Populate the Training Corpus folder, found in Products/, with a selection of mono .wav files on which the Markov chain shall be trained.
4. In main.cpp, change the variable Phase to AllPhase to run through the entire program.
dependencies. If you do not wish to see the multitude of console printouts, alter the macro in PrintfDebug.h

Several macros determine the type out output produced, all of which are found in SetupMacros.h 
#NSAMP - Audio segment size
#FNMELBANDS - Number of mel-filterbanks used during extraction. this increases the number of MFFCs
#MFCCGROUPING - Number of MFC frames grouped together
#MAXCLUSTITERATIONS - Maximum number of iterations of K-means algorithm. The larger the better.
#NUMCLUSTER - Maximum number of clusters. The majority of clusters do not have a unique centroid, so the number of active clusters will be significantly less than this.


Texture-synthesis like results:
NSAMP - 4096
FNMELBANDS - 14
MFCCGROUPING - 1

Automatic track remixer example:
NSAMP - 22050
FNMELBANDS - 14
MFCCGROUPING - 6

========
FFTW3,
libsndfile,
TNT,
LAPACK,
NanoFlann (C++ only subset of FLANN),
Accelerate framework (subset of within OSX veclib),
Dirent.h (included by default OSX 10.7 >  )
/*
A simple sound library for CSE 20211 by Douglas Thain
For course assignments, you should not change this file.
For complete documentation, see:
http://www.nd.edu/~dthain/courses/cse20211/fall2013/wavfile
*/

#ifndef WAVFILE_H
#define WAVFILE_H

#include <stdio.h>
//#include <inttypes.h>

FILE * wavfile_open( const char *filename );
void wavfile_write( FILE *file, short data[], int length );
void wavfile_close( FILE * file );

#define WAVFILE_SAMPLES_PER_SECOND 44100

#endif

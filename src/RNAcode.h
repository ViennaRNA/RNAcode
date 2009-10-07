/*  Copyright 2009, Stefan Washietl

    This file is part of RNAcode.

    RNAcode is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RNAcode is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RNAcode.  If not, see <http://www.gnu.org/licenses/>. */


#ifndef _RNACODE_H_
#define _RNACODE_H_

void usage(void);
void help(void);
void version(void);
void read_commandline(int argc, char *argv[]);

typedef struct _parameters parameters;

struct _parameters{
  float omega;
  float Omega;
  float Delta;
  float stopPenalty_0;
  float stopPenalty_k;
  FILE *inputFile;
  FILE *outputFile;
  FILE *debugFile;
  int sampleN;
  char limit[10000];
  float cutoff;
  int outputFormat;
  int stopEarly;
  int bestOnly;
  int bestRegion;
  int blosum;
  char debugFileName[1024];
  char inputFileName[1024];
  float printIfBelow; 
  float printIfAbove;
  int postscript;
  float postscript_cutoff;
  char postscriptDir[1024];

};

#endif

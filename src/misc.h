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

#ifndef _MISC_H_
#define _MISC_H_

#include "score.h"
#include "code.h"
#include "rnaz_utils.h"

void reintroduceGaps(const struct aln* origAln[], struct aln* sampledAln[]);

void sortAln(const struct aln* origAln[], struct aln* sampledAln[]);

void freeResults(segmentStats results[]);

//void getBlock(const char* seq_0, const char* seq_k, int i, char* block_0, char* block_k, int* z );
void getBlock(int i, const char* seq_0, const char* seq_k, const int* map_0, const int* map_k, char* block_0, char* block_k, int* z );

int pos2col(const char* seq, int pos);

int getSeqLength(char *seq);

int hDist(int a1, int a2, int a3, int b1, int b2, int b3);

float avg(float* data, int N);

float stddev(float* data, int N);

float gaussian (const float sigma);

void copyAln(struct aln *src[],struct aln *dest[]);

void printAlnClustal(FILE *out, const struct aln* AS[]);

void printResults(FILE* outfile, int outputFormat, segmentStats results[]);

#endif

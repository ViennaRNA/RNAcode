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

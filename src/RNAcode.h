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
  char debugFileName[1024];
  char inputFileName[1024];
  float printIfBelow; 
  float printIfAbove;
};

#endif

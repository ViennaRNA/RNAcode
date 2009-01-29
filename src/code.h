#ifndef _CODE_H_
#define _CODE_H_


#include "rnaz_utils.h"
#include "tree.h"
#include "treeSimulate.h"

typedef struct _bgModel bgModel;

struct _bgModel {

  double scores[4];
  double probs[4];
  double kappa;
  double dist;
  double freqs[4];
  int** matrix;

};

typedef struct _segmentStats segmentStats;

struct _segmentStats {

  int startSite;
  int endSite;
  int strand;
  int frame;
  int start;
  int end;
  char *name;
  double score;
  double pvalue;
  
};

/* Values set in main() */
int ntMap[256]; 

extern int transcode[4][4][4];

double probHKY(int i, int j, double d, double freqs[4], double kappa);

int** getScoringMatrix();

void freeScoringMatrix(int** matrix);

void calculateBG(bgModel* models);

int* encodeSeq(char* seq);

int encodeAA(char aa);

char decodeAA(int encodedAA);

char* translateSeq(char* seq);

int hDist(int a1, int a2, int a3, int b1, int b2, int b3);

double avg(double* data, int N);

double stddev(double* data, int N);

double gaussian (const double sigma);

int compareScores(const void * a, const void * b);

void printAlnClustal(FILE *out, const struct aln* AS[]);

void stripGaps(struct aln* AS[]);

void printResults(FILE* outfile, int outputFormat, segmentStats results[]);

void copyAln(struct aln *src[],struct aln *dest[]);

void countFreqsMono(const struct aln *alignment[], double freqs[]);

double* sumOfPairScore(bgModel* models, const struct aln *alignment[],int from, int to);

double* getCumSum(double* scores, int N);

bgModel* getModelMatrix(TTree* tree, struct aln *alignment[], double kappa);

void freeModelMatrix(bgModel* models, int N);

void getExtremeValuePars(TTree* tree, bgModel* models, const struct aln *alignment[],                          
                         int sampleN, int sampleMode, double* parMu, double* parLambda);

segmentStats* getHSS(bgModel* models, const struct aln** inputAln, double parMu, double parLambda, double cutoff);

void reintroduceGaps(const struct aln* origAln[], struct aln* sampledAln[]);

void randomlyPlaceGaps(const struct aln* origAln[], struct aln* sampledAln[]);

void sortAln(const struct aln* origAln[], struct aln* sampledAln[]);

void freeResults(segmentStats results[]);


#endif

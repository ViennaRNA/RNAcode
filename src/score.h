#ifndef _SCORE_H_
#define _SCORE_H_

#include "rnaz_utils.h"
#include "tree.h"
#include "treeSimulate.h"

#define UNDEF -9.0
#define MINUS_INF -99.0

#define MAX(x,y)       (((x)>(y)) ? (x) : (y))
#define MAX3(x,y,z)    (MAX(  (MAX((x),(y))) ,(z)))

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

int compareScores(const void * a, const void * b);

void stripGaps(struct aln* AS[]);

void countFreqsMono(const struct aln *alignment[], double freqs[]);

double* sumOfPairScore(bgModel* models, const struct aln *alignment[],int from, int to);

double* getCumSum(double* scores, int N);

bgModel* getModelMatrix(TTree* tree, struct aln *alignment[], double kappa);

void freeModelMatrix(bgModel* models, int N);

void getExtremeValuePars(TTree* tree, const struct aln *alignment[], 
                         int sampleN, double* parMu, double* parLambda);

segmentStats* getHSS(double** S, const struct aln** inputAln,  char strand, double parMu, double parLambda, double cutoff);


double**** getPairwiseScoreMatrix(bgModel* models, const struct aln *alignment[]);
double** getMultipleScoreMatrix(double**** Sk, bgModel* models, const struct aln *alignment[]);

double* backtrack(double**** S, int k, int from, int to, const struct aln *alignment[]);

void freeSk (double**** S, const struct aln *alignment[]);
void freeS (double** S, const struct aln *alignment[]);

segmentStats* scoreAln(const struct aln *alignment[], TTree* tree, double kappa, double parMu, double parLambda);


#endif

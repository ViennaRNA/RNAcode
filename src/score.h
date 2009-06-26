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

  float scores[4];
  float probs[4];
  float kappa;
  float dist;
  float freqs[4];
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
  float score;
  float pvalue;
  int hide;
  
};

/* Values set in main() */
int ntMap[256]; 

extern int transcode[4][4][4];

float probHKY(int i, int j, float d, float freqs[4], float kappa);

int** getScoringMatrix();

void freeScoringMatrix(int** matrix);

void calculateBG(bgModel* models);

void stripGaps(struct aln* AS[]);

void countFreqsMono(const struct aln *alignment[], float freqs[]);

float* sumOfPairScore(bgModel* models, const struct aln *alignment[],int from, int to);

float* getCumSum(float* scores, int N);

bgModel* getModelMatrix(TTree* tree, struct aln *alignment[], float kappa);

void freeModelMatrix(bgModel* models, int N);

int getExtremeValuePars(TTree* tree, const struct aln *alignment[], 
                        int sampleN, float maxNativeScore, float* parMu, float* parLambda);

segmentStats* getHSS(float** S, const struct aln** inputAln,  char strand);

void getPairwiseScoreMatrix(bgModel* models, const struct aln *alignment[]);
float** getMultipleScoreMatrix(float**** Sk, bgModel* models, const struct aln *alignment[]);

float* backtrack(float**** S, int k, int from, int to, const struct aln *alignment[]);

void freeSk (float**** S, const struct aln *alignment[]);
void freeS (float** S, const struct aln *alignment[]);

segmentStats* scoreAln(const struct aln *alignment[], TTree* tree, float kappa);

bgModel* getModels(TTree* tree, struct aln *alignment[], float kappa);


#endif

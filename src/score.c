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


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "score.h"
#include "RNAcode.h"
#include "code.h"
#include "misc.h"
#include "extreme_fit.h"

extern parameters pars;

extern int transcode[4][4][4];
extern int BLOSUM62[24][24];
extern int BLOSUM90[24][24];
extern bgModel *models;
extern bgModel *modelsRev;

extern float**** Sk;
extern float**** Sk_native;   
extern float**** Sk_native_rev;  

/*********************************************************************
  getScoringMatrix

  Returns matrix as int**, currently hard-coded for BLOSUM62

*********************************************************************/ 

int** getScoringMatrix(){

  int** matrix;
  int i,j;

  matrix=(int**)malloc(sizeof(int*)*24);

  for (i=0;i<24;i++){
    matrix[i]=(int*)malloc(sizeof(int)*24);
  }
  
  for (i=0;i<24;i++){
    for (j=0;j<24;j++){

      if (pars.blosum == 62){
        matrix[i][j]=BLOSUM62[i][j];
      }

      if (pars.blosum == 90){
        matrix[i][j]=BLOSUM90[i][j];
      }
    }
  }

  return matrix;

}

/*********************************************************************
  freeScoringMatrix

  Frees scoring matrix

*********************************************************************/ 

void freeScoringMatrix(int** matrix){

  int i;

  for (i=0;i<24;i++){
    free(matrix[i]);
  }
  free(matrix);
}

/*********************************************************************
  calculateBG

  model ... pointer to bgModel structure with all necessary information

  Fills the fields "scores[]" and "probs[]" in model based on kappa,
  dist, freqs and matrix. "probs[h]" holds the probabilities for pairs
  with hamming distance h=0,1,2,3 while scores[h] holds the expected
  score for pairs with hamming distance h.

*********************************************************************/ 

void calculateBG(bgModel* model){

  int a1,a2,a3, b1,b2,b3;
  int pepA, pepB;
  float f, prob, score, probStop;
  float counts[4];
  float probs[4][4];
  int h, i,j;

  for (i=0;i<4;i++){
    for (j=0;j<4;j++){
      probs[i][j]=probHKY(i,j,model->dist,model->freqs,model->kappa);
    }
  }
 
  counts[0]=counts[1]=counts[2]=counts[3]=0.0;
  model->scores[0]=model->scores[1]=model->scores[2]=model->scores[3]=0.0;

  // First calculate probability of pairs with stop
  probStop=0;
  for (a1=0;a1<4;a1++){
    for (a2=0;a2<4;a2++){
      for (a3=0;a3<4;a3++){
        for (b1=0;b1<4;b1++){
          for (b2=0;b2<4;b2++){
            for (b3=0;b3<4;b3++){
              pepA=transcode[a1][a2][a3];
              pepB=transcode[b1][b2][b3];
              if (pepA!=-1 && pepB!=-1) continue;

              f=(model->freqs[a1])*(model->freqs[a2])*(model->freqs[a3]);

              prob=probs[a1][b1]*probs[a2][b2]*probs[a3][b3];

              prob*=f;
              probStop+=prob;
            }
          }
        }
      }
    }
  }

  // Then calculate actual probabilities and expected scores for h=0,1,2,3
  for (a1=0;a1<4;a1++){
    for (a2=0;a2<4;a2++){
      for (a3=0;a3<4;a3++){

        pepA=transcode[a1][a2][a3];
        if (pepA==-1) continue;

        for (b1=0;b1<4;b1++){
          for (b2=0;b2<4;b2++){
            for (b3=0;b3<4;b3++){

              pepB=transcode[b1][b2][b3];

              if (pepB==-1) continue;

              h=hDist(a1,a2,a3,b1,b2,b3);

              f=(model->freqs[a1])*(model->freqs[a2])*(model->freqs[a3]);
              
              prob=probs[a1][b1]*probs[a2][b2]*probs[a3][b3];

              prob*=f;
              prob/=(1-probStop); //Correct probability for the condition that no stop exists
              
              if (pepA==-1) pepA=23;
              if (pepB==-1) pepB=23;
              
              score=model->matrix[pepA][pepB];
              counts[h]+=prob;
              model->scores[h]+=score*prob;

            }
          }
        }
      }
    }
  }

  for (i=0;i<4;i++){
    model->scores[i]/=counts[i];
    model->probs[i]=counts[i];
  }
}

/*********************************************************************
  probHKY

  Calculates probability of change i -> j given distance d, stationary
  frequencies freqs[] and ts/tv rate ration kappa under the HKY
  substitution model

*********************************************************************/ 

float probHKY(int i, int j, float d, float freqs[4], float kappa){

  float piA, piC, piG, piT;
  float piR, piY, r, l, k1, k2, exp1, exp2, exp22, exp21;
  float result[4][4];

  piA=freqs[0];
  piC=freqs[1];
  piG=freqs[2];
  piT=freqs[3];
  
  piR=piA+piG;
  piY=piT+piC;
  r=1./(2.*(piA*piC+piC*piG+piA*piT+piG*piT+kappa*(piC*piT+piA*piG)));
  l=r*d;
  k1=kappa*piY+piR;
  k2=kappa*piR+piY;
  exp1=exp(-l);
  exp22=exp(-k2*l);
  exp21=exp(-k1*l);

  result[0][0]=piA*(1.+(piY/piR)*exp1)+(piG/piR)*exp22;
  result[0][1]=piC*(1.-            exp1);
  result[0][2]=piG*(1.+(piY/piR)*exp1)-(piG/piR)*exp22;
  result[0][3]=piT*(1.-            exp1);
  result[1][0]=piA*(1.-            exp1);
  result[1][1]=piC*(1.+(piR/piY)*exp1)+(piT/piY)*exp21;
  result[1][2]=piG*(1.-            exp1);
  result[1][3]=piT*(1.+(piR/piY)*exp1)-(piT/piY)*exp21;
  result[2][0]=piA*(1.+(piY/piR)*exp1)-(piA/piR)*exp22;
  result[2][1]=piC*(1.-            exp1);
  result[2][2]=piG*(1.+(piY/piR)*exp1)+(piA/piR)*exp22;
  result[2][3]=piT*(1.-            exp1);
  result[3][0]=piA*(1.-            exp1);
  result[3][1]=piC*(1.+(piR/piY)*exp1)-(piC/piY)*exp21;
  result[3][2]=piG*(1.-            exp1);
  result[3][3]=piT*(1.+(piR/piY)*exp1)+(piC/piY)*exp21;

  return(result[i][j]);

}


/*********************************************************************
  countFreqsMono

  Counts mononucleotide frequencies in alignment aln and writes it to
  freqs.

*********************************************************************/ 

void countFreqsMono(const struct aln *alignment[], float freqs[]){

  int i,k;
  char* currSeq;
  char c;
  unsigned long counter;
  float sum;

  for (i=0;i<4;i++) freqs[i]=0.0;

  counter=0;

  for (i=0; alignment[i]!=NULL; i++){

    currSeq=alignment[i]->seq;

    k=0;
    while ((c=currSeq[k++])!='\0'){
      if (c=='-') continue;
      freqs[ntMap[c]]++;
      counter++;
    }
  }

  for (i=0;i<4;i++) freqs[i]/=(float)counter;
}


/*********************************************************************
  getModels

  Calculates background models for a tree, alignment and ts/tv rate
  ration kappa. Returns results as pointer to bgModel structure.

*********************************************************************/ 

bgModel* getModels(TTree* tree, struct aln *alignment[], float kappa){

  int i,j,N;
  int **scoringMatrix;
  float** distanceMatrix;
  float freqsMono[4];
  bgModel* models;
  float* w;

  for (N=0; alignment[N]!=NULL; N++);   
  
  distanceMatrix=getDistanceMatrix(tree,(struct aln**)alignment);

  //w = weights(distanceMatrix,N);

  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      //printf("%.2f ",w[i]);
      //printf("%.2f ",distanceMatrix[i][j]);
    }
    //printf("\n");
  }
  //printf("\n");

  for (j=0;j<N;j++){
    //printf("%.2f ",w[j]);
  }
  //printf("\n");


  countFreqsMono((const struct aln**)alignment, (float *) freqsMono);

  models=(bgModel*)malloc(sizeof(bgModel)*N);

  for (j=0;j<N;j++){
    scoringMatrix=getScoringMatrix();
    models[j].dist=distanceMatrix[0][j];
    models[j].kappa=kappa;
    models[j].freqs[0]=freqsMono[0];
    models[j].freqs[1]=freqsMono[1];
    models[j].freqs[2]=freqsMono[2];
    models[j].freqs[3]=freqsMono[3];
    models[j].matrix=scoringMatrix;
    calculateBG(&models[j]);
  }

  for (i=0;i<N;i++){
    free(distanceMatrix[i]);
  }
  free(distanceMatrix);

  return models;
  
}

/*********************************************************************
  freeModels

  Frees background model, N is the number of sequences in the
  for which the model was calculated

*********************************************************************/ 

void freeModels(bgModel* models, int N){

  int i, j;

  for (j=0;j<N;j++){
    freeScoringMatrix(models[j].matrix);
  }

  free(models);
}

/*********************************************************************
  calculateSigma

  Calculates sigma scores for two alignment blocks between the 
  reference sequence and sequence k. The required background model data
  is given in models. 

*********************************************************************/ 


float calculateSigma(char* block_0, char* block_k, int k, bgModel* models){

  char codonA[4]="XXX";
  char codonB[4]="XXX";
  int pepA, pepB;
  int i,j,h;
  float currScore,expectedScore, observedScore;

  i=j=0;
  while (block_0[i] != '\0') {
    if (block_0[i] != '-') {
      codonA[j]=block_0[i];
      codonB[j]=block_k[i];
      j++;
    }
    i++;
  }

  // Second sequence was gap only 
  if (strcmp(codonB,"XXX")==0){
    return 0.0;
  }

  codonA[3]=codonB[3]='\0';

  if (codonA[0] == 'N' || codonA[1] == 'N' || codonA[2] == 'N' ||
      codonB[0] == 'N' || codonB[1] == 'N' || codonB[2] == 'N'
      ) {
    return 0.0;
  }

  h=hDist(ntMap[codonA[0]], ntMap[codonA[1]], ntMap[codonA[2]],
          ntMap[codonB[0]], ntMap[codonB[1]], ntMap[codonB[2]]);
  
  if (h==0) return 0.0;

  pepA=transcode[ntMap[codonA[0]]][ntMap[codonA[1]]][ntMap[codonA[2]]];
  pepB=transcode[ntMap[codonB[0]]][ntMap[codonB[1]]][ntMap[codonB[2]]];

  if (pepA==-1) {
    return pars.stopPenalty_0;
  }

  if (pepB==-1) {
    return pars.stopPenalty_k;
  }

  expectedScore=models[k].scores[h];
  observedScore=(float)models[k].matrix[pepA][pepB];

  return observedScore-expectedScore;
}


/********************************************************************
  getPairwiseScoreMatrix

  Calculates the main score matrix for all pairs of the reference
  sequence with all other sequences. 

  k   ... number of sequence
  x   ... State (0 ... '0', 1 ... '+1', 2 ... '-1')
  b,i ... subsequence from b to i

*********************************************************************/ 

void getPairwiseScoreMatrix(bgModel* models, const struct aln *alignment[]){
  
  int L,N, colsN, pos,z,b,i,x, k,l;
  char *block_0, *block_k, *seq_0, *seq_k;
  char **blocks_0, **blocks_k;
  int** zs;
  int *map_0, *map_k;
  float** sigmas;

  seq_0=alignment[0]->seq;
  L=getSeqLength(seq_0);

  colsN=strlen(seq_0);
  for (N=0; alignment[N]!=NULL; N++);

  // Do this only once
  if (Sk == NULL){
    Sk = allocateSk(N, L);
  }
  
  // Allocate conservatively for full length of sequence + gaps;
  block_0 = (char*) malloc(sizeof(char)*(colsN+1));
  block_k = (char*) malloc(sizeof(char)*(colsN+1));
  
  // Precalculate blocks, z's and sigmas for all i to avoid doing this in the inner loop
  blocks_0 = (char**) malloc(sizeof(char*)*(colsN+1));
  blocks_k = (char**) malloc(sizeof(char*)*(colsN+1));

  zs = (int**) malloc(sizeof(int*)*N);
  sigmas = (float**) malloc(sizeof(float*)*N);

  map_0=(int*)malloc(sizeof(int)*(colsN+1));
  map_k=(int*)malloc(sizeof(int)*(colsN+1));

  for (l=1;l<=L;l++){
    map_0[l]=pos2col(seq_0,l);
  }

  for (k=1;k<N;k++){

    seq_k=alignment[k]->seq;
    
    zs[k]=(int*) malloc(sizeof(int)*(L+1));
    sigmas[k]=(float*) malloc(sizeof(float)*(L+1));
    for (l=1;l<=L;l++){
      map_k[l]=pos2col(seq_k,l);
    }
    for (x=3;x<L+1;x+=1){
      getBlock(x, seq_0, seq_k, map_0, map_k, block_0, block_k, &z );
      blocks_0[x]=strdup(block_0);
      blocks_k[x]=strdup(block_k);
      zs[k][x]=z;
      sigmas[k][x]= calculateSigma(block_0,block_k,k,models);
    }

    for (b=1;b<L+1;b++){
      for (i=b+2;i<L+1;i+=3){
        z=zs[k][i];
        
        if (i-3 < b) {
          Sk[k][0][b][i-3]=0.0;
          Sk[k][1][b][i-3]=0.0;
          Sk[k][2][b][i-3]=0.0;
        }

        if (z==0){
          Sk[k][0][b][i]=Sk[k][0][b][i-3]+sigmas[k][i];
          Sk[k][1][b][i]=Sk[k][1][b][i-3]+pars.omega;
          Sk[k][2][b][i]=Sk[k][2][b][i-3]+pars.omega;
        }

        if (z==+1){
          Sk[k][0][b][i]=MAX(Sk[k][0][b][i-3]+pars.Delta,
                             Sk[k][2][b][i-3]+pars.Omega);
          
          Sk[k][1][b][i]=MAX(Sk[k][0][b][i-3]+pars.Omega,
                             Sk[k][1][b][i-3]+pars.Delta);
        
          Sk[k][2][b][i]=MAX(Sk[k][1][b][i-3]+pars.Omega,
                             Sk[k][2][b][i-3]+pars.Delta);
        }
        
        if (z==-1){
          Sk[k][0][b][i]=MAX(Sk[k][0][b][i-3]+pars.Delta,
                             Sk[k][1][b][i-3]+pars.Omega);
        
          Sk[k][1][b][i]=MAX(Sk[k][1][b][i-3]+pars.Delta,
                             Sk[k][2][b][i-3]+pars.Omega);
        
          Sk[k][2][b][i]=MAX(Sk[k][2][b][i-3]+pars.Delta,
                             Sk[k][0][b][i-3]+pars.Omega);
          
        }
      }
    }

    for (x=3;x<L+1;x+=1){
      free(blocks_0[x]);
      free(blocks_k[x]);
    }

    free(sigmas[k]);
    free(zs[k]);

  }

  free(sigmas);
  free(zs);
  free(blocks_0);
  free(blocks_k);
  free(block_0);
  free(block_k);
  free(map_0);
  free(map_k);

}

backtrackData* backtrack(int opt_b, int opt_i, float**** SSk , const struct aln *alignment[]){

  float opt_score;
  int opt_state; 
  int curr_state, prev_state;
  char *display_line1;
  char *display_line2;
  char *display_line3;
  char *display_line4;

  char* block_0;
  char* block_k;
  int* states;
  int *map_0, *map_k;
  
  char* seq_0;
  char* seq_k;

  int pos,x,l,z,L, i, b,k, N, colsN;
  backtrackData* output;
  int transition;
  int* transitions;

  for (N=0; alignment[N]!=NULL; N++);
  
  seq_0=alignment[0]->seq;
  
  L=getSeqLength(seq_0);
  colsN=strlen(seq_0);
  
  output = (backtrackData*) malloc(sizeof(backtrackData)*N);

  for (k = 1; k < N; ++k){

    output[k].states = (int*)  malloc(sizeof(int) * (colsN+1));
    output[k].z =      (int*)  malloc(sizeof(int) * (colsN+1));
    output[k].transitions = (int*)  malloc(sizeof(int) * (colsN+1));
    output[k].scores = (float*)  malloc(sizeof(float) * (colsN+1));
    
    seq_k=alignment[k]->seq;
  
    char string[1000];

    opt_score=MINUS_INF;
    opt_state=-1;

    //printf("INHERE: %i %i\n", opt_b, opt_i);

    for (x=0;x<3;x++){
      //printf("%.2f\n",SSk[k][x][opt_b][opt_i]);
      if (SSk[k][x][opt_b][opt_i]>opt_score){
        opt_score=SSk[k][x][opt_b][opt_i];
        opt_state=x;
      }
    }

    //return NULL;

    //opt_score=SSk[k][0][b][i];
  
    //printf("Max score: %.1f at b=%i, i=%i at state %i\n", opt_score, opt_b, opt_i, opt_state);

    states=(int*)malloc(sizeof(int)*(L+3));
      
    for (i=0;i<=L+1;i++) states[i]=-1;
  
    // Allocate conservatively for full length of sequence;
    block_0 = (char*) malloc(sizeof(char)*(colsN+1));
    block_k = (char*) malloc(sizeof(char)*(colsN+1));

    map_0=(int*)malloc(sizeof(int)*(colsN+1));
    map_k=(int*)malloc(sizeof(int)*(colsN+1));

    for (l=1;l<=L;l++){
      map_0[l]=pos2col(seq_0,l);
      map_k[l]=pos2col(seq_k,l);
    }

    curr_state=opt_state;
    b=opt_b;
  
    for (i=opt_i;i>=opt_b+2;i-=3){
      getBlock(i, seq_0, seq_k, map_0, map_k, block_0, block_k, &z );

      if (z==0){
        prev_state = curr_state;
        transition=0;
      }

      if (z==+1){
        
        if (curr_state == 0){
          if (CMP(SSk[k][0][b][i],SSk[k][0][b][i-3]+pars.Delta)){
            transition=2;
            prev_state=0;
          }
          
          if (CMP(SSk[k][0][b][i],SSk[k][2][b][i-3]+pars.Omega)){
            transition=1;
            prev_state=2;
          }
        }
          
        if (curr_state == 1){
          if (CMP(SSk[k][1][b][i],SSk[k][0][b][i-3]+pars.Omega)){
            transition=1;
            prev_state=0;
          }
          if (CMP(SSk[k][1][b][i],SSk[k][1][b][i-3]+pars.Delta)){
            transition=1;
            prev_state=1;
          }
        }

        if (curr_state == 2){
          if (CMP(SSk[k][2][b][i],SSk[k][1][b][i-3]+pars.Omega)){
            transition=1;
            prev_state=1;
          }
          if (CMP(SSk[k][2][b][i],SSk[k][2][b][i-3]+pars.Delta)){
            transition=2;
            prev_state=2;
          }
        }
      }

        
      if (z==-1){
        if (curr_state == 0){
          if (CMP(SSk[k][0][b][i],SSk[k][0][b][i-3]+pars.Delta)){
            transition=2;
            prev_state=0;
          }
          if (CMP(SSk[k][0][b][i],SSk[k][1][b][i-3]+pars.Omega)){
            transition=1;
            prev_state=1;
          }
        }
          
        if (curr_state == 1){
          if (CMP(SSk[k][1][b][i],SSk[k][1][b][i-3]+pars.Delta)){
            transition=2;
            prev_state=1;
          }
          if (CMP(SSk[k][1][b][i],SSk[k][2][b][i-3]+pars.Omega)){
            transition=1;
            prev_state=2;
          }
        }

        if (curr_state == 2){
          if (CMP(SSk[k][2][b][i],SSk[k][2][b][i-3]+pars.Delta)){
            transition=2;
            prev_state=2;
          }
          if (CMP(SSk[k][2][b][i],SSk[k][0][b][i-3]+pars.Omega)){
            transition=1;
            prev_state=0;
          }
        }
      }

      states[i]=curr_state;

      
      output[k].states[i]=curr_state;
      output[k].transitions[i]=transition;
      output[k].z[i]=z;

      curr_state=prev_state;

    }

    /*

    
    for (i=opt_b+2;i<=opt_i;i+=3){
      getBlock(i, seq_0, seq_k, map_0, map_k, block_0, block_k, &z );
      //printf("%s ",block_0);
    }

    //printf("\n");

    for (i=opt_b+2;i<=opt_i;i+=3){
      getBlock(i, seq_0, seq_k, map_0, map_k, block_0, block_k, &z );
      //printf("%s ",block_k);
    }

    display_line1 = (char*) malloc(sizeof(char)*(L+1)*3);
    display_line2 = (char*) malloc(sizeof(char)*(L+1)*3);
    display_line3 = (char*) malloc(sizeof(char)*(L+1)*3);
    display_line4 = (char*) malloc(sizeof(char)*(L+1)*3);
  
    for (i=0;i<(L+1)*3;i++){
      display_line1[i]=' ';
      display_line2[i]=' ';
      display_line3[i]=' ';
      display_line4[i]=' ';
    }
  
    display_line1[(L+1)*3-1]='\0';
    display_line2[(L+1)*3-1]='\0';
    display_line3[(L+1)*3-1]='\0';
    display_line4[(L+1)*3-1]='\0';

    pos=0;
    for (i=opt_b+2;i<=opt_i;i+=3){
      getBlock(i, seq_0, seq_k, map_0, map_k, block_0, block_k, &z );
      //sprintf(string, "%i  ", states[i]);
      strncpy(display_line1+pos,string,strlen(string));
      //sprintf(string, "%+.0f  ", SSk[k][0][opt_b][i]);
      strncpy(display_line2+pos,string,strlen(string));
      //sprintf(string, "%+.0f  ", SSk[k][1][opt_b][i]);
      strncpy(display_line3+pos,string,strlen(string));
      //sprintf(string, "%+.0f  ", SSk[k][2][opt_b][i]);
      strncpy(display_line4+pos,string,strlen(string));
        
      pos+=strlen(block_0)+1;
    }
  
    //printf("%s\n",display_line1);
    //printf("%s\n",display_line2);
    //printf("%s\n",display_line3);
    //printf("%s\n",display_line4);

    //printf("%i\n", L);

    */
    free(states);
    free(block_0);
    free(block_k);
    free(map_0);
    free(map_k);

  }


  return output;

}




/*********************************************************************
  getMultipleScoreMatrix

  Takes the pairwise score matrix Sk and calculates a two dimensional
  matrix [b,i] which represents the average over all pairs in the
  alignment.
  
*********************************************************************/ 

float** getMultipleScoreMatrix(float**** Sk, bgModel* models, const struct aln *alignment[]){
  
  int L,N,b,i,j,k;
  float sum, max;
  float** S;
  float s0,s1,s2;

  for (N=0; alignment[N]!=NULL; N++);

  L=getSeqLength(alignment[0]->seq);

  S=(float**)malloc(sizeof(float*)*(L+1));
  for (i=0;i<L+1;i++){
    S[i]=(float*)malloc(sizeof(float)*(L+1));
    for (j=0;j<L+1;j++){
      S[i][j]=0.0;
    }
  }

  for (b=1;b<L+1;b++){
    for (i=b+2;i<L+1;i+=3){

      sum=0;
      for (k=1;k<N;k++){
        sum+=MAX3(Sk[k][0][b][i],
                  Sk[k][1][b][i],
                  Sk[k][2][b][i]);
      }

      
      S[b][i]=MAX3(sum,
                   S[b][i-1]+pars.Delta,
                   S[b][i-2]+pars.Delta)/(N-1);
    }
  }

  return(S);
}

/*********************************************************************
  getHSS

  Takes score matrix S as calculated by getMultipleScore matrix and 
  searches for high scoring segments. 

  S .. score matrix 
  inputAln ... alignment 
  strand ... whether inputAln is forward or reverse complement ('+' or '-')
  parMu, parLambda ... if these parameters describing the extreme value distribution are given
  
  
*********************************************************************/ 

segmentStats* getHSS(float** S, const struct aln** inputAln, 
                     char strand){

  segmentStats* results;
  int sites,L, frame;
  int segmentStart, segmentEnd, currPos, hssCount;
  float currMax, currEntry;
  float pvalue;
  int i,j,k,l;
  int minSegmentLength=2;

  L=getSeqLength(inputAln[0]->seq);

  results=NULL;
  hssCount=0;

  for (frame=0;frame<=2;frame++){

    sites=((int)(L-frame)/3);

    currMax=0.0;
    segmentStart=-1;
    segmentEnd= -1;
      
    for (i=0; i<sites; i++){
      for (j=i;j<sites;j++){
        currEntry=S[i*3+1+frame][j*3+3+frame];

        if (currEntry>0.0 ||           /* only consider positive entries */
            (i==sites-1 && j==sites-1)){ /* enter on last entry in
                                            any case to report last found
                                            maximum */ 
          
          /* if maximum has been found previously and we have reached a point
             not overlapping with this segment, it is reported as
             local maximum (we report on last entry in any case...) */
          if ((currMax>0.0 && segmentEnd<i ) || (i==sites-1 && j==sites-1)){

            if (segmentEnd-segmentStart>=minSegmentLength){
            
              /* re-allocate for each new result, leave room for last entry */
              results=(segmentStats*)realloc(results,sizeof(segmentStats)*(hssCount+2));

              if (results==NULL){
                exit(1);
              }

              /* chromosome name is also stored in results, note that we
                 use statically allocated memory there */
            
              results[hssCount].name=strdup(inputAln[0]->name);
              results[hssCount].strand=strand;
              results[hssCount].frame=frame;
              results[hssCount].startSite=segmentStart;
              results[hssCount].endSite=segmentEnd;
              results[hssCount].score=currMax;

              results[hssCount].start=segmentStart*3+frame+1;
              results[hssCount].end=segmentEnd*3+frame+3;

              /* If CLUSTAL W input without coordinates */
              if ((inputAln[0]->start==0) && (inputAln[0]->length==0)){
              
                results[hssCount].startGenomic=results[hssCount].start;
                results[hssCount].endGenomic=results[hssCount].end;

              } else {
                if (strand=='+'){
                  results[hssCount].startGenomic=inputAln[0]->start+segmentStart*3+frame;
                  results[hssCount].endGenomic=inputAln[0]->start+segmentEnd*3+frame+2;
                } else {
                  results[hssCount].endGenomic=(inputAln[0]->start+inputAln[0]->length-1)-segmentStart*3-frame;
                  results[hssCount].startGenomic=(inputAln[0]->start+inputAln[0]->length-1)-segmentEnd*3-frame-2;
                }
              }

              results[hssCount].pvalue=pvalue;
              
              hssCount++;
            }
              
            currMax=currEntry; 
            segmentStart=i; 
            segmentEnd=j; 
          } 
          /* current potential maximal segment overlaps with
             previously found segment */
          else { 
            /* set new maximum if curr entry has better score or if the score is the same and it is longer */
            if (currEntry>currMax ||
                ((fabs(currEntry-currMax) < 0.0001)  && ((j-i)>=(segmentEnd-segmentStart)))){ 
              currMax=currEntry; 
              segmentStart=i; 
              segmentEnd=j; 
            }
          } 
        }
      } /* j */
    } /* i */
  } /* frame */
    
  if (hssCount==0){
    results=(segmentStats*)malloc(sizeof(segmentStats));
    results[0].pvalue=1.0;
  }

  results[hssCount].score=-1.0; /* mark end of list */
  
  return results;
 
}

int getExtremeValuePars(TTree* tree, const struct aln *alignment[], 
                         int sampleN, float maxNativeScore, float* parMu, float* parLambda){
  
  int L, i, j, betterThanNative, stopCutoff;
  float* freqsMono;
  float kappa;
  float sum, maxSum;
  double* maxScores;
  double mu,lambda;
  
  struct aln *sampledAln[MAX_NUM_NAMES];
  segmentStats *results, *resultsRev, *allResults;
  int sampleWritten=0;
  int hssCount;
  float** S;

  stopCutoff = (int)(pars.cutoff*pars.sampleN);

  L=strlen(alignment[0]->seq);

  freqsMono=models[0].freqs;

  kappa=models[0].kappa;

  maxScores=(double*)malloc(sizeof(double)*sampleN);

  betterThanNative=0;
  
  for (i=0;i<sampleN;i++){

    simulateTree(tree,freqsMono,kappa,L);
      
    tree2aln(tree,sampledAln);

    sortAln(alignment, sampledAln);
      
    reintroduceGaps(alignment, sampledAln);
      
    /* Write one random sampling, debugging */
    /*
    if (!sampleWritten){
      int x=0;
      FILE *fp = fopen("samples.maf", "a");
      for (x=1; sampledAln[x]!=NULL; ++x){
        sampledAln[x]->strand='+';
      }
      printAlnMAF(fp,(const struct aln**)sampledAln,0);
      fprintf(fp,"\n");
      fclose(fp);
      sampleWritten=1;
    }
    */

    results=scoreAln((const struct aln**)sampledAln, tree, kappa, 0);

    hssCount=0;
    while (results[hssCount].score>=0) hssCount++;

    qsort((segmentStats*) results, hssCount,sizeof(segmentStats),compareScores);

    if (results[0].score > maxNativeScore) {
      betterThanNative++;
    }

    if ((pars.stopEarly) && (betterThanNative > stopCutoff )) {
      return -1;
    }

    maxScores[i]=results[0].score;

    freeAln((struct aln**)sampledAln);
    freeResults(results);
  }

  if (EVDMaxLikelyFit(maxScores, NULL, sampleN, &mu, &lambda) == 1){
    *parMu=mu;
    *parLambda=lambda;
    free(maxScores);
    return 1;
  } else {
    /* On failure of fit return -1, which means that downstream all
       hits in this alignment are considered insignificant, should be
       enough error handling as the fit usually only fails for
       pathological cases that do not contain real hits.*/
    free(maxScores);
    return -1;
  }

}


segmentStats* scoreAln(const struct aln *inputAln[], TTree* tree, float kappa, int backtrack){
  
  struct aln *inputAlnRev[MAX_NUM_NAMES];
  segmentStats *results, *resultsRev, *allResults;
  float** S;
  int hssCount, i, N, L;
  
  for (N=0; inputAln[N]!=NULL; N++);   

  L = getSeqLength(inputAln[0]->seq);

  copyAln((struct aln**)inputAln,(struct aln**)inputAlnRev);
  revAln((struct aln**)inputAlnRev);

  getPairwiseScoreMatrix(models,(const struct aln**)inputAln);
  S=getMultipleScoreMatrix(Sk,models,(const struct aln**)inputAln);

  if (backtrack){
    if (Sk_native == NULL){
      Sk_native = allocateSk(N,L);
      Sk_native_rev = allocateSk(N,L);
    }
    copySk(Sk, Sk_native, N, L);
  }

  results=getHSS(S, (const struct aln**)inputAln, '+');

  freeS(S, (const struct aln **)inputAln);

  getPairwiseScoreMatrix(modelsRev,(const struct aln**)inputAlnRev);
  S=getMultipleScoreMatrix(Sk,modelsRev,(const struct aln**)inputAlnRev);

  if (backtrack){
    copySk(Sk, Sk_native_rev, N, L);
  }

  resultsRev=getHSS(S, (const struct aln**)inputAlnRev, '-');

  freeS(S, (const struct aln **)inputAln);
    
  hssCount=0;

  allResults=NULL;

  i=0;
  while (results[i].score > 0.0){
    allResults=(segmentStats*)realloc(allResults,sizeof(segmentStats)*(hssCount+2));
    allResults[hssCount]=results[i];
    allResults[hssCount].name=strdup(results[i].name);
    hssCount++;
    i++;
  }

  i=0;
  while (resultsRev[i].score > 0.0){
    allResults=(segmentStats*)realloc(allResults,sizeof(segmentStats)*(hssCount+2));
    allResults[hssCount]=resultsRev[i];
    allResults[hssCount].name=strdup(resultsRev[i].name);
    hssCount++;
    i++;
  }
  
  if (hssCount==0){
    allResults=(segmentStats*)malloc(sizeof(segmentStats));
    allResults[0].pvalue=1.0;
  }
    
  allResults[hssCount].score=-1.0; 



  freeResults(results);
  freeResults(resultsRev);

  //qsort((segmentStats*) allResults, hssCount,sizeof(segmentStats),compareScores);

  freeAln((struct aln**)inputAlnRev);

  return(allResults);

}


void freeS (float** S, const struct aln *alignment[]){

  int i,L;

  L=getSeqLength(alignment[0]->seq);

  for (i=0;i<L+1;i++){
    free(S[i]);
  }

  free(S);

}

void freeSk (float**** S, const struct aln *alignment[]){

  int L,N,k,i,x;

  L=getSeqLength(alignment[0]->seq);
  for (N=0; alignment[N]!=NULL; N++);

  for (k=0;k<N;k++){
    for (x=0;x<3;x++){
      for (i=0;i<L+1;i++){
        free(S[k][x][i]);
      }
      free(S[k][x]);
    }
    free(S[k]);
  }
  free(S);
}



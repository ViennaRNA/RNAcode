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
#include <sys/stat.h>
#include <errno.h>
#include "RNAcode.h"
#include "misc.h"
#include "postscript.h"

extern parameters pars;

extern long int hitCounter;

float**** allocateSk(int N, int L){

  float**** S;
  int i, k, x;
  
  // We have one entry for each pair
  S=(float****)malloc(sizeof(float***)*(N+1));

  for (k=0;k<N;k++){
    // We have three states
    S[k]=(float***)malloc(sizeof(float**)*(3));
    
    for (x=0;x<3;x++){
      // indices are 1 based and we mark end with NULL, so we need L+2
      S[k][x]=(float**)malloc(sizeof(float*)*(L+1));
      
      for (i=0;i<L+1;i++){
        S[k][x][i]=(float*)malloc(sizeof(float)*(L+1));
      }
    }
  }
  return S;
}

void copySk(float**** from, float**** to, int N, int L){

  int i, j, k, x;
  
  for (k=0;k<N;k++){
    for (x=0;x<3;x++){
      for (i=0;i<L+1;i++){
        for (j=0;j<L+1;j++){
          to[k][x][i][j] = from[k][x][i][j];
        }
      }
    }
  }

}




/*********************************************************************
  compareScores

  compare function for qsort that compares two scores given as pointer
  to a segmentStats structure.

*********************************************************************/ 

int compareScores(const void * a, const void * b){

  segmentStats* statsA;
  segmentStats* statsB;

  statsA=(segmentStats*)a;
  statsB=(segmentStats*)b;

  if  ( statsA->score < statsB->score) {
    return 1;
  } else {
    return -1;
  }

  return 0;
}

/*********************************************************************
  compareLocation

  compare function for qsort that compares start positions given as
  pointer to a segmentStats structure.

*********************************************************************/ 

int compareLocation(const void * a, const void * b){

  segmentStats* statsA;
  segmentStats* statsB;

  statsA=(segmentStats*)a;
  statsB=(segmentStats*)b;

  if  ( statsA->startSite > statsB->startSite) {
    return +1;
  } else {
    return -1;
  }

  return 0;
}


void reintroduceGaps(const struct aln* origAln[], struct aln* sampledAln[]){

  int i,j,k;

  char* tmpSeq;
  char* tmpName;
  char* origSeq;
  char* sampledSeq;


  for (i=0;origAln[i]!=NULL;i++){

    origSeq=origAln[i]->seq;
    sampledSeq=sampledAln[i]->seq;
    
    for (k=0;k<strlen(origSeq);k++){
      if (origSeq[k]=='-'){
        sampledSeq[k]='-';
      }
    }
  }
}

void sortAln(const struct aln* origAln[], struct aln* sampledAln[]){

  int i,j,k;

  char* tmpSeq;
  char* tmpName;
  char* origSeq;
  char* sampledSeq;

  for (i=0;origAln[i]!=NULL;i++){
    for (j=0;sampledAln[j]!=NULL;j++){
      if (strcmp(origAln[i]->name,sampledAln[j]->name)==0){
        tmpSeq=sampledAln[j]->seq;
        tmpName=sampledAln[j]->name;
        sampledAln[j]->seq=sampledAln[i]->seq;
        sampledAln[j]->name=sampledAln[i]->name;
        sampledAln[i]->seq=tmpSeq;
        sampledAln[i]->name=tmpName;
      }
    }
  }
}

/*

returns block of potential codon that ends in i; the block contains leading gaps

seq_0    ... whole reference sequence
seq_k    ... whole sequence k
i        ... sequence position in reference sequence (1-based)
block_0  ... Block in reference sequence that is to be returned
block_k  ... Block in sequence k that is to be returned
z        ... difference of gaps between reference and k in the block

*/

void getBlock(int i, const char* seq_0, const char* seq_k, const int* map_0, const int* map_k, char* block_0, char* block_k, int* z ){

  int start,end;
  int gap_0, gap_k;
  int diff;
  char c;

  if (i<3){
    fprintf(stderr, "Fatal error in getBlock. i needs to be >=3");
    exit(0);
  }

  if (i>3){
    //start=pos2col(seq_0,i-3)+1;
    start=map_0[i-3]+1;
  }

  if (i==3){
    start=1;
  }
  
  //end=pos2col(seq_0,i);
  end=map_0[i];

  strncpy(block_0, seq_0+start-1,end-start+1);
  strncpy(block_k, seq_k+start-1,end-start+1);
  
  block_0[end-start+1]='\0';
  block_k[end-start+1]='\0';

  gap_0 = gap_k = 0;

  i=0;
  while (block_0[i] != '\0'){
    if (block_0[i] == '-') gap_0++;
    i++;
  }

  i=0;
  while (block_k[i] != '\0'){
    if (block_k[i] == '-') gap_k++;
    i++;
  }

  diff=gap_k-gap_0;
  if (diff<0) diff*=-1;


  if ( diff % 3 == 0 ){
    *z=0;
  }

  if ( diff % 3 == 1 ){
    *z=+1;
  }

  if ( diff % 3 == 2 ){
    *z=-1;
  }

}

// both pos and col are 1 based!

int pos2col(const char* seq, int pos){

  int i=0;
  int currPos=0;
    
  while (seq[i] != '\0'){
  
    if (seq[i]=='-'){
      i++;
      continue;
    } else {
      currPos++;
    }

    if (currPos==pos){
      return i+1;
    }
    i++;
  }
}


int getSeqLength(char* seq){

  int i=0;
  int counter=0;
  
  while (seq[i] != '\0'){
    if (seq[i]=='-'){
      i++;
      continue;
    } else {
      counter++;
      i++;
    }
  }

  return counter;

}


void freeResults(segmentStats results[]){
  
  int i=0;

  while (results[i].score>0){
    free(results[i].name);
    i++;
  }
  free(results);
}

int hDist(int a1, int a2, int a3, int b1, int b2, int b3){

  int dist=0;

  if (a1 != b1) dist++;
  if (a2 != b2) dist++;
  if (a3 != b3) dist++;
 
  return dist;

}


float avg(float* data, int N){

  int i;
  float sum;

  sum=0.0;

  for (i=0;i<N;i++){
    sum+=data[i];
  }
  
  return sum/(float)N;

}


float stddev(float* data, int N){

  int i;
  float sum;
  float mean;

  sum=0.0;

  for (i=0;i<N;i++){
    sum+=data[i];
  }

  mean=sum/(float)N;

  sum=0.0;

  for (i=0;i<N;i++){
    sum+=(mean-data[i])*(mean-data[i]);
  }
  
  return sqrt(sum/(float)(N-1));
}

void copyAln(struct aln *src[],struct aln *dest[]){

  int i,j,L;
  char* seq;
  
  L=strlen(src[0]->seq);

  for (i=0;src[i]!=NULL;i++){
    
    dest[i]=createAlnEntry(strdup(src[i]->name),
                           strdup(src[i]->seq),
                           src[i]->start, 
                           src[i]->length,
                           src[i]->fullLength,
                           src[i]->strand);
  }

  dest[i]=NULL;
}



void printAlnClustal(FILE *out, const struct aln* AS[]){

  int i;
  
  fprintf(out,"CLUSTAL W\n\n");

  for (i=0;AS[i]!=NULL;i++){
    fprintf(out, "%s %s\n",AS[i]->name,AS[i]->seq);
  }

  fprintf(out, "\n");

}


void printResults(FILE* outfile, int outputFormat, const struct aln* inputAln[], segmentStats results[]){

  int i,k, hssCount, currHSS, nextHSS;
  char c;
  char name[1024]="";
  char fileName[2048]="";
  char prefix[1024]="";
  char suffix[1024]="";

  hssCount = 0;
  while (results[hssCount].score > 0.0){
    results[hssCount].hide=0;
    //printf("Inner %.2f (%i)\n", results[hssCount].score, results[hssCount].hide);
    hssCount++;
  }
    
  if (pars.bestRegion){
  
    /* First sort by start position */
    qsort((segmentStats*) results, hssCount,sizeof(segmentStats),compareLocation);
    currHSS=0;
    nextHSS=1;
    while (nextHSS<=hssCount) {
    
      /* the two HSS overlap */
      if (!(results[currHSS].endSite <= results[nextHSS].startSite)){
      
        /* Hide the HSS with lower score */
        if (results[currHSS].score > results[nextHSS].score){
          results[nextHSS].hide=1;
          nextHSS++;
        } else {
          results[currHSS].hide=1;
          currHSS=nextHSS;
          nextHSS++;
        }
      } else {
        currHSS=nextHSS;
        nextHSS++;
      }
    }
  }

  qsort((segmentStats*) results, hssCount,sizeof(segmentStats),compareScores);

  if (results[0].score<0.0 || results[0].pvalue > pars.cutoff){
    if (outputFormat==0){
      fprintf(outfile,"\nNo significant coding regions found.\n");
    }
    return;
  } 

  if (outputFormat==0){
    fprintf(outfile, "\n%6s%5s%7s%6s%6s%12s%12s%12s%9s%9s\n",
            " HSS # ", "Frame","Length","From","To","Name","Start","End", "Score","P");
    fprintf(outfile, "======================================================================================\n");
  }

  
  i=0;

  while (results[i].score>0.0 && results[i].pvalue < pars.cutoff){
    
    if (results[i].hide) {
      i++;
      continue;
    }


    if (pars.postscript){
      if (results[i].pvalue < pars.postscript_cutoff){
        struct stat stat_p;	
        
        if (stat (pars.postscriptDir, &stat_p) != 0){
          if (mkdir(pars.postscriptDir, S_IRWXU|S_IROTH|S_IRGRP ) !=0){
            fprintf(stderr, "WARNING: Could not create directory: %s", pars.postscriptDir);
          }
        }
        
        sprintf(fileName,"%s/hss-%li.eps", pars.postscriptDir,hitCounter); 
        colorAln(fileName,(const struct aln**)inputAln, results[i]);
      }
    }


    if (outputFormat==0){

      fprintf(outfile, "%6li %4c%i%7i%6i%6i%12s%12i%12i%9.2f",
              hitCounter,
              results[i].strand, results[i].frame+1,
              results[i].endSite-results[i].startSite+1,
              results[i].startSite+1,results[i].endSite+1,
              results[i].name,
              results[i].startGenomic,results[i].endGenomic,
              results[i].score);

        if (results[i].pvalue < 0.001){
          /*Seems to be the minimum number I can get, don't know why
            we don't get down to 1e-37 which should be the limit for
            floats.
          */
          if (results[i].pvalue < 10e-16){
            fprintf(outfile, "   <1e-16\n");
          } else {
            fprintf(outfile, "% 9.1e\n",results[i].pvalue);
          }

        } else {
          fprintf(outfile, "% 9.3f\n",results[i].pvalue);
        }
    }

    if (outputFormat==1){
      /* if name is of the form hg18.chromX than only display chromX */
      k=0;
      while (1){
        if (results[i].name[k]=='\0' || results[i].name[k]=='.'){
          break;
        }
        k++;
      }

      if (k==strlen((char*)results[i].name)){
        strcpy(name,results[i].name);
      } else {
        strcpy(name,results[i].name+k+1);
      }

      fprintf(outfile,"%s\t%s\t%s\t%i\t%i\t%.2f|%.2e\t%c\t%c\t%s%li%s\n",
              name, "RNAcode","CDS",
              results[i].startGenomic+1,results[i].endGenomic+1,
              results[i].score,
              results[i].pvalue,
              results[i].strand, '.',"gene_id \"Gene", hitCounter ,"\"; transcript_id \"transcript 0\";");
    }
  
    if (outputFormat==2){

      fprintf(outfile, "%li\t%c\t%i\t%i\t%i\t%i\t%s\t%i\t%i\t%7.3f\t",
              hitCounter,
              results[i].strand, results[i].frame+1,
              results[i].endSite-results[i].startSite+1,
              results[i].startSite+1,results[i].endSite+1,
              results[i].name,
              results[i].startGenomic,results[i].endGenomic,
              results[i].score);
      
      if (results[i].pvalue < 0.001){
        fprintf(outfile, "% 9.3e\n",results[i].pvalue);
      } else {
        fprintf(outfile, "% 9.3f\n",results[i].pvalue);
      }
    }
    i++;

    if (pars.bestOnly) break;

    hitCounter++;

  }
}


int extendRegion(const struct aln* alignment[], int pos, int direction){

  int *map_0, *map_k;
  int N, k, x, z, L, l, colsN;
  char *block_0, *block_k, *seq;
  int pepA;
  int ii, jj;
  char codonA[4];
  
  seq = alignment[0]->seq;

  L=getSeqLength(seq);
  colsN=strlen(seq);

  block_0 = (char*) malloc(sizeof(char)*(colsN+1));
  block_k = (char*) malloc(sizeof(char)*(colsN+1));
  
  map_0=(int*)malloc(sizeof(int)*(colsN+1));
  map_k=(int*)malloc(sizeof(int)*(colsN+1));

  //printf("Pos %i L %i\n", pos, L);

  for (l=1;l<=L;l++){
    map_0[l]=pos2col(seq,l);
    map_k[l]=pos2col(alignment[1]->seq,l);
  }
  
  if (direction == 0 ){
    x=pos+2;
  } else {
    x=pos;
  }

  while (1){
    
    getBlock(x, seq, alignment[1]->seq, map_0, map_k, block_0, block_k, &z );

    ii=jj=0;
      
    while (block_0[ii] != '\0') {
      if (block_0[ii] != '-') {
        codonA[jj]=block_0[ii];
        jj++;
      }
      ii++;
    }

    pepA=transcode[ntMap[codonA[0]]][ntMap[codonA[1]]][ntMap[codonA[2]]];

    if (pepA == -1) break;

    if (direction == 0){
      if (x-3<3) break;
      x-=3;

    } else {
      if (x+3 > L) break;
      x+=3;
    }
  }

  free(block_0);
  free(block_k);
  free(map_0);
  free(map_k);




  if (direction == 0){
    return x-2;
  } else {
    return x;
  }


}

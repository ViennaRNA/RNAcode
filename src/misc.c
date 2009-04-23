#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "misc.h"

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

void getBlock(const char* seq_0, const char* seq_k, int i, char* block_0, char* block_k, int* z ){

  int start,end;
  int gap_0, gap_k;
  int diff;
  char c;

  if (i<3){
    fprintf(stderr, "Fatal error in getBlock. i needs to be >=3");
    exit(0);
  }

  if (i>3){
    start=pos2col(seq_0,i-3)+1;
  }

  if (i==3){
    start=1;
  }
  
  end=pos2col(seq_0,i);

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


double avg(double* data, int N){

  int i;
  double sum;

  sum=0.0;

  for (i=0;i<N;i++){
    sum+=data[i];
  }
  
  return sum/(double)N;

}


double stddev(double* data, int N){

  int i;
  double sum;
  double mean;

  sum=0.0;

  for (i=0;i<N;i++){
    sum+=data[i];
  }

  mean=sum/(double)N;

  sum=0.0;

  for (i=0;i<N;i++){
    sum+=(mean-data[i])*(mean-data[i]);
  }
  
  return sqrt(sum/(double)(N-1));
}


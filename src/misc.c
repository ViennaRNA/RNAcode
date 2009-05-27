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


void printResults(FILE* outfile, int outputFormat, segmentStats results[]){

  int i,k;
  char c;
  char name[1024]="";
  char prefix[1024]="";
  char suffix[1024]="";

  if (outputFormat==0){

    if (results[0].score<0.0){
      fprintf(outfile,"No significant coding regions found.\n");
    } else {

      fprintf(outfile, "\n%5s%7s%6s%6s%12s%12s%12s%9s%9s\n",
              "Frame","Length","From","To","Name","Start","End", "Score","P");
  
      fprintf(outfile, "==============================================================================\n");

      i=0;

      while (results[i].score>0){

        fprintf(outfile, "%4c%i%7i%6i%6i%12s%12i%12i%9.2f",
                results[i].strand, results[i].frame,
                results[i].endSite-results[i].startSite+1,
                results[i].startSite,results[i].endSite,
                results[i].name,
                results[i].start,results[i].end,
                results[i].score);

        if (results[i].pvalue < 0.001){
          fprintf(outfile, "% 9.1e\n",results[i].pvalue);
        } else {
          fprintf(outfile, "% 9.3f\n",results[i].pvalue);
        }

        i++;
      }
    }
  }

  
  if (outputFormat==1){
    i=0;
    while (results[i].score>0){
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

      fprintf(outfile,"%s\t%s\t%s\t%i\t%i\t%.2f|%.2e\t%c\t%c\t%s\n",
              name, "RNAcode","CDS",
              results[i].start+1,results[i].end+1,
              results[i].score,
              results[i].pvalue,
              results[i].strand, '.',"gene_id \"Gene 0\"; transcript_id \"transcript 0\";");
      i++;

      /* GTF, currently only outputs highest scoring hit */
      if (i > 0) break;
    }
  }

  if (outputFormat==2){

    i=0;

    while (results[i].score>0){
      fprintf(outfile, "%c\t%i\t%i\t%i\t%i\t%s\t%i\t%i\t%7.3f\t",
              results[i].strand, results[i].frame,
              results[i].endSite-results[i].startSite+1,
              results[i].startSite,results[i].endSite,
              results[i].name,
              results[i].start,results[i].end,
              results[i].score);

      if (results[i].pvalue < 0.001){
        fprintf(outfile, "% 9.3e\n",results[i].pvalue);
      } else {
        fprintf(outfile, "% 9.3f\n",results[i].pvalue);
      }
      break;
      i++;
    }
  }
}

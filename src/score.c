#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "score.h"
#include "code.h"
#include "misc.h"

double sigma=+4.0;
double omega=-2.0;
double Omega=-4.0;
double Delta=-10.0;

extern int transcode[4][4][4];
extern int BLOSUM62[24][24];


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
      matrix[i][j]=BLOSUM62[i][j];
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
  double f, prob, score;
  double counts[4];
  int h, i;

  counts[0]=counts[1]=counts[2]=counts[3]=0.0;
  model->scores[0]=model->scores[1]=model->scores[2]=model->scores[3]=0.0;

  for (a1=0;a1<4;a1++){
    for (a2=0;a2<4;a2++){
      for (a3=0;a3<4;a3++){
        for (b1=0;b1<4;b1++){
          for (b2=0;b2<4;b2++){
            for (b3=0;b3<4;b3++){
              
              h=hDist(a1,a2,a3,b1,b2,b3);

              f=(model->freqs[a1])*(model->freqs[a2])*(model->freqs[a3]);
              
              prob=probHKY(a1,b1,model->dist,model->freqs,model->kappa);

              prob*=probHKY(a2,b2,model->dist,model->freqs,model->kappa);
              prob*=probHKY(a3,b3,model->dist,model->freqs,model->kappa);

              prob*=f;
              
              pepA=transcode[a1][a2][a3];
              pepB=transcode[b1][b2][b3];
              
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


double probHKY(int i, int j, double d, double freqs[4], double kappa){

  double piA, piC, piG, piT;
  double piR, piY, r, l, k1, k2, exp1, exp2, exp22, exp21;
  double result[4][4];


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

void countFreqsMono(const struct aln *alignment[], double freqs[]){

  int i,k;
  char* currSeq;
  char c;
  unsigned long counter;
  double sum;

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

  for (i=0;i<4;i++) freqs[i]/=(double)counter;

}


bgModel* getModels(TTree* tree, struct aln *alignment[], double kappa){

  int i,j,N;
  int **scoringMatrix;
  double** distanceMatrix;
  double freqsMono[4];
  bgModel* models;

  for (N=0; alignment[N]!=NULL; N++);   
  
  distanceMatrix=getDistanceMatrix(tree,(struct aln**)alignment);
  
  countFreqsMono((const struct aln**)alignment, (double *) freqsMono);

  models=(bgModel*)malloc(sizeof(bgModel)*N);

  /*for (i=0;i<N;i++){
    models[i]=(bgModel*)malloc(sizeof(bgModel)*N);
  }
  */

  //i=0;
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

void freeModels(bgModel* models, int N){

  int i, j;

  for (j=0;j<N;j++){
    freeScoringMatrix(models[j].matrix);
  }

  free(models);
}



double calculateSigma(char* block_0, char* block_k, int k, bgModel* models){

  char codonA[4];
  char codonB[4];
  int pepA, pepB;
  int i,j,h;
  double currScore,expectedScore, observedScore;

  i=j=0;
  while (block_0[i] != '\0') {
    if (block_0[i]!='-') codonA[j++]=block_0[i++];
  }
  i=j=0;
  while (block_k[i] != '\0') {
    if (block_k[i]!='-') codonB[j++]=block_k[i++];
  }
  
  codonA[3]=codonB[3]='\0';

  h=hDist(ntMap[codonA[0]], ntMap[codonA[1]], ntMap[codonA[2]],
          ntMap[codonB[0]], ntMap[codonB[1]], ntMap[codonB[2]]);

  if (h==0) return 0.0;

  pepA=transcode[ntMap[codonA[0]]][ntMap[codonA[1]]][ntMap[codonA[2]]];
  pepB=transcode[ntMap[codonB[0]]][ntMap[codonB[1]]][ntMap[codonB[2]]];

  if (pepA==-1) pepA=23;
  if (pepB==-1) pepB=23;
  
  expectedScore=models[k].scores[h];
  observedScore=(double)models[k].matrix[pepA][pepB];

  //printf("codonA: %s %s %i %i %.2f %.2f %.2f\n",codonA, codonB, pepA, pepB, observedScore, expectedScore, observedScore-expectedScore );

  //return 4.0;

  return observedScore-expectedScore;
}

double** getMultipleScoreMatrix(double**** Sk, bgModel* models, const struct aln *alignment[]){
  
  int L,N,b,i,j,k;
  double sum, max;
  double** S;

  for (N=0; alignment[N]!=NULL; N++);

  L=getSeqLength(alignment[0]->seq);

  S=(double**)malloc(sizeof(double*)*(L+1));
  for (i=0;i<L+1;i++){
    S[i]=(double*)malloc(sizeof(double)*(L+1));
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

       
        S[b][i]=MAX3(sum,
                     S[b][i-1]+Delta,
                     S[b][i-2]+Delta);

        //printf("%.2f ",S[b][i]);
      }
    }
    //printf("\n");
  }
  
  return(S);
}

segmentStats* getHSS(double** S, const struct aln** inputAln, double parMu, double parLambda, double cutoff){

  segmentStats* results;
  int sites,L, frame, strand;
  int segmentStart, segmentEnd, currPos, hssCount;
  double currMax, currEntry;
  double pvalue;
  int i,j,k,l;
  int minSegmentLength=2;

  L=getSeqLength(inputAln[0]->seq);

  results=NULL;
  hssCount=0;

  strand=0;

  for (frame=0;frame<=2;frame++){

    sites=((int)(L-frame)/3);

    currMax=0.0;
    segmentStart=-1;
    segmentEnd= -1;
      
    for (i=0; i<sites; i++){
      for (j=i;j<sites;j++){
        currEntry=S[i*3+1+frame][j*3+3+frame];

        //printf("-------------------> Frame %i: %i %i %.2f (%i) \n", frame, i,j, currEntry, sites);

        if (currEntry>0.0 ||           /* only consider positive entries */
            (i==sites-1 && j==sites-1)){ /* enter on last entry in
                                            any case to report last found
                                            maximum */ 
          
          /* if maximum has been found previously and we have reached a point
             not overlapping with this segment, it is reported as
             local maximum (we report on last entry in any case...) */
          if ((currMax>0.0 && segmentEnd<i ) || (i==sites-1 && j==sites-1)){

            //printf("%u,%u:%.2f\n",segmentStart,segmentEnd,currMax);
          
            /* If both parameters are set to zero, no p-values are
               calculated and all positive scores are reported */
            if (parLambda==0.0 && parMu==0.0){
              pvalue=0.0;
            } else {
              pvalue=1-exp((-1)*exp((-1)*parLambda*(currMax-parMu)));
            }

            //printf("------------------->  %i: %i %i %.2f %.2f \n", frame, segmentStart, segmentEnd, currMax, cutoff);
            
            if (pvalue<cutoff && segmentEnd-segmentStart>=minSegmentLength){
             
              /* re-allocate for each new result, leave room for last entry */
              results=(segmentStats*)realloc(results,sizeof(segmentStats)*(hssCount+2));


              if (results==NULL){
                exit(1);
              }

              /* chromosome name is also stored in results, note that we
                 use statically allocated memory there */
            
              results[hssCount].name=strdup(inputAln[0]->name);
              results[hssCount].strand='+';
              results[hssCount].frame=frame;
              results[hssCount].startSite=segmentStart;
              results[hssCount].endSite=segmentEnd;
              results[hssCount].score=currMax;

              /* If CLUSTAL W input without coordinates */
              if ((inputAln[0]->start==0) && (inputAln[0]->length==0)){
              
                results[hssCount].start=segmentStart*3+frame;
                results[hssCount].end=segmentEnd*3+frame+2;

              } else {
                
                if (strand==0){
                  results[hssCount].start=inputAln[0]->start+segmentStart*3+frame;
                  results[hssCount].end=inputAln[0]->start+segmentEnd*3+frame+2;
                  
                } else {
                    
                  results[hssCount].end=(inputAln[0]->start+inputAln[0]->length-1)-segmentStart*3-frame;
                  results[hssCount].start=(inputAln[0]->start+inputAln[0]->length-1)-segmentEnd*3-frame-2;
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
            /* set new maximum if curr entry is larger */
            if (currEntry>currMax){ 
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

  results[hssCount].score=-1; /* mark end of list */
  
  return results;
 

}


void printResults(FILE* outfile, int outputFormat, segmentStats results[]){

  char direction;
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

        direction='+';
        if (results[i].strand==1)  direction='-';

   
        fprintf(outfile, "%4c%i%7i%6i%6i%12s%12i%12i%9.2f",
                direction, results[i].frame,
                results[i].endSite-results[i].startSite+1,
                results[i].startSite,results[i].endSite,
                results[i].name,
                results[i].start,results[i].end,
                results[i].score);

        if (results[i].pvalue < 0.001){
          fprintf(outfile, "% 9.1e\n",results[i].pvalue);
        } else {
          fprintf(outfile, "% 9.2f\n",results[i].pvalue);
        }


        i++;
      }
    }
  }

  
  if (outputFormat==1){
    i=0;
    while (results[i].score>0){
      direction='+';
      if (results[i].strand==1)  direction='-';
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
              direction, '.',"gene_id \"Gene 0\"; transcript_id \"transcript 0\";");
      i++;

      /* GTF, currently only outputs highest scoring hit */
      if (i > 0) break;
    }
  }

  if (outputFormat==2){

    i=0;

    while (results[i].score>0){
      direction='+';
      if (results[i].strand==1)  direction='-';
   
      fprintf(outfile, "%c\t%i\t%i\t%i\t%i\t%s\t%i\t%i\t%7.3f",
              direction, results[i].frame,
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


void getExtremeValuePars(TTree* tree,bgModel* models, const struct aln *alignment[], 
                         int sampleN, double* parMu, double* parLambda){

  int tmpCounter, L, i, j;
  double* freqsMono;
  double kappa;
  double sum, maxSum;
  double* maxScores;
  double mu, sigma;
  struct aln *sampledAln[MAX_NUM_NAMES];
  segmentStats *results;
  int sampleWritten=0;
  int hssCount;
  double**** Sk;
  double** S;

  L=strlen(alignment[0]->seq);

  freqsMono=models[0].freqs;
  kappa=models[0].kappa;

  maxScores=(double*)malloc(sizeof(double)*sampleN);

  //FILE *fp = fopen("scores.dat", "a");

  tmpCounter=0;

  for (i=0;i<sampleN;i++){
  
    simulateTree(tree,freqsMono,kappa,L);
      
    tree2aln(tree,sampledAln);

    sortAln(alignment, sampledAln);
      
    reintroduceGaps(alignment, sampledAln);
      
    /* Write one random sampling, debugging */
    /* if (!sampleWritten){ */
    /*   int x=0; */
    /*   FILE *fp = fopen("samples.maf", "a"); */
    /*   for (x=1; sampledAln[x]!=NULL; ++x){ */
    /*     sampledAln[x]->strand='+'; */
    /*   } */
    /*   printAlnMAF(fp,(const struct aln**)sampledAln,0); */
    /*   fprintf(fp,"\n"); */
    /*   fclose(fp); */
    /*   sampleWritten=1; */
    /* } */

    //printAlnMAF(stdout,(const struct aln**)sampledAln,0);

    Sk=getPairwiseScoreMatrix(models,(const struct aln**)sampledAln);
    S=getMultipleScoreMatrix(Sk,models,(const struct aln**)sampledAln);

    results=getHSS(S, (const struct aln**)sampledAln, 0.0, 0.0, 1.0);

    hssCount=0;
    while (results[hssCount].score>=0) hssCount++;

    //printf("%.2f\n", results[0].score);

    qsort((segmentStats*) results, hssCount,sizeof(segmentStats),compareScores);

    maxScores[i]=results[0].score;
        
    freeAln((struct aln**)sampledAln);
    free(results);
  }

  EVDMaxLikelyFit(maxScores, NULL, sampleN, parMu, parLambda);

  free(maxScores);
  //fclose(fp);

}




// Calculates matrix S[k][x][b][i]
// k ... number of sequence
// x ... State (0 ... '0', 1 ... '+1', 2 ... '-1')
// b,i ... subsequence from b to i

double**** getPairwiseScoreMatrix(bgModel* models, const struct aln *alignment[]){
  
  int L,N, pos,z,b,i,x, k;
  char *block_0, *block_k, *seq_0, *seq_k;
  double**** S;

  seq_0=alignment[0]->seq;

  L=getSeqLength(seq_0);

  for (N=0; alignment[N]!=NULL; N++);

  // We have one entry for each pair
  S=(double****)malloc(sizeof(double***)*(N+1));
  
  for (k=0;k<N;k++){
    // We have three states
    S[k]=(double***)malloc(sizeof(double**)*(3));

    for (x=0;x<3;x++){
      // indices are 1 based and we mark end with NULL, so we need L+2
      S[k][x]=(double**)malloc(sizeof(double*)*(L+2));
        
      for (i=0;i<L+2;i++){
        S[k][x][i]=(double*)malloc(sizeof(double)*(L+1));
      }
      S[k][x][L+1]=NULL;

      for (b=0;b<L+1;b++){
        for (i=0;i<L+1;i++){
          if (abs(b-i)<3 ){
            S[k][x][b][i]=0.0;
          } else {
            S[k][x][b][i]=0.0;
          }
        }
      }
    }
  }

  // Allocate conservatively for full length of sequence;
  block_0 = (char*) malloc(sizeof(char)*(L+1));
  block_k = (char*) malloc(sizeof(char)*(L+1));

  //printf("%s <-REF\n",seq_0);

  for (k=1;k<N;k++){

    seq_k=alignment[k]->seq;

    //printf("%s\n",seq_k);
      
    for (b=1;b<L+1;b++){
      for (i=b+2;i<L+1;i+=3){
        getBlock(seq_0, seq_k, i, block_0, block_k, &z );
        
        if (z==0){
          S[k][0][b][i]=S[k][0][b][i-3]+calculateSigma(block_0,block_k,k,models);
          S[k][1][b][i]=S[k][1][b][i-3]+omega;
          S[k][2][b][i]=S[k][2][b][i-3]+omega;
        }

        if (z==+1){
          S[k][0][b][i]=MAX(S[k][0][b][i-3]+Delta,
                            S[k][2][b][i-3]+Omega);
          
          S[k][1][b][i]=MAX(S[k][0][b][i-3]+Omega,
                            S[k][1][b][i-3]+Delta);
        
          S[k][2][b][i]=MAX(S[k][1][b][i-3]+Omega,
                            S[k][2][b][i-3]+Delta);
        }
        
        if (z==-1){
          S[k][0][b][i]=MAX(S[k][0][b][i-3]+Delta,
                            S[k][1][b][i-3]+Omega);
        
          S[k][1][b][i]=MAX(S[k][1][b][i-3]+Delta,
                            S[k][2][b][i-3]+Omega);
        
          S[k][2][b][i]=MAX(S[k][2][b][i-3]+Delta,
                            S[k][0][b][i-3]+Omega);
        
        }
      }
    }
  }

  /*
    for (x=0; x<3;x++){
    for (i=0;i<L+1;i++){
    free(S[x][i]);
    }
    free(S[x]);
    }
  */
  
  return(S);

}

double* backtrack(double**** S, int k, const struct aln *alignment[]){

  double opt_score=0.0;
  int opt_state, opt_b, opt_i; 
  int curr_state, prev_state;
  char *display_line1;
  char *display_line2;
  char *display_line3;
  char *display_line4;

  char* block_0;
  char* block_k;
  int* states;
  
  char* seq_0;
  char* seq_k;

  seq_0=alignment[0]->seq;
  seq_k=alignment[k]->seq;
  
  char string[1000];

  int pos,i,b,x,z,L;

  L=getSeqLength(seq_0);

  //for (i=1; S[k][0][i]!=NULL; i++) L++;

  for (b=1;b<L+1;b++){
    for (i=1;i<L+1;i++){
      if (i>=b){
        for (x=0;x<3;x++){
          if (S[k][x][b][i] > opt_score){
            opt_score=S[k][0][b][i];
            opt_i=i;
            opt_b=b;
            opt_state=x;
          }
        }
        //printf("%+.0f/%+.0f/%+.0f\t",S[k][0][b][i], S[k][1][b][i], S[k][2][b][i]);
      } else {
        //printf("    -   \t");
      }
    }
    //printf("\n\n");
  }


  printf("Max score: %.1f at b=%i, i=%i at state %i\n", opt_score, opt_b, opt_i, opt_state);

   states=(int*)malloc(sizeof(int)*(L+3));
      
  for (i=0;i<=L+1;i++) states[i]=-1;
  
  // Allocate conservatively for full length of sequence;
  block_0 = (char*) malloc(sizeof(char)*(L+1));
  block_k = (char*) malloc(sizeof(char)*(L+1));
      
  curr_state=opt_state;
  b=opt_b;
  
  for (i=opt_i;i>=opt_b+2;i-=3){
    getBlock(seq_0, seq_k, i, block_0, block_k, &z );
    
    printf("%i %i %i %i\n",b, i, z, curr_state);

    if (z==0){
      prev_state = curr_state;
    }
        
    if (z==+1){
      
      if (curr_state == 0){
        if (S[k][0][b][i]==S[k][0][b][i-3]+Delta) prev_state=0;
        if (S[k][0][b][i]==S[k][2][b][i-3]+Omega) prev_state=2;
      }
          
      if (curr_state == 1){
        if (S[k][1][b][i]==S[k][0][b][i-3]+Omega) prev_state=0;
        if (S[k][1][b][i]==S[k][1][b][i-3]+Delta) prev_state=1;
      }

      if (curr_state == 2){
        if (S[k][1][b][i]==S[k][1][b][i-3]+Omega) prev_state=1;
        if (S[k][1][b][i]==S[k][2][b][i-3]+Delta) prev_state=2;
      }
    }

        
    if (z==-1){

      if (curr_state == 0){
        if (S[k][0][b][i]==S[k][0][b][i-3]+Delta) prev_state=0;
        if (S[k][0][b][i]==S[k][1][b][i-3]+Omega) prev_state=1;
      }
          
      if (curr_state == 1){
        if (S[k][1][b][i]==S[k][1][b][i-3]+Delta) prev_state=1;
        if (S[k][1][b][i]==S[k][2][b][i-3]+Omega) prev_state=2;
      }

      if (curr_state == 2){
        if (S[k][1][b][i]==S[k][2][b][i-3]+Delta) prev_state=2;
        if (S[k][1][b][i]==S[k][0][b][i-3]+Omega) prev_state=0;
      }
    }
        
    states[i]=curr_state;

    curr_state=prev_state;
                
  }

  for (i=opt_b+2;i<=opt_i;i+=3){
    getBlock(seq_0, seq_k, i, block_0, block_k, &z );
    printf("%s ",block_0);
  }

  printf("\n");

  for (i=opt_b+2;i<=opt_i;i+=3){
    getBlock(seq_0, seq_k, i, block_0, block_k, &z );
    printf("%s ",block_k);
  }
  printf("\n");
      
  /* allocate a string to display annotations 3 times longer than
     the sequence to allow for spaces between codon in display */
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
    getBlock(seq_0, seq_k, i, block_0, block_k, &z );
    sprintf(string, "%i  ", states[i]);
    strncpy(display_line1+pos,string,strlen(string));
    sprintf(string, "%+.0f  ", S[k][0][opt_b][i]);
    strncpy(display_line2+pos,string,strlen(string));
    sprintf(string, "%+.0f  ", S[k][1][opt_b][i]);
    strncpy(display_line3+pos,string,strlen(string));
    sprintf(string, "%+.0f  ", S[k][2][opt_b][i]);
    strncpy(display_line4+pos,string,strlen(string));
        
    pos+=strlen(block_0)+1;
  }
  
  printf("%s\n",display_line1);
  printf("%s\n",display_line2);
  printf("%s\n",display_line3);
  printf("%s\n",display_line4);

  printf("%i\n", L);


}




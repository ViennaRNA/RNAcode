#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "code.h"

/*                        AAA AAC AAG AAT   ACA ACC ACG ACT   AGA AGC AGG AGT   ATA ATC ATG ATT   */
/*                         K   N   K   N     T   T   T   T     R   S   R   S     I   I   M   I    */
int transcode[4][4][4]={{{ 11, 2 , 11, 2 },{ 16, 16, 16, 16},{ 1 , 15, 1 , 15},{ 9 , 9 , 12, 9 }},
/*                        CAA CAC CAG CAT   CCA CCC CCG CCT   CGA CGC CGG CGT   CTA CTC CTG CTT   */  
/*                         Q   H   Q   H     P   P   P   P     R   R   R   R     L   L   L   L    */
                        {{ 5 , 8 , 5 , 8 },{ 14, 14, 14, 14},{ 1 , 1 , 1 , 1 },{ 10, 10, 10, 10 }},
/*                        GAA GAC GAG GAT   GCA GCC GCG GCT   GGA GGC GGG GGT   GTA GTC GTG GTT   */  
/*                         E   D   E   D     A   A   A   A     G   G   G   G     V   V   V   V    */
                        {{ 6 , 3 , 6 , 3 },{ 0 , 0 , 0 , 0 },{ 7 , 7 , 7 , 7 },{ 19, 19, 19, 19}},
/*                         *   Y   *   Y     S   S   S   S     *   C   W   C     L   F   L   F    */
/*                        TAA TAC TAG TAT   TCA TCC TCG TCT   TGA TGC TGG TGT   TTA TTC TTG TTT   */  
                        {{ -1, 18, -1, 18},{ 15, 15, 15, 15},{ -1, 4 , 17, 4 },{ 10, 13, 10, 13}}};


int BLOSUM62[24][24]={{ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4},
                      {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4},
                      {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4},
                      {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4},
                      { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},
                      {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4},
                      {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},
                      { 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4},
                      {-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4},
                      {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4},
                      {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4},
                      {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4},
                      {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4},
                      {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4},
                      {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4},
                      { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4},
                      { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4},
                      {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2, 11, 2,-3,-4,-3,-2,-4},
                      {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4},
                      { 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4},
                      {-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4},
                      {-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4},
                      { 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4},
                      {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1}};


int** getScoringMatrix(){

  int** matrix;
  int i,j;

  matrix=(int**)space(sizeof(int*)*24);

  for (i=0;i<24;i++){
    matrix[i]=(int*)space(sizeof(int)*24);
  }
  
  for (i=0;i<24;i++){
    for (j=0;j<24;j++){
      matrix[i][j]=BLOSUM62[i][j];
    }
  }

  return matrix;

}

void freeScoringMatrix(int** matrix){

  int i;

  for (i=0;i<24;i++){
    free(matrix[i]);
  }
  free(matrix);
}


void calculateBG(bgModel* model){


  int a1,a2,a3;
  int b1,b2,b3;
  int pepA, pepB;
  double f;
  double prob;
  double score;
  double counts[4];
  int h;
  int i;

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

char* translateSeq(char* seq){

  int i,j;
  int* encodedSeq;
  char* peptideSeq;

  encodedSeq=encodeSeq(seq);

  peptideSeq=(char*)space(sizeof(char)*strlen(seq));
  
  
  j=0;
  for (i=0;i+3<=strlen(seq);i+=3){
    if (encodedSeq[i]>3 || encodedSeq[i+1]>3 || encodedSeq[i+2]>3){
      peptideSeq[j++]='?';
    } else {
      peptideSeq[j++]=decodeAA(transcode[encodedSeq[i]][encodedSeq[i+1]][encodedSeq[i+2]]);
  
    }
  }
  
  peptideSeq[j]='\0';

  return peptideSeq;

}


int* encodeSeq(char* seq){

  int* encoded;
  int i;
  char c;

  encoded=(int*)space(sizeof(int)*(strlen(seq)+1));

  i=0;

  while (c=seq[i]){

    switch(c){

    case 'A': case 'a':
      encoded[i++]=0;
      break;

    case 'C': case 'c':
      encoded[i++]=1;
      break;

    case 'G': case 'g':
      encoded[i++]=2;
      break;

    case 'T': case 't': case 'U': case 'u':
      encoded[i++]=3;
      break;

    case '-': case '.':
      encoded[i++]=4;
      break;

    
    default:
      encoded[i++]=5;
    }
  }

  encoded[i]=-1;
  return encoded;

}


int encodeAA(char aa){

  switch(aa){

  case 'A': case 'a': return 0;
  case 'R': case 'r': return 1;
  case 'N': case 'n': return 2;
  case 'D': case 'd': return 3;
  case 'C': case 'c': return 4;
  case 'Q': case 'q': return 5;
  case 'E': case 'e': return 6;
  case 'G': case 'g': return 7;
  case 'H': case 'h': return 8;
  case 'I': case 'i': return 9;
  case 'L': case 'l': return 10;
  case 'K': case 'k': return 11;
  case 'M': case 'm': return 12;
  case 'F': case 'f': return 13;
  case 'P': case 'p': return 14;
  case 'S': case 's': return 15;
  case 'T': case 't': return 16;
  case 'W': case 'w': return 17;
  case 'Y': case 'y': return 18;
  case 'V': case 'v': return 19;
  case 'B': case 'b': return 20;
  case 'Z': case 'z': return 21;
  case 'X': case 'x': return 99;
  case '*': return -1;
  default: return 99;

  }

}

char decodeAA(int encodedAA){

  switch(encodedAA){

  case 0:   return 'A'; 
  case 1:   return 'R';
  case 2:   return 'N';
  case 3:   return 'D';
  case 4:   return 'C';
  case 5:   return 'Q';
  case 6:   return 'E';
  case 7:   return 'G';
  case 8:   return 'H';
  case 9:   return 'I';
  case 10:  return 'L';
  case 11:  return 'K';
  case 12:  return 'M';
  case 13:  return 'F';
  case 14:  return 'P';
  case 15:  return 'S';
  case 16:  return 'T';
  case 17:  return 'W';
  case 18:  return 'Y';
  case 19:  return 'V';
  case 20:  return 'B';
  case 21:  return 'Z';
  case 99:  return 'X';
  case  -1: return '*';
  default: return '?';
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


int hDist(int a1, int a2, int a3, int b1, int b2, int b3){

  int dist=0;

  if (a1 != b1) dist++;
  if (a2 != b2) dist++;
  if (a3 != b3) dist++;
 
  return dist;

}

/* from GSL */

double gaussian (const double sigma){
  
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
      x = -1 + 2 * (double)rand()/RAND_MAX;
      y = -1 + 2 * (double)rand()/RAND_MAX;

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
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
    
    //    seq=(char*)space((sizeof(char))*(numSites+1));

    dest[i]=createAlnEntry(strdup(src[i]->name),
                           strdup(src[i]->seq),
                           src[i]->start, 
                           src[i]->length,
                           src[i]->fullLength,
                           src[i]->strand);
  }


  dest[i]=NULL;
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

      fprintf(outfile, "\n% 5s% 7s % 6s % 6s % 15s % 9s % 9s % 7s % 9s\n",
              "Frame","Length","From","To","Name","Start","End", "Score","P");
  
      fprintf(outfile, "================================================================================\n");

      i=0;

      while (results[i].score>0){

        direction='+';
        if (results[i].strand==1)  direction='-';

   
        fprintf(outfile, "% 4c%i% 7i % 6i % 6i % 15s % 9i % 9i % 7.3f",
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


      fprintf(outfile,"%s\t%s\t%s\t%i\t%i\t%.2e\t%c\t%c\t%s\n",
              name, "RNAcode","CDS",
              results[i].start+1,results[i].end+1,
              results[i].pvalue,
              direction, '.',"gene_id \"Gene 0\"; transcript_id \"transcript 0\";");
      i++;

      /* GTF, currently only outputs highest scoring hit */
      if (i > 0) break;
    }
  }
}


void stripGaps(struct aln* alignment[]){

  int i,j, k, N;

  char** newSeqs;

  N=0;
  while (alignment[N]!=NULL) N++;
  
  newSeqs=(char**)malloc(sizeof(char*)*N);

  i=0;
  while (alignment[i]!=NULL){
    newSeqs[i]=(char*)malloc(sizeof(char)*(strlen(alignment[i]->seq)+1));
    i++;
  }
    
  j=0;
  k=0;

  while (alignment[0]->seq[j]!='\0'){
    if (alignment[0]->seq[j]!='-'){
      i=0;
      while (alignment[i]!=NULL){
        newSeqs[i][k]=alignment[i]->seq[j];
        i++;
      }
      k++;
    }
    j++;
  }

  i=0;
  while (alignment[i]!=NULL){
    newSeqs[i][k]='\0';
    /* free(alignment[i]->seq); */

    alignment[i]->fullSeq=alignment[i]->seq;

    alignment[i]->seq=newSeqs[i];
    i++;
  }

  free(newSeqs);

}

void freeResults(segmentStats results[]){
  
  int i=0;

  while (results[i].score>0){
    free(results[i].name);
    i++;
  }
  free(results);
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


double* sumOfPairScore(bgModel** modelMatrix, const struct aln *alignment[],int from, int to){

  double* scores;
  int i,j,k,L,N,h;
  char codonA[4];
  char codonB[4];
  int pepA, pepB;
  int gapsA, gapsB;
  double currScore,expectedScore, observedScore;
  int **scoringMatrix;
  int counts;
  int site=0;
  double sum;

  L=strlen(alignment[0]->seq); 
  for (N=0; alignment[N]!=NULL; N++);

  if (to>L-1) to=L-1;

  scores=(double*)malloc(sizeof(double)*L);

  site=0;

  for (k=from;k+2<=to;k+=3){
    
    counts=0;
    sum=0;
    
    for (i=0;i<N;i++){
      
      strncpy(codonA,alignment[i]->seq+k,3);
      codonA[3]='\0';

      gapsA=0;
      if (codonA[0]=='-') gapsA++;
      if (codonA[1]=='-') gapsA++;
      if (codonA[2]=='-') gapsA++;

      for (j=i+1;j<N;j++){
        strncpy(codonB,alignment[j]->seq+k,3);
        codonB[3]='\0';

        gapsB=0;
        if (codonB[0]=='-') gapsB++;
        if (codonB[1]=='-') gapsB++;
        if (codonB[2]=='-') gapsB++;


        h=hDist(ntMap[codonA[0]], ntMap[codonA[1]], ntMap[codonA[2]],
                ntMap[codonB[0]], ntMap[codonB[1]], ntMap[codonB[2]]);

        if (gapsA==0 && gapsB==0 && h>0){

          pepA=transcode[ntMap[codonA[0]]][ntMap[codonA[1]]][ntMap[codonA[2]]];
          pepB=transcode[ntMap[codonB[0]]][ntMap[codonB[1]]][ntMap[codonB[2]]];

          if (pepA==-1) pepA=23;
          if (pepB==-1) pepB=23;

          expectedScore=modelMatrix[i][j].scores[h];
          observedScore=(double)modelMatrix[i][j].matrix[pepA][pepB];

          /*
          printf("%i %s (%c) : %i %s (%c); observed: %.2f, expected: %.2f\n",
                 i, codonA, decodeAA(pepA), j, codonB, decodeAA(pepB),observedScore, expectedScore);
          */
          
          sum+=observedScore-expectedScore;
          
          counts++;
        
        }
      }
    }
    //printf("%.2f\n\n",sum/(double)counts);

    if (counts>0){
      scores[site++]=sum/(double)counts;
    } else {
      scores[site++]=0.0;
    }
  }

  return scores;
}

double* getCumSum(double* scores, int N){
  
  int i;
  double sum;
  double* cum;

  cum=(double*)malloc(sizeof(double)*N);

  sum=0.0;
  
  for (i=0; i<N; i++){
    sum+=scores[i];
    if (sum<0){
      sum=0.0;
    }
    cum[i]=sum;
  }
  return cum;
}




bgModel** getModelMatrix(TTree* tree, struct aln *alignment[], double kappa){

  int i,j,N;
  int **scoringMatrix;
  double** distanceMatrix;
  double freqsMono[4];
  bgModel** modelMatrix;

  for (N=0; alignment[N]!=NULL; N++);   
  
  distanceMatrix=getDistanceMatrix(tree,(struct aln**)alignment);
  
  countFreqsMono((const struct aln**)alignment, (double *) freqsMono);

  modelMatrix=(bgModel**)malloc(sizeof(bgModel*)*N);

  for (i=0;i<N;i++){
    modelMatrix[i]=(bgModel*)malloc(sizeof(bgModel)*N);
  }

  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      scoringMatrix=getScoringMatrix();
      modelMatrix[i][j].dist=distanceMatrix[i][j];
      modelMatrix[i][j].kappa=kappa;
      modelMatrix[i][j].freqs[0]=freqsMono[0];
      modelMatrix[i][j].freqs[1]=freqsMono[1];
      modelMatrix[i][j].freqs[2]=freqsMono[2];
      modelMatrix[i][j].freqs[3]=freqsMono[3];
      modelMatrix[i][j].matrix=scoringMatrix;
      calculateBG(modelMatrix[i]+j);
    }
  }

  for (i=0;i<N;i++){
    free(distanceMatrix[i]);
  }

  free(distanceMatrix);

  return modelMatrix;
  
}

void freeModelMatrix(bgModel** modelMatrix, int N){

  int i, j;

  for (i=0;i<N;i++){
    for (j=0;j<N;j++){
      freeScoringMatrix(modelMatrix[i][j].matrix);
    }
  }

  for (i=0;i<N;i++){
    free(modelMatrix[i]);
  }
  free(modelMatrix);
}

void getExtremeValuePars(TTree* tree, bgModel** modelMatrix, const struct aln *alignment[], 
                         int sampleN, int sampleMode, double* parMu, double* parLambda){

  int tmpCounter, L, i, j;
  double* freqsMono;
  double kappa;
  double sum, maxSum;
  double* scores;
  double* maxScores;
  double mu, sigma;
  struct aln *sampledAln[MAX_NUM_NAMES];

  L=strlen(alignment[0]->seq);

  freqsMono=modelMatrix[0][1].freqs;
  kappa=modelMatrix[0][1].kappa;

  maxScores=(double*)malloc(sizeof(double)*sampleN);

  if (sampleMode==1){

    tmpCounter=0;

    for (i=0;i<sampleN;i++){
  
      simulateTree(tree,freqsMono,kappa,L);
      
      tree2aln(tree,sampledAln);

      sortAln(alignment, sampledAln);
      
      reintroduceGaps(alignment, sampledAln);
      
      scores=sumOfPairScore(modelMatrix,(const struct aln**)sampledAln,0, L-1);

      freeAln((struct aln**)sampledAln);

      sum=0.0;
      maxSum=0;

      for (j=0; j<L/3; j++){
        sum+=scores[j];
        if (sum<0){
          sum=0.0;
        }
        maxSum=(sum>maxSum)?sum:maxSum;
      }
      
      maxScores[i]=maxSum;
      free(scores);
    }

  } else {

    simulateTree(tree,freqsMono,kappa,sampleN*3);

    tree2aln(tree,sampledAln);

    sortAln(alignment, sampledAln);

    randomlyPlaceGaps(alignment, sampledAln);

    //    printAlnClustal(stdout, (const struct aln**)sampledAln); 

    scores=sumOfPairScore(modelMatrix,(const struct aln**)sampledAln,0, sampleN*3);

    freeAln((struct aln**)sampledAln);

    mu=avg(scores,sampleN);
    sigma=stddev(scores,sampleN);

    free(scores);

    for (i=0;i<sampleN;i++){
      sum=0.0;
      maxSum=0;
      for (j=0; j<L/3; j++){
        sum+=gaussian(sigma)+mu;
        if (sum<0){
          sum=0.0;
        }
        maxSum=(sum>maxSum)?sum:maxSum;
      }
      maxScores[i]=maxSum;
    }
    
    for (i=0;i<sampleN;i++){
      //fprintf(stderr, "%.2f\n",maxScores[i]);
    }
  }

  EVDMaxLikelyFit(maxScores, NULL, sampleN, parMu, parLambda);

  free(maxScores);

}


segmentStats* getHSS(bgModel** modelMatrix, const struct aln** inputAln, double parMu, double parLambda,double cutoff){

  segmentStats* results;
  int strand, frame, sites,L;
  double* scores;
  double* cumSum;
  double* maxScores;
  int segmentStart, segmentEnd, currPos, hssCount;
  double currMax;
  double pvalue;
  int i;

  struct aln *inputAlnRev[MAX_NUM_NAMES];
  struct aln **currAln;
  
  L=strlen(inputAln[0]->seq);

  copyAln((struct aln**)inputAln,(struct aln**)inputAlnRev);

  revAln((struct aln**)inputAlnRev);
  
  results=NULL;


  hssCount=0;
  
  for (strand=0;strand<=1;strand++){
    
    if (strand==0){
      currAln=(struct aln**)inputAln;
    } else {
      currAln=(struct aln**)inputAlnRev;
    }

    for (frame=0;frame<3;frame++){

      sites=((int)(L-frame)/3);

      //printf("Frame %i, strand %i, sites: %i:\n",frame, strand, sites);

      scores=sumOfPairScore(modelMatrix,(const struct aln**)currAln,frame, L-1);

      cumSum=getCumSum(scores,sites);

      
      /*
        for (i=0;i<sites;i++){
        printf("%i %.2f %.2f\n",i, scores[i],cumSum[i]);
        }
      */

      /*

      for (i=0;i<sites;i++){
      printf("%i\t",i);
      }

      printf("\n");

      for (i=0;i<sites;i++){
         printf("%c%c%c\t",currAln[0]->seq[i*3+frame],currAln[0]->seq[i*3+1+frame],currAln[0]->seq[i*3+2+frame]);
      }
     
      printf("\n");

      for (i=0;i<sites;i++){
         printf("%.2f\t",cumSum[i]);
      }
      printf("\n");
      */

      segmentStart=0;
      segmentEnd=1;
      currMax=cumSum[0];
      for (currPos=0;currPos<sites;currPos++){

        if (cumSum[currPos]>=currMax){
          segmentEnd=currPos;
          currMax=cumSum[currPos];
        }

        if ((cumSum[currPos]<0.0001) || (currPos==sites-1) ){
          if (segmentEnd-segmentStart>2){

            if (scores[segmentStart]<0){
              segmentStart++;
            }
            
            pvalue=1-exp((-1)*exp((-1)*parLambda*(currMax-parMu)));

            if (pvalue<cutoff){

              /* re-allocate for each new result, leave room for last entry */
              results=(segmentStats*)realloc(results,sizeof(segmentStats)*(hssCount+2));

              if (results==NULL){
                exit(1);
              }

              /* chromosome name is also stored in results, note that we
                 use statically allocated memory there */
            
              //strncpy((char*)results[hssCount].name,inputAln[0]->name,256);
              //results[hssCount].name[strlen(inputAln[0]->seq)]='\0';

              results[hssCount].name=strdup(inputAln[0]->name);
              results[hssCount].strand=strand;
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

                        
              //printf("from: %i to %i\n", segmentStart, segmentEnd);
              hssCount++;
            }

          }
          segmentStart=currPos;
          segmentEnd=currPos;
          currMax=0.0;
          
        }
      }
      free(cumSum);
      free(scores);
    }
  }

  if (hssCount==0){
    results=(segmentStats*)malloc(sizeof(segmentStats));
  }

  results[hssCount].score=-1; /* mark end of list */

  freeAln((struct aln**)inputAlnRev);

  return results;


}


segmentStats* getHSSnew(bgModel** modelMatrix, const struct aln** inputAln, double parMu, double parLambda,double cutoff){

  segmentStats* results;
  int strand, frame, sites,L;
  double* scores;
  double* maxScores;
  int segmentStart, segmentEnd, currPos, hssCount;
  double currMax, cumSum;
  double pvalue;
  int i,j;

  struct aln *inputAlnRev[MAX_NUM_NAMES];
  struct aln **currAln;
  
  L=strlen(inputAln[0]->seq);

  copyAln((struct aln**)inputAln,(struct aln**)inputAlnRev);

  revAln((struct aln**)inputAlnRev);
  
  results=NULL;

  hssCount=0;
  
  for (strand=0;strand<=0;strand++){
    
    if (strand==0){
      currAln=(struct aln**)inputAln;
    } else {
      currAln=(struct aln**)inputAlnRev;
    }

    for (frame=0;frame<1;frame++){

      sites=((int)(L-frame)/3);

      //printf("Frame %i, strand %i, sites: %i:\n",frame, strand, sites);

      scores=sumOfPairScore(modelMatrix,(const struct aln**)currAln,frame, L-1);

      cumSum=0.0;
      segmentStart=0;
      segmentEnd=1;
      currMax=cumSum;
      for (currPos=0;currPos<sites;currPos++){

        cumSum+=scores[currPos];

        if (cumSum<0.0){
          cumSum=0.0;
        }

        if (cumSum>=currMax){
          segmentEnd=currPos;
          currMax=cumSum;
        }

        printf("%i-%i %.2f %.2f %.2f\n",segmentStart, segmentEnd,scores[currPos],cumSum, currMax);

        for (j=0;inputAln[j]!=NULL;j++){
          for (i=segmentStart*3;i<=segmentEnd*3;i+=3){
            printf("%c%c%c ",inputAln[j]->seq[i],inputAln[j]->seq[i+1],inputAln[j]->seq[i+2]);
          }
          printf("\n");
        }

        printf("======================\n");

        for (j=0;inputAln[j]!=NULL;j++){
          for (i=segmentStart*3;i<=segmentEnd*3;i+=3){
            printf("%c%c%c ",inputAln[j]->fullSeq[i],inputAln[j]->fullSeq[i+1],inputAln[j]->fullSeq[i+2]);
          }
          printf("\n");
        }


        

        printf("\n");

        if ((cumSum<0.0001) || (currPos==sites-1) ){
          if (segmentEnd-segmentStart>2){

            if (scores[segmentStart]<0){
              segmentStart++;
            }
            
            pvalue=1-exp((-1)*exp((-1)*parLambda*(currMax-parMu)));

            if (pvalue<cutoff){

              /* re-allocate for each new result, leave room for last entry */
              results=(segmentStats*)realloc(results,sizeof(segmentStats)*(hssCount+2));

              if (results==NULL){
                exit(1);
              }

              /* chromosome name is also stored in results, note that we
                 use statically allocated memory there */
            
              //strncpy((char*)results[hssCount].name,inputAln[0]->name,256);
              //results[hssCount].name[strlen(inputAln[0]->seq)]='\0';

              results[hssCount].name=strdup(inputAln[0]->name);
              results[hssCount].strand=strand;
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

                        
              //printf("from: %i to %i\n", segmentStart, segmentEnd);
              hssCount++;
            }

          }
          segmentStart=currPos;
          segmentEnd=currPos;
          currMax=0.0;
          
        }
      }
      //free(cumSum);
      free(scores);
    }
  }

  if (hssCount==0){
    results=(segmentStats*)malloc(sizeof(segmentStats));
  }

  results[hssCount].score=-1; /* mark end of list */

  freeAln((struct aln**)inputAlnRev);

  return results;


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


void randomlyPlaceGaps(const struct aln* origAln[], struct aln* sampledAln[]){

  int i,j,k,site,randomSite;
  int sampledL, origL;

  sampledL=strlen(sampledAln[0]->seq);
  origL=strlen(origAln[0]->seq);

  for (site=0;site<sampledL-2;site+=3){

    /* choose random site in  [0,L-2)   */
    randomSite=(int)((double)rand()/((double)RAND_MAX+(double)1)*(origL-2));

    for (i=0;origAln[i]!=NULL;i++){
      for (k=0;k<2;k++){
        if (origAln[i]->seq[randomSite+k]=='-'){
          sampledAln[i]->seq[site+k]='-';
        }
      }
    }
  }
}

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "score.h"
#include "code.h"
#include "misc.h"
#include "postscript.h"

float fontWidth, fontHeight, imageHeight, imageWidth,tmpColumns;
int length, maxName, maxNum, currPos, columnWidth;
float lineStep,blockStep,consStep,ssStep,rulerStep,nameStep,numberStep, maxConsBar;
float startY,namesX,seqsX, currY;

FILE *outfile;
float**** Sk_native;   
float**** Sk_native_rev;  

void setParameters() {

  columnWidth=60;            /* Display long alignments in blocks of this size */
  fontWidth=6;               /* Font metrics */
  fontHeight=6.5;
  lineStep=fontHeight+2;     /* distance between lines */
  blockStep=3.5*fontHeight;  /* distance between blocks */
  consStep=fontHeight*0.5;   /* distance between alignment and conservation curve */
  ssStep=12;                 /* distance between secondary structure line and sequences */
  rulerStep=2;               /* distance between sequences and ruler */
  nameStep=3*fontWidth;	     /* distance between names and sequences */
  numberStep=fontWidth;      /* distance between sequeces and numbers */
  maxConsBar=2.5*fontHeight; /* Height of conservation curve */
  startY=2;		             /* "y origin" */
  namesX=fontWidth;	         /* "x origin" */

}

int colorAln(const char *filename, const struct aln *alignment[],segmentStats region) {

  int N,i,j,k,x,tmp; 
  char *tmpBuffer,*ruler, *cons; 
  int *syn;
  int *non_syn;
  char c;
  float score,barHeight,xx,yy;
  int match,block;
  int hssCount, start, end;
  backtrackData* bt;
  char label[1024];
  char pString[1024];


  char *codeplot_header =
    "%%!PS-Adobe-3.0 EPSF-3.0\n"
    "%%%%BoundingBox: %i %i %i %i\n"
    "%%%%EndComments\n"
    "%%Created by RNAcode; visit wash.github.com/rnacode\n"
    "%% draws box in color given by hue and saturation\n"
    "/box { %% x1 y1 x2 y2 hue saturation\n"
    "  gsave\n"
    "  dup 0.3 mul 1 exch sub sethsbcolor\n"
    "  exch 3 index sub exch 2 index sub rectfill\n"
    "  grestore\n"
    "} def\n"
    "%% draws a box in current color\n"
    "/box2 { %% x1 y1 x2 y2\n"
    "  exch 3 index sub exch 2 index sub rectfill\n"
    "} def\n"
    "/string { %% (Text) x y\n"
    " 6 add\n"
    " moveto\n"
    "  show\n"
    "} def\n"
    "0 %i translate\n"
    "1 -1 scale\n"
    "/Courier findfont\n"
    "[10 0 0 -10 0 0] makefont setfont\n";
  

  int** scoringMatrix;

  int countMatrix[6][6][6];
  int countSyn;
  int countNonSyn;
  int pepA, pepB;
  char codonA[4];
  char codonB[4];

  struct aln *alignmentRev[MAX_NUM_NAMES];
  struct aln *currAln[MAX_NUM_NAMES];
    
  setParameters();
  
  outfile = fopen(filename, "w");

  if (outfile == NULL) {
    fprintf(stderr, "WARNING: Cannot write to file %s - not doing alignment plot.\n",
            filename);
    return 0;
  }
  
  /* Number of columns of the alignment */
  length=strlen(alignment[0]->seq);
  
  for (N=0; alignment[N]!=NULL; N++);   

  /* Allocate memory for various strings, length*2 is (more than)
     enough for all of them */
  tmpBuffer = (char *) malloc((unsigned) length*2);
  ruler=(char *) malloc((unsigned) length*2);
  //syn=(int *) malloc(sizeof(int)*length*2);
  //non_syn=(int *)malloc(sizeof(int)*length*2);

  /* Get length of longest name and count sequences in alignment*/

  for (i=maxName=N=0; alignment[i] != NULL; i++) {
    N++;
    tmp=strlen(alignment[i]->name);
    if (tmp>maxName)  maxName=tmp;
  }
  
  /* x-coord. where sequences start */
  seqsX=namesX+maxName*fontWidth+nameStep;

  /* calculate number of digits of the alignment length */
  //snprintf(tmpBuffer,length, "%i",length);
  //sprintf(tmpBuffer, "%i",length);
  //printf("tmpBuffer: %s\n",tmpBuffer);
  //maxNum=strlen(tmpBuffer);

  maxNum=10;

  /* Calculate bounding box */
  tmpColumns=columnWidth;
  if (length<columnWidth){
    columnWidth=length;
    tmpColumns=length;
  }
  
  imageWidth=ceil(namesX+(maxName+tmpColumns+maxNum)*fontWidth+2*nameStep+fontWidth+numberStep);
  imageHeight=startY+ceil((float)length/columnWidth)*((N+2)*lineStep+blockStep+consStep+ssStep+rulerStep);

  /* Write postscript header including correct bounding box */
  fprintf(outfile,codeplot_header,0,0,(int)imageWidth,(int)imageHeight,(int)imageHeight);

  /* Create ruler */
  i=0;
  /* Init all with dots */
  for (i=0;i<(length);i++){
    ruler[i]='.';
  }
  i=0;
  for (i=0;i<length;i++){
    /* Write number every 10th position, leave out block breaks */
    if ((i+1)%10==0 && (i+1)%columnWidth!=0){
      snprintf(tmpBuffer,length,"%i",i+1);
      memcpy(ruler+i, tmpBuffer, strlen(tmpBuffer));
    }
  }
  
  ruler[length]='\0';
  
  currY=startY;
  currPos=0;

  cons =  consensus(alignment);

  while (currPos<length) {

    fprintf(outfile,"0 setgray\n");
      
    currY+=ssStep+lineStep;
    
    // Display names, sequences and numbers

    for (i=0; i<N; i++) {
      
      strncpy(tmpBuffer,alignment[i]->seq+currPos,columnWidth);
      tmpBuffer[columnWidth]='\0';
      match=0;
           
      for (j=0;j<(currPos+strlen(tmpBuffer));j++){
        if (alignment[i]->seq[j] != '-') match++;
      }

      if (region.strand == '+'){
        match+=alignment[i]->start;
      } else {
        match=getSeqLength(alignment[i]->seq)-match+1;
      }
      
      fprintf(outfile, "(%s) %.1f %.1f string\n", alignment[i]->name,namesX,currY);

      //fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer,seqsX,currY);

      //for (j = 0; j < columnWidth; j++){
      //fprintf(outfile, "(%c) %.1f %.1f string\n", tmpBuffer[j],seqsX+j*fontWidth,currY);
      //}
      
      fprintf(outfile, "(%i) %.1f %.1f string\n", match,seqsX+fontWidth*(strlen(tmpBuffer))+numberStep,currY);
      currY+=lineStep;
    }
    currY+=rulerStep;
    strncpy(tmpBuffer,ruler+currPos,columnWidth);
    tmpBuffer[columnWidth]='\0';
    
    fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer,seqsX,currY);
    fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer,seqsX,currY);
    
    currY+=lineStep;
    currY+=consStep;
    
    //Display conservation bar
       
    fprintf(outfile,"0.6 setgray\n");
    for (i=currPos;(i<currPos+columnWidth && i<length);i++){
      match=0;

      for (j=0;j<N;j++){
        if (cons[i] == alignment[j]->seq[i]) match++;
        if (cons[i]=='U' && alignment[j]->seq[i]=='T') match++;
        if (cons[i]=='T' && alignment[j]->seq[i]=='U') match++;
      }
      
      score=(float)(match-1)/(N-1);
      
      if (cons[i] == '-' ||
          cons[i] == '_' ||
          cons[i] == '.'){
        score=0;
      }
      
      barHeight=maxConsBar*score;
      if (barHeight==0){
        barHeight=1;
      }
      
      xx=seqsX+(i-(columnWidth*currPos/columnWidth))*fontWidth;
      
      fprintf(outfile,"%.1f %.1f %.1f %.1f box2\n",
              xx,
              currY+maxConsBar-barHeight,
              xx+fontWidth,
              currY+maxConsBar);
    }
    
    currY+=blockStep;
    currPos+=columnWidth;
  }
  

  fprintf(outfile,"0.0 setgray\n");

  copyAln((struct aln**)alignment,(struct aln**)alignmentRev);
  revAln((struct aln**)alignmentRev);

  if (region.strand == '+'){
    copyAln((struct aln**)alignment,(struct aln**)currAln);
  } else {
    copyAln((struct aln**)alignmentRev,(struct aln**)currAln);
  }
    
  for (x=0;x<3;x++){

    // left
    if (x == 0){
      start = extendRegion((const struct aln**)currAln, region.start, 0);
      end = region.start-1;
      //printf("====> Extend left %i to %i", start, end);
      strcpy(label, "");
    }

    // top scoring hit
    if (x==1){
      start = region.start;
      end  = region.end;

      if (region.pvalue < 0.001){
        if (region.pvalue < 10e-16){
          sprintf(pString, "<1e-16\n");
        } else {
          sprintf(pString, "%9.1e\n",region.pvalue);
        }
      } else {
        sprintf(pString, "%9.3f\n",region.pvalue);
      }
      
      sprintf(label, "Frame %c%i p =%s", region.strand,region.frame+1,pString );
    }

    // right
    if (x==2){
      start = region.end+1;
      end = extendRegion((const struct aln**)currAln, region.end, 1);

      if (start >= end){
        break;
      }
      //printf("====> Extend right %i to %i", start, end);
      strcpy(label, "");
    }

    if (region.strand == '+'){
      bt = backtrack(start, end, Sk_native, (const struct aln**)currAln);
    } else {
      bt = backtrack(start, end, Sk_native_rev, (const struct aln**)currAln);
    }

    colorHSS(filename, (const struct aln**)currAln, bt, label, start, end);

    for (k=1; k<N; k++){
      free(bt[k].states);
      free(bt[k].z);
      free(bt[k].transitions);
      free(bt[k].scores);
    }

    free(bt);
    
    
  }

  free(cons);
  
  fprintf(outfile,"showpage\n");
  
  fclose(outfile);

  free(tmpBuffer);   
  free(ruler);
  
  freeAln((struct aln**)alignmentRev);
  freeAln((struct aln**)currAln);
 
  return 0;

}

void colorHSS(const char *filename, const struct aln *alignment[], backtrackData *bt, const char* label, int b, int i) {
  
  int N, k, x, z, L, l, colsN, currSum, currN;
  char *block_0, *block_k, *seq_0, *seq_k;
  char **blocks_0, **blocks_k;
  char *translated;
  int *map_0, *map_k;
  char codonA[4]="XXX";
  char codonB[4]="XXX";
  int ii, jj, kk;
  int syn_count, non_syn_count, currScore; 
  int alreadySeen[4][4][4];
  int** scoringMatrix;
  int pepA, pepB;
  int *syn;
  int *non_syn;
  float *scores;

  char * colorMatrix[2][6] = {
    //{"0.0 0.15", "0.0 0.2",  "0.0 0.4", "0.0 0.6", "0.0 0.8",  "0.0 1"},  /* red    */
    {"0.0 0.0", "0.0 0.2",  "0.0 0.4", "0.0 0.6", "0.0 0.8",  "0.0 1"},  /* red    */
    {"0.32 0.1", "0.32 0.2",  "0.32 0.4", "0.32 0.6", "0.32 0.8",  "0.32 1"},  /* green    */
  };

  scoringMatrix=getScoringMatrix();

  seq_0=alignment[0]->seq;
  L=getSeqLength(seq_0);

  syn=(int *) malloc(sizeof(int)*L*2);
  non_syn=(int *)malloc(sizeof(int)*L*2);
  scores=(float *) malloc(sizeof(float)*L*2);

  colsN=strlen(seq_0);
  for (N=0; alignment[N]!=NULL; N++);

  block_0 = (char*) malloc(sizeof(char)*(colsN+1));
  block_k = (char*) malloc(sizeof(char)*(colsN+1));
  blocks_0 = (char**) malloc(sizeof(char*)*(colsN+1));
  blocks_k = (char**) malloc(sizeof(char*)*(colsN+1));
  
  map_0=(int*)malloc(sizeof(int)*(colsN+1));
  map_k=(int*)malloc(sizeof(int)*(colsN+1));

  for (l=1;l<=L;l++){
    map_0[l]=pos2col(seq_0,l);
  }

  for (x=b+2;x<i+3;x+=3){
    
    syn_count=0;
    non_syn_count=0;
    currScore=0;

    for (ii=0;ii<4;ii++)
      for (jj=0;jj<4;jj++)
        for (kk=0;kk<4;kk++)
          alreadySeen[ii][jj][kk]=0;
  
    currSum=0;
    currN=0;
  
    for (k=1;k<N;k++){

      seq_k=alignment[k]->seq;
    
      for (l=1;l<=L;l++){
        map_k[l]=pos2col(seq_k,l);
      }
      
      getBlock(x, seq_0, seq_k, map_0, map_k, block_0, block_k, &z );

      ii=jj=0;
  
      while (block_0[ii] != '\0') {
        if (block_0[ii] != '-') {
          codonA[jj]=block_0[ii];
          codonB[jj]=block_k[ii];
          jj++;
        }
        ii++;
      }

      pepA=transcode[ntMap[codonA[0]]][ntMap[codonA[1]]][ntMap[codonA[2]]];

      if (codonB[0] == '-' || codonB[1] == '-' || codonB[2] == '-'){
        continue;
      }

      if (alreadySeen[ntMap[codonB[0]]][ntMap[codonB[1]]][ntMap[codonB[2]]]){
        continue;
      } else {
        alreadySeen[ntMap[codonB[0]]][ntMap[codonB[1]]][ntMap[codonB[2]]]=1;
      }

      pepB=transcode[ntMap[codonB[0]]][ntMap[codonB[1]]][ntMap[codonB[2]]];
      
      if (pepA != -1 && pepB != -1){

        currScore=scoringMatrix[pepA][pepB];

        if (strcmp(codonA, codonB) !=0){
          currSum+=currScore;
          currN++;
          //printf("  %s %s %i\n", codonA, codonB, currScore);
        
        
          if (pepA == pepB){
            if (strcmp(codonA, codonB)!=0){
              syn_count++;
            }
          } else {
            if (currScore<0){
              non_syn_count++;
            }
          }
        }
      } else {
        non_syn_count++;
      }
    }

    //printf("%i: %i vs %i\n", x, syn_count, non_syn_count);

    syn[x] = syn_count;
    non_syn[x] = non_syn_count;

    if (currN > 0){
      scores[x] = (float)currSum/currN;
    } else {
      scores[x] = 0.0;
    }
   
  }

  for (x=b+2;x<i+3;x+=3){

    char *syn_color;
    char *non_syn_color;
    int syn_count=0;
    int non_syn_count=0;
    int currScore=0;
    int currCol;
    int block;
    float xx, yy;
   
    syn_count=syn[x];
    non_syn_count=non_syn[x];

    if (non_syn_count>=5) {
      non_syn_count=5;
    }

    if (syn_count>=5) {
      syn_count=5;
    }

    for (k=0;k<N;k++){

      seq_k=alignment[k]->seq;
    
      for (l=1;l<=L;l++){
        map_k[l]=pos2col(seq_k,l);
      }
      
      getBlock(x, seq_0, seq_k, map_0, map_k, block_0, block_k, &z );

      syn_color=colorMatrix[1][syn_count];
      non_syn_color=colorMatrix[0][non_syn_count];

      ii=jj=0;
  
      while (block_0[ii] != '\0') {
        if (block_0[ii] != '-') {
          codonA[jj]=block_0[ii];
          codonB[jj]=block_k[ii];
          jj++;
        }
        ii++;
      }

      pepA=transcode[ntMap[codonA[0]]][ntMap[codonA[1]]][ntMap[codonA[2]]];

      if (codonB[0] == '-' || codonB[1] == '-' || codonB[2] == '-'){
        currScore=-1;
        pepB=-99; /* contains gap */
      } else {
        pepB=transcode[ntMap[codonB[0]]][ntMap[codonB[1]]][ntMap[codonB[2]]];
        if (pepA != -1 && pepB != -1){
          currScore=scoringMatrix[pepA][pepB];
        } else {
          currScore=-1;
          pepB=+99; /* Stop */;
        }
      }

      ii=jj=0;

      while (block_0[ii] != '\0') {
        
        currCol=pos2col(alignment[0]->seq,x)-strlen(block_0)+ii;
        block=ceil((float)(currCol+1)/columnWidth);
        xx=seqsX+(currCol-(block-1)*columnWidth)*fontWidth;

        yy=startY+(block-1)*(lineStep*(N+2)+blockStep+consStep+rulerStep)+ssStep*(block)+(k+1)*lineStep;

        if (k==0 && strcmp(label,"") !=0 && x==b+2 && ii == 0){
          fprintf(outfile,"0.15 0.5 0.6 sethsbcolor\n");
          fprintf(outfile,"/Helvetica findfont\n");
          fprintf(outfile,"[8 0 0 -8 0 0] makefont setfont\n");
          fprintf(outfile, "(%s) %.1f %.1f string\n", label,xx,yy-2*lineStep);
          fprintf(outfile,"0.0 setgray\n");
        }


        /* Print blocks of translation line */
        if (k==0 && strcmp(label,"") !=0){
          float offsetLeft=0.0;
          float offsetRight=0.0;

          if (ii==0){
            offsetLeft=0.5;
          }

          if (ii==strlen(block_0)-1){
            offsetRight=0.5;
          }
          
          fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                  xx+offsetLeft,yy-1,xx+fontWidth-offsetRight,yy-lineStep-1,"0.15 0.5");
        }

        /* Print translation in the mid of the codons*/
        if (k==0 && ii==ceil((float)strlen(block_0)/2.0)-1){
          fprintf(outfile,"/Courier findfont\n");
          fprintf(outfile,"[10 0 0 -10 0 0] makefont setfont\n");
          translated =  translateSeq(codonA);
          fprintf(outfile, "(%s) %.1f %.1f string\n",translated,xx,yy-lineStep);
          free(translated);
        }

        /* Print reference sequence */
        if (k == 0){
          fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                  xx,yy-1,xx+fontWidth,yy+fontHeight+1,"0.0 0.0");
          fprintf(outfile,"/Courier-Bold findfont\n");
          fprintf(outfile,"[10 0 0 -10 0 0] makefont setfont\n");
          fprintf(outfile, "(%c) %.1f %.1f string\n", block_k[ii],xx,yy);
        }

                
        /* In frame states */
        if (k > 0 && bt[k].states[x]==0 && bt[k].transitions[x] == 0) {

          if (currScore>=0){

            /* Synonymous mutations */
            if (strcmp(codonA,codonB)!=0){
              fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                      xx,yy-1,xx+fontWidth,yy+fontHeight+1,syn_color);
              /* Conservative mutations*/
            } else {
              fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                      xx,yy-1,xx+fontWidth,yy+fontHeight+1,"0.0 0.0");
            }

            //fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
            //        xx,yy-1,xx+fontWidth,yy+fontHeight+1,syn_color);


          } else {
            //printf("==> %i %i\n", currScore, pepB);
            // Stop
            if (pepB==+99){
              fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                      xx,yy-1,xx+fontWidth,yy+fontHeight+1,"0.6 1.0");
              //fprintf(outfile,"1 setgray\n");
            }

            // contains gap
            if (pepB == -99){
              fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                      xx,yy-1,xx+fontWidth,yy+fontHeight+1,"0.0 0.0");
            }
     
            if (pepB != 99 && pepB != -99){
              fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                      xx,yy-1,xx+fontWidth,yy+fontHeight+1,non_syn_color);
            }
          }
          
          if (pepA == pepB && strcmp(codonA,codonB)!=0){
            fprintf(outfile,"/Courier-Bold findfont\n");
            fprintf(outfile,"[10 0 0 -10 0 0] makefont setfont\n");
          } else {
            fprintf(outfile,"/Courier findfont\n");
            fprintf(outfile,"[10 0 0 -10 0 0] makefont setfont\n");
          }
          
        }

        if (k>0 && bt[k].transitions[x] == 2){
          //fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
          //        xx,yy-1,xx+fontWidth,yy+fontHeight+1,"0.5 0.1");
          fprintf(outfile,"/Courier-Bold findfont\n");
          fprintf(outfile,"[10 0 0 -10 0 0] makefont setfont\n");
          fprintf(outfile,"0.2 setgray\n");
          fprintf(outfile, "%.1f %.1f %.1f %.1f box2\n",
                  xx,yy-1,xx+fontWidth,yy+fontHeight+1);
          fprintf(outfile,"0.8 setgray\n");
        }

        if (k>0 && (bt[k].transitions[x] == 1 || ( bt[k].transitions[x] == 0 && bt[k].states[x] != 0))){
          fprintf(outfile,"/Courier findfont\n");
          fprintf(outfile,"[10 0 0 -10 0 0] makefont setfont\n");
          fprintf(outfile,"0.8 setgray\n");
          fprintf(outfile, "%.1f %.1f %.1f %.1f box2\n",
                  xx,yy-1,xx+fontWidth,yy+fontHeight+1);
          fprintf(outfile,"0 setgray\n");
        }
        
        fprintf(outfile, "(%c) %.1f %.1f string\n", block_k[ii],xx,yy);

        fprintf(outfile,"0 setgray\n");

        ii++;
      }
    }
  }

  free(blocks_0);
  free(blocks_k);
  free(block_0);
  free(block_k);
  free(map_0);
  free(map_k);
  free(syn);
  free(non_syn);
  free(scores);
  freeScoringMatrix(scoringMatrix);


}



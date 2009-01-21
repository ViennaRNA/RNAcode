#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "code.h"
#include "postscript.h"

extern char* codeplot_header;

int colorAln(const char *filename, const struct aln *alignment[],segmentStats region) {

  int N,i,j,k,x,y,tmp,columnWidth; 
  int cdsStart, cdsEnd;
  char *tmpBuffer,*ssEscaped,*ruler, *cds, *cons; 
  int *syn;
  int *non_syn;
  char c;
  float fontWidth, fontHeight, imageHeight, imageWidth,tmpColumns;
  int length, cdsLength, maxName, maxNum, currPos;
  float lineStep,blockStep,consStep,ssStep,rulerStep,nameStep,numberStep;
  float maxConsBar,startY,namesX,seqsX, currY;
  float score,barHeight,xx,yy;
  int match,block;
  FILE *outfile;

  char * colorMatrix[2][6] = {
    {"0.0 0.1", "0.0 0.2",  "0.0 0.4", "0.0 0.6", "0.0 0.8",  "0.0 1"},  /* red    */
    {"0.32 0.0", "0.32 0.2",  "0.32 0.4", "0.32 0.6", "0.32 0.8",  "0.32 1"},  /* green    */
  };

  char *codeplot_header =
    "%%!PS-Adobe-3.0 EPSF-3.0\n"
    "%%%%BoundingBox: %i %i %i %i\n"
    "%%%%EndComments\n"
    "%% draws Vienna RNA like colored boxes\n"
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

  printAlnClustal(stdout, alignment);

  scoringMatrix=getScoringMatrix();
  
  outfile = fopen(filename, "w");

  if (outfile == NULL) {
    fprintf(stderr, "can't open file %s - not doing alignment plot\n",
            filename);
    return 0;
  }
  
  columnWidth=60;            /* Display long alignments in blocks of this size */
  fontWidth=6;               /* Font metrics */
  fontHeight=6.5;
  lineStep=fontHeight+2;     /* distance between lines */
  blockStep=3.5*fontHeight;  /* distance between blocks */
  consStep=fontHeight*0.5;   /* distance between alignment and conservation curve */
  ssStep=2;                  /* distance between secondary structure line and sequences */
  rulerStep=2;               /* distance between sequences and ruler */
  nameStep=3*fontWidth;	     /* distance between names and sequences */
  numberStep=fontWidth;      /* distance between sequeces and numbers */
  maxConsBar=2.5*fontHeight; /* Height of conservation curve */
  startY=2;		               /* "y origin" */
  namesX=fontWidth;	         /* "x origin" */
  
  /* Number of columns of the alignment */
  length=strlen(alignment[0]->seq);
  
  for (N=0; alignment[N]!=NULL; N++);   

  /* Allocate memory for various strings, length*2 is (more than)
     enough for all of them */
  tmpBuffer = (char *) space((unsigned) length*2);
  ssEscaped=(char *) space((unsigned) length*2);
  ruler=(char *) space((unsigned) length*2);
  cds=(char *) space((unsigned) length*2);
  syn=(int *) space(sizeof(int)*length*2);
  non_syn=(int *)space(sizeof(int)*length*2);

  /* Get length of longest name and count sequences in alignment*/

  for (i=maxName=N=0; alignment[i] != NULL; i++) {
    N++;
    tmp=strlen(alignment[i]->name);
    if (tmp>maxName)  maxName=tmp;
  }
  
  /* x-coord. where sequences start */
  seqsX=namesX+maxName*fontWidth+nameStep;

  /* calculate number of digits of the alignment length */
  snprintf(tmpBuffer,length, "%i",length);
  maxNum=strlen(tmpBuffer);
 

  /* Calculate bounding box */
  tmpColumns=columnWidth;
  if (length<columnWidth){
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
      strncpy(ruler+i,tmpBuffer,strlen(tmpBuffer));
    }
  }
  ruler[length]='\0';
  
  //cdsStart=region.startSite*3;
  //cdsEnd=region.endSite*3;
  //cdsLength=cdsEnd-cdsStart;
  
  strcpy(tmpBuffer, translateSeq(alignment[0]->seq));

  for (i=0;i<length;i++){
    cds[i]='.';
  }
  
  if (region.strand == 0){
    cds[0]='>';
  } else {
    cds[cdsLength-1]='<';
  }

  j=0;
  for (i=1;i+3<length;i+=3){
    cds[i]=tmpBuffer[j++];
  }

  for (i=0;i+3<length;i+=3){
    int syn_count=0;
    int non_syn_count=0;
    int currScore=0;
    int alreadySeen[4][4][4];
    int ii,jj,kk;
    
    for (ii=0;ii<4;ii++)
      for (jj=0;jj<4;jj++)
        for (kk=0;kk<4;kk++)
          alreadySeen[ii][jj][kk]=0;
    
    //if (itart || i>cdsEnd) continue;
    strncpy(&codonA,alignment[0]->seq+i,3);
    pepA=transcode[ntMap[codonA[0]]][ntMap[codonA[1]]][ntMap[codonA[2]]];
    printf("%s (%i): ", codonA, pepA);
    for (j=1; j<N; j++){
      strncpy(&codonB,alignment[j]->seq+i,3);
      
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
        //printf("%s(%i) %i ", codonB, pepB, currScore);

        if (pepA == pepB){
          if (strcmp(codonA, codonB)!=0){
            syn_count++;
          }
        } else {
          if (currScore<0){
            non_syn_count++;
          }
        }
      } else {
        non_syn_count++;
      }
    }

    syn[i]=syn_count;
    non_syn[i]=non_syn_count;
    
    //printf("%u %u\n",syn_count, non_syn_count);
    
    for (ii=0;ii<4;ii++)
      for (jj=0;jj<4;jj++)
        for (kk=0;kk<4;kk++)
          alreadySeen[ii][jj][kk]=0;
  }

  currY=startY;
  currPos=0;

  cons =  consensus(alignment);

  while (currPos<length) {

    /* Display secondary structure line */
    fprintf(outfile,"0 setgray\n");

    //strncpy(tmpBuffer,cds+currPos,columnWidth); 
    //tmpBuffer[columnWidth]='\0'; 
    
    //for (i = 0; i < columnWidth; i++){
    //  fprintf(outfile, "(%c) %.1f %.1f string\n", tmpBuffer[i],seqsX+i*fontWidth,currY);
    //}
      
    currY+=ssStep+lineStep;
    
    /* Display names, sequences and numbers */

    for (i=0; i<N; i++) {
      strncpy(tmpBuffer,alignment[i]->seq+currPos,columnWidth);
      tmpBuffer[columnWidth]='\0';
      match=0;
      for (j=0;j<(currPos+strlen(tmpBuffer));j++){
        if (alignment[i]->seq[0] != '-') match++;
      }
      
      fprintf(outfile, "(%s) %.1f %.1f string\n", alignment[i]->name,namesX,currY);

      fprintf(outfile, "(%s) %.1f %.1f string\n", tmpBuffer,seqsX,currY);

      //for (j = 0; j < columnWidth; j++){
      //  fprintf(outfile, "(%c) %.1f %.1f string\n", tmpBuffer[j],seqsX+j*fontWidth,currY);
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
    
    /*Display conservation bar*/
    
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

  for (i=0;i+3<length;i+=3){
    char *syn_color;
    char *non_syn_color;
    
    int syn_count=0;
    int non_syn_count=0;
    int currScore=0;

    syn_count=syn[i];
    non_syn_count=non_syn[i];

    if (non_syn_count>=6) {
      non_syn_count=6;
    }

    syn_color=colorMatrix[1][syn_count];
    non_syn_color=colorMatrix[0][non_syn_count];

    strncpy(&codonA,alignment[0]->seq+i,3);
    pepA=transcode[ntMap[codonA[0]]][ntMap[codonA[1]]][ntMap[codonA[2]]];
    
    block=ceil((float)(i+1)/columnWidth);
    xx=seqsX+(i-(block-1)*columnWidth)*fontWidth;

    for (j=0; j<N; j++){
      strncpy(&codonB,alignment[j]->seq+i,3);
      
      if (codonB[0] == '-' || codonB[1] == '-' || codonB[2] == '-'){
        currScore=-1;
        pepB=-99;
      } else {
        pepB=transcode[ntMap[codonB[0]]][ntMap[codonB[1]]][ntMap[codonB[2]]];

        if (pepA != -1 && pepB != -1){
          currScore=scoringMatrix[pepA][pepB];
        } else {
          currScore=-1;
        }
      }

      yy=startY+(block-1)*(lineStep*(N+2)+blockStep+consStep+rulerStep)+ssStep*(block)+(j+1)*lineStep;

      if (currScore>0){
        fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                xx,yy-1,xx+fontWidth*3,yy+fontHeight+1,syn_color);
      } else {
        fprintf(outfile, "%.1f %.1f %.1f %.1f %s box\n",
                xx,yy-1,xx+fontWidth*3,yy+fontHeight+1,non_syn_color);
      }
      
      //strncpy(tmpBuffer,cds+currPos,columnWidth); 
      //tmpBuffer[columnWidth]='\0'; 

      if (pepA == pepB){
        fprintf(outfile,"/Courier-Bold findfont\n");
        fprintf(outfile,"[10 0 0 -10 0 0] makefont setfont\n");
      } else {
        fprintf(outfile,"/Courier findfont\n");
        fprintf(outfile,"[10 0 0 -10 0 0] makefont setfont\n");
      }

      fprintf(outfile, "(%s) %.1f %.1f string\n", codonB,xx,yy);
    }
  }
    

  free(cons);

  fprintf(outfile,"showpage\n");

  fprintf(outfile," (Inhere) 360.0 723.0 string\n showpage\n");


  fclose(outfile);

  /* free(tmpBuffer); */
  /* free(ssEscaped);free(ruler); */
  
//  return 0;

}

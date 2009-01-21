 /*******************************************************************\                
 *                                                                   *
 *                        rnaz_utils.c                               *
 *                                                                   *
 *                  Helper functions for RNAz                        *
 *                                                                   *
 *                                                                   *
 *	                    Stefan Washietl                              *
 *                                                                   *
 *	   $Id: rnaz_utils.c,v 1.3 2006-10-12 13:17:42 wash Exp $        *
 *                                                                   *
 \*******************************************************************/



#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "pair_mat.h"
#include "alifold.h"
#include "rnaz_utils.h"

/********************************************************************
 *                                                                  *
 * read_clustal -- read CLUSTAL W formatted file                    *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * clust ... filehandle pointing to the file to be read             *
 * alignedSeqs ... array of strings where read sequences are stored *
 * names ... array of sequence names                                *
 *                                                                  *
 * Returns number of sequences read                                 *
 *                                                                  *
 ********************************************************************/

int read_clustal(FILE *clust, struct aln *alignedSeqs[]) {

  char *line, name[100]={'\0'}, *seq;
  int  n, nn=0, num_seq = 0;

  if (feof(clust)){
	return 0;
  }
  
  line = get_line(clust);

  while (line!=NULL) {

	if (strncmp(line,"CLUSTAL", 7)==0) {
	  break;
	}
	
	if (((n=strlen(line))<4) || isspace((int)line[0])) {
	  /* skip non-sequence line */
	  free(line); line = get_line(clust);
	  nn=0; /* reset seqence number */
	  continue;
	} 
     
	seq = (char *) space( (n+1)*sizeof(char) );
	sscanf(line,"%99s %s", name, seq);
	if (nn == num_seq) { /* first time */
	  //names[nn] = strdup(name);
	  //alignedSeqs[nn] = strdup(seq);

	  alignedSeqs[nn]=createAlnEntry(strdup(name),strdup(seq),0,0,0,'?');
	  
	  
	}
	else {
	  if (strcmp(name, alignedSeqs[nn]->name)!=0) {
		/* name doesn't match */
		free(line); free(seq);
		nrerror("ERROR: Inconsistent sequence names in CLUSTAL file");
		return 0;
	  }
	  alignedSeqs[nn]->seq = (char *)
		xrealloc(alignedSeqs[nn]->seq, strlen(seq)+strlen(alignedSeqs[nn]->seq)+1);
	  strcat(alignedSeqs[nn]->seq, seq);
	}
	nn++;
	if (nn>num_seq) num_seq = nn;
	free(seq);
	free(line);
	if (num_seq>=MAX_NUM_NAMES) {
	  nrerror("ERROR: Too many sequences in CLUSTAL file");
	  return 0;
	}
	line = get_line(clust);
  }

  alignedSeqs[num_seq] = NULL;

  if (num_seq == 0) {
	return 0;
  }
  n = strlen(alignedSeqs[0]->seq); 
  for (nn=1; nn<num_seq; nn++) {
	if (strlen(alignedSeqs[nn]->seq)!=n) {
	  fprintf(stderr, "ERROR: Sequences are of unequal length.\n");
	  return 0;
	}
  }
  return num_seq;
}



/********************************************************************
 *                                                                  *
 * read_maf -- read MAF formatted file                              *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * clust ... filehandle pointing to the file to be read             *
 * alignedSeqs ... array of strings where read sequences are stored *
 * names ... array of sequence names                                *
 *                                                                  *
 * Returns number of sequences read                                 *
 *                                                                  *
 ********************************************************************/


int read_maf(FILE *clust, struct aln *alignedSeqs[]) {

  char *line;
  int num_seq = 0;
  char** fields;
  char *name,*seq;
  int start,end,length,fullLength;
  char strand;
  int n,nn;

  if (feof(clust)){
	return 0;
  }
  
  while ((line=get_line(clust))!=NULL) {

    fields=splitFields(line);

    /* Skip empty (=only whitespace) lines */
    if (fields==NULL){
      free(line);
      continue;
    }

    /* Skip comment (#) lines */
    if (fields[0][0]=='#'){
      free(line);
      freeFields(fields);
      continue;
    }

    /* Skip any other MAF lines */
    if ((fields[0][0]=='i' || fields[0][0]=='e' || fields[0][0]=='q') && fields[0][1]=='\0'){
      free(line);
      freeFields(fields);
      continue;
    }

   
    if (fields[0][0]=='s' && fields[0][1]=='\0'){

      n=0;
      while (fields[n]!=NULL){
        n++;
      }
      if (n!=7){
        nrerror("ERROR: Invalid MAF format (number of fields in 's' line not correct)");
      }
      
      name=strdup(fields[1]);
      seq=strdup(fields[6]);
      
      if (sscanf(fields[2],"%d",&start)!=1){
        fprintf(stderr,"ERROR: Invalid MAF format"
                " (start position '%s' is not an integer)\n",fields[2]);
        exit(EXIT_FAILURE);
      }
      
      if (sscanf(fields[3],"%d",&length)!=1){
        fprintf(stderr,"ERROR: Invalid MAF format"
                " (length '%s' is not an integer)\n",fields[3]);
        exit(EXIT_FAILURE);
      }
      
      if (sscanf(fields[5],"%d",&fullLength)!=1){
        fprintf(stderr,"ERROR: Invalid MAF format"
                " (source sequence length '%s' is not an integer)\n",fields[5]);
        exit(EXIT_FAILURE);
      }
		  
      
      strand=fields[4][0];
      
      if (strand !='+' && strand !='-'){
        fprintf(stderr,"ERROR: Invalid MAF format"
                " (strand field '%s' is not '+' or '-')\n",fields[4]);
        exit(EXIT_FAILURE);
      }
      
      alignedSeqs[num_seq++]=createAlnEntry(name,seq,start,length,fullLength,strand);
      free(line);
      freeFields(fields);
      continue;
    }

    if (fields[0][0]=='a' && fields[0][1]=='\0'){
      free(line);
      freeFields(fields);
      break;
    }
  }
  
  alignedSeqs[num_seq] = NULL;

  n = strlen(alignedSeqs[0]->seq); 
  for (nn=1; nn<num_seq; nn++) {
    if (strlen(alignedSeqs[nn]->seq)!=n) {
      nrerror("ERROR: Sequences are of unequal length.");
      return 0;
    }
  }
  return num_seq;
}


/********************************************************************
 *                                                                  *
 * consensus -- Calculates consensus of alignment                   *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * AS ... array with sequences                                      *
 *                                                                  *
 * Returns string with consensus sequence                           *
 *                                                                  *
 ********************************************************************/

 char *consensus(const struct aln *AS[]) {
  char *string;
  int i,n;
  n = strlen(AS[0]->seq);
  string = (char *) space((n+1)*sizeof(char));
  for (i=0; i<n; i++) {
    int s,c,fm, freq[8] = {0,0,0,0,0,0,0,0};
    for (s=0; AS[s]!=NULL; s++) 
      freq[encode_char(AS[s]->seq[i])]++;
    for (s=c=fm=0; s<8; s++) /* find the most frequent char */
      if (freq[s]>fm) {c=s, fm=freq[c];}
    if (s>4) s++; /* skip T */
    string[i]=Law_and_Order[c];
  }
  return string;
}


/********************************************************************
 *                                                                  *
 * meanPairID -- Calculates mean pairwise identity of alignment     *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * AS ... array with sequences                                      *
 *                                                                  *
 * Returns mean pair ID in percent                                  *
 *                                                                  *
 ********************************************************************/


 double meanPairID(const struct aln *AS[]) {

  int i,j,k,matches,pairs,length;

  matches=0;
  pairs=0;

  length=strlen(AS[0]->seq);
  
  for (i=0;AS[i]!=NULL;i++){
	for (j=i+1;AS[j]!=NULL;j++){
	  for (k=0;k<length;k++){
		if ((AS[i]->seq[k]!='-') || (AS[j]->seq[k]!='-')){
		  if (AS[i]->seq[k]==AS[j]->seq[k]){
			matches++;
		  }
		  pairs++;
		}
	  }
	}
  }

  return (double)(matches)/pairs*100;
  
}

/********************************************************************
 *                                                                  *
 * revAln -- Reverse complements sequences in an alignment          *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * AS ... array with sequences which is rev-complemented in place   *
 *                                                                  *
 ********************************************************************/

void revAln(struct aln *AS[]) {

  int i,j,length;
  char *tmp;
  char letter;
  length=strlen(AS[0]->seq);
  
  for (i=0;AS[i]!=NULL;i++){
    tmp = (char *) space((unsigned) length+1);
    for (j=length-1;j>=0;j--){
      letter=AS[i]->seq[j];
      switch(letter){
      case 'T': letter='A'; break;
      case 'U': letter='A'; break;
      case 'C': letter='G'; break;
      case 'G': letter='C'; break;
      case 'A': letter='T'; break;
      }
	  tmp[length-j-1]=letter;
	}
	tmp[length]='\0';
	strcpy(AS[i]->seq,tmp);

	if (AS[i]->strand =='+') {
	  AS[i]->strand ='-';
	} else if (AS[i]->strand =='-'){
	  AS[i]->strand ='+';
	}
		
	free(tmp);
	tmp=NULL;
  }
}

/********************************************************************
 *                                                                  *
 * sliceAln -- Gets a slice of an alignment                         *
 *                                                                  *
 ********************************************************************
 *                                                                  *
 * sourceAln ... array with sequences of source alignment           *
 * destAln   ... pointer to array where slice is stored             *
 * from, to  ... specifies slice, first column is column 1          *
 *                                                                  *
 ********************************************************************/

void sliceAln(const struct aln *sourceAln[], struct aln *destAln[],
			  int from, int to){

  int i;
  char *slice;
 
  for (i=0;sourceAln[i]!=NULL;i++){
	slice=(char *) space((unsigned) (to-from+2));
	strncpy(slice,(sourceAln[i]->seq)+from-1,(to-from+1));
	/* coordinates are NOT changed*/
	destAln[i]=createAlnEntry(strdup(sourceAln[i]->name),
							  slice,
							  sourceAln[i]->start, 
							  sourceAln[i]->length,
							  sourceAln[i]->fullLength,
							  sourceAln[i]->strand);
  }
  destAln[i]=NULL;
}

/********************************************************************
 *                                                                  *
 * freeAln -- Frees memory of alignment array                       *
 *                                                                  *
 ********************************************************************/

 void freeAln(struct aln *AS[]){

  int i;
  for (i=0;AS[i]!=NULL;i++){
	freeAlnEntry(AS[i]);
  }
  AS=NULL;
}


struct aln* createAlnEntry(char* name, char* seq, int start, int length, int fullLength, char strand){
  struct aln* entry;

  entry=(struct aln*)space(sizeof(struct aln));

  entry->name=name;
  entry->seq=seq;
  entry->start=start;
  entry->length=length;
  entry->fullLength=fullLength;
  entry->strand=strand;

  return entry;
}

 void freeAlnEntry(struct aln* entry){

   free(entry->name);
   free(entry->seq);
   if (entry->fullSeq != NULL) free(entry->fullSeq);
   free(entry);
   
}


void printAln(const struct aln* AS[]){
  int i;
  for (i=0;AS[i]!=NULL;i++){
    printf("%s %s\n",AS[i]->name,AS[i]->seq);
  }
}

 int checkFormat(FILE *file){

  char *line; 
  char** fields;
  
  while ((line=get_line(file)) != NULL){

	fields=splitFields(line);

	/* Skip empty (=only whitespace) and comments */

	if (fields==NULL){
	  free(line);
	  continue;
	}
	
	if (fields[0][0]=='#'){
	  free(line);
	  freeFields(fields);
	  continue;
	}

	/* Identify "CLUSTAL" header => CLUSTAL*/
	if (strcmp(fields[0],"CLUSTAL")==0){
	  free(line);
	  freeFields(fields);
	  return(CLUSTAL);
	}

	/* Identitfy "a" header => MAF */
	if (fields[0][0]=='a' && fields[0][1]=='\0'){
	  free(line);
	  freeFields(fields);
	  return(MAF);
	}

	/* Unknown format if the first non-empty, non-command data in the
	   file is non of the above */
	free(line);
	freeFields(fields);
	return(0);
  }

  free(line);
  nrerror("ERROR: Empty alignment file\n");
   
}



char** splitFields(char* string){

  char c;
  char* currField;
  char** output=NULL;
  int* seps;
  int nSep;
  int nField=0;
  int i=0;

  if (strlen(string)==0 || string==NULL){
	return NULL;
  }

  /* First find all characters which are whitespaces and store the
	 positions in the array seps */
    
  seps=(int *)space(sizeof(int));
  seps[0]=-1;
  nSep=1;
  
  while ((c=string[i])!='\0' && (c!='\n')){
	if (isspace(c)){
	  seps=(int*)xrealloc(seps,sizeof(int)*(nSep+1));
	  seps[nSep++]=i;
	}
	i++;
  }

  seps=(int*)xrealloc(seps,sizeof(int)*(nSep+1));
  seps[nSep]=strlen(string);


  /* Then go through all intervals in between of two whitespaces (or
	 end or start of string) and store the fields in the array
	 "output"; if there are two adjacent whitespaces this is ignored
	 resulting in a behaviour like "split /\s+/" in perl */
  
  for (i=0;i<nSep;i++){

	int start=seps[i];
	int stop=seps[i+1];
	int length=(stop-start);
	int notSpace,j;

	
	currField=(char *)space(sizeof(char)*(length+1));
	strncpy(currField,string+start+1,length-1);
	currField[length]='\0';

	/* check if field is not only whitespace */
	notSpace=0;
	j=0;
	while (c=currField[j]!='\0'){
	  if (!isspace(c)){
		notSpace=1;
		break;
	  }
	}

	if (notSpace){
	  output=(char**)xrealloc(output,sizeof(char**)*(nField+1));
	  output[nField++]=currField;
	  currField=NULL;
	} else {
	  free(currField);
	  currField=NULL;
	}

	//printf("%s|\n",output[nField-1]);
  }

  if (nField==0){
	return NULL;
  }

  
  output=(char**)xrealloc(output,sizeof(char**)*(nField+1));
  output[nField]=NULL;
  
  free(seps);
  return output;
  
}

 char** splitLines(char* string){

  char c;
  char* currLine=NULL;
  char** output=NULL;
  int i=0;
  int currLength=0;
  int lineN=0;

  while ((c=string[i])!='\0'){

	if (c=='\n'){
	  output=(char**)xrealloc(output,sizeof(char**)*(lineN+1));
	  currLine=(char*)xrealloc(currLine,sizeof(char)*(currLength+1));
	  currLine[currLength]='\0';
	  output[lineN]=currLine;
	  currLength=0;
	  currLine=NULL;
	  lineN++;
	} else {

	  currLine=(char*)xrealloc(currLine,sizeof(char)*(currLength+1));
	  currLine[currLength]=c;
	  currLength++;
	}
	i++;
  }

  output=(char**)xrealloc(output,sizeof(char**)*(lineN+1));
  output[lineN]=NULL;
  
  return output;

}

// for both splitLines and splitFields
 void freeFields(char** fields){

  int i=0;
  while (fields[i]!=NULL){
	free(fields[i++]);
  }
  free(fields);
}

double combPerPair(struct aln *AS[],char* structure){

  int* stack;
  int stackN;
  int i,j,k,l,x,y;
  int nPairs, nCombs;
  char c;
  int base1,base2;
  
  int pairMatrix[4][4]={{0,0,0,1},
						{0,0,1,1},
						{0,1,0,0},
						{1,1,0,0}};

  int seenMatrix[4][4];

  nPairs=0;
  nCombs=0;
  
  stack=(int*)space(sizeof(int)*strlen(structure));
  stackN=0;
  i=0;
 
  while ((c=structure[i])!='\0'){
	if (c=='('){
	  stack[stackN++]=i;
	}
	if (c==')'){
	  j=stack[--stackN];
	  //printf("Base pair: %d %d\n",j,i);
	  for (x=0;x<4;x++){
		for (y=0;y<4;y++){
		  seenMatrix[x][y]=0;
		}
	  }
	  k=0;
	  while (AS[k]!=NULL){
		base1=encodeBase(AS[k]->seq[j]);
		base2=encodeBase(AS[k]->seq[i]);

		if (base1==-1 || base2==-1){
		  k++;
		  continue;
		}
		
		if (pairMatrix[base1][base2]){
		  if (!seenMatrix[base1][base2]){
			nCombs++;
			seenMatrix[base1][base2]=1;
		  }
		}
		//printf("Seq %d: %d %d\n",k,base1,base2);
		k++;
	  }
	  nPairs++;
	}
	i++;
  }

  if (nPairs>0){
	return((double)nCombs/nPairs);
  } else {
	return 0.0;
  }
}

int encodeBase(char base){

  char b;

  b=toupper(base);

  if (b == 'A') return 0;
  if (b == 'G') return 1;
  if (b == 'C') return 2;
  if (b == 'T') return 3;
  if (b == 'U') return 3;

  return -1;
}

void printAlnMAF(FILE *out, const struct aln* AS[],int printU){

  int i,j,N;
  int L;
  char* tmpString;

  L=strlen(AS[0]->seq);

  fprintf(out, "a score=0\n");

  for (i=0;AS[i]!=NULL;i++){

    tmpString=AS[i]->seq;

    j=0;
    
    while (tmpString[j]){
      if (!printU){
        if (tmpString[j]=='U'){
          tmpString[j]='T';
        }
      }
      j++;
    }
  }
  
  for (i=0;AS[i]!=NULL;i++){
    fprintf(out, "s %s %i %i %c %i %s\n",AS[i]->name,AS[i]->start,AS[i]->length,AS[i]->strand, AS[i]->fullLength, AS[i]->seq);
  }
  fprintf(out, "\n");
  
}

void pruneAln(char* species, struct aln* alignment[]){

  int i,j,N,x, counter;
  int seen;
  struct aln *newAln[MAX_NUM_NAMES];
  char** list;

  list=splitString(species,",");

  counter=0;

  for (i=0;alignment[i]!=NULL;i++){
    seen=0;
    for (x=0;list[x]!=NULL;x++){
      if (strncmp(list[x],alignment[i]->name,strlen(list[x]))==0){
        //printf("matching %s\n",alignment[i]->name);
        seen=1;
      }
    }
    if (seen==0){
      //printf("deleting %s\n",alignment[i]->name);
      for (j=i;alignment[j]!=NULL;j++){
        alignment[j]=alignment[j+1];
      }
      i--;
    }
    //printAlnMAF(stdout,(const struct aln**)alignment,0); 
  }
}

char** splitString(char* string, char* separators){

  char sp[10000]; /* max length hardcoded */
  char* pch;
  char** output;
  char* currField;
  int counter=0;

  /* strtok only functions that way, cannot pass char pointer
     directly */
  strcpy(sp,string);

  pch=strtok(sp,",");

  counter=0;

  output=(char**)malloc(sizeof(char*));
  currField=(char*)malloc(sizeof(char)*(strlen(pch)+1));
  strcpy(currField, pch);
  output[counter++]=currField;
  
  while (pch != NULL){
    //printf ("-------> %s\n",pch);
    pch=strtok(NULL,",");
    if (pch == NULL){
      output=(char**)realloc(output,sizeof(char*)*(counter+1));
      output[counter]=NULL;
      break;
    }
    output=(char**)realloc(output,sizeof(char*)*(counter+1));
    currField=(char*)malloc(sizeof(char)*(strlen(pch)+1));
    strcpy(currField, pch);
    output[counter++]=currField;
  }
  
  return output;

}

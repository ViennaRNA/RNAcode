#include "code.h"

/* Data structure holds encoded amino acid for codon of three encoded nucleic acids  */ 

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

/* BLOSUM substitution matrix. Holds score for two encoded amino acids */ 

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


/*********************************************************************
  translateSeq

  seq ... Nucleic acid sequence (unencoded)

  Translates a nucleic acid sequence to a protein. Return protein
  as string in plain text (unencoded).

*********************************************************************/ 

char* translateSeq(char* seq){

  int i,j;
  int* encodedSeq;
  char* peptideSeq;

  encodedSeq=encodeSeq(seq);

  peptideSeq=(char*)malloc(sizeof(char)*strlen(seq));
  
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

/*********************************************************************
  encodeSeq

  seq ... Nucleic acid sequence (unencoded)

  Encodes a sequence from a plain text string into a string of
  integers.

*********************************************************************/ 

int* encodeSeq(char* seq){

  int* encoded;
  int i;
  char c;

  encoded=(int*)malloc(sizeof(int)*(strlen(seq)+1));

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

/*********************************************************************
  encodeAA

  aa ... Amino acid sequence (unencoded)

  Encodes an amino acid sequence into a sequence of integers

*********************************************************************/ 

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

/*********************************************************************
  decodeAA

  encodedAA ... Amino acid sequence (encoded)

  Decodes an amino acid encoded as integers back to plain sequence.

*********************************************************************/ 

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

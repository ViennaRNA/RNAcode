#ifndef _RNAZ_UTILS_H_
#define _RNAZ_UTILS_H_



#define PRIVATE static
#define MAX_NUM_NAMES 500
#define MIN2(A, B) ((A) < (B) ? (A) : (B))

enum alnFormat {UNKNOWN=0, CLUSTAL=1, MAF=2};

struct aln {
  char *name;
  char *seq; 
  char *fullSeq; /* Remembers full sequence when gaps are stripped in seq*/
  int start;
  int length;
  int fullLength;
  char strand;
};


int read_clustal(FILE *clust,
						 struct aln *alignedSeqs[]);

int read_maf(FILE *clust,
						 struct aln *alignedSeqs[]);


char *consensus(const struct aln *AS[]);

double meanPairID(const struct aln *AS[]);

void revAln(struct aln *AS[]);

double combPerPair(struct aln *AS[],char* structure);

int encodeBase(char base);

void sliceAln(const struct aln *sourceAln[], struct aln *destAln[],
					  int from, int to);

void freeAln(struct aln *AS[]);


struct aln* createAlnEntry(char* name, char* seq, int start, int length, int fullLength, char strand); 
void freeAlnEntry(struct aln* entry);
void printAln(const struct aln* AS[]);

int checkFormat(FILE *file);

char** splitFields(char* string);
void freeFields(char** fields);

char** splitLines(char* string);

#endif

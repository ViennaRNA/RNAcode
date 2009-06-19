#include "tree.h"
#include "rnaz_utils.h"

TNode* getLCA(TTree* tree, TNode* nodeA, TNode* nodeB);
void tree2aln(TTree* tree, struct aln *alignment[]);
TTree* string2tree(char* treeString);

void simulateTree(TTree* tree, float freqs[], float kap, int L);

float** getDistanceMatrix(TTree* tree, struct aln *alignment[]);

void freeSeqgenTree(TTree* tree);

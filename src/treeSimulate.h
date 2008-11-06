#include "tree.h"
#include "rnaz_utils.h"

TNode* getLCA(TTree* tree, TNode* nodeA, TNode* nodeB);
void tree2aln(TTree* tree, struct aln *alignment[]);
TTree* string2tree(char* treeString);

void simulateTree(TTree* tree, double freqs[], double kap, int L);

double** getDistanceMatrix(TTree* tree, struct aln *alignment[]);

void freeSeqgenTree(TTree* tree);

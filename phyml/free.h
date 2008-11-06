/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef FREE_H
#define FREE_H

void Free_All_Nodes_Light(arbre *tree);
void Free_All_Edges_Light(arbre *tree);
void Free_Mat(matrix *mat);
void Free_Partial_Lk(phydbl ****p_lk, int len, int n_catg);
void Free_Tree(arbre *tree);
void Free_Edge(edge *b);
void Free_Node(node *n);
void Free_Cseq(allseq *alldata);
void Free_Seq(seq **d, int n_otu);
void Free_All(seq **d, allseq *alldata, arbre *tree);
void Free_SubTree(edge *b_fcus, node *a, node *d, arbre *tree);
void Free_Tree_Ins_Tar(arbre *tree);
void Free_Tree_Lk(arbre *tree);
void Free_NNI(arbre *tree);
void Free_Edge_P_Lk_Struct(edge *b, arbre *tree);
void Free_Node_Lk(node *n);
void Free_Edge_Lk(arbre *tree, edge *b);
void Free_Model(model *mod);
void Free(void *p);
void Free_Input(option *input);
void Free_Reachable(arbre *tree);
void Free_St(superarbre *st);
void Free_Eigen(eigen *eigen_struct);
void Free_Triplet(ttriplet *t);
void Free_Tree_Pars(arbre *tree);
void Free_Edge_Pars(edge *b, arbre *tree);
void Free_One_Spr(spr *this_spr);
void Free_Spr_List(arbre *tree);
void Free_Actual_CSeq(allseq *data);
void Free_Prefix_Tree(pnode *n, int size);
void Free_Pnode(pnode *n);
void Free_Edge_Labels(edge *b);
void Free_Rates(int n_otu, trate *r);

#endif

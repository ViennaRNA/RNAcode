/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef ML_H
#define ML_H


void Init_Tips_At_One_Site_Nucleotides_Float(char state,phydbl *p_lk);
void Init_Tips_At_One_Site_AA_Float(char aa,phydbl *p_lk);
void Get_All_Partial_Lk(arbre *tree,edge *b_fcus,node *a,node *d);
void Get_All_Partial_Lk_Scale(arbre *tree,edge *b_fcus,node *a,node *d);
void Post_Order_Lk(node *pere,node *fils, arbre *tree);
void Pre_Order_Lk(node *pere,node *fils, arbre *tree);
void Lk(arbre *tree);
void Site_Lk(arbre *tree);
phydbl Lk_At_Given_Edge(edge *b_fcus,arbre *tree);
phydbl Return_Lk(arbre *tree);
phydbl Return_Abs_Lk(arbre *tree);
matrix *ML_Dist(allseq *data,model *mod);
phydbl Lk_Given_Two_Seq(allseq *data,int numseq1,int numseq2,phydbl dist,model *mod,phydbl *loglk);
void Unconstraint_Lk(arbre *tree);
void Update_P_Lk(arbre *tree,edge *b_fcus,node *n);
void Make_Tree_4_Lk(arbre *tree,allseq *alldata,int n_site);
void Init_P_Lk_Tips_Double(arbre *tree);
void Init_P_Lk_Tips_Int(arbre *tree);
void Init_P_Lk_At_One_Node(node *a, arbre *tree);
void Update_PMat_At_Given_Edge(edge *b_fcus, arbre *tree);
void Sort_Sites_Based_On_Lk(arbre *tree);
void Get_Partial_Lk_Scale(arbre *tree, edge *b_fcus, node *a, node *d);
void Get_Partial_Lk(arbre *tree, edge *b_fcus, node *a, node *d);
void Init_Tips_At_One_Site_Nucleotides_Int(char state, short int *p_pars);
void Init_Tips_At_One_Site_AA_Int(char aa, short int *p_pars);
void Update_P_Lk_Along_A_Path(node **path, int path_length, arbre *tree);
phydbl Lk_Dist(phydbl *F, phydbl dist, model *mod);
phydbl Update_Lk_At_Given_Edge(edge *b_fcus, arbre *tree);
void Update_P_Lk_Greedy(arbre *tree, edge *b_fcus, node *n);
void Get_All_Partial_Lk_Scale_Greedy(arbre *tree, edge *b_fcus, node *a, node *d);
phydbl Lk_Core(edge *b, arbre *tree);
phydbl Lk_Triplet(node *a, node *d, arbre *tree);
void Print_Lk_Given_Edge_Recurr(node *a, node *d, edge *b, arbre *tree);
phydbl *Post_Prob_Rates_At_Given_Edge(edge *b, phydbl *post_prob, arbre *tree);
phydbl Lk_With_MAP_Branch_Rates(arbre *tree);



#endif







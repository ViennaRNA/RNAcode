/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef PARS_H
#define PARS_H

void Make_Tree_4_Pars(arbre *tree, allseq *alldata, int n_site);
int  Pars(arbre *tree);
void Post_Order_Pars(node *a, node *d, arbre *tree);
void Pre_Order_Pars(node *a, node *d, arbre *tree);
void Get_Partial_Pars(arbre *tree, edge *b_fcus, node *a, node *d);
void Site_Pars(arbre *tree);
void Init_Ui_Tips(arbre *tree);
void Update_P_Pars(arbre *tree, edge *b_fcus, node *n);
int Pars_At_Given_Edge(edge *b, arbre *tree);
void Get_All_Partial_Pars(arbre *tree, edge *b_fcus, node *a, node *d);
int Update_Pars_At_Given_Edge(edge *b_fcus, arbre *tree);
void Init_P_Pars_Tips(arbre *tree);
void Get_Step_Mat(arbre *tree);
int Pars_Core(edge *b, arbre *tree);
int One_Pars_Step(edge *b,arbre *tree);

#endif

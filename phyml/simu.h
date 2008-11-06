/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef CURR_H
#define CURR_H

void Simu_Loop(arbre *tree);
void Simu(arbre *tree,int n_step_max);
void Select_Edges_To_Swap(arbre *tree,edge **sorted_b,int *n_neg);
void Update_Bl(arbre *tree,phydbl fact);
void Make_N_Swap(arbre *tree,edge **b,int beg,int end);
int Make_Best_Swap(arbre *tree);
int Mov_Backward_Topo_Bl(arbre *tree,phydbl lk_old,edge **tested_b,int n_tested);
void Unswap_N_Branch(arbre *tree,edge **b,int beg,int end);
void Swap_N_Branch(arbre *tree,edge **b,int beg,int end);
void Check_NNI_Scores_Around(node *a, node *d, edge *b, phydbl *best_score);
int Mov_Backward_Topo_Pars(arbre *tree, int pars_old, edge **tested_b, int n_tested);
void Simu_Pars(arbre *tree, int n_step_max);

#endif

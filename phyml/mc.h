/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef MC_H
#define MC_H

void MC_Least_Square_Node_Times_Pre(node *a, node *d, phydbl *A, phydbl *b, int n, arbre *tree);
int  MC_main(int argc, char **argv);
void MC_Bl_From_T_Post(node *a, node *d, edge *b, arbre *tree);
void MC_Bl_From_T(arbre *tree);
void MC_Optimize_Node_Times_Serie(node *a, node *d, edge *b, arbre *tree);
void MC_Round_Optimize(arbre *tree);
void MC_Print_Node_Times(node *a, node *d, arbre *tree);
edge *MC_Find_Best_Root_Position(arbre *tree);
void MC_Least_Square_Node_Times(edge *e_root, arbre *tree);
void MC_Mult_Time_Stamps(arbre *tree);
void MC_Div_Time_Stamps(arbre *tree);
void MC_Optimize_Tree_Height(arbre *tree);
void MC_Adjust_Node_Times(arbre *tree);
void MC_Adjust_Node_Times_Pre(node *a, node *d, arbre *tree);
void MC_Optimize_Root_Height(arbre *tree);
void MC_Estimate_Branch_Rates(arbre *tree);
edge *MC_Find_Best_Root_Position_Approx(arbre *tree);
void MC_Estimate_Branch_Rate_Parameter(arbre *tree);
phydbl MC_Classify_Branch_In_Rate_Class(arbre *tree);
void MC_Compute_Rates_And_Times_Least_Square_Adjustments(arbre *tree);
void MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(node *a, node *d, edge *b, arbre *tree);
void MC_Classify_Branch_Rates(arbre *tree);

#endif

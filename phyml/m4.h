/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef M4_H
#define M4_H

int M4_main(int argc, char **argv);
void M4_Make_Complete(int n_h, int n_o, m4 *m4mod);
m4 *M4_Make_Light(int n_o);
void M4_Free_M4_Model(m4 *m4mod);
void M4_Init_Qmat(m4 *m4mod, allseq *data, model *mod);
void M4_Update_Qmat(m4 *m4mod, model *mod);
void M4_Init_Model(m4 *m4mod, allseq *data, model *mod);
void M4_Init_P_Lk_Tips_Double(arbre *tree);
void M4_Init_P_Lk_Tips_Int(arbre *tree);
void M4_Post_Prob_H_Class_Edge_Site(edge *b, phydbl ****integral, phydbl *postprob, arbre *tree);
phydbl ****M4_Integral_Term_On_One_Edge(edge *b, arbre *tree);
phydbl ***M4_Compute_Proba_Hidden_States_On_Edges(arbre *tree);
void M4_Free_Integral_Term_On_One_Edge(phydbl ****integral, arbre *tree);
void M4_Compute_Posterior_Mean_Rates(phydbl ***post_probs, arbre *tree);
void M4_Scale_Br_Len(arbre *tree);
void M4_Detect_Site_Switches_Experiment(arbre *tree);
m4 *M4_Copy_M4_Model(model *ori_mod, m4 *ori_m4mod);
void M4_Posterior_Prediction_Experiment(arbre *tree);
void M4_Set_M4mod_Default(m4 *m4mod);
phydbl **M4_Site_Branch_Classification(phydbl ***post_probs, arbre *tree);
void M4_Site_Branch_Classification_Experiment(arbre *tree);

#endif

#ifndef MULTIGENE_H
#define MULTIGENE_H

#include "utilities.h"

void Menu_Supertree(option *input);
void Mg_Print_Nodes(node *a, node *d, superarbre *st);
superarbre *Mg_Make_Superarbre_Light(option *input);
void Mg_Make_Superarbre_Full(superarbre *st, option *input, allseq **data);
void Mg_Get_List_Of_Reachable_Tips(node *a, node *d, allseq **data, superarbre *st);
void Mg_Get_List_Of_Reachable_Tips_Pre(node *a, node *d, superarbre *st);
void Mg_Get_List_Of_Reachable_Tips_Post(node *a, node *d, superarbre *st);
void Mg_Prune_St_Topo(arbre *tree, allseq *data, superarbre *st);
void Mg_Match_St_Nodes_In_Gt_Recurr(node *a_gt, node *d_gt, node *a_st, node *d_st, arbre *gt, superarbre *st);
void Mg_Match_St_Nodes_In_Gt(arbre *tree, superarbre *st);
void Mg_Match_St_Edges_In_Gt(arbre *gt, superarbre *st);
void Mg_Match_St_Edges_In_Gt_Recurr(node *a, node *d, node *a_st, node *d_st, arbre *gt, superarbre *st);
void Mg_Simu(superarbre *tr);
int Mg_Mov_Backward_Topo_Bl(superarbre *st, phydbl lk_old, edge **tested_b, int n_tested);
int Mg_Get_Species_Found_In_St(superarbre *st, allseq *data);
void Mg_Map_St_Nodes_In_Gt_Pre(node *a_st, node *d_st, arbre *gt, superarbre *st);
void Mg_Map_St_Nodes_In_Gt_Post(node *a_st, node *d_st, arbre *gt, superarbre *st);
void Mg_Map_St_Nodes_In_Gt(arbre *gt, superarbre *st);
void Mg_Map_St_Nodes_In_Gt_One_Edge(node *a_st, node *d_st, edge *b_st, arbre *gt, superarbre *st);
void Mg_Map_St_Edges_In_Gt(arbre *gt, superarbre *st);
phydbl Mg_Lk(superarbre *st);
int Mg_Pars(superarbre *st);
int Mg_Spr(phydbl init_lnL, superarbre *st);
void Mg_Speed_Spr(superarbre *st);
int Map_Spr_Move(edge *st_pruned, edge *st_target, node *st_link, arbre *gt, superarbre *st);
void Mg_Test_All_Spr_Targets(edge *pruned, node *n_link, superarbre *st);
void Mg_Test_One_Spr_Target_Recur(node *a, node *d, edge *target, edge *pruned, node *n_link, superarbre *st);
void Mg_Test_One_Spr_Target(edge *st_p, edge *st_t, node *n_link, superarbre *st);
int Mg_Test_List_Of_Regraft_Pos(spr **st_spr_list, int list_size, superarbre *st);
int Mg_Try_One_Spr_Move(spr *st_move, superarbre *st);
void Mg_Map_Gt_Edges_In_St(arbre *gt, superarbre *st);
void Mg_NNI(edge *st_b, superarbre *st);
void Mg_Swap(node *st_a, node *st_b, node *st_c, node *st_d, superarbre *st);
void Mg_Set_Bl(double **bl, superarbre *st);
void Mg_Restore_Br_Len(superarbre *st);
void Mg_Record_Br_Len(superarbre *st);
phydbl Mg_Lk_At_Given_Edge(edge *st_b, superarbre *st);
phydbl Mg_Update_Lk_At_Given_Edge(edge *st_b, superarbre *st);
void Mg_Fill_Model_Partitions_Table(superarbre *st);
phydbl Mg_Br_Len_Brent(edge *st_b, superarbre *tree);
void Mg_Initialise_Bl_Partition(superarbre *st);
void Mg_Update_P_Lk(edge *st_b, node *st_n, superarbre *st);
void Mg_Optimize_Br_Len_Serie(node *st_a, node *st_d, edge *st_b, superarbre *st);
void Mg_Update_Bl_Swaped(edge **st_b, int n, superarbre *st);
void Mg_Update_Bl(phydbl fact, superarbre *st);
void Mg_Make_N_Swap(edge **st_b, int beg, int end, superarbre *st);
void Mg_Do_Mapping(superarbre *st);
void Mg_Update_PMat(edge *st_b, superarbre *st);
void Mg_Print_Bl(superarbre *st);
void Mg_Check_Extra_Taxa(superarbre *st);
int MG_main(int argc, char **argv);
#endif


/*
** spr.h: Header file for the SPR routines.
**
** Wim Hordijk   Last modified: 28 August 2006
*/

#ifndef _SPR_H_
#define _SPR_H_

#include "utilities.h"

#define ALL   1
#define BEST  2
#define ONE   3

/*
** _move_: Structure for holding the relevant information for candidate SPR moves.
*/

typedef struct
{
  node   *v_prune, *u_prune, *v_n, *v_nx1, *u_n, **path;
  edge   *e_prune, *e_regraft;
  phydbl  l_connect, l_est[3], delta_lk, d_L, d_up_v, d_un_v;
  int     dist, rgrft_rank, optim_rank, globl_rank;
} _move_;



void Init_SPR          (arbre *tree);
void Clean_SPR         (arbre *tree);
void Optim_SPR         (arbre *tree, int max_size, int method);
int  Perform_SPR_Moves (arbre *tree, int max_size);
int  Perform_Best_SPR  (arbre *tree, int max_size);
int  Perform_One_SPR   (arbre *tree, int max_size);

void Calc_Tree_Length (edge *e_prune, node *v_prune, arbre *tree);
void Tree_Length      (node *v_prune, node *u_prune, node *v_n, node *v_n_1,
		       node *v_nx1, node *v_0, node *u_n, phydbl d_up_v_1,
		       phydbl d_uu, phydbl d_L_1, int n, arbre *tree);
int  Est_Lk_Change    (edge *e_prune, node *v_prune, arbre *tree);
int  Best_Lk_Change   (edge *e_prune, node *v_prune, arbre *tree);
void Make_Move        (_move_ *move, int type, arbre *tree);
int  Find_Optim_Local (arbre *tree);
int  Find_Optim_Globl (arbre *tree);
void Prune            (edge *e, node *v, edge **e_connect, edge **e_avail,
		       arbre *tree);
void Regraft          (edge *e, node *v, edge *avail, arbre *tree);
void PostOrder_v      (arbre *tree, node *v, edge *e);
void PostOrder_w      (arbre *tree, node *v, edge *v_e, node *w, edge *e);


#endif  /* _SPR_H_ */


/*
** EOF: spr.h
*/

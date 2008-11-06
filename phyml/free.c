/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "utilities.h"
#include "free.h"
#include "m4.h"

/*********************************************************/

void Free_All_Nodes_Light(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) 
    Free_Node(tree->noeud[i]);
}

/*********************************************************/

void Free_All_Edges_Light(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-3) 
    if(tree->t_edges[i])
      Free_Edge(tree->t_edges[i]);
}


/*********************************************************/

void Free_Mat(matrix *mat)
{
  int i;

  For(i,mat->n_otu)
    {
      Free(mat->P[i]);
      Free(mat->Q[i]);
      Free(mat->dist[i]);
      Free(mat->name[i]);
    }

  Free(mat->P);
  Free(mat->Q);
  Free(mat->dist);
  Free(mat->name);
  Free(mat->tip_node);
      
  Free(mat->on_off);
  Free(mat);
}

/*********************************************************/

void Free_Partial_Lk(phydbl ****p_lk, int len, int n_catg)
{
  int i,j;

  For(i,len)
    {
      For(j,n_catg) Free((*p_lk)[i][j]);
      Free((*p_lk)[i]);
    }
  Free((*p_lk));
  (*p_lk) = NULL;
}

/*********************************************************/

void Free_Tree(arbre *tree)
{
  int i,j,k;
  edge *b;
  node *n;

  For(i,2*tree->n_otu-2) Free(tree->t_dir[i]);
  Free(tree->t_dir);

  if(tree->has_bip)
    {
      For(i,2*tree->n_otu-2)
	{
	  Free(tree->noeud[i]->bip_size);
	  For(j,3)
	    {
	      Free(tree->noeud[i]->bip_node[j]);
	      For(k,tree->n_otu) Free(tree->noeud[i]->bip_name[j][k]);
	      Free(tree->noeud[i]->bip_name[j]);
	    }
	  Free(tree->noeud[i]->bip_node);
	  Free(tree->noeud[i]->bip_name);
	}
    }
  
  Free(tree->curr_path);

  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      Free_Edge(b);
    }
  Free(tree->t_edges);


  For(i,2*tree->n_otu-2)
    {
      n = tree->noeud[i];
      Free_Node(n);
    }
  Free(tree->noeud);


  For(i,tree->n_dead_edges)
    Free_Edge(tree->t_dead_edges[i]);

  For(i,tree->n_dead_nodes)
    Free_Node(tree->t_dead_nodes[i]);
  
  Free(tree->t_dead_edges);
  Free(tree->t_dead_nodes);


  Free(tree);
}

/*********************************************************/

void Free_Edge_Labels(edge *b)
{
  int i;
  For(i,b->n_labels+b->n_labels%BLOCK_LABELS) Free(b->labels[i]);
  Free(b->labels);
  b->labels = NULL;
}

/*********************************************************/

void Free_Edge(edge *b)
{
  Free_Edge_Labels(b);
  Free(b);
}

/*********************************************************/

void Free_Node(node *n)
{
  int i;

  Free(n->b);
  Free(n->v);
  Free(n->l);
  Free(n->score);
  Free(n->name);

  if(n->list_of_reachable_tips)
    {
      For(i,3) Free(n->list_of_reachable_tips[i]);
      Free(n->list_of_reachable_tips);
      Free(n->n_of_reachable_tips);
    }

  Free(n);
}

/*********************************************************/

void Free_Cseq(allseq *data)
{
  int i;
  
  Free(data->invar);
  Free(data->wght);
  Free(data->ambigu);
  Free(data->b_frq);
  Free(data->sitepatt);
  For(i,data->n_otu)
    {
      Free(data->c_seq[i]->name);
      if(data->c_seq[i]->state) 
	{
	  Free(data->c_seq[i]->state);
	  if(data->c_seq[i]->is_ambigu) Free(data->c_seq[i]->is_ambigu);
	}
      Free(data->c_seq[i]);
    }
  Free(data->c_seq);
  Free(data);
}

/*********************************************************/

void Free_Seq(seq **d, int n_otu)
{
  int i;
  For(i,n_otu)
    {
      Free(d[i]->name);
      Free(d[i]->state);
      if(d[i]->is_ambigu) Free(d[i]->is_ambigu);
      Free(d[i]);
    }
  Free(d);
}


/*********************************************************/

void Free_All(seq **d, allseq *alldata, arbre *tree)
{
  Free_Cseq(alldata);
  Free_Seq(d,tree->n_otu);
  Free_Tree(tree);
}      

/*********************************************************/
void Free_SubTree(edge *b_fcus, node *a, node *d, arbre *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Free_SubTree(d->b[i],d,d->v[i],tree);
	      Free_Edge(d->b[i]);
	      Free_Node(d->v[i]);
	    }
	}
    }
}

/*********************************************************/
void Free_Tree_Ins_Tar(arbre *tree)
{
  return;
}

/*********************************************************/

void Free_Tree_Pars(arbre *tree)
{
  int i;
  
  Free(tree->step_mat);
  Free(tree->site_pars);
  For(i,2*tree->n_otu-3)
    Free_Edge_Pars(tree->t_edges[i],tree);
}

/*********************************************************/

void Free_Edge_Pars(edge *b, arbre *tree)
{
  int i;

  Free(b->pars_l);
  Free(b->pars_r);
  
  For(i,tree->data->crunch_len) 
    {
      Free(b->p_pars_l[i]);
      Free(b->p_pars_r[i]);
    }
  
  Free(b->ui_l);
  Free(b->ui_r);
  Free(b->p_pars_l);
  Free(b->p_pars_r);
}

/*********************************************************/

void Free_Tree_Lk(arbre *tree)
{
  int i;
  edge *b;
  node *n;

  b = NULL;
  n = NULL;

  For(i,3) Free(tree->log_lks_aLRT[i]);
  Free(tree->log_lks_aLRT);

  Free(tree->c_lnL_sorted);
  Free(tree->site_lk);

  For(i,tree->mod->n_catg) Free(tree->log_site_lk_cat[i]);
  Free(tree->log_site_lk_cat);
				
  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];      
      Free_Edge_Lk(tree,b);
    }
}


/*********************************************************/


void Free_Node_Lk(node *n)
{
/*   Free(n->n_ex_nodes); */
}

/*********************************************************/

void Free_Edge_Lk(arbre *tree, edge *b)
{
  int i,j;

  Free(b->nni);

  Free(b->div_post_pred_left);
  Free(b->div_post_pred_rght);

  if(b->p_lk_left)
    {
      For(i,tree->data->crunch_len)
	{
	  For(j,tree->mod->n_catg)
	    {
	      Free(b->p_lk_left[i][j]);
	    }
	  Free(b->p_lk_left[i]);
	}
      Free(b->p_lk_left);
      if(b->sum_scale_f_left) Free(b->sum_scale_f_left);
    }

  if(b->p_lk_tip_l)
    {
      For(i,tree->data->crunch_len)
	{
	  Free(b->p_lk_tip_l[i]);
	}
      Free(b->p_lk_tip_l);
    }


  if(b->p_lk_rght)
    {
      For(i,tree->data->crunch_len)
	{
	  For(j,tree->mod->n_catg)
	    {
	      Free(b->p_lk_rght[i][j]);
	    }
	  Free(b->p_lk_rght[i]);
	}
      Free(b->p_lk_rght);
      if(b->sum_scale_f_rght) Free(b->sum_scale_f_rght);
    }

  if(b->p_lk_tip_r)
    {
      For(i,tree->data->crunch_len)
	{
	  Free(b->p_lk_tip_r[i]);
	}
      Free(b->p_lk_tip_r);
    }


  For(i,tree->mod->n_catg)
    {
      For(j,tree->mod->ns)
	{
	  Free(b->Pij_rr[i][j]);
	}
      Free(b->Pij_rr[i]);
    }
  Free(b->Pij_rr);
}

/*********************************************************/

void Free_Model(model *mod)
{
  int i,j;

  Free(mod->modelname);
  Free(mod->custom_mod_string);
  Free(mod->user_b_freq);
  Free(mod->rr_num);
  Free(mod->rr);
  Free(mod->rr_val);
  Free(mod->n_rr_per_cat);
  Free(mod->s_opt);
  Free(mod->pi);
  Free(mod->pi_unscaled);
  Free(mod->gamma_r_proba);
  Free(mod->gamma_rr);
  Free(mod->qmat);
  Free(mod->qmat_buff);
  Free_Eigen(mod->eigen);

  For(i,mod->n_catg)
    {
      For(j,mod->ns) Free(mod->Pij_rr[i][j]);
      Free(mod->Pij_rr[i]);
    }
  Free(mod->Pij_rr);

  if(mod->n_rr_branch) 
    {
      Free(mod->rr_branch);
      Free(mod->p_rr_branch);
    }

  #ifndef PHYML
  /* Commented out by wash, we don't use M4 mode, setting PHYML does
     not help */
  /*  M4_Free_M4_Model(mod->m4mod); */
  #endif 
  Free(mod);
}

/*********************************************************/

void Free(void *p)
{
  free(p);
}

/*********************************************************/

void Free_Input(option *io)
{
  Free(io->in_seq_file);
  Free(io->in_tree_file);
  Free(io->out_tree_file);
  Free(io->out_best_tree_file);
  Free(io->out_boot_tree_file);
  Free(io->out_boot_stats_file);
  Free(io->out_stats_file);
  Free(io->out_lk_file); 
  Free(io->out_ps_file);
  Free(io->out_trace_file);
  Free(io->nt_or_cd);
  Free(io);
}

/*********************************************************/

void Free_St(superarbre *st)
{
  int i;

  For(i,2*st->tree->n_otu-3) 
    Free(st->tree->t_edges[i]->nni);

  For(i,st->n_gt) Free(st->match_st_node_in_gt[i]);

  Free(st->match_st_node_in_gt);

  Free_Tree(st->tree);
  
  Free(st);
}

/*********************************************************/

void Free_Eigen(eigen *eigen_struct)
{
  Free(eigen_struct->space_int);
  Free(eigen_struct->space);
  Free(eigen_struct->e_val);
  Free(eigen_struct->e_val_im);
  Free(eigen_struct->r_e_vect);
  Free(eigen_struct->r_e_vect_im);
  Free(eigen_struct->l_e_vect);
  Free(eigen_struct->q);
  Free(eigen_struct);
}

/*********************************************************/

void Free_One_Spr(spr *this_spr)
{
  Free(this_spr->path);
  Free(this_spr);
}

/*********************************************************/

void Free_Spr_List(arbre *tree)
{
  int i;

  For(i,tree->size_spr_list+1)
    {
      Free_One_Spr(tree->spr_list[i]);
    }
  Free(tree->spr_list);

}

/*********************************************************/


void Free_Triplet(ttriplet *t)
{
  int i,j,k;

  Free(t->F_bc);
  Free(t->F_cd);
  Free(t->F_bd);
  Free(t->pi_bc);
  Free(t->pi_cd);
  Free(t->pi_bd);

  For(k,t->mod->n_catg) 
    {
      For(i,t->size) 
	{
	  For(j,t->size) Free(t->core[k][i][j]);  
	  Free(t->core[k][i]);
	}
      Free(t->core[k]);	  
    }
  Free(t->core);

  For(i,t->size) 
    {
      For(j,t->size) Free(t->p_one_site[i][j]);  
      Free(t->p_one_site[i]);
    }
  Free(t->p_one_site);

  For(i,t->size) 
    {
      For(j,t->size) Free(t->sum_p_one_site[i][j]);  
      Free(t->sum_p_one_site[i]);
    }
  Free(t->sum_p_one_site);

  Free_Eigen(t->eigen_struct);
  
  Free(t);
}

/*********************************************************/

void Free_Actual_CSeq(allseq *data)
{
  int i;
  For(i,data->n_otu)
    {
      Free(data->c_seq[i]->state);
      data->c_seq[i]->state = NULL;
    }
}

/*********************************************************/

void Free_Prefix_Tree(pnode *n, int size)
{
  int i;
  
  For(i,size)
    {
      if(n->next[i])
	{
	  Free_Prefix_Tree(n->next[i],size);
	}
    }
  Free_Pnode(n);
}

/*********************************************************/

void Free_Pnode(pnode *n)
{
  Free(n->next);
  Free(n);
}

/*********************************************************/

void Free_Rates(int n_otu, trate *r)
{
  int i;

  Free(r->br_r);
  Free(r->lexp);

  For(i,n_otu) Free(r->mc_mr[i]);
  Free(r->mc_mr);

  Free(r);
}

/*********************************************************/
/*********************************************************/
/*********************************************************/

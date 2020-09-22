/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "mixt.h"

int n_sec1 = 0;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_All(t_tree *mixt_tree)
{
  t_tree *curr, *next;
  int i;
  
  curr = mixt_tree;
  next = mixt_tree->next;

  do
    {
      MIXT_Chain_String(curr->mod->aa_rate_mat_file,next->mod->aa_rate_mat_file);
      MIXT_Chain_String(curr->mod->modelname,next->mod->modelname);
      MIXT_Chain_String(curr->mod->custom_mod_string,next->mod->custom_mod_string);
      MIXT_Chain_Scalar_Dbl(curr->mod->kappa,next->mod->kappa);
      MIXT_Chain_Scalar_Dbl(curr->mod->lambda,next->mod->lambda);
      MIXT_Chain_Scalar_Dbl(curr->mod->br_len_mult,next->mod->br_len_mult);
      MIXT_Chain_Scalar_Dbl(curr->mod->br_len_mult_unscaled,next->mod->br_len_mult_unscaled);
      MIXT_Chain_Scalar_Dbl(curr->mod->mr,next->mod->mr);
      MIXT_Chain_Vector_Dbl(curr->mod->Pij_rr,next->mod->Pij_rr);
      MIXT_Chain_Vector_Dbl(curr->mod->e_frq->user_b_freq,next->mod->e_frq->user_b_freq);
      for(i=0;i<2*mixt_tree->n_otu-1;++i) MIXT_Chain_Scalar_Dbl(curr->a_edges[i]->l,next->a_edges[i]->l);
      for(i=0;i<2*mixt_tree->n_otu-1;++i) MIXT_Chain_Scalar_Dbl(curr->a_edges[i]->l_old,next->a_edges[i]->l_old);
      for(i=0;i<2*mixt_tree->n_otu-1;++i) MIXT_Chain_Scalar_Dbl(curr->a_edges[i]->l_var,next->a_edges[i]->l_var);
      for(i=0;i<2*mixt_tree->n_otu-1;++i) MIXT_Chain_Scalar_Dbl(curr->a_edges[i]->l_var_old,next->a_edges[i]->l_var_old);
      MIXT_Chain_Rmat(curr->mod->r_mat,next->mod->r_mat);
      MIXT_Chain_RAS(curr->mod->ras,next->mod->ras);
      MIXT_Chain_Efrq(curr->mod->e_frq,next->mod->e_frq);
      MIXT_Chain_Eigen(curr->mod->eigen,next->mod->eigen);

      curr = next;
      next = next->next;
    }
  while(next);

  MIXT_Chain_Edges(mixt_tree);
  MIXT_Chain_Nodes(mixt_tree);
  MIXT_Chain_Sprs(mixt_tree);
  Make_Rmat_Weight(mixt_tree);
  Make_Efrq_Weight(mixt_tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_Edges(t_tree *mixt_tree)
{
  int i;
  t_edge *b;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      for(i=0;i<2*tree->n_otu-1;++i)
        {
          b = tree->a_edges[i];
          
          if(tree->next)      b->next       = tree->next->a_edges[i];
          if(tree->prev)      b->prev       = tree->prev->a_edges[i];
          if(tree->next_mixt) b->next_mixt  = tree->next_mixt->a_edges[i];
          if(tree->prev_mixt) b->prev_mixt  = tree->prev_mixt->a_edges[i];
        }
      
      if(tree->e_root != NULL)
        {
          b = tree->e_root;
          
          if(tree->next)      b->next       = tree->next->e_root;
          if(tree->prev)      b->prev       = tree->prev->e_root;
          if(tree->next_mixt) b->next_mixt  = tree->next_mixt->e_root;
          if(tree->prev_mixt) b->prev_mixt  = tree->prev_mixt->e_root;
        }
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_Nodes(t_tree *mixt_tree)
{
  int i;
  t_node *n;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      for(i=0;i<2*tree->n_otu-2;++i)
        {
          n = tree->a_nodes[i];
          
          if(tree->next)      n->next       = tree->next->a_nodes[i];
          if(tree->prev)      n->prev       = tree->prev->a_nodes[i];
          if(tree->next_mixt) n->next_mixt  = tree->next_mixt->a_nodes[i];
          if(tree->prev_mixt) n->prev_mixt  = tree->prev_mixt->a_nodes[i];
        }
      
      if(tree->n_root != NULL)
        {
          n = tree->n_root;
          if(tree->next)      n->next       = tree->next->n_root;
          if(tree->prev)      n->prev       = tree->prev->n_root;
          if(tree->next_mixt) n->next_mixt  = tree->next_mixt->n_root;
          if(tree->prev_mixt) n->prev_mixt  = tree->prev_mixt->n_root;
        }
      tree = tree->next;
    }
  while(tree);
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_Sprs(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      if(tree->next)      tree->best_spr->next      = tree->next->best_spr;
      if(tree->prev)      tree->best_spr->prev      = tree->prev->best_spr;
      if(tree->next_mixt) tree->best_spr->next_mixt = tree->next_mixt->best_spr;
      if(tree->prev_mixt) tree->best_spr->prev_mixt = tree->prev_mixt->best_spr;
      
      for(i=0;i<2*tree->n_otu-2;++i)
        {
          if(tree->next)      tree->spr_list_one_edge[i]->next      = tree->next->spr_list_one_edge[i];
          if(tree->prev)      tree->spr_list_one_edge[i]->prev      = tree->prev->spr_list_one_edge[i];
          if(tree->next_mixt) tree->spr_list_one_edge[i]->next_mixt = tree->next_mixt->spr_list_one_edge[i];
          if(tree->prev_mixt) tree->spr_list_one_edge[i]->prev_mixt = tree->prev_mixt->spr_list_one_edge[i];

          if(tree->next)      tree->spr_list_all_edge[i]->next      = tree->next->spr_list_all_edge[i];
          if(tree->prev)      tree->spr_list_all_edge[i]->prev      = tree->prev->spr_list_all_edge[i];
          if(tree->next_mixt) tree->spr_list_all_edge[i]->next_mixt = tree->next_mixt->spr_list_all_edge[i];
          if(tree->prev_mixt) tree->spr_list_all_edge[i]->prev_mixt = tree->prev_mixt->spr_list_all_edge[i];

        }
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_String(t_string *curr, t_string *next)
{
  if(!next)
    {
      return;
    }
  else
    {
      t_string *buff,*last;

      last = NULL;

      /*! Search backward */
      buff = curr;
      while(buff)
        {
          if(buff == next) break;
          buff = buff->prev;
        }

      /*! Search forward */
      if(!buff)
        {
          buff = curr;
          while(buff)
            {
              if(buff == next) break;
              buff = buff->next;
            }
        }


      if(!buff)
        {
          last = curr;
          while(last->next) { last = last->next; }

          last->next = next;
          next->prev = last;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_Vector_Dbl(vect_dbl *curr, vect_dbl *next)
{
  if(!next)
    {
      return;
    }
  else
    {
      vect_dbl *buff,*last;

      last = NULL;

      buff = curr;
      while(buff)
        {
          if(buff == next) break;
          buff = buff->prev;
        }

      /*! Search forward */
      if(!buff)
        {
          buff = curr;
          while(buff)
            {
              if(buff == next) break;
              buff = buff->next;
            }
        }

      if(!buff)
        {
          last = curr;
          while(last->next) { last = last->next; }

          last->next = next;
          next->prev = last;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_Scalar_Dbl(scalar_dbl *curr, scalar_dbl *next)
{
  if(!next)
    {
      return;
    }
  else
    {
      scalar_dbl *buff, *last;

      last = NULL;

      /*! Search backward */
      buff = curr;
      while(buff)
        {
          if(buff == next) break;
          buff = buff->prev;
        }

      /*! Search forward */
      if(!buff)
        {
          buff = curr;
          while(buff)
            {
              if(buff == next) break;
              buff = buff->next;
            }
        }

      /*! Not chained yet. Add next at the end of chained list */
      if(!buff)
        {
          last = curr;
          while(last->next) { last = last->next; }

          last->next = next;
          next->prev = last;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_Rmat(t_rmat *curr, t_rmat *next)
{
  if(!next)
    {
      return;
    }
  else
    {
      t_rmat *buff, *last;

      last = NULL;

      buff = curr;
      while(buff)
        {
          if(buff == next) break;
          buff = buff->prev;
        }

      /*! Search forward */
      if(!buff)
        {
          buff = curr;
          while(buff)
            {
              if(buff == next) break;
              buff = buff->next;
            }
        }

      if(!buff)
        {
          last = curr;
          while(last->next) { last = last->next; }

          last->next = next;
          next->prev = last;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_Efrq(t_efrq *curr, t_efrq *next)
{
  if(!next)
    {
      return;
    }
  else
    {
      t_efrq *buff,*last;

      last = NULL;

      buff = curr;
      while(buff)
        {
          if(buff == next) break;
          buff = buff->prev;
        }

      /*! Search forward */
      if(!buff)
        {
          buff = curr;
          while(buff)
            {
              if(buff == next) break;
              buff = buff->next;
            }
        }

      if(!buff)
        {
          last = curr;
          while(last->next) { last = last->next; }

          last->next = next;
          next->prev = last;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_Eigen(eigen *curr, eigen *next)
{
  if(!next)
    {
      return;
    }
  else
    {
      eigen *buff,*last;

      last = NULL;

      buff = curr;
      while(buff)
        {
          if(buff == next) break;
          buff = buff->prev;
        }

      /*! Search forward */
      if(!buff)
        {
          buff = curr;
          while(buff)
            {
              if(buff == next) break;
              buff = buff->next;
            }
        }

      if(!buff)
        {
          last = curr;
          while(last->next) { last = last->next; }

          last->next = next;
          next->prev = last;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_RAS(t_ras *curr, t_ras *next)
{
  if(!next) return;
  else
    {
      t_ras *buff,*last;

      last = NULL;

      buff = curr;
      while(buff)
        {
          if(buff == next) break;
          buff = buff->prev;
        }

      /*! Search forward */
      if(!buff)
        {
          buff = curr;
          while(buff)
            {
              if(buff == next) break;
              buff = buff->next;
            }
        }

      if(!buff)
        {
          last = curr;
          while(last->next) { last = last->next; }

          last->next = next;
          next->prev = last;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_Rates(t_rate *curr, t_rate *next)
{
  if(!next) return;
  else
    {
      t_rate *buff,*last;

      last = NULL;

      buff = curr;
      while(buff)
        {
          if(buff == next) break;
          buff = buff->prev;
        }

      /*! Search forward */
      if(!buff)
        {
          buff = curr;
          while(buff)
            {
              if(buff == next) break;
              buff = buff->next;
            }
        }

      if(!buff)
        {
          last = curr;
          while(last->next) { last = last->next; }

          last->next = next;
          next->prev = last;
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Chain_Cal(t_tree *mixt_tree)
{
  int i;

  for(i=0;i<mixt_tree->rates->n_cal-1;i++) 
    {
      mixt_tree->rates->a_cal[i]->next   = mixt_tree->rates->a_cal[i+1];
      mixt_tree->rates->a_cal[i+1]->prev = mixt_tree->rates->a_cal[i];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Turn_Branches_OnOff_In_All_Elem(int onoff, t_tree *mixt_tree)
{
  t_tree *tree;

  /*! Turn all branches to ON state */
  tree = mixt_tree;
  do
    {
      MIXT_Turn_Branches_OnOff_In_One_Elem(ON,tree);
      tree = tree->next_mixt;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Turn_Branches_OnOff_In_One_Elem(int onoff, t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  if(mixt_tree->is_mixt_tree == NO)
    {
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
      Exit("\n");
    }

  tree = mixt_tree;

  do
    {
      for(i=0;i<2*tree->n_otu-1;++i) tree->a_edges[i]->l->onoff = onoff;
      tree = tree->next;
    }
  while(tree && tree->is_mixt_tree == NO);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Post_Order_Lk(t_node *mixt_a, t_node *mixt_d, t_tree *mixt_tree)
{
  t_tree *tree;
  t_node *a,*d;

  
  tree = mixt_tree;
  a    = mixt_a;
  d    = mixt_d;

  assert(a);
  assert(d);
  assert(tree);

  do
    {
      if(tree->is_mixt_tree)
        {
          tree = tree->next;
          a    = a->next;
          d    = d->next;
        }

      assert(a);
      assert(d);
      assert(tree);

      if(tree->mod->ras->invar == NO) Post_Order_Lk(a,d,tree);

      tree = tree->next;
      a    = a->next;
      d    = d->next;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Pre_Order_Lk(t_node *mixt_a, t_node *mixt_d, t_tree *mixt_tree)
{
  t_tree *tree;
  t_node *a,*d;

  tree = mixt_tree;
  a    = mixt_a;
  d    = mixt_d;

  assert(a);
  assert(d);
  assert(tree);

  do
    {
      if(tree->is_mixt_tree)
        {
          tree = tree->next;
          a    = a->next;
          d    = d->next;
        }


      assert(a);
      assert(d);
      assert(tree);

      if(tree->mod->ras->invar == NO) Pre_Order_Lk(a,d,tree);

      tree = tree->next;
      a    = a->next;
      d    = d->next;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIXT_Lk(t_edge *mixt_b, t_tree *mixt_tree)
{
  t_tree *tree,*cpy_mixt_tree;
  t_edge *b,*cpy_mixt_b;
  phydbl sum_lnL;
  unsigned int site, br, ns;
  phydbl *sum_scale_left_cat,*sum_scale_rght_cat;
  phydbl sum;
  phydbl site_lk_cat,site_lk,log_site_lk,inv_site_lk;
  int num_prec_issue;
  int ambiguity_check,state;
  int l;
  phydbl r_mat_weight_sum, e_frq_weight_sum, sum_probas;
  phydbl len;
  phydbl *expl,*dot_prod;
  
  tree          = NULL;
  b             = NULL;
  expl          = NULL;
  cpy_mixt_tree = mixt_tree;
  cpy_mixt_b    = mixt_b;
  len           = -1.;
  ns            = 0;
  
  /* MIXT_Update_Br_Len_Multipliers(mixt_tree->mod); */

#if (defined PHYTIME || defined INVITEE || defined PHYREX || defined DATE)
  if((mixt_tree->rates) && (mixt_tree->rates->bl_from_rt)) RATES_Update_Cur_Bl(mixt_tree);
#endif

  
  /* Update RAS structure (mixt_tree level) */
  tree = mixt_tree;
  do
    {
      Update_RAS(tree->mod);
      tree = tree->next_mixt;
    }
  while(tree);

  /* Update other model structure (tree level) */
  tree = mixt_tree->next;
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;
      
      if(!cpy_mixt_b)
        {
          if(!Update_Boundaries(tree->mod) ||
             !Update_Efrq(tree->mod) ||
             !Update_Eigen(tree->mod))
            {
              PhyML_Fprintf(stderr,"\n. move: %s",mixt_tree->mcmc->move_name[mixt_tree->mcmc->move_idx]);
              Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
            }
        }
      
      tree = tree->next;
    }
  while(tree);

  
  if(!cpy_mixt_b)
    {
      for(br=0;br<2*mixt_tree->n_otu-3;++br) MIXT_Update_PMat_At_Given_Edge(mixt_tree->a_edges[br],mixt_tree);
      
      if(mixt_tree->n_root && mixt_tree->ignore_root == NO)
        {
          MIXT_Update_PMat_At_Given_Edge(mixt_tree->n_root->b[1],mixt_tree);
          MIXT_Update_PMat_At_Given_Edge(mixt_tree->n_root->b[2],mixt_tree);
        }
    }
  else
    {
      MIXT_Update_PMat_At_Given_Edge(mixt_b,mixt_tree);
    }

  if(!cpy_mixt_b)
    {
      if(mixt_tree->n_root != NULL)
        {
          if(mixt_tree->ignore_root == NO)
            {                  
              MIXT_Post_Order_Lk(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree);
              MIXT_Post_Order_Lk(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree);
              
              MIXT_Update_Partial_Lk(mixt_tree,mixt_tree->n_root->b[1],mixt_tree->n_root);
              MIXT_Update_Partial_Lk(mixt_tree,mixt_tree->n_root->b[2],mixt_tree->n_root);
              
              if(mixt_tree->both_sides == YES)
                {
                  MIXT_Pre_Order_Lk(mixt_tree->n_root,mixt_tree->n_root->v[1],mixt_tree);              
                  MIXT_Pre_Order_Lk(mixt_tree->n_root,mixt_tree->n_root->v[2],mixt_tree);
                }
            }
          else
            {
              MIXT_Post_Order_Lk(mixt_tree->e_root->rght,mixt_tree->e_root->left,mixt_tree);              
              MIXT_Post_Order_Lk(mixt_tree->e_root->left,mixt_tree->e_root->rght,mixt_tree);
              
              if(mixt_tree->both_sides == YES)
                {
                  MIXT_Pre_Order_Lk(mixt_tree->e_root->rght,mixt_tree->e_root->left,mixt_tree);
                  MIXT_Pre_Order_Lk(mixt_tree->e_root->left,mixt_tree->e_root->rght,mixt_tree);
                }
            }
        }
      else
        {
          MIXT_Post_Order_Lk(mixt_tree->a_nodes[0],mixt_tree->a_nodes[0]->v[0],mixt_tree);
          if(mixt_tree->both_sides == YES)
            {
              MIXT_Pre_Order_Lk(mixt_tree->a_nodes[0],mixt_tree->a_nodes[0]->v[0],mixt_tree);
            }
        }
    }
  
  do /*! Consider each element of the data partition */
    {
      mixt_tree->c_lnL = 0.0;
      tree = mixt_tree->next;
      do
        {
          tree->c_lnL = 0.0;
          tree = tree->next;
        }
      while(tree && tree->is_mixt_tree == NO);
            
      Set_Br_Len_Var(mixt_b,mixt_tree);

      if(!cpy_mixt_b) 
        {
          if(mixt_tree->n_root) 
            {
              if(mixt_tree->ignore_root == NO)
                mixt_b = (mixt_tree->n_root->v[1]->tax == NO)?(mixt_tree->n_root->b[2]):(mixt_tree->n_root->b[1]);
              else
                mixt_b = mixt_tree->e_root;
            }          
          else                  
            {
              mixt_b = mixt_tree->a_nodes[0]->b[0];
            }
        }
      
      sum_scale_left_cat = (phydbl *)mCalloc(MAX(mixt_tree->mod->ras->n_catg,mixt_tree->mod->n_mixt_classes),sizeof(phydbl));
      sum_scale_rght_cat = (phydbl *)mCalloc(MAX(mixt_tree->mod->ras->n_catg,mixt_tree->mod->n_mixt_classes),sizeof(phydbl));

      r_mat_weight_sum = MIXT_Get_Sum_Chained_Scalar_Dbl(mixt_tree->next->mod->r_mat_weight);
      e_frq_weight_sum = MIXT_Get_Sum_Chained_Scalar_Dbl(mixt_tree->next->mod->e_frq_weight);
      sum_probas       = MIXT_Get_Sum_Of_Probas_Across_Mixtures(r_mat_weight_sum,e_frq_weight_sum,mixt_tree);


      b    = mixt_b->next;
      tree = mixt_tree->next;
      do
        {
          while(tree->mod->ras->invar == YES)
            {
              tree = tree->next;
              b    = b->next;

              if(tree == NULL || tree->is_mixt_tree == YES)
                {
                  PhyML_Fprintf(stderr,"\n. %p",(void *)tree);
                  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
                  Exit("\n");
                }
            }

          if(tree->update_eigen_lr == YES) Update_Eigen_Lr(b,tree);
          
          if(tree->use_eigen_lr == YES)
            {
              len = b->l->v;
              len *= tree->mod->br_len_mult->v;
              len *= tree->mixt_tree->mod->ras->gamma_rr->v[tree->mod->ras->parent_class_number];
              if(len < tree->mod->l_min)      len = tree->mod->l_min;
              else if(len > tree->mod->l_max) len = tree->mod->l_max;
              expl = tree->expl;
              for(l=0;l<tree->mod->ns;l++) expl[l] = (phydbl)POW(tree->mod->eigen->e_val[l],len);
            }

          tree = tree->next;
          b    = b->next;

        }
      while(tree && tree->is_mixt_tree == NO);

      tree = mixt_tree;
      do
        {
          tree->numerical_warning = NO;
          tree = tree->next;
        }
      while(tree);

      mixt_tree->c_lnL = .0;
      
      for(site=0;site<mixt_tree->n_pattern;++site)
        {
          b    = mixt_b->next;
          tree = mixt_tree->next;
          
          /*! Skip calculations if model has zero rate */
          while(tree->mod->ras->invar == YES)
            {
              tree = tree->next;
              b    = b->next;

              if(tree == NULL || tree->is_mixt_tree == YES)
                {
                  PhyML_Fprintf(stderr,"\n. %p",(void *)tree);
                  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d (function '%s') \n",__FILE__,__LINE__,__FUNCTION__);
                  Exit("\n");
                }
            }

          ambiguity_check = -1;
          state           = -1;

          if((b->rght->tax)                   &&
             (!mixt_tree->mod->s_opt->greedy) &&
             (mixt_tree->data->wght[site] > SMALL))
            {
              ambiguity_check = b->rght->c_seq->is_ambigu[site];
              if(!ambiguity_check) state = b->rght->c_seq->d_state[site];
            }
          
          /*! For all classes in the mixture */
          do
            {
              if(tree->is_mixt_tree == YES)
                {
                  tree = tree->next;
                  b    = b->next;
                }

              ns                     = tree->mod->ns;
              tree->curr_site        = site;
              tree->apply_lk_scaling = NO;
              site_lk_cat            = 0.0;
              expl                   = tree->expl;
              dot_prod               = tree->dot_prod;
              
              if(!(tree->mod->ras->invar == YES && mixt_tree->is_mixt_tree == YES) && (tree->data->wght[tree->curr_site] > SMALL))
                {
                  if((b->rght->tax) && (tree->mod->s_opt->greedy == NO))
                    {
                      ambiguity_check = b->rght->c_seq->is_ambigu[tree->curr_site];
                      if(ambiguity_check == NO) state = b->rght->c_seq->d_state[tree->curr_site];                        
                    }
                  
                  if(tree->use_eigen_lr == YES)
                    {
                      site_lk_cat = Lk_Core_Eigen_Lr(expl,dot_prod + site*ns,b,tree);
                    }
                  else
                    {
                      if(b->rght->tax == YES)
                        site_lk_cat = Lk_Core(state,ambiguity_check,
                                              b->p_lk_left + site*ns,
                                              b->p_lk_tip_r + site*ns,
                                              b->Pij_rr,b->tPij_rr,b,tree);
                      else
                        site_lk_cat = Lk_Core(state,ambiguity_check,
                                              b->p_lk_left + site*ns,
                                              b->p_lk_rght  + site*ns,
                                              b->Pij_rr,b->tPij_rr,b,tree);
                    }
                }
              
              
              tree->apply_lk_scaling = YES;

              sum_scale_left_cat[tree->mod->ras->parent_class_number] =
                (b->sum_scale_left)?
                (b->sum_scale_left[site]):
                (0.0);

              sum_scale_rght_cat[tree->mod->ras->parent_class_number] =
                (b->sum_scale_rght)?
                (b->sum_scale_rght[site]):
                (0.0);

              sum =
                sum_scale_left_cat[tree->mod->ras->parent_class_number] +
                sum_scale_rght_cat[tree->mod->ras->parent_class_number];

              if(sum > 1024.)
                {
                  /* PhyML_Fprintf(stderr,"\n. Numerical precision issue detected (sum = %g)!!!",sum); */
                  sum = 1023.;
                  tree->mixt_tree->numerical_warning = YES;
                  tree->numerical_warning = YES;
                }
              
              tree->unscaled_site_lk_cat[site] = site_lk_cat;
              
              site_lk_cat /= pow(2,sum);
              
              tree->site_lk_cat[0] = site_lk_cat;
                      
              tree = tree->next;
              b    = b->next;          
            }
          while(tree && tree->is_mixt_tree == NO); 
          // done with all trees in the mixture for this partition element.
          // All likelihood values are in site_lk_cat[0]
          

          tree    = mixt_tree->next;
          b       = mixt_b->next;
          site_lk = .0;

          do
            {
              if(tree->mod->ras->invar == YES)
                {
                  tree = tree->next;
                  b    = b->next;
                  if(!(tree && tree->is_mixt_tree == NO)) break;
                }

              site_lk +=
                tree->site_lk_cat[0] *
                mixt_tree->mod->ras->gamma_r_proba->v[tree->mod->ras->parent_class_number] *
                tree->mod->r_mat_weight->v / r_mat_weight_sum *
                tree->mod->e_frq_weight->v / e_frq_weight_sum /
                sum_probas;
                            
              tree = tree->next;
              b    = b->next;
            }
          while(tree && tree->is_mixt_tree == NO);


          /* Scaling for invariants */
          if(mixt_tree->mod->ras->invar == YES)
            {
              num_prec_issue = NO;

              tree = mixt_tree->next;
              while(tree->mod->ras->invar == NO)
                {
                  tree = tree->next;
                  if(tree == NULL || tree->is_mixt_tree == YES)
                    {
                      PhyML_Fprintf(stderr,"\n. tree: %p",tree);
                      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d",__FILE__,__LINE__);
                      Exit("\n");
                    }
                }

              /*! 'tree' will give the correct state frequencies (as opposed to mixt_tree) */
              inv_site_lk = Invariant_Lk(0,site,&num_prec_issue,tree);

              if(num_prec_issue == YES) // inv_site_lk >> site_lk
                {
                  site_lk = inv_site_lk * mixt_tree->mod->ras->pinvar->v;
                }
              else
                {
                  site_lk = site_lk * (1. - mixt_tree->mod->ras->pinvar->v) + inv_site_lk * mixt_tree->mod->ras->pinvar->v;
                }
            }

          if(site_lk < SMALL)
            {
              /* PhyML_Fprintf(stderr,"\n. site = %d",site); */
              /* PhyML_Fprintf(stderr,"\n. invar = %d",mixt_tree->data->invar[site]); */
              /* PhyML_Fprintf(stderr,"\n. mixt = %d",mixt_tree->is_mixt_tree); */
              /* PhyML_Fprintf(stderr,"\n. lk = %G log(lk) = %f < %G",site_lk,log_site_lk,-BIG); */
              /* for(class=0;class<mixt_tree->mod->ras->n_catg;class++) PhyML_Fprintf(stderr,"\n. rr=%f p=%f",mixt_tree->mod->ras->gamma_rr->v[class],mixt_tree->mod->ras->gamma_r_proba->v[class]); */
              /* PhyML_Fprintf(stderr,"\n. pinv = %G",mixt_tree->mod->ras->pinvar->v); */
              /* PhyML_Fprintf(stderr,"\n. bl mult = %G",mixt_tree->mod->br_len_mult->v); */
              /* PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d.\n",__FILE__,__LINE__); */
              /* Exit("\n"); */
              
              site_lk = SMALL;
              mixt_tree->numerical_warning = YES;
            }

          log_site_lk = log(site_lk);

          
                    
          // ... or using the log-likelihood
          if(isinf(site_lk) || isnan(site_lk))
            {
              mixt_tree->cur_site_lk[site] = exp(log_site_lk);
            }
          else
            {
              mixt_tree->cur_site_lk[site] = site_lk;
            }
          

          /* Multiply log likelihood by the number of times this site pattern is found in the data */
          mixt_tree->c_lnL_sorted[site] = mixt_tree->data->wght[site]*log_site_lk;

          
          mixt_tree->c_lnL += mixt_tree->data->wght[site]*log_site_lk;

        }

      Free(sum_scale_left_cat);
      Free(sum_scale_rght_cat);

      mixt_tree = mixt_tree->next_mixt;
      mixt_b    = mixt_b->next_mixt;
    }
  while(mixt_tree);

  mixt_tree = cpy_mixt_tree;
  mixt_b    = cpy_mixt_b;

  sum_lnL = .0;
  do
    {
      sum_lnL += mixt_tree->c_lnL;
      mixt_tree = mixt_tree->next_mixt;
    }
  while(mixt_tree);

  mixt_tree = cpy_mixt_tree;
  do
    {
      mixt_tree->c_lnL = sum_lnL;
      mixt_tree = mixt_tree->next_mixt;
    }
  while(mixt_tree);

  mixt_tree = cpy_mixt_tree;

  return mixt_tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Update_Partial_Lk(t_tree *mixt_tree, t_edge *mixt_b, t_node *mixt_d)
{
  t_tree *tree;
  t_edge *b;
  t_node *d;

  tree = mixt_tree;
  b    = mixt_b;
  d    = mixt_d;

  do
    {
      if(tree->is_mixt_tree)
        {
          tree = tree->next;
          b    = b->next;
          d    = d->next;
        }

      if(tree->mod->ras->invar == NO) Update_Partial_Lk(tree,b,d);

      tree = tree->next;
      b    = b->next;
      d    = d->next;

    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Update_Partial_Pars(t_tree *mixt_tree, t_edge *mixt_b, t_node *mixt_d)
{
  t_tree *tree;
  t_edge *b;
  t_node *d;

  tree = mixt_tree;
  b    = mixt_b;
  d    = mixt_d;

  do
    {
      if(tree->is_mixt_tree)
        {
          tree = tree->next;
          b    = b->next;
          d    = d->next;
        }

      Update_Partial_Pars(tree,b,d);

      tree = tree->next;
      b    = b->next;
      d    = d->next;

    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Update_PMat_At_Given_Edge(t_edge *mixt_b, t_tree *mixt_tree)
{
  t_tree *tree;
  t_edge *b;

  tree = mixt_tree;
  b    = mixt_b;

  do
    {
      if(tree->is_mixt_tree)
        {
          tree = tree->next;
          b    = b->next;
        }

      if(tree->mod->ras->invar == NO) Update_PMat_At_Given_Edge(b,tree);

      tree = tree->next;
      b    = b->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int *MIXT_Get_Number_Of_Classes_In_All_Mixtures(t_tree *mixt_tree)
{
  int *n_catg;
  t_tree *tree;
  int class;

  if(mixt_tree->is_mixt_tree == YES)
    {
      n_catg = NULL;
      tree = mixt_tree;
      class = 0;
      do
        {
          if(!class) n_catg = (int *)mCalloc(1,sizeof(int));
          else       n_catg = (int *)realloc(n_catg,(class+1)*sizeof(int));
          
          tree = tree->next;
          
          n_catg[class]=0;
          do
            {
              n_catg[class]++;
              tree = tree->next;
            }
          while(tree && tree->is_mixt_tree == NO);
          
      class++;
    }
  while(tree);
    }
  else
    {
      n_catg = (int *)mCalloc(1,sizeof(int));
      n_catg[0] = mixt_tree->mod->ras->n_catg;
      if(mixt_tree->mod->ras->invar == YES) n_catg[0]++;
    }
  return(n_catg);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree **MIXT_Record_All_Mixtures(t_tree *mixt_tree)
{
  t_tree **tree_list;
  int n_trees;
  t_tree *tree;

  tree_list = NULL;
  n_trees   = 0;
  tree      = mixt_tree;
  do
    {
      if(!tree_list) tree_list = (t_tree **)mCalloc(1,sizeof(t_tree *));
      else           tree_list = (t_tree **)realloc(tree_list,(n_trees+1)*sizeof(t_tree *));

      tree_list[n_trees] = tree;
      n_trees++;
      tree = tree->next;
    }
  while(tree);

  tree_list = (t_tree **)realloc(tree_list,(n_trees+1)*sizeof(t_tree *));
  tree_list[n_trees] = NULL;

  return(tree_list);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Break_All_Mixtures(int *c_max, t_tree *mixt_tree)
{
  t_tree *tree;
  int c,i,n;

  if(mixt_tree->is_mixt_tree == NO) return;

  c = 0;
  n = -1;
  tree = mixt_tree;
  do
    {
      if(tree->is_mixt_tree == YES) 
        {
          c = 0;
          n++;
          tree = tree->next;
        }

      if(c == (c_max[n]-1)           && 
         tree->next != NULL          &&
         tree->next->is_mixt_tree == NO)
        {
          if(tree->mixt_tree->next_mixt == NULL)
            {
              tree->next = NULL;
              for(i=0;i<2*tree->n_otu-1;++i) tree->a_edges[i]->next  = NULL;
              for(i=0;i<2*tree->n_otu-1;++i) tree->a_nodes[i]->next  = NULL;
              for(i=0;i<2*tree->n_otu-2;++i) tree->spr_list_one_edge[i]->next = NULL;
            }
          else
            {
              tree->next = tree->mixt_tree->next_mixt;
              for(i=0;i<2*tree->n_otu-1;++i) tree->a_edges[i]->next  = tree->mixt_tree->next_mixt->a_edges[i];
              for(i=0;i<2*tree->n_otu-1;++i) tree->a_nodes[i]->next  = tree->mixt_tree->next_mixt->a_nodes[i];
              for(i=0;i<2*tree->n_otu-2;++i) tree->spr_list_one_edge[i]->next = tree->mixt_tree->next_mixt->spr_list_one_edge[i];
            }
        }

      tree = tree->next;
      c++;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Reconnect_All_Mixtures(t_tree **tree_list, t_tree *mixt_tree)
{
  t_tree *tree;
  int n_trees;

  if(mixt_tree->is_mixt_tree == NO) return;

  tree = mixt_tree;
  n_trees = 0;
  do
    {
      tree = tree_list[n_trees];
      if(tree->is_mixt_tree == NO) tree->next = tree_list[n_trees+1];
      n_trees++;
      tree = tree->next;
    }
  while(tree);

  MIXT_Chain_All(mixt_tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int *MIXT_Record_Has_Invariants(t_tree *mixt_tree)
{
  int *has_invariants;
  t_tree *tree;
  int n_trees;

  has_invariants = NULL;
  tree = mixt_tree;
  n_trees = 0;
  do
    {
      if(!n_trees) has_invariants = (int *)mCalloc(1,sizeof(int));
      else         has_invariants = (int *)realloc(has_invariants,(n_trees+1)*sizeof(int));
      has_invariants[n_trees] = (tree->mod->ras->invar == YES)?1:0;
      n_trees++;
      tree = tree->next;
    }
  while(tree);

  return(has_invariants);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Reset_Has_Invariants(int *has_invariants, t_tree *mixt_tree)
{
  t_tree *tree;
  int n_trees;

  tree = mixt_tree;
  n_trees = 0;
  do
    {
      tree->mod->ras->invar = has_invariants[n_trees];
      n_trees++;
      tree = tree->next;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Check_Invar_Struct_In_Each_Partition_Elem(t_tree *mixt_tree)
{
  if(mixt_tree->is_mixt_tree == NO) return;
  else
    {
      t_tree *tree;
      int n_inv;

      n_inv = 0;
      tree = mixt_tree;
      do
        {
          if(tree->is_mixt_tree)
            {
              tree  = tree->next;
              n_inv = 0;
            }

          if(tree->mod->ras->invar == YES) n_inv++;

          if(n_inv > 1)
            {
              PhyML_Fprintf(stderr,"\n. Found %d classes of the mixture for file '%s' set to",n_inv,tree->mixt_tree->io->in_align_file);
              PhyML_Fprintf(stderr,"\n. invariable. Only one such class per mixture is allowed.");
              PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n\n",__FILE__,__LINE__);
              Warn_And_Exit("\n");
            }

          if(tree->mixt_tree->mod->ras->invar == NO &&
             tree->mod->ras->invar == YES)
            {
              PhyML_Fprintf(stderr,"\n. Unexpected settings for 'siterates' in a partition element (file '%s')",tree->mixt_tree->io->in_align_file);
              PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n\n",__FILE__,__LINE__);
              Warn_And_Exit("\n");
            }

          tree = tree->next;
        }
      while(tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Check_RAS_Struct_In_Each_Partition_Elem(t_tree *mixt_tree)
{
  if(mixt_tree->is_mixt_tree == NO) return;
  else
    {
      t_tree *tree;
      int n_classes;

      n_classes = 0;
      tree = mixt_tree;
      do
        {
          if(tree->is_mixt_tree)
            {
              if(tree->mod->ras->invar == YES)
                {
                  if(tree->next->mod->ras->invar == NO)
                    {
                      PhyML_Fprintf(stderr,"\n. The invariant site class has to be the first element in");
                      PhyML_Fprintf(stderr,"\n. each <mixtureelem> component. Please amend you XML");
                      PhyML_Fprintf(stderr,"\n. file accordingly.\n");
                      Exit("\n.");
                    }
                }
              tree  = tree->next;
              n_classes = 0;
            }

          if(tree && tree->mod->ras->invar == NO) n_classes++;

          if((tree->next && tree->next->is_mixt_tree == YES) || (!tree->next)) /*! current tree is the last element of this mixture */
            {
              if(n_classes < tree->mixt_tree->mod->ras->n_catg)
                {
                  PhyML_Fprintf(stderr,"\n. %d class%s found in 'partitionelem' for file '%s' while",
                                n_classes,
                                (n_classes>1)?"es\0":"\0",
                                tree->mixt_tree->io->in_align_file);
                  PhyML_Fprintf(stderr,"\n. the corresponding 'siterates' element defined %d class%s.",
                                tree->mixt_tree->mod->ras->n_catg,
                                (tree->mixt_tree->mod->ras->n_catg>1)?"es\0":"\0");
                  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n\n",__FILE__,__LINE__);
                  Warn_And_Exit("\n");
                }
            }

          tree = tree->next;
        }
      while(tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Prune_Subtree(t_node *mixt_a, t_node *mixt_d, t_edge **mixt_target, t_edge **mixt_residual, t_tree *mixt_tree)
{
  t_node *a,*d;
  t_edge *target, *residual;
  t_tree *tree;

  MIXT_Turn_Branches_OnOff_In_One_Elem(OFF,mixt_tree);

  tree     = mixt_tree;
  a        = mixt_a;
  d        = mixt_d;  
  target   = (mixt_target) ? (*mixt_target) : NULL;
  residual = (mixt_residual) ? (*mixt_residual) : NULL;

  do
    {
      if(tree->is_mixt_tree == YES)
        {
          tree     = tree->next;
          a        = a->next;
          d        = d->next;
          target   = target ? target->next : NULL;
          residual = residual ? residual->next : NULL;
        }

      Prune_Subtree(a,d,&target,&residual,tree);

      tree     = tree->next;
      a        = a->next;
      d        = d->next;
      target   = target ? target->next : NULL;
      residual = residual ? residual->next : NULL;
    }
  while(tree && tree->is_mixt_tree == NO);

  if(tree) Prune_Subtree(a,d,&target,&residual,tree);

  /*! Turn branches of this mixt_tree to ON after recursive call
    to Prune_Subtree such that, if branches of mixt_tree->next
    point to those of mixt_tree, they are set to OFF when calling
    Prune */
  MIXT_Turn_Branches_OnOff_In_One_Elem(ON,mixt_tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Graft_Subtree(t_edge *mixt_target, t_node *mixt_link, t_node *mixt_link_daughter, t_edge *mixt_residual, t_node *mixt_target_nd, t_tree *mixt_tree)
{
  t_edge *target,*residual;
  t_node *link,*link_daughter,*target_nd;
  t_tree *tree;

  MIXT_Turn_Branches_OnOff_In_One_Elem(OFF,mixt_tree);

  tree          = mixt_tree;
  target        = mixt_target;
  residual      = mixt_residual;
  link          = mixt_link;
  link_daughter = mixt_link_daughter;
  target_nd     = mixt_target_nd;

  do
    {
      if(tree->is_mixt_tree == YES)
        {
          tree          = tree->next;
          target        = target->next;
          residual      = residual->next;
          link          = link->next;
          link_daughter = link_daughter ? link_daughter->next : NULL;
          target_nd     = target_nd ? target_nd->next : NULL;
        }

      Graft_Subtree(target,link,link_daughter,residual,target_nd,tree);

      tree          = tree->next;
      target        = target->next;
      residual      = residual->next;
      link          = link->next;
      link_daughter = link_daughter ? link_daughter->next : NULL;
      target_nd     = target_nd ? target_nd->next : NULL;
    }
  while(tree && tree->is_mixt_tree == NO);

  if(tree) Graft_Subtree(target,link,link_daughter,residual,target_nd,tree);

  /*! Turn branches of this mixt_tree to ON after recursive call
    to Graft_Subtree such that, if branches of mixt_tree->next
    point to those of mixt_tree, they are set to OFF when calling
    Graft */
  MIXT_Turn_Branches_OnOff_In_One_Elem(ON,mixt_tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Multiply_Scalar_Dbl(scalar_dbl *this, phydbl scalar)
{
  if(this == NULL) return;
  else
    {
      scalar_dbl *buff;
      buff = this;
      do
        {
          buff->v *= scalar;
          buff = buff->next;
        }
      while(buff);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Br_Len_Opt(t_edge *mixt_b, t_tree *mixt_tree)
{
  t_edge *b;
  t_tree *tree;
  scalar_dbl *l;

  MIXT_Turn_Branches_OnOff_In_All_Elem(ON,mixt_tree);

  mixt_tree->ignore_mixt_info = YES;
  
  b = mixt_b;
  tree = mixt_tree;
  l = NULL;
  
  do
    {
      while(tree->is_mixt_tree == YES)
        {
          tree = tree->next;
          b = b->next;
        }

      l = b->l;
      do
        {
          if(l->onoff == ON)
            {
              Br_Len_Opt(&(l->v),mixt_b,mixt_tree);      
              l->onoff = OFF;
            }
          l = l->next;
        }
      while(l);
      
      tree = tree->next;
      b    = b->next;
    }
  while(tree);

  mixt_tree->ignore_mixt_info = NO;
}


/*////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////*/

void MIXT_Br_Len_Involving_Invar(t_tree *mixt_tree)
{
  int i;
  scalar_dbl *l;

  for(i=0;i<2*mixt_tree->n_otu-1;++i) 
    {
      l = mixt_tree->a_edges[i]->l;
      do
        {
          l->v *= (1.-mixt_tree->mod->ras->pinvar->v);
          l = l->next;
        }
      while(l);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Br_Len_Not_Involving_Invar(t_tree *mixt_tree)
{
  int i;
  scalar_dbl *l;

  for(i=0;i<2*mixt_tree->n_otu-1;++i) 
    {
      l = mixt_tree->a_edges[i]->l;
      do
        {
          l->v /= (1.-mixt_tree->mod->ras->pinvar->v);
          l = l->next;
        }
      while(l);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIXT_Unscale_Br_Len_Multiplier_Tree(t_tree *mixt_tree)
{
  int i;
  scalar_dbl *l;

  for(i=0;i<2*mixt_tree->n_otu-1;++i) 
    {
      l = mixt_tree->a_edges[i]->l;
      do
        {
          l->v /= mixt_tree->mod->br_len_mult->v;
          l = l->next;
        }
      while(l);
    }
  return(-1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIXT_Rescale_Br_Len_Multiplier_Tree(t_tree *mixt_tree)
{
  int i;
  scalar_dbl *l;

  for(i=0;i<2*mixt_tree->n_otu-1;++i) 
    {
      l = mixt_tree->a_edges[i]->l;
      do
        {
          l->v *= mixt_tree->mod->br_len_mult->v;
          l = l->next;
        }
      while(l);
    }
  return(-1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIXT_Rescale_Free_Rate_Tree(t_tree *mixt_tree)
{
  int i,side_effect,at_boundary;
  t_edge *b;

  side_effect = NO;
  for(i=0;i<2*mixt_tree->n_otu-1;++i)
    {
      b = mixt_tree->a_edges[i]->next;

      at_boundary = NO;
      if(b->l->v > mixt_tree->mod->l_max-1.E-100 && b->l->v < mixt_tree->mod->l_max+1.E-100) at_boundary = YES;
      if(b->l->v > mixt_tree->mod->l_min-1.E-100 && b->l->v < mixt_tree->mod->l_min+1.E-100) at_boundary = YES;
      
      b->l->v *= mixt_tree->mod->ras->free_rate_mr->v;
      
      if(b->l->v > mixt_tree->mod->l_max && at_boundary == NO) side_effect = YES;
      if(b->l->v < mixt_tree->mod->l_min && at_boundary == NO) side_effect = YES;
    }
  
  return side_effect;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Ignore_Root(int yesno, t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      tree->ignore_root = yesno;
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Alias_Subpatt(int onoff, t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      tree->update_alias_subpatt = onoff;
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Check_Edge_Lens_In_All_Elem(t_tree *mixt_tree)
{
  t_tree *tree;

  /*! Check that all the edges in a mixt_tree at pointing
    to a single set of lengths
  */
  tree = mixt_tree;
  do
    {
      MIXT_Check_Edge_Lens_In_One_Elem(tree);
      tree = tree->next_mixt;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Check_Edge_Lens_In_One_Elem(t_tree *mixt_tree)
{
  t_tree *tree;
  int i;

  tree = mixt_tree->next;
  do
    {
      if(tree->next && tree->next->is_mixt_tree == NO)
        {
          for(i=0;i<2*tree->n_otu-1;++i)
            {
              if(tree->a_edges[i]->l != tree->next->a_edges[i]->l)
                {
                  PhyML_Fprintf(stderr,"\n. %p %p",tree->a_edges[i]->l,tree->next->a_edges[i]->l);
                  PhyML_Fprintf(stderr,"\n. Only one set of edge lengths is allowed ");
                  PhyML_Fprintf(stderr,"\n. in a 'partitionelem'. Please fix your XML file.");
                  Exit("\n");
                }
            }
        }
      tree = tree->next;
    }
  while(tree && tree->next && tree->next->is_mixt_tree == NO);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int MIXT_Pars(t_edge *mixt_b, t_tree *mixt_tree)
{
  t_edge *b;
  t_tree *tree;

  b                 = mixt_b;
  tree              = mixt_tree;
  mixt_tree->c_pars = 0;

  do
    {
      if(tree->next)
        {
          Pars(b?b->next:NULL,tree->next);
          mixt_tree->c_pars += tree->next->c_pars;      
        }
      
      if(mixt_b != NULL) b = b->next_mixt;
      tree = tree->next_mixt;
    }
  while(tree);

  tree = mixt_tree;
  do
    {
      tree->c_pars = mixt_tree->c_pars;
      tree = tree->next_mixt;
    }
  while(tree);

  return(mixt_tree->c_pars);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Bootstrap(char *best_tree, xml_node *root)
{
  xml_node *n,*p_elem;
  char *bootstrap;

  assert(root);
  
  n = XML_Search_Node_Name("phyml",NO,root);

  bootstrap = XML_Get_Attribute_Value(n,"bootstrap");

  if(!bootstrap) return;
  else
    {
      int n_boot,i,j,k;
      xml_attr *boot_attr,*seqfile_attr,*out_attr,*boot_out_attr;
      char *orig_align,*boot_out_file_name,*xml_boot_file_name,*buff;
      FILE *boot_fp_in_align,*xml_boot_file_fp;
      option *io,*dum;
      align **boot_data,**orig_data;
      int position,elem;
      xml_node *boot_root;
      int pid;
      char *s;

      orig_align = (char *)mCalloc(T_MAX_NAME,sizeof(char));

      xml_boot_file_name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(xml_boot_file_name,"phyml_boot_config.");
      pid = (int)getpid();
      sprintf(xml_boot_file_name+strlen(xml_boot_file_name),"%d",pid);
      strcat(xml_boot_file_name,".xml");

      out_attr = XML_Search_Attribute(root,"output.file");
      assert(out_attr);
      boot_out_file_name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(boot_out_file_name,out_attr->value);
      s = XML_Get_Attribute_Value(root,"run.id");
      if(s)
        {
          strcat(boot_out_file_name,"_");
          strcat(boot_out_file_name,s);
        }

      n_boot = atoi(bootstrap);

      io = NULL;
      for(i=0;i<n_boot;++i)
        {
          boot_root = XML_Copy_XML_Graph(root);

          /*! Set the number of bootstrap repeats to 0
            in each generated XML file */
          boot_attr = XML_Search_Attribute(boot_root,"bootstrap");
          assert(boot_attr);
          strcpy(boot_attr->value,"0");

          /*! Set the output file name for each bootstrap analysis */
          boot_out_attr = XML_Search_Attribute(boot_root,"output.file");
          assert(boot_out_attr);
          buff = (char *)mCalloc(T_MAX_NAME,sizeof(char));
          strcpy(buff,boot_out_attr->value);
          Free(boot_out_attr->value);
          boot_out_attr->value = buff;
          sprintf(boot_out_attr->value+strlen(boot_out_attr->value),"_boot.%d",pid);

          p_elem = boot_root;
          elem   = 0;
          do
            {
              p_elem = XML_Search_Node_Name("partitionelem",YES,p_elem);
              if(!p_elem) break;

              io = (option *)Make_Input();
              Set_Defaults_Input(io);

              /*! Get the original sequence file name and the corresponding
                attribute in the XML graph
              */
              seqfile_attr = NULL;
              seqfile_attr = XML_Search_Attribute(p_elem,"file.name");
              assert(seqfile_attr);

              strcpy(orig_align,seqfile_attr->value);

              /*! Open the original sequence file */
              io->fp_in_align = Openfile(orig_align,0);

              /*! Read in the original sequence file */
              orig_data = Get_Seq(io);
              rewind(io->fp_in_align);

              /*! Read in the original sequence file and put
               it in 'boot_data' structure */
              boot_data = Get_Seq(io);

              fclose(io->fp_in_align);

              /*! Bootstrap resampling: sample from original and put in boot */
              for(j=0;j<boot_data[0]->len;++j)
                {
                  position = Rand_Int(0,(int)(boot_data[0]->len-1.0));
                  for(k=0;k<io->n_otu;++k) boot_data[k]->state[j] = orig_data[k]->state[position];
                }

              /*! Modify the sequence file attribute in the original XML
                graph */
              buff = (char *)mCalloc(T_MAX_NAME,sizeof(char));
              Free(seqfile_attr->value);
              seqfile_attr->value = buff;
              sprintf(seqfile_attr->value,"%s_%d_%d",orig_align,elem,i);

              /*! Open a new sequence file with the modified attribute name */
              boot_fp_in_align = Openfile(seqfile_attr->value,1);

              /*! Print the bootstrap data set in it */
              Print_Seq(boot_fp_in_align,boot_data,io->n_otu);
              fclose(boot_fp_in_align);

              Free_Seq(orig_data,io->n_otu);
              Free_Seq(boot_data,io->n_otu);

              Free_Input(io);
              elem++;
            }
          while(p_elem);
          
          /*! Open bootstrap XML file in writing mode */
          xml_boot_file_fp = Openfile(xml_boot_file_name,1);

          /*! Write the bootstrap XML graph */
          XML_Write_XML_Graph(xml_boot_file_fp,boot_root);
          fclose(xml_boot_file_fp);

          /*! Reconstruct the tree */
          dum = PhyML_XML(xml_boot_file_name);
          Free(dum);

          /*! Remove the bootstrap alignment files */
          p_elem = boot_root;
          do
            {
              p_elem = XML_Search_Node_Name("partitionelem",YES,p_elem);
              if(!p_elem) break;
              seqfile_attr = XML_Search_Attribute(p_elem,"file.name");
              unlink(seqfile_attr->value);
            }
          while(p_elem);

          XML_Free_XML_Tree(boot_root);
        }

      Free(xml_boot_file_name);
      Free(orig_align);
      Free(boot_out_file_name);
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Pars_Thresh(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      tree->mod->s_opt->pars_thresh = (tree->io->datatype == AA)?(15):(5);
      tree = tree->next_mixt;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIXT_Get_Mean_Edge_Len(t_edge *mixt_b, t_tree *mixt_tree)
{
  phydbl sum;
  int n;
  t_tree *tree;
  t_edge *b;

  if(mixt_tree->is_mixt_tree == NO) return mixt_b->l->v;
  
  b    = mixt_b;
  tree = mixt_tree;
  sum  = .0;
  n    = 0 ;
  do
    {
      if(tree->is_mixt_tree == YES)
        {
          tree = tree->next;
          b    = b->next;
        }

      sum += b->l->v * (tree->mixt_tree ? tree->mixt_tree->mod->br_len_mult->v : 1.0);
      n++;
      b    = b->next;
      tree = tree->next;
    }
  while(b);

  return(sum / (phydbl)n);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIXT_Get_Sum_Chained_Scalar_Dbl(scalar_dbl *s)
{
  scalar_dbl *s_buff;
  phydbl sum;

  s_buff = s;
  sum = .0;
  do
    {
      sum += s_buff->v;
      s_buff = s_buff->next;
    }
  while(s_buff);

  return sum;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl MIXT_Get_Sum_Of_Probas_Across_Mixtures(phydbl r_mat_weight_sum, phydbl e_frq_weight_sum, t_tree *mixt_tree)
{
  t_tree *tree;
  phydbl sum;

  sum = .0;
  tree = mixt_tree->next;
  do
    {
      // e.g., if mixture has two classes, one of these 
      // corresponding to invariable sites. We need to skip it.
      if(tree->mod->ras->invar == YES) tree = tree->next;
          
      sum +=
        mixt_tree->mod->ras->gamma_r_proba->v[tree->mod->ras->parent_class_number] *
        tree->mod->r_mat_weight->v / r_mat_weight_sum *
        tree->mod->e_frq_weight->v / e_frq_weight_sum;
      

      tree = tree->next;

    }
  while(tree && tree->is_mixt_tree == NO);

  return(sum);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Br_Len_Var(t_edge *mixt_b, t_tree *mixt_tree)
{
  t_tree *tree;
  
  if(mixt_b != NULL)
    {
      t_edge *b;
  
      tree = mixt_tree->next;
      b    = mixt_b->next;

      do
        {
          Set_Br_Len_Var(b,tree);
          tree = tree->next;
          b    = b->next;
        }
      while(tree);
    }
  else
    {
      tree = mixt_tree->next;

      do
        {
          Set_Br_Len_Var(NULL,tree);
          tree = tree->next;
        }
      while(tree);

    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Br_Len(phydbl val, t_edge *mixt_b, t_tree *mixt_tree)
{
  scalar_dbl *l;

  l = mixt_b->l;
  do
    {
      l->v = val;
      l = l->next;
    }
  while(l);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Add_Root(t_edge *mixt_b, t_tree *mixt_tree)
{
  t_tree *tree;
  t_edge *b;

  tree = mixt_tree;
  b    = mixt_b;
  do
    {
      if(tree->is_mixt_tree)
        {
          tree = tree->next;
          b    = b->next;
        }

      // Condition is true when tree is not chained
      if(b == NULL) break; 

      Add_Root(b,tree);

      tree = tree->next;
      b    = b->next;
    }
  while(tree);

  tree = mixt_tree;
  do
    {
      assert(tree->n_root != NULL);
      
      if(tree->next)      tree->n_root->next       = tree->next->n_root;
      if(tree->prev)      tree->n_root->prev       = tree->prev->n_root;
      if(tree->next_mixt) tree->n_root->next_mixt  = tree->next_mixt->n_root;
      if(tree->prev_mixt) tree->n_root->prev_mixt  = tree->prev_mixt->n_root;

      tree = tree->next;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_RATES_Update_Cur_Bl(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree->next;
  do
    {
      RATES_Update_Cur_Bl(tree);
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Update_Br_Len_Multipliers(t_mod *mod)
{
  phydbl sum;
  t_mod *loc;
  int n_mixt;

  loc = mod;
  sum = 0.0;
  n_mixt = 0;
  do
    {
      /* if(loc->s_opt->opt_br_len_mult == YES) */
      /*   { */
          sum += loc->br_len_mult_unscaled->v;
          n_mixt++;
        /* } */
      loc = loc->next_mixt;
    }
  while(loc);

  loc = mod;
  do
    {
      if(loc->s_opt->opt_br_len_mult == YES)
        {
          loc->br_len_mult->v = loc->br_len_mult_unscaled->v / sum;
          loc->br_len_mult->v *= (phydbl)(n_mixt);
          /* printf("\n. HERE %f %f\n",loc->br_len_mult_unscaled->v,loc->br_len_mult->v); */
        }
      loc = loc->next_mixt;
    }
  while(loc);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Init_Model(t_tree *mixt_tree)
{
  t_mod *mod,*mixt_mod;
  option *io;
  t_tree *tree;

  assert(mixt_tree);

  mixt_mod = mixt_tree->mod;
  io       = mixt_tree->io;

  mod = mixt_mod;
  do
    {
      Init_Model(mod->io->cdata,mod,io);
      mod = mod->next;
    }
  while(mod);

  tree = mixt_tree;
  do
    {
      if(tree->next_mixt != NULL)
        {
          tree->mod->next_mixt = tree->next_mixt->mod;
          tree->next_mixt->mod->prev_mixt = tree->mod;
        }
      tree = tree->next_mixt;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Check_Model_Validity(t_tree *mixt_tree)
{
  // Verify that models associated to distinct data partition elements do not 
  // share the same empirical character frequencies.
  
  t_mod *mod_in, *mod_out;

  mod_out = mixt_tree->mod;

  do
    {
      mod_in = mod_out;
      do
        {
          if(mod_in->io->cdata != mod_out->io->cdata)
            {
              if(mod_in->e_frq == mod_out->e_frq)
                {
                  if(mod_in->io->datatype == NT &&
                     mod_in->e_frq->user_state_freq == NO &&
                     mod_in->whichmodel != JC69 &&
                     mod_in->whichmodel != K80)
                    {
                      PhyML_Fprintf(stderr,"\n. A vector of observed nucleotide frequencies should correspond ");
                      PhyML_Fprintf(stderr,"\n. to one data set only. If you are using the XML interface, ");
                      PhyML_Fprintf(stderr,"\n. please amend your file accordingly.");          
                      Exit("\n");
                    }
                  else if(mod_in->io->datatype == AA && mod_in->e_frq->empirical_state_freq == YES)
                    {
                      PhyML_Fprintf(stderr,"\n. A vector of observed amino-acid frequencies should correspond ");
                      PhyML_Fprintf(stderr,"\n. to one data set only. If you are using the XML interface, ");
                      PhyML_Fprintf(stderr,"\n. please amend your file accordingly.");          
                      Exit("\n");
                    }
                }
            }
          mod_in = mod_in->next;
        }
      while(mod_in);      
      mod_out = mod_out->next;
    }
  while(mod_out);


}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *MIXT_Starting_Tree(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = NULL;
  
  if(mixt_tree->io->mod->s_opt->random_input_tree == NO)
    {
      switch(mixt_tree->io->in_tree)
        {
        case 2: // user-defined input tree 
          {
            assert(mixt_tree->io->fp_in_tree);
            
            tree = Read_User_Tree(mixt_tree->io->cdata,
                                  mixt_tree->mod,
                                  mixt_tree->io);
            
            break;
          }
        case 1: case 0:
          {
            // Build a BioNJ tree from the analysis of
            // the first partition element
            tree = Dist_And_BioNJ(mixt_tree->data,
                                  mixt_tree->mod,
                                  mixt_tree->io);
            break;
          }
        default : assert(FALSE);
        }
    }
  else
    {
      tree = Make_Tree_From_Scratch(mixt_tree->n_otu,mixt_tree->data);
      Random_Tree(tree);
    }
  
  return tree;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Connect_Cseqs_To_Nodes(t_tree *mixt_tree)
{
  t_tree *tree;

  Copy_Tree(mixt_tree,mixt_tree->next);
  
  tree = mixt_tree;
  do
    {
      Connect_CSeqs_To_Nodes(tree->data,mixt_tree->io,tree);
      tree = tree->next;
    }
  while(tree);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Init_T_Beg(t_tree *mixt_tree)
{
  t_tree *tree;

  /*! Initialize t_beg in each mixture tree */
  tree = mixt_tree;
  do
    {
      time(&(tree->t_beg));
      tree = tree->next_mixt;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Init_T_End(t_tree *mixt_tree)
{
  t_tree *tree;

  /*! Initialize t_beg in each mixture tree */
  tree = mixt_tree;
  do
    {
      time(&(tree->t_current));
      tree = tree->next_mixt;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Prepare_All(int num_rand_tree, t_tree *mixt_tree)
{
  t_tree *tree;

  MIXT_Check_Model_Validity(mixt_tree);
  MIXT_Init_Model(mixt_tree);
  Print_Data_Structure(NO,stdout,mixt_tree);
  tree = MIXT_Starting_Tree(mixt_tree);
  Copy_Tree(tree,mixt_tree);
  Free_Tree(tree);
  
  if(mixt_tree->io->mod->s_opt->random_input_tree)
    {
      PhyML_Printf("\n\n. [%3d/%3d]",num_rand_tree+1,mixt_tree->io->mod->s_opt->n_rand_starts);
      Random_Tree(mixt_tree);
    }
  
  MIXT_Connect_Cseqs_To_Nodes(mixt_tree);
  MIXT_Init_T_Beg(mixt_tree);
  
  MIXT_Make_Tree_For_Pars(mixt_tree);
  MIXT_Make_Tree_For_Lk(mixt_tree);
  MIXT_Make_Spr(mixt_tree);
  
  MIXT_Chain_All(mixt_tree);
  MIXT_Check_Edge_Lens_In_All_Elem(mixt_tree);
  MIXT_Turn_Branches_OnOff_In_All_Elem(ON,mixt_tree);
  MIXT_Check_Invar_Struct_In_Each_Partition_Elem(mixt_tree);
  MIXT_Check_RAS_Struct_In_Each_Partition_Elem(mixt_tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Make_Spr(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      Make_Spr(tree);
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Make_Tree_For_Pars(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      if(tree->is_mixt_tree == NO) Make_Tree_For_Pars(tree);
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Make_Tree_For_Lk(t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      if(tree->is_mixt_tree == NO) Make_Tree_For_Lk(tree);
      else
        {
          tree->c_lnL_sorted         = (phydbl *)mCalloc(tree->n_pattern,sizeof(phydbl));
          tree->cur_site_lk          = (phydbl *)mCalloc(tree->n_pattern,sizeof(phydbl));
          tree->old_site_lk          = (phydbl *)mCalloc(tree->n_pattern,sizeof(phydbl));
          tree->site_lk_cat          = (phydbl *)mCalloc(MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes),sizeof(phydbl));
          tree->unscaled_site_lk_cat = (phydbl *)mCalloc(MAX(tree->mod->ras->n_catg,tree->mod->n_mixt_classes)*tree->n_pattern,sizeof(phydbl));
          tree->fact_sum_scale       = (int *)mCalloc(tree->n_pattern,sizeof(int));


#if (defined(__AVX__) || defined(__SSE3__))
#ifndef WIN32
          if(posix_memalign((void **)&tree->expl,BYTE_ALIGN,(size_t)2*tree->mod->n_mixt_classes*tree->mod->ns*sizeof(phydbl))) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
#else
          tree->expl = _aligned_malloc(2*tree->mod->n_mixt_classes*tree->mod->ns*sizeof(phydbl),BYTE_ALIGN);
#endif
#else
          tree->expl = (phydbl *)mCalloc(2*tree->mod->n_mixt_classes*tree->mod->ns,sizeof(phydbl));
#endif
          
          for(int i=0;i<2*tree->n_otu-1;++i) Make_Edge_NNI(tree->a_edges[i]);

          tree->log_lks_aLRT = (phydbl **)mCalloc(3,sizeof(phydbl *));
          for(int i=0;i<3;i++) tree->log_lks_aLRT[i] = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));
        }
      
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Ancestral_Sequences_One_Node(t_node *mixt_d, t_tree *mixt_tree, int print)
{
  if(mixt_d->tax) return;
  else
    {
      t_node *v0,*v1,*v2; // three neighbours of d
      t_edge *b0,*b1,*b2;
      int i,j;
      int catg;
      phydbl p0, p1, p2;
      phydbl *p;
      t_node *d,*curr_mixt_d;
      t_tree *tree, *curr_mixt_tree;
      int site,csite;
      phydbl *p_lk0, *p_lk1, *p_lk2;
      int *sum_scale0, *sum_scale1, *sum_scale2;
      phydbl r_mat_weight_sum, e_frq_weight_sum, sum_probas;
      phydbl *Pij0, *Pij1, *Pij2;
      int NsNs, Ns, NsNg;
      FILE *fp;

      if(!mixt_d) return;

      curr_mixt_tree = mixt_tree;
      curr_mixt_d    = mixt_d;
      fp             = mixt_tree->io->fp_out_ancestral_seq;


      do /* For each partition element */
        {
          if(curr_mixt_tree->next)
            {
              r_mat_weight_sum = MIXT_Get_Sum_Chained_Scalar_Dbl(curr_mixt_tree->next->mod->r_mat_weight);
              e_frq_weight_sum = MIXT_Get_Sum_Chained_Scalar_Dbl(curr_mixt_tree->next->mod->e_frq_weight);
              sum_probas       = MIXT_Get_Sum_Of_Probas_Across_Mixtures(r_mat_weight_sum, e_frq_weight_sum, curr_mixt_tree);
            }
          else
            {
              r_mat_weight_sum = 1.;
              e_frq_weight_sum = 1.;
              sum_probas       = 1.;
            }

          Ns   = curr_mixt_tree->next ? curr_mixt_tree->next->mod->ns : curr_mixt_tree->mod->ns;
          NsNs = Ns*Ns;
          NsNg = Ns*curr_mixt_tree->mod->ras->n_catg;

          p = (phydbl *)mCalloc(Ns,sizeof(phydbl));

          /* for(site=0;site<curr_mixt_tree->n_pattern;site++) // For each site in the current partition element */
          for(site=0;site<curr_mixt_tree->data->init_len;site++) // For each site in the current partition element
            {
              csite = curr_mixt_tree->data->sitepatt[site];

              d    = curr_mixt_d->next ? curr_mixt_d->next : curr_mixt_d;
              tree = curr_mixt_tree->next ? curr_mixt_tree->next : curr_mixt_tree;

              for(i=0;i<tree->mod->ns;i++) p[i] = .0;

              do // For each class of the mixture model that applies to the current partition element
                {
                  if(tree->is_mixt_tree == YES)
                    {
                      tree = tree->next;
                      d    = d->next;
                    }

                  v0 = d->v[0];
                  v1 = d->v[1];
                  v2 = d->v[2];

                  b0 = d->b[0];
                  b1 = d->b[1];
                  b2 = d->b[2];

                  Pij0 = b0->Pij_rr;
                  Pij1 = b1->Pij_rr;
                  Pij2 = b2->Pij_rr;

                  if(v0 == b0->left)
                    {
                      p_lk0 = b0->p_lk_left;
                      sum_scale0 = b0->sum_scale_left;
                    }
                  else
                    {
                      p_lk0 = b0->p_lk_rght;
                      sum_scale0 = b0->sum_scale_rght;
                    }

                  if(v1 == b1->left)
                    {
                      p_lk1 = b1->p_lk_left;
                      sum_scale1 = b1->sum_scale_left;
                    }
                  else
                    {
                      p_lk1 = b1->p_lk_rght;
                      sum_scale1 = b1->sum_scale_rght;
                    }

                  if(v2 == b2->left)
                    {
                      p_lk2 = b2->p_lk_left;
                      sum_scale2 = b2->sum_scale_left;
                    }
                  else
                    {
                      p_lk2 = b2->p_lk_rght;
                      sum_scale2 = b2->sum_scale_rght;
                    }


                  for(catg=0;catg<tree->mod->ras->n_catg;catg++)
                    {
                      for(i=0;i<Ns;i++)
                        {
                          p0 = .0;
                          if(v0->tax)
                            for(j=0;j<tree->mod->ns;j++)
                              {
                                p0 += v0->b[0]->p_lk_tip_r[csite*Ns+j] * Pij0[catg*NsNs+i*Ns+j];

                                /* printf("\n. p0 %d %f", */
                                /*        v0->b[0]->p_lk_tip_r[site*Ns+j], */
                                /*        Pij0[catg*NsNs+i*Ns+j]); */
                              }
                          else
                            for(j=0;j<tree->mod->ns;j++)
                              {
                                p0 += p_lk0[csite*NsNg+catg*Ns+j] * Pij0[catg*NsNs+i*Ns+j] / (phydbl)POW(2,sum_scale0[catg*curr_mixt_tree->n_pattern+csite]);

                                /* p0 += p_lk0[site*NsNg+catg*Ns+j] * Pij0[catg*NsNs+i*Ns+j]; */

                                /* printf("\n. p0 %f %f", */
                                /*        p_lk0[site*NsNg+catg*Ns+j], */
                                /*        Pij0[catg*NsNs+i*Ns+j]); */
                              }
                          p1 = .0;
                          if(v1->tax)
                            for(j=0;j<tree->mod->ns;j++)
                              {
                                p1 += v1->b[0]->p_lk_tip_r[csite*Ns+j] * Pij1[catg*NsNs+i*Ns+j];

                                /* printf("\n. p1 %d %f", */
                                /*        v1->b[0]->p_lk_tip_r[site*Ns+j], */
                                /*        Pij1[catg*NsNs+i*Ns+j]); */
                              }

                          else
                            for(j=0;j<tree->mod->ns;j++)
                              {
                                p1 += p_lk1[csite*NsNg+catg*Ns+j] * Pij1[catg*NsNs+i*Ns+j] / (phydbl)POW(2,sum_scale1[catg*curr_mixt_tree->n_pattern+csite]);

                                /* p1 += p_lk1[site*NsNg+catg*Ns+j] * Pij1[catg*NsNs+i*Ns+j];  */

                                /* printf("\n. p1 %f %f", */
                                /*        p_lk1[site*NsNg+catg*Ns+j], */
                                /*        Pij1[catg*NsNs+i*Ns+j]); */
                             }


                          p2 = .0;
                          if(v2->tax)
                            for(j=0;j<tree->mod->ns;j++)
                              {
                                p2 += v2->b[0]->p_lk_tip_r[csite*Ns+j] * Pij2[catg*NsNs+i*Ns+j];
                                /* printf("\n. p2 %d %f", */
                                /*        v2->b[0]->p_lk_tip_r[site*Ns+j], */
                                /*        Pij2[catg*NsNs+i*Ns+j]); */
                              }
                          else
                            for(j=0;j<tree->mod->ns;j++)
                              {
                                p2 += p_lk2[csite*NsNg+catg*Ns+j] * Pij2[catg*NsNs+i*Ns+j] / (phydbl)POW(2,sum_scale2[catg*curr_mixt_tree->n_pattern+csite]);

                                /* p2 += p_lk2[site*NsNg+catg*Ns+j] * Pij2[catg*NsNs+i*Ns+j]; */

                                /* printf("\n. p2 %f %f", */
                                /*        p_lk2[site*NsNg+catg*Ns+j], */
                                /*        Pij2[catg*NsNs+i*Ns+j]);  */
                             }

                          p[i] +=
                            p0*p1*p2*
                            tree->mod->e_frq->pi->v[i] /
                            tree->cur_site_lk[csite] *
                            curr_mixt_tree->mod->ras->gamma_r_proba->v[tree->mod->ras->parent_class_number] *
                            tree->mod->r_mat_weight->v / r_mat_weight_sum *
                            tree->mod->e_frq_weight->v / e_frq_weight_sum /
                            sum_probas;

                          if(print == YES)
                            printf("\n class: %d prob: %f",
                                   tree->mod->ras->parent_class_number,
                                   curr_mixt_tree->mod->ras->gamma_r_proba->v[tree->mod->ras->parent_class_number]);
                        }
                    }
                  
                  if(print == YES)
                    {
                      PhyML_Fprintf(fp,"%4d\t%4d\t",site+1,d->num);
                      for(i=0;i<Ns;i++)
                        {
                          PhyML_Fprintf(fp,"%.4f\t",p[i]);
                        }
                      PhyML_Fprintf(fp,"\n");
                      fflush(NULL);
                    }

                  tree = tree->next;
                  d    = d->next;


                }
              while(tree && d && tree->is_mixt_tree == NO);
            }

          Free(p);
          curr_mixt_tree = curr_mixt_tree->next_mixt;
          curr_mixt_d    = curr_mixt_d->next_mixt;
        }
      while(curr_mixt_tree != NULL);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// First and second derivative of the log-likelihood with respect
// to the length of edge b
phydbl MIXT_dLk(phydbl *l, t_edge *mixt_b, t_tree *mixt_tree)
{
  t_tree *tree,*cpy_mixt_tree;
  t_edge *b,*cpy_mixt_b;
  unsigned int site, class, nclasses, ns, state;
  phydbl *sum_scale_left_cat,*sum_scale_rght_cat;
  phydbl sum,mult;
  phydbl *lk,*dlk;
  phydbl site_lk,site_dlk;
  phydbl log_site_lk,inv_site_lk;
  int num_prec_issue;
  phydbl r_mat_weight_sum, e_frq_weight_sum, sum_probas;
  phydbl len;
  phydbl *expl,*dot_prod;
  phydbl rr;
  phydbl ev,expevlen;
  phydbl one_m_pinv;
  short int length_found;
  
  tree          = NULL;
  b             = NULL;
  expl          = NULL;
  lk            = NULL;
  dlk           = NULL;
  cpy_mixt_tree = mixt_tree;
  cpy_mixt_b    = mixt_b;
  len           = -1.;
  
  if(mixt_tree->update_eigen_lr == YES)
    {
      MIXT_Update_Eigen_Lr(mixt_b,mixt_tree);
    }

  /*! Make sure that l is one of the lengths of mixt_b */
  b = mixt_b;
  while(b) { if(&(b->l->v) == l) break; b = b->next; }
  assert(b != NULL);

  
  do /*! Consider each element of the data partition */
    {
      b    = mixt_b;
      tree = mixt_tree;
      length_found = NO;
      do
        {
          if(&(b->l->v) == l)
            {
              length_found = YES;
              break;
            }

          b = b->next;
          tree = tree->next;
        }
      while(tree && tree->is_mixt_tree == NO);

      /*! For computational efficiency, the dlk and lk computation for
        data partition elements corresponding to trees in which l is not
        found could be skipped. We compute dlk and lk nonetheless so that
        these two values are the derivative and lielihood computed for the
        whole data set, not juste the data partitions that "contain" l
      */
      
      tree = mixt_tree;
      do
        {
          tree->c_lnL  = .0;
          tree->c_dlnL = .0;
          tree = tree->next;
        }
      while(tree);
      
      ns = mixt_tree->mod->ns;
      nclasses = MIXT_Mixt_Size(mixt_tree);
      
      lk   = (phydbl *)mCalloc(nclasses,sizeof(phydbl));
      dlk  = (phydbl *)mCalloc(nclasses,sizeof(phydbl));
      
      sum_scale_left_cat = (phydbl *)mCalloc(nclasses,sizeof(phydbl));
      sum_scale_rght_cat = (phydbl *)mCalloc(nclasses,sizeof(phydbl));
      
      r_mat_weight_sum = MIXT_Get_Sum_Chained_Scalar_Dbl(mixt_tree->next->mod->r_mat_weight);
      e_frq_weight_sum = MIXT_Get_Sum_Chained_Scalar_Dbl(mixt_tree->next->mod->e_frq_weight);
      sum_probas       = MIXT_Get_Sum_Of_Probas_Across_Mixtures(r_mat_weight_sum,e_frq_weight_sum,mixt_tree);
      
      // Fill in expl vector for each rate class
      tree = mixt_tree->next;
      b = mixt_b->next;
      class = 0;
      do
        {          
          if(tree->mod->ras->invar == YES) rr = 0.0;
          else                             
            {
              rr = 1.0;
              rr *= tree->mod->br_len_mult->v;
              rr *= mixt_tree->mod->ras->gamma_rr->v[tree->mod->ras->parent_class_number];
            }
          

          if(length_found == YES) len = (*l) * rr;  
          else len = b->l->v * rr;
          
          if(isinf(len) || isnan(len)) 
            {
              PhyML_Fprintf(stderr,"\n. len=%f rr=%f l=%f",len,rr,*l);
              assert(FALSE);
            }
          
          if(tree->mod->ras->invar == NO) 
            {
              if(len < tree->mod->l_min)      len = tree->mod->l_min;
              else if(len > tree->mod->l_max) len = tree->mod->l_max;
            }
          else
            {
              len = 0.0;
            }
          
          
          expl = mixt_tree->expl;
          
          for(state=0;state<tree->mod->ns;++state) 
            {
              ev = tree->mod->eigen->e_val[state];
              expevlen = exp(ev*len);
              
              expl[class*2*ns + 2*state]     = expevlen;
              expl[class*2*ns + 2*state + 1] = expevlen*ev*rr;
            }
          
          class++;
          
          tree = tree->next;
          b = b->next;
        }
      while(tree && tree->is_mixt_tree == NO); 
      
      assert(class == nclasses);
      
      mixt_tree->c_lnL  = .0;
      mixt_tree->c_dlnL = .0;
      
      for(site=0;site<mixt_tree->n_pattern;++site)
        {
          for(class=0;class<nclasses;++class)
            {
              dlk[class] = 0.0;
              lk[class]  = 0.0;
            }
          
          
          b     = mixt_b->next;
          tree  = mixt_tree->next;          
          class = 0;
          /*! For all classes in the mixture */
          do
            {
              if(tree->mod->ras->invar == NO && tree->data->wght[tree->curr_site] > SMALL)
                {
                  tree->curr_site = site;
                  dot_prod        = tree->dot_prod + site * ns;
                  expl            = mixt_tree->expl + 2*ns*class;
                  
                  if(tree->mod->io->datatype == NT || tree->mod->io->datatype == AA)
                    {
#if (defined(__AVX__))
                      AVX_Lk_dLk_Core_One_Class_Eigen_Lr(dot_prod,
                                                         expl ? expl : NULL,
                                                         ns,lk+class,dlk+class);
#elif (defined(__SSE3__))
                      SSE_Lk_dLk_Core_One_Class_Eigen_Lr(dot_prod,
                                                         expl ? expl : NULL,
                                                         ns,lk+class,dlk+class);
#else
                      Lk_dLk_Core_One_Class_Eigen_Lr(dot_prod,
                                                     expl ? expl : NULL,
                                                     ns,lk+class,dlk+class);
#endif
                    }
                  else
                    {
                      Lk_dLk_Core_One_Class_Eigen_Lr(dot_prod,
                                                     expl ? expl : NULL,
                                                     ns,lk+class,dlk+class);
                    }
                  
                  tree->apply_lk_scaling = YES;

                  
                  sum_scale_left_cat[tree->mod->ras->parent_class_number] =
                    (b->sum_scale_left)?
                    (b->sum_scale_left[site]):
                    (0.0);
                  
                  sum_scale_rght_cat[tree->mod->ras->parent_class_number] =
                    (b->sum_scale_rght)?
                    (b->sum_scale_rght[site]):
                    (0.0);
                  
                  sum =
                    sum_scale_left_cat[tree->mod->ras->parent_class_number] +
                    sum_scale_rght_cat[tree->mod->ras->parent_class_number];
                  
                  if(sum > 1024.)
                    {
                      /* PhyML_Fprintf(stderr,"\n. Numerical precision issue detected (sum = %g)!!!",sum); */
                      sum = 1023.;
                      tree->mixt_tree->numerical_warning = YES;
                      tree->numerical_warning = YES;
                    }
                  
                  mult = pow(2,sum);
                  
                  /* PhyML_Printf("\n> tree: %p class: %d lk: %f dlk: %g sum: %g site: %d dot_prod[0]: %g expl[0]: %g expl[1]: %g expl[2]: %g",tree,class,log(lk[class]),dlk[class],sum,site,dot_prod[0],expl[0],expl[1],expl[2]); */
                  
                  lk[class]  /= mult;
                  dlk[class] /= mult;
                  
                }
              
              tree = tree->next;
              b    = b->next;
              class++;
              
            }
          while(tree && tree->is_mixt_tree == NO); 
          // done with all trees in the mixture for this partition element.
          
          tree      = mixt_tree->next;
          class     = 0;
          site_lk   = .0;
          site_dlk  = .0;
          
          if(mixt_tree->mod->ras->invar == YES) one_m_pinv = 1. - mixt_tree->mod->ras->pinvar->v;
          else one_m_pinv = 1.;
          
          
          do
            {              
              if(tree->mod->ras->invar == NO && mixt_tree->data->wght[tree->curr_site] > SMALL)
                {
                  site_lk +=
                    lk[class] *
                    mixt_tree->mod->ras->gamma_r_proba->v[tree->mod->ras->parent_class_number] *
                    tree->mod->r_mat_weight->v / r_mat_weight_sum *
                    tree->mod->e_frq_weight->v / e_frq_weight_sum /
                    sum_probas;
                  
                  
                  site_dlk +=
                    dlk[class] *
                    one_m_pinv *
                    mixt_tree->mod->ras->gamma_r_proba->v[tree->mod->ras->parent_class_number] *
                    tree->mod->r_mat_weight->v / r_mat_weight_sum *
                    tree->mod->e_frq_weight->v / e_frq_weight_sum /
                    sum_probas;
                                    
                  /* PhyML_Printf("\n. MIXT_DLK site: %d mixt_tree: %p tree: %p class: %d l: %p lk: %g dlk: %g", */
                  /*              site, */
                  /*              mixt_tree, */
                  /*              tree, */
                  /*              class, */
                  /*              l, */
                  /*              lk[class], */
                  /*              dlk[class]); */
                  
                }
              
              class++;
              tree = tree->next;
            }
          while(tree && tree->is_mixt_tree == NO);

          
          
          /* Scaling for invariants */
          if(mixt_tree->mod->ras->invar == YES)
            {
              num_prec_issue = NO;
              
              tree = mixt_tree->next;
              while(tree->mod->ras->invar == NO)
                {
                  tree = tree->next;
                  if(!tree || tree->is_mixt_tree == YES)
                    {
                      PhyML_Fprintf(stderr,"\n. tree: %p",tree);
                      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d",__FILE__,__LINE__);
                      Exit("\n");
                    }
                }
              
              tree->apply_lk_scaling = YES;
              
              /*! 'tree' will give the correct state frequencies (as opposed to mixt_tree) */              
              inv_site_lk = Invariant_Lk(0,site,&num_prec_issue,tree);
              
              
              if(num_prec_issue == YES) // inv_site_lk >> site_lk
                {
                  site_lk = inv_site_lk * mixt_tree->mod->ras->pinvar->v;
                }
              else
                {
                  site_lk = site_lk * (1. - mixt_tree->mod->ras->pinvar->v) + inv_site_lk * mixt_tree->mod->ras->pinvar->v;
                }
            }
          
          
          log_site_lk = log(site_lk);
          
          if(isinf(log_site_lk) || isnan(log_site_lk))
            {
              PhyML_Fprintf(stderr,"\n. site = %d",site);
              PhyML_Fprintf(stderr,"\n. invar = %d",mixt_tree->data->invar[site]);
              PhyML_Fprintf(stderr,"\n. mixt = %d",mixt_tree->is_mixt_tree);
              PhyML_Fprintf(stderr,"\n. lk = %G log(lk) = %f < %G",site_lk,log_site_lk,-BIG);
              for(class=0;class<mixt_tree->mod->ras->n_catg;++class) PhyML_Printf("\n. rr=%f p=%f",mixt_tree->mod->ras->gamma_rr->v[class],mixt_tree->mod->ras->gamma_r_proba->v[class]);
              PhyML_Fprintf(stderr,"\n. pinv = %G",mixt_tree->mod->ras->pinvar->v);
              PhyML_Fprintf(stderr,"\n. bl mult = %G",mixt_tree->mod->br_len_mult->v);
              PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d.\n",__FILE__,__LINE__);
              Exit("\n");
            }
          
          
          mixt_tree->c_lnL  += mixt_tree->data->wght[site] * log_site_lk;
          mixt_tree->c_dlnL += mixt_tree->data->wght[site] * ( site_dlk / site_lk );

          /* PhyML_Printf("\n. site: %4d lnL: %15f dlnL: %15f",site,site_lk,site_dlk); */
        }
      
      Free(sum_scale_left_cat);
      Free(sum_scale_rght_cat);
      
      Free(lk);
      Free(dlk);
    
      mixt_tree = mixt_tree->next_mixt;
      mixt_b    = mixt_b->next_mixt;
    }
  while(mixt_tree);


  // Sum c_dlnL across all partition element and set each c_dlnL equal to
  // this sum  

  sum = .0;
  mixt_tree = cpy_mixt_tree;
  do
    {
      sum += mixt_tree->c_dlnL;
      mixt_tree = mixt_tree->next_mixt;
    }
  while(mixt_tree);

  mixt_tree = cpy_mixt_tree;
  do
    {
      mixt_tree->c_dlnL = sum;
      mixt_tree = mixt_tree->next_mixt;
    }
  while(mixt_tree);



  // Do the same with likelihood values
  sum = .0;
  mixt_tree = cpy_mixt_tree;
  do
    {
      sum += mixt_tree->c_lnL;
      mixt_tree = mixt_tree->next_mixt;
    }
  while(mixt_tree);

  mixt_tree = cpy_mixt_tree;
  do
    {
      mixt_tree->c_lnL = sum;
      mixt_tree = mixt_tree->next_mixt;
    }
  while(mixt_tree);

  mixt_tree = cpy_mixt_tree;
  mixt_b    = cpy_mixt_b;

  return mixt_tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Returns the number of trees (including mixt_trees) in the
// whole model.
int MIXT_Part_Mixt_Size(t_tree *mixt_tree)
{
  if(mixt_tree->is_mixt_tree == NO) return 1;
  else
    {
      int num;
      t_tree *tree;
      
      num = 0;
      tree = mixt_tree;
      do
        {
          num++;
          tree = tree->next;
        }
      while(tree);
      return num;
    }
  return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Returns the number of trees in the partition element corresponding to mixt_tree.
int MIXT_Mixt_Size(t_tree *mixt_tree)
{
  if(mixt_tree->is_mixt_tree == NO) return 1;
  else
    {
      int num;
      t_tree *tree;
      
      num = 0;
      tree = mixt_tree->next;
      do
        {
          num++;
          tree = tree->next;
        }
      while(tree && tree->is_mixt_tree == NO);
 
     return num;
    }
  return -1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Both_Sides(int yesno, t_tree *mixt_tree)
{
  t_tree *tree;
  
  assert(mixt_tree->is_mixt_tree == YES);

  tree = mixt_tree;
  
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;        
      Set_Both_Sides(yesno,tree);
      tree = tree->next;
    }
  while(tree);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Model_Parameters(t_mod *mixt_mod)
{
  t_mod *mod = mixt_mod->next;

  if(mod != NULL)
    {
      do
        {
          Set_Model_Parameters(mod);
          mod = mod->next;
        }
      while(mod);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Update_Eigen(t_mod *mixt_mod)
{
  t_mod *mod;

  mod = mixt_mod;

  do
    {
      if(mod->is_mixt_mod) mod = mod->next;
      if(!Update_Boundaries(mod)) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      if(!Update_Eigen(mod)) Generic_Exit(__FILE__,__LINE__,__FUNCTION__);
      mod = mod->next;
    }
  while(mod);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Update_Eigen_Lr(t_edge *mixt_b, t_tree *mixt_tree)
{
  t_edge *b;
  t_tree *tree;

  tree = mixt_tree;
  b    = mixt_b;

  do
    {
      if(tree->is_mixt_tree == YES)
        {
          tree = tree->next;
          b    = b->next;
        }
      
      Update_Eigen_Lr(b,tree);

      tree = tree->next;
      b    = b->next;
    }
  while(tree);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Update_Eigen(int yn, t_mod *mixt_mod)
{
  t_mod *mod;

  mod = mixt_mod;
  do
    {      
      mod->update_eigen = yn;
      mod = mod->next;
    }
  while(mod); 
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Update_Eigen_Lr(int yn, t_tree *mixt_tree)
{
  t_tree *tree;
  int dum;
  
  tree = mixt_tree->next;
  do
    {
      dum = mixt_tree->is_mixt_tree;
      mixt_tree->is_mixt_tree = NO;

      Set_Update_Eigen_Lr(yn,tree);

      mixt_tree->is_mixt_tree = dum;
      
      tree = tree->next;
    }
  while(tree); 
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Use_Eigen_Lr(int yn, t_tree *mixt_tree)
{
  t_tree *tree;
  int dum;

  tree = mixt_tree->next;
  do
    {
      dum = mixt_tree->is_mixt_tree;
      mixt_tree->is_mixt_tree = NO;
      
      Set_Use_Eigen_Lr(yn,tree);

      mixt_tree->is_mixt_tree = dum;

      tree = tree->next;
    }
  while(tree); 
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Sample_Ancestral_Seq(int fullmutmap, int fromprior, t_tree *mixt_tree)
{
  t_tree *tree,*loc_mixt_tree;
  phydbl r_mat_weight_sum, e_frq_weight_sum, sum_probas, *class_proba;
  int i;
  
  r_mat_weight_sum = MIXT_Get_Sum_Chained_Scalar_Dbl(mixt_tree->next->mod->r_mat_weight);
  e_frq_weight_sum = MIXT_Get_Sum_Chained_Scalar_Dbl(mixt_tree->next->mod->e_frq_weight);
  sum_probas       = MIXT_Get_Sum_Of_Probas_Across_Mixtures(r_mat_weight_sum,e_frq_weight_sum,mixt_tree);


  tree = mixt_tree;
  loc_mixt_tree = tree;

  // For each element in the partition (i.e., each gene), select a
  // class of the mixture according to its post prob.
  do
    {
      tree = loc_mixt_tree->next;

      class_proba = NULL;
      i = 0;
      do
        {
          if(i == 0) class_proba = (phydbl *)mCalloc(1,sizeof(phydbl));
          else class_proba = (phydbl *)mRealloc(class_proba,i+1,sizeof(phydbl));
          
          class_proba[i] =
            tree->site_lk_cat[0] *
            loc_mixt_tree->mod->ras->gamma_r_proba->v[tree->mod->ras->parent_class_number] *
            tree->mod->r_mat_weight->v / r_mat_weight_sum *
            tree->mod->e_frq_weight->v / e_frq_weight_sum /
            sum_probas;
          
          tree = tree->next;
          i++;
        }
      while(tree && tree->is_mixt_tree == NO);

      tree = loc_mixt_tree->next;
      i = Sample_i_With_Proba_pi(class_proba,i);

      do
        {
          i--;
          if(i < 0) break;
          tree = tree->next;
          assert(tree);
        }
      while(1);

      assert(tree->is_mixt_tree == NO);
      
      Sample_Ancestral_Seq(fullmutmap,fromprior,tree);
      
      Free(class_proba);

      loc_mixt_tree = loc_mixt_tree->next_mixt;
    }
  while(loc_mixt_tree);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// Propagate changes in the "master" tree to other trees. Partial
// likelihood vectors (and scaling stuff) is not update here so
// requires full tree traversal to be updated properly.
void MIXT_Propagate_Tree_Update(t_tree *mixt_tree)
{
  t_tree *tree;
  int i;
  
  assert(!mixt_tree->prev);

  tree = mixt_tree->next;
  if(!tree) return;
  
  do
    {
      for(i=0;i<2*mixt_tree->n_otu-1;++i)
        {          
          tree->a_nodes[i]->v[0] = mixt_tree->a_nodes[i]->v[0] ? tree->a_nodes[mixt_tree->a_nodes[i]->v[0]->num] : NULL;
          tree->a_nodes[i]->v[1] = mixt_tree->a_nodes[i]->v[1] ? tree->a_nodes[mixt_tree->a_nodes[i]->v[1]->num] : NULL;
          tree->a_nodes[i]->v[2] = mixt_tree->a_nodes[i]->v[2] ? tree->a_nodes[mixt_tree->a_nodes[i]->v[2]->num] : NULL;

          if(tree->rates != NULL) tree->rates->nd_t[i] = mixt_tree->rates->nd_t[i];
        }

      Connect_Edges_To_Nodes_Serial(tree);          

      Reorganize_Edges_Given_Lk_Struct(tree);
      Init_Partial_Lk_Tips_Double(tree);
      
      if(mixt_tree->n_root != NULL)
        {
          assert(mixt_tree->e_root);
          assert(mixt_tree->n_root->v[1]);
          assert(mixt_tree->n_root->v[2]);
          assert(mixt_tree->n_root->b[1]);
          assert(mixt_tree->n_root->b[2]);
          
          tree->n_root = tree->a_nodes[mixt_tree->n_root->num];
          tree->e_root = tree->a_edges[mixt_tree->e_root->num];
          
          tree->n_root->v[1] = tree->a_nodes[mixt_tree->n_root->v[1]->num];
          tree->n_root->v[2] = tree->a_nodes[mixt_tree->n_root->v[2]->num];

          tree->n_root->b[1] = tree->a_edges[mixt_tree->n_root->b[1]->num];
          tree->n_root->b[2] = tree->a_edges[mixt_tree->n_root->b[2]->num];

          tree->n_root->b[1]->left = tree->n_root;
          tree->n_root->b[2]->left = tree->n_root;
          tree->n_root->b[1]->rght = tree->n_root->v[1];
          tree->n_root->b[2]->rght = tree->n_root->v[2];
        }
      
      Update_Ancestors(tree->n_root,tree->n_root->v[2],tree);
      Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);

      tree = tree->next;
      
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Set_Bl_From_Rt(int yn, t_tree *mixt_tree)
{
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      assert(tree->rates);
      tree->rates->bl_from_rt = yn;
      tree = tree->next;
    }
  while(tree);  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Copy_Tree(t_tree *ori, t_tree *cpy)
{
  int ori_is_mixt_tree,cpy_is_mixt_tree;

  assert(!((cpy && !ori) || (!cpy && ori)));
  
  if(cpy == NULL || ori == NULL) return;
  
  if(ori->is_mixt_tree == YES && cpy->is_mixt_tree == YES)
    {
      do
        {
          ori_is_mixt_tree = ori->is_mixt_tree;
          ori->is_mixt_tree = NO;
          cpy_is_mixt_tree = cpy->is_mixt_tree;      
          cpy->is_mixt_tree = NO;
          
          Copy_Tree(ori,cpy);
          
          cpy->is_mixt_tree = cpy_is_mixt_tree;
          ori->is_mixt_tree = ori_is_mixt_tree;
          
          ori = ori->next;
          cpy = cpy->next;
        }
      while(cpy);
    }
  else if(ori->is_mixt_tree == NO && cpy->is_mixt_tree == YES)
    {
      do
        {
          cpy_is_mixt_tree = cpy->is_mixt_tree;      
          cpy->is_mixt_tree = NO;
          
          Copy_Tree(ori,cpy);
          
          cpy->is_mixt_tree = cpy_is_mixt_tree;
          
          cpy = cpy->next;
        }
      while(cpy);      
    }
  else if(ori->is_mixt_tree == YES && cpy->is_mixt_tree == NO)
    {
      ori_is_mixt_tree = ori->is_mixt_tree;      
      ori->is_mixt_tree = NO;
      Copy_Tree(ori,cpy);
      ori->is_mixt_tree = ori_is_mixt_tree;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Init_NNI_Score(phydbl val, t_edge *mixt_b, t_tree *mixt_tree)
{
  t_edge *b;
  t_tree *tree;

  tree = mixt_tree->next;
  b = mixt_b->next;
  do
    {
      if(tree->is_mixt_tree)
        {
          tree = tree->next;
          b = b->next;
        }
      
      Init_NNI_Score(val,b,tree);
      tree = tree->next;
      b = b->next;
    }
  while(tree);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_tree *MIXT_Duplicate_Tree(t_tree *ori)
{
  t_tree *cpy,*first,*prev;
  int dum;
  
  ori->is_mixt_tree = NO;
  first = Duplicate_Tree(ori);
  first->prev = NULL;
  first->is_mixt_tree = YES;
  ori->is_mixt_tree = YES;

  ori = ori->next;
  prev = cpy = first;
  do
    {
      dum = ori->is_mixt_tree;
      ori->is_mixt_tree = NO;

      prev = cpy;

      cpy = Duplicate_Tree(ori);

      ori->is_mixt_tree = dum;      
      cpy->is_mixt_tree = ori->is_mixt_tree; 

      cpy->prev = prev;
      prev->next = cpy;
      
      ori = ori->next;
    }
  while(ori);

  cpy->next = NULL;
  
  return(first);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MIXT_Print_Site_Lk(t_tree *mixt_tree, FILE *fp)
{
  t_tree *tree;
  int site,catg;
  char *s;
  phydbl postmean,sum;
  
  assert(mixt_tree->is_mixt_tree == YES);
  assert(mixt_tree->io->print_site_lnl == YES);

  
  if(!mixt_tree->io->print_trace)
    {
      tree = mixt_tree;
      do
        {
          PhyML_Fprintf(fp,"Note : P(D|M) is the probability of site D given the model M (i.e., the site likelihood)\n");
          if(tree->mod->ras->n_catg > 1 || tree->mod->ras->invar)
            {
              PhyML_Fprintf(fp,"P*(D|M,rr[x]) is the scaled probability of site D given the model M and the relative rate\n");
              PhyML_Fprintf(fp,"of evolution rr[x], where x is the class of rate to be considered.\n");
              PhyML_Fprintf(fp,"The actual conditional probability is given by P*(D|M,rr[x])/2^F, where\n");
              PhyML_Fprintf(fp,"F is the scaling factor (see column 'Scaler').\n");
              PhyML_Fprintf(fp,"For invariant sites, P(D|M,rr[0]=0) is the actual conditional probability\n");
              PhyML_Fprintf(fp,"(i.e., it is not scaled).\n");
              break;
            }
          tree = tree->next_mixt;
        }
      while(tree);
      PhyML_Fprintf(fp,"\n");
    

      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));

      do
        {
          PhyML_Fprintf(fp,"Alignment file name: %s\n\n",mixt_tree->io->in_align_file);

          sprintf(s,"Site");
          PhyML_Fprintf(fp, "%-12s",s);
          
          sprintf(s,"P(D|M)");
          PhyML_Fprintf(fp,"%-15s",s);
          
          sprintf(s,"Scaler");
          PhyML_Fprintf(fp,"%-7s",s);
          
          sprintf(s,"Pattern");
          PhyML_Fprintf(fp, "%-9s",s);
          
          if(mixt_tree->mod->ras->n_catg > 1)
            {
              for(catg=0;catg<mixt_tree->mod->ras->n_catg;catg++)
                {
                  sprintf(s,"P*(D|M,rr[%d]=%5.4f)",catg+1,mixt_tree->mod->ras->gamma_rr->v[catg]);
                  PhyML_Fprintf(fp,"%-23s",s);
                }
              
              sprintf(s,"Posterior mean");
              PhyML_Fprintf(fp,"%-22s",s);
            }
          
          
          if(mixt_tree->mod->ras->invar)
            {
              sprintf(s,"P(D|M,rr[0]=0)");
              PhyML_Fprintf(fp,"%-16s",s);
            }
          
          sprintf(s,"NDistinctStates");
          PhyML_Fprintf(fp,"%-16s",s);
          
          PhyML_Fprintf(fp,"\n");
          
          assert(mixt_tree->next->is_mixt_tree == NO);
          Init_Ui_Tips(mixt_tree->next);
          
          for(site=0;site<mixt_tree->data->init_len;site++)
            {                    
              PhyML_Fprintf(fp,"%-12d",site+1);
              
              PhyML_Fprintf(fp,"%-15g",mixt_tree->cur_site_lk[mixt_tree->data->sitepatt[site]]);      
              
              PhyML_Fprintf(fp,"%-7d",mixt_tree->fact_sum_scale[mixt_tree->data->sitepatt[site]]);      
              
              PhyML_Fprintf(fp,"%-9d",mixt_tree->data->sitepatt[site]);
              
              if(mixt_tree->mod->ras->n_catg > 1)
                {
                  tree = mixt_tree->next;
                  do
                    {
                      PhyML_Fprintf(fp,"%-23g",tree->unscaled_site_lk_cat[tree->data->sitepatt[site]]);
                      tree = tree->next;
                    }
                  while(tree && tree->is_mixt_tree == NO);

                  
                  tree = mixt_tree->next;
                  postmean = .0;
                  sum = .0;
                  do
                    {
                      postmean +=
                        mixt_tree->mod->ras->gamma_rr->v[tree->mod->ras->parent_class_number] *
                        tree->unscaled_site_lk_cat[tree->data->sitepatt[site]] *
                        mixt_tree->mod->ras->gamma_r_proba->v[tree->mod->ras->parent_class_number];

                      sum +=
                        tree->unscaled_site_lk_cat[tree->data->sitepatt[site]] *
                        mixt_tree->mod->ras->gamma_r_proba->v[tree->mod->ras->parent_class_number];

                      tree = tree->next;
                    }
                  while(tree && tree->is_mixt_tree == NO);
                  
                  postmean /= sum;
                  
                  PhyML_Fprintf(fp,"%-22g",postmean);
                }
              
              if(mixt_tree->mod->ras->invar)
                {
                  if((phydbl)mixt_tree->data->invar[mixt_tree->data->sitepatt[site]] > -0.5)
                    PhyML_Fprintf(fp,"%-16g",mixt_tree->mod->e_frq->pi->v[mixt_tree->data->invar[mixt_tree->data->sitepatt[site]]]);
                  else
                    PhyML_Fprintf(fp,"%-16g",0.0);
                }
              
              assert(mixt_tree->next != NULL);
              PhyML_Fprintf(fp,"%-16d",Number_Of_Diff_States_One_Site(mixt_tree->data->sitepatt[site],mixt_tree->next));
              
              PhyML_Fprintf(fp,"\n");
            }
          Free(s);
          
          mixt_tree = mixt_tree->next_mixt;
        }
      while(mixt_tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////




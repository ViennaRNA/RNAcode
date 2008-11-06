/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */

#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "options.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "mc.h"
#include "m4.h"
#include "draw.h"
#include "rates.h"


/*********************************************************/

int MC_main(int argc, char **argv)
{
  seq **data;
  allseq *alldata;
  option *io;
  arbre *tree;
  int n_otu, num_data_set;
  int num_tree,tree_line_number,num_rand_tree;
  matrix *mat;
  model *mod;
  m4 *m4mod;
  time_t t_beg,t_end;
  div_t hour,min;
  phydbl best_lnL;
  int bootstrap_this_tree;
  int r_seed;


#ifdef QUIET
  setvbuf(stdout,NULL,_IOFBF,2048);
#endif

  tree = NULL;
  mod  = NULL;
  data = NULL;
  bootstrap_this_tree = 1;
  best_lnL = UNLIKELY;

  io = (option *)Get_Input(argc,argv);
  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
  srand(r_seed);
  Make_Model_Complete(io->mod);
  mod = io->mod;
  m4mod = mod->m4mod;
  if(io->in_tree) Test_Multiple_Data_Set_Format(io);
  else io->n_trees = 1;

  if(io->mod->s_opt->random_input_tree) bootstrap_this_tree = 0;

  mat = NULL;
  tree_line_number = 0;

  if((io->n_data_sets > 1) && (io->n_trees > 1))
    {
      io->n_data_sets = MIN(io->n_trees,io->n_data_sets);
      io->n_trees     = MIN(io->n_trees,io->n_data_sets);
    }

  For(num_data_set,io->n_data_sets)
    {
      n_otu = 0;
      best_lnL = UNLIKELY;
      data = Get_Seq(io,0);

      if(data)
	{
	  if(io->n_data_sets > 1) printf("\n. Data set [#%d]\n",num_data_set+1);
	  printf("\n. Compressing sequences...\n");
	  alldata = Compact_Seq(data,io);
	  Free_Seq(data,alldata->n_otu);
	  Check_Ambiguities(alldata,io->mod->datatype,io->mod->stepsize);

	  for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++)
	    {
	      if(!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;

	      For(num_rand_tree,io->mod->s_opt->n_rand_starts)
		{
		  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search == SPR_MOVE))
		    printf("\n. [Random start %3d/%3d]\n",num_rand_tree+1,io->mod->s_opt->n_rand_starts);

		  Init_Model(alldata,mod);
		  if(io->m4_model) M4_Init_Model(m4mod,alldata,mod);

		  if(!io->in_tree)
		    {
		      printf("\n. Computing pairwise distances...\n");
		      mat = ML_Dist(alldata,mod);
		      Fill_Missing_Dist(mat);
		      printf("\n. Building BIONJ tree...\n");
		      mat->tree = Make_Tree_From_Scratch(alldata->n_otu,alldata);
		      Bionj(mat);
		      tree      = mat->tree;
		      tree->mat = mat;

		    }
		  else
		    {
		      if((io->n_trees == 1) || (!num_tree))
			{
			  rewind(io->fp_in_tree);
			  tree_line_number = 0;
			}

		      if(io->n_trees > 1) printf("\n. Reading tree [#%d]\n",tree_line_number+1);
		      else printf("\n. Reading tree...\n");
		      fflush(NULL);

		      tree = Read_Tree_File(io->fp_in_tree);
		      tree_line_number++;

		      if(!tree)
			{
			  printf("\n. Input tree not found...\n");
			  Exit("\n\n");
			}

		      if(!tree->has_branch_lengths)
			{
			  printf("\n. Computing branch length estimates...\n");
			  Order_Tree_CSeq(tree,alldata);
			  mat = ML_Dist(alldata,mod);
			  mat->tree = tree;
			  mat->method = 0;
			  Bionj_Br_Length(mat);
			  tree->mat = mat;
			}

		      tree->mod        = mod;
		      tree->io         = io;
		      tree->data       = alldata;
		      tree->both_sides = 1;
		      tree->n_pattern  = tree->data->crunch_len/tree->mod->stepsize;
		    }

		  if(!tree) continue;

		  time(&t_beg);
		  time(&(tree->t_beg));

		  tree->mod         = mod;
		  tree->io          = io;
		  tree->data        = alldata;
		  tree->both_sides  = 1;
		  tree->n_pattern   = tree->data->crunch_len/tree->mod->stepsize;

		  if((!num_data_set) && (!num_tree) && (!num_rand_tree))
		    {
#ifndef BATCH
		      Check_Memory_Amount(tree);
#endif
		    }

		  Order_Tree_CSeq(tree,alldata);

		  if((tree->mod->s_opt->random_input_tree) && (tree->mod->s_opt->topo_search == SPR_MOVE))
		    {
		      printf("\n. Randomising the tree...\n");
		      Random_Tree(tree);
		    }

		  Fill_Dir_Table(tree);
		  Update_Dirs(tree);
		  Make_Tree_4_Pars(tree,alldata,alldata->init_len);
		  Make_Tree_4_Lk(tree,alldata,alldata->init_len);
		  tree->triplet_struct = Make_Triplet_Struct(mod);
		  Br_Len_Not_Involving_Invar(tree);

 		  if((tree->mod->s_opt->topo_search == SPR_MOVE) ||
		     (tree->mod->s_opt->topo_search == NNI_MOVE  &&
		      tree->mod->s_opt->spr_step_after_nnis))
		    {
		      Make_Spr_List(tree);
		      Make_Best_Spr(tree);
		    }

/* 		  M4_Site_Branch_Classification_Experiment(tree); */
/* 		  Exit("\n"); */

/* 		  M4_Detect_Site_Switches_Experiment(tree); */
/* 		  Site_Diversity(tree); */
/* 		  M4_Posterior_Prediction_Experiment(tree); */
/* 		  Exit("\n"); */

		  if(tree->mod->s_opt->opt_topo)
		    {
		      if(tree->mod->s_opt->topo_search == NNI_MOVE)
			{
			  Simu_Loop(tree);
			}
		      else
			{
			  if(tree->mod->s_opt->steph_spr)
			    {
			      Speed_Spr_Loop(tree);
			    }
			  else
			    {
			      Init_SPR(tree);
			      Optim_SPR(tree,0,ALL);
			      Clean_SPR(tree);
			    }
			}
		    }
		  else
		    {
		      if(tree->mod->s_opt->opt_num_param || tree->mod->s_opt->opt_bl)
			{
			  Round_Optimize(tree,tree->data);
			}
		      else
			{
			  Lk(tree);
			  Print_Lk(tree,"");
			}
		    }


		  Lk(tree);
		  printf("\n. Best unconstrained lnL = %f",tree->c_lnL);
		  fprintf(tree->io->fp_out_stats,"\n. Best unconstrained lnL = %f",tree->c_lnL);
		  fprintf(tree->io->fp_out_tree,"\n[unconstrained tree] = %s\n",Write_Tree(tree));
		  fflush(NULL);

		  Record_Br_Len(tree);

		  edge *best_edge;
		  tree->bl_from_node_stamps = 1;
		  tree->rates = RATES_Make_Rate_Struct(tree);
		  RATES_Init_Rate_Struct(tree->rates,tree);
		  best_edge = MC_Find_Best_Root_Position_Approx(tree);
		  MC_Least_Square_Node_Times(best_edge,tree);      
		  MC_Adjust_Node_Times(tree);
		  MC_Round_Optimize(tree);
		  MC_Estimate_Branch_Rate_Parameter(tree);
		  MC_Classify_Branch_Rates(tree);
/* 		  MC_Compute_Rates_And_Times_Least_Square_Adjustments(tree); */
		  RATES_Print_Rates(tree);
		  For(tree->rates->curr_mc_run,tree->rates->n_mc_runs) RATES_Monte_Carlo_Mean_Rates(tree);
		  
		  int i,j,k;
		  phydbl density;
		  phydbl **x, *where;

		  where = (phydbl *)mCalloc(1,sizeof(phydbl));
		  x = (phydbl **)mCalloc(1,sizeof(phydbl));
		  x[0] = (phydbl *)mCalloc(tree->rates->n_mc_runs,sizeof(phydbl));
		  
		  printf("\n. Kernel values \n");
		  For(i,2*tree->n_otu-1)
		    {
		      For(j,tree->rates->n_mc_runs)
			{
			  density = Univariate_Kernel_Density_Estimate(tree->rates->mc_mr[i][j],
								       tree->rates->mc_mr[i],
								       tree->rates->n_mc_runs);

			  printf("\nS %d %f %f",i,tree->rates->mc_mr[i][j],density);
			  
			  where[0] = tree->rates->mc_mr[i][j];
			  For(k,tree->rates->n_mc_runs) 
			    {
			      x[0][k] = tree->rates->mc_mr[i][k];
			    }
			  density = Multivariate_Kernel_Density_Estimate(where,
									 x,
									 tree->rates->n_mc_runs,
									 1);
			  	
			  printf("\nM %d %f %f",i,tree->rates->mc_mr[i][j],density);
			}
		    }

		  tree->bl_from_node_stamps = 0;


/* 		  /\* Constrain the tree to be molecular-clock like and find best root position *\/ */
/* 		  tree->bl_from_node_stamps = 1; */
/* 		  edge *best_edge; */
/* 		  best_edge = MC_Find_Best_Root_Position(tree); */
/* /\* 		  best_edge = Find_Edge_With_Label("ROOT",tree); *\/ */
/* 		  if(!best_edge) */
/* 		    { */
/* 		      printf("\n. Could not find the root !"); */
/* 		      Warn_And_Exit("\n"); */
/* 		    } */


		  if(tree->io->ratio_test) aLRT(tree);

		  Lk(tree);
		  printf("\n\n. Final log likelihood : %f",tree->c_lnL);


		  if((tree->c_lnL > best_lnL) && (io->mod->s_opt->n_rand_starts > 1))
		    {
		      best_lnL = tree->c_lnL;
		      io->fp_out_best_tree = (FILE *)fopen(io->out_best_tree_file,"w");
		      Print_Tree(io->fp_out_best_tree,tree);
		      fflush(NULL);
		      fclose(io->fp_out_best_tree);
		    }

		  if((tree->mod->bootstrap) && (bootstrap_this_tree))
		    {
		      if(num_rand_tree > 0) io->in_tree = 0;
		      Bootstrap(tree);
		      tree->mod->bootstrap = 0;
		    }

		  Br_Len_Involving_Invar(tree);
		  Print_Tree(io->fp_out_tree,tree);

		  Unconstraint_Lk(tree);
		  time(&t_end);
		  hour = div(t_end-t_beg,3600);
		  min  = div(t_end-t_beg,60  );
		  min.quot -= hour.quot*60;
		  printf("\n\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
		  printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
		  Print_Fp_Out(io->fp_out_stats,t_beg,t_end,tree,
			       io,num_data_set+1,
			       (tree->mod->s_opt->n_rand_starts > 1)?(num_rand_tree):(num_tree));

		  if(tree->io->print_site_lnl) Print_Site_Lk(tree,io->fp_out_lk);

		  /* Start from BioNJ tree */
		  if((num_rand_tree == io->mod->s_opt->n_rand_starts-1)
		     && (io->mod->s_opt->n_rand_starts > 1)
		     && (tree->mod->s_opt->random_input_tree))
		    {
		      num_rand_tree--;
		      tree->mod->s_opt->random_input_tree = 0;
		    }

		  if((num_rand_tree == io->mod->s_opt->n_rand_starts - 1) &&
		     (!tree->mod->s_opt->random_input_tree) &&
		     (io->mod->s_opt->n_rand_starts > 1))
		    {
		      if(tree->mod->bootstrap)
			{
			  num_rand_tree--;
			  io->in_tree = 1;
			  io->fp_in_tree = io->fp_out_best_tree;
			  bootstrap_this_tree  = 1;
			  io->fp_in_tree = (FILE *)fopen(io->out_best_tree_file,"r");
			}
		      else
			{
			  io->fp_out_best_tree = (FILE *)fopen(io->out_best_tree_file,"w");
			  Print_Tree(io->fp_out_best_tree,tree);
			  fflush(NULL);
			  fclose(io->fp_out_best_tree);
			}
		    }

 		  if((tree->mod->s_opt->topo_search == SPR_MOVE) ||
		     (tree->mod->s_opt->topo_search == NNI_MOVE  &&
		      tree->mod->s_opt->spr_step_after_nnis))
		    {
		      Free_Spr_List(tree);
		      Free_One_Spr(tree->best_spr);
		    }

		  if(tree->mat) Free_Mat(tree->mat);
		  Free_Triplet(tree->triplet_struct);
		  Free_Tree_Pars(tree);
		  Free_Tree_Lk(tree);
		  Free_Tree(tree);
		}
	      if(io->n_trees > 1 && io->n_data_sets > 1) break;
	    }
	  Free_Cseq(alldata);
	}
    }

  if(io->mod->s_opt->n_rand_starts > 1) printf("\n\n. Best log likelihood : %f\n",best_lnL);

  Free_Model(mod);

  if(io->fp_in_seq)    fclose(io->fp_in_seq);
  if(io->fp_in_tree)   fclose(io->fp_in_tree);
  if(io->fp_out_lk)    fclose(io->fp_out_lk);
  if(io->fp_out_tree)  fclose(io->fp_out_tree);
  if(io->fp_out_stats) fclose(io->fp_out_stats);

  Free_Input(io);
  return 0;
}

/*********************************************************/

void MC_Least_Square_Node_Times(edge *e_root, arbre *tree)
{

  /* Solve A.x = b, where x are the node time estimated
     under the least square criterion.

     A is a n x n matrix, with n being the number of
     nodes in a rooted tree (i.e. 2*n_otu-1).
   */

  phydbl *A, *b, *x;
  int n;
  int i,j;
  node *root;


  printf("\n. Making the tree molecular clock like.");

  n = 2*tree->n_otu-1;

  A = (phydbl *)mCalloc(n*n,sizeof(phydbl));
  b = (phydbl *)mCalloc(n,  sizeof(phydbl));
  x = (phydbl *)mCalloc(n,  sizeof(phydbl));

  (e_root)?(Add_Root(e_root,tree)):(Add_Root(tree->t_edges[0],tree));

  root = tree->n_root;

  MC_Least_Square_Node_Times_Pre(root,root->v[0],A,b,n,tree);
  MC_Least_Square_Node_Times_Pre(root,root->v[1],A,b,n,tree);

  b[root->num] = tree->e_root->l/2.;

  A[root->num * n + root->num]       = 1.0;
  A[root->num * n + root->v[0]->num] = -.5;
  A[root->num * n + root->v[1]->num] = -.5;

  Matinv(A, n, n, NULL);

  For(i,n) x[i] = .0;
  For(i,n) For(j,n) x[i] += A[i*n+j] * b[j];

  For(i,n-1) tree->noeud[i]->t = x[i] / tree->rates->mean_r;
  root->t = x[n-1] / tree->rates->mean_r;  

  Free(A);
  Free(b);
  Free(x);

}

/*********************************************************/

void MC_Least_Square_Node_Times_Pre(node *a, node *d, phydbl *A, phydbl *b, int n, arbre *tree)
{
  if(d->tax)
    {
      A[d->num * n + d->num] = 1.;
      
      /* Set the time stamp at tip nodes to 0.0 */
/*       printf("\n. Tip node date set to 0"); */
      b[d->num] = 0.0;
      return;
    }
  else
    {
      int i;
      
      For(i,3) 
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  MC_Least_Square_Node_Times_Pre(d,d->v[i],A,b,n,tree);
      
      A[d->num * n + d->num] = 1.;
      b[d->num] = .0;
      For(i,3)
	{
	  A[d->num * n + d->v[i]->num] = -1./3.;
	  if(d->v[i] != a) b[d->num] += d->b[i]->l;
	  else             b[d->num] -= d->b[i]->l;
	}
      b[d->num] /= 3.;
    }
}

/*********************************************************/

/* Adjust node times in order to have correct time stamp ranking with
 respect to the tree topology */

void MC_Adjust_Node_Times(arbre *tree)
{
  MC_Adjust_Node_Times_Pre(tree->n_root->v[0],tree->n_root->v[1],tree);
  MC_Adjust_Node_Times_Pre(tree->n_root->v[1],tree->n_root->v[0],tree);
  if(tree->n_root->t < MAX(tree->n_root->v[0]->t,tree->n_root->v[1]->t))
/*     tree->n_root->t = MAX(tree->n_root->v[0]->t,tree->n_root->v[1]->t) + BL_MIN; */
    tree->n_root->t = MAX(tree->n_root->v[0]->t,tree->n_root->v[1]->t);
}

/*********************************************************/

void MC_Adjust_Node_Times_Pre(node *a, node *d, arbre *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      phydbl max_height;

      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    MC_Adjust_Node_Times_Pre(d,d->v[i],tree);
	  }

      max_height = -1.0;
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      if(d->v[i]->t > max_height)
		{
		  max_height = d->v[i]->t;
		}
	    }
	}
      
      if(d->t < max_height) 
	{
/* 	  d->t = max_height + BL_MIN; */
	  d->t = max_height;
	}
    }
}
/*********************************************************/

  /* Multiply each time stamp at each internal 
     node by  'tree->time_stamp_mult'.
   */

void MC_Mult_Time_Stamps(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->noeud[i]->t *= fabs(tree->mod->s_opt->tree_size_mult);
  tree->n_root->t *= fabs(tree->mod->s_opt->tree_size_mult);
}

/*********************************************************/

/* Divide each time stamp at each internal 
   node by  'tree->time_stamp_mult'.
*/
void MC_Div_Time_Stamps(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->noeud[i]->t /= fabs(tree->mod->s_opt->tree_size_mult);
  tree->n_root->t /= fabs(tree->mod->s_opt->tree_size_mult);
}

/*********************************************************/

void MC_Bl_From_T(arbre *tree)
{
  phydbl mean_rate, branch_rate;

  /* Branch lengths are deduced from time stamps */
  MC_Bl_From_T_Post(tree->n_root,tree->n_root->v[0],NULL,tree);
  MC_Bl_From_T_Post(tree->n_root,tree->n_root->v[1],NULL,tree);

  mean_rate   = tree->rates->mean_r;
  branch_rate = tree->rates->br_r[tree->e_root->num];

  tree->e_root->l = 
    mean_rate * branch_rate * (tree->n_root->t - tree->e_root->left->t) + 
    mean_rate * branch_rate * (tree->n_root->t - tree->e_root->rght->t);

  /* -tree->e_root->left->t - tree->e_root->rght->t + 2*tree->n_root->t; */

  /* Actual formula =>  tree->e_root->l = 
     (tree->n_root->t - tree->e_root->left->t) + 
     (tree->n_root->t - tree->e_root->rght->t); */
  
  tree->n_root_pos = (tree->n_root->t - tree->e_root->left->t)/tree->e_root->l;

}

/*********************************************************/

void MC_Bl_From_T_Post(node *a, node *d, edge *b, arbre *tree)
{

  if(b)
    {
      phydbl mean_rate, branch_rate;

      mean_rate   = tree->rates->mean_r;
      branch_rate = tree->rates->br_r[b->num];

      b->l = (a->t - d->t) * mean_rate * branch_rate;
      if(b->l < 0.0)
	{
	  printf("\n. Correction failed.");
	  printf("\n. d->t = %f a->t = %f",d->t,a->t);
	  printf("\n. a->num=%d d->num=%d",a->num,d->num);
	  Warn_And_Exit("\n");
	}
    }

  if(d->tax) return;
  else
    {
      int i;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) MC_Bl_From_T_Post(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/

void MC_Round_Optimize(arbre *tree)
{
  int n_round,each;
  phydbl lk_old, lk_new, tol;

  lk_new = UNLIKELY;
  lk_old = UNLIKELY;
  n_round = 0;
  each = 5;
  tol = 1.e-2;

  tree->both_sides = 1;
  Lk(tree);

  while(n_round < ROUND_MAX)
    {
      if(tree->mod->s_opt->opt_bl)
	{
	  MC_Optimize_Node_Times_Serie(tree->n_root,tree->n_root->v[0],NULL,tree);
	  MC_Optimize_Node_Times_Serie(tree->n_root,tree->n_root->v[1],NULL,tree);
	  MC_Optimize_Root_Height(tree);
	}

      tree->both_sides = 1;
      Lk(tree);

      if(tree->mod->s_opt->print) Print_Lk(tree,"[Node times         ]");

      lk_new = tree->c_lnL;
      if((!each) ||
	 (fabs(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_global))
	{
	  each = 1;
	  MC_Optimize_Tree_Height(tree);
	  Optimiz_All_Free_Param(tree,tree->mod->s_opt->print);
	  tree->both_sides = 1;
	  Lk(tree);
	}

      lk_new = tree->c_lnL;

      if(lk_new < lk_old - tree->mod->s_opt->min_diff_lk_global*10.) 
	{
	  printf("\n. lk_new = %f lk_old = %f",lk_new,lk_old);
	  Warn_And_Exit("\n. Optimisation failed ! (Round_Optimize)\n");
	}
      if(fabs(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_global)  break;
/*       if(fabs(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_local)  break; */
      else lk_old  = lk_new;
      n_round++;
      each--;
    }
}



/*********************************************************/

void MC_Optimize_Node_Times_Serie(node *a, node *d, edge *b, arbre *tree)
{
  int i;      

  if(d->tax) return;
  else
    {
      node *v1, *v2; /* the two sons of d */
      phydbl t_sup, t_inf;
      phydbl lk_init;
      
      lk_init = tree->c_lnL;
      
      v1 = v2 = NULL;
      For(i,3) if(d->v[i] != a) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
	  
      t_inf = MAX(v1->t,v2->t);
      t_sup = a->t;
	  
      if(t_sup < t_inf - MDBL_MAX)
	{
	  printf("\n. t_sup = %f t_inf = %f",t_sup,t_inf);
	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      else
	{
	  Node_Time_Brent(t_inf,d->t,t_sup,
			  tree->mod->s_opt->min_diff_lk_local,
			  a,d,tree,
			  tree->mod->s_opt->brent_it_max);  
	}

      if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local*10.)
/*       if(tree->c_lnL < lk_init - 1.E-03) */
	{
	  printf("\n. t-inf= %f t-sup=%f t-est=%f",t_inf,t_sup,d->t);
	  printf("\n. %f -- %f",lk_init,tree->c_lnL);
	  printf("\n. a->num = %d, d->num = %d",a->num,d->num);
	  Warn_And_Exit("\n. Err. in MC_Optimize_Node_Times_Serie.");
	}

/*       printf("\n. init_lnL = %f c_lnL = %f",lk_init,tree->c_lnL); */
      
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    Update_P_Lk(tree,d->b[i],d);
	    MC_Optimize_Node_Times_Serie(d,d->v[i],d->b[i],tree);
	  }

      For(i,3) 
	if((d->v[i] == a) || (d->b[i] == tree->e_root)) 
	  {
	    Update_P_Lk(tree,d->b[i],d);
	    break;
	  }
    }
}

/*********************************************************/

void MC_Print_Node_Times(node *a, node *d, arbre *tree)
{
  edge *b;
  int i;
  
  b = NULL;
  For(i,3) if((d->v[i]) && (d->v[i] == a)) {b = d->b[i]; break;}

  printf("\n. (%3d %3d) a->t = %f d->t = %f (#=%f) b->l = %f",a->num,d->num,a->t,d->t,a->t-d->t,(b)?(b->l):(-1.0));
  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  MC_Print_Node_Times(d,d->v[i],tree);
    }
}

/*********************************************************/

edge *MC_Find_Best_Root_Position(arbre *tree)
{
  int i;
  edge *best_edge;
  phydbl best_lnL,best_pos;

  Record_Br_Len(tree);
  best_pos = -1.;
  best_edge = NULL;
  best_lnL = UNLIKELY;
  For(i,2*tree->n_otu-3)
    {
      Restore_Br_Len(tree);
      printf("\n. Root positioned on edge %3d",i);
      MC_Least_Square_Node_Times(tree->t_edges[i],tree);      
      MC_Adjust_Node_Times(tree);
      MC_Round_Optimize(tree);
      if(tree->c_lnL > best_lnL)
	{
	  best_lnL = tree->c_lnL;
	  best_edge = tree->t_edges[i];
	  best_pos = tree->n_root_pos;
	}
    }
  tree->n_root_pos = best_pos; /* Set the root node to its best position */ 
  printf("\n. Best root position: edge %3d",best_edge->num);
  printf("\n. Best constrained lnL = %f",best_lnL);
  fprintf(tree->io->fp_out_stats,"\n. Best constrained lnL = %f",best_lnL);
  return best_edge;
}

/*********************************************************/

edge *MC_Find_Best_Root_Position_Approx(arbre *tree)
{
  int i;
  edge *best_edge;
  phydbl best_lnL,best_pos;

  Record_Br_Len(tree);
  best_pos  = -1.;
  best_edge = NULL;
  best_lnL  = UNLIKELY;
  For(i,2*tree->n_otu-3)
    {
      Restore_Br_Len(tree);
      printf("\n. Root positioned on edge %3d",i);
      MC_Least_Square_Node_Times(tree->t_edges[i],tree);      
      MC_Adjust_Node_Times(tree);
      Lk(tree);
      printf("\n. LnL = %f",tree->c_lnL);
      if(tree->c_lnL > best_lnL)
	{
	  best_lnL  = tree->c_lnL;
	  best_edge = tree->t_edges[i];
	  best_pos  = tree->n_root_pos;
	}
    }
  tree->n_root_pos = best_pos; /* Set the root node to its best position */ 
  printf("\n. Best root position: edge %3d",best_edge->num);
  printf("\n. Best constrained lnL = %f",best_lnL);
  fprintf(tree->io->fp_out_stats,"\n. Best constrained lnL = %f",best_lnL);
  return best_edge;
}

/*********************************************************/

void MC_Optimize_Tree_Height(arbre *tree)
{
  tree->mod->s_opt->tree_size_mult = 1.0;
  Time_Stamps_Mult_Brent(0.1,1.0,10.0,
			 tree->mod->s_opt->min_diff_lk_global,
			 tree,100);
}

/*********************************************************/

void MC_Optimize_Root_Height(arbre *tree)
{
  /*
           t_root
           --x---
l_2 / mu  |      | l_2 / mu
          |      |
          x      x t_r
l_1 / mu  |
          |
          x t_l

 l* = 2 l_2 + l_1 is the ML estimate of the root branch length.
 l_1 = (Max(t_l,t_r) - Min(t_l,t_r)) * mu 
 l_2 = (l* - l1)/2
 t_root = Max(t_l,t_r) + l_2 / mu 

  */
  
  phydbl l_1, l_2;
  phydbl mean_rate, branch_rate;

  mean_rate   = tree->rates->mean_r;
  branch_rate = tree->rates->br_r[tree->e_root->num];

  Br_Len_Brent_Default(tree->e_root,tree);

  l_1 = 
    (MAX(tree->e_root->left->t,tree->e_root->rght->t) -
     MIN(tree->e_root->left->t,tree->e_root->rght->t)) *
    mean_rate * branch_rate;
  
  l_2 = (tree->e_root->l - l_1) / 2.;

  if(l_2 < 0.0) 
    {
      l_2 = 0.0;
      tree->e_root->l = l_1;
      Lk_At_Given_Edge(tree->e_root,tree);
    }

  tree->n_root->t = 
    MAX(tree->e_root->left->t,tree->e_root->rght->t) +
    l_2 / (mean_rate * branch_rate);

/*  /\* Check that the optimal 'root branch' length is longer than the  */
/*      lower bound determined by time stamps on the left and */
/*      right handside of the root branch  */
/*   *\/ */
/*   if(tree->e_root->l <  */
/*      MAX(tree->e_root->left->t,tree->e_root->rght->t) -  */
/*      MIN(tree->e_root->left->t,tree->e_root->rght->t)) */
/*     { */
/*       tree->e_root->l =  */
/* 	MAX(tree->e_root->left->t,tree->e_root->rght->t) -  */
/* 	MIN(tree->e_root->left->t,tree->e_root->rght->t); */
/*       Lk_At_Given_Edge(tree->e_root,tree); */
/*     } */

/*   tree->n_root->t = .5 * (tree->e_root->l + tree->e_root->left->t + tree->e_root->rght->t); */
 
/* /\*   printf("\n. Root->t = %f left->t=%f rght->t=%f e_root->l=%f", *\/ */
/* /\* 		 tree->n_root->t, *\/ */
/* /\* 		 tree->e_root->left->t, *\/ */
/* /\* 		 tree->e_root->rght->t, *\/ */
/* /\* 		 tree->e_root->l); *\/ */

  if((tree->n_root->t < tree->e_root->left->t-1.E-4) ||
     (tree->n_root->t < tree->e_root->rght->t-1.E-4))
    {
      printf("\n. t_root = %f t_left = %f t_rght = %f",
	     tree->n_root->t,
	     tree->e_root->left->t,
	     tree->e_root->rght->t);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
}

/*********************************************************/

void MC_Estimate_Branch_Rate_Parameter(arbre *tree)
{

  /* The tree should be clock-like already */
  Branch_Rate_Shape_Brent(0.3, tree->mod->rr_branch_alpha, 100.,
			  tree->mod->s_opt->min_diff_lk_global, 
			  &(tree->mod->rr_branch_alpha), 
			  tree,tree->mod->s_opt->brent_it_max);  
}

/*********************************************************/

void MC_Compute_Rates_And_Times_Least_Square_Adjustments(arbre *tree)
{
  MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(tree->n_root,tree->n_root->v[0],NULL,tree);
  MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(tree->n_root,tree->n_root->v[1],NULL,tree);
}

/*********************************************************/

void MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(node *a, node *d, edge *b, arbre *tree)
{
  int i;

  if(d->tax) return;

  if(b)
    {
      phydbl t0, t1, t2, t3;
      phydbl mu1, mu2, mu3;
      phydbl K;

      t0 = a->t;
      t1 = d->t;

      mu1 = tree->rates->br_r[b->num];
      mu2 = -1.;
      mu3 = -1.;

      t2 = t3 = -1.;
      For(i,3)
	if(d->v[i] != a)
	  {
	    if(t2 < 0) 
	      {
		t2  = d->v[i]->t;
		mu2 = tree->rates->br_r[d->b[i]->num];
	      }
	    if(t3 < 0) 
	      {
		t3  = d->v[i]->t;
		mu3 = tree->rates->br_r[d->b[i]->num];
	      }
	  }

      printf("\n. t0=%f, t1=%f, t2=%f, t3=%f, mu1=%f, mu2=%f, mu3=%f",
	     t0,t1,t2,t3,
	     mu1,mu2,mu3);

      if(t2 < 1.E-10 && t3 < 1.E-10)
	{
	  K = t0 / t1;
	}
      else
	{
	  K = 
	    (pow(mu2,2)*t2)/(pow(t1-t2,2)) +
	    (pow(mu1,2)*t0)/(pow(t0-t1,2)) +
	    (pow(mu3,2)*t3)/(pow(t1-t3,2)) +
	    (pow(mu1,2)*t0)/(pow(t0-t1,2)) +
	    mu2 * mu1 * (t0 + t2) / ((t1-t2)*(t0-t1)) +
	    mu3 * mu1 * (t0 + t3) / ((t1-t3)*(t0-t1)) ;

	  K /= 
	    t1 *
	    (
	    (pow(mu2,2))/(pow(t1-t2,2))  +
	    (pow(mu1,2))/(pow(t0-t1,2))  +
	    (pow(mu3,2))/(pow(t1-t3,2))  +
	    (pow(mu1,2))/(pow(t0-t1,2))  +
	    2.*mu2*mu1/((t1-t2)*(t0-t1)) +
	    2.*mu3*mu1/((t1-t3)*(t0-t1)) 
	    );
	}
      printf("\n. K = %f",K); 
    }

  For(i,3) 
    if(d->v[i] != a)
      MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(d,d->v[i],d->b[i],tree);

}

/*********************************************************/

void MC_Classify_Branch_Rates(arbre *tree)
{
  int br,rcat,best_post_prob_cat;
  phydbl *post_prob;
  phydbl best_post_prob;
  edge *b;
  
  post_prob = (phydbl *)mCalloc(tree->mod->n_rr_branch,sizeof(phydbl));

  Lk(tree);

  Record_Br_Len(tree);

  For(br,2*tree->n_otu-3)
    {
      b = tree->t_edges[br];
      post_prob = (phydbl *)Post_Prob_Rates_At_Given_Edge(b,post_prob,tree);

      best_post_prob = UNLIKELY;
      best_post_prob_cat = -1;
      For(rcat,tree->mod->n_rr_branch)
	{
	  if(post_prob[rcat] > best_post_prob)
	    {
	      best_post_prob = post_prob[rcat];
	      best_post_prob_cat = rcat;
	    }
	}
      tree->rates->br_r[br] = tree->mod->rr_branch[best_post_prob_cat];
    }

  Free(post_prob);  
}

/*********************************************************/
/*********************************************************/
/*********************************************************/

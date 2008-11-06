/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

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


#ifdef PHYML

int main(int argc, char **argv)
{
  seq **data;
  allseq *alldata;
  option *io;
  arbre *tree;
  int n_otu, num_data_set;
  int num_tree,tree_line_number,num_rand_tree;
  matrix *mat;
  model *mod;
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

#ifndef BATCH
		  if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);
#endif

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

		  if(io->mod->s_opt->spr_step_after_nnis)
		    {
		      printf("\n. The NNI+SPR option is not available yet, sorry.");
		      Warn_And_Exit("");
		    }

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
		  printf("\n\n. Final log likelihood : %f",tree->c_lnL);

		  if(tree->io->ratio_test) aLRT(tree);

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

/*  		  //printing loglk for each site, to compute SH-like tests *\/ */
/*  		 phydbl sum=0.0; */
/*  		  printf("\n\nSITES LKS:\n"); */
/*  		  int n_patterns = (int)floor(tree->n_pattern*tree->prop_of_sites_to_consider); */
/*  		  int site=0; */
/*  		  For(site,n_patterns) { */
/*  		    int wei=0; */
/*  		    For(wei,tree->data->wght[site]) { */
/*  		      printf("%f\n",tree->c_lnL_sorted[site] / tree->data->wght[site]); */
/*  		      sum+=tree->c_lnL_sorted[site] / tree->data->wght[site]; */
/*  		    } */
/*  		  } */

/*  		  printf("\n\nsum=%f\n\n",sum); */
/*  		  int i=0; */
/*  		  For(i,2*tree->n_otu-3) */
/*  		    { */
/*  		      if((!tree->t_edges[i]->left->tax) && (!tree->t_edges[i]->rght->tax)) */
/*  			{ */
/*  			  printf("%3d %f %f %f\n", */
/*  				 tree->t_edges[i]->bip_score,tree->t_edges[i]->alrt_statistic, tree->t_edges[i]->ratio_test,tree->t_edges[i]->l); */
/*  			} */
/*  		    } */


/* 		  //printing loglk for each site, to compute SH-like tests */
/* 		  phydbl sum=0.0; */
/* 		  printf("\n\nSITES LKS:\n"); */
/* 		  int n_patterns = (int)floor(tree->n_pattern*tree->prop_of_sites_to_consider); */
/* 		  int site=0; */
/* 		  For(site,n_patterns) { */
/* 		    int wei=0; */
/* 		    For(wei,tree->data->wght[site]) { */
/* 		      printf("%f\n",tree->c_lnL_sorted[site] / tree->data->wght[site]); */
/* 		      sum+=tree->c_lnL_sorted[site] / tree->data->wght[site]; */
/* 		    } */
/* 		  } */

/* 		  printf("\n\nsum=%f\n\n",sum); */

/* 		  int i=0; */
/* 		  For(i,2*tree->n_otu-3) */
/* 		    { */
/* 		      if((!tree->t_edges[i]->left->tax) && (!tree->t_edges[i]->rght->tax)) */
/* 			{ */
/* 			  printf("%3d %f %f %f\n", */
/* 				 tree->t_edges[i]->bip_score,tree->t_edges[i]->alrt_statistic, tree->t_edges[i]->ratio_test,tree->t_edges[i]->l); */
/* 			} */
/* 		    } */


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

  if(io->fp_in_seq)     fclose(io->fp_in_seq);
  if(io->fp_in_tree)    fclose(io->fp_in_tree);
  if(io->fp_out_lk)     fclose(io->fp_out_lk);
  if(io->fp_out_tree)   fclose(io->fp_out_tree);
  if(io->fp_out_stats)  fclose(io->fp_out_stats);


  Free_Input(io);
  return 0;
}

#elif(MG)
#include "mg.h"
int main(int argc, char **argv)
{
  MC_main(argc, argv);
  return 1;
}

#elif(M4)
#include "m4.h"
int main(int argc, char **argv)
{
  MC_main(argc, argv);
  return 1;
}

#elif(MC)
#include "mc.h"
int main(int argc, char **argv)
{
  MC_main(argc, argv);
  return 1;
}

#endif

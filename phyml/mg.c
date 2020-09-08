/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/


#include "mg.h"
#include "free.h"
#include "options.h"
#include "utilities.h"
#include "optimiz.h"
#include "models.h"
#include "simu.h"
#include "lk.h"
#include "pars.h"
#include "interface.h"

/*********************************************************/

int MG_main(int argc, char **argv)
{
  option *io;
  char *s_tree;
  FILE *fp_phyml_tree,*fp_phyml_stats,*fp_phyml_lk;
  int set;
  time_t t_beg,t_end;
  div_t hour,min;
  int r_seed;
  int i;

  fflush(NULL);
  io = (option *)Get_Input(argc,argv);
  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
  srand(r_seed);
  Make_Model_Complete(io->mod);
  fp_phyml_stats = Openfile(io->out_stats_file,io->out_stats_file_open_mode);
  //fprintf(fp_phyml_stats,"\n- PHYML %s -\n\n", VERSION);
  fp_phyml_tree = Openfile(io->out_tree_file,io->out_tree_file_open_mode);
  fp_phyml_lk = fopen(io->out_lk_file,"w");

  time(&t_beg);

  if(io->multigene)
    {
      seq ***data;
      allseq **alldata;
      model **mod;
      matrix **mat;
      arbrelist *treelist;
      superarbre *st;

      data      = (seq ***)  mCalloc(io->n_gt,sizeof(seq **));
      alldata   = (allseq **)mCalloc(io->n_gt,sizeof(allseq *));
      mod       = (model **) mCalloc(io->n_gt,sizeof(model *));
      mat       = (matrix **)mCalloc(io->n_gt,sizeof(matrix *));



      /* Read the sequences (for each subdataset) */
      For(set,io->n_gt)
	{
	  Make_Model_Complete(io->st->optionlist[set]->mod); /* Complete model for each data set */
	  mod[set]  = io->st->optionlist[set]->mod;
	  data[set] = Get_Seq(io->st->optionlist[set],0);
	  printf("\n. Data set [#%d]\n",set+1);
	  printf("\n. Compressing sequences...\n");
	  alldata[set] = Compact_Seq(data[set],io->st->optionlist[set]);
	  fclose(io->st->optionlist[set]->fp_in_seq);
	  Free_Seq(data[set],alldata[set]->n_otu);
	  Init_Model(alldata[set],mod[set]);
	  Check_Ambiguities(alldata[set],
			    io->st->optionlist[set]->mod->datatype,
			    io->st->optionlist[set]->mod->stepsize);
	}

      Mg_Make_Superarbre_Full(io->st,io,alldata);
      st = io->st;
      treelist = st->treelist;
      Fill_Dir_Table(st->tree);
      Update_Dirs(st->tree);

      For(set,io->n_gt)
	{
	  st->data_of_interest = alldata[set];
	  if(!Mg_Get_Species_Found_In_St(st,alldata[set])) break;
	  treelist->tree[set] = Make_Tree_From_Scratch(st->tree->n_otu,NULL);
	  Copy_Tree_Topology_With_Labels(st->tree,treelist->tree[set]);
 	  treelist->tree[set]->num_curr_branch_available = 0;
	  Connect_Edges_To_Nodes_Recur(treelist->tree[set]->noeud[0],
				       treelist->tree[set]->noeud[0]->v[0],
				       treelist->tree[set]);
	  Mg_Prune_St_Topo(treelist->tree[set],alldata[set],st);

	  if(treelist->tree[set]->n_otu != alldata[set]->n_otu)
	    {
	      printf("\n. Problem with sequence file '%s'\n",io->st->optionlist[set]->in_seq_file);
	      printf("\n. # taxa found in supertree restricted to '%s' taxa = %d\n",
		     io->st->optionlist[set]->in_seq_file,
		     treelist->tree[set]->n_otu);
	      printf("\n. # sequences in '%s' = %d\n",
		     io->st->optionlist[set]->in_seq_file,
		     alldata[set]->n_otu);
	      Exit("\n");
	    }

	  treelist->tree[set]->dp         = set;
	  treelist->tree[set]->n_otu      = alldata[set]->n_otu;
	  treelist->tree[set]->mod        = mod[set];
	  treelist->tree[set]->io         = io->st->optionlist[set];
	  treelist->tree[set]->data       = alldata[set];
	  treelist->tree[set]->both_sides = 1;
	  treelist->tree[set]->n_pattern  = treelist->tree[set]->data->crunch_len/
	                                     treelist->tree[set]->mod->stepsize;

	  Order_Tree_CSeq(treelist->tree[set],alldata[set]);
	  Fill_Dir_Table(treelist->tree[set]);
	  Update_Dirs(treelist->tree[set]);
	  Make_Tree_4_Lk(treelist->tree[set],alldata[set],alldata[set]->init_len);
	  Make_Tree_4_Pars(treelist->tree[set],alldata[set],alldata[set]->init_len);
	  treelist->tree[set]->triplet_struct = Make_Triplet_Struct(treelist->tree[set]->mod);
	}

      if(set != io->n_gt)
	{
	  printf("\n. Sequence data set found in '%s' has one or more taxa not found in the '%s' tree file\n",
		 io->st->optionlist[set]->in_seq_file,
		 io->in_tree_file);
	  Exit("\n");
	}

      Mg_Check_Extra_Taxa(st);
	
      st->tree->c_lnL = .0;
      For(set,io->n_gt)
	{
	  Lk(treelist->tree[set]);
	  Get_List_Of_Reachable_Tips(treelist->tree[set]);
	  Mg_Match_St_Nodes_In_Gt(treelist->tree[set],st);
	  Mg_Match_St_Edges_In_Gt(treelist->tree[set],st);
	  Mg_Map_St_Nodes_In_Gt(treelist->tree[set],st);
	  Mg_Map_St_Edges_In_Gt(treelist->tree[set],st);
	  Mg_Map_Gt_Edges_In_St(treelist->tree[set],st);
	  st->tree->c_lnL += treelist->tree[set]->c_lnL;
	}

      Mg_Initialise_Bl_Partition(st);
      Mg_Set_Bl(st->bl,st);

      time(&(st->tree->t_beg));
      time(&(st->tree->t_current));
      printf("\n. (%5d sec) [00] [%10.2f] [%5d]\n",
	     (int)(st->tree->t_current-st->tree->t_beg),
	     Mg_Lk(st),
	     Mg_Pars(st));

      
/*       int n_iter=0; */
/*       do */
/* 	{ */
/* 	  Mg_Optimize_Br_Len_Serie(st->tree->noeud[0], */
/* 				   st->tree->noeud[0]->v[0], */
/* 				   st->tree->noeud[0]->b[0], */
/* 				   st); */

/* 	  st->tree->both_sides = 1; */
/* 	  Mg_Lk(st); */
/* 	  printf("\n. %f",st->tree->c_lnL); */
/* /\* 	  For(set,st->n_gt) printf("\n. %s",Write_Tree(st->treelist->tree[set])); *\/ */
/* 	  n_iter++; */
/* 	}while(n_iter < 5); */


/*       Mg_Lk(st); */
/*       printf("\n> %f",st->tree->c_lnL); */
/*       For(i,2*st->tree->n_otu-3) */
/* 	{ */
/* 	  if((!st->tree->t_edges[i]->left->tax) && (!st->tree->t_edges[i]->rght->tax)) */
/* 	    { */
/* 	      Mg_NNI(st->tree->t_edges[i],st); */
/* 	    } */
/* 	} */

      
/*       if(io->mod->s_opt->topo_search == NNI_MOVE) */
      Mg_Simu(st);
/*       else */
/*       Mg_Speed_Spr(st); */


      time(&t_end);

      hour = div(t_end-t_beg,3600);
      min  = div(t_end-t_beg,60  );

      min.quot -= hour.quot*60;


      fprintf(fp_phyml_stats,"\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
      fprintf(fp_phyml_stats,"\n. Number of partitions = %d\n\n",st->n_gt);
      fprintf(fp_phyml_stats,"\n. Full data set -- lnL = %f\n\n",st->tree->c_lnL);
      fprintf(fp_phyml_stats,"\n. Tree search algorithm : %s\n\n",(io->mod->s_opt->topo_search == NNI_MOVE)?("NNIs"):("SPRs"));
      fprintf(fp_phyml_stats,"\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n\n");
      For(set,io->n_gt)
	{
	  Print_Fp_Out(fp_phyml_stats,t_beg,t_end,st->treelist->tree[set],
		       io->st->optionlist[set],1,1);
	}


      printf("\n\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
      printf("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");

      For(i,2*st->tree->n_otu-3) st->tree->t_edges[i]->l = 0.1;
      s_tree = Write_Tree(st->tree);
      fprintf(fp_phyml_tree,"Supertree\n");
      fprintf(fp_phyml_tree,"%s\n",s_tree);
      Free(s_tree);
      For(set,st->n_gt)
	{
	  fprintf(fp_phyml_tree,"Gene tree number %d\n",set+1);
	  s_tree = Write_Tree(st->treelist->tree[set]);
	  fprintf(fp_phyml_tree,"%s\n",s_tree);
	  Free(s_tree);
	}


      For(set,st->n_gt)
	{
	  if(io->mod->s_opt->topo_search == SPR_MOVE) Free_Spr_List(treelist->tree[set]);
	  Free_Tree_Lk(treelist->tree[set]);
	  Free_Tree_Pars(treelist->tree[set]);
	  Free_Triplet(treelist->tree[set]->triplet_struct);
	  Free_Tree(treelist->tree[set]);
	  Free_Cseq(alldata[set]);
	  Free_Model(mod[set]);
	  Free_Input(io->st->optionlist[set]);

	}

      if(io->mod->s_opt->topo_search == SPR_MOVE) Free_Spr_List(st->tree);
      Free(mat);
      Free(mod);
      Free(data);
      Free(alldata);
      Free(treelist->tree);
      Free(treelist);
      Free_St(st);
    }


  if(io->fp_in_seq ) fclose(io->fp_in_seq);
  if(io->fp_in_tree) fclose(io->fp_in_tree);


  Free_Model(io->mod);
  Free_Input(io);

  fclose(fp_phyml_lk);
  fclose(fp_phyml_tree);
  fclose(fp_phyml_stats);


  return 0;
}

/*********************************************************/

void Mg_Print_Nodes(node *a, node *d, superarbre *st)
{
  int i;
  printf(">>>>>>>>>>>>>>>>>>>>\n");
  printf("num node = %d\n",d->num);
  if(d->tax) printf("name='%s'\n",d->name);
  else
    {
      printf("n_of_reachable_tips : \n");
      For(i,3)
	{
/* 	  printf("dir%d=%d; ",i,st->n_of_reachable_tips[st->num_data_of_interest][d->num][i]); */
/* 	  For(j,st->n_of_reachable_tips[st->num_data_of_interest][d->num][i]) */
/* 	    { */
/* 	      printf("%s ", */
/* 		     st->list_of_reachable_tips[st->num_data_of_interest][d->num][i][j]->name); */
/* 	    } */
	  
	  printf("\n");
	}
    }
  printf("<<<<<<<<<<<<<<<<<<<\n\n");
  if(d->tax) return;
  else
    {
      For(i,3) if(d->v[i] != a) Mg_Print_Nodes(d,d->v[i],st);
    }
}

/*********************************************************/

superarbre *Mg_Make_Superarbre_Light(option *input)
{
  superarbre *st;

  st = (superarbre *)mCalloc(1,sizeof(superarbre));
  st->optionlist = (option **)mCalloc(input->n_gt,sizeof(option *));
  st->bl_partition = (int *)mCalloc(input->n_gt,sizeof(int ));
  return st;
}

/*********************************************************/

void Mg_Make_Superarbre_Full(superarbre *st, option *io, allseq **data)
{
  int i,j,k;

  if(io->in_tree)
    {
      printf("\n. Reading user tree...\n");
      rewind(io->fp_in_tree);
      
      st->tree = Read_Tree_File(io->fp_in_tree);
      
      if(!st->tree->has_branch_lengths)
	{
	  printf("\n. Branch lengths are all set to 0.1...\n");
	  For(i,2*st->tree->n_otu-3) st->tree->t_edges[i]->l = 0.1;
	}
    }
  else
    {
      Warn_And_Exit("\n. A user-defined input tree is needed\n");
    }

  st->tree->io      = io;
  st->treelist      = (arbrelist *)Make_Treelist(io->n_gt);
  st->n_gt          = io->n_gt;
  st->tree->mod     = io->mod;
  st->lock_br_len   = 0;

  st->map_st_node_in_gt = (node *****)mCalloc(st->n_gt,sizeof(node ****));
  For(i,st->n_gt) 
    {
      st->map_st_node_in_gt[i] = (node ****)mCalloc(2*st->tree->n_otu-2,sizeof(node ***));
      For(j,2*st->tree->n_otu-2) 
	{
	  st->map_st_node_in_gt[i][j] = (node ***)mCalloc(3,sizeof(node **));
	  For(k,3) st->map_st_node_in_gt[i][j][k] = (node **)mCalloc(2,sizeof(node *));
	}
    }

  st->map_st_edge_in_gt = (edge ***)mCalloc(st->n_gt,sizeof(edge **));
  For(i,st->n_gt) st->map_st_edge_in_gt[i] = (edge **)mCalloc(2*st->tree->n_otu-3,sizeof(edge *));

  st->map_gt_edge_in_st = (edge ****)mCalloc(st->n_gt,sizeof(edge ***));
  For(i,st->n_gt)
    {
      st->map_gt_edge_in_st[i] = (edge ***)mCalloc(2*st->tree->n_otu-3,sizeof(edge **));
      For(j,2*st->tree->n_otu-3) st->map_gt_edge_in_st[i][j] = (edge **)mCalloc(2*st->tree->n_otu-3,sizeof(edge *));
    }

  st->size_map_gt_edge_in_st = (int **)mCalloc(st->n_gt,sizeof(int *));
  For(i,st->n_gt) st->size_map_gt_edge_in_st[i] = (int *)mCalloc(2*st->tree->n_otu-3,sizeof(int));


  st->match_st_edge_in_gt = (edge ***)mCalloc(st->n_gt,sizeof(edge **));
  For(i,st->n_gt) st->match_st_edge_in_gt[i] = (edge **)mCalloc(2*st->tree->n_otu-3,sizeof(edge *));

  st->match_gt_edge_in_st = (edge ***)mCalloc(st->n_gt,sizeof(edge **));
  For(i,st->n_gt) st->match_gt_edge_in_st[i] = (edge **)mCalloc(2*st->tree->n_otu-3,sizeof(edge *));

  st->bl = (phydbl **)mCalloc(st->n_gt,sizeof(phydbl *));
  For(i,st->n_gt) st->bl[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl));

  st->bl_cpy = (phydbl **)mCalloc(st->n_gt,sizeof(phydbl *));
  For(i,st->n_gt) st->bl_cpy[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl));

  st->bl0 = (phydbl **)mCalloc(st->n_gt,sizeof(phydbl *));
  For(i,st->n_gt) st->bl0[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl));

  st->bl1 = (phydbl **)mCalloc(st->n_gt,sizeof(phydbl *));
  For(i,st->n_gt) st->bl1[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl));

  st->bl2 = (phydbl **)mCalloc(st->n_gt,sizeof(phydbl *));
  For(i,st->n_gt) st->bl2[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl));

  st->s_mod = (model **)mCalloc(st->n_gt,sizeof(model *));

  For(i,2*st->tree->n_otu-3) Make_Edge_NNI(st->tree->t_edges[i]);

  Get_List_Of_Reachable_Tips(st->tree);

  st->match_st_node_in_gt = (node ***)mCalloc(io->n_gt,sizeof(node **));
  For(i,io->n_gt) st->match_st_node_in_gt[i] = (node **)mCalloc(2*st->tree->n_otu-2,sizeof(node *));
}

/*********************************************************/

void Mg_Prune_St_Topo(arbre *tree, allseq *data, superarbre *st)
{
  int i,j,not_found;
  int curr_ext_node, curr_int_node, curr_br, n_pruned_nodes;;
  node **pruned_nodes;
  edge **residual_edges;

  pruned_nodes   = (node **)mCalloc(st->tree->n_otu,sizeof(node *));
  residual_edges = (edge **)mCalloc(st->tree->n_otu,sizeof(edge *));

  n_pruned_nodes = 0;
  For(i,st->tree->n_otu)
    {
      For(j,data->n_otu)
	{
	  if(!strcmp(data->c_seq[j]->name,st->tree->noeud[i]->name))
	    break;
	}

      not_found = 1;
      if(j == data->n_otu)
	{
	  For(j,tree->n_otu)
	    {
	      if(!strcmp(tree->noeud[j]->name,st->tree->noeud[i]->name))
		{
		  Prune_Subtree(tree->noeud[j]->v[0],
				tree->noeud[j],
				NULL,&(residual_edges[n_pruned_nodes]),
				tree);

		  pruned_nodes[n_pruned_nodes] = tree->noeud[j];
		  n_pruned_nodes++;
		  not_found = 0;
		  break;
		}	      
	    }


	  if(not_found)	    
	    {
	      printf("\n. Taxon '%s'",st->tree->noeud[i]->name);
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	}
    }

  For(i,2*tree->n_otu-2) Free(tree->t_dir[i]); Free(tree->t_dir);

  tree->n_otu -= n_pruned_nodes;

  curr_ext_node = 0;
  curr_int_node = tree->n_otu;  
  curr_br = 0;
  For(i,st->tree->n_otu)
    {
      For(j,n_pruned_nodes)
	{
	  if(!strcmp(pruned_nodes[j]->name,st->tree->noeud[i]->name))
	    break;
	}
      if(j == n_pruned_nodes) /* That node still belongs to the tree */
	{
	  Reassign_Node_Nums(tree->noeud[i],tree->noeud[i]->v[0], 
			     &curr_ext_node, &curr_int_node,tree);
	  break;
	}
    }
  
  Reassign_Edge_Nums(tree->noeud[0],tree->noeud[0]->v[0],&curr_br,tree);

  tree->t_dir = (int **)mCalloc(2*tree->n_otu-2,sizeof(int *));
  For(i,2*tree->n_otu-2) tree->t_dir[i] = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));

  For(i,n_pruned_nodes) 
    {
      Free_Edge(residual_edges[i]);
      Free_Edge(pruned_nodes[i]->b[0]);
      Free_Node(pruned_nodes[i]->v[0]);
      Free_Node(pruned_nodes[i]);
    }

  Free(pruned_nodes);
  Free(residual_edges);

}

/*********************************************************/

void Mg_Match_St_Nodes_In_Gt(arbre *gt, superarbre *st)
{
  int i,j;

  For(i,2*st->tree->n_otu-2) st->match_st_node_in_gt[gt->dp][i] = NULL; /* don't forget that step ! */

  /* Map tips */
  For(i,st->tree->n_otu)
    {
      For(j,gt->n_otu)
	{
	  if(!strcmp(st->tree->noeud[i]->name,gt->noeud[j]->name))
	    {
	      st->match_st_node_in_gt[gt->dp][st->tree->noeud[i]->num] = gt->noeud[j];
	      break;
	    }
	}
    }

#ifdef DEBUG
  /* Checking that the results are correct so far */
  int n_matches;
  n_matches = 0;
  For(i,2*st->tree->n_otu-2)
    if(st->match_st_node_in_gt[gt->dp][i])
      n_matches++;

  if(n_matches != gt->n_otu)
    {
      printf("\n");
      printf("\n. n_matches = %d 2*gt->n_otu-2 = %d\n",n_matches,2*gt->n_otu-2);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif


  /* Map internal nodes */
  For(i,st->tree->n_otu)
    {
      if(st->match_st_node_in_gt[gt->dp][st->tree->noeud[i]->num])
	{
	  Mg_Match_St_Nodes_In_Gt_Recurr(st->match_st_node_in_gt[gt->dp][st->tree->noeud[i]->num],
					 st->match_st_node_in_gt[gt->dp][st->tree->noeud[i]->num]->v[0],
					 st->tree->noeud[i],
					 st->tree->noeud[i]->v[0],
					 gt,
					 st);
	  break;
	}
    }
  


#ifdef DEBUG
  /* Checking that the results are correct */
  n_matches = 0;
  For(i,2*st->tree->n_otu-2) 
    if(st->match_st_node_in_gt[gt->dp][st->tree->noeud[i]->num])
	n_matches++;

  if(n_matches != 2*gt->n_otu-2)
    {
      int j;
      printf("\n");
      printf("\n. n_matches = %d 2*gt->n_otu-2 = %d\n",n_matches,2*gt->n_otu-2);
      For(j,2*gt->n_otu-2)
	{
	  For(i,2*st->tree->n_otu-2) 
	    if(st->match_st_node_in_gt[gt->dp][i] == gt->noeud[j])
	      break;

 	  if(i == 2*st->tree->n_otu-2)
	    {
	      printf("\n. Gt %3d node %3d (%3d %3d %3d) (%s %s %s) (%f %f %f) does not match\n",
		     gt->dp,
		     gt->noeud[j]->num,
		     gt->noeud[j]->v[0] ? gt->noeud[j]->v[0]->num : -1,
		     gt->noeud[j]->v[1] ? gt->noeud[j]->v[1]->num : -1,
		     gt->noeud[j]->v[2] ? gt->noeud[j]->v[2]->num : -1,
		     gt->noeud[j]->v[0]->tax ? gt->noeud[j]->v[0]->name : NULL,
		     gt->noeud[j]->v[1]->tax ? gt->noeud[j]->v[1]->name : NULL,
		     gt->noeud[j]->v[2]->tax ? gt->noeud[j]->v[2]->name : NULL,
		     gt->noeud[j]->v[0] ? gt->noeud[j]->b[0]->l : -1.,
		     gt->noeud[j]->v[1] ? gt->noeud[j]->b[1]->l : -1.,
		     gt->noeud[j]->v[2] ? gt->noeud[j]->b[2]->l : -1.);
	    }
	}

      printf("oooooooo\n");
      Print_Node(st->tree->noeud[0],
		 st->tree->noeud[0]->v[0],
		 st->tree);
      printf(">>>>>>>\n");
      For(i,st->n_gt)
	{
	  Print_Node(st->treelist->tree[i]->noeud[0],
		     st->treelist->tree[i]->noeud[0]->v[0],
		     st->treelist->tree[i]);
	  printf("<<<<<<<\n");
	}
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif
}

/*********************************************************/

void Mg_Match_St_Nodes_In_Gt_Recurr(node *a_gt, node *d_gt, node *a_st, node *d_st, arbre *gt, superarbre *st)
{
  int i,j,k;
  int *score_d_st;


  if((d_gt->tax) || (d_st->tax)) return;
  else
    {
      score_d_st = (int *)mCalloc(3,sizeof(int));
      
      For(i,3)
	{
	  For(j,3)
	    {
	      For(k,d_st->n_of_reachable_tips[j])
		{
		  if(!strcmp(d_gt->list_of_reachable_tips[i][0]->name,d_st->list_of_reachable_tips[j][k]->name))
		    break;
		}
	      if(k != d_st->n_of_reachable_tips[j]) score_d_st[j] += 1;
	    }
	}


      if((score_d_st[0] == 1) && (score_d_st[1] == 1) && (score_d_st[2] == 1))
	{
	  st->match_st_node_in_gt[gt->dp][d_st->num] = d_gt;

	  For(i,3)
	    {
	      if(d_gt->v[i] != a_gt)
		{
		  For(j,3)
		    {
		      For(k,d_st->n_of_reachable_tips[j])
			if(!strcmp(d_gt->list_of_reachable_tips[i][0]->name,
				   d_st->list_of_reachable_tips[j][k]->name))
			  {
			    Mg_Match_St_Nodes_In_Gt_Recurr(d_gt,d_gt->v[i],d_st,d_st->v[j],gt,st);
			    break;
			  }
		      if(k != d_st->n_of_reachable_tips[j]) break;
		    }
		}
	    }
	}
      else
	{
	  For(i,3)
	    if(d_st->v[i] != a_st)
	      Mg_Match_St_Nodes_In_Gt_Recurr(a_gt,d_gt,d_st,d_st->v[i],gt,st);
	}
      Free(score_d_st);	
    }
}

/*********************************************************/

void Mg_Match_St_Edges_In_Gt(arbre *gt, superarbre *st)
{
  int i;

  For(i,2*st->tree->n_otu-3) 
    {
      st->match_st_edge_in_gt[gt->dp][i] = NULL; 
      st->match_gt_edge_in_st[gt->dp][i] = NULL;
    }

  For(i,st->tree->n_otu) 
    if(st->match_st_node_in_gt[gt->dp][i])
      {
	Mg_Match_St_Edges_In_Gt_Recurr(st->match_st_node_in_gt[gt->dp][i],
				       st->match_st_node_in_gt[gt->dp][i]->v[0],
				       st->tree->noeud[i],
				       st->tree->noeud[i]->v[0],
				       gt,st);
	break;
      }

}

/*********************************************************/

void Mg_Match_St_Edges_In_Gt_Recurr(node *a_gt, node *d_gt, node *a_st, node *d_st, arbre *gt, superarbre *st)
{
  edge *b_gt, *b_st;
  int i,j,k;

  b_gt = b_st = NULL;

  if((st->match_st_node_in_gt[gt->dp][a_st->num] == a_gt) &&
     (st->match_st_node_in_gt[gt->dp][d_st->num] == d_gt))
    {
      For(i,3) if((a_st->v[i]) && (a_st->v[i] == d_st)) {b_st = a_st->b[i]; break;}
      For(i,3) if((a_gt->v[i]) && (a_gt->v[i] == d_gt)) {b_gt = a_gt->b[i]; break;}

      st->match_st_edge_in_gt[gt->dp][b_st->num] = b_gt;
      st->match_gt_edge_in_st[gt->dp][b_gt->num] = b_st;
    }


  if(!d_gt)
    {
      printf("\n");
      printf("\n. a_gt->num = %d\n",a_gt->num);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(d_gt->tax || d_st->tax) return;
  else
    {
      if(st->match_st_node_in_gt[gt->dp][d_st->num] == d_gt)
	{
	  For(i,3)
	    {
	      if(d_gt->v[i] != a_gt)
		{
		  For(j,3)
		    {
		      For(k,d_st->n_of_reachable_tips[j])
			if(!strcmp(d_gt->list_of_reachable_tips[i][0]->name,d_st->list_of_reachable_tips[j][k]->name))
			  {
			    Mg_Match_St_Edges_In_Gt_Recurr(d_gt,d_gt->v[i],d_st,d_st->v[j],gt,st);
			    break;
			  }
		      if(k != d_st->n_of_reachable_tips[j]) break;
		    }
		}
	    }
	}
      else
	{
	  For(i,3)
	    if(d_st->v[i] != a_st)
	      Mg_Match_St_Edges_In_Gt_Recurr(a_gt,d_gt,d_st,d_st->v[i],gt,st);
	}
    }
}

/*********************************************************/

void Mg_Simu(superarbre *st)
{
  int i,j,step,n_without_swap,it_lim_without_swap;
  edge **sorted_b,*st_b,**tested_b;
  int n_neg,n_tested,each;
  phydbl lambda,old_loglk;

  sorted_b = (edge **)mCalloc(st->tree->n_otu-3,sizeof(edge *));
  tested_b = (edge **)mCalloc(st->tree->n_otu-3,sizeof(edge *));

  For(i,st->n_gt) Update_Dirs(st->treelist->tree[i]);
  Update_Dirs(st->tree);

  each                = 4;
  step                = 0;
  lambda              = .75;
  n_tested            = 0;
  n_without_swap      = 0;
  old_loglk           = UNLIKELY; 
  it_lim_without_swap = 2;
  st_b                = NULL; 
  do
    {
      For(i,st->n_gt) Check_Dirs(st->treelist->tree[i]);

      ++step;      
      each--;

      Mg_Print_Bl(st);

      /* Compute the likelihood of the supertreee */
      st->tree->c_lnL  = Mg_Lk(st);
      st->tree->c_pars = Mg_Pars(st);
/*       For(i,st->n_gt) printf("\n. %s",Write_Tree(st->treelist->tree[i])); */
/*       printf("\n"); */

      time(&(st->tree->t_current));
      printf("\n. (%5d sec) [tot lnL=%15.5f] [# swaps=%3d]",
	     (int)(st->tree->t_current-st->tree->t_beg),
	     st->tree->c_lnL,n_tested);
/*       For(i,st->n_gt) printf("\n[gt %3d lnL=%15.5f]",i,st->treelist->tree[i]->c_lnL); */
      
      if((fabs(old_loglk-st->tree->c_lnL) < st->tree->mod->s_opt->min_diff_lk_global) || 
	 (n_without_swap > it_lim_without_swap)) break;

      if(st->tree->c_lnL < old_loglk)
	{
	  printf("\n. Moving backward (topology + branch lengths) \n");
	  
	  if(!Mg_Mov_Backward_Topo_Bl(st,old_loglk,tested_b,n_tested))
	    Warn_And_Exit("\n. Err: mov_back failed\n");

	  if(!st->tree->n_swap) n_neg = 0;
	    
	  Mg_Record_Br_Len(st);
	  For(i,st->n_gt) Optimiz_All_Free_Param(st->treelist->tree[i],0);
	}
      else 
	{
	  if(!each)
	    {
	      each = 4;
	      /* Markov model parameters are free to vary across data partitions */
	      For(i,st->n_gt) Optimiz_All_Free_Param(st->treelist->tree[i],0);	      
	      For(i,st->n_gt) st->treelist->tree[i]->both_sides = 1;
	      st->tree->c_lnL  = Mg_Lk(st);
	      st->tree->c_pars = Mg_Pars(st);
	    }
	  
	  old_loglk = st->tree->c_lnL;
	  

	  For(i,2*st->tree->n_otu-3) Init_NNI(st->tree->t_edges[i]->nni);

	  /* Test NNIs */
	  For(i,2*st->tree->n_otu-3)
	    {
	      st_b = st->tree->t_edges[i];
	      if((!st_b->left->tax) && (!st_b->rght->tax)) Mg_NNI(st_b,st);
	    }
	  
	  /* Optimise external branch lengths */
	  For(i,2*st->tree->n_otu-3)
	    {
	      st_b = st->tree->t_edges[i];
	      if((st_b->left->tax) || (st_b->rght->tax))
		{
		  Mg_Record_Br_Len(st);
		  Mg_Br_Len_Brent(st_b,st);
		  For(j,st->n_gt) st->bl0[st->bl_partition[j]][st_b->num] = st->bl[st->bl_partition[j]][st_b->num];
		  st_b->nni->score     = .0;
		  st_b->nni->best_conf =  0;
		  Mg_Restore_Br_Len(st);
		  Mg_Lk_At_Given_Edge(st_b,st);
		}
	    }

	  /* Select and sort swaps */
	  n_neg = 0;
	  Select_Edges_To_Swap(st->tree,sorted_b,&n_neg);
	  Sort_Edges_NNI_Score(st->tree,sorted_b,n_neg);	  

	  n_tested = 0;
	  For(i,(int)ceil((phydbl)n_neg*(lambda)))
	    tested_b[n_tested++] = sorted_b[i];

	  if(n_tested > 0) n_without_swap = 0;
	  else             n_without_swap++;

	  Mg_Record_Br_Len(st);
	  
	  /* Apply swaps */
	  Mg_Make_N_Swap(tested_b,0,n_tested,st);

	  /* Update branch lengths (all edges first and then swaped edges) */
	  Mg_Update_Bl(lambda,st);
	  Mg_Update_Bl_Swaped(tested_b,n_tested,st);

	    
	}
    }
  while(1);
  
  printf("\n\n. End of Mg_Simu \n");
  Free(sorted_b);
  Free(tested_b);

  if((n_without_swap > it_lim_without_swap))
    {
      printf("\n. Last optimization step...\n");
      For(i,st->n_gt) Round_Optimize(st->treelist->tree[i],st->treelist->tree[i]->data);
      st->tree->c_lnL = Mg_Lk(st);
      Mg_Simu(st);
    }
}


/*********************************************************/

int Mg_Mov_Backward_Topo_Bl(superarbre *st, phydbl lk_old, edge **tested_b, int n_tested)
{
  int i,j,step,beg,end;
  edge *st_b;
  phydbl **l_init;

  l_init = (phydbl **)mCalloc(st->n_gt,sizeof(phydbl *));
  For(i,st->n_gt) l_init[i] = (phydbl *)mCalloc(2*st->tree->n_otu-3,sizeof(phydbl ));

  For(i,2*st->tree->n_otu-3) For(j,st->n_gt) l_init[st->bl_partition[j]][i] = st->bl[st->bl_partition[j]][i];

  step = 2;
  do
    {
      For(i,2*st->tree->n_otu-3)
	{
	  For(j,st->n_gt)
	    {
	      st->bl[st->bl_partition[j]][i] = st->bl_cpy[st->bl_partition[j]][i] + 
		(1./step) * (l_init[st->bl_partition[j]][i] - st->bl_cpy[st->bl_partition[j]][i]);
/* 	      st->bl[st->bl_partition[j]][i] = st->bl_cpy[st->bl_partition[j]][i]; */
	    }
	}
      
      beg = (int)floor((phydbl)n_tested/(step-1));
      end = 0;
      st_b = NULL;

      for(i=beg-1;i>=end;i--)
	{
	  st_b = tested_b[i];

	  Mg_Swap(st_b->nni->swap_node_v2->v[st->tree->t_dir[st_b->nni->swap_node_v2->num][st_b->nni->swap_node_v1->num]],
		  st_b->nni->swap_node_v2,
		  st_b->nni->swap_node_v3,
		  st_b->nni->swap_node_v3->v[st->tree->t_dir[st_b->nni->swap_node_v3->num][st_b->nni->swap_node_v4->num]],
		  st);

	  Swap(st_b->nni->swap_node_v2->v[st->tree->t_dir[st_b->nni->swap_node_v2->num][st_b->nni->swap_node_v1->num]],
	       st_b->nni->swap_node_v2,
	       st_b->nni->swap_node_v3,
	       st_b->nni->swap_node_v3->v[st->tree->t_dir[st_b->nni->swap_node_v3->num][st_b->nni->swap_node_v4->num]],
	       st->tree);

	  Mg_Do_Mapping(st);

	}

      beg = 0;
      end = (int)floor((phydbl)n_tested/step);
      st_b = NULL;

      for(i=beg;i<end;i++)
	{
	  st_b = tested_b[i];
	  
	  Mg_Swap(st_b->nni->swap_node_v2->v[st->tree->t_dir[st_b->nni->swap_node_v2->num][st_b->nni->swap_node_v1->num]],
		  st_b->nni->swap_node_v2,
		  st_b->nni->swap_node_v3,
		  st_b->nni->swap_node_v3->v[st->tree->t_dir[st_b->nni->swap_node_v3->num][st_b->nni->swap_node_v4->num]],
		  st);

	  Swap(st_b->nni->swap_node_v2->v[st->tree->t_dir[st_b->nni->swap_node_v2->num][st_b->nni->swap_node_v1->num]],
	       st_b->nni->swap_node_v2,
	       st_b->nni->swap_node_v3,
	       st_b->nni->swap_node_v3->v[st->tree->t_dir[st_b->nni->swap_node_v3->num][st_b->nni->swap_node_v4->num]],
	       st->tree);

	  Mg_Do_Mapping(st);

	}
     
      if(!end) st->tree->n_swap = 0;      

      Mg_Lk(st);

      printf("\n. lnL = %15.5f",st->tree->c_lnL);
      step++;
    }
  while((st->tree->c_lnL < lk_old) && (step < 100));

  if(step == 100)
    {
      For(i,2*st->tree->n_otu-3) For(j,st->n_gt)
	st->bl[st->bl_partition[j]][i] = st->bl_cpy[st->bl_partition[j]][i];
    }

  st->tree->n_swap = 0;
  For(i,2*st->tree->n_otu-3) 
    {
      if(st->tree->t_edges[i]->nni->score < 0.0) st->tree->n_swap++;
      st->tree->t_edges[i]->nni->score = +1.0;
    }

  Mg_Lk(st);

  if(st->tree->c_lnL > lk_old)                                                    return  1;
  else if(fabs(st->tree->c_lnL-lk_old) < st->tree->mod->s_opt->min_diff_lk_local) return -1;
  else                                                                            return  0;

  For(i,st->n_gt) Free(l_init[i]);
  Free(l_init);
}

/*********************************************************/

void Mg_Check_Extra_Taxa(superarbre *st)
{
  int i,j,k;
  int sum;
  int *st_taxa;

  st_taxa = (int *)mCalloc(st->tree->n_otu,sizeof(int));

  For(i,st->tree->n_otu)
    {
      For(j,st->n_gt)
	{
	  For(k,st->treelist->tree[j]->n_otu) 
	    if(!strcmp(st->treelist->tree[j]->noeud[k]->name,st->tree->noeud[i]->name)) break;
	  if(k != st->treelist->tree[j]->n_otu) { st_taxa[i] = 1; break; }
	}
    }

  sum = 0;
  For(i,st->tree->n_otu) if(st_taxa[i]) sum++;
  if(sum != st->tree->n_otu) 
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  Free(st_taxa);
}

/*********************************************************/

int Mg_Get_Species_Found_In_St(superarbre *st, allseq *data)
{
  int i,j;
  
  For(i,data->n_otu)
    {
      For(j,st->tree->n_otu)
	{
	  if(!strcmp(data->c_seq[i]->name,st->tree->noeud[j]->name))
	    {
	      break;
	    }
	}
      if(j == st->tree->n_otu) return 0;
    }
  return 1;

}

/*********************************************************/

void Mg_Map_St_Nodes_In_Gt(arbre *gt, superarbre *st)
{
  int i;

  For(i,2*st->tree->n_otu-2)
    {
      st->map_st_node_in_gt[gt->dp][i][0][0] = NULL;
      st->map_st_node_in_gt[gt->dp][i][1][0] = NULL;
      st->map_st_node_in_gt[gt->dp][i][2][0] = NULL;

      st->map_st_node_in_gt[gt->dp][i][0][1] = NULL;
      st->map_st_node_in_gt[gt->dp][i][1][1] = NULL;
      st->map_st_node_in_gt[gt->dp][i][2][1] = NULL;
    }

  
  /* Root */
  Mg_Map_St_Nodes_In_Gt_One_Edge(st->tree->noeud[0]->v[0],
				 st->tree->noeud[0],
				 st->tree->noeud[0]->b[0],
				 gt,st);

  /* Internal nodes */
  Mg_Map_St_Nodes_In_Gt_Post(st->tree->noeud[0],st->tree->noeud[0]->v[0],gt,st);
  Mg_Map_St_Nodes_In_Gt_Pre (st->tree->noeud[0],st->tree->noeud[0]->v[0],gt,st);
  
  /* Root */
  Mg_Map_St_Nodes_In_Gt_One_Edge(st->tree->noeud[0],
				 st->tree->noeud[0]->v[0],
				 st->tree->noeud[0]->b[0],
				 gt,st);
  
}

/*********************************************************/

void Mg_Map_St_Nodes_In_Gt_Post(node *a_st, node *d_st, arbre *gt, superarbre *st)
{
  int i,dir;

  dir = -1;

  if(d_st->tax) return;
  else
    {
      For(i,3)
	if(d_st->v[i] != a_st)
	  Mg_Map_St_Nodes_In_Gt_Post(d_st,d_st->v[i],gt,st);
      
      For(i,3)
	if(d_st->v[i] != a_st)
	  {
	    Mg_Map_St_Nodes_In_Gt_One_Edge(d_st,d_st->v[i],d_st->b[i],gt,st);
	  }
    }
}

/*********************************************************/

void Mg_Map_St_Nodes_In_Gt_Pre(node *a_st, node *d_st, arbre *gt, superarbre *st)
{
  int i;

  if(d_st->tax) return;
  else
    {
      For(i,3)
	{
	  if(d_st->v[i] != a_st)
	    {	      
	      Mg_Map_St_Nodes_In_Gt_One_Edge(d_st->v[i],d_st,d_st->b[i],gt,st);
	      Mg_Map_St_Nodes_In_Gt_Pre(d_st,d_st->v[i],gt,st);
	    }
	}
    }
}

/*********************************************************/

void Mg_Map_St_Nodes_In_Gt_One_Edge(node *a_st, node *d_st, edge *b_st, arbre *gt, superarbre *st)
{
  if(d_st->tax)
    {
#ifdef DEBUG
      if(b_st->rght != d_st)
	{
	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
#endif
      
      st->map_st_node_in_gt[gt->dp][d_st->num][0][0] = st->match_st_node_in_gt[gt->dp][d_st->num];
    }
  else
    {
      node **list_of_nodes_d, **list_of_nodes_v1, **list_of_nodes_v2;
      int dir1, dir2;
      int i;

      list_of_nodes_d  = NULL;
      list_of_nodes_v1 = NULL;
      list_of_nodes_v2 = NULL;

      dir1 = dir2 = -1;
      For(i,3) 
	{
	  if(d_st->v[i] != a_st) (dir1 < 0)?(dir1 = i):(dir2 = i);
	  else list_of_nodes_d = st->map_st_node_in_gt[gt->dp][d_st->num][i];
	}

      For(i,3) 
	if((d_st->v[dir1]->v[i]) && (d_st->v[dir1]->v[i] == d_st)) 
	  {
	    list_of_nodes_v1 = st->map_st_node_in_gt[gt->dp][d_st->v[dir1]->num][i];
	    break;
	  }

      For(i,3) 
	if((d_st->v[dir2]->v[i]) && (d_st->v[dir2]->v[i] == d_st)) 
	  {
	    list_of_nodes_v2 = st->map_st_node_in_gt[gt->dp][d_st->v[dir2]->num][i];
	    break;
	  }

      /* d_st matches one node in gt */
      if(st->match_st_node_in_gt[gt->dp][d_st->num])
	{
	  list_of_nodes_d[0] = st->match_st_node_in_gt[gt->dp][d_st->num];
	  list_of_nodes_d[1] = NULL;
	}
      else
	{
	  /* list_of_nodes = union of  list_of_nodes_v1  &  list_of_nodes_v2 */
	  
	  if(!list_of_nodes_v1[0])
	    {
	      list_of_nodes_d[0] = list_of_nodes_v2[0];
	      list_of_nodes_d[1] = list_of_nodes_v2[1];
	    }
	  else if(!list_of_nodes_v2[0])
	    {
	      list_of_nodes_d[0] = list_of_nodes_v1[0];
	      list_of_nodes_d[1] = list_of_nodes_v1[1];
	    }
	  else
	    {
	      list_of_nodes_d[0] = list_of_nodes_v1[0];
	      list_of_nodes_d[1] = list_of_nodes_v2[0];

	      if(list_of_nodes_v1[1] || list_of_nodes_v2[1])
		{
		  Print_Node(st->tree->noeud[0],
			     st->tree->noeud[0]->v[0],
			     st->tree);
		  
		  printf("\n\n--------------------------\n\n");
		  Print_Node(gt->noeud[0],
			     gt->noeud[0]->v[0],
			     gt);

		  printf("\n\n--------------------------\n\n");
		  
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Warn_And_Exit("");
		}
	    }
	}
    }
}

/*********************************************************/

void Mg_Map_St_Edges_In_Gt(arbre *gt, superarbre *st)
{
  int i,j;
  edge *st_b;
  node *gt_a, *gt_d;
  
  gt_a = NULL;
  gt_d = NULL;

  For(i,2*st->tree->n_otu-3) st->map_st_edge_in_gt[gt->dp][i] = NULL;

  For(i,2*st->tree->n_otu-3)
    {
      st_b = st->tree->t_edges[i];

      if(!st->map_st_node_in_gt[gt->dp][st_b->left->num][st_b->l_r][0])
	{
	  gt_a = st->map_st_node_in_gt[gt->dp][st_b->rght->num][st_b->r_l][0];
	  gt_d = st->map_st_node_in_gt[gt->dp][st_b->rght->num][st_b->r_l][1];

	  For(j,3)
	    {
	      if((gt_a->v[j]) && (gt_a->v[j] == gt_d))
		{
		  st->map_st_edge_in_gt[gt->dp][st_b->num] = gt_a->b[j];
		  break;
		}
	    }
	}
      else if(!st->map_st_node_in_gt[gt->dp][st_b->rght->num][st_b->r_l][0])
	{
	  gt_a = st->map_st_node_in_gt[gt->dp][st_b->left->num][st_b->l_r][0];
	  gt_d = st->map_st_node_in_gt[gt->dp][st_b->left->num][st_b->l_r][1];

	  For(j,3)
	    {
	      if((gt_a->v[j]) && (gt_a->v[j] == gt_d))
		{
		  st->map_st_edge_in_gt[gt->dp][st_b->num] = gt_a->b[j];
		  break;
		}
	    }
	}
      else
	{	  
	  gt_a = st->map_st_node_in_gt[gt->dp][st_b->left->num][st_b->l_r][0];
	  gt_d = st->map_st_node_in_gt[gt->dp][st_b->rght->num][st_b->r_l][0];

	  #ifdef DEBUG
	  if((!gt_a) || (!gt_d))
	    {
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  #endif

	  For(j,3)
	    {
	      if((gt_a->v[j]) && (gt_a->v[j] == gt_d))
		{
		  st->map_st_edge_in_gt[gt->dp][st_b->num] = gt_a->b[j];
		  break;
		}
	    }
	}
    }
}

/*********************************************************/

void Mg_Map_Gt_Edges_In_St(arbre *gt, superarbre *st)
{
  int i;
  edge *st_b, *gt_b;
  
  For(i,2*st->tree->n_otu-3) st->size_map_gt_edge_in_st[gt->dp][i] = 0;

  st_b = NULL;
  gt_b = NULL;
  For(i,2*st->tree->n_otu-3)
    {
      st_b = st->tree->t_edges[i];
      gt_b = st->map_st_edge_in_gt[gt->dp][st_b->num];

      if(gt_b)
	{
	  st->map_gt_edge_in_st[gt->dp][gt_b->num][st->size_map_gt_edge_in_st[gt->dp][gt_b->num]] = st_b;
	  st->size_map_gt_edge_in_st[gt->dp][gt_b->num]++;
	}
    }
}

/*********************************************************/

int Mg_Pars(superarbre *st)
{
  int i;
  
  st->tree->c_pars = 0;
  For(i,st->n_gt) 
    {
      st->treelist->tree[i]->both_sides = 1;	  
      Pars(st->treelist->tree[i]);
/*       printf("\n. Tree %3d pars = %d",i+1,st->treelist->tree[i]->c_pars); */
      st->tree->c_pars += st->treelist->tree[i]->c_pars;
    }

  return st->tree->c_pars;
}

/*********************************************************/

int Mg_Spr(phydbl init_lnL, superarbre *st)
{
  int gt;
  int i;
  edge *target, *pruned;
  int move,best_move;
  node *gt_a, *gt_d;

  //st->tree->root       = st->tree->noeud[0];
  st->tree->both_sides = 1;
  pruned               = NULL;
  target               = NULL;
  move                 = -2;
  gt_a                 = NULL;
  gt_d                 = NULL;


  For(i,2*st->tree->n_otu-3)
    {
      pruned            = st->tree->t_edges[i];
      st->tree->n_moves = 0;

      Reset_Spr_List(st->tree);
      For(gt,st->n_gt) Reset_Spr_List(st->treelist->tree[gt]);
      
      if(!pruned->rght->tax)
	{
	  For(gt,st->n_gt)
	    {
	      /* Check constraints at prune site on gt tree */
 	      gt_a = st->map_st_node_in_gt[gt][pruned->rght->num][pruned->r_l][0];
 	      gt_d = st->map_st_node_in_gt[gt][pruned->left->num][pruned->l_r][0];

	      if((gt_a) && (gt_d) && (!st->map_st_edge_in_gt[gt][pruned->num]->rght->tax))
		{
		  Test_All_Spr_Targets(st->map_st_edge_in_gt[gt][pruned->num],
				       st->map_st_edge_in_gt[gt][pruned->num]->rght,
				       st->treelist->tree[gt]);
		}
	    }
	}

      if(!pruned->left->tax)
	{
	  For(gt,st->n_gt)
	    {
	      /* Check constraints at prune site on gt tree */	      
 	      gt_a = st->map_st_node_in_gt[gt][pruned->rght->num][pruned->r_l][0];
 	      gt_d = st->map_st_node_in_gt[gt][pruned->left->num][pruned->l_r][0];

	      if((gt_a) && (gt_d) && (!st->map_st_edge_in_gt[gt][pruned->num]->left->tax))
		{
		  Test_All_Spr_Targets(st->map_st_edge_in_gt[gt][pruned->num],
				       st->map_st_edge_in_gt[gt][pruned->num]->left,
				       st->treelist->tree[gt]);
		}
	    }
	}

      
      if(!pruned->left->tax)
	{
	  Mg_Test_All_Spr_Targets(st->tree->t_edges[i],
				  st->tree->t_edges[i]->left,
				  st);      
	}
      
      if(!pruned->rght->tax)
	{
	  Mg_Test_All_Spr_Targets(st->tree->t_edges[i],
				  st->tree->t_edges[i]->rght,
				  st);      
	}


      if(st->tree->n_moves)
	{
	  best_move = Mg_Test_List_Of_Regraft_Pos(st->tree->spr_list,
						  (int)ceil(0.1*(st->tree->n_moves)),
						  st);	  
	  
	  if(st->tree->spr_list[best_move]->lnL > init_lnL)
	    {
	      Mg_Try_One_Spr_Move(st->tree->spr_list[best_move],st);
	    }
	  else
	    {
	      st->tree->both_sides = 1;
	      st->tree->c_lnL      = Mg_Lk(st);
	      st->tree->c_pars     = Mg_Pars(st);	      
	    }
	}
    }
  return 1;
}

/*********************************************************/

void Mg_Speed_Spr(superarbre *st)
{
  int step;
  int gt;
  phydbl old_lnL;

  Make_Spr_List(st->tree);
  For(gt,st->n_gt) Make_Spr_List(st->treelist->tree[gt]);

  st->tree->both_sides = 1; 
  For(gt,st->n_gt) 
    {
      st->treelist->tree[gt]->both_sides = 1;
      Record_Br_Len(st->treelist->tree[gt]);
    }
  
  st->tree->c_pars = Mg_Pars(st);
  st->tree->c_lnL  = Mg_Lk(st);
  

  st->tree->best_lnL = st->tree->c_lnL;
  old_lnL            = st->tree->c_lnL;
  step               = 0;
  do
    {
      ++step;

      printf("\n. Starting a SPR cycle... \n");

      old_lnL = st->tree->c_lnL;

      st->tree->n_improvements         = 0;
      st->tree->perform_spr_right_away = 1;
      Mg_Spr(UNLIKELY,st);

 
      time(&(st->tree->t_current));      
      printf("\n. (%5d sec) [00] [%10.2f] [%5d]\n",
	     (int)(st->tree->t_current-st->tree->t_beg),
	     Mg_Lk(st),Mg_Pars(st));

      /* Optimise parameters of the Markov model */
      For(gt,st->n_gt) Optimiz_All_Free_Param(st->treelist->tree[gt],
					      st->treelist->tree[gt]->mod->s_opt->print);

      time(&(st->tree->t_current));      
      printf("\n. (%5d sec) [ 0] [%10.2f] [%5d]\n",
	     (int)(st->tree->t_current-st->tree->t_beg),
	     Mg_Lk(st),Mg_Pars(st));

      /* Optimise branch lengths */
      For(gt,st->n_gt)
	{
	  Optimize_Br_Len_Serie(st->treelist->tree[gt]->noeud[0],
				st->treelist->tree[gt]->noeud[0]->v[0],
				st->treelist->tree[gt]->noeud[0]->b[0],
				st->treelist->tree[gt],
				st->treelist->tree[gt]->data);
	}


      /* Update partial likelihoods & parsimony */
      st->tree->both_sides = 1; 
      st->tree->c_pars = Mg_Pars(st);
      st->tree->c_lnL  = Mg_Lk(st);
      
      
      time(&(st->tree->t_current));      
      printf("\n. (%5d sec) [**] [%10.2f] [%5d]\n",
	     (int)(st->tree->t_current-st->tree->t_beg),
	     st->tree->c_lnL,st->tree->c_pars);

      /* Record the current best log-likleihood  */
      st->tree->best_lnL = st->tree->c_lnL;

      if(st->tree->c_lnL < old_lnL)
	{
	  printf("\n. old_lnL = %f c_lnL = %f\n",old_lnL,st->tree->c_lnL); 
	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      /* Record the current best branch lengths  */
      For(gt,st->n_gt) Record_Br_Len(st->treelist->tree[gt]);

      /* Exit if no improvements after complete optimization */
      if((!st->tree->n_improvements) || 
	 (fabs(old_lnL-st->tree->c_lnL) < st->tree->mod->s_opt->min_diff_lk_global)) break;
            
    }while(1);
  
}

/*********************************************************/


void Mg_Test_All_Spr_Targets(edge *pruned, node *n_link, superarbre *st)
{
  int i,j;

  For(i,3)
    {
      if(n_link->b[i] != pruned)
	{
	  For(j,3)
	    {
	      if((n_link->v[i]->v[j]) && (n_link->v[i]->v[j] != n_link))
		{
		  Mg_Test_One_Spr_Target_Recur(n_link->v[i],n_link->v[i]->v[j],n_link->v[i]->b[j],pruned,n_link,st);
		}
	    }
	}
    }
}

/*********************************************************/

void Mg_Test_One_Spr_Target_Recur(node *a, node *d, edge *target, edge *pruned, node *n_link, superarbre *st)
{

  Mg_Test_One_Spr_Target(pruned,target,n_link,st);

  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	if(d->v[i] != a)
	  {
	    Mg_Test_One_Spr_Target_Recur(d,d->v[i],d->b[i],pruned,n_link,st);
	  }
    }
}

/*********************************************************/

void Mg_Test_One_Spr_Target(edge *st_p, edge *st_t, node *n_link, superarbre *st) 
{
  int gt, move;

  st->tree->n_moves++;
  st->tree->spr_list[st->tree->size_spr_list]->b_target      = st_t;
  st->tree->spr_list[st->tree->size_spr_list]->n_link        = n_link;
  st->tree->spr_list[st->tree->size_spr_list]->n_opp_to_link = (n_link == st_p->left)?(st_p->rght):(st_p->left);
  st->tree->spr_list[st->tree->size_spr_list]->b_opp_to_link = st_p;
  st->tree->spr_list[st->tree->size_spr_list]->pars          = 0;

  For(gt,st->n_gt)
    {
      move = Map_Spr_Move(st_p,st_t,n_link,st->treelist->tree[gt],st);
      
      if(move > -1)
	st->tree->spr_list[st->tree->size_spr_list]->pars += st->treelist->tree[gt]->spr_list[move]->pars;
      else if(move == -1 || move == -2)
	st->tree->spr_list[st->tree->size_spr_list]->pars += st->treelist->tree[gt]->c_pars;
    }

  Include_One_Spr_To_List_Of_Spr(st->tree->spr_list[st->tree->size_spr_list],st->tree);
}

/*********************************************************/

int Map_Spr_Move(edge *st_pruned, edge *st_target, node *st_link, arbre *gt, superarbre *st)
{
  int i;
  edge *gt_pruned, *gt_target;
  node *gt_link, *gt_a, *gt_d;

  gt_pruned = NULL;
  gt_target = NULL;
  gt_link   = NULL;
  gt_a      = NULL;
  gt_d      = NULL;

  /* Check contraints at prune and regraft sites on the gt tree */

  /* st_pruned is not on a path that matches a branch in gt */
  gt_a = st->map_st_node_in_gt[gt->dp][st_pruned->left->num][st_pruned->l_r][0];
  gt_d = st->map_st_node_in_gt[gt->dp][st_pruned->rght->num][st_pruned->r_l][0];

  if((!gt_a) || (!gt_d)) return -1;
  else
    {
      /* which gt nodes matches st_link ? */
      gt_link = (st_pruned->left == st_link)?(gt_a):(gt_d);
      
      if(gt_link->tax) return -1;
      else
	{
	  gt_pruned = st->map_st_edge_in_gt[gt->dp][st_pruned->num];
	  gt_target = st->map_st_edge_in_gt[gt->dp][st_target->num];
	   
	  if((gt_pruned->left == gt_target->left) ||
	     (gt_pruned->left == gt_target->rght) ||
	     (gt_pruned->rght == gt_target->left) ||
	     (gt_pruned->rght == gt_target->rght)) return -1;
	  else
	    {
	      For(i,gt->size_spr_list)
		{
		  if((gt_pruned == gt->spr_list[i]->b_opp_to_link) && (gt_target == gt->spr_list[i]->b_target))
		    return i;
		}
	    }
	}
    }
  return -2;
}

/*********************************************************/

int Mg_Test_List_Of_Regraft_Pos(spr **st_spr_list, int list_size, superarbre *st)
{

  int i,j,best_move;
  spr *move;
  edge *init_target, *b_residual;
  phydbl best_lnL, init_lnL;
  int dir_v0, dir_v1, dir_v2;
  int gt;
  int move_num;
  

  best_lnL = UNLIKELY;
  init_target = b_residual = NULL;
  best_move = -1;

#ifdef DEBUG
  if(!list_size)
    {
      printf("\n\n. List size is 0 !");
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit(""); 
    }
#endif
  
  init_lnL = UNLIKELY;

  For(i,list_size)
    {
      st->tree->spr_list[i]->lnL = .0;

      For(gt,st->n_gt) 
	{
	  move_num = Map_Spr_Move(st->tree->spr_list[i]->b_opp_to_link,
				  st->tree->spr_list[i]->b_target,
				  st->tree->spr_list[i]->n_link,
				  st->treelist->tree[gt],st);

	  if(move_num > -1)
	    {
	      move = st->treelist->tree[gt]->spr_list[move_num];

	      if(move->b_target)
		{
		  init_lnL = st->treelist->tree[gt]->c_lnL;

		  /* Record edge lengths */
		  Record_Br_Len(st->treelist->tree[gt]);
		  
		  /* Prune subtree */
		  Prune_Subtree(move->n_link,
				move->n_opp_to_link,			    
				&init_target,
				&b_residual,
				st->treelist->tree[gt]);

		  /* Rough optimisation of the branch length */
		  Fast_Br_Len(init_target,st->treelist->tree[gt]);
		  
		  /* Update the change proba matrix at prune position */
		  Update_PMat_At_Given_Edge(init_target,st->treelist->tree[gt]);
	      
		  /* Update partial likelihood along the path from the prune to
		     the regraft position */
		  Update_P_Lk_Along_A_Path(move->path,move->depth_path,st->treelist->tree[gt]);

		  /* Regraft subtree */
		  Graft_Subtree(move->b_target,move->n_link,b_residual,st->treelist->tree[gt]);
	      
		  /* Estimate the three edge lengths at the regraft site */
		  Triple_Dist(move->n_link,st->treelist->tree[gt]);
	      
		  /* Update the transition proba matrices along edges defining 
		     the regraft site */
		  For(j,3)
		    if(move->n_link->v[j] != move->n_opp_to_link)
		      Update_PMat_At_Given_Edge(move->n_link->b[j],st->treelist->tree[gt]);
	      
		  /* Compute the likelihood */
		  Update_P_Lk(st->treelist->tree[gt],
			      move->b_opp_to_link,
			      move->n_link);
		  
		  move->lnL = Lk_At_Given_Edge(move->b_opp_to_link,st->treelist->tree[gt]);

		  
		  st->tree->spr_list[i]->lnL += move->lnL;

		  /* Record branch lengths */
		  dir_v1 = dir_v2 = dir_v0 = -1;
		  For(j,3)
		    {
		      if(move->n_link->v[j] == move->n_opp_to_link) dir_v0 = j;
		      else if(dir_v1 < 0)                           dir_v1 = j;
		      else                                          dir_v2 = j;
		    }
		  
		  move->l0 = move->n_link->b[dir_v0]->l;
		  
		  if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num)
		    {
		      move->l1 = move->n_link->b[dir_v2]->l;
		      move->l2 = move->n_link->b[dir_v1]->l;
		    }
		  else
		    {
		      move->l1 = move->n_link->b[dir_v1]->l;
		      move->l2 = move->n_link->b[dir_v2]->l;
		    }
		  	  
		  /* Regraft the subtree at its original position */
		  Prune_Subtree(move->n_link,
				move->n_opp_to_link,
				&move->b_target,
				&b_residual,
				st->treelist->tree[gt]);
		  
		  Graft_Subtree(init_target,
				move->n_link,
				b_residual,
				st->treelist->tree[gt]);
		  
		  /* Restore branch lengths */
		  Restore_Br_Len(st->treelist->tree[gt]);
	      
		  /* Update relevant change proba matrices */
		  Update_PMat_At_Given_Edge(move->b_target,st->treelist->tree[gt]);
		  For(j,3) Update_PMat_At_Given_Edge(move->n_link->b[j],st->treelist->tree[gt]);
		  
		  /* Update relevant partial likelihoods */
		  For(j,3) Update_P_Lk(st->treelist->tree[gt],move->n_link->b[j],move->n_link);
		  
		  st->treelist->tree[gt]->c_lnL = init_lnL;
		}
	    }
	  else
	    {
	      st->tree->spr_list[i]->lnL += st->treelist->tree[gt]->c_lnL;
	    }
	}

      if(st->tree->spr_list[i]->lnL > best_lnL)
	{
	  best_lnL  = st->tree->spr_list[i]->lnL;
	  best_move = i;
	}
    }

  return best_move;  

}


/*********************************************************/

int Mg_Try_One_Spr_Move(spr *st_move, superarbre *st)
{
  int j;
  spr **gt_move;
  edge **init_target, **b_residual;
  int dir_v0, dir_v1, dir_v2;
  int gt;
  int gt_move_num;
  int n_moves;

  
  init_target = (edge **)mCalloc(st->n_gt,sizeof(edge *));
  b_residual  = (edge **)mCalloc(st->n_gt,sizeof(edge *));
  gt_move     = (spr **)mCalloc(st->n_gt,sizeof(spr *));


  n_moves = 0;
  For(gt,st->n_gt) 
    {
      gt_move_num = Map_Spr_Move(st_move->b_opp_to_link,
				 st_move->b_target,
				 st_move->n_link,
				 st->treelist->tree[gt],st);
      
      if(gt_move_num > -1)
	{
	  n_moves++;

	  gt_move[gt] = st->treelist->tree[gt]->spr_list[gt_move_num];
	  
	  if(gt_move[gt]->b_target)
	    {
	      /* Record edge lengths */
	      Record_Br_Len(st->treelist->tree[gt]);

	      /* Prune subtree */
	      Prune_Subtree(gt_move[gt]->n_link,
			    gt_move[gt]->n_opp_to_link,
			    &(init_target[gt]),
			    &(b_residual[gt]),
			    st->treelist->tree[gt]);
	      
	      /* Rough optimisation of the branch length */
	      Fast_Br_Len(init_target[gt],st->treelist->tree[gt]);
	      
	      /* Update the change proba matrix at prune position */
	      Update_PMat_At_Given_Edge(init_target[gt],st->treelist->tree[gt]); /* TO DO : NECESSARY ?? */
	      
	      /* Update partial likelihood along the path from the prune to
		 the regraft position */
	      Update_P_Lk_Along_A_Path(gt_move[gt]->path,gt_move[gt]->depth_path,st->treelist->tree[gt]); /* TO DO : NECESSARY ?? */
	      
	      /* Regraft subtree */
	      Graft_Subtree(gt_move[gt]->b_target,gt_move[gt]->n_link,b_residual[gt],st->treelist->tree[gt]);
	      
	      dir_v1 = dir_v2 = dir_v0 = -1;
	      For(j,3)
		{
		  if(gt_move[gt]->n_link->v[j] == gt_move[gt]->n_opp_to_link) dir_v0 = j;
		  else if(dir_v1 < 0)                                         dir_v1 = j;
		  else                                                        dir_v2 = j;
		}
	      
	      gt_move[gt]->n_link->b[dir_v0]->l = gt_move[gt]->l0;
		  
	      if(gt_move[gt]->n_link->v[dir_v1]->num > gt_move[gt]->n_link->v[dir_v2]->num)
		{
		  gt_move[gt]->n_link->b[dir_v2]->l = gt_move[gt]->l1;
		  gt_move[gt]->n_link->b[dir_v1]->l = gt_move[gt]->l2;
		}
	      else
		{
		  gt_move[gt]->n_link->b[dir_v1]->l = gt_move[gt]->l1;
		  gt_move[gt]->n_link->b[dir_v2]->l = gt_move[gt]->l2;
		}
	    }
	}
    }
  

  if(n_moves)
    {
      if(st_move->lnL > st->tree->best_lnL)
	{
	  edge *st_target, *st_residual;

	  /* Apply the move on the super-tree */
	  Prune_Subtree(st_move->n_link,
			st_move->n_opp_to_link,
			&st_target,
			&st_residual,
			st->tree);
	  
	  Graft_Subtree(st_move->b_target,
			st_move->n_link,
			st_residual,
			st->tree);
	  
	  	  
	  /* Map gt and st nodes and edges */
	  Mg_Do_Mapping(st);

	  time(&(st->tree->t_current));
	  st->tree->both_sides = 1;
	  
	  st->tree->both_sides = 1;
	  st->tree->c_lnL      = Mg_Lk(st);
	  st->tree->c_pars     = Mg_Pars(st);
	  
	  
	  if(fabs(st->tree->c_lnL - st_move->lnL) > st->tree->mod->s_opt->min_diff_lk_local)
	    {
	      printf("\n. st->tree->c_lnL = %f st_move->lnL = %f\n",
		     st->tree->c_lnL,st_move->lnL);

	      For(gt,st->n_gt)
		{
		  printf("\n. truth -> %f ; move -> %f",
			 Return_Lk(st->treelist->tree[gt]),
			 gt_move[gt] ? gt_move[gt]->lnL : -1.);
		}
	    }
	  
	  printf("\n. (%5d sec) [+ ] [%10.2f] [%5d] -- ",
		 (int)(st->tree->t_current - st->tree->t_beg),
		 st->tree->c_lnL,st->tree->c_pars);	  
	  
	  For(gt,st->n_gt)
	    printf("[%10.2f] ",st->treelist->tree[gt]->c_lnL);
	  
	  
	  st->tree->n_improvements++;
	  st->tree->best_lnL = st->tree->c_lnL;
	  For(gt,st->n_gt) Record_Br_Len(st->treelist->tree[gt]);
	  
	  Free(init_target);
	  Free(b_residual);
	  Free(gt_move);
	  
	  return 1;
	}
      else
	{
	  For(gt,st->n_gt) 
	    {
	      if(gt_move[gt])
		{
		  Return_Lk(st->treelist->tree[gt]);
		  Fast_Br_Len_Recur(st->treelist->tree[gt]->noeud[0],
				    st->treelist->tree[gt]->noeud[0]->v[0],
				    st->treelist->tree[gt]->noeud[0]->b[0],
				    st->treelist->tree[gt]);
		}
	    }
	  
	  time(&(st->tree->t_current));
	  st->tree->both_sides = 1;
	  st->tree->c_lnL      = Mg_Lk(st);
	  
	  if(st->tree->c_lnL > st->tree->best_lnL)
	    {
	      edge *st_target, *st_residual;
	      
	      /* Apply the move on the super-tree */
	      Prune_Subtree(st_move->n_link,
			    st_move->n_opp_to_link,			    
			    &st_target,
			    &st_residual,
			    st->tree);
	      
	      Graft_Subtree(st_move->b_target,
			    st_move->n_link,
			    st_residual,
			    st->tree);
	      
	      
	      /* Map gt and st nodes and edges */
	      Mg_Do_Mapping(st);


	      st->tree->c_pars = Mg_Pars(st);
	      printf("\n. (%5d sec) [++] [%10.2f] [%5d] -- ",
		     (int)(st->tree->t_current-st->tree->t_beg),
		     st->tree->c_lnL,
		     st->tree->c_pars);
	      For(gt,st->n_gt)
		printf("[%10.2f] ",st->treelist->tree[gt]->c_lnL);
	      
	      st->tree->n_improvements++;
	      st->tree->best_lnL = st->tree->c_lnL;
	      For(gt,st->n_gt) Record_Br_Len(st->treelist->tree[gt]);

	      Free(init_target);
	      Free(b_residual);
	      Free(gt_move);

	      return 1;
	    }
	}
    }
  
  For(gt,st->n_gt) 
    {
      if(gt_move[gt])
	{	  
	  /* Regraft the subtree at its original position */
	  Prune_Subtree(gt_move[gt]->n_link,
			gt_move[gt]->n_opp_to_link,
			&(gt_move[gt]->b_target),
			&(b_residual[gt]),
			st->treelist->tree[gt]);

	  Graft_Subtree(init_target[gt],
			gt_move[gt]->n_link,
			b_residual[gt],
			st->treelist->tree[gt]);	  

	  /* Restore branch lengths */
	  Restore_Br_Len(st->treelist->tree[gt]);
	}
    }
  
  st->tree->both_sides = 1;
  st->tree->c_lnL      = Mg_Lk(st);
  st->tree->c_pars     = Mg_Pars(st);

  time(&(st->tree->t_current));
  
  printf("\n. (%5d sec) [--] [%10.2f] [%5d] -- ",
	 (int)(st->tree->t_current - st->tree->t_beg),
	 st->tree->c_lnL,st->tree->c_pars);	  
  
  For(gt,st->n_gt) printf("[%10.2f] ",st->treelist->tree[gt]->c_lnL);

  Free(init_target);
  Free(b_residual);
  Free(gt_move);

  return 0;

}

/*********************************************************/

void Mg_NNI(edge *st_b, superarbre *st)
{  
  node *v1, *v2, *v3, *v4;
  phydbl lk_init;
  phydbl lk0_init, lk1_init, lk2_init;
  phydbl lk0_opt, lk1_opt, lk2_opt;
  int i,j;
  double *init_bl;
  edge **map_edge_bef_swap, **map_edge_aft_swap;


  init_bl = (double *)mCalloc(st->n_bl_partition,sizeof(double));
  map_edge_bef_swap = (edge **)mCalloc(st->n_gt,sizeof(edge *));
  map_edge_aft_swap = (edge **)mCalloc(st->n_gt,sizeof(edge *));


  v1 = st_b->left->v[st_b->l_v1];
  v2 = st_b->left->v[st_b->l_v2];
  v3 = st_b->rght->v[st_b->r_v1];
  v4 = st_b->rght->v[st_b->r_v2];

  lk0_init = lk1_init = lk2_init = UNLIKELY;
  lk0_opt  = lk1_opt  = lk2_opt  = UNLIKELY;

  if(v1->num < v2->num)
    {
      Check_Dirs(st->tree);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(v3->num < v4->num)
    {
      Check_Dirs(st->tree);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  lk_init = st->tree->c_lnL;
  
/*   printf("oooooooo\n"); */
/*   Print_Node(st->tree->noeud[0], */
/* 	     st->tree->noeud[0]->v[0], */
/* 	     st->tree); */
/*   printf(">>>>>>>\n"); */
/*   For(i,st->n_gt) */
/*     { */
/*       Print_Node(st->treelist->tree[i]->noeud[0], */
/* 		 st->treelist->tree[i]->noeud[0]->v[0], */
/* 		 st->treelist->tree[i]); */
/*       printf("<<<<<<<\n"); */
/*     } */

  
  Mg_Record_Br_Len(st);

  For(i,st->n_gt) map_edge_bef_swap[i] = NULL;
  For(i,st->n_gt) if(st->map_st_edge_in_gt[i][st_b->num]) map_edge_bef_swap[i] = st->map_st_edge_in_gt[i][st_b->num];

  /* First alternative topological configuration */
  /* Swap */
  Mg_Swap(v2,st_b->left,st_b->rght,v3,st);
  Swap(v2,st_b->left,st_b->rght,v3,st->tree);
  Mg_Do_Mapping(st);
  Mg_Set_Bl(st->bl,st);
  For(i,st->n_gt) map_edge_aft_swap[i] = NULL;
  For(i,st->n_gt) if(st->map_st_edge_in_gt[i][st_b->num]) map_edge_aft_swap[i] = st->map_st_edge_in_gt[i][st_b->num];
  For(i,st->n_gt) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_gt) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  For(i,st->n_gt) if(map_edge_bef_swap[i] && map_edge_aft_swap[i]) 
    {
      For(j,3) if((map_edge_aft_swap[i]->left->v[j]) && 
		  (map_edge_aft_swap[i]->left->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->left);
      For(j,3) if((map_edge_aft_swap[i]->rght->v[j]) && 
		  (map_edge_aft_swap[i]->rght->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->rght);
    }
  lk1_init = Mg_Update_Lk_At_Given_Edge(st_b,st);
  lk1_opt  = Mg_Br_Len_Brent(st_b,st);
  For(i,st->n_gt) st->bl1[st->bl_partition[i]][st_b->num] = st->bl[st->bl_partition[i]][st_b->num];
  /* Unswap */
  Mg_Swap(v3,st_b->left,st_b->rght,v2,st);
  Swap(v3,st_b->left,st_b->rght,v2,st->tree);
  Mg_Do_Mapping(st);
  Mg_Restore_Br_Len(st);
  Mg_Set_Bl(st->bl,st);
  For(i,st->n_gt) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_gt) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  For(i,st->n_gt) if(map_edge_bef_swap[i] && map_edge_aft_swap[i]) 
    {
      For(j,3) if((map_edge_aft_swap[i]->left->v[j]) && 
		  (map_edge_aft_swap[i]->left->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->left);
      For(j,3) if((map_edge_aft_swap[i]->rght->v[j]) && 
		  (map_edge_aft_swap[i]->rght->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->rght);
    }




  /* Second alternative topological configuration */
  /* Swap */
  Mg_Swap(v2,st_b->left,st_b->rght,v4,st);
  Swap(v2,st_b->left,st_b->rght,v4,st->tree);
  Mg_Do_Mapping(st);
  Mg_Set_Bl(st->bl,st);
  For(i,st->n_gt) map_edge_aft_swap[i] = NULL;
  For(i,st->n_gt) if(st->map_st_edge_in_gt[i][st_b->num]) map_edge_aft_swap[i] = st->map_st_edge_in_gt[i][st_b->num];
  For(i,st->n_gt) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_gt) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  For(i,st->n_gt) if(map_edge_bef_swap[i] && map_edge_aft_swap[i]) 
    {
      For(j,3) if((map_edge_aft_swap[i]->left->v[j]) && 
		  (map_edge_aft_swap[i]->left->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->left);
      For(j,3) if((map_edge_aft_swap[i]->rght->v[j]) && 
		  (map_edge_aft_swap[i]->rght->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->rght);
    }

  lk2_init = Mg_Update_Lk_At_Given_Edge(st_b,st);
  lk2_opt  = Mg_Br_Len_Brent(st_b,st);
  For(i,st->n_gt) st->bl2[st->bl_partition[i]][st_b->num] = st->bl[st->bl_partition[i]][st_b->num];
  /*   printf("\n. lk2_init = %f lk2_opt = %f",lk2_init,lk2_opt); */
  /* Unswap */
  Mg_Swap(v4,st_b->left,st_b->rght,v2,st);
  Swap(v4,st_b->left,st_b->rght,v2,st->tree);
  Mg_Do_Mapping(st);
  Mg_Restore_Br_Len(st);
  Mg_Set_Bl(st->bl,st);
  For(i,st->n_gt) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_gt) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  For(i,st->n_gt) if(map_edge_bef_swap[i] && map_edge_aft_swap[i]) 
    {
      For(j,3) if((map_edge_aft_swap[i]->left->v[j]) && 
		  (map_edge_aft_swap[i]->left->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->left);
      For(j,3) if((map_edge_aft_swap[i]->rght->v[j]) && 
		  (map_edge_aft_swap[i]->rght->b[j] == map_edge_bef_swap[i]) &&
		  (map_edge_aft_swap[i] != map_edge_bef_swap[i])) 
	Update_P_Lk(st->treelist->tree[i],
		    map_edge_aft_swap[i],
		    map_edge_aft_swap[i]->rght);
    }


  /* Back to the initial topological configuration 
   * and branch lengths.
   */
  Mg_Do_Mapping(st);
  Mg_Set_Bl(st->bl,st);
  Mg_Restore_Br_Len(st);
  For(i,st->n_gt) map_edge_aft_swap[i] = NULL;
  For(i,st->n_gt) if(st->map_st_edge_in_gt[i][st_b->num]) map_edge_aft_swap[i] = st->map_st_edge_in_gt[i][st_b->num];
  For(i,st->n_gt) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_gt) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  lk0_init = Mg_Update_Lk_At_Given_Edge(st_b,st);
  lk0_opt  = Mg_Br_Len_Brent(st_b,st);
  For(i,st->n_gt) st->bl0[st->bl_partition[i]][st_b->num] = st->bl[st->bl_partition[i]][st_b->num];

  Mg_Restore_Br_Len(st);
  Mg_Set_Bl(st->bl,st);
  For(i,st->n_gt) if(map_edge_bef_swap[i]) Update_PMat_At_Given_Edge(map_edge_bef_swap[i],st->treelist->tree[i]);
  For(i,st->n_gt) if(map_edge_aft_swap[i]) Update_PMat_At_Given_Edge(map_edge_aft_swap[i],st->treelist->tree[i]);
  Mg_Update_Lk_At_Given_Edge(st_b,st);



/*   For(i,2*st->tree->n_otu-3) */
/*     printf("\n. 3 Edge %3d --> lnL=%f",i,Mg_Lk_At_Given_Edge(st->tree->t_edges[i],st)); */



  st_b->nni->lk0 = lk0_opt;
  st_b->nni->lk1 = lk1_opt;
  st_b->nni->lk2 = lk2_opt;

  st_b->nni->score = lk0_opt - MAX(lk1_opt,lk2_opt);

  if((st_b->nni->score <  st->tree->mod->s_opt->min_diff_lk_local) &&
     (st_b->nni->score > -st->tree->mod->s_opt->min_diff_lk_local))
    {
      st_b->nni->score = .0;
      st_b->nni->lk1   = st_b->nni->lk0;
      st_b->nni->lk2   = st_b->nni->lk0;
     }

  Mg_Restore_Br_Len(st);
  Mg_Update_Lk_At_Given_Edge(st_b,st); /* to replace by Mg_Update_PMat_At_Given_Edge(st_b,st); */
/*   printf("\n. lk_end = %f",st->tree->c_lnL); */
/*   For(i,2*st->tree->n_otu-3) printf("\n. %f",Mg_Lk_At_Given_Edge(st->tree->t_edges[i],st)); */
/*   printf("\n. lk_end = %f",Mg_Lk(st)); */
/*   printf("\n"); */

/*   printf("\n. Edge %3d, score = %20f",st_b->num,st_b->nni->score); */

  if(st_b->num == 90)
    printf("\n. v1=%d v2=%d v3=%d v4=%d left-%d right-%d",
	   v1->num,
	   v2->num,
	   v3->num,
	   v4->num,
	   st_b->left->num,
	   st_b->rght->num);



  if(lk0_opt > MAX(lk1_opt,lk2_opt))
    {
      st_b->nni->best_conf = 0;
      st_b->nni->swap_node_v1 = NULL;
      st_b->nni->swap_node_v2 = NULL;
      st_b->nni->swap_node_v3 = NULL;
      st_b->nni->swap_node_v4 = NULL;
    }
  else if(lk1_opt > MAX(lk0_opt,lk2_opt))
    {
      st_b->nni->best_conf    = 1;
      st_b->nni->swap_node_v1 = v2;
      st_b->nni->swap_node_v2 = st_b->left;
      st_b->nni->swap_node_v3 = st_b->rght;
      st_b->nni->swap_node_v4 = v3;
    }
  else if(lk2_opt > MAX(lk0_opt,lk1_opt))
    {
      st_b->nni->best_conf    = 2;
      st_b->nni->swap_node_v1 = v2;
      st_b->nni->swap_node_v2 = st_b->left;
      st_b->nni->swap_node_v3 = st_b->rght;
      st_b->nni->swap_node_v4 = v4;
    }
  else
    {
      st_b->nni->score        = +1.0;
      st_b->nni->best_conf    = 0;
      st_b->nni->swap_node_v1 = NULL;
      st_b->nni->swap_node_v2 = NULL;
      st_b->nni->swap_node_v3 = NULL;
      st_b->nni->swap_node_v4 = NULL;
    }
  
  Free(init_bl);
  Free(map_edge_aft_swap);
  Free(map_edge_bef_swap);
}

/*********************************************************/

void Mg_Swap(node *st_a, node *st_b, node *st_c, node *st_d, superarbre *st)
{
  int i,j;
  node *gt_a, *gt_b, *gt_c, *gt_d;
  int ab, ba, cd, dc, bc;

  ab = ba = cd = dc = bc = -1;

  For(i,3) if(st_a->v[i] == st_b) { ab = i; break; }
  For(i,3) if(st_b->v[i] == st_a) { ba = i; break; }
  For(i,3) if(st_c->v[i] == st_d) { cd = i; break; }
  For(i,3) if(st_d->v[i] == st_c) { dc = i; break; }
  For(i,3) if(st_b->v[i] == st_c) { bc = i; break; }

  if(ab < 0 || ba < 0 || cd < 0 || dc < 0)
    {
      printf("\n. Nodes %d %d %d %d\n",st_a->num,st_b->num,st_c->num,st_d->num);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  gt_a = gt_b = gt_c = gt_d = NULL;
  
  For(i,st->n_gt)
    {
      gt_b = st->match_st_node_in_gt[i][st_b->num];
      gt_c = st->match_st_node_in_gt[i][st_c->num];
      
      if(gt_b && gt_c) /* The st edge with st_b and st_c at its extremities
			* matches an edge in gt 
		        */
	{
#ifdef DEBUG
	  For(j,3) if((gt_b->v[j]) && (gt_b->v[j] == gt_c)) break;
	  if(j == 3)
	    {
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
#endif
	  gt_a = st->map_st_node_in_gt[i][st_a->num][ab][0];
	  gt_d = st->map_st_node_in_gt[i][st_d->num][dc][0];
	  Swap(gt_a,gt_b,gt_c,gt_d,st->treelist->tree[i]);
	}
    }
}

/*********************************************************/

void Mg_Set_Bl(double **bl, superarbre *st)
{
  int i,j;
  edge *gt_b;

  gt_b = NULL;
						
  /* Set all the actual branch lengths to 0.0 
   */
  For(i,st->n_gt)
    {
      For(j,2*st->treelist->tree[i]->n_otu-3)
	{
	  gt_b = st->treelist->tree[i]->t_edges[j];
	  gt_b->l = .0;
	}
    }

  /* Update every branch length 
   */  
  For(i,2*st->tree->n_otu-3)
    {
      For(j,st->n_gt)
	{
	  gt_b = st->map_st_edge_in_gt[j][i];	
	  
	  /* Need to make sure that st->tree->t_edges[i] is on an existing path in gt */
	  if((st->map_st_node_in_gt[j][st->tree->t_edges[i]->left->num][st->tree->t_edges[i]->l_r][0]) &&
	     (st->map_st_node_in_gt[j][st->tree->t_edges[i]->rght->num][st->tree->t_edges[i]->r_l][0]))
	    {
	      gt_b->l += bl[st->bl_partition[j]][i];
	    }
	}
    }
}

/*********************************************************/

void Mg_Record_Br_Len(superarbre *st)
{
  int i,j;
  For(i,st->n_gt) For(j,2*st->tree->n_otu-3) st->bl_cpy[i][j] = st->bl[i][j];
}

/*********************************************************/

void Mg_Restore_Br_Len(superarbre *st)
{
  int i,j;
  For(i,st->n_gt) For(j,2*st->tree->n_otu-3) st->bl[i][j] = st->bl_cpy[i][j];
}

/*********************************************************/

phydbl Mg_Lk(superarbre *st)
{
  int i;

  Mg_Do_Mapping(st);
  Mg_Set_Bl(st->bl,st);  

  st->tree->c_lnL = .0;
  For(i,st->n_gt) 
    {
      st->treelist->tree[i]->both_sides = 1;	  
      Lk(st->treelist->tree[i]);
/*       printf("\n. Tree %3d lnL = %f",i+1,st->treelist->tree[i]->c_lnL); */
      st->tree->c_lnL += st->treelist->tree[i]->c_lnL;
    }
  return st->tree->c_lnL;
}

/*********************************************************/

phydbl Mg_Lk_At_Given_Edge(edge *st_b, superarbre *st)
{
  int i;
  edge *gt_b;
  phydbl lnL;

  Mg_Set_Bl(st->bl,st);

  gt_b = NULL;
  st->tree->c_lnL = .0;
  lnL = .0;
  For(i,st->n_gt)
    {      
      gt_b = st->map_st_edge_in_gt[i][st_b->num];
      lnL = Lk_At_Given_Edge(gt_b,st->treelist->tree[i]);
      st->tree->c_lnL += lnL;
/*       printf("\n. gt %d st edge %d gt edge %d lnL=%f l=%f ",i,st_b->num,gt_b->num,lnL,gt_b->l); */
    }
  return st->tree->c_lnL;
}

/*********************************************************/

phydbl Mg_Update_Lk_At_Given_Edge(edge *st_b, superarbre *st)
{
  int i;
  edge *gt_b;

  Mg_Set_Bl(st->bl,st);
  
  gt_b = NULL;
  st->tree->c_lnL = .0;
  For(i,st->n_gt)
    {
      gt_b = st->map_st_edge_in_gt[i][st_b->num];
      if(gt_b) st->tree->c_lnL += Update_Lk_At_Given_Edge(gt_b,st->treelist->tree[i]);
      else     st->tree->c_lnL += st->treelist->tree[i]->c_lnL;
    }
  return st->tree->c_lnL;
}


/*********************************************************/

void Mg_Fill_Model_Partitions_Table(superarbre *st)
{
  int i,j;
  char *c;
  char *abc;
  int lig, col;
  int n_groups;
  int *encountered_vals;

  c = (char *)mCalloc(10,sizeof(char));
  abc = (char *)mCalloc(20,sizeof(char));
  encountered_vals = (int *)mCalloc(st->n_gt,sizeof(int));

  strcpy(abc,"ABCDEFGHIJKLMNOP\0");

  printf("\n\n\n");
  lig = col = 0;
  while(1)
    {
      printf("\n\n");
      For(i,st->n_gt)
	printf(". Data set %3d : %s\n",i+1,st->optionlist[i]->in_seq_file);

      printf("\n. Data set             ");
      For(i,st->n_gt) printf("%3d ",i+1);
      printf("\n. -A- edge lengths     ");
      For(i,st->n_gt) printf("%3d ",st->bl_partition[i]);

      if(lig == 1) break;

      printf("\n. (%c-%2d)> ",abc[lig],col+1);
      Getstring_Stdin(c);
      
      switch(lig)
	{
	case 0 :
	  {
	    st->bl_partition[col] = atoi(c);
	    break;
	  }
	default :
	  {
	    break;
	  }
	}

      col++;

      if(col == st->n_gt)
	{
	  col = 0;
	  lig++;
	}
    }
    
  n_groups = 0;
  For(i,st->n_gt) 
    {
      For(j,n_groups)
	if(encountered_vals[j] == st->bl_partition[i])
	  break;

      if(j == n_groups) 
	{
	  encountered_vals[n_groups] = st->bl_partition[i];
	  n_groups++;
	}
      st->bl_partition[i] = j;
    }
  
  st->n_bl_partition = n_groups;

  Free(encountered_vals);
  Free(c);
  Free(abc);
}

/*********************************************************/

phydbl Mg_Br_Len_Brent(edge *st_b, superarbre *st)
{
  int iter, n_iter_max;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL, init_lnL;
  phydbl ax, bx, cx;
  int bl_part;
  phydbl *param;
  double tol;

  tol = st->tree->mod->s_opt->min_diff_lk_local;
  n_iter_max = 1000;

  For(bl_part,st->n_bl_partition)
    {
      param = &(st->bl[bl_part][st_b->num]);

      ax = 10.*(*param);
      bx = *param;
      cx = 0.01*(*param);

      d=0.0;
      a=((ax < cx) ? ax : cx);
      b=((ax > cx) ? ax : cx);
      x=w=v=bx;
      old_lnL = UNLIKELY;
      (*param) = fabs(bx);
      fw=fv=fx=fu=-Mg_Lk_At_Given_Edge(st_b,st);
      init_lnL = -fw;
      
      for(iter=1;iter<=BRENT_ITMAX;iter++)
	{
	  xm=0.5*(a+b);
	  tol2=2.0*(tol1=tol*fabs(x)+BRENT_ZEPS);
	  if(
	     ((fabs(st->tree->c_lnL-old_lnL) < tol) && 
	      (st->tree->c_lnL > init_lnL - tol)) ||	 
	     (iter > n_iter_max - 1))
	    {
	      (*param)=x;
	      Mg_Lk_At_Given_Edge(st_b,st);
	      break;
	    }
      
	  if(fabs(e) > tol1)
	    {
	      r=(x-w)*(fx-fv);
	      q=(x-v)*(fx-fw);
	      p=(x-v)*q-(x-w)*r;
	      q=2.0*(q-r);
	      if(q > 0.0) p = -p;
	      q=fabs(q);
	      etemp=e;
	      e=d;
	      if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
		d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	      else
		{
		  d=p/q;
		  u=x+d;
		  if (u-a < tol2 || b-u < tol2)
		    d=SIGN(tol1,xm-x);
		}
	    }
	  else
	    {
	      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	    }
	  u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
	  if(u<BL_MIN) u = BL_MIN;
	  (*param)=fabs(u);
	  old_lnL = st->tree->c_lnL;
	  fu=-Mg_Lk_At_Given_Edge(st_b,st);
	  
	  /*       printf("\n. BRENT edge %3d l=%f lnL=%20f iter=%3d",b_fcus->num,(*param),fu,iter); */

	  if(fu <= fx)
	    {
	      if(u >= x) a=x; else b=x;
	      SHFT(v,w,x,u)
	      SHFT(fv,fw,fx,fu)
	    }
	  else
	    {
	      if (u < x) a=u; else b=u;
	      if (fu <= fw || w == x)
		{
		  v=w;
		  w=u;
		  fv=fw;
		  fw=fu;
		}
	      else if (fu <= fv || v == x || v == w) 
		{
		  v=u;
		  fv=fu;
		}
	    }
	}
      if(iter > BRENT_ITMAX) printf("\n. Too many iterations in BRENT (%d) (%f)",iter,(*param));
      /* Not Reached ??  *xmin=x;   */
      /* Not Reached ??  return fx; */
    }
  return st->tree->c_lnL;
}

/*********************************************************/

void Mg_Initialise_Bl_Partition(superarbre *st)
{
  int i,j;
  
  For(i,st->n_bl_partition)
    {
      For(j,2*st->tree->n_otu-3)
	{
	  st->bl[i][j] = .1;
	}
    }
}

/*********************************************************/

void Mg_Optimize_Br_Len_Serie(node *st_a, node *st_d, edge *st_b, superarbre *st)
{
  phydbl lk_init;
  int i;

  lk_init = st->tree->c_lnL;
  
  Mg_Br_Len_Brent(st_b,st);

  if(st->tree->c_lnL < lk_init - st->tree->mod->s_opt->min_diff_lk_local)
    { 
      printf("\n. %f -- %f",lk_init,st->tree->c_lnL);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
    
  if(st_d->tax) return;
  else For(i,3) if(st_d->v[i] != st_a)
    {
      Mg_Update_P_Lk(st_d->b[i],st_d,st);
      Mg_Optimize_Br_Len_Serie(st_d,st_d->v[i],st_d->b[i],st);
    }

  For(i,3) if((st_d->v[i] == st_a) && (!st_d->v[i]->tax)) Mg_Update_P_Lk(st_d->b[i],st_d,st);

}

/*********************************************************/

void Mg_Update_P_Lk(edge *st_b, node *st_n, superarbre *st)
{
  int i,dir;

  dir = -1;
  For(i,3) if((st_n->b[i]) && (st_n->b[i] == st_b)) {dir = i; break;}
  
  if(dir < 0)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  For(i,st->n_gt)
    {
      if((st->map_st_node_in_gt[i][st_n->num][dir][0]) && (!st->map_st_node_in_gt[i][st_n->num][dir][0]->tax))
	{	  
	  Update_P_Lk(st->treelist->tree[i],
		      st->map_st_edge_in_gt[i][st_b->num],
		      st->map_st_node_in_gt[i][st_n->num][dir][0]);
	}
    }
}

/*********************************************************/

void Mg_Make_N_Swap(edge **st_b, int beg, int end, superarbre *st)
{
  int i;

  st->tree->n_swap = 0;
  for(i=beg;i<end;i++)
    {
      if(st_b[i]->left->tax || st_b[i]->rght->tax)
	{
	  printf("\n. Edge %d is external.",st_b[i]->num);
	  printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      Mg_Swap(st_b[i]->nni->swap_node_v2->v[st->tree->t_dir[st_b[i]->nni->swap_node_v2->num][st_b[i]->nni->swap_node_v1->num]],
	      st_b[i]->nni->swap_node_v2,
	      st_b[i]->nni->swap_node_v3,
	      st_b[i]->nni->swap_node_v3->v[st->tree->t_dir[st_b[i]->nni->swap_node_v3->num][st_b[i]->nni->swap_node_v4->num]],
	      st);

      Swap(st_b[i]->nni->swap_node_v2->v[st->tree->t_dir[st_b[i]->nni->swap_node_v2->num][st_b[i]->nni->swap_node_v1->num]],
	   st_b[i]->nni->swap_node_v2,
	   st_b[i]->nni->swap_node_v3,
	   st_b[i]->nni->swap_node_v3->v[st->tree->t_dir[st_b[i]->nni->swap_node_v3->num][st_b[i]->nni->swap_node_v4->num]],
	   st->tree);

      Mg_Do_Mapping(st);

      st->tree->n_swap++;
    }
}

/*********************************************************/

void Mg_Update_Bl(phydbl fact, superarbre *st)
{
  int i,j;
  
  For(i,2*st->tree->n_otu-3)
    {
      For(j,st->n_gt)
	st->bl[st->bl_partition[j]][i] = 
	st->bl_cpy[st->bl_partition[j]][i] + 
	(st->bl0[st->bl_partition[j]][i] - st->bl_cpy[st->bl_partition[j]][i]) * fact;
    }
  Mg_Set_Bl(st->bl,st);
}

/*********************************************************/

void Mg_Update_Bl_Swaped(edge **st_b, int n, superarbre *st)
{
  int i,j;
  
  For(i,n)
    {
      For(j,st->n_gt)
	{
	  st->bl[st->bl_partition[j]][st_b[i]->num] = 
	    (st_b[i]->nni->best_conf == 1)?
	    (st->bl1[st->bl_partition[j]][st_b[i]->num]):
	    (st->bl2[st->bl_partition[j]][st_b[i]->num]);
	}
    }
  Mg_Set_Bl(st->bl,st);
}

/*********************************************************/

void Mg_Do_Mapping(superarbre *st)
{
  int k;

  Fill_Dir_Table(st->tree);
  For(k,st->n_gt)
    {
      Fill_Dir_Table(st->treelist->tree[k]);
      Mg_Match_St_Nodes_In_Gt(st->treelist->tree[k],st);	      
      Mg_Match_St_Edges_In_Gt(st->treelist->tree[k],st);	      
      Mg_Map_St_Nodes_In_Gt(st->treelist->tree[k],st);
      Mg_Map_St_Edges_In_Gt(st->treelist->tree[k],st);
      Mg_Map_Gt_Edges_In_St(st->treelist->tree[k],st);
    }
}

/*********************************************************/

void Mg_Print_Bl(superarbre *st)
{
  int i,j;
  
  For(j,2*st->tree->n_otu-3)
    { 
      printf("\n. edge %4d ",j);
      For(i,st->n_bl_partition)
	{
	  printf("%f ",st->bl[i][j]);
	}
    }
}

/*********************************************************/
/*********************************************************/

/*  Copyright 2009, Stefan Washietl

    This file is part of RNAcode.

    RNAcode is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RNAcode is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RNAcode.  If not, see <http://www.gnu.org/licenses/>. */


#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
//#include "options.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "rnaz_utils.h"
#include "make.h"

align **Get_Seq_local(align**, option *io,  int rw);

int treeML(const struct aln *alignment[], char** treeString, float* kappa) {

  calign *cdata;
  option *io;
  t_tree *tree;
  int num_data_set;
  int num_tree, num_rand_tree;
  t_mod *mod;
  time_t t_beg, t_end;
  phydbl best_lnL;
  int L,r_seed;
  char *most_likely_tree = NULL;
  int orig_random_input_tree;
  t_opt *s_opt;

  align **data;
  int n_otu;


  /* Initialize data structures */
  tree = NULL;
  mod  = NULL;
  best_lnL = UNLIKELY;


  io    = (option *)Make_Input();
  mod   = (t_mod *)Make_Model_Basic();
  s_opt = (t_opt *)Make_Optimiz();

  Set_Defaults_Input(io);
  Set_Defaults_Model(mod);
  Set_Defaults_Optimiz(s_opt);
  io->mod        = mod;
  mod->io        = io;
  mod->s_opt     = s_opt;


  if (!io) return (0);
  else {
    if (io->use_xml == YES) {
      Free(io);
      return (0);
    }
  }
  r_seed = (io->r_seed < 0) ? (time(NULL)) : (io->r_seed);
  srand(r_seed);
  io->r_seed = r_seed;

  if (io->in_tree == 2) Test_Multiple_Data_Set_Format(io);
  else io->n_trees = 1;

  if (io->n_trees == 0 && io->in_tree == 2)
  {
    PhyML_Printf("\n== Err.: the input tree file does not provide a tree in valid format.");
    Exit("\n");
  }

  if ((io->n_data_sets > 1) && (io->n_trees > 1))
  {
    io->n_data_sets = MIN(io->n_trees, io->n_data_sets);
    io->n_trees     = MIN(io->n_trees, io->n_data_sets);
  }

  for (num_data_set = 0; num_data_set < io->n_data_sets; num_data_set++)
  {
    best_lnL = UNLIKELY;


    n_otu = 0;

    L = strlen(alignment[0]->seq);
    for (n_otu = 0; alignment[n_otu] != NULL; n_otu++);

    io->n_otu = n_otu;

    data = (align **)mCalloc(n_otu, sizeof(align *));

    for (int i = 0; i < n_otu; i++) {
      data[i] = (align *)mCalloc(1, sizeof(align));
      data[i]->len = L;
      data[i]->name = (char *)mCalloc(T_MAX_NAME, sizeof(char));
      strcpy(data[i]->name, alignment[i]->name);
      data[i]->state = (char *)mCalloc(T_MAX_SEQ, sizeof(char));
      strcpy(data[i]->state, alignment[i]->seq);
      data[i]->is_ambigu = NULL;
    }

    /* Call modified version of Get_Seq, that does some processing */
    data = Get_Seq_local(data, io, 0);
    io->data = data;
    //Get_Seq(io);


    Make_Model_Complete(io->mod);
    Set_Model_Name(io->mod);
    Print_Settings(io);
    mod = io->mod;
    orig_random_input_tree = io->mod->s_opt->random_input_tree;


    if (io->data) {
      if (io->n_data_sets > 1) PhyML_Printf("\n. Data set [#%d]\n", num_data_set + 1);
      cdata = Compact_Data(io->data, io);

      Free_Seq(io->data, cdata->n_otu);

      for (num_tree = (io->n_trees == 1) ? (0) : (num_data_set); num_tree < io->n_trees; num_tree++)
      {
        if (io->mod->s_opt->random_input_tree == NO) io->mod->s_opt->n_rand_starts = 1;

        if (orig_random_input_tree == YES && io->n_trees > 1)
        {
          PhyML_Printf("\n== Cannot combine random starting trees with multiple input trees.");
          Exit("\n");
        }

        for (num_rand_tree = 0; num_rand_tree < io->mod->s_opt->n_rand_starts; num_rand_tree++)
        {
          if ((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE))
            if (!io->quiet) PhyML_Printf("\n\n. [Random start %3d/%3d]", num_rand_tree + 1, io->mod->s_opt->n_rand_starts);

          Init_Model(cdata, mod, io);
          switch (io->in_tree)
          {
          case 0 : case 1 : { tree = Dist_And_BioNJ(cdata, mod, io); break; }
          case 2 :          { tree = Read_User_Tree(cdata, mod, io); break; }
          }

          if (io->mod->s_opt->opt_topo == YES) Remove_Duplicates(cdata, io, tree);

          if (io->fp_in_constraint_tree != NULL)
          {
            char *s;

            PhyML_Printf("\n. Reading constraint tree file...");

            io->cstr_tree = Read_Tree_File_Phylip(io->fp_in_constraint_tree);

            if (io->cstr_tree->n_root != NULL)
            {
              PhyML_Printf("\n== The constraint tree file must be unrooted");
              Exit("\n");
            }

            s = Add_Taxa_To_Constraint_Tree(io->fp_in_constraint_tree, cdata);
            fflush(NULL);
            Free_Tree(tree);
            tree = Read_Tree(&s);
            io->in_tree = 2;
            Free(s);
            Check_Constraint_Tree_Taxa_Names(io->cstr_tree, cdata);
            Alloc_Bip(io->cstr_tree);
            Get_Bip(io->cstr_tree->a_nodes[0],
                    io->cstr_tree->a_nodes[0]->v[0],
                    io->cstr_tree);
            if (tree->has_branch_lengths == NO) Add_BioNJ_Branch_Lengths(tree, cdata, mod, NULL);
          }


          if (!tree) continue;

          time(&t_beg);
          time(&(tree->t_beg));

          tree->mod          = mod;
          tree->io           = io;
          tree->data         = cdata;
          tree->n_pattern    = tree->data->crunch_len;
          tree->n_root       = NULL;
          tree->e_root       = NULL;
          tree->n_tot_bl_opt = 0;

          Set_Both_Sides(YES, tree);

          if ((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);

          if (io->cstr_tree && !Check_Topo_Constraints(tree, io->cstr_tree))
          {
            PhyML_Printf("\n\n== The initial tree does not satisfy the topological constraint.");
            PhyML_Printf("\n== Please use the user input tree option with an adequate tree topology.");
            Exit("\n");
          }

          Connect_CSeqs_To_Nodes(tree->data, tree->io, tree);
          Make_Tree_For_Pars(tree);
          Make_Tree_For_Lk(tree);
          Make_Spr(tree);
          Br_Len_Not_Involving_Invar(tree);
          Unscale_Br_Len_Multiplier_Tree(tree);
          if (tree->io->print_json_trace == YES) JSON_Tree_Io(tree, tree->io->fp_out_json_trace);



          Set_Update_Eigen(YES, tree->mod);
          Lk(NULL, tree);
          Set_Update_Eigen(NO, tree->mod);


          if (tree->mod->s_opt->opt_topo) {
            Global_Spr_Search(tree);
            if (tree->n_root) Add_Root(tree->a_edges[0], tree);
          }
          else
          {
            if (tree->mod->s_opt->opt_subst_param || tree->mod->s_opt->opt_bl) Round_Optimize(tree, ROUND_MAX);
          }


          if (tree->mod->gamma_mgf_bl) Best_Root_Position_IL_Model(tree);

          Set_Both_Sides(YES, tree);
          Lk(NULL, tree);
          Pars(NULL, tree);
          Get_Tree_Size(tree);
          PhyML_Printf("\n\n. Log likelihood of the current tree: %.*f.", DECIMAL_DIG, tree->c_lnL);

          if (tree->io->ancestral == YES) Ancestral_Sequences(tree, YES);

          Check_Br_Lens(tree);
          Br_Len_Involving_Invar(tree);
          Rescale_Br_Len_Multiplier_Tree(tree);
          if (!tree->n_root) Get_Best_Root_Position(tree);

          /* Print the tree estimated using the current random (or BioNJ) starting tree */
          /* if(io->mod->s_opt->n_rand_starts > 1) */
          if (orig_random_input_tree == YES)
          {
            Print_Tree(io->fp_out_trees, tree);
            fflush(NULL);
          }

          /* Record the most likely tree in a string of characters */
          if (tree->c_lnL > best_lnL)
          {
            best_lnL = tree->c_lnL;
            if (most_likely_tree) Free(most_likely_tree);
            most_likely_tree = Write_Tree(tree);

            time(&t_end);
            /*
            Print_Fp_Out(io->fp_out_stats, t_beg, t_end, tree,
                         io, num_data_set + 1,
                         (orig_random_input_tree == YES) ? (num_rand_tree) : (num_tree),
                         (num_rand_tree == io->mod->s_opt->n_rand_starts - 1) ? (YES) : (NO), io->precision);
            */
            if (tree->io->print_site_lnl) Print_Site_Lk(tree, io->fp_out_lk);
          }



          /* Start from BioNJ tree */
          if ((num_rand_tree == io->mod->s_opt->n_rand_starts - 1) && (tree->mod->s_opt->random_input_tree))
          {
            /* Do one more iteration in the loop, but don't randomize the tree */
            tree->mod->s_opt->n_rand_starts++;
            tree->mod->s_opt->random_input_tree = NO;
          }
          Free_One_Spr(tree->best_spr);
          Free_Spr_List_One_Edge(tree);
          Free_Spr_List_All_Edge(tree);
          Free_Tree_Pars(tree);
          Free_Tree_Lk(tree);
          Free_Tree(tree);
        } //Tree done

        if (io->n_data_sets == 1) rewind(io->fp_out_tree);
        if (most_likely_tree) PhyML_Fprintf(io->fp_out_tree, "%s\n", most_likely_tree);


        /* Launch bootstrap analysis */
        if (io->do_boot || io->do_tbe)
        {
          if (!io->quiet) PhyML_Printf("\n\n. Launch bootstrap analysis on the most likely tree...");
          most_likely_tree = Bootstrap_From_String(most_likely_tree, cdata, mod, io);

          PhyML_Printf("\n\n. Completed the bootstrap analysis succesfully."); fflush(NULL);
        }
        else if (io->ratio_test != NO)
        {
          /* Launch aLRT */
          most_likely_tree = aLRT_From_String(most_likely_tree, cdata, mod, io);
        }


        /* Print the most likely tree in the output file */
        if (!io->quiet) PhyML_Printf("\n\n. Printing the most likely tree in file '%s'.", Basename(io->out_tree_file));
        if (io->n_data_sets == 1) rewind(io->fp_out_tree);

        t_tree *dum;
        dum = Read_Tree(&most_likely_tree);
        dum->data = cdata;
        dum->mod  = mod;
        dum->io   = io;
        Connect_CSeqs_To_Nodes(cdata, io, dum);
        Insert_Duplicates(dum);
        Free(most_likely_tree);
        most_likely_tree = Write_Tree(dum);
        Free_Tree(dum);

        PhyML_Fprintf(io->fp_out_tree, "%s\n", most_likely_tree);

        if (io->n_trees > 1 && io->n_data_sets > 1) break;
      }
      Free_Calign(cdata);
    }
    else
    {
      PhyML_Printf("\n== No data was found.\n");
      PhyML_Printf("\n== Err. in file %s at line %d\n", __FILE__, __LINE__);
      Exit("\n");
    }
    Free_Model_Complete(mod);
  }

  if (most_likely_tree) Free(most_likely_tree);

  if (mod->s_opt->n_rand_starts > 1) PhyML_Printf("\n. Best log likelihood: %f\n", best_lnL);

  Free_Optimiz(mod->s_opt);
  Free_Model_Basic(mod);

  if (io->fp_in_constraint_tree) fclose(io->fp_in_constraint_tree);
  if (io->fp_in_align)           fclose(io->fp_in_align);
  if (io->fp_in_tree)            fclose(io->fp_in_tree);
  if (io->fp_out_lk)             fclose(io->fp_out_lk);
  if (io->fp_out_tree)           fclose(io->fp_out_tree);
  if (io->fp_out_trees)          fclose(io->fp_out_trees);
  if (io->fp_out_stats)          fclose(io->fp_out_stats);
  if (io->fp_out_trace)          fclose(io->fp_out_trace);
  if (io->fp_out_json_trace)     fclose(io->fp_out_json_trace);

  if (io->fp_in_constraint_tree != NULL) Free_Tree(io->cstr_tree);
  Free_Input(io);
  return 1;

}

/* Function from Seqgen (c)  Andrew Rambaut & Nick Grassly */

align **Get_Seq_local(align** data, option *io,  int rw)
{
  int i, j;
  char **buff;
  int n_unkn, n_removed, pos;
  int *remove;


  /*   rewind(fp_seq); */


  if (data)
  {
    buff = (char **)mCalloc(io->n_otu, sizeof(char *));
    For(i, io->n_otu) buff[i] = (char *)mCalloc(data[0]->len, sizeof(char));
    remove = (int *)mCalloc(data[0]->len, sizeof(int));

    n_removed = 0;

    For(i, data[0]->len)
    {
      For(j, io->n_otu)
      {
        if ((data[j]->state[i] == '?') || (data[j]->state[i] == '-')) data[j]->state[i] = 'X';
        if ((io->datatype == NT) && (data[j]->state[i] == 'N')) data[j]->state[i] = 'X';
        if (data[j]->state[i] == 'U') data[j]->state[i] = 'T';
      }

      n_unkn = 0;
      For(j, io->n_otu) if (data[j]->state[i] == 'X') n_unkn++;

      if (n_unkn == io->n_otu)
      {
        remove[i] = 1;
        n_removed++;
      }

      For(j, io->n_otu) buff[j][i] = data[j]->state[i];
    }

    if (n_removed > 0)
    {
      if (io->datatype == NT)
        printf("\n. %d sites are made from completely undetermined states ('X', '-', '?' or 'N')...\n", n_removed);
      else
        printf("\n. %d sites are made from completely undetermined states ('X', '-', '?')...\n", n_removed);
    }

    pos = 0;
    For(i, data[0]->len)
    {
      /*    if(!remove[i]) */
      /*      { */
      For(j, io->n_otu) data[j]->state[pos] = buff[j][i];
      pos++;
      /*      } */
    }

    For(i, io->n_otu) data[i]->len = pos;
    For(i, io->n_otu) Free(buff[i]);
    Free(buff);
    Free(remove);
  }
  return data;
}


/*



  Make_Model_Complete(io->mod);
  mod = io->mod;


  //Set options

  io->datatype = NT; // Nucleotides not amino acids
  io->mod->s_opt->opt_topo = 0; // Do not optimize topology
  io->mod->s_opt->opt_bl   = 1; // Optimize branch lengths
  io->mod->s_opt->opt_subst_param = 1; // Optimize kappa, initialize with 4.0
  io->mod->s_opt->opt_kappa = 1;
  io->mod->kappa->v = 4.0;

  // Manually fill data structure with the alignment

  n_otu = 0;

  L = strlen(alignment[0]->seq);
  for (n_otu = 0; alignment[n_otu] != NULL; n_otu++);

  io->n_otu = n_otu;

  data = (align **)mCalloc(n_otu, sizeof(align *));

  for (i = 0; i < n_otu; i++) {
    data[i] = (align *)mCalloc(1, sizeof(align));
    data[i]->len = L;
    data[i]->name = (char *)mCalloc(T_MAX_NAME, sizeof(char));
    strcpy(data[i]->name, alignment[i]->name);
    data[i]->state = (char *)mCalloc(T_MAX_SEQ, sizeof(char));
    strcpy(data[i]->state, alignment[i]->seq);
    data[i]->is_ambigu = NULL;
  }

  // Call modified version of Get_Seq, that does some processing
  data = Get_Seq_local(data, io, 0);
  alldata = Compact_Data(data, io);
  Free_Seq(data, alldata->n_otu);
  Check_Ambiguities(alldata, io->datatype, io->state_len);
  Init_Model(alldata, mod, io);


  // Calculate pairwise distances and make BIONJ tree
  mat = ML_Dist(alldata, mod);
  Fill_Missing_Dist(mat);
  mat->tree = Make_Tree_From_Scratch(alldata->n_otu, alldata);


  // bionj.h

  Bionj(mat);
  tree = mat->tree;
  tree->verbose = 0; // Shut of verbose output
  tree->mod = mod;
  tree->io = io;
  tree->data = alldata;
  tree->both_sides = 1;
  tree->n_pattern = tree->data->crunch_len / io->state_len;
  //time(&t_beg);
  //time(&(tree->t_beg));

  // Prepare for optimization and optimize

  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  Connect_CSeqs_To_Nodes(tree->data, tree->io, tree);
  Make_Tree_For_Pars(tree);
  Make_Tree_For_Lk(tree);
  //tree->triplet_struct = Make_Triplet_Struct(mod);
  Br_Len_Not_Involving_Invar(tree);
  //Connect_CSeqs_To_Nodes(alldata,io,tree);
  Round_Optimize(tree, ROUND_MAX);

  //treeString=(char *)mCalloc(T_MAX_LINE,sizeof(char));

  *treeString = Write_Tree(tree);
  //what that for??
  // *kappa=io->mod->kappa;

  Free_Mat(tree->mat);
  //Free_Triplet(tree->triplet_struct);
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Tree(tree);

  Free_Calign(alldata);

  Free_Model(mod);

  Free_Input(io);

  return (1);
  */
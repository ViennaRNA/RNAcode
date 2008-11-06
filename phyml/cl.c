/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <unistd.h>
#include <getopt.h>
#include "utilities.h"
#include "options.h"
#include "cl.h"
#include "models.h"
#include "free.h"
#include "interface.h"
#include "mg.h"
#include "m4.h"


/*********************************************************/
/**
* Fill the Option fields, with the argc array
*/
void Read_Command_Line(option *io, int argc, char **argv)
{
  int c;
  int open_ps_file;
  int use_gamma;
  struct option longopts[] =
    {
      {"n_rgrft",           required_argument,NULL,0},
      {"n_globl",           required_argument,NULL,1},
      {"max_dist",          required_argument,NULL,2},
      {"n_optim",           required_argument,NULL,3},
      {"n_best",            required_argument,NULL,4},
      {"model",             required_argument,NULL,5},
      {"search",            required_argument,NULL,6},
      {"datatype",          required_argument,NULL,7},
      {"multiple",          required_argument,NULL,8},
      {"input",             required_argument,NULL,9},
      {"bootstrap",         required_argument,NULL,10},
      {"ts/tv",             required_argument,NULL,11},
      {"nclasses",          required_argument,NULL,12},
      {"pinv",              required_argument,NULL,13},
      {"alpha",             required_argument,NULL,14},
      {"inputtree",         required_argument,NULL,15},
      {"min_diff_lk_local", required_argument,NULL,16},
      {"min_diff_lk_global",required_argument,NULL,17},
      {"steph_spr",         no_argument,NULL,18},
      {"brent_it_max",      required_argument,NULL,19},
      {"rand_start",        no_argument,NULL,20},
      {"n_rand_starts",     required_argument,NULL,21},
      {"sequential",        no_argument,NULL,22},
      {"inside_opt",        no_argument,NULL,23},
      {"p_moves",           required_argument,NULL,24},
      {"fast_nni",          no_argument,NULL,25},
      {"g_pars",            no_argument,NULL,26},
      {"r_seed",            required_argument,NULL,27},
      {"collapse_boot",     required_argument,NULL,28},
      {"random_boot",       required_argument,NULL,29},
      {"print_trace",       no_argument,NULL,30},
      {"print_site_lnl",    no_argument,NULL,31},
      {"cov",               no_argument,NULL,32},
      {"cov_delta",         required_argument,NULL,33},
      {"cov_alpha",         required_argument,NULL,34},
      {"cov_ncats",         required_argument,NULL,35},
      {"ps",                no_argument,NULL,36},
      {"cov_free",          no_argument,NULL,37},
      {"no_gap",            no_argument,NULL,38},
      {"n_rr_branch",       required_argument,NULL,39},
      {0,0,0,0}
    };

  open_ps_file = 0;
  use_gamma = 0;
  while((c = getopt_long(argc,argv,"qi:d:m:b:n:t:f:zk:v:c:a:u:ho:s:x:g:l:e",longopts,NULL)) != -1)
      {
	switch(c)
	  {
	  case 39 :
	    {
	      io->mod->n_rr_branch = (int)atoi(optarg);
	      if(io->mod->n_rr_branch < 1)
		{
		  printf("\n. The number of classes must be an integer greater than 0.\n");
		  printf("\n. Type any key to exit.\n");
		  Exit("\n");
		}
	      break;
	    }
	  case 38 :
	    {
	      io->rm_ambigu = 1;
	      break;
	    }
	  case 37 :
	    {
	      io->mod->s_opt->opt_cov_free_rates = 1;
	      break;
	    }
	  case 36 :
	    {
	      open_ps_file = 1;
	      break;
	    }
	  case 35 :
	    {
	      io->m4_model = YES;
	      if(!io->mod->m4mod) io->mod->m4mod = (m4 *)M4_Make_Light((io->mod->datatype == NT)?(4):(20));
	      io->mod->m4mod->n_h = (int)atoi(optarg);
	      
	      if(io->mod->m4mod->n_h < 1)
		{
		  char choix;
		  printf("\n. The number of classes must be greater than 0.\n");
		  printf("\n. Type any key to exit.\n");
		  scanf("%c",&choix);
		  Exit("\n");
		}
	      break;
	    }
	  case 34 :
	    {
	      io->m4_model = YES;
	      if(!io->mod->m4mod) io->mod->m4mod = (m4 *)M4_Make_Light((io->mod->datatype == NT)?(4):(20));
	      io->mod->m4mod->use_cov_alpha = 1;
	      io->mod->m4mod->use_cov_free  = 0;
	      
	      if(!strcmp(optarg,"e") || !strcmp(optarg,"E") ||
		 !strcmp(optarg,"estimated") || !strcmp(optarg,"ESTIMATED"))
		{
		  io->mod->s_opt->opt_cov_alpha = 1;
		  io->mod->m4mod->alpha         = 1.0;
		}
	      else
		{
		  io->mod->m4mod->alpha = (phydbl)atof(optarg);

		  if(io->mod->m4mod->alpha < 1.E-5)
		    {
		      char choix;
		      printf("\n. The value of alpha must be greater than 1.E-5.\n");
		      printf("\n. Type any key to exit.\n");
		      scanf("%c",&choix);
		      Exit("\n");
		    }
		}
	      break;
	    }

	  case 33 :
	    {
	      io->m4_model = YES;
	      if(!io->mod->m4mod) io->mod->m4mod = (m4 *)M4_Make_Light((io->mod->datatype == NT)?(4):(20));

	      if(!strcmp(optarg,"e") || !strcmp(optarg,"E") ||
		 !strcmp(optarg,"estimated") || !strcmp(optarg,"ESTIMATED"))
		{
		  io->mod->s_opt->opt_cov_delta = 1;
		  io->mod->m4mod->delta         = 1.0;
		}
	      else
		{
		  io->mod->m4mod->delta = (phydbl)atof(optarg);

		  if(atof(optarg) < 1.E-10)
		    {
		      char choix;
		      printf("\n. The value of delta must be larger than 1.E-10.\n");
		      printf("\n. Type any key to exit.\n");
		      scanf("%c",&choix);
		      Exit("\n");
		    }
		}
	      break;
	    }
	  case 32 :
	    {
	      io->m4_model = YES;
	      break;
	    }
	  case 31 :
	    {
	      io->print_site_lnl = 1;
	      break;
	    }
	  case 30 :
	    {
	      io->print_trace = 1;
	      break;
	    }
	  case 29 :
	    {
	      io->random_boot_seq_order = (int)atoi(optarg);
	      break;
	    }
	  case 28 :
	    {
	      io->collapse_boot = (int)atoi(optarg);
	      break;
	    }
	  case 27 :
	    {
	      io->r_seed = (int)atoi(optarg);
	      break;
	    }
	  case 26 :
	    {
	      io->mod->s_opt->general_pars = 1;
	      break;
	    }
	  case 25 :
	    {
	      io->mod->s_opt->fast_nni = 1;
	      break;
	    }
	  case 24 :
	    {
	      io->mod->s_opt->p_moves_to_examine = (double)atof(optarg);
	      break;
	    }
	  case 23 :
	    {
	      io->mod->s_opt->wim_inside_opt = 1;
	      break;
	    }
	  case 0 :
	    {
	      io->mod->s_opt->wim_n_rgrft = atoi(optarg);
	      break;
	    }
	  case 1 :
	    {
	      io->mod->s_opt->wim_n_globl = atoi(optarg);
	      break;
	    }
	  case 2 :
	    {
	      io->mod->s_opt->wim_max_dist = atoi(optarg);
	      break;
	    }
	  case 3 :
	    {
	      io->mod->s_opt->wim_n_optim = atoi(optarg);
	      break;
	    }
	  case 4 :
	    {
	      io->mod->s_opt->wim_n_best = atoi(optarg);
	      break;
	    }
	  case 16 :
	    {
	      io->mod->s_opt->min_diff_lk_local = atof(optarg);
	      break;
	    }
	    
	  case 17 :
	    {
	      io->mod->s_opt->min_diff_lk_global = atof(optarg);
	      break;
	    }
	  case 18 :
	    {
	      io->mod->s_opt->steph_spr = 0;
	      io->mod->s_opt->greedy    = 1;
	      break;
	    }
	  case 19 :
	    {
	      io->mod->s_opt->brent_it_max = atoi(optarg);
	      break;
	    }
	  case 20 :
	    {
	      io->mod->s_opt->random_input_tree = 1;
	      break;
	    }
	  case 21 :
	    {
	      io->mod->s_opt->random_input_tree = 1;
	      io->mod->s_opt->n_rand_starts = atoi(optarg);
	    }
	  case 's':case 6:
	    {
	      if((!strcmp(optarg,"spr")) || (!strcmp(optarg,"SPR")))
		{
		  io->mod->s_opt->topo_search = SPR_MOVE;
		  io->mod->s_opt->greedy      = (io->mod->s_opt->steph_spr)?(0):(1);
		}
	      else if((!strcmp(optarg,"nni")) || (!strcmp(optarg,"NNI")))
		{
		  io->mod->s_opt->topo_search         = NNI_MOVE;
		  io->mod->s_opt->random_input_tree   = 0;
		}
	      else if((!strcmp(optarg,"nni+spr")) || (!strcmp(optarg,"NNI+SPR")))
		{
		  io->mod->s_opt->topo_search         = NNI_MOVE;
		  io->mod->s_opt->random_input_tree   = 0;
		  io->mod->s_opt->spr_step_after_nnis = 1;
		  io->mod->s_opt->greedy              = (io->mod->s_opt->steph_spr)?(0):(1);
		}
	      break;
	    }
	    
	  case 'd':case 7:
	    if(!strcmp(optarg,"nt"))
	      /*         if(atoi(optarg) == NT) */
	      {
		io->mod->datatype = NT;
		io->mod->stepsize = 1;
		io->mod->ns = 4;
		
		if(
		   (io->mod->whichmodel == LG)       ||
		   (io->mod->whichmodel == WAG)       ||
		   (io->mod->whichmodel == DAYHOFF)   ||
		   (io->mod->whichmodel == JTT)       ||
		   (io->mod->whichmodel == BLOSUM62)  ||
		   (io->mod->whichmodel == MTREV)     ||
		   (io->mod->whichmodel == RTREV)     ||
		   (io->mod->whichmodel == CPREV)     ||
		   (io->mod->whichmodel == DCMUT)     ||
		   (io->mod->whichmodel == VT)        ||
		   (io->mod->whichmodel == MTMAM)     ||
		   (io->mod->whichmodel == MTART)     ||
		   (io->mod->whichmodel == HIVW)      ||
		   (io->mod->whichmodel == HIVB)      ||
		   (io->mod->whichmodel == CUSTOMAA)
		   )
		  {
		    io->mod->whichmodel = HKY85;
		    strcpy(io->mod->modelname, "HKY85\0");
		  }
	      }
	    else if (!strcmp(optarg,"aa"))
	      /*         else if (atoi(optarg) == AA) */
	      {
		io->mod->datatype = AA;
		io->mod->stepsize = 1;
		io->mod->ns = 20;
		if(
		   (io->mod->whichmodel == JC69)   ||
		   (io->mod->whichmodel == K80)    ||
		   (io->mod->whichmodel == F81)    ||
		   (io->mod->whichmodel == HKY85)  ||
		   (io->mod->whichmodel == F84)    ||
		   (io->mod->whichmodel == TN93)   ||
		   (io->mod->whichmodel == GTR)    ||
		   (io->mod->whichmodel == CUSTOM)
		   )
		  {
		    io->mod->whichmodel = LG;
		    strcpy(io->mod->modelname, "LG\0");
		  }
	      }
	    else
	      {
		char choix;
		printf("\n. Unknown argument to -d option: please use `nt' for DNA or `aa' for Amino-Acids\n");
		printf("\n. Type any key to exit.\n");
		scanf("%c",&choix);
		Exit("\n");
	      }
	    
	    break;
	    
	  case 'm': case 5 :
	    {
	      if(!isalpha(optarg[0]))
		{
		  strcpy(io->mod->custom_mod_string,optarg);

		  if(strlen(io->mod->custom_mod_string) != 6)
		    {
		      Warn_And_Exit("\n. The string should be of length 6.\n");
		    }
		  else
		    {
		      Translate_Custom_Mod_String(io->mod);
		    }

		  io->mod->datatype = NT;
		  io->mod->whichmodel = CUSTOM;
		  strcpy(io->mod->modelname, "custom");
		  io->mod->s_opt->opt_kappa     = 0;
		  io->mod->s_opt->opt_rr        = 1;
		  io->mod->s_opt->opt_num_param = 1;
		}

	      if (strcmp(optarg, "JC69") == 0)
		{
		  io->mod->datatype = NT;
		  io->mod->whichmodel = JC69;
		  strcpy(io->mod->modelname, "JC69");
		  io->mod->s_opt->opt_kappa  = 0;
		}
	      else if(strcmp(optarg, "K80") == 0)
		{
		  io->mod->datatype = NT;
		  io->mod->whichmodel = K80;
		  strcpy(io->mod->modelname, "K80");
		}
	      else if(strcmp(optarg, "F81") == 0)
		{
		  io->mod->datatype = NT;
		  io->mod->whichmodel = F81;
		  strcpy(io->mod->modelname, "F81");
		  io->mod->s_opt->opt_kappa  = 0;
		}
	      else if (strcmp(optarg, "HKY85") == 0)
		{
		  io->mod->datatype = NT;
		  io->mod->whichmodel = HKY85;
		  strcpy(io->mod->modelname, "HKY85");
		}
	      else if(strcmp(optarg, "F84") == 0)
		{
		  io->mod->datatype = NT;
		  io->mod->whichmodel = F84;
		  strcpy(io->mod->modelname, "F84");
		}
	      else if (strcmp (optarg,"TN93") == 0)
		{
		  io->mod->datatype = NT;
		  io->mod->whichmodel = TN93;
		  strcpy(io->mod->modelname, "TN93");
		  if(io->mod->s_opt->opt_kappa)
		    {
		      io->mod->s_opt->opt_lambda = 1;
		    }
		}
	      else if(strcmp (optarg, "GTR") == 0)
		{
		  io->mod->datatype = NT;
		  io->mod->whichmodel = GTR;
		  strcpy(io->mod->modelname, "GTR");
		  io->mod->s_opt->opt_kappa = 0;
		}
	      else if(strcmp(optarg, "Dayhoff") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = DAYHOFF;
		  strcpy(io->mod->modelname, "Dayhoff");
		}
	      else if(strcmp (optarg, "JTT") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = JTT;
		  strcpy(io->mod->modelname, "JTT");
		}
	      else if(strcmp(optarg, "MtREV") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = MTREV;
		  strcpy(io->mod->modelname,"MtREV");
		}
	      else if(strcmp (optarg, "LG") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = LG;
		  strcpy(io->mod->modelname, "LG");
		}
	      else if(strcmp (optarg, "WAG") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = WAG;
		  strcpy(io->mod->modelname, "WAG");
		}
	      else if(strcmp(optarg, "DCMut") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = DCMUT;
		  strcpy(io->mod->modelname, "DCMut");
		}
	      else if(strcmp (optarg, "RtREV") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = RTREV;
		  strcpy(io->mod->modelname, "RtREV");
		}
	      else if(strcmp(optarg, "CpREV") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = CPREV;
		  strcpy(io->mod->modelname, "CpREV");
		}
	      else if(strcmp(optarg, "VT") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = VT;
		  strcpy(io->mod->modelname, "VT");
		}
	      else if(strcmp(optarg, "Blosum62") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = BLOSUM62;
		  strcpy(io->mod->modelname,"Blosum62");
		}
	      else if(strcmp(optarg, "MtMam") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = MTMAM;
		  strcpy(io->mod->modelname, "MtMam");
		}
	      else if (strcmp(optarg,"MtArt") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = MTART;
		  strcpy(io->mod->modelname, "MtArt");
		}
	      else if (strcmp(optarg,"HIVw") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = HIVW;
		  strcpy(io->mod->modelname, "HIVw");
		}
	      else if(strcmp(optarg, "HIVb") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = HIVB;
		  strcpy(io->mod->modelname, "HIVb");
		}
	      else if (strcmp(optarg, "custom") == 0)
		{
		  io->mod->datatype = AA;
		  io->mod->whichmodel = CUSTOMAA;
		  strcpy(io->mod->modelname, "Read from file");
		}
	      break;
	    }

	  case 'a':case 14 :
	    {
	      use_gamma = 1;
	      if ((strcmp (optarg, "e") == 0) ||
		  (strcmp (optarg, "E") == 0) ||
		  (strcmp (optarg, "estimated") == 0) ||
		  (strcmp (optarg, "ESTIMATED") == 0))
		{
		  io->mod->s_opt->opt_alpha     = 1;
		  io->mod->s_opt->opt_num_param = 1;
		}
	      else if ((!atof(optarg)) || (atof(optarg) < .0))
		{
		  char choix;
		  printf("\n. Alpha must be a positive number\n");
		  printf("\n. Type any key to exit.\n");
		  scanf("%c",&choix);
		  Exit("\n");
		}
	      else
		{
		  io->mod->alpha = (phydbl)atof(optarg);
		  io->mod->s_opt->opt_alpha  = 0;
		}
	      break;
	    }
	  case 'b':case 10:
	  {
	    if (atoi(optarg) < -4)
	      {
		char choix;
		printf("\n. Branch test value must be a positive integer for bootstrap, or between -1 and -4 for aLRT branch test\n");
		printf("\n. Type any key to exit.\n");
		scanf("%c",&choix);
		Exit("\n");
	      }
	    else
	      {
		if((int)atoi(optarg) > 0)
		  {
		    io->ratio_test       = 0;
		    io->mod->bootstrap   = (int)atoi(optarg);
		    io->print_boot_trees = 1;

		    if(io->n_data_sets > 1)
		      {
			char choix;
			printf("\n. Bootstrap option is not allowed with multiple data sets\n");
			printf("\n. Type any key to exit.\n");
			scanf("%c",&choix);
			Exit("\n");
		      }
		  }
		else if (atoi(optarg)==0)
		  {
		    io->mod->bootstrap = 0;
		    io->ratio_test     = 0;
		  }
		else
		  {
		    io->mod->bootstrap = 0;
		    io->ratio_test     = -(int)atoi(optarg);
		  }
	      }
	    break;
	  }
	  case 'c':case 12:
	    {
	      if ((!atoi(optarg)) || (atoi(optarg) < 0))
		{
		  char choix;
		  printf("\n. Unknown argument to -c option: the number of categories must be a positive integer\n");
		  printf("\n. Type any key to exit.\n");
		  scanf("%c",&choix);
		  Exit("\n");
		}
	      else io->mod->n_catg = atoi(optarg);
	      break;
	    }
	case 'f':
	  {
	    if(!strcmp(optarg,"e"))
	      {
		io->mod->s_opt->opt_state_freq = 1;
		
		if((io->mod->whichmodel == JC69) ||
		   (io->mod->whichmodel == K80))
		  {
		    char choix;
		    printf("\n. Invalid model settings (option '-f').\n");
		    printf("\n. Type any key to exit.\n");
		    scanf("%c",&choix);
		    Exit("\n");
		  }
	      }
	    else if(!strcmp(optarg,"d"))
	      {
		io->mod->s_opt->opt_state_freq = 0;
	      }
	    else if(!isalpha(optarg[0]))
	      {
		phydbl sum;

		io->mod->s_opt->opt_state_freq  = 0;
		io->mod->s_opt->user_state_freq = 1;
		
		sum =
		  (io->mod->user_b_freq[0] +
		   io->mod->user_b_freq[1] +
		   io->mod->user_b_freq[2] +
		   io->mod->user_b_freq[3]);

		io->mod->user_b_freq[0] /= sum;
		io->mod->user_b_freq[1] /= sum;
		io->mod->user_b_freq[2] /= sum;
		io->mod->user_b_freq[3] /= sum;

		if(io->mod->user_b_freq[0] < .0 ||
		   io->mod->user_b_freq[1] < .0 ||
		   io->mod->user_b_freq[2] < .0 ||
		   io->mod->user_b_freq[3] < .0 ||
		   io->mod->user_b_freq[0] > 1. ||
		   io->mod->user_b_freq[1] > 1. ||
		   io->mod->user_b_freq[2] > 1. ||
		   io->mod->user_b_freq[3] > 1.)
		  {
		    Warn_And_Exit("\n. Invalid base frequencies.\n");
		  }
	      }
	    break;
	  }
	  
	case 'h':
	  {
	    Usage();
	    break;
	  }

	  case 'i':case 9:
	    {
	      char *tmp;
	      tmp = (char *) mCalloc (T_MAX_FILE, sizeof(char));
	      if (strlen (optarg) > T_MAX_FILE -16)
		{
		  char choix;
		  strcpy (tmp, "\n. The file name'");
		  strcat (tmp, optarg);
		  strcat (tmp, "' is too long.\n");
		  printf("%s",tmp);
		  printf("\n. Type any key to exit.\n");
		  scanf("%c",&choix);
		  Exit("\n");
		}
	  
	      else if (!Filexists (optarg))
		{
		  char choix;
		  strcpy (tmp, "\n. The file '");
		  strcat (tmp, optarg);
		  strcat (tmp, "' does not exist.\n");
		  printf("%s",tmp);
		  printf("\n. Type any key to exit.\n");
		  scanf("%c",&choix);
		  Exit("\n");
		}
	      else
		{
		  strcpy(io->in_seq_file, optarg);
		  io->fp_in_seq = Openfile(io->in_seq_file,0);
		  strcpy(io->out_tree_file,optarg);
#ifdef PHYML
		  strcat(io->out_tree_file,"_phyml_tree.txt");
#else
		  strcat(io->out_tree_file,"_mc_tree.txt");
#endif
		  io->fp_out_tree = Openfile(io->out_tree_file,1);
		  strcpy(io->out_stats_file,optarg);
#ifdef PHYML
		  strcat(io->out_stats_file,"_phyml_stats.txt");
#else
		  strcat(io->out_stats_file,"_mc_stats.txt");
#endif
		  io->fp_out_stats = Openfile(io->out_stats_file,1);
		}
	      Free (tmp);
	      break;
	    }
	  
	  case 't':case 11:
	    {
	      if ((io->mod->whichmodel != JC69) &&
		  (io->mod->whichmodel != F81)  &&
		  (io->mod->whichmodel != GTR))
		{
		  if ((strcmp(optarg, "e") == 0) ||
		      (strcmp(optarg, "E") == 0) ||
		      (strcmp(optarg, "estimated") == 0) ||
		      (strcmp(optarg, "ESTIMATED") == 0))
		    {
		      io->mod->kappa                 = 4.0;
		      io->mod->s_opt->opt_num_param  = 1;
		      io->mod->s_opt->opt_kappa      = 1;
		      if (io->mod->whichmodel == TN93)
			io->mod->s_opt->opt_lambda   = 1;
		    }
		  else
		    {
		      if ((!atof(optarg)) || (atof(optarg) < .0))
			{
			  char choix;
			  printf("\n. The ts/tv ratio must be a positive number\n");
			  printf("\n. Type any key to exit.\n");
			  scanf("%c",&choix);
			  Exit("\n");
			}
		      else
			{
			  io->mod->kappa = (phydbl)atof(optarg);
			  io->mod->s_opt->opt_kappa  = 0;
			  io->mod->s_opt->opt_lambda = 0;
			}
		    }
		}
	      break;
	    }
	  case 'n':case 8:
	    {
	      if ((!atoi(optarg)) || (atoi(optarg) < 0))
		{
		  char choix;
		  printf("\n. The number of alignments must be a positive integer\n");
		  printf("\n. Type any key to exit.\n");
		  scanf("%c",&choix);
		  Exit("\n");
		}
	      else io->n_data_sets = atoi (optarg);
	      break;
	    }
	  case 'q':case 22:
	    {
	      io->interleaved = 0;
	      break;
	    }
	  case 'u':case 15:
	    {
	      char *tmp;
	      tmp = (char *)mCalloc(T_MAX_FILE, sizeof(char));
	      if (strlen (optarg) > T_MAX_FILE -11)
		{
		  char choix;
		  strcpy (tmp, "\n. The file name'");
		  strcat (tmp, optarg);
		  strcat (tmp, "' is too long.\n");
		  printf("%s",tmp);
		  printf("\n. Type any key to exit.\n");
		  scanf("%c",&choix);
		  Exit("\n");
		}
	      else if (! Filexists (optarg))
		{
		  char choix;
		  strcpy (tmp, "\n. The file '");
		  strcat (tmp, optarg);
		  strcat (tmp, "' doesn't exist.\n");
		  printf("%s",tmp);
		  printf("\n. Type any key to exit.\n");
		  scanf("%c",&choix);
		  Exit("\n");
		}
	      else
		{
		  strcpy(io->in_tree_file, optarg);
		  io->in_tree = 1;
		  io->fp_in_tree = Openfile(io->in_tree_file,0);
		}
	      Free(tmp);
	      break;
	    }

	  case 'v':case 13:
	    {
	      if ((strcmp (optarg, "e") == 0) ||
		  (strcmp (optarg, "E") == 0) ||
		  (strcmp (optarg, "estimated") == 0) ||
		  (strcmp (optarg, "ESTIMATED") == 0)) {
		io->mod->s_opt->opt_num_param = 1;
		io->mod->s_opt->opt_pinvar    = 1;
		io->mod->pinvar               = 0.2;
		io->mod->invar                = 1;
	      }
	      else if ((atof(optarg) < 0.0) || (atof(optarg) > 1.0))
		{
		  char choix;
		  printf("\n. The proportion of invariable site must be a number between 0.0 and 1.0\n");
		  printf("\n. Type any key to exit.");
		  scanf("%c",&choix);
		  Exit("\n");
		}
	      else
		{
		  io->mod->pinvar = (phydbl)atof(optarg);
		  if (io->mod->pinvar > 0.0+MDBL_MIN)
		    io->mod->invar = 1;
		  else
		    io->mod->invar = 0;
		  io->mod->s_opt->opt_pinvar = 0;
		}
	      break;
	    }
	case 'o':
	  {
	    if(!strcmp(optarg,"tlr"))
	      {
		io->mod->s_opt->opt_topo = 1;
		io->mod->s_opt->opt_bl   = 1;
	      }
	    else if(!strcmp(optarg,"tl"))
	      {
		io->mod->s_opt->opt_topo = 1;
		io->mod->s_opt->opt_bl   = 1;
	      }
	    else if(!strcmp(optarg,"t"))
	      {
		Warn_And_Exit("\n. You can't optimize the topology without adjusting branch length too...\n");
	      }
	    else if(!strcmp(optarg,"lr"))
	      {
		io->mod->s_opt->opt_topo = 0;
		io->mod->s_opt->opt_bl   = 1;
	      }
	    else if(!strcmp(optarg,"l"))
	      {
		io->mod->s_opt->opt_topo = 0;
		io->mod->s_opt->opt_bl   = 1;
	      }
	    else if(!strcmp(optarg,"r"))
	      {
		io->mod->s_opt->opt_topo = 0;
		io->mod->s_opt->opt_bl   = 0;
	      }
	    else if(!strcmp(optarg,"none") || !strcmp(optarg,"n"))
	      {
		io->mod->s_opt->opt_topo = 0;
		io->mod->s_opt->opt_bl   = 0;
	      }
	    else
	      {
		char choix;
		printf ("\n. The optimization parameter must be 'tlr' or 'tl' or 'lr' or 'l' or 'r' or ''.");
		printf("\n. Type any key to exit.\n");
		scanf("%c",&choix);
		Exit("\n");
	      }
	    break;
	  }
	case '?':
	  {
	    char choix;
	    if (isprint (optopt))
	      printf ("\n. Unknown option `-%c'.\n", optopt);
	    else
	      printf ("\n. Unknown option character `\\x%x'.\n", optopt);
	    printf("\n. Type any key to exit.\n");
	    scanf("%c",&choix);
	    break;
	  }
	  
	default:
	  Usage();
      }
}
    
    
    if(io->print_site_lnl)
      {
	strcpy(io->out_lk_file,io->in_seq_file);
	strcat(io->out_lk_file, "_phyml_lk.txt");
	io->fp_out_lk = Openfile(io->out_lk_file,1);
      }

    if(io->print_trace)
      {
	strcpy(io->out_trace_file,io->in_seq_file);
	strcat(io->out_trace_file,"_phyml_trace.txt");
	io->fp_out_trace = Openfile(io->out_trace_file,1);
      }

    if(io->mod->s_opt->random_input_tree)
      {
	strcpy(io->out_best_tree_file,io->in_seq_file);
	strcat(io->out_best_tree_file,"_phyml_best_tree.txt");
      }

    if((io->print_boot_trees) && (io->mod->bootstrap > 0))
      {
	strcpy(io->out_boot_tree_file,io->in_seq_file);
	strcat(io->out_boot_tree_file,"_phyml_boot_trees.txt");
	io->fp_out_boot_tree = Openfile(io->out_boot_tree_file,1);
	
	strcpy(io->out_boot_stats_file,io->in_seq_file);
	strcat(io->out_boot_stats_file,"_phyml_boot_stats.txt");
	io->fp_out_boot_stats = Openfile(io->out_boot_stats_file,1);
      }

#ifndef PHYML
    if((open_ps_file) || (io->m4_model == YES))
      {
	strcpy(io->out_ps_file,io->in_seq_file);
	strcat(io->out_ps_file, "_mc_tree.ps");
	io->fp_out_ps = Openfile(io->out_ps_file,1);
      }
#endif 


    if((use_gamma) && (io->mod->n_catg == 1))
      {
	io->mod->n_catg = 4;
      }
 
    if((io->mod->s_opt->n_rand_starts)           && 
       (io->mod->s_opt->topo_search != SPR_MOVE) && 
       (io->mod->s_opt->random_input_tree))
      {
	Warn_And_Exit("\n. Random starting tree option is only compatible with a full-SPR search.\n"); 
      }

    if ((io->mod->datatype == NT) && (io->mod->whichmodel > 10))
      {
	char choix;
	printf("\n. Err: model incompatible with the data type. Please use JC69, K80, F81, HKY, F84, TN93 or GTR\n");
	printf("\n. Type any key to exit.\n");
	scanf("%c",&choix);
	Warn_And_Exit("\n");
      }
    else if ((io->mod->datatype == AA) && (io->mod->whichmodel < 11))
      {
	char choix;
	printf("\n. Err: model incompatible with the data type. Please use LG, Dayhoff, JTT, MtREV, WAG, DCMut, RtREV, CpREV, VT, Blosum62, MtMam, MtArt, HIVw or HIVb.\n");
	printf("\n. Type any key to exit.\n");
	scanf("%c",&choix);
	Exit("\n");
      }

    if(io->m4_model == YES)
      {
	io->mod->ns *= io->mod->m4mod->n_h;
	io->mod->use_m4mod = 1;
	M4_Make_Complete(io->mod->m4mod->n_h,
			 io->mod->m4mod->n_o,
			 io->mod->m4mod);
      }
    else
      {
	io->mod->s_opt->opt_cov_delta      = 0;
	io->mod->s_opt->opt_cov_alpha      = 0;
	io->mod->s_opt->opt_cov_free_rates = 0;
      }

    if((io->mod->s_opt->opt_cov_free_rates) && (io->mod->s_opt->opt_cov_alpha))
      {
	io->mod->s_opt->opt_cov_free_rates = 0;
	io->mod->m4mod->use_cov_alpha      = 0;
	io->mod->m4mod->use_cov_free       = 1;
      }
    
    return;
}

/*********************************************************/

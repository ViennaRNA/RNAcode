/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "utilities.h"
#include "options.h"
#include "models.h"
#include "free.h"
#include "options.h"
#include "interface.h"
#ifdef MG
#include "multigene.h"
#endif


void Launch_Interface(option *io)
{
  Launch_Interface_Input(io);
  io->ready_to_go = 0;
  do
    {
      switch(io->curr_interface)
	{
	case INTERFACE_DATA_TYPE :
	  {
	    Launch_Interface_Data_Type(io);
	    break;
	  }
	case INTERFACE_MULTIGENE :
	  {
	    Launch_Interface_Multigene(io);
	    break;
	  }
	case INTERFACE_MODEL :
	  {
	    Launch_Interface_Model(io);
	    break;
	  }
	case INTERFACE_TOPO_SEARCH :
	  {
	    Launch_Interface_Topo_Search(io);
	    break;
	  }
	case INTERFACE_BRANCH_SUPPORT :
	  {
	    Launch_Interface_Branch_Support(io);
	    break;
	  }
	default :
	  {
	    printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Exit("");
	    break;
	  }
	}

    }while(!io->ready_to_go);
}

/*********************************************************/

void Clear()
{
#ifdef WIN32
  system("cls");
#elif UNIX
  printf("\033[2J\033[H");
#endif
}

/*********************************************************/

void Launch_Interface_Input(option *io)
{
  char choix;
  int n_trial;

  Clear();
  Print_Banner(stdout);


#ifdef EVOLVE

  char *n_data_sets;

  printf("\n\n");
  printf("\n. Enter the tree file name > "); fflush(NULL);
  Getstring_Stdin(io->in_tree_file);
  io->fp_in_tree = Openfile(io->in_tree_file,0);
  printf("\n");

  printf("\n. Enter the reference sequence file name > "); fflush(NULL);
  Getstring_Stdin(io->in_seq_file);
  io->fp_in_seq = Openfile(io->in_seq_file,0);
  printf("\n");

  printf("\n. Number of data sets > ");
  n_data_sets = (char *)mCalloc(10000,sizeof(char));
  Getstring_Stdin(n_data_sets);
  n_trial = 0;
  while((!atoi(n_data_sets)) || (atoi(n_data_sets) < 0))
    {
      if(++n_trial > 10) Exit("\n. Err : the number of sets must be a positive integer");
      printf("\n. The number of sets must be a positive integer");
      printf("\n. Enter a new value > ");
      Getstring_Stdin(n_data_sets);
    }
  io->n_data_set_asked = atoi(n_data_sets);
  Free(n_data_sets);

#elif OPTIMIZ

  printf("\n. Enter the tree file name > "); fflush(NULL);
  Getstring_Stdin(io->in_tree_file);
  io->fp_in_tree = Openfile(io->in_tree_file,0);
  printf("\n");

  printf("\n. Enter the reference sequence file name > "); fflush(NULL);
  Getstring_Stdin(io->in_seq_file);
  io->fp_in_seq = Openfile(io->in_seq_file,0);
  printf("\n");

#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)

  printf("\n. Enter the sequence file name > "); fflush(NULL);
  Getstring_Stdin(io->in_seq_file);
  io->fp_in_seq = Openfile(io->in_seq_file,0);

#endif


#if defined(PHYML) || defined(MG) || defined(PHYML_INSERT)

  strcpy(io->out_stats_file,io->in_seq_file);
  strcat(io->out_stats_file,"_phyml_stats.txt");
  
  strcpy(io->out_tree_file,io->in_seq_file);
  strcat(io->out_tree_file,"_phyml_tree.txt");

  strcpy(io->out_lk_file,io->in_seq_file);
  strcat(io->out_lk_file,"_phyml_lk.txt");


#endif


#ifdef WIN32
#ifdef EVOLVE
  if(Filexists("evolve_out.txt"));
#elif OPTIMIZ
  if(Filexists("optimiz_out.txt"))
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
  if(Filexists(io->out_stat_file))
#endif
#elif UNIX
#ifdef EVOLVE
  if(Filexists("evolve_out"));
#elif OPTIMIZ
  if(Filexists("optimiz_out"))
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
   if(Filexists(io->out_stats_file))
#endif
#endif
    {
      printf("\n");
#ifdef EVOLVE
      printf("\n. A file 'evolve_out' already exists");
#elif OPTIMIZ
      printf("\n. A file 'optimiz_out' already exists");
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
      printf("\n. A file '%s' already exists",io->out_stats_file);
#endif
      printf("\n. Do you want to Replace it or Append to it ? ");
      n_trial = 0;
      do
	{
	  printf("\n. Please type R or A > ");
	  scanf("%c",&choix);
	  if(choix == '\n') choix = 'r';
	  else getchar();
	  if(++n_trial>10) Exit("\n");
	  Uppercase(&choix);
	}
      while((choix != 'R') && (choix != 'A'));
      if(choix == 'R') io->out_stats_file_open_mode = 1;
      else             io->out_stats_file_open_mode = 2;
    }

  io->fp_out_stats = Openfile(io->out_stats_file,io->out_stats_file_open_mode);

#ifdef WIN32
#ifdef EVOLVE
  if(Filexists("evolve_seq.txt"))
#elif OPTIMIZ
  if(Filexists("optimiz_tree.txt"))
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
  if(Filexists(io->out_tree_file))
#endif
#elif UNIX
#ifdef EVOLVE
  if(Filexists("evolve_seq"))
#elif OPTIMIZ
  if(Filexists("optimiz_tree"))
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
  if(Filexists(io->out_tree_file))
#endif
#endif
    {
      printf("\n");
#ifdef EVOLVE
      printf("\n. A file 'evolve_seq' already exists");
#elif OPTIMIZ
      printf("\n. A file 'optimiz_tree' already exists");
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
      printf("\n. A file '%s' already exists",io->out_tree_file);
#endif
      printf("\n. Do you want to Replace it or Append to it ? ");
      n_trial = 0;
      do
	{
	  printf("\n. Please type R or A > ");
	  scanf("%c",&choix);
	  if(choix == '\n') choix = 'r';
	  else getchar();
	  Uppercase(&choix);
	  if(++n_trial>10) Exit("\n");
	}
      while((choix != 'R') && (choix != 'A'));
      if(choix == 'R') io->out_tree_file_open_mode = 1;
      else             io->out_tree_file_open_mode = 2;
    }
  io->fp_out_tree = Openfile(io->out_tree_file,io->out_tree_file_open_mode);
}

/*********************************************************/

void Launch_Interface_Data_Type(option *io)
{
  char choix;
  char *s;
  char *buff;

  Clear();
  Print_Banner(stdout);
  if(io->config_multigene) Print_Data_Set_Number(io,stdout);

  s    = (char *)mCalloc(100,sizeof(char));
  buff = (char *)mCalloc(100,sizeof(char));


  printf("\n\n");

  printf("                                   ...................                                              \n");
  printf("                                    Menu : Input Data                                               \n");
  printf("                                .........................                                           \n");

  printf("\n\n");

  printf("                [+] "
	 ".................................... Next sub-menu\n");

  printf("                [-] "
	 "................................ Previous sub-menu\n");

  printf("                [Y] "
	 ".............................. Launch the analysis\n");

  printf("\n");

  printf("                [D] "
	 "............................... Data type (DNA/AA) "
	 " %-15s \n",
	 (io->mod->datatype)?("AA"):("DNA"));

  printf("                [I] "
	 "...... Input sequences interleaved (or sequential) "
	 " %-15s \n",
	 (io->interleaved)?("interleaved"):("sequential"));

  strcpy(s,"");
  sprintf(s," (%d sets)",io->n_data_sets);
  strcpy(buff,(io->n_data_sets > 1)?("yes"):("no"));
  buff=strcat(buff,(io->n_data_sets > 1)?(s):("\0"));
  printf("                [M] "
	 "....................... Analyze multiple data sets "
	 " %-15s \n",buff);


  printf("\n\n. Are these settings correct ? "
	 "(type '+', '-', 'Y' or other letter for one to change)  ");

 
  scanf("%c",&choix);
  if(choix != '\n') getchar(); /* \n */
  fflush(NULL);

  Uppercase(&choix);

  switch(choix)
    {
/*     case '\n': */
/*       { */
/* 	io->curr_interface++; */
/* 	break; */
/*       } */
    case 'M' :
      {
	char *c;
	int n_trial;

	printf("\n. How many data sets > ");
	c = (char *)mCalloc(100,sizeof(char));
	Getstring_Stdin(c);
	n_trial = 0;
	while((!atoi(c)) || (atoi(c) < 0))
	  {
	    if(++n_trial > 10) Exit("\n. Err : The number of data sets must be a positive integer");
	    printf("\n. The number of data sets must be a positive integer");
	    printf("\n. Enter a new value > ");
	    Getstring_Stdin(c);
	  }
	io->n_data_sets = atoi(c);

	#ifdef MG
	if(io->n_data_sets > 1) 
	  {
	    io->multigene = 1;
	  }
	#endif

	if((io->mod->bootstrap > 1) && (io->n_data_sets > 1))
	  {
	    printf("\n. Bootstrap option is not allowed with multiple data sets !\n");
	    printf("\n. Type any key to exit.");
	    scanf("%c",&choix);
	    Exit("\n");
	  }

	Free(c);
	break;
      }
    case 'I' :
      {
	if(io->interleaved)
	  io->interleaved = 0;
	else io->interleaved = 1;
	break;
      }
    case 'D' :
      {
	if(io->mod->datatype == NT)
	  {
	    io->mod->datatype   = 1;
	    io->mod->stepsize   = 1;
	    io->mod->ns         = 20;
	    io->mod->whichmodel = LG;
	    strcpy(io->mod->modelname,"LG");
	  }
	else
	  {
	    io->mod->datatype   = 0;
	    io->mod->stepsize   = 1;
	    io->mod->ns         = 4;
	    io->mod->whichmodel = HKY85;
	    strcpy(io->mod->modelname,"HKY85");
	    strcpy(io->nt_or_cd,"nucleotides");
	  }
	break;
      }
    case '-' :
      {
	io->curr_interface = (io->config_multigene)?(INTERFACE_MODEL):(INTERFACE_BRANCH_SUPPORT);
	break;
      }
    case '+' :
      {
	io->curr_interface = (io->multigene)?(INTERFACE_MULTIGENE):(INTERFACE_MODEL);
	break;
      }
    case 'Y' :
      {
	io->ready_to_go = 1;
	break;
      }
    default :
      {
	break;
      }
    }

  Free(s);
  Free(buff);
}
/*********************************************************/

void Launch_Interface_Model(option *io)
{
  char choix;
  char *s;

  s = (char *)mCalloc(100,sizeof(char));

  Clear();
  Print_Banner(stdout);
  if(io->config_multigene) Print_Data_Set_Number(io,stdout);


  printf("\n\n");

  printf("                                ...........................                                      \n");
  printf("                                 Menu : Substitution Model                                       \n");
  printf("                             .................................                                   \n");

  printf("\n\n");

  printf("                [+] "
	 ".................................... Next sub-menu\n");

  printf("                [-] "
	 "................................ Previous sub-menu\n");

  printf("                [Y] "
	 ".............................. Launch the analysis\n");

  printf("\n");

  if (io->mod->datatype == NT)
    {
      if(!strcmp(io->nt_or_cd,"nucleotides"))
	{
	  printf("                [M] "
		 "................. Model of nucleotide substitution "
		 " %-15s \n", io->mod->modelname);

	  if((io->mod->whichmodel == F81)   ||
	     (io->mod->whichmodel == HKY85) ||
	     (io->mod->whichmodel == F84)   ||
	     (io->mod->whichmodel == TN93)  ||
	     (io->mod->whichmodel == GTR)   ||
	     (io->mod->whichmodel == CUSTOM))
	    {
/* 	      printf("                [F] " */
/* 		     ".......... Base frequency estimates (empirical/ML) " */
/* 		     " %-15s \n", */
/* 		     (io->mod->s_opt->opt_state_freq)?("ML"):("empirical")); */
/* 	    } */
/* 	  else if(io->mod->whichmodel == CUSTOM) */
/* 	    { */
	      printf("                [F] "
		     "................. Optimise equilibrium frequencies "
		     " %-15s \n",
		     (io->mod->s_opt->opt_state_freq)?("yes"):("no"));
	    }


	  if(io->mod->whichmodel == CUSTOM)
	    {
	      if(!io->mod->s_opt->opt_state_freq)
		{
		  printf("                [E] "
			 "......... Equilibrium frequencies (empirical/user) "
			 " %-15s \n",
			 (io->mod->s_opt->user_state_freq)?("user defined"):("empirical"));
		}
	      printf("                [K] "
		     "............................. Current custom model "
		     " %-15s \n", io->mod->custom_mod_string);

	      printf("                [O] "
		     "................ Optimise relative rate parameters "
		     " %-15s \n",(io->mod->s_opt->opt_rr)?("yes"):("no"));
	    }
	}
    }
  else
    {
      printf("                [M] "
	     "................ Model of amino-acids substitution "
	     " %-15s \n", io->mod->modelname);

      printf("                [F] "
	     ". Amino acid frequencies (empirical/model defined) "
	     " %-15s \n",
	     (io->mod->s_opt->opt_state_freq)?("empirical"):("model"));
    }


  if ((io->mod->datatype    == NT)   &&
      ((io->mod->whichmodel == K80)  ||
       (io->mod->whichmodel == HKY85)||
       (io->mod->whichmodel == F84)  ||
       (io->mod->whichmodel == TN93)))
    {
      strcpy(s,(io->mod->s_opt->opt_kappa)?("estimated"):("fixed"));
      (io->mod->s_opt->opt_kappa)?((char *)strcat(s,"")):((char *)strcat(s," (ts/tv = "));
/*       (io->mod->s_opt->opt_kappa)?((char *)strcat(s,"")):((char *)sprintf(s+(int)strlen((char *)s),"%3.2f)",io->mod->kappa)); */
      if(io->mod->s_opt->opt_kappa)
	{
	  strcat((char *)s,"");
	}
      else
	{
	  sprintf((char *)(s+(int)strlen(s)),"%3.2f)",io->mod->kappa);
	}

      printf("                [T] "
	     ".................... Ts/tv ratio (fixed/estimated) "
	     " %-15s \n",s);
    }


  (io->mod->s_opt->opt_pinvar)?(strcpy(s,"estimated")):(strcpy(s,"fixed"));
  (io->mod->s_opt->opt_pinvar)?((char *)strcat(s,"")):((char *)strcat(s," (p-invar = "));

  if(io->mod->s_opt->opt_pinvar)
    {
      strcat(s,"");
    }
  else
    {
      sprintf((char *)(s+strlen((char *)s)),"%3.2f)",io->mod->pinvar);
    }

  printf("                [V] "
	 ". Proportion of invariable sites (fixed/estimated)"
	 "  %-15s \n",s);

  printf("                [R] "
	 "....... One category of substitution rate (yes/no) "
	 " %-15s \n",
	 (io->mod->n_catg > 1)?("no"):("yes"));

  if(io->mod->n_catg > 1)
    {
      printf("                [C] "
	     "........... Number of substitution rate categories "
	     " %-15d \n",
	     io->mod->n_catg);
    }


  if(io->mod->n_catg > 1)
    {
      strcpy(s,(io->mod->s_opt->opt_alpha)?("estimated"):("fixed"));
      (io->mod->s_opt->opt_alpha)?(strcat(s, "")):(strcat(s," (alpha = "));
      
      if(io->mod->s_opt->opt_alpha)
	{
	  strcat(s, "");
	}
      else
	{
	  sprintf(s+strlen(s),"%3.2f)",io->mod->alpha);
	}

      printf("                [A] "
	     "... Gamma distribution parameter (fixed/estimated) "
	     " %-15s \n",s);
    }

  printf("\n\n. Are these settings correct ? "
	 "(type '+', '-', 'Y' or other letter for one to change)  ");


  scanf("%c",&choix);
  if(choix != '\n') getchar(); /* \n */

  Uppercase(&choix);

  switch(choix)
    {
/*     case '\n': */
/*       { */
/* 	io->curr_interface++; */
/* 	break; */
/*       } */

    case 'O' :
      {
	io->mod->s_opt->opt_rr = (io->mod->s_opt->opt_rr)?(0):(1);
	break;
      }

    case 'K' :
      {
	int i,j;
	char **rr_param,*rr;
	model *mod;
	int curr_param;
	int n_trial;

	if(io->mod->whichmodel == CUSTOM)
	  {
	    rr_param = (char **)mCalloc(6,sizeof(char *));
	    For(i,6) rr_param[i] = (char *)mCalloc(10,sizeof(char));
	    rr = (char *)mCalloc(50,sizeof(char));

	    mod = io->mod;

	    n_trial = 0;
	    do
	      {
		printf("\n. Enter a new custom model > ");
		Getstring_Stdin(io->mod->custom_mod_string);
		if(strlen(io->mod->custom_mod_string) == 6)
		  {
		    For(i,6)
		      {
			while(!isdigit((int)io->mod->custom_mod_string[i]))
			  {
			    if(++n_trial > 10) Exit("\n. Err : this string is not valid !\n");
			    printf("\n. This string is not valid\n");
			    printf("\n. Enter a new model > ");
			    Getstring_Stdin(io->mod->custom_mod_string);
			  }
		      }
		    if(i == 6) break;
		  }
		else
		  {
		    printf("\n. The string should be of length 6\n");
		    n_trial++;
		  }
	      }while(n_trial < 10);
	    if(n_trial == 10) Exit("");

	    Translate_Custom_Mod_String(io->mod);

	    strcpy(rr_param[0],"A<->C");
	    strcpy(rr_param[1],"A<->G");
	    strcpy(rr_param[2],"A<->T");
	    strcpy(rr_param[3],"C<->G");
	    strcpy(rr_param[4],"C<->T");
	    strcpy(rr_param[5],"G<->T");

	    printf("\n. Set the relative rate values\n");
	    curr_param = 0;
	    For(i,mod->n_diff_rr)
	      {
		sprintf(rr,"\n. [");
		For(j,6)
		  {
		    if(mod->rr_num[j] == i) 
		      {
			sprintf(rr+strlen(rr),"%s = ",rr_param[j]);
		      }
		  }
		sprintf(rr+strlen(rr)-3,"]");
		printf("%s",rr);

		printf("  (current=%.2f) > ",mod->rr_val[i]);
		
		Getstring_Stdin(rr);
		
		if(rr[0] != '\0')
		  {
		    n_trial = 0;
		    while((atof(rr) < .0))
		      {
			if(++n_trial > 10)
			  Exit("\n. Err : the value of this parameter must be a positive number\n");
			printf("\n. The value of this parameter must be a positive number\n");
			printf("\n. Enter a new value > ");
			Getstring_Stdin(rr);
		      }
		    io->mod->rr_val[i] = (phydbl)atof(rr);
		  }
	      }

	    For(i,5) Free(rr_param[i]);
	    Free(rr_param);
	    Free(rr);
	  }
	break;
      }
    case 'E' :
      {
	int i;

	if(io->mod->whichmodel == CUSTOM)
	  {
	    io->mod->s_opt->user_state_freq = (io->mod->s_opt->user_state_freq)?(0):(1);

	    if(io->mod->s_opt->user_state_freq)
	      {
		if(!io->mod->s_opt->opt_state_freq)
		  {
		    char **bases;
		    char *bs;
		    phydbl sum;
		    int n_trial;
		    
		    bases = (char **)mCalloc(4,sizeof(char *));
		    For(i,4) bases[i] = (char *)mCalloc(50,sizeof(char));
		    bs = (char *)mCalloc(100,sizeof(char));
		    
		    strcpy(bases[0],". f(A)> ");
		    strcpy(bases[1],". f(C)> ");
		    strcpy(bases[2],". f(G)> ");
		    strcpy(bases[3],". f(T)> ");
		    
		    printf("\n. Set nucleotide frequencies \n");
		    sum = .0;
		    For(i,4)
		      {
			printf("%s",bases[i]);
			Getstring_Stdin(bs);
			n_trial = 0;
			
			while((atof(bs) < .0001) || (bs[0] == '\0'))
			  {
			    if(++n_trial > 10)
			      Exit("\n. Err : the value of this parameter must be a positive number\n");
			    printf("\n. The value of this parameter must be a positive number\n");
			    printf("\n. Enter a new value > ");
			    Getstring_Stdin(bs);
			  }
			io->mod->user_b_freq[i] = (phydbl)atof(bs);
			sum += io->mod->user_b_freq[i];
		      }
		
		    For(i,4) io->mod->user_b_freq[i] /= sum;

		    if(sum > 1.0 || sum < 1.0)
		      {
			printf("\n. The nucleotide frequencies have to be normalised in order to sum to 1.0.\n");
			printf("\n. The frequencies are now : f(A)=%f, f(C)=%f, f(G)=%f, f(T)=%f.\n",
			       io->mod->user_b_freq[0],
			       io->mod->user_b_freq[1],
			       io->mod->user_b_freq[2],
			       io->mod->user_b_freq[3]);			  
			printf("\n. Enter any key to continue.\n");
			scanf("%c",bs);
		      }

		    For(i,4) Free(bases[i]);
		    Free(bases);
		    Free(bs);
		  }
		else
		  {
		    Warn_And_Exit("\n. 'E' is not a valid option with these model settings.\n");
		  }
	      }	    
	  }
	break;
      }
    case 'A' :
      {
	char answer;
	int n_trial;

	switch(io->mod->s_opt->opt_alpha)
	  {
	  case 0 :
	    {
	      printf("\n. Optimise alpha ? [Y/n] ");
	      scanf("%c",&answer);
	      if(answer == '\n') answer = 'Y';
	      else getchar();
	      break;
	    }
	  case 1 :
	    {
	      printf("\n. Optimise alpha ? [N/y] ");
	      scanf("%c",&answer);
	      if(answer == '\n') answer = 'N';
	      else getchar();
	      break;
	    }
	  default : Exit("\n");
	  }

	n_trial = 0;
	while((answer != 'Y') && (answer != 'y') &&
	      (answer != 'N') && (answer != 'n'))
	  {
	    if(++n_trial > 10) Exit("\n. Err : wrong answers !");
	    printf("\n. Optimise alpha ? [N/y] ");
	    scanf("%c",&answer);
	    if(answer == '\n') answer = 'N';
	    else getchar();
	  }

	switch(answer)
	  {
	  case 'Y' : case 'y' :
	    {
	      io->mod->s_opt->opt_alpha = 1;
	      io->mod->s_opt->opt_num_param = 1;
	      break;
	    }
	  case 'N' : case 'n' :
	    {
	      char *a;
	      a = (char *)mCalloc(100,sizeof(char));
	      io->mod->alpha = 10.0;
	      io->mod->s_opt->opt_alpha = 0;
	      printf("\n. Enter your value of alpha > ");
	      Getstring_Stdin(a);
	      n_trial = 0;
	      while((!atof(a)) || (atof(a) < .0))
		{
		  if(++n_trial > 10) Exit("\n. Err : alpha must be a positive number\n");
		  printf("\n. Alpha must be a positive number\n");
		  printf("\n. Enter a new value > ");
		  Getstring_Stdin(a);
		}
	      io->mod->alpha = (phydbl)atof(a);
	      Free(a);
	      io->mod->s_opt->opt_alpha  = 0;
	      break;
	    }
	  }
	break;
      }

    case 'C' :
      {
	char *c;
	int n_trial;

	printf("\n. Enter your number of categories > ");
	c = (char *)mCalloc(100,sizeof(char));
	Getstring_Stdin(c);
	n_trial = 0;
	while((!atoi(c)) || (atoi(c) < 0))
	  {
	    if(++n_trial > 10) Exit("\n. Err : the number of categories must be a positive integer\n");
	    printf("\n. The number of categories must be a positive integer\n");
	    printf("\n. Enter a new value > ");
	    Getstring_Stdin(c);
	  }
	io->mod->n_catg = atoi(c);
	Free(c);
	break;
      }

    case 'R' :
      {
	(io->mod->n_catg == 1)?(io->mod->n_catg = 4):(io->mod->n_catg = 1);
	break;
      }

    case 'V' :
      {
	char answer;
	int n_trial;

	switch(io->mod->s_opt->opt_pinvar)
	  {
	  case 0 :
	    {
	      printf("\n. Optimise p-invar ? [Y/n] ");
	      scanf("%c", &answer);
	      if(answer == '\n') answer = 'Y';
	      else getchar();
	      break;
	    }
	  case 1 :
	    {
	      printf("\n. Optimise p-invar ? [N/y] ");
	      scanf("%c", &answer);
	      if(answer == '\n') answer = 'N';
	      else getchar();
	      break;
	    }
	  default : Exit("\n");
	  }

	n_trial = 0;
	while((answer != 'Y') && (answer != 'y') &&
	      (answer != 'N') && (answer != 'n'))
	  {
	    if(++n_trial > 10) Exit("\n. Err : wrong answers !");
	    printf("\n. Optimise p-invar ? [N/y] ");
	    scanf("%c", &answer);
	    if(answer == '\n') answer = 'N';
	    else getchar();
	  }

	switch(answer)
	  {
	  case 'Y' : case 'y' :
	    {
	      io->mod->s_opt->opt_num_param = 1;
	      io->mod->s_opt->opt_pinvar = 1;
	      io->mod->pinvar = 0.2;
	      io->mod->invar  = 1;
	      break;
	    }
	  case 'N' : case 'n' :
	    {
	      char *p;
	      p = (char *)mCalloc(100,sizeof(char));
	      printf("\n. Enter your value of p-invar > ");
	      Getstring_Stdin(p);
	      n_trial = 0;
	      while((atof(p) < 0.0) || (atof(p) > 1.0))
		{
		  if(++n_trial > 10)
		    Exit("\n. Err : the proportion of invariable sites must be a positive number between 0.0 and 1.0");
		  printf("\n. The proportion must be a positive number between 0.0 and 1.0\n");
		  printf("\n. Enter a new value > ");
		  Getstring_Stdin(p);
		}
	      io->mod->pinvar = (phydbl)atof(p);

	      if(io->mod->pinvar > 0.0+MDBL_MIN) io->mod->invar = 1;
	      else                             io->mod->invar = 0;

	      Free(p);

	      io->mod->s_opt->opt_pinvar = 0;
	      break;
	    }
	  }
	break;
      }

    case 'T' :
      {
	char answer;
	int n_trial;

	if((io->mod->datatype   == AA)  ||
	   (io->mod->whichmodel == JC69)||
	   (io->mod->whichmodel == F81) ||
	   (io->mod->whichmodel == GTR) ||
	   (io->mod->whichmodel == CUSTOM))
	  {
	    printf("\n. 'K' is not a valid choice for this model\n");
	    printf("\n. Type any key to exit.\n");
	    scanf("%c",&choix);
	    Exit("\n");
	  }

	switch(io->mod->s_opt->opt_kappa)
	  {
	  case 0 :
	    {
	      printf("\n. Optimise ts/tv ratio ? [Y/n] ");
	      scanf("%c", &answer);
	      if(answer == '\n') answer = 'Y';
	      else getchar();
	      break;
	    }
	  case 1 :
	    {
	      printf("\n. Optimise ts/tv ratio ? [N/y] ");
	      scanf("%c", &answer);
	      if(answer == '\n') answer = 'N';
	      else getchar();
	      break;
	    }
	  default : Exit("\n");
	  }

	n_trial = 0;
	while((answer != 'Y') && (answer != 'y') &&
	      (answer != 'N') && (answer != 'n'))
	  {
	    if(++n_trial > 10) Exit("\n. Err : wrong answers !");
	    printf("\n. Optimise ts/tv ratio ? [N/y] ");
	    scanf("%c", &answer);
	    if(answer == '\n') answer = 'N';
	    else getchar();
	  }

	switch(answer)
	  {
	  case 'Y' : case 'y' :
	    {
	      io->mod->kappa = 4.0;
	      io->mod->s_opt->opt_num_param = 1;
	      io->mod->s_opt->opt_kappa = 1;
	      io->mod->s_opt->opt_kappa = 1;
	      if(io->mod->whichmodel == TN93)
		io->mod->s_opt->opt_lambda = 1;
	      break;
	    }
	  case 'N' : case 'n' :
	    {
	      char *t;
	      t = (char *)mCalloc(100,sizeof(char));
	      io->mod->s_opt->opt_kappa = 0;
	      printf("\n. Enter your value of the ts/tv ratio > ");
	      Getstring_Stdin(t);
	      n_trial = 0;
	      while((!atof(t)) || (atof(t) < .0))
		{
		  if(++n_trial > 10) Exit("\n. Err : the ts/tv ratio must be a positive number\n");
		  printf("\n. The ratio must be a positive number");
		  printf("\n. Enter a new value > ");
		  Getstring_Stdin(t);
		}
	      io->mod->kappa = (phydbl)atof(t);
	      io->mod->s_opt->opt_kappa  = 0;
	      io->mod->s_opt->opt_lambda = 0;
	      Free(t);
	      break;
	    }
	  }
	break;
      }

    case 'F' :
      {
	if((io->mod->whichmodel == JC69) ||
	   (io->mod->whichmodel == K80)) 
	  {
	    Warn_And_Exit("\n. 'F' is not a valid choice with these model settings.\n");
	  }
	io->mod->s_opt->opt_state_freq = (io->mod->s_opt->opt_state_freq)?(0):(1);
	break;
      }

    case 'M' :
      {
	if(io->mod->datatype == NT)
	  {
	    if(!strcmp(io->nt_or_cd,"nucleotides"))
	      {
		if(io->mod->whichmodel == JC69)
		  {
		    io->mod->whichmodel = K80;
		    strcpy(io->mod->modelname,"K80");
		  }
		else if(io->mod->whichmodel == K80)
		  {
		    io->mod->whichmodel = F81;
		    strcpy(io->mod->modelname,"F81");
		    io->mod->s_opt->opt_kappa = 0;
		  }
		else if(io->mod->whichmodel == F81)
		  {
		    io->mod->whichmodel = HKY85;
		    strcpy(io->mod->modelname,"HKY85");
		  }
		else if(io->mod->whichmodel == HKY85)
		  {
		    io->mod->whichmodel = F84;
		    strcpy(io->mod->modelname,"F84");
		  }
		else if(io->mod->whichmodel == F84)
		  {
		    io->mod->whichmodel = TN93;
		    strcpy(io->mod->modelname,"TN93");
		    if(io->mod->s_opt->opt_kappa) io->mod->s_opt->opt_lambda = 1;
		  }
		else if(io->mod->whichmodel == TN93)
		  {
		    io->mod->whichmodel = GTR;
		    strcpy(io->mod->modelname,"GTR");
		    io->mod->s_opt->opt_kappa = 0;
		  }
		else if(io->mod->whichmodel == GTR)
		  {
		    io->mod->whichmodel = CUSTOM;
		    strcpy(io->mod->modelname,"custom");
		    io->mod->s_opt->opt_kappa = 0;
		  }

		else if(io->mod->whichmodel == CUSTOM)
		  {
		    io->mod->whichmodel = JC69;
		    strcpy(io->mod->modelname,"JC69");
		    io->mod->s_opt->opt_kappa = 0;
		  }
	      }
	  }
	else
	  {
	    if(io->mod->whichmodel == LG)
	      {
		io->mod->whichmodel = WAG;
		strcpy(io->mod->modelname,"WAG");
	      }
	    else if(io->mod->whichmodel == WAG)
	      {
		io->mod->whichmodel = DAYHOFF;
		strcpy(io->mod->modelname,"Dayhoff");
	      }
	    else if(io->mod->whichmodel == DAYHOFF)
	      {
		io->mod->whichmodel = JTT;
		strcpy(io->mod->modelname,"JTT");
	      }
	    else if(io->mod->whichmodel == JTT)
	      {
		io->mod->whichmodel = BLOSUM62;
		strcpy(io->mod->modelname,"Blossum62");
	      }
	    else if(io->mod->whichmodel == BLOSUM62)
	      {
		io->mod->whichmodel = MTREV;
		strcpy(io->mod->modelname,"MtREV");
	      }
	    else if(io->mod->whichmodel == MTREV)
	      {
		io->mod->whichmodel = RTREV;
		strcpy(io->mod->modelname,"RtREV");
	      }
	    else if(io->mod->whichmodel == RTREV)
	      {
		io->mod->whichmodel = CPREV;
		strcpy(io->mod->modelname,"CpREV");
	      }
	    else if(io->mod->whichmodel == CPREV)
	      {
		io->mod->whichmodel = DCMUT;
		strcpy(io->mod->modelname,"DCMut");
	      }
	    else if(io->mod->whichmodel == DCMUT)
	      {
		io->mod->whichmodel = VT;
		strcpy(io->mod->modelname,"VT");
	      }
	    else if(io->mod->whichmodel == VT)
	      {
		io->mod->whichmodel = MTMAM;
		strcpy(io->mod->modelname,"MtMam");
	      }
	    else if(io->mod->whichmodel == MTMAM)
	      {
		io->mod->whichmodel = CUSTOMAA;
		strcpy(io->mod->modelname,"Read from file");
	      }
	    else if(io->mod->whichmodel == CUSTOMAA)
	      {
		io->mod->whichmodel = LG;
		strcpy(io->mod->modelname,"LG");
	      }
	  }
	break;
      }
    case '-' :
      {
	io->curr_interface = INTERFACE_DATA_TYPE;
	break;
      }
    case '+' :
      {
	io->curr_interface = (io->config_multigene)?(INTERFACE_DATA_TYPE):(INTERFACE_TOPO_SEARCH);
	break;
      }
    case 'Y' :
      {
 	io->ready_to_go = 1;
	break;
      }
    default :
      {
	break;
      }
    }

  if(io->mod->s_opt->opt_alpha  ||
     io->mod->s_opt->opt_kappa  ||
     io->mod->s_opt->opt_lambda ||
     io->mod->s_opt->opt_pinvar ||
     io->mod->s_opt->opt_rr) io->mod->s_opt->opt_num_param = 1;
  else                       io->mod->s_opt->opt_num_param = 0;





  Free(s);
}

/*********************************************************/

void Launch_Interface_Topo_Search(option *io)
{
  char choix;

  Clear();
  Print_Banner(stdout);

  printf("\n\n");


  printf("                                   .......................                                     \n");
  printf("                                    Menu : Tree Searching                                        \n");
  printf("                                .............................                                  \n");


  printf("\n\n");

  printf("                [+] "
	 ".................................... Next sub-menu\n");

  printf("                [-] "
	 "................................ Previous sub-menu\n");

  printf("                [Y] "
	 ".............................. Launch the analysis\n");

  printf("\n");

  printf("                [O] "
	 "........................... Optimise tree topology "
	 " %-15s \n",
	 (io->mod->s_opt->opt_topo)?("yes"):("no"));


  if(io->mod->s_opt->opt_topo)
    {
      char *s;
      
      s = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
      if(io->mod->s_opt->topo_search == NNI_MOVE)
	{
	  if(!io->mod->s_opt->spr_step_after_nnis)
	    {
	      strcpy(s,"NNI moves (default)\0");
	    }
	  else
	    {
	      strcpy(s,"NNI moves plus one SPR step\0");
	    }
	}
      else if(io->mod->s_opt->topo_search == SPR_MOVE)
	strcpy(s,"SPR moves\0");

      printf("                [S] "
	     ".................. Tree topology search operations "
	     " %-15s \n",s);

      Free(s);

      if(io->mod->s_opt->topo_search == SPR_MOVE)
	{
	  printf("                [R] "
		 "......................... Use random starting tree "
		 " %-15s \n",
		 (io->mod->s_opt->random_input_tree)?("yes"):("no"));

	  if(io->mod->s_opt->random_input_tree)
	    {
	      printf("                [N] "
		     ".................. Number of random starting trees "
		     " %-15d \n",io->mod->s_opt->n_rand_starts);	      
	    }
	}
    }
  else
    {
      printf("                [L] "
	     ".......................... Optimise branch lengths "
	     " %-15s \n",
	     (io->mod->s_opt->opt_bl)?("yes"):("no"));
    }
  
  if(!io->mod->s_opt->random_input_tree)
    {
      printf("                [U] "
	     "..................... Input tree (BioNJ/user tree) "
	     " %-15s \n",
	     (!io->in_tree)?("BioNJ"):("user tree"));
    }

  printf("\n\n. Are these settings correct ? "
	 "(type '+', '-', 'Y' or other letter for one to change)  ");


  scanf("%c",&choix);
  if(choix != '\n') getchar(); /* \n */

  Uppercase(&choix);

  switch(choix)
    {
    case '-' :
      {
	io->curr_interface = INTERFACE_MODEL;
	break;
      }
    case '+' :
      {
	io->curr_interface = INTERFACE_BRANCH_SUPPORT;
	break;
      }
    case 'Y' :
      {
 	io->ready_to_go = 1;
	break;
      }

    case 'A' :
      {
	io->mod->s_opt->spr_step_after_nnis = 
	(io->mod->s_opt->spr_step_after_nnis)?(0):(1);
	break;
      }
    case 'U' :
      {
	if(!io->in_tree)
	  {
	    io->in_tree = 1;
	    printf("\n. Enter the name of the tree file > ");
	    Getstring_Stdin(io->in_tree_file);
	    io->fp_in_tree = Openfile(io->in_tree_file,0);
	  }
	else io->in_tree = 0;
	break;
      }

    case 'N' : 
      {	
	char *n;
	int n_trial;

	printf("\n. Enter your number of starting trees > ");
	n = (char *)mCalloc(100,sizeof(char));
	Getstring_Stdin(n);
	n_trial = 0;
	while((!atoi(n)) || (atoi(n) < 0))
	  {
	    if(++n_trial > 10) Exit("\n. Err : the number of starting trees must be a positive integer\n");
	    printf("\n. The number of starting trees must be a positive integer\n");
	    printf("\n. Enter a new value > ");
	    Getstring_Stdin(n);
	  }
	io->mod->s_opt->n_rand_starts = atoi(n);
	io->n_trees = 1;
	Free(n);
	break;
      }
    case 'O' :
      {
	io->mod->s_opt->opt_topo = (io->mod->s_opt->opt_topo)?(0):(1);
	break;
      }

    case 'L' :
      {
	if(!io->mod->s_opt->opt_topo)
	  {
	    io->mod->s_opt->opt_bl = (io->mod->s_opt->opt_bl)?(0):(1);
	  }
	break;
      }

    case 'S' :
      {
	if(io->mod->s_opt->topo_search == NNI_MOVE)
	  {
	    if(!io->mod->s_opt->spr_step_after_nnis)
	      io->mod->s_opt->spr_step_after_nnis = 1;
	    else 
	      {
		io->mod->s_opt->topo_search         = SPR_MOVE;
		io->mod->s_opt->n_rand_starts       = 1;
		io->mod->s_opt->random_input_tree   = 0;
		io->mod->s_opt->greedy              = 0;
		io->mod->s_opt->spr_step_after_nnis = 0;
	      }
	  }
	else if(io->mod->s_opt->topo_search == SPR_MOVE)
	  {
	    io->mod->s_opt->n_rand_starts       = 1;
	    io->mod->s_opt->random_input_tree   = 0;
	    io->mod->s_opt->topo_search         = NNI_MOVE;
	    io->mod->s_opt->greedy              = 0;
	    io->mod->s_opt->spr_step_after_nnis = 0;
	  }
	break;
      }
    case 'R' :
      {
	io->mod->s_opt->random_input_tree = (io->mod->s_opt->random_input_tree)?(0):(1);

	if(io->mod->s_opt->random_input_tree)
	  {
	    if(io->fp_in_tree) fclose(io->fp_in_tree);
	    io->in_tree                   = 0;
	    io->n_trees                   = 1;
	    io->mod->s_opt->n_rand_starts = 5;

	    strcpy(io->out_best_tree_file,io->in_seq_file);
	    strcat(io->out_best_tree_file,"_phyml_best_tree.txt");
	  }
	break;
      }
    default :
      {
	break;
      }
    }
}

/*********************************************************/

void Launch_Interface_Branch_Support(option *io)
{
  char choix;
  char *s;

  s = (char *)mCalloc(100,sizeof(char));

  Clear();
  Print_Banner(stdout);

  printf("\n\n");

  printf("                                  .......................                                          \n");
  printf("                                   Menu : Branch Support                                         \n");
  printf("                               .............................                                       \n");

  printf("\n\n");

  printf("                [+] "
	 ".................................... Next sub-menu\n");

  printf("                [-] "
	 "................................ Previous sub-menu\n");

  printf("                [Y] "
	 ".............................. Launch the analysis\n");

  printf("\n");


  strcpy(s,(io->mod->bootstrap > 0)?("yes"):("no"));
  if(io->mod->bootstrap > 0) sprintf(s+strlen(s)," (%d replicate%s)",
					io->mod->bootstrap,
					(io->mod->bootstrap>1)?("s"):(""));
  
  /*   printf("                [+] " */
  printf("                [B] "
	 "................ Non parametric bootstrap analysis "
	 " %-15s \n",s);

  if(io->ratio_test == 0)
    {
      strcpy(s,"no");
    }
  else if(io->ratio_test == 1)
    {
      strcpy(s,"yes / aLRT statistics");
    }
  else if(io->ratio_test == 2)
    {
      strcpy(s,"yes / Chi2-based supports");
    }
  else if(io->ratio_test == 3)
    {
      strcpy(s,"yes / min of SH-like & Chi2-based supports");
    }
  else
    {
      strcpy(s,"yes / SH-like supports");
    }
  
  
  /*   printf("                [+] " */
  printf("                [A] "
	 "................ Approximate likelihood ratio test "
	 " %-15s \n",s);      

  printf("\n. Are these settings correct ? "
	 "(type '+', '-', 'Y' or other letter for one to change)  ");


  scanf("%c",&choix);
  if(choix != '\n') getchar(); /* \n */

  Uppercase(&choix);
  
  Free(s);

  switch(choix)
    {
/*     case '\n': */
/*       { */
/* 	io->curr_interface++; */
/* 	break; */
/*       } */

    case 'B' :
      {
	if(io->mod->bootstrap > 0) io->mod->bootstrap = 0;
	else
	  {
	    char *r;
	    char answer;
	    int n_trial;

	    io->ratio_test = 0;

	    if(io->n_data_sets > 1)
	      {
		printf("\n. Bootstrap option is not allowed with multiple data sets.\n");
		printf("\n. Type any key to exit.\n");
		scanf("%c",&choix);
		Exit("\n");
	      }


	    printf("\n. Number of replicates > ");
	    r = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
	    Getstring_Stdin(r);
	    n_trial = 0;
	    while((!atoi(r)) || (atoi(r) < 0))
	      {
		if(++n_trial > 10) Exit("\n. Err : the number of replicates must be a positive integer\n");
		printf("\n. The number of replicates must be a positive integer");
		printf("\n. Enter a new value > ");
		Getstring_Stdin(r);
	      }
	    io->mod->bootstrap = atoi(r);

	    printf("\n. Print bootstrap trees (and statistics) ? (%s) > ",
		   (io->print_boot_trees)?("Y/n"):("y/N"));

	    scanf("%c",&answer);
	    if(answer == '\n') answer = (io->print_boot_trees)?('Y'):('N');
	    else getchar();

	    switch(answer)
	      {
	      case 'Y' : case 'y' :
		{
		  io->print_boot_trees = 1;

		  strcpy(io->out_boot_tree_file,io->in_seq_file);
		  strcat(io->out_boot_tree_file,"_phyml_boot_trees.txt");
		  io->fp_out_boot_tree = Openfile(io->out_boot_tree_file,1);

		  strcpy(io->out_boot_stats_file,io->in_seq_file);
		  strcat(io->out_boot_stats_file,"_phyml_boot_stats.txt");
		  io->fp_out_boot_stats = Openfile(io->out_boot_stats_file,1);

		  break;
		}
	      case 'N' : case 'n' :
		{
		  io->print_boot_trees  = 0;
		  io->fp_out_boot_tree  = NULL;
		  io->fp_out_boot_stats = NULL;
		  break;
		}
	      }
	    Free(r);
	  }
	break;
      }
    case 'A' :
      {
	io->mod->bootstrap = 0;

	switch(io->ratio_test)
	  {
	  case 0 :
	    {
	      io->ratio_test = 1;
	      break;
	    }
	  case 1 :
	    {
	      io->ratio_test = 2;
	      break;
	    }
	  case 2 :
	    {
	      io->ratio_test = 3;
	      break;
	    }
	  case 3 :
	    {
	      io->ratio_test = 4;
	      break;
	    }
	  case 4 :
	    {
	      io->ratio_test = 0;
	      break;
	    }
	  }
	break;
      }

    case '-' :
      {
	io->curr_interface = INTERFACE_TOPO_SEARCH;
	break;
      }
    case '+' :
      {
	io->curr_interface = INTERFACE_DATA_TYPE;
	break;
      }
    case 'Y' :
      {
 	io->ready_to_go = 1;
	break;
      }
    default : break;
    }
}

/*********************************************************/

void Launch_Interface_Multigene(option *io)
{

#ifdef MG

  if((io->n_data_sets > 1) && (io->multigene))
    {
      int set, n_trial;
      char choix;
      
      io->n_gt     = io->n_data_sets;
      io->st       = (superarbre *)Mg_Make_Superarbre_Light(io);
      io->st->n_gt = io->n_data_sets;

      For(set,io->n_gt)
	{
	  io->st->optionlist[set] = Make_Input();
	  Set_Defaults_Input(io->st->optionlist[set]);
	  Set_Defaults_Model(io->st->optionlist[set]->mod);
	  Set_Defaults_Optimiz(io->st->optionlist[set]->mod->s_opt);
	  io->st->optionlist[set]->curr_gt = set;
	  printf("\n. Enter the sequence file name [data set %2d] > ",set+1); fflush(NULL);
	  Getstring_Stdin(io->st->optionlist[set]->in_seq_file);
	  io->st->optionlist[set]->fp_in_seq = Openfile(io->st->optionlist[set]->in_seq_file,0);
	}
      
      
      printf("\n. Do you want to use your own initial tree ? [N/y] ");
      scanf("%c", &choix);
      if(choix != '\n') getchar();
      
      n_trial = 0;
      while((choix != 'Y') && (choix != 'y') &&
	    (choix != 'N') && (choix != 'n'))  
	{
	  if(++n_trial > 10) Exit("\n. Err : wrong answers !");
	  printf("\n. Do you want to use your own initial tree ? [N/y] ? ");
	  scanf("%c", &choix);
	  if(choix == '\n') choix = 'N';
	  else getchar();
	}
      
      switch(choix)
	{
	case '\n' : break;
	case 'Y' : case 'y' : 
	  {
	    io->in_tree = 1;
	    printf("\n. Enter the name of the tree file > ");
	    Getstring_Stdin(io->in_tree_file);
	    io->fp_in_tree = Openfile(io->in_tree_file,0);
	    break;
	  }
	case 'N' : case 'n' : 
	  {
	    io->in_tree = 0;
	    break;
	  }
	default : break;
	}
      
      io->curr_gt = 0;
      
      do
	{
	  io->st->optionlist[io->curr_gt]->config_multigene = 1;
	  do
	    {
	      switch(io->st->optionlist[io->curr_gt]->curr_interface)
		{
		case INTERFACE_DATA_TYPE :
		  {
		    Launch_Interface_Data_Type(io->st->optionlist[io->curr_gt]);
		    break;
		  }
		case INTERFACE_MODEL : 
		  {
		    Launch_Interface_Model(io->st->optionlist[io->curr_gt]);
		    break;
		  }
		default : 
		  {
		    printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		    Exit("");
		    break;		  
		  }
		}
	    }while(!io->st->optionlist[io->curr_gt]->ready_to_go);
	  io->curr_gt++;
	}while(io->curr_gt < io->n_gt);      
    }
  io->ready_to_go = 1;
  Mg_Fill_Model_Partitions_Table(io->st);
#endif
}

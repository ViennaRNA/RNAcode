/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <stdlib.h>
#include "cl.h"
#include "utilities.h"
#include "options.h"
#include "models.h"
#include "free.h"
#include "interface.h"
#ifdef MG
#include "mg.h"
#endif

/* int  T_MAX_FILE; */
/* phydbl MDBL_MIN; */
/* phydbl UNLIKELY; */

/*********************************************************/

void Usage()
{

  char *BOLD=(char *)mCalloc(10,sizeof(char));
  char *FLAT=(char *)mCalloc(10,sizeof(char));
  char *LINE=(char *)mCalloc(10,sizeof(char));
  char *cha;


  cha =getenv("OS");

  if(cha!=NULL) 
    {
      strcpy(BOLD, "");
      strcpy(FLAT, "");
      strcpy(LINE, "");
    } 
  else 
    {
      strcpy(BOLD, "\033[00;01m");
      strcpy(FLAT, "\033[00;00m");
      strcpy(LINE, "\033[00;04m");
    }

      printf("%sNAME\n"
   "%s\t- PhyML %s - \n\n"
	 "%s\t\''A simple, fast, and accurate algorithm to estimate\n"
	 "%s\tlarge phylogenies by maximum likelihood.\''\n\n"
	 "%s\tStephane Guindon and Olivier Gascuel,\n"
	 "%s\tSystematic Biology 52(5):696-704, 2003.\n\n"
	 "%s\tPlease cite this paper if you use this software in your publications.\n",BOLD,FLAT,VVERSION,FLAT,FLAT,FLAT,FLAT,FLAT);
  
  

  printf("%s\nSYNOPSIS:\n\n"
	 "%s\tphyml %s[command args]\n",BOLD,BOLD,BOLD);
  
  printf("%s\nCommand options:\n%s",BOLD,FLAT);

  printf("\n\t%s-i (or --input) %sseq_file_name%s\n",BOLD,LINE,FLAT); 
  printf("\t\t%sseq_file_name%s is the name of the nucleotide or amino-acid sequence file in PHYLIP format.\n",LINE,FLAT);
  
  printf("\n");

  printf("%s\n\t-d (or --datatype) ""%sdata_type%s\n",BOLD,LINE,FLAT);
  printf("%s\t\tdata_type%s is 'nt' for nucleotide (default) and 'aa' for amino-acid sequences.\n",LINE,FLAT);

  printf("\n");

  printf("%s\n\t-q (or --sequential)\n",BOLD);
  printf("%s\t\tChanges interleaved format (default) to sequential format.\n",FLAT);

  printf("\n");

  printf("%s\n\t-n (or --multiple) ""%snb_data_sets%s\n",BOLD,LINE,FLAT);
  printf("%s\t\tnb_data_sets%s is an integer corresponding to the number of data sets to analyse.\n",LINE,FLAT);

  printf("\n");

  printf("%s\n\t-b (or --bootstrap) %sint%s\n",BOLD,LINE,FLAT);
  printf("\t\t%sint%s >  0 : %sint%s is the number of bootstrap replicates.\n",LINE,FLAT,LINE,FLAT);
  printf("\t\t%sint%s =  0 : neither approximate likelihood ratio test nor bootstrap values are computed.\n",LINE,FLAT);
  printf("\t\t%sint%s = -1 : approximate likelihood ratio test returning aLRT statistics.\n",LINE,FLAT);
  printf("\t\t%sint%s = -2 : approximate likelihood ratio test returning Chi2-based parametric branch supports.\n",LINE,FLAT);
  printf("\t\t%sint%s = -3 : minimum of Chi2-based parametric and SH-like branch supports.\n",LINE,FLAT);
  printf("\t\t%sint%s = -4 : SH-like branch supports alone.\n",LINE,FLAT);

  printf("\n");
  
  printf("%s\n\t-m (or --model) %smodel%s\n",BOLD,LINE,FLAT);
  printf("\t\tmodel%s : substitution model name.\n",FLAT);
  printf("\t\t%s- %sNucleotide%s-based models : %sHKY85%s (default) | %sJC69%s | %sK80%s | %sF81%s | %sF84%s | %sTN93%s | %sGTR%s | %scustom (*)%s\n",
	 FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT);
  printf("\t\t(*) : for the custom option, a string of six digits identifies the model. For instance, 000000\n");
  printf("\t\t corresponds to F81 (or JC69 provided the distribution of nucleotide frequencies is uniform).\n");
  printf("\t\t 012345 corresponds to GTR. This option can be used for encoding any model that is a nested within GTR.\n");
  printf("\n");
  printf("\t\t%s- %sAmino-acid%s based models : %sWAG%s (default) | %sJTT%s | %sMtREV%s | %sDayhoff%s | %sDCMut%s | %sRtREV%s | %sCpREV%s | %sVT%s\n",	 
	 FLAT,LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT);
  printf("\t\t %sBlosum62%s | %sMtMam%s | %sMtArt%s | %sHIVw%s |  %sHIVb%s | %scustom%s\n",
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT);

  printf("\n");

  printf("%s\n\t-f %se%s, %sd%s, or %s\"fA fC fG fT\"%s\n",BOLD,LINE,BOLD,LINE,BOLD,LINE,FLAT);
  printf("\t\t%se%s : the character frequencies are determined as follows : \n",LINE,FLAT);
  printf("%s\t\t- %sNucleotide%s sequences: the equilibrium base frequencies are estimated using maximum likelihood \n",FLAT,LINE,FLAT);
  printf("%s\t\t- %sAmino-acid%s sequences: the equilibrium amino-acid frequencies are estimated by counting the\n"
"\t\t occurence of the different amino-acids in the data.\n",FLAT,LINE,FLAT);
  printf("\n");
  printf("\t\t%sd%s : the character frequencies are determined as follows : \n",LINE,FLAT);
  printf("%s\t\t- %sNucleotide%s sequences: the equilibrium base frequencies are estimated by counting the occurence\n"
	 "\t\t of the different bases in the alignment.\n",FLAT,LINE,FLAT);
  printf("%s\t\t- %sAmino-acid%s sequences: the equilibrium amino-acid frequencies are estimated using the frequencies\n"
"\t\t defined by the substitution model.\n",FLAT,LINE,FLAT);
  printf("\n");
  printf("\t\t%s\"fA fC fG fT\"%s : only valid for nucleotide-based models. fA, fC, fG and fT are floating numbers that \n",LINE,FLAT);
  printf("\t\t correspond to the frequencies of A, C, G and T respectively.\n");

  printf("\n");

  printf("%s\n\t-t (or --ts/tv) %sts/tv_ratio%s\n",BOLD,LINE,FLAT);
  printf("\t\tts/tv_ratio%s : transition/transversion ratio. DNA sequences only.\n",FLAT);
  printf("\t\tCan be a fixed positive value (ex:4.0) or %se%s to get the maximum likelihood estimate.\n",LINE,FLAT);

  printf("\n");

  printf("%s\n\t-v (or --pinv) %sprop_invar%s\n",BOLD,LINE,FLAT);
  printf("\t\tprop_invar%s : proportion of invariable sites.\n",FLAT);
  printf("\t\tCan be a fixed value in the [0,1] range or %se%s to get the maximum likelihood estimate.\n",LINE,FLAT);
  
  printf("\n");

  printf("%s\n\t-c (or --nclasses) %snb_subst_cat%s\n",BOLD,LINE,FLAT);
  printf("\t\tnb_subst_cat%s : number of relative substitution rate categories. Default : %snb_subst_cat%s=1.\n",
	 FLAT,LINE,FLAT);
  printf("\t\tMust be a positive integer.\n");

  printf("\n");

  printf("%s\n\t-a (or --alpha) %sgamma%s\n",BOLD,LINE,FLAT);
  printf("\t\tgamma%s : distribution of the gamma distribution shape parameter.\n",FLAT);
  printf("\t\tCan be a fixed positive value or %se%s to get the maximum likelihood estimate.\n",LINE,FLAT);

  printf("\n");

  printf("%s\n\t-s (or --search) %smove%s\n",BOLD,LINE,FLAT);
  printf("\t\tTree topology search operation option.\n");
  printf("\t\tCan be either %sNNI%s (default) or %sNNI+SPR%s or %sSPR%s.\n",LINE,FLAT,LINE,FLAT,LINE,FLAT);
  
  printf("\n");

  printf("%s\n\t-u (or --inputtree) %suser_tree_file%s\n",BOLD,LINE,FLAT);
  printf("\t\tuser_tree_file%s : starting tree filename. The tree must be in Newick format.\n",FLAT);

  printf("\n");

  printf("%s\n\t-o %sparams%s\n",BOLD,LINE,FLAT);
  printf("\t\tThis option focuses on specific parameter optimisation.\n");
  printf("\t\t%sparams%s=tlr : tree topology (t), branch length (l) and rate parameters (r) are optimised.\n",LINE,FLAT);
  printf("\t\t%sparams%s=tl  : tree topology and branch length are optimised.\n",LINE,FLAT);
  printf("\t\t%sparams%s=lr  : branch length and rate parameters are optimised.\n",LINE,FLAT);
  printf("\t\t%sparams%s=l   : branch length are optimised.\n",LINE,FLAT);
  printf("\t\t%sparams%s=r   : rate parameters are optimised.\n",LINE,FLAT);
  printf("\t\t%sparams%s=n   : no parameter is optimised.\n",LINE,FLAT);
  
  printf("\n");

  printf("%s\n\t--rand_start%s\n",BOLD,FLAT);
  printf("\t\tThis option sets the initial tree to random.\n");
  printf("\t\tIt is only valid if SPR searches are to be performed.\n");

  printf("\n");

  printf("%s\n\t--n_rand_starts %snum%s\n",BOLD,LINE,FLAT);
  printf("\t\tnum%s is the number of initial random trees to be used.\n",FLAT);
  printf("\t\tIt is only valid if SPR searches are to be performed.\n");

  printf("\n");

  printf("%s\n\t--r_seed %snum%s\n",BOLD,LINE,FLAT);
  printf("\t\tnum%s is the seed used to initiate the random number generator.\n",FLAT);
  printf("\t\tMust be an integer.\n");

  printf("\n");

  printf("%s\n\t--print_site_lnl%s\n",BOLD,FLAT);
  printf("\t\t%sPrint the likelihood for each site in file *_phyml_lk.txt.\n",FLAT);

  printf("\n");

  printf("%s\n\t--print_trace%s\n",BOLD,FLAT);
  printf("\t\t%sPrint each phylogeny explored during the tree search process\n",FLAT);
  printf("\t\t%sin file *_phyml_trace.txt.\n",FLAT);

  printf("\n");

  printf("%sPHYLIP-LIKE INTERFACE\n""%s\n\tYou can use phyml with no arguments, in this case change the value of\n",BOLD,FLAT);
  printf("%s\ta parameter by typing its corresponding character as shown on screen.\n\n",FLAT);
 
  printf("%sEXAMPLES\n\n"
	 "%s\tDNA interleaved sequence file, default parameters : ""%s  ./phyml -i seqs1"
	 "%s\n\tAA interleaved sequence file, default parameters :  ""%s  ./phyml -i seqs2 -d aa"
	 "%s\n\tAA sequential sequence file, with customization :   ""%s  ./phyml -i seqs3 -q -d aa -m JTT -c 4 -a e\n",BOLD,FLAT,BOLD,FLAT,BOLD,FLAT,BOLD);
  Exit("");
}

/*********************************************************/

#define N_SEQUENCEFILE 1
#define N_DATATYPE 2
#define N_FORMAT 3
#define N_DATASETS 4
#define N_BOOTSTRAPSETS 5
#define N_MODELNAME 6
#define N_KAPPA 7
#define N_PROPORTIONINVAR 7 /*same as kappa*/
#define N_NBCATG 8
#define N_ALPHA 9
#define N_STARTINGTREE 10
#define N_OPT_TOPO 11
#define N_OPT_LENGTHSRATES 12

#define N_NB_PARAMS_DNA 13
#define N_NB_PARAMS_AA 12

option *Get_Input(int argc, char **argv)
{

  option *io;

  io = (option *)Make_Input();
  Set_Defaults_Input(io);
  Set_Defaults_Model(io->mod);
  Set_Defaults_Optimiz(io->mod->s_opt);

  putchar('\n');

  switch (argc)
    {
    case 1:
/* #if defined(MG) || defined(USE_OLD_INTERFACE) */
/*       Get_Input_Interactive(io); */
/* #else */
      Launch_Interface(io);
/* #endif */

      break;
    case 2:
      Usage();
      break;
    default:
      Read_Command_Line(io,argc,argv);
    }

  Print_Settings(io);
  return io;
}

/*********************************************************/



#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "RNAcode.h"
#include "rnaz_utils.h"
#include "treeML.h"
#include "cmdline.h"
#include "code.h"
#include "score.h"
#include "tree.h"
#include "treeSimulate.h"
#include "extreme_fit.h"
#include "postscript.h"
#include "misc.h"

/* Global variables shared troughout the application */

parameters pars; /* user options */
bgModel *models, *modelsRev; /* Background model data for current alignment*/
float**** Sk; /* Main score matrix which is allocated only once and re-used*/

int main(int argc, char *argv[]){

  int i,j,k,x,L,N,hssCount;
  char *tmpSeq, *treeString;
  float kappa;
  TTree* tree;
  segmentStats *results;
  float parMu, parLambda;

  int (*readFunction)(FILE *clust,struct aln *alignedSeqs[]);

  struct aln *inputAln[MAX_NUM_NAMES];
  struct aln *inputAlnRev[MAX_NUM_NAMES];

  pars.omega=-2.0;
  pars.Omega=-4.0;
  pars.Delta=-10.0;
  pars.stopPenalty_0=-9999.0;
  pars.stopPenalty_k=-8.0;
  pars.inputFile=stdin;
  pars.outputFile=stdout;
  pars.debugFile=stdout;
  pars.sampleN=100;
  strcpy(pars.limit,"");
  pars.cutoff=1.0;
  pars.outputFormat=0; /* 0: normal list; 1: GTF; 2:compact list (debugging) */
  strcpy(pars.debugFileName,"");
  strcpy(pars.inputFileName,"STDIN");

  read_commandline(argc, argv);

  srand(time(NULL));

  ntMap['A']=ntMap['a']=0;
  ntMap['C']=ntMap['c']=1;
  ntMap['G']=ntMap['g']=2;
  ntMap['T']=ntMap['t']=3;
  ntMap['U']=ntMap['u']=3;

  switch(checkFormat(pars.inputFile)){
  case CLUSTAL:
    readFunction=&read_clustal;
    break;
  case MAF:
    readFunction=&read_maf;
    break;
  case 0:
    nrerror("ERROR: Unknown alignment file format. Use Clustal W or MAF format.\n");
  }

  while (readFunction(pars.inputFile, inputAln)!=0){
    /* Currently repeat masked regions are ignored and Ns are converted to gaps */
    /* Fix this */
    for (i=0; inputAln[i]!=NULL; i++){
      tmpSeq=inputAln[i]->seq;
      j=0;
      while (tmpSeq[j]){
        tmpSeq[j]=toupper(tmpSeq[j]);
        if (tmpSeq[j]=='N'){
          tmpSeq[j]='-';
        }
        j++;
      }
    }
    
    if (!strcmp(pars.limit,"")==0){
      pruneAln(pars.limit,(struct aln**)inputAln);
    }
    
    //printAlnMAF(stdout,(const struct aln**)inputAln,0); 
    
    //L=strlen(inputAln[0]->seq);

    L=getSeqLength(inputAln[0]->seq);
    
    for (N=0; inputAln[N]!=NULL; N++);
    
    /* Currently minimum number of sequences is 3 because BIONJ seg-faults with 2 */
    /* Fix this that it works with two sequences*/
    if (N<=2){
      fprintf(stderr,"Skipping alignment. There must be at least three sequences in the alignment.\n");
      continue;
    }

    if (L<3){
      fprintf(stderr,"Skipping alignment. Too short.\n");
      continue;
    }
   
    
    if (treeML((const struct aln**)inputAln,&treeString,&kappa)==0){
      fprintf(stderr,"\nSkipping alignment. Failed to build ML tree.\n");
      continue;
    }

    tree=string2tree(treeString);
    
    copyAln((struct aln**)inputAln,(struct aln**)inputAlnRev);
    revAln((struct aln**)inputAlnRev);

    models=getModels(tree,(struct aln**)inputAln,kappa);
    modelsRev=getModels(tree,inputAlnRev,kappa);

    Sk=NULL;

    getExtremeValuePars(tree, (const struct aln**)inputAln, pars.sampleN, &parMu, &parLambda);
    
    results=scoreAln((const struct aln**)inputAln, tree, kappa, parMu, parLambda);
    printResults(pars.outputFile,pars.outputFormat,results);

    //getPairwiseScoreMatrix(models,(const struct aln**)inputAln);
    //getMultipleScoreMatrix(Sk,models,(const struct aln**)inputAln);
    //results=scoreAln((const struct aln**)inputAln, tree, kappa, 0.0, 0.0);

    
    for (k=0;k<N;k++){
      for (x=0;x<3;x++){
        for (i=0;i<L+1;i++){
          free(Sk[k][x][i]);
        }
        free(Sk[k][x]);
      }
      free(Sk[k]);
      } 
    free(Sk);
    Sk=NULL;

    freeSeqgenTree(tree);
    freeResults(results);
    freeModels(models,N);
    freeModels(modelsRev,N);
    freeAln((struct aln**)inputAln);
    freeAln((struct aln**)inputAlnRev);

  }

  exit(EXIT_SUCCESS);

}


void usage(void){
  help();
}

void help(void){

  cmdline_parser_print_version ();

  printf("\nUsage: %s [OPTIONS]... [FILES]\n\n", CMDLINE_PARSER_PACKAGE);
  printf("%s\n","  -h, --help                       Help screen");
  printf("%s\n\n","  -V, --version                    Print version");
}

void version(void){
  //printf("RNAcode v Wed Dec 10 14:53:47 2008" PACKAGE_VERSION "\n");
  printf("RNAcode v Fri Jun 19 12:05:51 2009\n");
  exit(EXIT_SUCCESS);
}

void read_commandline(int argc, char *argv[]){

  struct gengetopt_args_info args;

  if (cmdline_parser (argc, argv, &args) != 0){
    usage();
    exit(EXIT_FAILURE);
  }
  
  if (args.inputs_num >= 1){
    pars.inputFile = fopen(args.inputs[0], "r"); 
    if (pars.inputFile == NULL){
      fprintf(stderr, "ERROR: Can't open input file %s\n", args.inputs[0]);
      exit(1);
    }
    strcpy(pars.inputFileName,args.inputs[0]);
  }
  
  if (args.outfile_given){
    pars.outputFile = fopen(args.outfile_arg, "w");
    if (pars.outputFile == NULL){
      fprintf(stderr, "ERROR: Can't open output file %s\n", args.outfile_arg);
      exit(1);
    }
  }

  if (args.num_samples_given){
    pars.sampleN=args.num_samples_arg;
  }

  if (args.gtf_given){
    pars.outputFormat=1;
  }

  if (args.concise_given){
    pars.outputFormat=2;
  }

  if (args.cutoff_given){
    pars.cutoff=args.cutoff_arg;
  }

  if (args.debug_file_given){
    strcpy(pars.debugFileName,args.debug_file_arg);
  }

  if (args.print_if_below_given){
    pars.printIfBelow=args.print_if_below_arg; 
    pars.debugFile = fopen(pars.debugFileName, "w");
  }

  if (args.print_if_above_given){
    pars.printIfAbove=args.print_if_above_arg; 
    pars.debugFile = fopen(pars.debugFileName, "w");
  }

  if (args.limit_given){
    strcpy(pars.limit,args.limit_arg);
  }

  if (args.help_given){
    help();
    exit(EXIT_SUCCESS);
  }

  if (args.version_given){
    version();
    exit(EXIT_SUCCESS);
  }

  cmdline_parser_free(&args);

}

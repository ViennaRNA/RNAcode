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
#include "utils.h"

/* Global variables shared troughout the application */

parameters pars; /* user options */
bgModel *models, *modelsRev; /* Background model data for current alignment*/
float**** Sk;  /* Main score matrix which is allocated only once and re-used*/

/* Save matrix for native alignment for backtracking afterwards */ 
float**** Sk_native;   
float**** Sk_native_rev;  

long int hitCounter;

void freeModels(bgModel* models, int N);

int main(int argc, char *argv[]){

  int i,j,k,x,L,N,hssCount, alnCounter;
  char *tmpSeq, *treeString;
  float kappa, maxScore;
  TTree* tree;
  segmentStats *results;
  float parMu, parLambda;
  clock_t startTime;	
  float runtime;

  int (*readFunction)(FILE *clust,struct aln *alignedSeqs[]);

  struct aln *inputAln[MAX_NUM_NAMES];
  struct aln *inputAlnRev[MAX_NUM_NAMES];

  pars.Delta=-10.0;
  pars.Omega=-4.0;
  pars.omega=-2.0;
  pars.stopPenalty_k=-8.0;
  pars.stopPenalty_0=-9999.0;
  pars.inputFile=stdin;
  pars.outputFile=stdout;
  pars.debugFile=stdout;
  pars.bestOnly=0;
  pars.bestRegion=0;
  pars.stopEarly=0;
  pars.postscript=0;
  pars.postscript_cutoff=0.05;
  strcpy(pars.postscriptDir, "eps");
  pars.sampleN=100;
  pars.blosum=62;
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

  alnCounter=0;
  hitCounter=0;
  startTime=clock();
  
  while (readFunction(pars.inputFile, inputAln)!=0){
  
    alnCounter++;

    /* Currently repeat masked regions are ignored*/
    /* Fix this */
    for (i=0; inputAln[i]!=NULL; i++){
      tmpSeq=inputAln[i]->seq;
      j=0;
      while (tmpSeq[j]){
        tmpSeq[j]=toupper(tmpSeq[j]);
        j++;
      }
    }
    
    if (!strcmp(pars.limit,"")==0){
      pruneAln(pars.limit,(struct aln**)inputAln);
    }
    
    //printAlnMAF(stdout,(const struct aln**)inputAln,0); 
   
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
    free(treeString);
    
    copyAln((struct aln**)inputAln,(struct aln**)inputAlnRev);
    revAln((struct aln**)inputAlnRev);

    models=getModels(tree,(struct aln**)inputAln,kappa);
    modelsRev=getModels(tree,inputAlnRev,kappa);

    Sk=NULL;
    Sk_native=NULL;
    Sk_native_rev=NULL;

    results=scoreAln((const struct aln**)inputAln, tree, kappa, 1);
    
    hssCount = 0;
    while (results[hssCount++].score > 0.0);

    qsort((segmentStats*) results, hssCount,sizeof(segmentStats),compareScores);

    maxScore=results[0].score;

    if (getExtremeValuePars(tree, (const struct aln**)inputAln, pars.sampleN, maxScore, &parMu, &parLambda) == 1){
      for (i=0;i<hssCount;i++){
        results[i].pvalue=1-exp((-1)*exp((-1)*parLambda*(results[i].score-parMu)));
      }
    } else {
      for (i=0;i<hssCount;i++){
        results[i].pvalue=99.0;
      }
    }

    printResults(pars.outputFile,pars.outputFormat,(const struct aln**)inputAln, results);

    for (k=0;k<N;k++){
      for (x=0;x<3;x++){
        for (i=0;i<L+1;i++){
          free(Sk[k][x][i]);
          free(Sk_native[k][x][i]);
          free(Sk_native_rev[k][x][i]);
        }
        free(Sk[k][x]);
        free(Sk_native[k][x]);
        free(Sk_native_rev[k][x]);
      }
      free(Sk[k]);
      free(Sk_native[k]);
      free(Sk_native_rev[k]);
    } 
    free(Sk);
    free(Sk_native);
    free(Sk_native_rev);
    Sk=NULL;
    Sk_native=NULL;
    Sk_native_rev=NULL;

    freeSeqgenTree(tree);
    freeResults(results);
    freeModels(models,N);
    freeModels(modelsRev,N);
    freeAln((struct aln**)inputAln);
    freeAln((struct aln**)inputAlnRev);

  }

  if (pars.outputFormat==0){

    runtime = (float)(clock() - startTime) / CLOCKS_PER_SEC;

    fprintf(pars.outputFile,
            "\n%i alignment(s) scored in %.2f seconds. Parameters used:\nN=%i, Delta=%.2f, Omega=%.2f, omega=%.2f, stop penalty=%.2f\n\n", alnCounter,runtime,pars.sampleN, pars.Delta,pars.Omega,pars.omega,pars.stopPenalty_k);
  }

  exit(EXIT_SUCCESS);

}


void usage(void){
  help();
}

void help(void){

  cmdline_parser_print_version ();

  printf("\nUsage: %s [OPTIONS]... [FILE]\n\n", CMDLINE_PARSER_PACKAGE);
  printf("--outfile     -o  Output file  (default: stdout)\n");
  printf("--gtf         -g  Format output as GTF\n");
  printf("--tabular     -t  Format output as tab delimited fields\n");
  printf("--best-only   -b  Show only best hit\n");
  printf("--best-region -r  Show only best non-overlapping hits\n");
  printf("--stop-early  -s  Don't calculate p-values for hits likely to be above cutoff\n");
  printf("--num-samples -n  Number of samples to calculate p-value (default: 100)\n");
  printf("--cutoff      -p  p-value cutoff (default: 1.0)\n");
  printf("--pars        -c  Parameters as comma separated string (see README for details)\n");
  printf("--eps         -e  Create colored plots in EPS format\n");
  printf("--eps-cutoff  -i  Create plots only if p better than this cutoff (default: 0.05)\n");
  printf("--eps-dir     -d  Directory to put eps-files (default: eps)\n");
  printf("--help        -h  Print this help screen\n");
  printf("--version     -v  Print version\n\n");
}

void version(void){
  printf("RNAcode v " PACKAGE_VERSION " Wed Jan  5 13:14:57 EST 2011\n");
  exit(EXIT_SUCCESS);
}

void read_commandline(int argc, char *argv[]){

  int x;

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

  if (args.blosum_given){
    if (args.blosum_arg !=62 && args.blosum_arg !=90){
      fprintf(stderr, "ERROR: Currently only BLOSUM62 and BLOSUM90 are supported.\n");
      exit(1);
    } else {
      pars.blosum=args.blosum_arg;
    }
  }

  if (args.gtf_given){
    pars.outputFormat=1;
  }

  if (args.tabular_given){
    pars.outputFormat=2;
  }

  if (args.pars_given){

    x=sscanf(args.pars_arg,"%f,%f,%f,%f",&pars.Delta,&pars.Omega,&pars.omega, &pars.stopPenalty_0);

    if (!x){
      fprintf(stderr, "ERROR: Format error in parameter string. Refer to README how to use --pars.\n");
      exit(1);
    }
  }
  
  if (args.best_only_given){
    pars.bestOnly=1;
  }

  if (args.best_region_given){
    pars.bestRegion=1;
  }

  if (args.stop_early_given){
    pars.stopEarly=1;    
  }

  if (args.cutoff_given){
    pars.cutoff=args.cutoff_arg;
  }

  if (args.eps_cutoff_given){
    pars.postscript_cutoff = args.eps_cutoff_arg;
  }

  if (args.eps_given){
    pars.postscript = 1;
  }

  if (args.eps_dir_given){
    strcpy(pars.postscriptDir, args.eps_dir_arg);
  }

  if (args.debug_file_given){
    strcpy(pars.debugFileName,args.debug_file_arg);
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

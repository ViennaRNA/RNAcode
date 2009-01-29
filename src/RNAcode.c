#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "rnaz_utils.h"
#include "treeML.h"
#include "cmdline.h"
#include "code.h"
#include "tree.h"
#include "treeSimulate.h"
#include "extreme_fit.h"
#include "postscript.h"

void usage(void);
void help(void);
void version(void);

int main(int argc, char *argv[]){

  int i,j,k,L,N,hssCount;
 
  char inputFileName[1024]="STDIN";
  char *tmpSeq, *treeString;
  double kappa, parMu, parLambda;
  bgModel* models;
  segmentStats *results;
  TTree* tree;
  char debugFileName[1024]="";
  FILE *inputFile=stdin;
  FILE *outputFile=stdout;
  FILE *debugFile=stdout;

  int (*readFunction)(FILE *clust,struct aln *alignedSeqs[]);

  struct gengetopt_args_info args;
  struct aln *inputAln[MAX_NUM_NAMES];

  int sampleN=1000;
  int sampleMode=1;
  int outputFormat=0; /* 0: normal list; 1: GTF; 2:compact list (debugging) */
  double cutoff=1.0;

  /* For debugging purposes, set these variable to print out
     false-positive or false-negatives on test sets of noncoding or
     coding regions of known  annotation */

  double printIfBelow=-1.0; 
  double printIfAbove=-1.0;

  srand(time(NULL));

  ntMap['A']=ntMap['a']=0;
  ntMap['C']=ntMap['c']=1;
  ntMap['G']=ntMap['g']=2;
  ntMap['T']=ntMap['t']=3;
  ntMap['U']=ntMap['u']=3;

  /* Read command line arguments */
  
  if (cmdline_parser (argc, argv, &args) != 0){
    usage();
    exit(EXIT_FAILURE);
  }
  
  if (args.inputs_num>=1){
    inputFile = fopen(args.inputs[0], "r"); 
    if (inputFile == NULL){
      fprintf(stderr, "ERROR: Can't open input file %s\n", args.inputs[0]);
      exit(1);
    }
    strcpy(inputFileName,args.inputs[0]);
  }
  
  if (args.outfile_given){
    outputFile = fopen(args.outfile_arg, "w");
    if (outputFile == NULL){
      fprintf(stderr, "ERROR: Can't open output file %s\n", args.outfile_arg);
      exit(1);
    }
  }

  if (args.num_samples_given){
    sampleN=args.num_samples_arg;
  }


  if (args.gtf_given){
    outputFormat=1;
  }

  if (args.concise_given){
    outputFormat=2;
  }

  if (args.cutoff_given){
    cutoff=args.cutoff_arg;
  }

  if (args.debug_file_given){
    strcpy(debugFileName,args.debug_file_arg);
  }

  if (args.print_if_below_given){
      printIfBelow=args.print_if_below_arg; 
      debugFile = fopen(debugFileName, "w");
  }

  if (args.print_if_above_given){
    printIfAbove=args.print_if_above_arg; 
    debugFile = fopen(debugFileName, "w");
  }

  if (args.help_given){
    help();
    exit(EXIT_SUCCESS);
  }

  if (args.version_given){
    version();
    exit(EXIT_SUCCESS);
  }

  switch(checkFormat(inputFile)){
  case CLUSTAL:
    readFunction=&read_clustal;
    break;
  case MAF:
    readFunction=&read_maf;
    break;
  case 0:
    nrerror("ERROR: Unknown alignment file format. Use Clustal W or MAF format.\n");
  }

  while (readFunction(inputFile, inputAln)!=0){
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

    if (args.limit_given){
      pruneAln(args.limit_arg,(struct aln**)inputAln);
    }
    
    
    //printAlnMAF(stdout,(const struct aln**)inputAln,0); 

    stripGaps((struct aln**)inputAln);

    L=strlen(inputAln[0]->seq);

    if (L<3){
      continue;
    }

    for (N=0; inputAln[N]!=NULL; N++);   


    /* Currently minimum number of sequences is 3 because BIONJ seg-faults with 2 */
    /* Fix this that it works with two sequences*/

    if (N<=2){
      //fprintf(stderr,"There must be at least three sequences in the alignment.\n");
      continue;
      //exit(1);
    }

    
    if (treeML((const struct aln**)inputAln,&treeString,&kappa)==0){
      fprintf(stderr,"\nFailed to build ML tree.\n");
      continue; 
    }
  
    tree=string2tree(treeString);

    models=getModels(tree,inputAln,kappa);

    getExtremeValuePars(tree, models, (const struct aln**)inputAln, sampleN, sampleMode, &parMu, &parLambda);
    results=getHSS(models, (const struct aln**)inputAln, parMu, parLambda,cutoff);
    
    hssCount=0;
    
    while (results[hssCount].score>=0) hssCount++;

    qsort((segmentStats*) results, hssCount,sizeof(segmentStats),compareScores);

    if (printIfAbove > -1.0){
      if (results[0].pvalue > printIfAbove){
        printAlnMAF(debugFile,(const struct aln**)inputAln,0);
      }
    }
      
    if (printIfBelow > -1.0){
      if (results[0].pvalue < printIfBelow){
        printAlnMAF(debugFile,(const struct aln**)inputAln,0);
      }
    }
     
    printResults(outputFile,outputFormat,results);

    if (args.gfx_given){
      colorAln("color.ps", (const struct aln**)inputAln, results[0]);
    }

    freeResults(results);

    freeModels(models,N);

    free(treeString);

    freeSeqgenTree(tree);

    freeAln((struct aln**)inputAln);

  }
  
  if ((printIfAbove > -1.0) || (printIfBelow > -1.0)){
    fclose(debugFile);
  }


  fclose(inputFile);

  cmdline_parser_free(&args);
  
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
  printf("RNAcode v Thu Jan 22 15:07:09 2009\n");
  exit(EXIT_SUCCESS);
}


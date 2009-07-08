#include <stdio.h>
#include <strings.h>
#include <malloc.h>
#include <math.h>

float *weights(float **D, int N) { 
  float **DD, s;
  float *aans;
  float MINDIST = 0.000001;
  int i,j;

  /** working arrays **/
  DD = (float **) malloc(sizeof(float*)*N);
  for(i=0;i<N;i++) {
    DD[i] = (float *) malloc(sizeof(float)*N);
    for(j=0;j<N;j++) {
      DD[i][j] = D[i][j];
      if( (i!=j)&&(D[i][j]<MINDIST) ) DD[i][j] = MINDIST;      
    }
  }  
  aans = (float *) malloc(sizeof(float)*N);
  for(i=0;i<N;i++) { aans[i]= 1.0; } 

  /** call Gauss elimination **/
  (void)  linsys(DD,aans,N);

  /** normalize **/
  //for(s=0.0,i=1;i<N;i++) s += aans[i];
  //for(i=0;i<N;i++) aans[i] /= s;
 
  for(i=0;i<N;i++) free(DD[i]);
  free(DD);
  return aans;
}
  

/* ------------------------------------------------------------------ */

/*
  imported from:
  
  http://www.fizik.itu.edu.tr/guvenh/CODES/gauss.c

  Programname: lin.c
  K. Atkinson, Elementary Numerical Analysis, p.221, C version by H. Guven

  interface modified.

*/

/*------------------------
 THIS FUNCTION SOLVES A SYSTEM OF LINEAR EQUATIONS

 A*X =  B

 THE METHOD USED IS GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING

 INPUT:
 THE COEFFICIENT MATRIX A IS STORED IN THE ARRAY a[][].
 THE RIGHT SIDE CONSTANTS ARE IN THE ARRAY b[].
 THE ORDER OF THE LINEAR SYSTEM IS Norder.
 THE SIZE OF MAXPIV, GIVEN BELOW, MUST BE GREATER THEN
 Norder. IF NOT, IT IS A FATAL ERROR.

 OUTPUT:
 THE ARRAY b[] CONSTANTS THE SOLUTION X
 THE MATRIX a[][] MADE UPPER TRIANGULAR MATRIX U OBTAINED
 BY ELIMINATION. THE ROW MULTIPLIERS USED IN THE ELIMINATION
 ARE STORED IN THE LOWER TRIANGULAR PART OF a[][].
 ----------------     */

int linsys(float **a, float *b,  int N )
{
  float absa, amax, temp, mult,sum;
  int i, j, k;
  int *pivot;

  pivot = (int *) malloc(sizeof(int)*N);
 
  
/*     BEGIN ELIMINATION STEPS */

 for (k = 0; k < N; k++) {

/*        CHOOSE PIVOT ROW */

 pivot[k] = k;
 amax= fabs(a[k][k]);

 for (i = k + 1; i < N; i++) {
 absa = fabs(a[i][k]);
 if (absa > amax) {
 pivot[k] = i;
 amax = absa;
 }
 }

 if (amax == 0.) goto endf;


//     COEFFICIENTS MATRIX IS SINGULAR


 if (pivot[k] != k) {

//      SWITCH ROWS K AND PIVOT(K)

 i    = pivot[k];
 temp = b[k];
 b[k] = b[i];
 b[i] = temp;

 for (j = k; j < N; j++) {
 temp = a[k][j];
 a[k][j] = a[i][j];
 a[i][j] = temp;
 }
 }

/*      PERFORM STEP #K OF ELIMINATION */

 for (i = k + 1; i < N; i++) {
 mult = a[i][k] / a[k][ k];
 a[i][k] = mult;
 b[i] =b[i] - mult * b[k];

 for (j = k + 1; j < N; j++)
 a[i][j] = a[i][j] - mult * a[k][j];
 }

 }
 if (a[N-1][N-1] == 0.) goto endf; // COEFFICIENT MATRIX IS SINGULAR


/*  SOLVE FOR SOLUTION X USING BACK SUBSTITUTION */

 for (i= N-1; i>= 0; i--) {
 sum = 0.;
 for (j = i+1; j < N; j++)
 sum =sum+ a[i][j] * b[j];

 b[i] = (b[i] - sum)/a[i][i];
 }

 endf:
return 1;
}


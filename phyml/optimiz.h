/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef OPTIMIZ_H
#define OPTIMIZ_H

void      Optimiz_Ext_Br(arbre *tree);
void      Optimize_Alpha(arbre *tree);
void      Optimize_Kappa(arbre *tree);
void      Optimize_Lambda(arbre *tree);
void      Optimize_Param_Parall(arbre *tree);
phydbl    Optimize_Branch_Quad(arbre *tree, allseq *alldata, edge *b_fcus);
void      Optimize_After_Hide(arbre *tree, allseq *alldata, node *h);
void      Round_Optimize(arbre *tree, allseq *data);
int       Dist_Seq_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
			phydbl *fa, phydbl *fb, phydbl *fc, 
			allseq *data, int num1, int num2, model *mod);
phydbl    Dist_Seq_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			 phydbl *xmin, allseq *data, 
			 int num1, int num2, model *mod);
phydbl    Kappa_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		       phydbl *xmin, arbre *tree, allseq *alldata);
phydbl    Lambda_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			phydbl *xmin, arbre *tree, allseq *alldata);
phydbl    Alpha_Golden_Br_Opt(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			      phydbl *xmin, arbre *tree, allseq *alldata, 
			      int n_opt, phydbl *init_l);
phydbl    Alpha_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol,phydbl *xmin, 
		       arbre *tree, allseq *alldata);
phydbl    Br_Len_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			phydbl *xmin, edge *b_fcus, arbre *tree);
phydbl    Br_Len_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		       edge *b_fcus, arbre *tree, int n_iter_max);
int       Br_Len_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		      phydbl *fa, phydbl *fb, phydbl *fc, 
		      edge *b_fcus, arbre *tree);
phydbl    Optimize_Path_Length(model *mod, allseq *alldata, edge *a, 
			       int lra, edge *b, int lrb, phydbl i_len);
void      Optimize_Param_Serie(node *a, node *d, edge *b_fcus, arbre *tree, 
			       allseq *alldata, int n_passes);
phydbl    Optimize_Dist(model *mod, phydbl init, allseq *twoseqs);
phydbl    Pinvar_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			phydbl *xmin, arbre *tree, allseq *alldata, int n_iter_max);
void      Optimize_Pinvar(arbre *tree);
int       Lambda_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		      phydbl *fa, phydbl *fb, phydbl *fc, 
		      arbre *tree);
int       Kappa_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		      phydbl *fa, phydbl *fb, phydbl *fc, 
		      arbre *tree);
int       Alpha_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		      phydbl *fa, phydbl *fb, phydbl *fc, 
		      arbre *tree);
int       Pinvar_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		      phydbl *fa, phydbl *fb, phydbl *fc, 
		      arbre *tree);
void Optimiz_All_Free_Param(arbre *tree, int verbose);
void      Optimiz_RRparam_GTR(arbre *tree, int num_param);
phydbl    RRparam_GTR_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		   	     phydbl *xmin, arbre *tree, allseq *alldata, phydbl *param, int n_iter_max);

int Powell_GTR_Param(arbre *tree, phydbl *p, int n, phydbl ftol);
phydbl Linmin_GTR_Param(arbre *tree,phydbl *p, phydbl *xi, int n);
phydbl F1dim(arbre *tree, phydbl x, phydbl *p, phydbl *xi, phydbl n);
int Mnbrak_1dim(phydbl *ax, phydbl *bx, phydbl *cx, 
		phydbl *fa, phydbl *fb, phydbl *fc,
		arbre *tree,
		phydbl *p,  phydbl *xi, phydbl n);
phydbl Brent_1dim(phydbl ax, phydbl bx, phydbl cx, 
		  phydbl tol, phydbl *xmin,
		  arbre *tree,
		  phydbl *p, phydbl *xi, phydbl n);

int Min_With_Derivatives(arbre *tree, phydbl *p, int n, phydbl ftol, phydbl step_size, 
			 phydbl (*func) (), void (*dfunc)(), phydbl (*linmin)());
void BFGS(arbre *tree, phydbl *p, int n, phydbl gtol, phydbl step_size,
	  phydbl(*func)(), void (*dfunc)(), void (*lnsrch)(),int *failed);
void Lnsrch_RR_Param(arbre *tree, int n, phydbl *xold, phydbl fold, phydbl *g, phydbl *p, phydbl *x,
		     phydbl *f, phydbl stpmax, int *check);
void Optimize_Single_Param_Generic(arbre *tree, phydbl *param, phydbl lim_inf, phydbl lim_sup, phydbl tol, int n_max_iter);
int Generic_Brak(phydbl *param,
		 phydbl *ax, phydbl *bx, phydbl *cx, 
		 phydbl *fa, phydbl *fb, phydbl *fc,
		 phydbl lim_inf, phydbl lim_sup,
		 arbre *tree);
phydbl Generic_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		     phydbl *xmin, arbre *tree, int n_iter_max);
void Optimize_Br_Len_Serie(node *a, node *d, edge *b_fcus, arbre *tree,allseq *alldata);
void Lnsrch_Nucleotide_Frequencies(arbre *tree, int n, phydbl *xold, 
				   phydbl fold, phydbl *g, phydbl *p, phydbl *x,
				   phydbl *f, phydbl stpmax, int *check);

void Optimize_Global_Rate(arbre *tree);
phydbl Br_Len_Brent_Default(edge *b_fcus, arbre *tree);

void EM_Dist(model *mod, allseq *data);
phydbl Dist_F_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
		    phydbl *param, phydbl *F, model *mod);
int Dist_F_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *F, phydbl *param, model *mod);
void Opt_Dist_F(phydbl *dist, phydbl *F, model *mod);
phydbl Missing_Dist_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
			  int x, int y, matrix *mat);
int Missing_Dist_Brak(phydbl *ax, phydbl *bx, phydbl *cx, int x, int y, matrix *mat);
void Opt_Missing_Dist(int x, int y, matrix *mat);
int Optimiz_Alpha_And_Pinv(arbre *tree);
void Lnsrch_RR_Cov_Param(arbre *tree, int n, phydbl *xold, phydbl fold, 
			 phydbl *g, phydbl *p, phydbl *x,
			 phydbl *f, phydbl stpmax, int *check);
phydbl Node_Time_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
		       node *anc, node *des, arbre *tree, int n_iter_max);
phydbl Time_Stamps_Mult_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
			      arbre *tree, int n_iter_max);
phydbl Branch_Rate_Shape_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			       phydbl *xmin, arbre *tree, int n_iter_max);


#endif


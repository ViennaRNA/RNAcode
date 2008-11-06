/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "bionj.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "m4.h"
#include "mc.h"

#ifdef MG
#include "mg.h"
#endif


/*********************************************************/

double Uni()
{
  double r; 
  r=rand();
  r/=RAND_MAX;
  return r;
}

/*********************************************************************/
/********************* random Gamma generator ************************
* Properties:
* (1) X = Gamma(alpha,lambda) = Gamma(alpha,1)/lambda
* (2) X1 = Gamma(alpha1,1), X2 = Gamma(alpha2,1) independent
*     then X = X1+X2 = Gamma(alpha1+alpha2,1)
* (3) alpha = k = integer then
*     X = Gamma(k,1) = Erlang(k,1) = -sum(log(Ui)) = -log(prod(Ui))
*     where U1,...Uk iid uniform(0,1)
*
* Decompose alpha = k+delta with k = [alpha], and 0<delta<1
* Apply (3) for Gamma(k,1)
* Apply Ahrens-Dieter algorithm for Gamma(delta,1)
*/
 
double Ahrensdietergamma(double alpha)
{
  double x = 0.;

  if (alpha>0.) 
    {
      double y = 0.;
      double b = (alpha+exp(1.))/exp(1.);
      double p = 1./alpha;
      int go = 0;
      while (go==0) 
	{
	  double u = Uni();
	  double w = Uni();
	  double v = b*u;
	  if (v<=1.) 
	    {
	      x = pow(v,p);
	      y = exp(-x);
	    }
	  else 
	    {
	      x = -log(p*(b-v));
	      y = pow(x,alpha-1.);
	    }
	  go = (w<y); // x is accepted when go=1
	}
    }
  return x;
}

/*********************************************************/

double Rgamma(double shape, double scale)
{
  int i;
  double x1 = 0.;
  double delta = shape;
  if (shape>=1.) 
    {
      int k = (int) floor(shape);
      delta = shape - k;
      double u = 1.;
      for (i=0; i<k; i++)
	u *= Uni();
      x1 = -log(u);
    }
  double x2 = Ahrensdietergamma(delta);
  return (x1 + x2)*scale;
}

double Rexp(double lambda)
{
  return -log(Uni()+1.E-30)/lambda;
}

/*********************************************************/
/* NUMERICAL RECIPES ROUTINES FOR COMPUTING C(n,k)       */
/*********************************************************/

phydbl bico(int n, int k)
{
  return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

phydbl factln(int n)
{
  static phydbl a[101];
  
  if (n < 0){ Warn_And_Exit("Err: negative factorial in routine FACTLN"); }
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
  else return gammln(n+1.0);
}

/*********************************************************/

phydbl gammln(phydbl xx)
{
  phydbl x,tmp,ser;
  static phydbl cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;
  
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*(phydbl)log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) 
    {
      x += 1.0;
      ser += cof[j]/x;
    }
  return -tmp+(phydbl)log(2.50662827465*ser);
}

/*********************************************************/
/*          END OF NUMERICAL RECIPES ROUTINES            */
/*********************************************************/

phydbl Pbinom(int N, int ni, phydbl p)
{
  return bico(N,ni)*pow(p,ni)*pow(1-p,N-ni);
}

/*********************************************************/

void Plim_Binom(phydbl pH0, int N, phydbl *pinf, phydbl *psup)
{
  *pinf = pH0 - 1.64*sqrt(pH0*(1-pH0)/(phydbl)N);
  if(*pinf < 0) *pinf = .0;
  *psup = pH0 + 1.64*sqrt(pH0*(1-pH0)/(phydbl)N);
}

/*********************************************************/

/*wash. commented out following functions since they also exist in seq-gen gamma.c*/

/* phydbl LnGamma (phydbl alpha) */
/* { */
/* /\* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places. */
/*    Stirling's formula is used for the central polynomial part of the procedure. */
/*    Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function. */
/*    Communications of the Association for Computing Machinery, 9:684 */
/* *\/ */
/*    phydbl x=alpha, f=0, z; */

/*    if (x<7) { */
/*       f=1;  z=x-1; */
/*       while (++z<7)  f*=z; */
/*       x=z;   f=-(phydbl)log(f); */
/*    } */
/*    z = 1/(x*x); */
/*    return  f + (x-0.5)*(phydbl)log(x) - x + .918938533204673 */
/* 	  + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z */
/* 	       +.083333333333333)/x; */
/* } */


/*********************************************************/

/* phydbl IncompleteGamma(phydbl x, phydbl alpha, phydbl ln_gamma_alpha) */
/* { */
/* /\* returns the incomplete gamma ratio I(x,alpha) where x is the upper */
/* 	   limit of the integration and alpha is the shape parameter. */
/*    returns (-1) if in error */
/*    ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant. */
/*    (1) series expansion     if (alpha>x || x<=1) */
/*    (2) continued fraction   otherwise */
/*    RATNEST FORTRAN by */
/*    Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics, */
/*    19: 285-287 (AS32) */
/* *\/ */
/*    int i; */
/*    phydbl p=alpha, g=ln_gamma_alpha; */
/*    phydbl accurate=1e-8, overflow=1e30; */
/*    phydbl factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6]; */

/*    if (x==0) return (0); */
/*    if (x<0 || p<=0) return (-1); */

/*    factor=(phydbl)exp(p*(phydbl)log(x)-x-g); */
/*    if (x>1 && x>=p) goto l30; */
/*    /\* (1) series expansion *\/ */
/*    gin=1;  term=1;  rn=p; */
/*  l20: */
/*    rn++; */
/*    term*=x/rn;   gin+=term; */

/*    if (term > accurate) goto l20; */
/*    gin*=factor/p; */
/*    goto l50; */
/*  l30: */
/*    /\* (2) continued fraction *\/ */
/*    a=1-p;   b=a+x+1;  term=0; */
/*    pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b; */
/*    gin=pn[2]/pn[3]; */
/*  l32: */
/*    a++;  b+=2;  term++;   an=a*term; */
/*    for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i]; */
/*    if (pn[5] == 0) goto l35; */
/*    rn=pn[4]/pn[5];   dif=fabs(gin-rn); */
/*    if (dif>accurate) goto l34; */
/*    if (dif<=accurate*rn) goto l42; */
/*  l34: */
/*    gin=rn; */
/*  l35: */
/*    for (i=0; i<4; i++) pn[i]=pn[i+2]; */
/*    if (fabs(pn[4]) < overflow) goto l32; */
/*    for (i=0; i<4; i++) pn[i]/=overflow; */
/*    goto l32; */
/*  l42: */
/*    gin=1-factor*gin; */

/*  l50: */
/*    return (gin); */
/* } */


/*********************************************************/

/* phydbl PointChi2 (phydbl prob, phydbl v) */
/* { */
/* /\* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v */
/*    returns -1 if in error.   0.000002<prob<0.999998 */
/*    RATNEST FORTRAN by */
/*        Best DJ & Roberts DE (1975) The percentage points of the */
/*        Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91) */
/*    Converted into C by Ziheng Yang, Oct. 1993. */
/* *\/ */
/*    phydbl e=.5e-6, aa=.6931471805, p=prob, g; */
/*    phydbl xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6; */

/*    if (p<.000002 || p>.999998 || v<=0) return (-1); */

/*    g = LnGamma (v/2); */
/*    xx=v/2;   c=xx-1; */
/*    if (v >= -1.24*(phydbl)log(p)) goto l1; */

/*    ch=pow((p*xx*(phydbl)exp(g+xx*aa)), 1/xx); */
/*    if (ch-e<0) return (ch); */
/*    goto l4; */
/* l1: */
/*    if (v>.32) goto l3; */
/*    ch=0.4;   a=(phydbl)log(1-p); */
/* l2: */
/*    q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch)); */
/*    t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2; */
/*    ch-=(1-(phydbl)exp(a+g+.5*ch+c*aa)*p2/p1)/t; */
/*    if (fabs(q/ch-1)-.01 <= 0) goto l4; */
/*    else                       goto l2; */

/* l3: */
/*    x=PointNormal (p); */
/*    p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0); */
/*    if (ch>2.2*v+6)  ch=-2*((phydbl)log(1-p)-c*(phydbl)log(.5*ch)+g); */
/* l4: */
/*    q=ch;   p1=.5*ch; */
/*    if ((t=IncompleteGamma (p1, xx, g))<0) { */
/*       printf ("\nerr IncompleteGamma"); */
/*       return (-1); */
/*    } */
/*    p2=p-t; */
/*    t=p2*(phydbl)exp(xx*aa+g+p1-c*(phydbl)log(ch)); */
/*    b=t/ch;  a=0.5*t-b*c; */

/*    s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420; */
/*    s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520; */
/*    s3=(210+a*(462+a*(707+932*a)))/2520; */
/*    s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040; */
/*    s5=(84+264*a+c*(175+606*a))/2520; */
/*    s6=(120+c*(346+127*c))/5040; */
/*    ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6)))))); */
/*    if (fabs(q/ch-1) > e) goto l4; */

/*    return (ch); */
/* } */

/* /\*********************************************************\/ */

/* phydbl PointNormal (phydbl prob) */
/* { */
/* /\* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12) */
/*    returns (-9999) if in error */
/*    Odeh RE & Evans JO (1974) The percentage points of the normal distribution. */
/*    Applied Statistics 22: 96-97 (AS70) */

/*    Newer methods: */
/*      Wichura MJ (1988) Algorithm AS 241: the percentage points of the */
/*        normal distribution.  37: 477-484. */
/*      Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage */
/*        points of the normal distribution.  26: 118-121. */

/* *\/ */
/*    phydbl a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245; */
/*    phydbl a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495; */
/*    phydbl b2=.531103462366, b3=.103537752850, b4=.0038560700634; */
/*    phydbl y, z=0, p=prob, p1; */

/*    p1 = (p<0.5 ? p : 1-p); */
/*    if (p1<1e-20) return (-9999); */

/*    y = sqrt ((phydbl)log(1/(p1*p1))); */
/*    z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0); */
/*    return (p<0.5 ? -z : z); */
/* } */
/* /\*********************************************************\/ */

/* int DiscreteGamma (phydbl freqK[], phydbl rK[], */
/* 		   phydbl alfa, phydbl beta, int K, int median) */
/* { */
/*   /\* discretization of gamma distribution with equal proportions in each */
/*      category */
/*   *\/ */
   
/*   int i; */
/*   phydbl gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1; */

/*   if(K==1) */
/*     { */
/*       rK[0] = 1.0; */
/*       return 0; */
/*     } */

/*    if (median)  */
/*      { */
/*        for (i=0; i<K; i++)     rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta); */
/*        for (i=0,t=0; i<K; i++) t+=rK[i]; */
/*        for (i=0; i<K; i++)     rK[i]*=factor/t; */
/*      } */
/*    else { */
/*       lnga1=LnGamma(alfa+1); */
/*       for (i=0; i<K-1; i++) */
/* 	 freqK[i]=PointGamma((i+1.0)/K, alfa, beta); */
/*       for (i=0; i<K-1; i++) */
/* 	 freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1); */
/*       rK[0] = freqK[0]*factor; */
/*       rK[K-1] = (1-freqK[K-2])*factor; */
/*       for (i=1; i<K-1; i++)  rK[i] = (freqK[i]-freqK[i-1])*factor; */
/*    } */
/*    for (i=0; i<K; i++) freqK[i]=1.0/K; */
/*    return (0); */
/* } */

/*********************************************************/

phydbl CDF_Normal(phydbl x, phydbl mean, phydbl var)
{
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;
  
  x = (x-mean)/var;
  
  if(x >= 0.0) 
    {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
	      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
  else 
    {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
	       ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}

/*********************************************************/

phydbl CDF_Gamma(phydbl x, phydbl mean, phydbl var)
{
  phydbl scale,shape;

  scale = var / mean;
  shape = mean * mean / var;

/*   scale = mean; */
/*   shape = var; */

  printf("\n. shape=%f, scale=%f",shape,scale);
  return IncompleteGamma(x/scale,shape,LnGamma(shape));
}

/*********************************************************/

phydbl Rand_Normal_Deviate(phydbl mean, phydbl sd)
{
  int i;
  phydbl x=.0;
  For(i,12) x += (phydbl)rand()/(RAND_MAX);
  return sd * (x-6.0) + mean;
}

/*********************************************************/

arbre *Read_Tree(char *s_tree)
{
  char **subs;
  int i,n_ext,n_int,n_otu;
  arbre *tree;
  int degree;


  n_otu=0;
  For(i,(int)strlen(s_tree)) if(s_tree[i] == ',') n_otu++;
  n_otu+=1;

  tree = (arbre *)Make_Tree(n_otu);
  Init_Tree(tree,tree->n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);
  Make_Tree_Path(tree);
  Make_List_Of_Reachable_Tips(tree);

  tree->noeud[n_otu]->num = n_otu;
  tree->noeud[n_otu]->tax = 0;

  subs = Sub_Trees(s_tree,&degree);
  Clean_Multifurcation(subs,degree,3);
  if(degree == 2) Unroot_Tree(subs);
  degree = 3;

  tree->has_branch_lengths = 0;
  tree->num_curr_branch_available = 0;
  n_int = n_ext = 0;
  For(i,degree) R_rtree(s_tree,subs[i],tree->noeud[n_otu],tree,&n_int,&n_ext);

  For(i,NODE_DEG_MAX) Free(subs[i]);
  Free(subs);
  return tree;
}

/*********************************************************/

/* void Make_All_Edges_Light(node *a, node *d, int *curr_num_edge) */
/* { */
/*   int i; */

/*   Make_Edge_Light(a,d,*curr_num_edge); */
/*   (*curr_num_edge)++; */
/*   if(d->tax) return; */
/*   else */
/*     { */
/*       For(i,3) */
/* 	{ */
/* 	  if(d->v[i] != a) */
/* 	    Make_All_Edges_Light(d,d->v[i],curr_num_edge); */
/* 	} */
/*     } */
/* } */

/*********************************************************/

void Make_All_Edges_Lk(node *a, node *d, arbre *tree)
{
  int i;

  For(i,3) if((a->v[i]) && (a->v[i] == d)) Make_Edge_Lk(a->b[i],tree);
  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Make_All_Edges_Lk(d,d->v[i],tree);
	}
    }
}

/*********************************************************/
/* 'a' in node a stands for ancestor. 'd' stands for descendant */ 
void R_rtree(char *s_tree_a, char *s_tree_d, node *a, arbre *tree, int *n_int, int *n_ext)
{
  int i;
  node *d;
  int n_otu = tree->n_otu;

  if(strstr(s_tree_a," ")) Warn_And_Exit("\n Err : tree must not contain a ' ' character\n");

  if(s_tree_d[0] == '(')
    {
      char **subs;
      int degree;

      (*n_int)+=1;
      d      = tree->noeud[n_otu+*n_int];
      d->num = n_otu+*n_int;
      d->tax = 0;

      Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
      Read_Branch_Length(s_tree_d,s_tree_a,tree);

      For(i,3)
       {
	 if(!a->v[i])
	   {
	     a->v[i]=d;
	     d->l[0]=a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
	     break;
	   }
       }
      d->v[0]=a;

      Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
      tree->num_curr_branch_available++;

      subs=Sub_Trees(s_tree_d,&degree);
      Clean_Multifurcation(subs,degree,2);
      R_rtree(s_tree_d,subs[0],d,tree,n_int,n_ext);
      R_rtree(s_tree_d,subs[1],d,tree,n_int,n_ext);
      For(i,NODE_DEG_MAX) Free(subs[i]);
      Free(subs);
    }

  else
    {
      int i;

      d      = tree->noeud[*n_ext];
      d->tax = 1;

      Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]); 
      Read_Branch_Length(s_tree_d,s_tree_a,tree);
      Read_Node_Name(d,s_tree_d,tree);
      
      For(i,3)
	{
	 if(!a->v[i])
	   {
	     a->v[i]=d;
	     d->l[0]=a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
	     break;
	   }
	}
      d->v[0]=a;

      Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
      tree->num_curr_branch_available++;
      
      d->num=*n_ext;
      (*n_ext)+=1;
    }
}

/*********************************************************/

void Read_Branch_Label(char *s_d, char *s_a, edge *b)
{
  char *sub_tp;
  char *p;
  int i,pos;

  sub_tp = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  strcpy(sub_tp,s_d);
  strcat(sub_tp,"#");
  p = strstr(s_a,sub_tp);
  i = 0;
  b->n_labels = 0;
  if(p)
    {
      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
      b->n_labels++;
      
      pos = 0;
      do 
	{
	  b->labels[b->n_labels-1][pos] = p[i+strlen(s_d)+1];
	  i++;
	  pos++;
	  if(p[i+strlen(s_d)+1] == '#') 
	    { 
	      b->labels[b->n_labels-1][pos] = '\0';
	      b->n_labels++;
	      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
	      i++;
	      pos=0;
	    }
	}
      while((p[i+strlen(s_d)+1] != ':') && 
	    (p[i+strlen(s_d)+1] != ',') && 
	    (p[i+strlen(s_d)+1] != '('));

      b->labels[b->n_labels-1][pos] = '\0';
    }

  if(p) 
    {
      if(b->n_labels == 1)
	printf("\n. Found label '%s' on edge %3d.",b->labels[0],b->num);
      else
	{
	  printf("\n. Found labels ");
	  For(i,b->n_labels) printf("'%s' ",b->labels[i]);
	  printf("on edge %3d.",b->num);
	}
    }

  Free(sub_tp);
}

/*********************************************************/

void Read_Branch_Length(char *s_d, char *s_a, arbre *tree)
{
  char *sub_tp;
  char *p;
  edge *b;
  int i;

  b = tree->t_edges[tree->num_curr_branch_available];

  sub_tp = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  For(i,b->n_labels) 
    {
      strcat(s_d,"#");
      strcat(s_d,b->labels[i]);
    }

  strcpy(sub_tp,s_d);
  strcat(sub_tp,":");
  p = strstr(s_a,sub_tp);
  if(p) 
    {
      b->l = atof((char *)p+(int)strlen(sub_tp)+1);
      tree->has_branch_lengths = 1;
    }
  Free(sub_tp);
}

/*********************************************************/

void Read_Node_Name(node *d, char *s_tree_d, arbre *tree)
{
  int i;

  if(!tree->t_edges[tree->num_curr_branch_available]->n_labels)
    {
      strcpy(d->name,s_tree_d);
    }
  else
    {
      i = 0;
      do
	{
	  d->name[i] = s_tree_d[i];
	  i++;
	}
      while(s_tree_d[i] != '#');
      d->name[i] = '\0';
    }
}
/*********************************************************/

void Unroot_Tree(char **subtrees)
{
  char **tmp_sub;
  int degree,i,j;

  printf("\n. Removing the root...\n");
  
  tmp_sub = Sub_Trees(subtrees[0],&degree);
  if(degree >= 2)
    {
      strcpy(subtrees[2],subtrees[1]);
      Clean_Multifurcation(tmp_sub,degree,2);
      For(j,2) strcpy(subtrees[j],tmp_sub[j]);
    }
  else
    {
      tmp_sub = Sub_Trees(subtrees[1],&degree);
      strcpy(subtrees[2],subtrees[0]);
      Clean_Multifurcation(tmp_sub,degree,2);
      For(j,2) strcpy(subtrees[j],tmp_sub[j]);
    }

  For(i,degree) Free(tmp_sub[i]);
  Free(tmp_sub);
}

/*********************************************************/

void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg)
{

  if(current_deg <= end_deg) return;
  else
    {
      char *s_tmp;
      int i;

      s_tmp = (char *)mCalloc(T_MAX_LINE,sizeof(char));

      strcat(s_tmp,"(\0");
      strcat(s_tmp,subtrees[0]);
      strcat(s_tmp,",\0");
      strcat(s_tmp,subtrees[1]);
      strcat(s_tmp,")\0");
      Free(subtrees[0]);
      subtrees[0] = s_tmp;

      for(i=1;i<current_deg-1;i++) strcpy(subtrees[i],subtrees[i+1]);

      Clean_Multifurcation(subtrees,current_deg-1,end_deg);
    }
}

/*********************************************************/

char **Sub_Trees(char *tree, int *degree)
{
  char **subs;
  int posbeg,posend;
  int i;

  if(tree[0] != '(') {*degree = 1; return NULL;}

  subs=(char **)mCalloc(NODE_DEG_MAX,sizeof(char *));

  For(i,NODE_DEG_MAX) subs[i]=(char *)mCalloc(strlen(tree)+1,sizeof(char));

  
  posbeg=posend=1;
  (*degree)=0;
  do
    {
      posbeg = posend;
      if(tree[posend] != '(')
	{
	  while((tree[posend] != ',' ) &&
		(tree[posend] != ':' ) &&
		(tree[posend] != '#' ) &&
		(tree[posend] != ')' )) 
	    {
	      posend++ ;
	    }
	  posend -= 1;
	}
      else posend=Next_Par(tree,posend);

      while((tree[posend+1] != ',') &&
	    (tree[posend+1] != ':') &&
	    (tree[posend+1] != '#') &&
	    (tree[posend+1] != ')')) {posend++;}


      strncpy(subs[(*degree)],tree+posbeg,posend-posbeg+1);
      strcat(subs[(*degree)],"\0");

      posend += 1;
      while((tree[posend] != ',') &&
	    (tree[posend] != ')')) {posend++;}
      posend+=1;


      (*degree)++;
      if((*degree) == NODE_DEG_MAX)
	{
	  For(i,(*degree))
	    printf("\n. Subtree %d : %s\n",i+1,subs[i]);

	  printf("\n. The degree of a node cannot be greater than %d\n",NODE_DEG_MAX);
	  Warn_And_Exit("\n");
	}
    }
  while(tree[posend-1] != ')');

  return subs;
}


/*********************************************************/

int Next_Par(char *s, int pos)
{
  int curr;

  curr=pos+1;

  while(*(s+curr) != ')')
    {
      if(*(s+curr) == '(') curr=Next_Par(s,curr);
      curr++;
    }

  return curr;
}

/*********************************************************/

void Print_Tree(FILE *fp, arbre *tree)
{
  char *s_tree;
  int i;

  s_tree = (char *)Write_Tree(tree);

  if(OUTPUT_TREE_FORMAT == 0)
    fprintf(fp,"%s\n",s_tree);
  else if(OUTPUT_TREE_FORMAT == 1)
    {
      fprintf(fp,"#NEXUS\n");
      fprintf(fp,"BEGIN TREES;\n");
      fprintf(fp,"\tTRANSLATE\n");
      For(i,tree->n_otu) fprintf(fp,"\t%3d\t%s,\n",i+1,tree->noeud[i]->name);
      fprintf(fp,"\tUTREE PAUP_1=\n");
      fprintf(fp,"%s\n",s_tree);
      fprintf(fp,"ENDBLOCK;");
    }
  Free(s_tree);
}

/*********************************************************/

char *Write_Tree(arbre *tree)
{

  char *s;
  int i;

  s=(char *)mCalloc(T_MAX_LINE,sizeof(char));

  s[0]='(';
  
  #ifdef PHYML 
  tree->n_root = NULL;
  tree->e_root = NULL;
  #endif

/*   if(!tree->n_root) */
/*     { */
      i = 0;
      while((!tree->noeud[tree->n_otu+i]->v[0]) ||
	    (!tree->noeud[tree->n_otu+i]->v[1]) ||
	    (!tree->noeud[tree->n_otu+i]->v[2])) i++;
      
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[0],s,tree);
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[1],s,tree);
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[2],s,tree);
/*     } */
/*   else */
/*     { */
/*       R_wtree(tree->n_root,tree->n_root->v[0],s,tree); */
/*       R_wtree(tree->n_root,tree->n_root->v[1],s,tree); */
/*     } */

  s[(int)strlen(s)-1]=')';
  s[(int)strlen(s)]=';';

  return s;
}

/*********************************************************/

void R_wtree(node *pere, node *fils, char *s_tree, arbre *tree)
{
  int i,p;

  p = -1;
  if(fils->tax)
    {
      if(OUTPUT_TREE_FORMAT == 0)
	strcat(s_tree,fils->name);
      else
	sprintf(s_tree+(int)strlen(s_tree),"%d",fils->num+1);

      if((fils->b[0]) && (fils->b[0]->l != -1))
	{
	  if(tree->print_labels)
	    {
	      if(fils->b[0]->n_labels < 10)
		For(i,fils->b[0]->n_labels) sprintf(s_tree+(int)strlen(s_tree),"#%s",fils->b[0]->labels[i]);
	      else
		sprintf(s_tree+(int)strlen(s_tree),"#%d_labels",fils->b[0]->n_labels);
	    }

	  strcat(s_tree,":");
/* 	  sprintf(s_tree+(int)strlen(s_tree),"%.10f",fils->b[0]->l); */
 	  sprintf(s_tree+(int)strlen(s_tree),"%f",fils->b[0]->l);
	}
      sprintf(s_tree+(int)strlen(s_tree),",");
   }
  else
    {
      s_tree[(int)strlen(s_tree)]='(';
      For(i,3)
	{
/* 	  if((fils->v[i] != pere) && (fils->b[i] != tree->e_root)) */
	  if(fils->v[i] != pere)
	    R_wtree(fils,fils->v[i],s_tree,tree);
	  else p=i;
	}
      s_tree[(int)strlen(s_tree)-1]=')';
      if(fils->b[0]->l != -1)
	{
	  if(tree->print_boot_val)
	    sprintf(s_tree+(int)strlen(s_tree),"%d",fils->b[p]->bip_score);
	  else if(tree->print_alrt_val)
	    sprintf(s_tree+(int)strlen(s_tree),"%f",fils->b[p]->ratio_test);

	  if(tree->print_labels)
	    {
	      if(fils->b[p]->n_labels < 10)
		For(i,fils->b[p]->n_labels) sprintf(s_tree+(int)strlen(s_tree),"#%s",fils->b[p]->labels[i]);
	      else
		sprintf(s_tree+(int)strlen(s_tree),"#%d_labels",fils->b[p]->n_labels);
	    }

	  strcat(s_tree,":");
	  sprintf(s_tree+(int)strlen(s_tree),"%f",fils->b[p]->l);

	  strcat(s_tree,",");
	}
    }
}

/*********************************************************/

void Init_Tree(arbre *tree, int n_otu)
{
  tree->n_otu                     = n_otu;
  tree->best_tree                 = NULL;
  tree->old_tree                  = NULL;
  tree->mat                       = NULL;
  tree->n_root                    = NULL;
  tree->e_root                    = NULL;
  tree->ps_tree                   = NULL;

  tree->depth_curr_path           = 0;
  tree->has_bip                   = 0;
  tree->n_moves                   = 0;
  tree->n_improvements            = 0;
  tree->number_of_lk_calls        = 0;
  tree->number_of_branch_lk_calls = 0;
  tree->bl_from_node_stamps       = 0;
  tree->lock_topo                 = 0;
  tree->ps_page_number            = 0;
  tree->init_lnL                  = UNLIKELY;
  tree->best_lnL                  = UNLIKELY;
  tree->c_lnL                     = UNLIKELY;
  tree->n_swap                    = 0;

  tree->n_pattern                 = -1;
  tree->prop_of_sites_to_consider = 1.;
  tree->n_root_pos                = -1.;
  tree->print_labels              = 1;

  tree->print_boot_val            = 0;
  tree->print_alrt_val            = 0;
  tree->num_curr_branch_available = 0;
}

/*********************************************************/

void Make_New_Edge_Label(edge *b)
{
  int i;

  b->labels = (char **)realloc(b->labels,(b->n_labels+BLOCK_LABELS)*sizeof(char *));

  if(!b->labels)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else
    {
      for(i=b->n_labels;i<b->n_labels+100;i++) b->labels[i] = (char *)mCalloc(T_MAX_LABEL,sizeof(char));
    }
}

/*********************************************************/

edge *Make_Edge_Light(node *a, node *d, int num)
{
  edge *b;

  b = (edge *)mCalloc(1,sizeof(edge));


  Init_Edge_Light(b,num);

  if(a && b)
    {
      b->left = a;  b->rght = d;
      if(a->tax) {b->rght = a; b->left = d;} /* root */
      /* a tip is necessary on the right side of the edge */

      (b->left == a)?
	(Make_Edge_Dirs(b,a,d)):
	(Make_Edge_Dirs(b,d,a));

      b->l                    = a->l[b->l_r];
      if(a->tax) b->l         = a->l[b->r_l];
      if(b->l < BL_MIN)  b->l = BL_MIN;
      else if(b->l > BL_MAX) b->l = BL_MAX;
      b->l_old                = b->l;
    }
  else
    {
      b->left = NULL;
      b->rght = NULL;
    }

  return b;

}

/*********************************************************/

void Init_Edge_Light(edge *b, int num)
{
  b->num                  = num;
  b->bip_score            = 0;
  b->dist_btw_edges       = .0;
  b->topo_dist_btw_edges  = 0;
  b->has_zero_br_len      = 0;
  b->is_p_lk_l_u2d        = 0;
  b->is_p_lk_r_u2d        = 0;

  b->p_lk_left            = NULL;
  b->p_lk_rght            = NULL;
  b->Pij_rr               = NULL;
}

/*********************************************************/

void Init_Node_Light(node *n, int num)
{
  n->list_of_reachable_tips = NULL;
  n->num                    = num;
  n->tax                    = -1;
  n->dist_to_root           = .0;
}

/*********************************************************/

void Make_Edge_Dirs(edge *b, node *a, node *d)
{
  int i;

  if(a == b->rght)
    {
      printf("\n. a->num = %3d ; d->num = %3d",a->num,d->num);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  if(d == b->left)
    {
      printf("\n. a->num = %3d ; d->num = %3d",a->num,d->num);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  b->l_r = b->r_l = -1;
  For(i,3)
    {
      if((a->v[i]) && (a->v[i] == d))
	{
	  b->l_r  = i; /* we consider here that 'a' is on the left handside of 'b'*/
	  a->b[i] = b;
	}
      if((d->v[i]) && (d->v[i] == a))
	{
	  b->r_l  = i; /* we consider here that 'd' is on the right handside of 'b'*/
	  d->b[i] = b;
	}
    }

  if(a->tax) {b->r_l = 0; For(i,3) if(d->v[i]==a) {b->l_r = i; break;}}


  b->l_v1 = b->l_v2 = b->r_v1 = b->r_v2 = -1;
  For(i,3)
    {
      if(b->left->v[i] != b->rght)
	{
	  if(b->l_v1 < 0) b->l_v1 = i;
	  else            b->l_v2 = i;
	}

      if(b->rght->v[i] != b->left)
	{
	  if(b->r_v1 < 0) b->r_v1 = i;
	  else            b->r_v2 = i;
	}
    }
}

/*********************************************************/

void Make_Edge_Pars(edge *b, arbre *tree)
{
  int site;

  b->pars_l = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->pars_r = (int *)mCalloc(tree->data->crunch_len,sizeof(int));

  b->ui_l = (unsigned int *)mCalloc(tree->data->crunch_len,sizeof(unsigned int));
  b->ui_r = (unsigned int *)mCalloc(tree->data->crunch_len,sizeof(unsigned int));

  b->p_pars_l = (int **)mCalloc(tree->data->crunch_len,sizeof(int *));
  b->p_pars_r = (int **)mCalloc(tree->data->crunch_len,sizeof(int *));

  For(site,tree->data->crunch_len)
    {
      b->p_pars_l[site] = (int *)mCalloc(tree->mod->ns,sizeof(int));
      b->p_pars_r[site] = (int *)mCalloc(tree->mod->ns,sizeof(int));
    }
}

/*********************************************************/

void Make_Edge_Lk(edge *b, arbre *tree)
{
  int i,j,k;

  b->l_old = b->l;

  b->div_post_pred_left = (short int *)mCalloc((tree->mod->datatype == NT)?(4):(20),sizeof(short int));
  b->div_post_pred_rght = (short int *)mCalloc((tree->mod->datatype == NT)?(4):(20),sizeof(short int));

  b->Pij_rr   = (double ***)mCalloc(tree->mod->n_catg,sizeof(double **));
  For(i,tree->mod->n_catg)
    {
      b->Pij_rr[i] = (double **)mCalloc(tree->mod->ns,sizeof(double *));
      For(j,tree->mod->ns) b->Pij_rr[i][j] = (double *)mCalloc(tree->mod->ns,sizeof(double ));
    }
  
  b->scale_left = b->scale_rght = 0;
  
  if(!b->left->tax)
    b->sum_scale_f_left = (phydbl *)mCalloc(tree->data->crunch_len,sizeof(phydbl ));
  else
    b->sum_scale_f_left = NULL;
  
  if(!b->rght->tax)
    b->sum_scale_f_rght = (phydbl *)mCalloc(tree->data->crunch_len,sizeof(phydbl ));
  else
    b->sum_scale_f_rght = NULL;
  
  
  if((!b->left->tax) || (tree->mod->s_opt->greedy))
    {
      b->p_lk_left = (phydbl ***)mCalloc(tree->data->crunch_len,sizeof(phydbl **));
      For(j,tree->data->crunch_len)
	{
	  b->p_lk_left[j] = (phydbl **)mCalloc(tree->mod->n_catg,sizeof(phydbl *));	  
	  For(k,tree->mod->n_catg) b->p_lk_left[j][k] = (phydbl *)mCalloc(tree->mod->ns,sizeof(phydbl ));
	}
      b->p_lk_tip_l = NULL;
    }
  else if(b->left->tax)
    {
      b->p_lk_left   = NULL;
      
      b->p_lk_tip_l  = (short int **)mCalloc(tree->data->crunch_len,sizeof(short int *));
      For(j,tree->data->crunch_len) b->p_lk_tip_l[j] = (short int *)mCalloc(tree->mod->ns,sizeof(short int));
    }
  
  
  if((!b->rght->tax) || (tree->mod->s_opt->greedy))
    {
      b->p_lk_rght = (phydbl ***)mCalloc(tree->data->crunch_len,sizeof(phydbl **));
      
      For(j,tree->data->crunch_len)
	{
	  b->p_lk_rght[j] = (phydbl **)mCalloc(tree->mod->n_catg,sizeof(phydbl *));
	  
	  For(k,tree->mod->n_catg) b->p_lk_rght[j][k] = (phydbl *)mCalloc(tree->mod->ns,sizeof(phydbl ));
	}
      b->p_lk_tip_r = NULL;
    }
  else if(b->rght->tax)
    {
      b->p_lk_rght = NULL;
      
      b->p_lk_tip_r  = (short int **)mCalloc(tree->data->crunch_len,sizeof(short int *));
      For(j,tree->data->crunch_len) b->p_lk_tip_r[j] = (short int *)mCalloc(tree->mod->ns,sizeof(short int));
    }
}

/*********************************************************/

void Make_Edge_NNI(edge *b)
{
  b->nni    = Make_NNI();
  b->nni->b = b;
  b->nni->left = b->left;
  b->nni->rght = b->rght;
}

/*********************************************************/

nni *Make_NNI()
{
  nni *a_nni;
  a_nni = (nni *)mCalloc(1,sizeof(nni ));
  Init_NNI(a_nni);
  return a_nni;
}

/*********************************************************/

void Init_NNI(nni *a_nni)
{
  a_nni->left         = NULL;
  a_nni->rght         = NULL;
  a_nni->b            = NULL;
  a_nni->init_l       = -1.;
  a_nni->init_lk      = .0;
  a_nni->score        = +1.0;
  a_nni->best_l       = -1.;
  a_nni->swap_node_v1 = NULL;
  a_nni->swap_node_v2 = NULL;
  a_nni->swap_node_v3 = NULL;
  a_nni->swap_node_v4 = NULL;
  a_nni->lk0          = UNLIKELY;
  a_nni->lk1          = UNLIKELY;
  a_nni->lk2          = UNLIKELY;
  a_nni->l0           = -1.0;
  a_nni->l1           = -1.0;
  a_nni->l2           = -1.0;
}

/*********************************************************/

node *Make_Node_Light(int num)
{
  node *n;
  n        = (node *)mCalloc(1,sizeof(node));
  n->v     = (node **)mCalloc(3,sizeof(node *));
  n->l     = (phydbl *)mCalloc(3,sizeof(phydbl));
  n->b     = (edge **)mCalloc(3,sizeof(edge *));
  n->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  n->score = (phydbl *)mCalloc(3,sizeof(phydbl));
  Init_Node_Light(n,num);
  return n;
}

/*********************************************************/


void Make_Node_Lk(node *n)
{
/*   n->n_ex_nodes = (int *)mCalloc(2,sizeof(int)); */
  return;
}

/*********************************************************/

seq **Get_Seq(option *io,  int rw)
{
  seq **data;
  int i,j;
  char **buff;
  int n_unkn,n_removed,pos;
  int *remove;


/*   rewind(fp_seq); */

  if(io->interleaved) data = Read_Seq_Interleaved(io->fp_in_seq,&(io->mod->n_otu));
  else                data = Read_Seq_Sequential(io->fp_in_seq,&(io->mod->n_otu));

  if(data)
    {
      buff = (char **)mCalloc(io->mod->n_otu,sizeof(char *));
      For(i,io->mod->n_otu) buff[i] = (char *)mCalloc(data[0]->len,sizeof(char));
      remove = (int *)mCalloc(data[0]->len,sizeof(int));

      n_removed = 0;

      For(i,data[0]->len)
	{
	  For(j,io->mod->n_otu)
	    {
	      if((data[j]->state[i] == '?') || (data[j]->state[i] == '-')) data[j]->state[i] = 'X';
	      if((io->mod->datatype == NT) && (data[j]->state[i] == 'N')) data[j]->state[i] = 'X';
	      if(data[j]->state[i] == 'U') data[j]->state[i] = 'T';
	    }

	  n_unkn = 0;
	  For(j,io->mod->n_otu) if(data[j]->state[i] == 'X') n_unkn++;

	  if(n_unkn == io->mod->n_otu)
	    {
	      remove[i] = 1;
	      n_removed++;
	    }

	  For(j,io->mod->n_otu) buff[j][i] = data[j]->state[i];
	}

      if(n_removed > 0)
	{
	  if(io->mod->datatype == NT)
	    printf("\n. %d sites are made from completely undetermined states ('X', '-', '?' or 'N')...\n",n_removed);
	  else
	    printf("\n. %d sites are made from completely undetermined states ('X', '-', '?')...\n",n_removed);
	}

      pos = 0;
      For(i,data[0]->len)
	{
/* 	  if(!remove[i]) */
/* 	    { */
	      For(j,io->mod->n_otu) data[j]->state[pos] = buff[j][i];
	      pos++;
/* 	    } */
	}

      For(i,io->mod->n_otu) data[i]->len = pos;
      For(i,io->mod->n_otu) Free(buff[i]);
      Free(buff);
      Free(remove);
    }
  return data;
}

/*********************************************************/

seq **Read_Seq_Sequential(FILE *in, int *n_otu)
{
  int i;
  char *line;
  int len,readok;
  seq **data;
  char c;
  char *format = (char *)mCalloc(T_MAX_NAME, sizeof(char));

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  readok = len = 0;
  do
    {
      if(fscanf(in,"%s",line) == EOF)
	{
	  Free(line); return NULL;
	}
      else
	{
	  if(strcmp(line,"\n") && strcmp(line,"\n") && strcmp(line,"\t"))
	    {
	      *n_otu = atoi(line);
	      data = (seq **)mCalloc(*n_otu,sizeof(seq *));
	      if(*n_otu <= 0) Warn_And_Exit("\n. Problem with sequence format\n");
	      fscanf(in,"%s",line);
	      len = atoi(line);
	      if(len <= 0) Warn_And_Exit("\n. Problem with sequence format\n");
	      else readok = 1;
	    }
	}
    }while(!readok);


/*   while((c=fgetc(in))!='\n'); */
  while(((c=fgetc(in))!='\n') && (c != ' ') && (c != '\r') && (c != '\t'));

  For(i,*n_otu)
    {
      data[i] = (seq *)mCalloc(1,sizeof(seq));
      data[i]->len = 0;
      data[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      data[i]->state = (char *)mCalloc(T_MAX_SEQ,sizeof(char));
      data[i]->is_ambigu = NULL;
      sprintf(format, "%%%ds", T_MAX_NAME);
      fscanf(in, format, data[i]->name);

      while(data[i]->len < len)
	Read_One_Line_Seq(&data,i,in);

      if(data[i]->len != len)
	{
	  printf("\n. Err: Problem with species %s's sequence (check the format)\n",
		 data[i]->name);
	  Warn_And_Exit("");
	}
    }

  /*   fgets(line,T_MAX_LINE,in);  */
  /* inter data sets */

  Free(format);
  Free(line);
  return data;
}

/*********************************************************/

seq **Read_Seq_Interleaved(FILE *in, int *n_otu)
{
  int i,end,num_block;
  char *line;
  int len,readok;
  seq **data;
  char c;
  char *format;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  format = (char *)mCalloc(T_MAX_NAME, sizeof(char));

  readok = len = 0;
  do
    {
      if(fscanf(in,"%s",line) == EOF)
	{
	  Free(format);
	  Free(line); return NULL;
	}
      else
	{
	  if(strcmp(line,"\n") && strcmp(line,"\r") && strcmp(line,"\t"))
	    {
	      *n_otu = atoi(line);
	      data = (seq **)mCalloc(*n_otu,sizeof(seq *));
	      if(*n_otu <= 0) Warn_And_Exit("\n. Problem with sequence format\n");
	      fscanf(in,"%s",line);
	      len = atoi(line);
	      if(len <= 0) Warn_And_Exit("\n. Problem with sequence format\n");
	      else readok = 1;
	    }
	}
    }while(!readok);


  while(((c=fgetc(in))!='\n') && (c != ' ') && (c != '\r') && (c != '\t'));

  end = 0;
  For(i,*n_otu)
    {
      data[i] = (seq *)mCalloc(1,sizeof(seq));
      data[i]->len = 0;
      data[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      data[i]->state = (char *)mCalloc(T_MAX_SEQ,sizeof(char));
      data[i]->is_ambigu = NULL;
      sprintf(format, "%%%ds", T_MAX_NAME);
      fscanf(in, format, data[i]->name);
      if(!Read_One_Line_Seq(&data,i,in))
	{
	  end = 1;
	  if((i != *n_otu) && (i != *n_otu-1))
	    {
	      printf("\n. Err: Problem with species %s's sequence\n",data[i]->name);
	      Warn_And_Exit("");
	    }
	  break;
	}
    }

  if(data[0]->len == len) end = 1;

  if(!end)
    {
      end = 0;

      num_block = 1;
      do
	{
	  num_block++;

	  /* interblock */
	  if(!fgets(line,T_MAX_LINE,in)) break;

	  if(line[0] != 13 && line[0] != 10)
	    {
	      printf("\n. One or more missing sequences in block %d\n",num_block-1);
	      Warn_And_Exit("");
	    }
	  
	  For(i,*n_otu)
	    if(data[i]->len != len)
	      break;
	  
	  if(i == *n_otu) break;
	  
	  
	  For(i,*n_otu)
	    {
	      if(data[i]->len > len)
		{
		  printf("\n. Observed length=%d expected length=%d\n",data[i]->len,len);
		  printf("\n. Err: Problem with species %s's sequence\n",data[i]->name);
		  Warn_And_Exit("");
		}
	      else if(!Read_One_Line_Seq(&data,i,in))
		{
		  end = 1;
		  if((i != *n_otu) && (i != *n_otu-1))
		    {
		      printf("\n. Err: Problem with species %s's sequence\n",data[i]->name);
		      Warn_And_Exit("");
		    }
		  break;
		}
	    }
	}while(!end);
    }

  For(i,*n_otu)
    {
      if(data[i]->len != len)
	{
	  printf("\n. Check sequence '%s' length...\n",data[i]->name);
	  Warn_And_Exit("");
	}
    }

  Free(format);
  Free(line);
  return data;
}

/*********************************************************/

int Read_One_Line_Seq(seq ***data, int num_otu, FILE *in)
{
  char c;

  c=' ';
  while(1)
    {
/*       if((c == EOF) || (c == '\n') || (c == '\r')) break; */
        if((c == EOF) || (c == 13) || (c == 10)) break;
      else if((c==' ') || (c=='\t')) {c=(char)fgetc(in); continue;}
      Uppercase(&c);

      if (strchr("ABCDEFGHIKLMNOPQRSTUVWXYZ?-.", c) == NULL)
	{
	  printf("\n. Err: bad symbol: \"%c\" at position %d of species %s\n",
		 c,(*data)[num_otu]->len,(*data)[num_otu]->name);
	  Warn_And_Exit("");
	}

      if(c == '.')
	{
	  c = (*data)[0]->state[(*data)[num_otu]->len];
	  if(!num_otu)
	    Warn_And_Exit("\n. Err: Symbol \".\" should not appear in the first sequence\n");
	}
      (*data)[num_otu]->state[(*data)[num_otu]->len]=c;
      (*data)[num_otu]->len++;
      c = (char)fgetc(in);
    }
  if(c == EOF) return 0;
  else return 1;
}

/*********************************************************/

void Uppercase(char *ch)
{
  /* convert ch to upper case -- either ASCII or EBCDIC */
   *ch = isupper((int)*ch) ? *ch : toupper((int)*ch);
}

/*********************************************************/

allseq *Compact_Seq(seq **data, option *io)
{
  allseq *alldata_tmp,*alldata;
  int i,j,k,site;
  int n_patt,which_patt,n_invar;
  char **sp_names;
  int n_otu, n_sites;
  pnode *proot;
  int compress;
  int n_ambigu,is_ambigu;

  n_otu        = io->mod->n_otu;
  n_patt       = 0;
  which_patt   = 0;

  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  For(i,n_otu)
    {
      sp_names[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(sp_names[i],data[i]->name);
    }

  alldata_tmp = Make_Cseq(n_otu,data[0]->len,data[0]->len,sp_names);
  proot       = (pnode *)Create_Pnode(T_MAX_ALPHABET);
 
  For(i,n_otu) Free(sp_names[i]);
  Free(sp_names);


  if(data[0]->len%io->mod->stepsize)
    {
      printf("\n. Sequence length is not a multiple of %d\n",io->mod->stepsize);
      Warn_And_Exit("");
    }
  
  compress = 1;
/*   compress = 0; */
  n_ambigu = 0;
  is_ambigu = 0;

  Fors(site,data[0]->len,io->mod->stepsize)
    {
      if(io->rm_ambigu)
	{
	  is_ambigu = 0;
	  For(j,n_otu)
	    {
	      if(Is_Ambigu(data[j]->state+site,io->mod->datatype,io->mod->stepsize))
		{
		  break;
		}
	    }
	  if(j != n_otu)
	    {
	      is_ambigu = 1;
	      n_ambigu++;
	    }
	}

      if(!is_ambigu)
	{
	  if(compress)
	    {
	      which_patt = -1;
	      Traverse_Prefix_Tree(site,-1,&which_patt,&n_patt,data,io,proot);
	      if(which_patt == n_patt-1) /* New pattern found */
		{
		  n_patt--;
		  k=n_patt;
		}
	      else
		{
		  k = n_patt-10;
		}
	    }
	  else
	    {
	      printf("\n. WARNING: sequences are not compressed !");
	      k = n_patt;
	    }
	  
	  if(k == n_patt) /* add a new site pattern */
	    {
	      For(j,n_otu)
		Copy_One_State(data[j]->state+site,
			       alldata_tmp->c_seq[j]->state+n_patt,
			       io->mod->stepsize);
	      
	      
	      For(i,n_otu)
		{
		  For(j,n_otu)
		    {
		      if(!(Are_Compatible(alldata_tmp->c_seq[i]->state+n_patt,
					  alldata_tmp->c_seq[j]->state+n_patt,
					  io->mod->stepsize,
					  io->mod->datatype))) break;
		    }
		  if(j != n_otu) break;
		}
	      
	      if((j == n_otu) && (i == n_otu)) /* all characters at that site are compatible -> the site is invariant */
		{
		  For(j,n_otu)
		    {
		      alldata_tmp->invar[n_patt] = Assign_State(alldata_tmp->c_seq[j]->state+n_patt,
								io->mod->datatype,
								io->mod->stepsize);
		      if(alldata_tmp->invar[n_patt] > -1.) break;
		    }
		}
	      else alldata_tmp->invar[n_patt] = -1;
	      
	      alldata_tmp->sitepatt[site] = n_patt;
	      alldata_tmp->wght[n_patt]  += 1;
	      n_patt                     += io->mod->stepsize;
	    }
	  else
	    {
	      alldata_tmp->sitepatt[site]    = which_patt;
	      alldata_tmp->wght[which_patt] += 1;
	    }
	}
    }
  
  data[0]->len -= n_ambigu;
  
  alldata_tmp->init_len                   = data[0]->len;
  alldata_tmp->crunch_len                 = n_patt;
  For(i,n_otu) alldata_tmp->c_seq[i]->len = n_patt;

  /* wash. */
  /*  printf("\n. %d patterns found. (out of a total of %d sites) \n",n_patt,data[0]->len);*/

  if((io->rm_ambigu) && (n_ambigu))
    {
      printf("\n. Removed %d columns of the alignment as the contain ambiguous characters (e.g., gaps) \n",n_ambigu);
    }

  n_invar=0;
  For(i,alldata_tmp->crunch_len) if(alldata_tmp->invar[i] > -1.) n_invar+=(int)alldata_tmp->wght[i];

  /* wash. */
  /*  printf("\n. %d sites without polymorphism (%.2f%c).\n",n_invar,100.*(phydbl)n_invar/data[0]->len,'%');*/

  alldata_tmp->obs_pinvar = (phydbl)n_invar/data[0]->len;

  n_sites = 0;
  For(i,alldata_tmp->crunch_len) n_sites += alldata_tmp->wght[i];
  if(n_sites != data[0]->len)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(io->mod->datatype == NT) Get_Base_Freqs(alldata_tmp);
  else                        Get_AA_Freqs(alldata_tmp);

  alldata = Copy_Cseq(alldata_tmp, alldata_tmp->crunch_len, io->mod->ns);

  Free_Cseq(alldata_tmp);
  Free_Prefix_Tree(proot,T_MAX_ALPHABET);

  return alldata;
}

/*********************************************************/

allseq *Compact_CSeq(allseq *data, model *mod)
{
  allseq *alldata;
  int i,j,k,site;
  int n_patt,which_patt;
  int n_otu;

  n_otu = data->n_otu;

  alldata         = (allseq *)mCalloc(1,sizeof(allseq));
  alldata->n_otu  = n_otu;
  alldata->c_seq  = (seq **)mCalloc(n_otu,sizeof(seq *));
  alldata->wght   = (int *)mCalloc(data->crunch_len,sizeof(int));
  alldata->b_frq  = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  alldata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  alldata->invar  = (short int *)mCalloc(data->crunch_len,sizeof(short int));

  alldata->crunch_len = alldata->init_len = -1;
  For(j,n_otu)
    {
      alldata->c_seq[j]            = (seq *)mCalloc(1,sizeof(seq));
      alldata->c_seq[j]->name      = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(alldata->c_seq[j]->name,data->c_seq[j]->name);
      alldata->c_seq[j]->state     = (char *)mCalloc(data->crunch_len,sizeof(char));
      alldata->c_seq[j]->is_ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
      alldata->c_seq[j]->state[0]  = data->c_seq[j]->state[0];
    }

  n_patt = which_patt =  0;

  Fors(site,data->crunch_len,mod->stepsize)
    {
      if(data->wght[site])
	{
	  Fors(k,n_patt,mod->stepsize)
	    {
	      For(j,n_otu)
		{
		  if(strncmp(alldata->c_seq[j]->state+k,
			     data->c_seq[j]->state+site,
			     mod->stepsize))
		    break;
		}
	      
	      if(j == n_otu)
		{
		  which_patt = k;
		  break;
		}
	    }
	  
	  /*       /\* TO DO *\/ */
	  /*       k = n_patt; */
	  
	  if(k == n_patt)
	    {
	      For(j,n_otu) Copy_One_State(data->c_seq[j]->state+site,
					  alldata->c_seq[j]->state+n_patt,
					  mod->stepsize);
	      
	      For(i,n_otu)
		{
		  For(j,n_otu)
		    {
		      if(!(Are_Compatible(alldata->c_seq[i]->state+n_patt,
					  alldata->c_seq[j]->state+n_patt,
					  mod->stepsize,
					  mod->datatype))) break;
		    }
		  if(j != n_otu) break;
		}
	      
	      if((j == n_otu) && (i == n_otu)) 
		{
		  For(j,n_otu)
		    {
		      alldata->invar[n_patt] = Assign_State(alldata->c_seq[j]->state+n_patt,
							    mod->datatype,
							    mod->stepsize);
		      if(alldata->invar[n_patt] > -1.) break;
		    }
		}
	      else alldata->invar[n_patt] = -1;
	      
	      alldata->wght[n_patt] += data->wght[site];
	      n_patt+=mod->stepsize;
	    }
	  else alldata->wght[which_patt] += data->wght[site];
	  
	  /*       Print_Site(alldata,k,n_otu,"\n",mod->stepsize); */
	}
    }
  
  alldata->init_len   = data->crunch_len;
  alldata->crunch_len = n_patt;
  For(i,n_otu) alldata->c_seq[i]->len = n_patt;

  (mod->datatype == NT)?
    (Get_Base_Freqs(alldata)):
    (Get_AA_Freqs(alldata));

  return alldata;
}

/*********************************************************/

void Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt, seq **data, option *io, pnode *n)
{
  int ret_val;

  ret_val = -1;

  if(seqnum == io->mod->n_otu-1)
    {
      n->weight++;
      if(n->weight == 1)
	{
	  n->num = *n_patt;
	  (*n_patt) += 1;
	}
      (*patt_num) = n->num;
      return;
    }
  else
    {
      int next_state;

      next_state = -1;
      next_state = Assign_State_With_Ambiguity(data[seqnum+1]->state+site,
					       io->mod->datatype,
					       io->mod->stepsize);

      if(!n->next[next_state])
	{
	  n->next[next_state] = Create_Pnode(T_MAX_ALPHABET);
	}
      Traverse_Prefix_Tree(site,seqnum+1,patt_num,n_patt,data,io,n->next[next_state]);
    }
}

/*********************************************************/

pnode *Create_Pnode(int size)
{
  pnode *n;
  int i;

  n = (pnode *)mCalloc(1,sizeof(pnode ));
  n->next = (pnode **)mCalloc(size,sizeof(pnode *));
  For(i,size) n->next[i] = NULL;
  n->weight = 0;
  n->num = -1;
  return n;
}
/*********************************************************/
/*********************************************************/

void Get_Base_Freqs(allseq *data)
{
  int i,j,k;
  phydbl A,C,G,T;
  phydbl fA,fC,fG,fT;
  int w;

  fA = fC = fG = fT = .25;

  For(k,8)
    {
      A = C = G = T = .0;
      For(i,data->n_otu)
	{
	  For(j,data->crunch_len)
	    {
	      w = data->wght[j];
	      if(w)
		{
		  switch(data->c_seq[i]->state[j])
		    {
		    case 'A' : A+=w;
		      break;
		    case 'C' : C+=w;
		      break;
		    case 'G' : G+=w;
		      break;
		    case 'T' : T+=w;
		      break;
		    case 'U' : T+=w;
		      break;
		    case 'M' : C+=w*fC/(fC+fA); A+=w*fA/(fA+fC);
		      break;
		    case 'R' : G+=w*fG/(fA+fG); A+=w*fA/(fA+fG);
		      break;
		    case 'W' : T+=w*fT/(fA+fT); A+=w*fA/(fA+fT);
		      break;
		    case 'S' : C+=w*fC/(fC+fG); G+=w*fG/(fC+fG);
		      break;
		    case 'Y' : C+=w*fC/(fC+fT); T+=w*fT/(fT+fC);
		      break;
		    case 'K' : G+=w*fG/(fG+fT); T+=w*fT/(fT+fG);
		      break;
		    case 'B' : C+=w*fC/(fC+fG+fT); G+=w*fG/(fC+fG+fT); T+=w*fT/(fC+fG+fT);
		      break;
		    case 'D' : A+=w*fA/(fA+fG+fT); G+=w*fG/(fA+fG+fT); T+=w*fT/(fA+fG+fT);
		      break;
		    case 'H' : A+=w*fA/(fA+fC+fT); C+=w*fC/(fA+fC+fT); T+=w*fT/(fA+fC+fT);
		      break;
		    case 'V' : A+=w*fA/(fA+fC+fG); C+=w*fC/(fA+fC+fG); G+=w*fG/(fA+fC+fG);
		      break;
		    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
		      A+=w*fA; C+=w*fC; G+=w*fG; T+=w*fT; break;
		    default : break;
		    }
		}
	    }
	}
      fA = A/(A+C+G+T);
      fC = C/(A+C+G+T);
      fG = G/(A+C+G+T);
      fT = T/(A+C+G+T);
    }
  
  data->b_frq[0] = fA;
  data->b_frq[1] = fC;
  data->b_frq[2] = fG;
  data->b_frq[3] = fT;
}

/*********************************************************/

void Get_AA_Freqs(allseq *data)
{
  int i,j,k;
  phydbl A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y;
  phydbl fA,fC,fD,fE,fF,fG,fH,fI,fK,fL,fM,fN,fP,fQ,fR,fS,fT,fV,fW,fY;
  int w;
  phydbl sum;
  
  fA = fC = fD = fE = fF = fG = fH = fI = fK = fL =
  fM = fN = fP = fQ = fR = fS = fT = fV = fW = fY = 1./20.;
  
  For(k,8)
    {
      A = C = D = E = F = G = H = I = K = L =
      M = N = P = Q = R = S = T = V = W = Y = .0;

      For(i,data->n_otu)
	{
	  For(j,data->crunch_len)
	    {
	      w = data->wght[j];
	      if(w)
		{
		  switch(data->c_seq[i]->state[j])
		    {
		    case 'A' : A+=w;		break;
		    case 'C' : C+=w;		break;
		    case 'D' : D+=w;		break;
		    case 'E' : E+=w;		break;
		    case 'F' : F+=w;		break;
		    case 'G' : G+=w;		break;
		    case 'H' : H+=w;		break;
		    case 'I' : I+=w;		break;
		    case 'K' : K+=w;		break;
		    case 'L' : L+=w;		break;
		    case 'M' : M+=w;		break;
		    case 'N' : N+=w;		break;
		    case 'P' : P+=w;		break;
		    case 'Q' : Q+=w;		break;
		    case 'R' : R+=w;		break;
		    case 'S' : S+=w;		break;
		    case 'T' : T+=w;		break;
		    case 'V' : V+=w;		break;
		    case 'W' : W+=w;		break;
		    case 'Y' : Y+=w;		break;
		    case 'Z' : Q+=w;		break;
		    case 'X' : case '?' : case 'O' : case '-' :
		      A+=w*fA;
		      C+=w*fC;
		      D+=w*fD;
		      E+=w*fE;
		      F+=w*fF;
		      G+=w*fG;
		      H+=w*fH;
		      I+=w*fI;
		      K+=w*fK;
		      L+=w*fL;
		      M+=w*fM;
		      N+=w*fN;
		      P+=w*fP;
		      Q+=w*fQ;
		      R+=w*fR;
		      S+=w*fS;
		      T+=w*fT;
		      V+=w*fV;
		      W+=w*fW;
		      Y+=w*fY;
		      break;
		    default : break;
		    }
		}
	    }
	}
      sum = (A+C+D+E+F+G+H+I+K+L+M+N+P+Q+R+S+T+V+W+Y);
      fA = A/sum;      fC = C/sum;      fD = D/sum;      fE = E/sum;
      fF = F/sum;      fG = G/sum;      fH = H/sum;      fI = I/sum;
      fK = K/sum;      fL = L/sum;      fM = M/sum;      fN = N/sum;
      fP = P/sum;      fQ = Q/sum;      fR = R/sum;      fS = S/sum;
      fT = T/sum;      fV = V/sum;      fW = W/sum;      fY = Y/sum;
    }

  data->b_frq[0]  = fA;  data->b_frq[1]  = fR;  data->b_frq[2]  = fN;  data->b_frq[3]  = fD;
  data->b_frq[4]  = fC;  data->b_frq[5]  = fQ;  data->b_frq[6]  = fE;  data->b_frq[7]  = fG;
  data->b_frq[8]  = fH;  data->b_frq[9]  = fI;  data->b_frq[10] = fL;  data->b_frq[11] = fK;
  data->b_frq[12] = fM;  data->b_frq[13] = fF;  data->b_frq[14] = fP;  data->b_frq[15] = fS;
  data->b_frq[16] = fT;  data->b_frq[17] = fW;  data->b_frq[18] = fY;  data->b_frq[19] = fV;
}

/*********************************************************/

arbre *Read_Tree_File(FILE *fp_input_tree)
{
  char *line;
  arbre *tree;
  int i;
  char c;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  do
    c=fgetc(fp_input_tree);
  while((c != '(') && (c != EOF));

  if(c==EOF)
      {
          Free(line);
          return NULL;
      }

  i=0;
  for(;;)
    {
      if((c == ' ') || (c == '\n'))
	{
	  c=fgetc(fp_input_tree);
	  if(c==EOF) break;
	  else continue;
	}

      line[i]=c;
      i++;
      c=fgetc(fp_input_tree);
      if(c==EOF || c==';') break;
    }

  tree = Read_Tree(line);
  Free(line);
  return tree;
}

/*********************************************************/

void Connect_Edges_To_Nodes_Recur(node *a, node *d, arbre *tree)
{
  int i;

  Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
  tree->num_curr_branch_available += 1;

  if(d->tax) return;
  else For(i,3) if(d->v[i] != a) Connect_Edges_To_Nodes_Recur(d,d->v[i],tree);
}

/*********************************************************/

void Connect_One_Edge_To_Two_Nodes(node *a, node *d, edge *b, arbre *tree)
{
  int i,dir_a_d;

  dir_a_d = -1;
  For(i,3) if(a->v[i] == d) {dir_a_d = i; break;}


  a->b[dir_a_d] = b;
  b->num        = tree->num_curr_branch_available;
  b->left       = a;
  b->rght       = d;
  if(a->tax) {b->rght = a; b->left = d;} /* root */
  /* a tip is necessary on the right side of the edge */

  (b->left == a)?
    (Make_Edge_Dirs(b,a,d)):
    (Make_Edge_Dirs(b,d,a));

  b->l                    = a->l[b->l_r];
  if(a->tax) b->l         = a->l[b->r_l];
  if(b->l < BL_MIN)  b->l = BL_MIN;
  else if(b->l > BL_MAX) b->l = BL_MAX;
  b->l_old                = b->l;
}

/*********************************************************/

void Update_Dirs(arbre *tree)
{
  int i;
  int buff;
  edge *b;

  b = NULL;
  buff = -1;
  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      
      if((!b->left->tax) && (b->left->v[b->l_v1]->num < b->left->v[b->l_v2]->num))
	{
	  buff    = b->l_v1;
	  b->l_v1 = b->l_v2;
	  b->l_v2 = buff;
	}
      if((!b->rght->tax) && (b->rght->v[b->r_v1]->num < b->rght->v[b->r_v2]->num))
	{
	  buff    = b->r_v1;
	  b->r_v1 = b->r_v2;
	  b->r_v2 = buff;
	}
    }

}

/*********************************************************/

void Exit(char *message)
{
  fflush(NULL);
  fprintf(stderr,"%s",message);
  exit(1);
}

/*********************************************************/

void *mCalloc(int nb, size_t size)
{
  void *allocated;

  if((allocated = calloc((size_t)nb,(size_t)size)) != NULL)
    {
      return allocated;
    }
  else
    Warn_And_Exit("\n. Err: low memory\n");

  return NULL;
}

/*********************************************************/

void *mRealloc(void *p,int nb, size_t size)
{
  if((p = realloc(p,(size_t)nb*size)) != NULL)
	return p;
  else
    Warn_And_Exit("\n. Err: low memory\n");

  return NULL;
}

/*********************************************************/

/* arbre *Make_Light_Tree_Struct(int n_otu) */
/* { */
/*   arbre *tree; */
/*   int i; */

/*   tree          = (arbre *)mCalloc(1,sizeof(arbre )); */
/*   tree->t_edges = (edge **)mCalloc(2*n_otu-3,sizeof(edge *)); */
/*   tree->noeud   = (node **)mCalloc(2*n_otu-2,sizeof(node *)); */
/*   tree->n_otu   = n_otu; */

/*   For(i,2*n_otu-3) */
/*     tree->t_edges[i] = Make_Edge_Light(NULL,NULL,i); */

/*   For(i,2*n_otu-2) */
/*     tree->noeud[i] = Make_Node_Light(i); */

/*   return tree; */
/* } */

/*********************************************************/

int Sort_Phydbl_Decrease(const void *a, const void *b)
{
    if((*(phydbl *)(a)) >= (*(phydbl *)(b))) return -1;
    else return 1;
}

/*********************************************************/

void Qksort(phydbl* A, int ilo, int ihi)
{
    phydbl pivot;	// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    phydbl tempEntry;	// temporary entry used for swapping

    if (ilo >= ihi) {
	return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
      if (A[uhi] > pivot) {
	    // Here, we can reduce the size of the unpartitioned region and
	    // try again.
	    uhi--;
	} else {
	    // Here, A[uhi] <= pivot, so swap entries at indices ulo and
	    // uhi.
	    tempEntry = A[ulo];
	    A[ulo] = A[uhi];
	    A[uhi] = tempEntry;
	    // After the swap, A[ulo] <= pivot.
	    if (A[ulo] < pivot) {
		// Swap entries at indices ieq and ulo.
		tempEntry = A[ieq];
		A[ieq] = A[ulo];
		A[ulo] = tempEntry;
		// After the swap, A[ieq] < pivot, so we need to change
		// ieq.
		ieq++;
		// We also need to change ulo, but we also need to do
		// that when A[ulo] = pivot, so we do it after this if
		// statement.
	    }
	    // Once again, we can reduce the size of the unpartitioned
	    // region and try again.
	    ulo++;
	}
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    Qksort(A, ilo, ieq - 1);
    Qksort(A, uhi + 1, ihi);
}

/********************************************************/

void Qksort_Matrix(phydbl **A, int col, int ilo, int ihi)
{
    phydbl pivot;	// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    phydbl *tempEntry;	// temporary entry used for swapping

    tempEntry = NULL;

    if (ilo >= ihi) {
	return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2][col];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
	if (A[uhi][col] > pivot) {
	    // Here, we can reduce the size of the unpartitioned region and
	    // try again.
	    uhi--;
	} else {
	    // Here, A[uhi] <= pivot, so swap entries at indices ulo and
	    // uhi.
	    tempEntry = A[ulo];
	    A[ulo] = A[uhi];
	    A[uhi] = tempEntry;
	    // After the swap, A[ulo] <= pivot.
	    if (A[ulo][col] < pivot) {
		// Swap entries at indices ieq and ulo.
		tempEntry = A[ieq];
		A[ieq] = A[ulo];
		A[ulo] = tempEntry;
		// After the swap, A[ieq] < pivot, so we need to change
		// ieq.
		ieq++;
		// We also need to change ulo, but we also need to do
		// that when A[ulo] = pivot, so we do it after this if
		// statement.
	    }
	    // Once again, we can reduce the size of the unpartitioned
	    // region and try again.
	    ulo++;
	}
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    Qksort_Matrix(A, col, ilo, ieq - 1);
    Qksort_Matrix(A, col, uhi + 1, ihi);
}

/********************************************************/

void Print_Site(allseq *alldata, int num, int n_otu, char *sep, int stepsize)
{
  int i,j;
  For(i,n_otu)
    {
      printf("%s   ",alldata->c_seq[i]->name);
      For(j,stepsize)
	printf("%c",alldata->c_seq[i]->state[num+j]);
      printf("%s",sep);
    }
  fprintf(stderr,"%s",sep);
}

/*********************************************************/

void Print_Site_Lk(arbre *tree, FILE *fp)
{
  int site;
  int catg;
  char *s;
  phydbl postmean;

  if(!tree->io->print_site_lnl)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(!tree->io->print_trace)
    {
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      
      fprintf(fp,"Note : P(D|M) is the probability of site D given the model M (i.e., the site likelihood)\n");
      if(tree->mod->n_catg > 1 || tree->mod->invar)
	fprintf(fp,"P(D|M,rr[x]) is the probability of site D given the model M and the relative rate\nof evolution rr[x], where x is the class of rate to be considered.\nWe have P(D|M) = \\sum_x P(x) x P(D|M,rr[x]).\n");
      fprintf(fp,"\n\n");
      
      sprintf(s,"Site");
      fprintf(fp, "%-7s",s);
      
      sprintf(s,"P(D|M)");
      fprintf(fp,"%-16s",s);
      
      if(tree->mod->n_catg > 1)
	{
	  For(catg,tree->mod->n_catg)
	    {
	      sprintf(s,"P(D|M,rr[%d]=%5.4f)",catg+1,tree->mod->gamma_rr[catg]);
	      fprintf(fp,"%-22s",s);
	    }
	  
	  sprintf(s,"Posterior mean");
	  fprintf(fp,"%-22s",s);
	}
      
      
      if(tree->mod->invar)
	{
	  sprintf(s,"P(D|M,rr[0]=0)");
	  fprintf(fp,"%-16s",s);
	}
      fprintf(fp,"\n");
      
      For(site,tree->data->init_len)
	{
	  fprintf(fp,"%-7d",site+1);
	  fprintf(fp,"%-16g",(phydbl)exp(tree->site_lk[tree->data->sitepatt[site]]));      
	  if(tree->mod->n_catg > 1)
	    {
	      For(catg,tree->mod->n_catg)
		fprintf(fp,"%-22g",(phydbl)exp(tree->log_site_lk_cat[catg][tree->data->sitepatt[site]]));

	      postmean = .0;
	      For(catg,tree->mod->n_catg) 
		postmean += 
		tree->mod->gamma_rr[catg] * 
		exp(tree->log_site_lk_cat[catg][tree->data->sitepatt[site]]) * 
		tree->mod->gamma_r_proba[catg];
	      postmean /= exp(tree->site_lk[tree->data->sitepatt[site]]);

	      fprintf(fp,"%-22g",postmean);
	    }
	  if(tree->mod->invar)
	    {
	      if((phydbl)tree->data->invar[tree->data->sitepatt[site]] > -0.5)
		fprintf(fp,"%-16g",tree->mod->pi[tree->data->invar[tree->data->sitepatt[site]]]);
	      else
		fprintf(fp,"%-16g",0.0);
	    }
	  fprintf(fp,"\n");
	}
      Free(s);
    }
  else
    {
      For(site,tree->data->init_len)
	fprintf(fp,"%.2f\t",tree->site_lk[tree->data->sitepatt[site]]);
      fprintf(fp,"\n");
    }
}


/*********************************************************/

void Print_Seq(seq **data, int n_otu)
{
  int i,j;

  printf("%d\t%d\n",n_otu,data[0]->len);
  For(i,n_otu)
    {
/*       For(j,30) */
/* 	{ */
/* 	  if(j<(int)strlen(data[i]->name)) */
/* 	     putchar(data[i]->name[j]); */
/* 	  else putchar(' '); */
/* 	} */
      printf("%10d  ",i);
      For(j,data[i]->len)
	{
	  printf("%c",data[i]->state[j]);
	}
      printf("\n");
    }
}

/*********************************************************/

void Print_CSeq(FILE *fp, allseq *alldata)
{
    int i,j,k;
  int n_otu;

  n_otu = alldata->n_otu;
  fprintf(fp,"%d\t%d\n",n_otu,alldata->init_len);
  For(i,n_otu)
    {
      For(j,50)
	{
	  if(j<(int)strlen(alldata->c_seq[i]->name))
	     fputc(alldata->c_seq[i]->name[j],fp);
	  else fputc(' ',fp);
	}

      For(j,alldata->crunch_len)
	{
	  For(k,alldata->wght[j])
	    fprintf(fp,"%c",alldata->c_seq[i]->state[j]);
	}
      fprintf(fp,"\n");
    }
  fprintf(fp,"\n");

/*   printf("\t"); */
/*   For(j,alldata->crunch_len) */
/*     printf("%.0f ",alldata->wght[j]); */
/*   printf("\n"); */
}

/*********************************************************/

void Order_Tree_Seq(arbre *tree, seq **data)
{
    int i,j,n_otu;
    seq *buff;

    n_otu = tree->n_otu;

    For(i,n_otu)
      {
	For(j,n_otu)
	  {
	    if(!strcmp(tree->noeud[i]->name,data[j]->name))
	      break;
	  }
	buff = data[j];
	data[j] = data[i];
	data[i] = buff;
      }
}

/*********************************************************/

void Order_Tree_CSeq(arbre *tree, allseq *data)
{
    int i,j,n_otu_tree,n_otu_seq;
    seq *buff;


    n_otu_tree = tree->n_otu;
    n_otu_seq  = data->n_otu;


    if(n_otu_tree != n_otu_seq)
        {
            /*       printf("%d(tree) != %d(seq) \n",n_otu_tree,n_otu_seq); */
            Warn_And_Exit("\n. The number of tips in the tree is not the same as the number of sequences\n");
        }
    For(i,MAX(n_otu_tree,n_otu_seq))
        {
            For(j,MIN(n_otu_tree,n_otu_seq))
                {
                    if(!strcmp(tree->noeud[i]->name,data->c_seq[j]->name))
                        break;
                }

            if(j==MIN(n_otu_tree,n_otu_seq))
                {
                    printf("\n. Err: %s is not found in sequence data set\n",
                           tree->noeud[i]->name);
                    Warn_And_Exit("");
                }
            buff = data->c_seq[j];
            data->c_seq[j] = data->c_seq[i];
            data->c_seq[i] = buff;
        }
}

/*********************************************************/

matrix *Make_Mat(int n_otu)
{
  matrix *mat;
  int i;

  mat = (matrix *)mCalloc(1,sizeof(matrix));

  mat->n_otu = n_otu;

  mat->P        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->Q        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->dist     = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->on_off   = (int *)mCalloc(n_otu,sizeof(int));
  mat->name     = (char **)mCalloc(n_otu,sizeof(char *));
  mat->tip_node = (node **)mCalloc(n_otu,sizeof(node *));


  For(i,n_otu)
    {
      mat->P[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->Q[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->dist[i] = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->name[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

  return mat;
}

/*********************************************************/

void Init_Mat(matrix *mat, allseq *data)
{
  int i;

  mat->n_otu = data->n_otu;
  mat->r = mat->n_otu;
  mat->curr_int = mat->n_otu;
  mat->method = 1;

  For(i,data->n_otu)
    {
      strcpy(mat->name[i],data->c_seq[i]->name);
      mat->on_off[i] = 1;
    }
}

/*********************************************************/

arbre *Make_Tree_From_Scratch(int n_otu, allseq *data)
{
  arbre *tree;

  tree = Make_Tree(n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);
  Make_Tree_Path(tree);
  Make_List_Of_Reachable_Tips(tree);
  if(data)
    {
      Copy_Tax_Names_To_Tip_Labels(tree,data);
      tree->data = data;
    }
  return tree;
}

/*********************************************************/

arbre *Make_Tree(int n_otu)
{
  arbre *tree;
  int i;
  tree = (arbre *)mCalloc(1,sizeof(arbre ));
  Init_Tree(tree,n_otu);
  tree->t_dir = (int **)mCalloc(2*n_otu-2,sizeof(int *));
  For(i,2*n_otu-2) tree->t_dir[i] = (int *)mCalloc(2*n_otu-2,sizeof(int));
  return tree;
}

/*********************************************************/

void Make_Tree_Path(arbre *tree)
{
  tree->curr_path = (node **)mCalloc(tree->n_otu,sizeof(node *));
}

/*********************************************************/

void Make_All_Tree_Nodes(arbre *tree)
{
  int i;
  tree->noeud          = (node **)mCalloc(2*tree->n_otu-2,sizeof(node *));
  tree->t_dead_nodes   = (node **)mCalloc(2*tree->n_otu-2,sizeof(node *));

  For(i,2*tree->n_otu-2)
    {
      tree->noeud[i] = (node *)Make_Node_Light(i);
      if(i < tree->n_otu) tree->noeud[i]->tax = 1;
      else                tree->noeud[i]->tax = 0;
    }
}

/*********************************************************/

void Make_All_Tree_Edges(arbre *tree)
{
  int i;

  tree->t_edges      = (edge **)mCalloc(2*tree->n_otu-3,sizeof(edge *));
  tree->t_dead_edges = (edge **)mCalloc(2*tree->n_otu-3,sizeof(edge *));

  For(i,2*tree->n_otu-3) tree->t_edges[i] = (edge *)Make_Edge_Light(NULL,NULL,i);
}

/*********************************************************/

void Copy_Tax_Names_To_Tip_Labels(arbre *tree, allseq *data)
{
  int i;

  For(i,tree->n_otu)
    {
      strcpy(tree->noeud[i]->name,data->c_seq[i]->name);
      tree->noeud[i]->tax = 1;
      tree->noeud[i]->num = i;
    }
}

/*********************************************************/

void Print_Dist(matrix *mat)
{
  int i,j;

  For(i,mat->n_otu)
    {
      printf("%s ",mat->name[i]);

      For(j,mat->n_otu)
	printf("%9.6f ",mat->dist[i][j]);
      printf("\n");
    }
}

/*********************************************************/

void Print_Node(node *a, node *d, arbre *tree)
{
  int i;
  int dir;
  dir = -1;
  For(i,3) if(a->v[i] == d) {dir = i; break;}
  printf("Node nums = %3d %3d  (dir=%d);",a->num,d->num,dir);
  printf("Node names = '%s' '%s' ; ",a->name,d->name);
  For(i,3) if(a->v[i] == d)
    {
      printf("Branch num = %3d (%d %d) %f",
	     a->b[i]->num,a->b[i]->left->num,
	     a->b[i]->rght->num,a->b[i]->l);
      if(a->b[i]->left->tax) printf(" WARNING LEFT->TAX!");
      break;
    }
  printf("\n");

  if(d->tax) return;
  else
    For(i,3)
      if(d->v[i] != a) Print_Node(d,d->v[i],tree);
}

/*********************************************************/

void Share_Lk_Struct(arbre *t_full, arbre *t_empt)
{
  int i,j,n_otu;
  edge *b_e,*b_f;
  node *n_e, *n_f;

  n_otu                   = t_full->n_otu;
  t_empt->n_root          = t_full->n_root;
  t_empt->e_root          = t_full->e_root;
  t_empt->c_lnL_sorted    = t_full->c_lnL_sorted;
  t_empt->log_site_lk_cat = t_full->log_site_lk_cat;
  t_empt->site_lk         = t_full->site_lk;
  t_empt->triplet_struct  = t_full->triplet_struct;
  t_empt->log_lks_aLRT    = t_full->log_lks_aLRT;

  For(i,2*n_otu-3)
    {
      b_f = t_full->t_edges[i];
      b_e = t_empt->t_edges[i];

      b_e->Pij_rr = b_f->Pij_rr;

      b_e->nni = b_f->nni;
    }


  for(i=n_otu;i<2*n_otu-2;i++)
    {
      n_f = t_full->noeud[i];
      n_e = t_empt->noeud[i];
            
      For(j,3)
	{
	  if(n_f->b[j]->left == n_f)
	    {
	      if(n_e->b[j]->left == n_e)
		{
		  n_e->b[j]->p_lk_left        = n_f->b[j]->p_lk_left;
		  n_e->b[j]->sum_scale_f_left = n_f->b[j]->sum_scale_f_left;
		  n_e->b[j]->p_lk_tip_l       = n_f->b[j]->p_lk_tip_l;
		}
	      else
		{
		  n_e->b[j]->p_lk_rght        = n_f->b[j]->p_lk_left;
		  n_e->b[j]->sum_scale_f_rght = n_f->b[j]->sum_scale_f_left;
		  n_e->b[j]->p_lk_tip_r       = n_f->b[j]->p_lk_tip_l;
		}
	    }
	  else
	    {
	      if(n_e->b[j]->rght == n_e)
		{
		  n_e->b[j]->p_lk_rght        = n_f->b[j]->p_lk_rght;
		  n_e->b[j]->sum_scale_f_rght = n_f->b[j]->sum_scale_f_rght;
		  n_e->b[j]->p_lk_tip_r       = n_f->b[j]->p_lk_tip_r;
		}
	      else
		{
		  n_e->b[j]->p_lk_left        = n_f->b[j]->p_lk_rght;
		  n_e->b[j]->sum_scale_f_left = n_f->b[j]->sum_scale_f_rght;
		  n_e->b[j]->p_lk_tip_l       = n_f->b[j]->p_lk_tip_r;
		}
	    }
	}
    }

  For(i,n_otu)
    {
      n_f = t_full->noeud[i];
      n_e = t_empt->noeud[i];

      if(n_f->b[0]->rght == n_f)
	{
	  n_e->b[0]->p_lk_rght        = n_f->b[0]->p_lk_rght;
	  n_e->b[0]->sum_scale_f_rght = n_f->b[0]->sum_scale_f_rght;
	  n_e->b[0]->p_lk_tip_r       = n_f->b[0]->p_lk_tip_r;
	}
      else
	{
	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
    }
}

/*********************************************************/

void Share_Spr_Struct(arbre *t_full, arbre *t_empt)
{
  t_empt->size_spr_list = t_full->size_spr_list;
  t_empt->spr_list      = t_full->spr_list;
  t_empt->best_spr      = t_full->best_spr;
}

/*********************************************************/

void Share_Pars_Struct(arbre *t_full, arbre *t_empt)
{
  int i;

  t_empt->site_pars = t_full->site_pars;
  t_empt->step_mat  = t_full->step_mat;

  For(i,2*t_full->n_otu-3)
    {
      t_empt->t_edges[i]->ui_l     = t_full->t_edges[i]->ui_l;
      t_empt->t_edges[i]->ui_r     = t_full->t_edges[i]->ui_r;

      t_empt->t_edges[i]->pars_l   = t_full->t_edges[i]->pars_l;
      t_empt->t_edges[i]->pars_r   = t_full->t_edges[i]->pars_r;

      t_empt->t_edges[i]->p_pars_l = t_full->t_edges[i]->p_pars_l;
      t_empt->t_edges[i]->p_pars_r = t_full->t_edges[i]->p_pars_r;
    }
}

/*********************************************************/

void Share_List_Of_Reachable_Tips_Struct(arbre *t_full, arbre *t_empt)
{
  int i;

  For(i,2*t_full->n_otu-2)
    {
      t_empt->noeud[i]->list_of_reachable_tips = t_full->noeud[i]->list_of_reachable_tips;
      t_empt->noeud[i]->n_of_reachable_tips    = t_full->noeud[i]->n_of_reachable_tips;
    }
}

/*********************************************************/

void Print_Mat(matrix *mat)
{
  int i,j;

  printf("%d",mat->n_otu);
  printf("\n");

  For(i,mat->n_otu)
    {
      For(j,13)
	{
	  if(j>=(int)strlen(mat->name[i])) putchar(' ');
	  else putchar(mat->name[i][j]);
	}

      For(j,mat->n_otu)
	{
	  if(mat->dist[i][j] == -1)
	    printf("   -     ");
	  else
	    printf("%7.8f  ",mat->dist[i][j]);
	}
      printf("\n");
    }
}

/*********************************************************/

int Sort_Edges_NNI_Score(arbre *tree, edge **sorted_edges, int n_elem)
{
  int i,j;
  edge *buff;

  For(i,n_elem-1)
    {
      for(j=i+1;j<n_elem;j++)
	{
	  if(sorted_edges[j]->nni->score  < sorted_edges[i]->nni->score)
	    {
	      buff = sorted_edges[j];
	      sorted_edges[j] = sorted_edges[i];
	      sorted_edges[i] = buff;
	    }
	}
    }
  return 1;
}

/*********************************************************/

void NNI(arbre *tree, edge *b_fcus, int do_swap)
{
  int l_r, r_l, l_v1, l_v2, r_v3, r_v4;
  node *v1,*v2,*v3,*v4;
  phydbl lk0, lk1, lk2;
  phydbl lk0_init, lk1_init, lk2_init;
  phydbl bl_init;
  phydbl l0,l1,l2;
  phydbl l_infa, l_infb, l_max;
/*   phydbl lk_infa, lk_infb, lk_max; */
  phydbl lk_init;

  bl_init                = b_fcus->l;
  lk_init                = tree->c_lnL;

  b_fcus->nni->init_l    = b_fcus->l;
  b_fcus->nni->init_lk   = tree->c_lnL;;

  b_fcus->nni->best_conf = 0;
  b_fcus->nni->score     = +1.0;

  lk0 = lk1 = lk2        = UNLIKELY;
  v1 = v2 = v3 = v4      = NULL;

  l_r = r_l = l_v1 = l_v2 = r_v3 = r_v4 = -1;

  l_r                    = b_fcus->l_r;
  r_l                    = b_fcus->r_l;

  v1                     = b_fcus->left->v[b_fcus->l_v1];
  v2                     = b_fcus->left->v[b_fcus->l_v2];
  v3                     = b_fcus->rght->v[b_fcus->r_v1];
  v4                     = b_fcus->rght->v[b_fcus->r_v2];

  if(v1->num < v2->num)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  if(v3->num < v4->num)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  l0 = l1 = l2 = -1.;

  
  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
  tree->both_sides = 1;

  lk1_init = Update_Lk_At_Given_Edge(b_fcus,tree);

  l_infa = 10.*b_fcus->l;
  l_max  = b_fcus->l;
  l_infb = BL_MIN;

  if(tree->mod->s_opt->fast_nni)
    {
      Fast_Br_Len(b_fcus,tree);
      lk1 = Lk_At_Given_Edge(b_fcus,tree);
    }
  else
    {
      lk1 = Br_Len_Brent(l_infa,l_max,l_infb,
			 tree->mod->s_opt->min_diff_lk_local,
			 b_fcus,tree,
			 tree->mod->s_opt->brent_it_max);
    }

  if(lk1 < lk1_init - tree->mod->s_opt->min_diff_lk_local)
    {
      printf("%f %f %f %f\n",l_infa,l_max,l_infb,b_fcus->l);
      printf("%f -- %f \n",lk1_init,lk1);
      printf("\n. Err. in NNI (1)\n");
    }

  l1  = b_fcus->l;
  Swap(v3,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/


  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
  b_fcus->l = bl_init;
  tree->both_sides = 1;

  lk2_init = Update_Lk_At_Given_Edge(b_fcus,tree);

  l_infa = 10.*b_fcus->l;
  l_max  = b_fcus->l;
  l_infb = BL_MIN;

  if(tree->mod->s_opt->fast_nni)
    {
      Fast_Br_Len(b_fcus,tree);
      lk2 = Lk_At_Given_Edge(b_fcus,tree);
    }
  else
    {
      lk2 = Br_Len_Brent(l_infa,l_max,l_infb,
			 tree->mod->s_opt->min_diff_lk_local,
			 b_fcus,tree,
			 tree->mod->s_opt->brent_it_max);
    }

  if(lk2 < lk2_init - tree->mod->s_opt->min_diff_lk_local)
    {
      printf("%f %f %f %f\n",l_infa,l_max,l_infb,b_fcus->l);
      printf("%f -- %f \n",lk2_init,lk2);
      printf("\n. Err. in NNI (2)\n");
   }

  l2  = b_fcus->l;
  Swap(v4,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/



  /***********/
   b_fcus->l = bl_init;
   tree->both_sides = 1;

   lk0_init = Update_Lk_At_Given_Edge(b_fcus,tree);

   if(fabs(lk0_init - lk_init) > tree->mod->s_opt->min_diff_lk_local)
     {
       printf("\n. lk_init = %f; lk = %f diff = %f\n",
	      lk_init,
	      lk0_init,
	      lk_init-lk0_init);
       printf("\n. Curr_lnL = %f\n",Return_Lk(tree));
       Warn_And_Exit("\n. Err. in NNI (3)\n");
     }

   l_infa = 10.*b_fcus->l;
   l_max  = b_fcus->l;
   l_infb = BL_MIN;

   if(tree->mod->s_opt->fast_nni)
     {
       Fast_Br_Len(b_fcus,tree);
       lk0 = Lk_At_Given_Edge(b_fcus,tree);
     }
   else
     {
       lk0 = Br_Len_Brent(l_infa,l_max,l_infb,
			  tree->mod->s_opt->min_diff_lk_local,
			  b_fcus,tree,
			  tree->mod->s_opt->brent_it_max);
     }

   if(lk0 < lk_init - tree->mod->s_opt->min_diff_lk_local)
     {
       printf("\n\n%f %f %f %f\n",l_infa,l_max,l_infb,b_fcus->l);
       printf("%f -- %f \n",lk0_init,lk0);
       printf("\n. Err. in NNI (3)\n");
       Warn_And_Exit("\n");
     }

   l0  = b_fcus->l;
   /***********/

   b_fcus->nni->lk0 = lk0;
   b_fcus->nni->lk1 = lk1;
   b_fcus->nni->lk2 = lk2;

   b_fcus->nni->l0  = l0;
   b_fcus->nni->l1  = l1;
   b_fcus->nni->l2  = l2;

   b_fcus->nni->score = lk0 - MAX(lk1,lk2);

   if((b_fcus->nni->score <  tree->mod->s_opt->min_diff_lk_local) &&
      (b_fcus->nni->score > -tree->mod->s_opt->min_diff_lk_local))
     {
       b_fcus->nni->score = .0;
       b_fcus->nni->lk1 = b_fcus->nni->lk0;
       b_fcus->nni->lk2 = b_fcus->nni->lk0;
     }

   if(lk0 > MAX(lk1,lk2))
     {
       b_fcus->nni->best_conf    = 0;
       b_fcus->nni->best_l       = l0;
       b_fcus->nni->swap_node_v1 = NULL;
       b_fcus->nni->swap_node_v2 = NULL;
       b_fcus->nni->swap_node_v3 = NULL;
       b_fcus->nni->swap_node_v4 = NULL;
      }
   else if(lk1 > MAX(lk0,lk2))
     {
       b_fcus->nni->best_conf    = 1;
       b_fcus->nni->best_l       = l1;
       b_fcus->nni->swap_node_v1 = v2;
       b_fcus->nni->swap_node_v2 = b_fcus->left;
       b_fcus->nni->swap_node_v3 = b_fcus->rght;
       b_fcus->nni->swap_node_v4 = v3;
     }
   else if(lk2 > MAX(lk0,lk1))
     {
       b_fcus->nni->best_conf    = 2;
       b_fcus->nni->best_l       = l2;
       b_fcus->nni->swap_node_v1 = v2;
       b_fcus->nni->swap_node_v2 = b_fcus->left;
       b_fcus->nni->swap_node_v3 = b_fcus->rght;
       b_fcus->nni->swap_node_v4 = v4;
     }
   else
     {
       b_fcus->nni->score        = +1.0;
       b_fcus->nni->best_conf    = 0;
       b_fcus->nni->best_l       = l0;
       b_fcus->nni->swap_node_v1 = NULL;
       b_fcus->nni->swap_node_v2 = NULL;
       b_fcus->nni->swap_node_v3 = NULL;
       b_fcus->nni->swap_node_v4 = NULL;
     }

   if((do_swap) && ((lk1 > lk0) || (lk2 > lk0)))
     {
      tree->n_swap++;
      printf("Swap edge %d -> %f\n",b_fcus->num,MAX(lk1,lk2));

      if(lk1 > lk2)
	 {
	   tree->best_lnL = lk1;
	   Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
	   b_fcus->l = l1;
	   tree->both_sides = 1;
	   Lk(tree);
	 }
       else
	 {
	   tree->best_lnL = lk2;
	   Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
	   b_fcus->l = l2;
	   tree->both_sides = 1;
	   Lk(tree);
	 }
     }
   else
     {
       b_fcus->l = bl_init;
       Update_PMat_At_Given_Edge(b_fcus,tree);
       tree->c_lnL = lk_init;
     }
}

/*********************************************************/

void NNI_Pars(arbre *tree, edge *b_fcus, int do_swap)
{
  int l_r, r_l, l_v1, l_v2, r_v3, r_v4;
  node *v1,*v2,*v3,*v4;
  int pars0, pars1, pars2;
  int pars_init;

  pars_init              = tree->c_pars;
  b_fcus->nni->best_conf = 0;
  b_fcus->nni->score     = +1.0;

  pars0 = pars1 = pars2  = 0;
  v1 = v2 = v3 = v4      = NULL;

  l_r = r_l = l_v1 = l_v2 = r_v3 = r_v4 = -1;

  l_r                    = b_fcus->l_r;
  r_l                    = b_fcus->r_l;

  v1                     = b_fcus->left->v[b_fcus->l_v1];
  v2                     = b_fcus->left->v[b_fcus->l_v2];
  v3                     = b_fcus->rght->v[b_fcus->r_v1];
  v4                     = b_fcus->rght->v[b_fcus->r_v2];

  if(v1->num < v2->num)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  if(v3->num < v4->num)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  
  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
  tree->both_sides = 1;
  pars1 = Update_Pars_At_Given_Edge(b_fcus,tree);
  Swap(v3,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/

  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
  tree->both_sides = 1;
  pars2 = Update_Pars_At_Given_Edge(b_fcus,tree);
  Swap(v4,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/


  /***********/
   tree->both_sides = 1;
   pars0 = Update_Pars_At_Given_Edge(b_fcus,tree);
 
   if(pars0 != pars_init)
     {
       printf("\n. pars_init = %d; pars0 = %d\n",
	      pars_init,
	      pars0);
       Warn_And_Exit("\n. Err. in NNI (3)\n");
     }
   /***********/

   tree->c_pars = pars0;

   b_fcus->nni->score = MIN(pars1,pars2) - pars0;

   if(pars0 < MIN(pars1,pars2))
     {
       b_fcus->nni->best_conf    = 0;
       b_fcus->nni->swap_node_v1 = NULL;
       b_fcus->nni->swap_node_v2 = NULL;
       b_fcus->nni->swap_node_v3 = NULL;
       b_fcus->nni->swap_node_v4 = NULL;
      }
   else if(pars1 < MIN(pars0,pars2))
     {
       b_fcus->nni->best_conf    = 1;
       b_fcus->nni->swap_node_v1 = v2;
       b_fcus->nni->swap_node_v2 = b_fcus->left;
       b_fcus->nni->swap_node_v3 = b_fcus->rght;
       b_fcus->nni->swap_node_v4 = v3;
     }
   else if(pars2 > MIN(pars0,pars1))
     {
       b_fcus->nni->best_conf    = 2;
       b_fcus->nni->swap_node_v1 = v2;
       b_fcus->nni->swap_node_v2 = b_fcus->left;
       b_fcus->nni->swap_node_v3 = b_fcus->rght;
       b_fcus->nni->swap_node_v4 = v4;
     }
   else
     {
       b_fcus->nni->score        = +1.0;
       b_fcus->nni->swap_node_v1 = NULL;
       b_fcus->nni->swap_node_v2 = NULL;
       b_fcus->nni->swap_node_v3 = NULL;
       b_fcus->nni->swap_node_v4 = NULL;
     }
}

/*********************************************************/

void Swap(node *a, node *b, node *c, node *d, arbre *tree)
{
  int ab, ba, cd, dc, bc;
  int i;


  /* \             /d      \             /a
   *  \           /         \           /
   *   \b__...__c/    ->     \b__...__c/
   *   /         \	     /		\
   *  /           \	    /		 \
   * /a            \  	   /d             \ 
   *
   * nodes b and c are not necessarily on the same branch 
   */


#ifdef DEBUG
  if(!a || !b || !c || !d)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif


  ab = ba = cd = dc = bc = -1;

  For(i,3) if(a->v[i] == b) { ab = i; break; }
  For(i,3) if(b->v[i] == a) { ba = i; break; }
  For(i,3) if(c->v[i] == d) { cd = i; break; }
  For(i,3) if(d->v[i] == c) { dc = i; break; }
  For(i,3) if(b->v[i] == c) { bc = i; break; }

#ifdef DEBUG
  if(ab < 0 || ba < 0 || cd < 0 || dc < 0)
    {
      printf("\n. Nodes %d %d %d %d\n",a->num,b->num,c->num,d->num);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif

  a->v[ab] = c;
  d->v[dc] = b;
  b->v[ba] = d;
  c->v[cd] = a;
  b->b[ba] = d->b[dc];
  c->b[cd] = a->b[ab];

  (a->b[ab]->left == b)?
  (a->b[ab]->left = c):
  (a->b[ab]->rght = c);

  (d->b[dc]->left == c)?
  (d->b[dc]->left = b):
  (d->b[dc]->rght = b);

  For(i,3)
    {
      if(a->b[ab]->left->v[i] == a->b[ab]->rght) a->b[ab]->l_r = i;
      if(a->b[ab]->rght->v[i] == a->b[ab]->left) a->b[ab]->r_l = i;
      if(d->b[dc]->left->v[i] == d->b[dc]->rght) d->b[dc]->l_r = i;
      if(d->b[dc]->rght->v[i] == d->b[dc]->left) d->b[dc]->r_l = i;
    }


  a->b[ab]->l_v1 = a->b[ab]->l_v2 =
  a->b[ab]->r_v1 = a->b[ab]->r_v2 =
  d->b[dc]->l_v1 = d->b[dc]->l_v2 =
  d->b[dc]->r_v1 = d->b[dc]->r_v2 = -1;


  For(i,3)
    {
      if(i != a->b[ab]->l_r)
	{
	  if(a->b[ab]->l_v1 < 0) a->b[ab]->l_v1 = i;
	  else a->b[ab]->l_v2 = i;
	}
      if(i != a->b[ab]->r_l)
	{
	  if(a->b[ab]->r_v1 < 0) a->b[ab]->r_v1 = i;
	  else a->b[ab]->r_v2 = i;
	}
      if(i != d->b[dc]->l_r)
	{
	  if(d->b[dc]->l_v1 < 0) d->b[dc]->l_v1 = i;
	  else d->b[dc]->l_v2 = i;
	}
      if(i != d->b[dc]->r_l)
	{
	  if(d->b[dc]->r_v1 < 0) d->b[dc]->r_v1 = i;
	  else d->b[dc]->r_v2 = i;
	}
    }
  Update_Dirs(tree);
}

/*********************************************************/

void Update_All_Partial_Lk(edge *b_fcus, arbre *tree)
{

  Update_SubTree_Partial_Lk(b_fcus->left->b[b_fcus->l_v1],
			    b_fcus->left,
			    b_fcus->left->v[b_fcus->l_v1],
			    tree);

  Update_SubTree_Partial_Lk(b_fcus->left->b[b_fcus->l_v2],
			    b_fcus->left,
			    b_fcus->left->v[b_fcus->l_v2],
			    tree);

  Update_SubTree_Partial_Lk(b_fcus->rght->b[b_fcus->r_v1],
			    b_fcus->rght,
			    b_fcus->rght->v[b_fcus->r_v1],
			    tree);

  Update_SubTree_Partial_Lk(b_fcus->rght->b[b_fcus->r_v2],
			    b_fcus->rght,
			    b_fcus->rght->v[b_fcus->r_v2],
			    tree);

  tree->c_lnL = Lk_At_Given_Edge(b_fcus,tree);
}

/*********************************************************/

void Update_SubTree_Partial_Lk(edge *b_fcus, node *a, node *d, arbre *tree)
{
  int i;

  Update_P_Lk(tree,b_fcus,a);
  if(d->tax) return;
  else For(i,3) if(d->v[i] != a)
    Update_SubTree_Partial_Lk(d->b[i],d,d->v[i],tree);
}

/*********************************************************/

allseq *Make_Cseq(int n_otu, int crunch_len, int init_len, char **sp_names)
{
  allseq *alldata;
  int j;

  alldata                        = (allseq *)mCalloc(1,sizeof(allseq));
  alldata->n_otu                 = n_otu;
  alldata->c_seq                 = (seq **)mCalloc(n_otu,sizeof(seq *));
  alldata->b_frq                 = (phydbl *)mCalloc(T_MAX_ALPHABET,sizeof(phydbl));
  alldata->wght                  = (int *)mCalloc(crunch_len,sizeof(int));
  alldata->ambigu                = (short int *)mCalloc(crunch_len,sizeof(short int));
  alldata->invar                 = (short int *)mCalloc(crunch_len,sizeof(short int));
  alldata->sitepatt              = (int *)mCalloc(  init_len,sizeof(int ));

  alldata->crunch_len = crunch_len;
  alldata->init_len   = init_len;
  alldata->obs_pinvar = .0;

  For(j,n_otu)
    {
      alldata->c_seq[j]            = (seq *)mCalloc(1,sizeof(seq));
      alldata->c_seq[j]->name      = (char *)mCalloc((int)(strlen(sp_names[j])+1),sizeof(char));
      strcpy(alldata->c_seq[j]->name,sp_names[j]);
      alldata->c_seq[j]->state     = (char *)mCalloc(crunch_len,sizeof(char));
      alldata->c_seq[j]->is_ambigu = (short int *)mCalloc(crunch_len,sizeof(short int));
    }

  return alldata;
}

/*********************************************************/

arbrelist *Make_Treelist(int list_size)
{
  arbrelist *tlist;

  tlist = (arbrelist *)mCalloc(1,sizeof(arbrelist));
  tlist->list_size = list_size;
  tlist->tree = (arbre **)mCalloc(list_size,sizeof(arbre *));

  return tlist;
}


/*********************************************************/

void Copy_Seq_Names_To_Tip_Labels(arbre *tree, allseq *data)
{
  int i;
  For(i,tree->n_otu)
    {
      strcpy(tree->noeud[i]->name,data->c_seq[i]->name);
    }
}

/*********************************************************/

allseq *Copy_Cseq(allseq *ori, int len, int ns)
{
  allseq *new;
  int i,j,n_otu;
  char **sp_names;

  n_otu = ori->n_otu;

  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  For(i,n_otu)
    {
      sp_names[i] = (char *)mCalloc(strlen(ori->c_seq[i]->name)+1,sizeof(char));
      strcpy(sp_names[i],ori->c_seq[i]->name);
    }

  new = Make_Cseq(n_otu, len, ori->init_len, sp_names);

  new->obs_pinvar = ori->obs_pinvar;

  For(i,ori->init_len) new->sitepatt[i] = ori->sitepatt[i];

  For(j,ori->crunch_len)
    {
      For(i,ori->n_otu) 
	{
	  new->c_seq[i]->state[j]     = ori->c_seq[i]->state[j];
	  new->c_seq[i]->is_ambigu[j] = ori->c_seq[i]->is_ambigu[j];
	}

      new->wght[j]   = ori->wght[j];
      new->ambigu[j] = ori->ambigu[j];
      new->invar[j]  = ori->invar[j];
    }

  For(i,ori->n_otu)
    {
      new->c_seq[i]->len = ori->c_seq[i]->len;
      strcpy(new->c_seq[i]->name,ori->c_seq[i]->name);
    }

  new->init_len           = ori->init_len;
  new->clean_len          = ori->clean_len;
  new->crunch_len         = ori->crunch_len;
  For(i,ns) new->b_frq[i] = ori->b_frq[i];
  new->n_otu              = ori->n_otu;

  For(i,n_otu) Free(sp_names[i]);
  Free(sp_names);

  return new;
}

/*********************************************************/

optimiz *Alloc_Optimiz()
{
  optimiz *s_opt;
  s_opt = (optimiz *)mCalloc(1,sizeof(optimiz));
  return s_opt;
}

/*********************************************************/


int Filexists(char *filename)
{
  FILE *fp;
  fp =fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}

/*********************************************************/

FILE *Openfile(char *filename, int mode)
{
  /* mode = 0 -> read */
  /* mode = 1 -> write */
  /* mode = 2 -> append */

  FILE *fp;
  char *s;
  int open_test=0;

/*   s = (char *)mCalloc(T_MAX_FILE,sizeof(char)); */

/*   strcpy(s,filename); */

  s = filename;

  fp = NULL;

  switch(mode)
    {
    case 0 :
      {
	while(!(fp = (FILE *)fopen(s,"r")) && ++open_test<10)
	  {
	    printf("\n. Can't open file '%s', enter a new name : ",s);
	    Getstring_Stdin(s);
	  }
	break;
      }
    case 1 :
      {
	fp = (FILE *)fopen(s,"w");
	break;
      }
    case 2 :
      {
	fp = (FILE *)fopen(s,"a");
	break;
      }

    default : break;

    }

/*   Free(s); */

  return fp;
}

/*********************************************************/

void Print_Fp_Out(FILE *fp_out, time_t t_beg, time_t t_end, arbre *tree, option *io, int n_data_set, int num_tree)
{
  char *s;
  div_t hour,min;

  if((!n_data_set) || (!num_tree)) 
    {
      Print_Banner(fp_out);
    }

  fprintf(fp_out,"\n\n");
  fprintf(fp_out,". Sequence file : %s\n", io->in_seq_file);

  fprintf(fp_out,". Data set : #%d\n",n_data_set);

  if(io->mod->s_opt->random_input_tree)
    fprintf(fp_out,". Random init tree : #%d\n",num_tree+1);
  else if(io->n_trees > 1)
    fprintf(fp_out,". Starting tree number : #%d\n",num_tree+1);
  
  if(io->mod->s_opt->opt_topo)
    fprintf(fp_out,". Tree search : %s",(io->mod->s_opt->topo_search == SPR_MOVE)?("SPRs"):("NNIs"));



  /*was after Sequence file ; moved here FLT*/
  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  if(io->in_tree)
    {
      strcat(strcat(strcat(s,"user tree ("),io->in_tree_file),")");
    }
  else
    {
      if(!io->mod->s_opt->random_input_tree)
	{
	  strcat(s,"BIONJ");
	}
      else
	{
	  strcat(s,"random tree");
	}
    }

  fprintf(fp_out,". Initial tree : %s\n",s);
  Free(s);

  (tree->mod->datatype == NT)?
    (fprintf(fp_out,". Model of nucleotides substitution : %s\n",io->mod->modelname)):
    (fprintf(fp_out,". Model of amino acids substitution : %s\n",io->mod->modelname));



  fprintf(fp_out,". Number of taxa : %d\n",tree->n_otu);/*added FLT*/

  fprintf(fp_out,". Log-likelihood : %.5f\n",tree->c_lnL);/*was last ; moved here FLT*/

  fprintf(fp_out,". Discrete gamma model : %s\n",
	  (tree->mod->n_catg>1)?("Yes"):("No"));
  if(tree->mod->n_catg > 1)
    {
      fprintf(fp_out,"  - Number of categories : %d\n",tree->mod->n_catg);
      fprintf(fp_out,"  - Gamma shape parameter : %.3f\n",tree->mod->alpha);
    }

  if(tree->mod->invar) fprintf(fp_out,". Proportion of invariant : %.3f\n",tree->mod->pinvar);

  /*was before Discrete gamma model ; moved here FLT*/
  if((tree->mod->whichmodel == K80)   ||
     (tree->mod->whichmodel == HKY85) ||
     (tree->mod->whichmodel == F84))
    fprintf(fp_out,". Transition/transversion ratio : %.3f\n",tree->mod->kappa);
  else if(tree->mod->whichmodel == TN93)
    {
      fprintf(fp_out,". Transition/transversion ratio for purines :     %.3f\n",
	      tree->mod->kappa*2.*tree->mod->lambda/(1.+tree->mod->lambda));
      fprintf(fp_out,". Transition/transversion ratio for pyrimidines : %.3f\n",
	      tree->mod->kappa*2./(1.+tree->mod->lambda));
    }

  if(tree->mod->datatype == NT)
    {
      fprintf(fp_out,". Nucleotides frequencies :\n");
      fprintf(fp_out,"  - f(A)=%8.5f\n",tree->mod->pi[0]);
      fprintf(fp_out,"  - f(C)=%8.5f\n",tree->mod->pi[1]);
      fprintf(fp_out,"  - f(G)=%8.5f\n",tree->mod->pi[2]);
      fprintf(fp_out,"  - f(T)=%8.5f\n",tree->mod->pi[3]);
    }

  /*****************************************/
  if((tree->mod->whichmodel == GTR) ||
     (tree->mod->whichmodel == CUSTOM))
    {
      int i,j;

      Update_Qmat_GTR(tree->mod->rr, 
		      tree->mod->rr_val, 
		      tree->mod->rr_num, 
		      tree->mod->pi, 
		      tree->mod->qmat);

      printf("\n");
      fprintf(fp_out,". GTR relative rate parameters : \n\n");
      fprintf(fp_out,"A <-> C   %8.5f\n",  tree->mod->rr[0]);
      fprintf(fp_out,"A <-> G   %8.5f\n",  tree->mod->rr[1]);
      fprintf(fp_out,"A <-> T   %8.5f\n",  tree->mod->rr[2]);
      fprintf(fp_out,"C <-> G   %8.5f\n",  tree->mod->rr[3]);
      fprintf(fp_out,"C <-> T   %8.5f\n",  tree->mod->rr[4]);
      fprintf(fp_out,"G <-> T   %8.5f\n\n",tree->mod->rr[5]);


      fprintf(fp_out,"\n. Instantaneous rate matrix : \n");
      fprintf(fp_out,"\n[A---------C---------G---------T------]\n");
      For(i,4)
	{
	  For(j,4)
	    fprintf(fp_out,"%8.5f  ",tree->mod->qmat[i*4+j]);
	  fprintf(fp_out,"\n");
	}
      fprintf(fp_out,"\n");
    }
  /*****************************************/


  if(io->ratio_test == 1)
    {
      fprintf(fp_out,". aLRT statistics to test branches");
    }
  else if(io->ratio_test == 2)
    {
      /* TO DO
	 cubic approximation ??
      */
      fprintf(fp_out,". aLRT branch supports (cubic approximation, mixture of Chi2s distribution)");
    }


  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );

  min.quot -= hour.quot*60;

  fprintf(fp_out,". Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  if(t_end-t_beg > 60)
    fprintf(fp_out,". -> %d seconds\n",(int)(t_end-t_beg));

  fprintf(fp_out," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");

}

/*********************************************************/
/*FLT wrote this function*/
void Print_Fp_Out_Lines(FILE *fp_out, time_t t_beg, time_t t_end, arbre *tree, option *io, int n_data_set)
{
  char *s;
  /*div_t hour,min;*/

  if (n_data_set==1)
      {

	fprintf(fp_out,". Sequence file : [%s]\n\n", io->in_seq_file);

	(tree->mod->datatype == NT)?
	  (fprintf(fp_out,". Model of nucleotides substitution : %s\n\n",io->mod->modelname)):
	  (fprintf(fp_out,". Model of amino acids substitution : %s\n\n",io->mod->modelname));

	s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	fprintf(fp_out,". Initial tree : [%s]\n\n",
		(!io->in_tree)?("BIONJ"):
		(strcat(strcat(strcat(s,"user tree ("),io->in_tree_file),")")));
	Free(s);

	fprintf(fp_out,"\n");

	/*headline 1*/
	fprintf(fp_out, ". Data\t");

	fprintf(fp_out,"Nb of \t");

	fprintf(fp_out,"Likelihood\t");

	fprintf(fp_out, "Discrete   \t");

	if(tree->mod->n_catg > 1)
	  fprintf(fp_out, "Number of \tGamma shape\t");

	fprintf(fp_out,"Proportion of\t");

	if(tree->mod->whichmodel <= 6)
	  fprintf(fp_out,"Transition/ \t");

	fprintf(fp_out,"Nucleotides frequencies               \t");

	if((tree->mod->whichmodel == GTR) ||
	   (tree->mod->whichmodel == CUSTOM))
	  fprintf(fp_out,"Instantaneous rate matrix              \t");

	/*    fprintf(fp_out,"Time\t");*/

	fprintf(fp_out, "\n");


	/*headline 2*/
	fprintf(fp_out, "  set\t");

	fprintf(fp_out,"taxa\t");

	fprintf(fp_out,"loglk     \t");

	fprintf(fp_out, "gamma model\t");

	if(tree->mod->n_catg > 1)
	  fprintf(fp_out, "categories\tparameter  \t");

	fprintf(fp_out,"invariant    \t");

	if(tree->mod->whichmodel <= 6)
	  fprintf(fp_out,"transversion\t");

	fprintf(fp_out,"f(A)      f(C)      f(G)      f(T)    \t");

	if((tree->mod->whichmodel == GTR) ||
	   (tree->mod->whichmodel == CUSTOM))
	  fprintf(fp_out,"[A---------C---------G---------T------]\t");

	/*    fprintf(fp_out,"used\t");*/

	fprintf(fp_out, "\n");


	/*headline 3*/
	if(tree->mod->whichmodel == TN93)
	  {
	    fprintf(fp_out,"    \t      \t          \t           \t");
	    if(tree->mod->n_catg > 1) fprintf(fp_out,"         \t         \t");
	    fprintf(fp_out,"             \t");
	    fprintf(fp_out,"purines pyrimid.\t");
	    fprintf(fp_out, "\n");
          }

          fprintf(fp_out, "\n");
      }


  /*line items*/

  fprintf(fp_out,"  #%d\t",n_data_set);

  fprintf(fp_out,"%d   \t",tree->n_otu);

  fprintf(fp_out,"%.5f\t",tree->c_lnL);

  fprintf(fp_out,"%s        \t",
	  (tree->mod->n_catg>1)?("Yes"):("No "));
  if(tree->mod->n_catg > 1)
    {
      fprintf(fp_out,"%d        \t",tree->mod->n_catg);
      fprintf(fp_out,"%.3f    \t",tree->mod->alpha);
    }

  /*if(tree->mod->invar)*/
    fprintf(fp_out,"%.3f    \t",tree->mod->pinvar);

  if(tree->mod->whichmodel <= 5)
    {
      fprintf(fp_out,"%.3f     \t",tree->mod->kappa);
    }
  else if(tree->mod->whichmodel == TN93)
    {
      fprintf(fp_out,"%.3f   ",
	      tree->mod->kappa*2.*tree->mod->lambda/(1.+tree->mod->lambda));
      fprintf(fp_out,"%.3f\t",
	      tree->mod->kappa*2./(1.+tree->mod->lambda));
    }


  if(tree->mod->datatype == NT)
    {
      fprintf(fp_out,"%8.5f  ",tree->mod->pi[0]);
      fprintf(fp_out,"%8.5f  ",tree->mod->pi[1]);
      fprintf(fp_out,"%8.5f  ",tree->mod->pi[2]);
      fprintf(fp_out,"%8.5f\t",tree->mod->pi[3]);
    }
  /*
  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );

  min.quot -= hour.quot*60;

  fprintf(fp_out,"%dh%dm%ds\t", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  if(t_end-t_beg > 60)
    fprintf(fp_out,". -> %d seconds\t",(int)(t_end-t_beg));
  */

  /*****************************************/
  if((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
    {
      int i,j;

      For(i,4)
	{
	  if (i!=0) {
	    /*format*/
	    fprintf(fp_out,"      \t     \t          \t           \t");
	    if(tree->mod->n_catg > 1) fprintf(fp_out,"          \t           \t");
	    fprintf(fp_out,"             \t                                      \t");
	  }
	  For(j,4)
	    fprintf(fp_out,"%8.5f  ",tree->mod->qmat[i*4+j]);
	  if (i<3) fprintf(fp_out,"\n");
	}
    }
  /*****************************************/

  fprintf(fp_out, "\n\n");
}

/*********************************************************/

matrix *K80_dist(allseq *data, phydbl g_shape)
{
  int i,j,k;
  int diff;
  phydbl unc_len;
  matrix *mat;
  phydbl **len;

  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  For(i,data->n_otu)
    len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));

  unc_len = .0;

  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);


  For(i,data->c_seq[0]->len)
    {
      For(j,data->n_otu-1)
	{
	  for(k=j+1;k<data->n_otu;k++)
	    {
	      if(((data->c_seq[j]->state[i] == 'A' || data->c_seq[j]->state[i] == 'G') &&
		  (data->c_seq[k]->state[i] == 'C' || data->c_seq[k]->state[i] == 'T'))||
		 ((data->c_seq[j]->state[i] == 'C' || data->c_seq[j]->state[i] == 'T') &&
		  (data->c_seq[k]->state[i] == 'A' || data->c_seq[k]->state[i] == 'G')))
		{
		  diff++;
		  mat->Q[j][k]+=data->wght[i];
		  len[j][k]+=data->wght[i];
		  len[k][j]=len[j][k];
		}

	      else
		if(((data->c_seq[j]->state[i] == 'A' && data->c_seq[k]->state[i] == 'G') ||
		    (data->c_seq[j]->state[i] == 'G' && data->c_seq[k]->state[i] == 'A'))||
		   ((data->c_seq[j]->state[i] == 'C' && data->c_seq[k]->state[i] == 'T') ||
		    (data->c_seq[j]->state[i] == 'T' && data->c_seq[k]->state[i] == 'C')))
		  {
		    diff++;
		    mat->P[j][k]+=data->wght[i];
		    len[j][k]+=data->wght[i];
		    len[k][j]=len[j][k];
		  }
		else
		  if((data->c_seq[j]->state[i] == 'A' ||
		      data->c_seq[j]->state[i] == 'C' ||
		      data->c_seq[j]->state[i] == 'G' ||
		      data->c_seq[j]->state[i] == 'T')&&
		     (data->c_seq[k]->state[i] == 'A' ||
		      data->c_seq[k]->state[i] == 'C' ||
		      data->c_seq[k]->state[i] == 'G' ||
		      data->c_seq[k]->state[i] == 'T'))
		    {
		      len[j][k]+=data->wght[i];
		      len[k][j]=len[j][k];
		    }
	    }
	}
    }


  For(i,data->n_otu-1)
    for(j=i+1;j<data->n_otu;j++)
      {
	if(len[i][j])
	  {
	    mat->P[i][j] /= len[i][j];
	    mat->Q[i][j] /= len[i][j];
	  }
	else
	  {
	    mat->P[i][j] = .5;
	    mat->Q[i][j] = .5;
	  }

	mat->P[j][i] = mat->P[i][j];
	mat->Q[j][i] = mat->Q[i][j];


	if((1-2*mat->P[i][j]-mat->Q[i][j] <= .0) || (1-2*mat->Q[i][j] <= .0))
	  {
	    mat->dist[i][j] = -1.;
	    mat->dist[j][i] = -1.;
	    continue;
	  }

	mat->dist[i][j] = (g_shape/2)*
	  (pow(1-2*mat->P[i][j]-mat->Q[i][j],-1./g_shape) +
	   0.5*pow(1-2*mat->Q[i][j],-1./g_shape) - 1.5);


	if(mat->dist[i][j] > DIST_MAX)
	  {
	    mat->dist[i][j] = DIST_MAX;
	  }
	mat->dist[j][i] = mat->dist[i][j];
      }

  For(i,data->n_otu) free(len[i]);
  free(len);
  return mat;
}

/*********************************************************/

matrix *JC69_Dist(allseq *data, model *mod)
{
  int site,i,j,k;
  phydbl unc_len;
  matrix *mat;
  phydbl **len;


  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  For(i,data->n_otu)
    len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));

  unc_len = .0;

  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);

  Fors(site,data->c_seq[0]->len,mod->stepsize)
    {
      For(j,data->n_otu-1)
	{
	  for(k=j+1;k<data->n_otu;k++)
	    {
	      if((!Is_Ambigu(data->c_seq[j]->state+site,mod->datatype,mod->stepsize)) &&
		 (!Is_Ambigu(data->c_seq[k]->state+site,mod->datatype,mod->stepsize)))
		{
		  len[j][k]+=data->wght[site];
		  len[k][j]=len[j][k];
		  if(strncmp(data->c_seq[j]->state+site,
			     data->c_seq[k]->state+site,
			     mod->stepsize))
		    mat->P[j][k]+=data->wght[site];
		}
	    }
	}
    }
  

  For(i,data->n_otu-1)
    for(j=i+1;j<data->n_otu;j++)
      {
	if(len[i][j])
	  {
	    mat->P[i][j] /= len[i][j];
	  }
	else
	  {
	    mat->P[i][j] = 1.;
	  }

	mat->P[j][i] = mat->P[i][j];

	if((1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]) < .0)
	  {
	    mat->dist[i][j] = DIST_MAX;
	  }
	else
	  mat->dist[i][j] = -(mod->ns-1.)/(mod->ns)*(phydbl)log(1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]);


	if(mat->dist[i][j] > DIST_MAX)
	  {
	    mat->dist[i][j] = DIST_MAX;
	  }
	mat->dist[j][i] = mat->dist[i][j];
      }

  For(i,data->n_otu) free(len[i]);
  free(len);

  return mat;
}

/*********************************************************/

matrix *Hamming_Dist(allseq *data, model *mod)
{
  int i,j,k;
  phydbl unc_len;
  matrix *mat;
  phydbl **len;


  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  For(i,data->n_otu)
    len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));

  unc_len = .0;

  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);

  For(i,data->c_seq[0]->len)
    {
      For(j,data->n_otu-1)
	{
	  for(k=j+1;k<data->n_otu;k++)
	    {
	      if((!Is_Ambigu(data->c_seq[j]->state+i,mod->datatype,mod->stepsize)) &&
		 (!Is_Ambigu(data->c_seq[k]->state+i,mod->datatype,mod->stepsize)))
		{
		  len[j][k]+=data->wght[i];
		  len[k][j]=len[j][k];
		  if(data->c_seq[j]->state[i] != data->c_seq[k]->state[i])
		    mat->P[j][k]+=data->wght[i];
		}
	    }
	}
    }


  For(i,data->n_otu-1)
    for(j=i+1;j<data->n_otu;j++)
      {
	if(len[i][j])
	  {
	    mat->P[i][j] /= len[i][j];
	  }
	else
	  {
	    mat->P[i][j] = 1.;
	  }

	mat->P[j][i] = mat->P[i][j];

	mat->dist[i][j] = mat->P[i][j];


	if(mat->dist[i][j] > DIST_MAX)
	  {
	    mat->dist[i][j] = DIST_MAX;
	  }
	mat->dist[j][i] = mat->dist[i][j];
      }

  For(i,data->n_otu) free(len[i]);
  free(len);

  return mat;
}

/*********************************************************/
/* Test if the given site pattern is invariant. Does not handle ambiguities */

int Is_Invar(int patt_num, int stepsize, int datatype, allseq *data)
{
  int i, j;

  For(i,data->n_otu)
    {
      For(j,data->n_otu)
	{
	  if(!(Are_Compatible(data->c_seq[i]->state+patt_num,
			      data->c_seq[j]->state+patt_num,
			      stepsize,
			      datatype))) 
	    {
	      break;
	    }
	}
      if(j != data->n_otu) break;
    }
  
  if(i == data->n_otu) return 1;
  else                 return 0;
}


/*********************************************************/

int Is_Ambigu(char *state, int datatype, int stepsize)
{
  int i;

  if(datatype == NT)
    {
      For(i,stepsize)
	{
	  if(strchr("MRWSYKBDHVNXO?-.",state[i]))
	    return 1;
	}
    }
  else
    {
      if(strchr("X?-.",state[0])) return 1;
    }

  return 0;
}

/*********************************************************/

void Check_Ambiguities(allseq *data, int datatype, int stepsize)
{
  int i,j;

  Fors(j,data->crunch_len,stepsize) 
    {
      For(i,data->n_otu)
	{
	  data->ambigu[j]              = 0;
	  data->c_seq[i]->is_ambigu[j] = 0;
	}

      For(i,data->n_otu)
	{
	  if(Is_Ambigu(data->c_seq[i]->state+j,
		       datatype,
		       stepsize))
	    {
	      data->ambigu[j]              = 1;
	      data->c_seq[i]->is_ambigu[j] = 1;
	    }
	}
    }
}

/*********************************************************/

int Get_State_From_Ui(int ui, int datatype)
{
  if(datatype == NT)
    {
      switch(ui)
	{
	case 1 : {return 0; break;}
	case 2 : {return 1; break;}
	case 4 : {return 2; break;}
	case 8 : {return 3; break;}
	default : 
	  {
	    printf("\n. ui=%d",ui);
	    printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	    break;
	  }
	}
    }
  else if(datatype == AA)
    {
      switch(ui)
	{
	case 1 :      {return 0;  break;}
	case 2 :      {return 1;  break;} 
	case 4 :      {return 2;  break;} 
	case 8 :      {return 3;  break;} 
	case 16 :     {return 4;  break;} 
	case 32 :     {return 5;  break;} 
	case 64 :     {return 6;  break;} 
	case 128 :    {return 7;  break;} 
	case 256 :    {return 8;  break;} 
	case 512 :    {return 9;  break;} 
	case 1024 :   {return 10; break;} 
	case 2048 :   {return 11; break;} 
	case 4096 :   {return 12; break;} 
	case 8192 :   {return 13; break;} 
	case 16384 :  {return 14; break;} 
	case 32768 :  {return 15; break;} 
	case 65536 :  {return 16; break;} 
	case 131072 : {return 17; break;} 
	case 262144 : {return 18; break;} 
	case 524288 : {return 19; break;} 
	default : 
	  {
	    printf("\n. ui=%d",ui);
	    printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }
	}
    }
  else
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  return -1;
}

/*********************************************************/

int Assign_State(char *c, int datatype, int stepsize)
{
  int state[3];
  int i;

  state[0] = state[1] = state[2] = -1;
  if(datatype == NT)
    {
      For(i,stepsize)
	{
	  switch(c[i])
	    {
	    case 'A' : {state[i]=0;  break;}
	    case 'C' : {state[i]=1;  break;}
	    case 'G' : {state[i]=2;  break;}
	    case 'T' : {state[i]=3;  break;}
	    case 'U' : {state[i]=3;  break;}
	    default  : {state[i]=-1; break;}
	    }
	}
      return (stepsize>1)?(state[0]*16+state[1]*4+state[2]):(state[0]);
    }
  else
    {
      switch(c[0])
	{
	case 'A' : {state[0]=0 ; break;}
	case 'R' : {state[0]=1 ; break;}
	case 'N' : {state[0]=2 ; break;}
	case 'D' : {state[0]=3 ; break;}
	case 'C' : {state[0]=4 ; break;}
	case 'Q' : {state[0]=5 ; break;}
	case 'E' : {state[0]=6 ; break;}
	case 'G' : {state[0]=7 ; break;}
	case 'H' : {state[0]=8 ; break;}
	case 'I' : {state[0]=9 ; break;}
	case 'L' : {state[0]=10; break;}
	case 'K' : {state[0]=11; break;}
	case 'M' : {state[0]=12; break;}
	case 'F' : {state[0]=13; break;}
	case 'P' : {state[0]=14; break;}
	case 'S' : {state[0]=15; break;}
	case 'T' : {state[0]=16; break;}
	case 'W' : {state[0]=17; break;}
	case 'Y' : {state[0]=18; break;}
	case 'V' : {state[0]=19; break;}

	case 'B' : {state[0] = 2; break;}
	case 'Z' : {state[0] = 5; break;}
	default  : {state[0]=-1;  break;}
	}
      return state[0];
    }
  return -1;
}

/*********************************************************/

char Reciproc_Assign_State(int i_state, int datatype)
{
  
  if(datatype == NT)
    {
      i_state = i_state%4;
      switch(i_state)
	{
	case 0 :   {return 'A';  break;}
	case 1 :   {return 'C';  break;}
	case 2 :   {return 'G';  break;}
	case 3 :   {return 'T';  break;}
	default  : 
	  {
	    printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	    break;
	  }
	}
    }
  else
    {
      i_state = i_state%20;
      switch(i_state)
	{
	case 0  : {return 'A' ; break;}
	case 1  : {return 'R' ; break;}
	case 2  : {return 'N' ; break;}
	case 3  : {return 'D' ; break;}
	case 4  : {return 'C' ; break;}
	case 5  : {return 'Q' ; break;}
	case 6  : {return 'E' ; break;}
	case 7  : {return 'G' ; break;}
	case 8  : {return 'H' ; break;}
	case 9  : {return 'I' ; break;}
	case 10 : {return 'L';  break;}
	case 11 : {return 'K';  break;}
	case 12 : {return 'M';  break;}
	case 13 : {return 'F';  break;}
	case 14 : {return 'P';  break;}
	case 15 : {return 'S';  break;}
	case 16 : {return 'T';  break;}
	case 17 : {return 'W';  break;}
	case 18 : {return 'Y';  break;}
	case 19 : {return 'V';  break;}
	default  : 
	  {
	    printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	    break;
	  }
	}
    }
  return -1;
}

/*********************************************************/

int Assign_State_With_Ambiguity(char *c, int datatype, int stepsize)
{
  int state[3];
  int i;

  state[0] = state[1] = state[2] = -1;
  if(datatype == NT)
    {
      For(i,stepsize)
	{
	  switch(c[i])
	    {
	    case 'A' : {state[i]= 0;  break;}
	    case 'C' : {state[i]= 1;  break;}
	    case 'G' : {state[i]= 2;  break;}
	    case 'T' : {state[i]= 3;  break;}
	    case 'U' : {state[i]= 3;  break;}
	    case 'M' : {state[i]= 4;  break;}
	    case 'R' : {state[i]= 5;  break;}
	    case 'W' : {state[i]= 6;  break;}
	    case 'S' : {state[i]= 7;  break;}
	    case 'Y' : {state[i]= 8;  break;}
	    case 'K' : {state[i]= 9;  break;}
	    case 'B' : {state[i]=10;  break;}
	    case 'D' : {state[i]=11;  break;}
	    case 'H' : {state[i]=12;  break;}
	    case 'V' : {state[i]=13;  break;}
	    case 'N' : case 'X' : case '?' : case 'O' : case '-' : {state[i]=14;  break;}
	    default :
	      {
		printf("\n. Unknown character state : %c\n",state[i]);
		Warn_And_Exit("\n. Init failed (check the data type)\n");
		break;
	      }
	    }
	  return (stepsize>1)?(state[0]*16+state[1]*4+state[2]):(state[0]);
	}
    }
  else
    {
      switch(c[0])
	{
	case 'A' : {state[0]= 0; break;}
	case 'R' : {state[0]= 1; break;}
	case 'N' : {state[0]= 2; break;}
	case 'D' : {state[0]= 3; break;}
	case 'C' : {state[0]= 4; break;}
	case 'Q' : {state[0]= 5; break;}
	case 'E' : {state[0]= 6; break;}
	case 'G' : {state[0]= 7; break;}
	case 'H' : {state[0]= 8; break;}
	case 'I' : {state[0]= 9; break;}
	case 'L' : {state[0]=10; break;}
	case 'K' : {state[0]=11; break;}
	case 'M' : {state[0]=12; break;}
	case 'F' : {state[0]=13; break;}
	case 'P' : {state[0]=14; break;}
	case 'S' : {state[0]=15; break;}
	case 'T' : {state[0]=16; break;}
	case 'W' : {state[0]=17; break;}
	case 'Y' : {state[0]=18; break;}
	case 'V' : {state[0]=19; break;}
	case 'B' : {state[0]= 2; break;}
	case 'Z' : {state[0]= 5; break;}
	case 'X' : case '?' : case '-' : {state[0]=20; break;}
	default  : 
	  {
	    printf("\n. Unknown character state : %c\n",state[0]);
	    Warn_And_Exit("\n. Init failed (check the data type)\n");
	    break;
	  }
	}
      return state[0];
    }
  return -1;
}

/*********************************************************/

void Clean_Tree_Connections(arbre *tree)
{

  int i;
  For(i,2*tree->n_otu-2)
    {
      tree->noeud[i]->v[0] = NULL;
      tree->noeud[i]->v[1] = NULL;
      tree->noeud[i]->v[2] = NULL;
      tree->noeud[i]->b[0] = NULL;
      tree->noeud[i]->b[1] = NULL;
      tree->noeud[i]->b[2] = NULL;
    }
}

/*********************************************************/

void Bootstrap(arbre *tree)
{
  int *site_num, n_site;
  int replicate,j,k;
  int position,init_len;
  phydbl buff;
  allseq *boot_data;
  arbre *boot_tree;
  model *boot_mod;
  matrix *boot_mat;
  char *s;
/*   phydbl rf; */

  tree->print_boot_val       = 1;
  tree->print_alrt_val       = 0;
  boot_tree                  = NULL;

  site_num = (int *)mCalloc(tree->data->init_len,sizeof(int));

  Alloc_Bip(tree);
  Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);

  n_site = 0;
  For(j,tree->data->crunch_len) For(k,tree->data->wght[j])
    {
      site_num[n_site] = j;
      n_site++;
    }

  boot_data = Copy_Cseq(tree->data, tree->data->crunch_len, tree->mod->ns);

  printf("\n\n. Non parametric bootstrap analysis \n\n");
  printf("  ["); 


  For(replicate,tree->mod->bootstrap)
    {
      For(j,boot_data->crunch_len) boot_data->wght[j] = 0;

      init_len = 0;
      For(j,boot_data->init_len)
	{
	  buff  = rand();
	  buff /= RAND_MAX;
	  buff *= (phydbl)(tree->data->init_len-1.0);
	  if(buff-(int)(buff) > 0.5-MDBL_MAX) position = (int)(buff)+1;
	  else position = (int)(buff);
	  boot_data->wght[site_num[position]] += 1;
	  init_len++;
	}
      
      if(init_len != tree->data->init_len) Warn_And_Exit("\n. Pb when copying sequences\n");

      (tree->mod->datatype == NT)?
	(Get_Base_Freqs(boot_data)):
	(Get_AA_Freqs(boot_data));

      if(tree->io->random_boot_seq_order) Randomize_Sequence_Order(boot_data);

      boot_mod = Copy_Model(tree->mod);
      Init_Model(boot_data,boot_mod);

      if(tree->io->in_tree)
	{
	  rewind(tree->io->fp_in_tree);
	  boot_tree = Read_Tree_File(tree->io->fp_in_tree);
	}
      else
	{
	  boot_mat = ML_Dist(boot_data,boot_mod);
	  boot_mat->tree = Make_Tree_From_Scratch(boot_data->n_otu,boot_data);
	  Fill_Missing_Dist(boot_mat);
	  Bionj(boot_mat);
	  boot_tree = boot_mat->tree;
	  boot_tree->mat = boot_mat;
	}

      boot_tree->mod                = boot_mod;
      boot_tree->io                 = tree->io;
      boot_tree->data               = boot_data;
      boot_tree->both_sides         = 1;
      boot_tree->mod->s_opt->print  = 0;
      boot_tree->n_pattern          = boot_tree->data->crunch_len/
	                              boot_tree->mod->stepsize;
      boot_tree->io->print_site_lnl = 0;
      boot_tree->io->print_trace    = 0;

      if((boot_tree->mod->s_opt->random_input_tree) && (boot_tree->mod->s_opt->topo_search == SPR_MOVE)) Random_Tree(boot_tree);
      Order_Tree_CSeq(boot_tree,boot_data);
      Share_Lk_Struct(tree,boot_tree);
      Share_Spr_Struct(tree,boot_tree);
      Share_Pars_Struct(tree,boot_tree);
      Fill_Dir_Table(boot_tree);
      Update_Dirs(boot_tree);
      if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(boot_tree);
      else                         Init_P_Lk_Tips_Int(boot_tree);
      Init_Ui_Tips(boot_tree);
      Init_P_Pars_Tips(boot_tree); 
      Br_Len_Not_Involving_Invar(boot_tree);

      if(boot_tree->mod->s_opt->opt_topo)
	{
	  if(boot_tree->mod->s_opt->topo_search == NNI_MOVE) 
	    {
	      Simu_Loop(boot_tree);
	    }
	  else
	    {
	      if(tree->mod->s_opt->steph_spr)
		{
		  Speed_Spr_Loop(boot_tree);
		}
	      else
		{
		  Init_SPR (boot_tree);
		  Optim_SPR (boot_tree,0,ALL);
		  Clean_SPR (boot_tree);
		}
	    }
	}
      else
	{
	  if(boot_tree->mod->s_opt->opt_num_param || boot_tree->mod->s_opt->opt_bl)
	    Round_Optimize(boot_tree,boot_tree->data);
	  else
	    Lk(boot_tree);
	}

      Alloc_Bip(boot_tree);

      Get_Bip(boot_tree->noeud[0],
	      boot_tree->noeud[0]->v[0],
	      boot_tree);

      if(!tree->io->collapse_boot) 
	Compare_Bip(tree,boot_tree);
      else
	Compare_Bip_On_Existing_Edges(1,tree,boot_tree);

      Br_Len_Involving_Invar(boot_tree);

      if(tree->io->print_boot_trees)
	{
	  s = Write_Tree(boot_tree);
	  fprintf(tree->io->fp_out_boot_tree,"%s\n",s);
	  Free(s);
          Print_Fp_Out_Lines(tree->io->fp_out_boot_stats,0,0,boot_tree,tree->io,replicate+1);
	}


      /*       rf = .0; */
      /*       For(j,2*tree->n_otu-3)  */
      /* 	rf += tree->t_edges[j]->bip_score; */


      printf("."); 
#ifndef QUIET
fflush(stdout);
#endif
      if(!((replicate+1)%20))
	{
	  printf("] %4d/%4d\n  ",replicate+1,tree->mod->bootstrap);
	  if(replicate != tree->mod->bootstrap-1) printf("[");
	}

      if(boot_tree->mat) Free_Mat(boot_tree->mat);
      Free_Tree(boot_tree);      
      Free_Model(boot_mod);
    }

  if(((replicate)%20)) printf("] %4d/%4d\n ",replicate,tree->mod->bootstrap);

  tree->lock_topo = 1; /* Topology should not be modified afterwards */

  if(tree->io->print_boot_trees)
    {
      fclose(tree->io->fp_out_boot_tree);
      fclose(tree->io->fp_out_boot_stats);
    }

  Free_Cseq(boot_data);
  Free(site_num);
}

/*********************************************************/

void Br_Len_Involving_Invar(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l *= (1.0-tree->mod->pinvar);
}

/*********************************************************/

void Br_Len_Not_Involving_Invar(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l /= (1.0-tree->mod->pinvar);
}

/*********************************************************/

void Getstring_Stdin(char *file_name)
{
  fgets(file_name,T_MAX_LINE,stdin);
  if (strchr(file_name, '\n') != NULL)
    *strchr(file_name, '\n') = '\0';
}

/*********************************************************/

void Print_Freq(arbre *tree)
{

  switch(tree->mod->datatype)
    {
    case NT:
      {
	printf("A : %f\n",tree->mod->pi[0]);
	printf("C : %f\n",tree->mod->pi[1]);
	printf("G : %f\n",tree->mod->pi[2]);
	printf("T : %f\n",tree->mod->pi[3]);

	printf("U : %f\n",tree->mod->pi[4]);
	printf("M : %f\n",tree->mod->pi[5]);
	printf("R : %f\n",tree->mod->pi[6]);
	printf("W : %f\n",tree->mod->pi[7]);
	printf("S : %f\n",tree->mod->pi[8]);
	printf("Y : %f\n",tree->mod->pi[9]);
	printf("K : %f\n",tree->mod->pi[10]);
	printf("B : %f\n",tree->mod->pi[11]);
	printf("D : %f\n",tree->mod->pi[12]);
	printf("H : %f\n",tree->mod->pi[13]);
	printf("V : %f\n",tree->mod->pi[14]);
	printf("N : %f\n",tree->mod->pi[15]);
	break;
      }
    case AA:
      {
	printf("A : %f\n",tree->mod->pi[0]);
	printf("R : %f\n",tree->mod->pi[1]);
	printf("N : %f\n",tree->mod->pi[2]);
	printf("D : %f\n",tree->mod->pi[3]);
	printf("C : %f\n",tree->mod->pi[4]);
	printf("Q : %f\n",tree->mod->pi[5]);
	printf("E : %f\n",tree->mod->pi[6]);
	printf("G : %f\n",tree->mod->pi[7]);
	printf("H : %f\n",tree->mod->pi[8]);
	printf("I : %f\n",tree->mod->pi[9]);
	printf("L : %f\n",tree->mod->pi[10]);
	printf("K : %f\n",tree->mod->pi[11]);
	printf("M : %f\n",tree->mod->pi[12]);
	printf("F : %f\n",tree->mod->pi[13]);
	printf("P : %f\n",tree->mod->pi[14]);
	printf("S : %f\n",tree->mod->pi[15]);
	printf("T : %f\n",tree->mod->pi[16]);
	printf("W : %f\n",tree->mod->pi[17]);
	printf("Y : %f\n",tree->mod->pi[18]);
	printf("V : %f\n",tree->mod->pi[19]);

	printf("N : %f\n",tree->mod->pi[20]);
	break;
      }
    default : {break;}
    }
}

/*********************************************************/

phydbl Num_Derivatives_One_Param(phydbl (*func)(arbre *tree), arbre *tree,
				 phydbl f0, phydbl *param, phydbl stepsize,
				 phydbl *err, int precise)
{
  int i,j;
  phydbl errt,fac,hh,**a,ans;
  int n_iter;
  a = (phydbl **)mCalloc(11,sizeof(phydbl *));
  For(i,11) a[i] = (phydbl *)mCalloc(11,sizeof(phydbl));


  n_iter = 10; /* */

  ans  = .0;

  if (stepsize == 0.0) Warn_And_Exit("\n. h must be nonzero in Dfridr.");

  hh=stepsize;

  if(!precise)
    {

      *param   = *param+hh;
      a[0][0]  = (*func)(tree);
      a[0][0]  -= f0;
      a[0][0]  /= hh;
      *param   = *param-hh;

      ans =  a[0][0];
    }
  else
    {
      *param   = *param+hh;
      a[0][0]  = (*func)(tree);
      /*   *param   = *param-2*hh; */
      /*   a[0][0] -= (*func)(tree); */
      /*   a[0][0] /= (2.0*hh); */
      /*   *param   = *param+hh; */
      a[0][0]  -= f0;
      a[0][0]  /= hh;
      *param   = *param-hh;


      *err=1e30;
      for(i=1;i<n_iter;i++)
	{
	  hh /= 1.4;

	  /*       *param   = *param+hh; */
	  /*       a[0][i]  = (*func)(tree); */
	  /*       *param   = *param-2*hh; */
	  /*       a[0][i] -= (*func)(tree); */
	  /*       a[0][i] /= (2.0*hh); */
	  /*       *param   = *param+hh; */


	  *param   = *param+hh;
	  a[0][i]  = (*func)(tree);
	  /*   *param   = *param-2*hh; */
	  /*   a[0][i] -= (*func)(tree); */
	  /*   a[0][i] /= (2.0*hh); */
	  /*   *param   = *param+hh; */
	  a[0][i]  -= f0;
	  a[0][i]  /= hh;
	  *param   = *param-hh;


	  fac=1.4*1.4;
	  for (j=1;j<=i;j++)
	    {
	      a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
	      fac=1.4*1.4*fac;

	      errt=MAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));

	      if (errt <= *err)
		{
		  *err=errt;
		  ans=a[j][i];
		}
	    }

	  if(fabs(a[i][i]-a[i-1][i-1]) >= 2.0*(*err)) break;
	}
    }
  For(i,11) Free(a[i]);
  Free(a);

  return ans;
}

/*********************************************************/

void Num_Derivative_Several_Param(arbre *tree, phydbl *param, int n_param, phydbl stepsize,
				  phydbl (*func)(arbre *tree), phydbl *derivatives)
{
  int i;
  phydbl err,f0;

  f0 = (*func)(tree);

  For(i,n_param)
    {
      derivatives[i] = Num_Derivatives_One_Param(func,
						 tree,
						 f0,
						 param+i,
						 stepsize,
						 &err,
						 0
						 );
    }
}

/*********************************************************/

int Compare_Two_States(char *state1, char *state2, int state_size)
{
  /* 1 the two states are identical */
  /* 0 the two states are different */
  int i;

  For(i,state_size) if(state1[i] != state2[i]) break;

  return (i==state_size)?(1):(0);
}

/*********************************************************/

void Copy_One_State(char *from, char *to, int state_size)
{
  int i;
  For(i,state_size) to[i] = from[i];
}

/*********************************************************/

model *Make_Model_Basic()
{
  model *mod;

  mod                     = (model *)mCalloc(1,sizeof(model));

  mod->modelname          = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  mod->custom_mod_string  = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
  mod->user_b_freq        = (phydbl *)mCalloc(T_MAX_OPTION,sizeof(phydbl));

  mod->rr                 = (phydbl *)mCalloc(6,sizeof(phydbl));
  mod->rr_val             = (phydbl *)mCalloc(6,sizeof(phydbl));
  mod->rr_num             = (int *)mCalloc(6,sizeof(int *));
  mod->n_rr_per_cat       = (int *)mCalloc(6,sizeof(int));
  mod->s_opt              = (optimiz *)Alloc_Optimiz();

  return mod;
}

/*********************************************************/

void Make_Model_Complete(model *mod)
{
  int i,j;

  mod->pi             = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  mod->gamma_r_proba  = (phydbl *)mCalloc(mod->n_catg,sizeof(phydbl));
  mod->gamma_rr       = (phydbl *)mCalloc(mod->n_catg,sizeof(phydbl));
  mod->pi_unscaled    = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));

  mod->Pij_rr   = (double***)mCalloc(mod->n_catg,sizeof(double **));
  
  For(i,mod->n_catg)
    {
      mod->Pij_rr[i] = (double **)mCalloc(mod->ns,sizeof(double *));
      For(j,mod->ns) mod->Pij_rr[i][j] = (double *)mCalloc(mod->ns,sizeof(double));
    }

  mod->qmat      = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
  mod->qmat_buff = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
  mod->eigen     = (eigen *)Make_Eigen_Struct(mod);
  
  if(mod->n_rr_branch)
    {
      mod->rr_branch   = (phydbl *)mCalloc(mod->n_rr_branch,sizeof(phydbl));
      mod->p_rr_branch = (phydbl *)mCalloc(mod->n_rr_branch,sizeof(phydbl));
    }
}

/*********************************************************/

void Copy_Dist(phydbl **cpy, phydbl **orig, int n)
{
  int i,j;
  For(i,n) For(j,n) cpy[i][j] = orig[i][j];
}

/*********************************************************/

model *Copy_Model(model *ori)
{
  int i;
  model *cpy;

  cpy = Make_Model_Basic();

  Copy_Optimiz(ori->s_opt,cpy->s_opt);

  cpy->ns      = ori->ns;
  cpy->n_catg  = ori->n_catg;
  cpy->n_catg  = ori->n_catg;

  Make_Model_Complete(cpy);

  cpy->datatype     = ori->datatype;
  cpy->n_otu        = ori->n_otu;
  cpy->alpha_old    = ori->alpha_old;
  cpy->kappa_old    = ori->alpha_old;
  cpy->lambda_old   = ori->lambda_old;
  cpy->pinvar_old   = ori->pinvar_old;
  cpy->whichmodel   = ori->whichmodel;
  cpy->seq_len      = ori->seq_len;
  cpy->update_eigen = ori->update_eigen;
  cpy->kappa        = ori->kappa;
  cpy->alpha        = ori->alpha;
  cpy->lambda       = ori->lambda;
  cpy->bootstrap    = ori->bootstrap;
  cpy->invar        = ori->invar;
  cpy->pinvar       = ori->pinvar;
  cpy->stepsize     = ori->stepsize;
  cpy->n_diff_rr    = ori->n_diff_rr;

  For(i,6) cpy->rr_num[i] = ori->rr_num[i];

  For(i,6)
    {
      cpy->rr_val[i]  = ori->rr_val[i];
      cpy->rr[i] = cpy->rr[i];
    }
  
  For(i,cpy->ns)
      {
	cpy->pi[i]          = ori->pi[i];
	cpy->user_b_freq[i] = ori->user_b_freq[i];
      }
  
  For(i,cpy->ns*cpy->ns) cpy->qmat[i] = ori->qmat[i];

  For(i,cpy->n_catg)
    {
      cpy->gamma_r_proba[i] = ori->gamma_r_proba[i];
      cpy->gamma_rr[i]      = ori->gamma_rr[i];
    }
  
#ifndef PHYML
  if(ori->m4mod) 
    {
      cpy->m4mod     = M4_Copy_M4_Model(ori, ori->m4mod);
      cpy->use_m4mod = ori->use_m4mod;
    }
#endif 

  cpy->eigen->size = ori->eigen->size;
  For(i,2*ori->ns)       cpy->eigen->space[i]       = ori->eigen->space[i];
  For(i,2*ori->ns)       cpy->eigen->space_int[i]   = ori->eigen->space_int[i];
  For(i,ori->ns)         cpy->eigen->e_val[i]       = ori->eigen->e_val[i];
  For(i,ori->ns)         cpy->eigen->e_val_im[i]    = ori->eigen->e_val_im[i];
  For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect[i]    = ori->eigen->r_e_vect[i];
  For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect[i]    = ori->eigen->r_e_vect[i];
  For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect_im[i] = ori->eigen->r_e_vect_im[i];
  For(i,ori->ns*ori->ns) cpy->eigen->l_e_vect[i]    = ori->eigen->l_e_vect[i];
  For(i,ori->ns*ori->ns) cpy->eigen->q[i]           = ori->eigen->q[i];
 
  return cpy;
}

/*********************************************************/

option *Make_Input()
{
  option* io               = (option *)mCalloc(1,sizeof(option));
  io->mod                  = (model *)Make_Model_Basic();
  io->nt_or_cd             = (char *)mCalloc(T_MAX_FILE,sizeof(char));

  io->in_seq_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->in_tree_file         = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_best_tree_file   = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_boot_tree_file   = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_boot_stats_file  = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_stats_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_tree_file        = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_trace_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_lk_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_ps_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  return io;
}

/*********************************************************/

void Set_Defaults_Input(option* io)
{
  io->fp_in_seq                  = NULL;
  io->fp_in_tree                 = NULL;
  io->fp_out_tree                = NULL;
  io->fp_out_best_tree           = NULL;
  io->fp_out_boot_tree           = NULL;
  io->fp_out_boot_stats          = NULL;
  io->fp_out_stats               = NULL;

  io->tree                       = NULL;
  io->mod->datatype              = 0;
  strcpy(io->nt_or_cd,"nucleotides");
  io->n_data_sets                = 1;
  io->interleaved                = 1;
  io->in_tree                    = 0;
  io->out_tree_file_open_mode    = 1;
  io->out_stats_file_open_mode   = 1;
  io->seq_len                    = -1;
  io->n_data_set_asked           = -1;
  io->print_boot_trees           = 1;
  io->n_gt                       = 1;
  io->ratio_test		 = 4;
  io->multigene                  = 0;
  io->config_multigene           = 0;
  io->curr_interface             = 0;
  io->r_seed                     = -1;
  io->collapse_boot              = 0;
  io->random_boot_seq_order      = 1;
  io->print_trace                = 0;
  io->print_site_lnl             = 0;
  io->m4_model                   = NO;
  io->rm_ambigu                  = 0;
}

/*********************************************************/

void Set_Defaults_Model(model *mod)
{
  int i;

  strcpy(mod->modelname,"HKY85");
  strcpy(mod->custom_mod_string,"000000");
  mod->whichmodel              = HKY85;
  mod->n_catg                  = 1;
  mod->n_catg                  = 1;
  mod->kappa                   = 4.0;
  mod->alpha                   = 1.0;
  mod->lambda                  = 1.0;
  mod->bootstrap               = 0;
  mod->invar                   = 0;
  mod->pinvar                  = 0.0;
  mod->stepsize                = 1;
  mod->ns                      = 4;
  mod->n_diff_rr               = 0;
  For(i,6) mod->rr_val[i]      = 1.0;
  For(i,4) mod->user_b_freq[i] = -1.;
  mod->m4mod                   = NULL;
  mod->use_m4mod               = 0;
  mod->n_rr_branch             = 0;
  mod->rr_branch_alpha         = 0.1;
}

/*********************************************************/

void Set_Defaults_Optimiz(optimiz *s_opt)
{
  s_opt->print                = 1;
  s_opt->last_opt             = 1;
  s_opt->opt_alpha            = 0;
  s_opt->opt_kappa            = 0;
  s_opt->opt_bl               = 1;
  s_opt->opt_lambda           = 0;
  s_opt->opt_pinvar           = 0;
  s_opt->opt_num_param        = 0;
  s_opt->opt_cov_delta        = 0;
  s_opt->opt_cov_alpha        = 0;
  s_opt->opt_cov_free_rates   = 0;
  s_opt->opt_rr               = 0;
  s_opt->init_lk              = UNLIKELY;
  s_opt->n_it_max             = 1000;
  s_opt->opt_topo             = 1;
  s_opt->topo_search          = NNI_MOVE;
  s_opt->random_input_tree    = 0;
  s_opt->n_rand_starts        = 5;
  s_opt->brent_it_max         = 500;
  s_opt->steph_spr            = 1;
  s_opt->user_state_freq      = 0;
  s_opt->min_diff_lk_local    = 1.E-05;
  s_opt->min_diff_lk_global   = 1.E-03;
  s_opt->spr_step_after_nnis  = 0;
  s_opt->p_moves_to_examine   = 0.1;
  s_opt->fast_nni             = 0;
  s_opt->greedy               = 0;
  s_opt->general_pars         = 0;
  s_opt->tree_size_mult       = 1;

  s_opt->wim_n_rgrft          = -1;
  s_opt->wim_n_globl          = -1;
  s_opt->wim_max_dist         = -1;
  s_opt->wim_n_optim          = -1;
  s_opt->wim_n_best           = -1;
  s_opt->wim_inside_opt       =  0;
}

/*********************************************************/

void Copy_Optimiz(optimiz *ori, optimiz *cpy)
{
  cpy->print                = ori->print;
  cpy->last_opt             = ori->last_opt;
  cpy->opt_alpha            = ori->opt_alpha;
  cpy->opt_kappa            = ori->opt_kappa;
  cpy->opt_bl               = ori->opt_bl;
  cpy->opt_lambda           = ori->opt_lambda;
  cpy->opt_pinvar           = ori->opt_pinvar;
  cpy->opt_cov_delta        = ori->opt_cov_delta;
  cpy->opt_cov_alpha        = ori->opt_cov_alpha;
  cpy->opt_num_param        = ori->opt_num_param;
  cpy->opt_cov_free_rates   = ori->opt_cov_free_rates;
  cpy->opt_rr               = ori->opt_rr;
  cpy->init_lk              = ori->init_lk;
  cpy->n_it_max             = ori->n_it_max;
  cpy->opt_topo             = ori->opt_topo;
  cpy->topo_search          = ori->topo_search;
  cpy->random_input_tree    = ori->random_input_tree;
  cpy->n_rand_starts        = ori->n_rand_starts;
  cpy->brent_it_max         = ori->brent_it_max;
  cpy->steph_spr            = ori->steph_spr;
  cpy->user_state_freq      = ori->user_state_freq;
  cpy->min_diff_lk_local    = ori->min_diff_lk_local;
  cpy->min_diff_lk_global   = ori->min_diff_lk_global;
  cpy->spr_step_after_nnis  = ori->spr_step_after_nnis;
  cpy->p_moves_to_examine   = ori->p_moves_to_examine;
  cpy->fast_nni             = ori->fast_nni;
  cpy->greedy               = ori->greedy;
  cpy->general_pars         = ori->general_pars;
  cpy->tree_size_mult       = ori->tree_size_mult;

  cpy->wim_n_rgrft          = ori->wim_n_rgrft;
  cpy->wim_n_globl          = ori->wim_n_globl;
  cpy->wim_max_dist         = ori->wim_max_dist;
  cpy->wim_n_optim          = ori->wim_n_optim;
  cpy->wim_n_best           = ori->wim_n_best;
  cpy->wim_inside_opt       = ori->wim_inside_opt;
}

/*********************************************************/

void Get_Bip(node *a, node *d, arbre *tree)
{
  if(d->tax)
    {
      d->bip_node[0][0] = d;
      d->bip_size[0]    = 1;
      return;
    }
  else
    {
      int i,j,k;
      int d_a;


      d_a = -1;

      For(i,3)
	{
	  if(d->v[i] != a)
	    Get_Bip(d,d->v[i],tree);
	  else d_a = i;
	}

      d->bip_size[d_a] = 0;
      For(i,3)
	if(d->v[i] != a)
	  {
	    For(j,3)
	      {
		if(d->v[i]->v[j] == d)
		  {
		    For(k,d->v[i]->bip_size[j])
		      {
			d->bip_node[d_a][d->bip_size[d_a]] = d->v[i]->bip_node[j][k];
			strcpy(d->bip_name[d_a][d->bip_size[d_a]],d->v[i]->bip_node[j][k]->name);
			d->bip_size[d_a]++;
		      }
		    break;
		  }
	      }
	  }

      qsort(d->bip_name[d_a],d->bip_size[d_a],sizeof(char *),Sort_String);

      For(i,3)
	if(a->v[i] == d)
	  {
	    a->bip_size[i] = 0;
	    For(j,tree->n_otu)
	      {
		For(k,d->bip_size[d_a])
		  {
		    if(d->bip_node[d_a][k] == tree->noeud[j])
		      break;
		  }

		if(k == d->bip_size[d_a])
		  {
		    a->bip_node[i][a->bip_size[i]] = tree->noeud[j];
		    strcpy(a->bip_name[i][a->bip_size[i]],tree->noeud[j]->name);
		    a->bip_size[i]++;
		  }
	      }

	    qsort(a->bip_name[i],a->bip_size[i],sizeof(char *),Sort_String);

	    if(a->bip_size[i] != tree->n_otu - d->bip_size[d_a])
	      {
		printf("%d %d \n",a->bip_size[i],tree->n_otu - d->bip_size[d_a]);
		Warn_And_Exit("\n. Problem in counting bipartitions \n");
	      }
	    break;
	  }
    }
}

/*********************************************************/

void Alloc_Bip(arbre *tree)
{
  int i,j,k;

  tree->has_bip = 1;

  For(i,2*tree->n_otu-2)
    {
      tree->noeud[i]->bip_size = (int *)mCalloc(3,sizeof(int));
      tree->noeud[i]->bip_node = (node ***)mCalloc(3,sizeof(node **));
      tree->noeud[i]->bip_name = (char ***)mCalloc(3,sizeof(char **));
      For(j,3)
	{
	  tree->noeud[i]->bip_node[j] =
	    (node **)mCalloc(tree->n_otu,sizeof(node *));

	  tree->noeud[i]->bip_name[j] =
	    (char **)mCalloc(tree->n_otu,sizeof(char *));

	  For(k,tree->n_otu)
	    tree->noeud[i]->bip_name[j][k] =
	    (char *)mCalloc(T_MAX_NAME,sizeof(char ));
	}
    }
}

/*********************************************************/

int Sort_Phydbl_Increase(const void *a, const void *b)
{
  if((*(phydbl *)(a)) <= (*(phydbl *)(b))) return -1;
  else return 1;
}

/*********************************************************/

int Sort_String(const void *a, const void *b)
{
  return(strcmp((*(const char **)(a)), (*(const char **)(b))));
}

/*********************************************************/

void Compare_Bip_On_Existing_Edges(int discard, arbre *tree1, arbre *tree2)
{
  int i,j,k;
  edge *b1,*b2;
  char **bip1,**bip2;
  int bip_size;


  For(i,2*tree1->n_otu-3)
    {
      if((!tree1->t_edges[i]->left->tax) &&
	 (!tree1->t_edges[i]->rght->tax))
	{
	  b1 = tree1->t_edges[i];

	  For(j,2*tree2->n_otu-3)
	    {
	      if((!tree2->t_edges[j]->left->tax) &&
		 (!tree2->t_edges[j]->rght->tax))
		{
		  b2 = tree2->t_edges[j];

		  if(MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]) ==
		     MIN(b2->left->bip_size[b2->l_r],b2->rght->bip_size[b2->r_l]))
		    {
		      bip_size = MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]);

		      if(b1->left->bip_size[b1->l_r] == b1->rght->bip_size[b1->r_l])
			{
			  if(b1->left->bip_name[b1->l_r][0][0] < b1->rght->bip_name[b1->r_l][0][0])
			    {
			      bip1 = b1->left->bip_name[b1->l_r];
			    }
			  else
			    {
			      bip1 = b1->rght->bip_name[b1->r_l];
			    }
			}
		      else if(b1->left->bip_size[b1->l_r] < b1->rght->bip_size[b1->r_l])
			{
			  bip1 = b1->left->bip_name[b1->l_r];
			}
		      else
			{
			  bip1 = b1->rght->bip_name[b1->r_l];
			}

		      if(b2->left->bip_size[b2->l_r] == b2->rght->bip_size[b2->r_l])
			{
			  if(b2->left->bip_name[b2->l_r][0][0] < b2->rght->bip_name[b2->r_l][0][0])
			    {
			      bip2 = b2->left->bip_name[b2->l_r];
			    }
			  else
			    {
			      bip2 = b2->rght->bip_name[b2->r_l];
			    }
			}
		      else if(b2->left->bip_size[b2->l_r] < b2->rght->bip_size[b2->r_l])
			{
			  bip2 = b2->left->bip_name[b2->l_r];
			}
		      else
			{
			  bip2 = b2->rght->bip_name[b2->r_l];
			}

		      if(bip_size == 1) Warn_And_Exit("\n. Problem in Compare_Bip\n");


		      For(k,bip_size)
			{
			  if(strcmp(bip1[k],bip2[k])) break;
			}

		      if(k == bip_size)
			{
			  if(!((discard) && (b2->l < .5 / (phydbl)tree2->data->init_len)))
			    {
/* 			      One_Pars_Step(b2,tree2); */
			      b2->bip_score++;
			      b1->bip_score++;
			    }

/* 			  if((b1->l > .5 / (phydbl)tree1->data->init_len) || */
/* 			     (b2->l > .5 / (phydbl)tree2->data->init_len))  */
/* 			    { */
/* 			      b1->bip_score++; */
/* 			      b2->bip_score++; */
/* 			    } */
			  break;
			}
		    }
		}
	    }
	}
    }
}

/*********************************************************/

void Compare_Bip(arbre *tree1, arbre *tree2)
{
  int i,j,k;
  edge *b1,*b2;
  char **bip1,**bip2;
  int bip_size;


  For(i,2*tree1->n_otu-3)
    {
      if((!tree1->t_edges[i]->left->tax) &&
	 (!tree1->t_edges[i]->rght->tax))
	{

	  b1 = tree1->t_edges[i];

	  For(j,2*tree2->n_otu-3)
	    {
	      if((!tree2->t_edges[j]->left->tax) &&
		 (!tree2->t_edges[j]->rght->tax))
		{

		  b2 = tree2->t_edges[j];

		  if(MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]) ==
		     MIN(b2->left->bip_size[b2->l_r],b2->rght->bip_size[b2->r_l]))
		    {
		      bip_size = MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]);

		      if(b1->left->bip_size[b1->l_r] == b1->rght->bip_size[b1->r_l])
			{
			  if(b1->left->bip_name[b1->l_r][0][0] < b1->rght->bip_name[b1->r_l][0][0])
			    {
			      bip1 = b1->left->bip_name[b1->l_r];
			    }
			  else
			    {
			      bip1 = b1->rght->bip_name[b1->r_l];
			    }
			}
		      else if(b1->left->bip_size[b1->l_r] < b1->rght->bip_size[b1->r_l])
			{
			  bip1 = b1->left->bip_name[b1->l_r];
			}
		      else
			{
			  bip1 = b1->rght->bip_name[b1->r_l];
			}


		      if(b2->left->bip_size[b2->l_r] == b2->rght->bip_size[b2->r_l])
			{
			  if(b2->left->bip_name[b2->l_r][0][0] < b2->rght->bip_name[b2->r_l][0][0])
			    {
			      bip2 = b2->left->bip_name[b2->l_r];
			    }
			  else
			    {
			      bip2 = b2->rght->bip_name[b2->r_l];
			    }
			}
		      else if(b2->left->bip_size[b2->l_r] < b2->rght->bip_size[b2->r_l])
			{
			  bip2 = b2->left->bip_name[b2->l_r];
			}
		      else
			{
			  bip2 = b2->rght->bip_name[b2->r_l];
			}

		      if(bip_size == 1) Warn_And_Exit("\n. Problem in Compare_Bip\n");


		      For(k,bip_size)
			{
			  if(strcmp(bip1[k],bip2[k])) break;
			}

		      if(k == bip_size)
			{
			  b1->bip_score++;
			  b2->bip_score++;
			  break;
			}
		    }
		}
	    }
	}
    }
}

/*********************************************************/

void Test_Multiple_Data_Set_Format(option *io)
{
  char *line;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  io->n_trees = 0;

  while(fgets(line,T_MAX_LINE,io->fp_in_tree)) if(strstr(line,";")) io->n_trees++;

  Free(line);

  if((io->mod->bootstrap > 1) && (io->n_trees > 1))
    Warn_And_Exit("\n. Bootstrap option is not allowed with multiple input trees !\n");

  rewind(io->fp_in_tree);

  return;
}

/*********************************************************/

int Are_Compatible(char *statea, char *stateb, int stepsize, int datatype)
{
  int i,j;
  char a,b;


  if(datatype == NT)
    {
      For(i,stepsize)
	{
	  a = statea[i];
	  For(j,stepsize)
	    {
	      b = stateb[j];

	      switch(a)
		{
		case 'A':
		  {
		    switch(b)
		      {
		      case 'A' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'G':
		  {
		    switch(b)
		      {
		      case 'G' :
		      case 'R' :
		      case 'S' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'C':
		  {
		    switch(b)
		      {
		      case 'C' :
		      case 'M' :
		      case 'S' :
		      case 'Y' :
		      case 'B' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'T':
		  {
		    switch(b)
		      {
		      case 'T' :
		      case 'W' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'X' :
			{b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'M' :
		  {
		    switch(b)
		      {
		      case 'M' :
		      case 'A' :
		      case 'C' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' :
			{b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'R' :
		  {
		    switch(b)
		      {
		      case 'R' :
		      case 'A' :
		      case 'G' :
		      case 'M' :
		      case 'W' :
		      case 'S' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }

		case 'W' :
		  {
		    switch(b)
		      {
		      case 'W' :
		      case 'A' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }

		case 'S' :
		  {
		    switch(b)
		      {
		      case 'S' :
		      case 'C' :
		      case 'G' :
		      case 'M' :
		      case 'R' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }

		case 'Y' :
		  {
		    switch(b)
		      {
		      case 'Y' :
		      case 'C' :
		      case 'T' :
		      case 'M' :
		      case 'W' :
		      case 'S' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }

		case 'K' :
		  {
		    switch(b)
		      {
		      case 'K' :
		      case 'G' :
		      case 'T' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'B' :
		  {
		    switch(b)
		      {
		      case 'B' :
		      case 'C' :
		      case 'G' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'D' :
		  {
		    switch(b)
		      {
		      case 'D' :
		      case 'A' :
		      case 'G' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'H' :
		  {
		    switch(b)
		      {
		      case 'H' :
		      case 'A' :
		      case 'C' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'V' :
		  {
		    switch(b)
		      {
		      case 'V' :
		      case 'A' :
		      case 'C' :
		      case 'G' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'X' :
		  {
		    switch(b)
		      {
		      case 'X' :
		      case 'A' :
		      case 'C' :
		      case 'G' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		default :
		  {
                      printf("\n. Err. in Are_Compatible\n");
                      printf("\n. Please check that characters `%c` and `%c`\n",a,b);
                      printf("  correspond to existing amino-acids.\n");
                      Warn_And_Exit("\n");
                      return 0;
		  }
		}
	    }
	}
    }
  else
    {
      a = statea[0]; b = stateb[0];
      switch(a)
	{
	case 'A' :
	  {
	    switch(b)
	      {
	      case 'A' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'R' :
	  {
	    switch(b)
	      {
	      case 'R' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'N' :
	  {
	    switch(b)
	      {
	      case 'N' :
	      case 'B' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'B' :
	  {
	    switch(b)
	      {
	      case 'N' :
	      case 'B' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'D' :
	  {
	    switch(b)
	      {
	      case 'D' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'C' :
	  {
	    switch(b)
	      {
	      case 'C' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Q' :
	  {
	    switch(b)
	      {
	      case 'Q' :
	      case 'Z' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Z' :
	  {
	    switch(b)
	      {
	      case 'Q' :
	      case 'Z' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'E' :
	  {
	    switch(b)
	      {
	      case 'E' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'G' :
	  {
	    switch(b)
	      {
	      case 'G' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'H' :
	  {
	    switch(b)
	      {
	      case 'H' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'I' :
	  {
	    switch(b)
	      {
	      case 'I' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'L' :
	  {
	    switch(b)
	      {
	      case 'L' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'K' :
	  {
	    switch(b)
	      {
	      case 'K' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'M' :
	  {
	    switch(b)
	      {
	      case 'M' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'F' :
	  {
	    switch(b)
	      {
	      case 'F' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'P' :
	  {
	    switch(b)
	      {
	      case 'P' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'S' :
	  {
	    switch(b)
	      {
	      case 'S' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'T' :
	  {
	    switch(b)
	      {
	      case 'T' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'W' :
	  {
	    switch(b)
	      {
	      case 'W' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Y' :
	  {
	    switch(b)
	      {
	      case 'Y' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'V' :
	  {
	    switch(b)
	      {
	      case 'V' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'X' :
	  {
	    switch(b)
	      {
	      case 'A':case 'R':case 'N' :case 'B' :case 'D' :
	      case 'C':case 'Q':case 'Z' :case 'E' :case 'G' :
	      case 'H':case 'I':case 'L' :case 'K' :case 'M' :
	      case 'F':case 'P':case 'S' :case 'T' :case 'W' :
	      case 'Y':case 'V': case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	default :
	  {
	    printf("\n. Err. in Are_Compatible\n");
            printf("\n. Please check that characters `%c` and `%c`\n",a,b);
            printf("  correspond to existing amino-acids.\n");
            Warn_And_Exit("\n");
	    return 0;
	  }
	}
    }
  return 1;
}

/*********************************************************/

void Hide_Ambiguities(allseq *data)
{
  int i;
  For(i,data->crunch_len) if(data->ambigu[i]) data->wght[i] = 0;
}

/*********************************************************/

void Copy_Tree_Topology_With_Labels(arbre *ori, arbre *cpy)
{
  int i,j;

  For(i,2*ori->n_otu-2)
    {
      For(j,3)
	{
	  if(ori->noeud[i]->v[j])
	    {
	      cpy->noeud[i]->v[j] = cpy->noeud[ori->noeud[i]->v[j]->num];
	      cpy->noeud[i]->l[j] = ori->noeud[i]->l[j];
	    }
	  else
	    cpy->noeud[i]->v[j] = NULL;
	}
      cpy->noeud[i]->num = ori->noeud[i]->num;
      cpy->noeud[i]->tax = 0;
    }

  For(i,2*ori->n_otu-3)
    {
      cpy->t_edges[i]->l = ori->t_edges[i]->l;
    }

  For(i,ori->n_otu)
    {
      cpy->noeud[i]->tax = 1;
      strcpy(cpy->noeud[i]->name,ori->noeud[i]->name);
    }

}

/*********************************************************/

void Prune_Subtree(node *a, node *d, edge **target, edge **residual, arbre *tree)
{
  node *v1, *v2;
  edge *b1, *b2;
  int dir_v1, dir_v2;
  int i;
  phydbl ***buff_p_lk,*buff_scale;
  int **buff_p_pars, *buff_pars;
  unsigned int *buff_ui;
  short int **buff_p_lk_tip;

  if(a->tax)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  dir_v1 = dir_v2 = -1;
  For(i,3)
    {
      if(a->v[i] != d)
	{
	  if(dir_v1 < 0) dir_v1 = i;
	  else           dir_v2 = i;
	}
    }

  if(a->v[dir_v1]->num < a->v[dir_v2]->num)
    {
      v1 = a->v[dir_v1];
      v2 = a->v[dir_v2];
      b1 = a->b[dir_v1];
      b2 = a->b[dir_v2];
    }
  else
    {
      v1 = a->v[dir_v2];
      v2 = a->v[dir_v1];
      b1 = a->b[dir_v2];
      b2 = a->b[dir_v1];
    }


  a->v[dir_v1] = NULL;
  a->v[dir_v2] = NULL;
  a->b[dir_v1] = NULL;
  a->b[dir_v2] = NULL;


  if(v1 == b1->left)
    {
      b1->rght = v2;

      if(v2 == b2->left)
	{
	  buff_p_lk            = b1->p_lk_rght;
	  b1->p_lk_rght        = b2->p_lk_left;
	  b2->p_lk_left        = buff_p_lk;

	  buff_p_lk_tip        = b1->p_lk_tip_r;
	  b1->p_lk_tip_r       = b2->p_lk_tip_l;
	  b2->p_lk_tip_l       = buff_p_lk_tip;

	  buff_scale           = b1->sum_scale_f_rght;
	  b1->sum_scale_f_rght = b2->sum_scale_f_left;
	  b2->sum_scale_f_left = buff_scale;

	  buff_pars            = b1->pars_r;
	  b1->pars_r           = b2->pars_l;
	  b2->pars_l           = buff_pars;

	  buff_ui              = b1->ui_r;
	  b1->ui_r             = b2->ui_l;
	  b2->ui_l             = buff_ui;

	  buff_p_pars          = b1->p_pars_r;
	  b1->p_pars_r         = b2->p_pars_l;
	  b2->p_pars_l         = buff_p_pars;
	}
      else
	{
	  buff_p_lk            = b1->p_lk_rght; /* b1->p_lk_rght = NULL if b1->rght->tax */
	  b1->p_lk_rght        = b2->p_lk_rght; /* b2->p_lk_rght = NULL if b2->rght->tax */ 
	  b2->p_lk_rght        = buff_p_lk;

	  buff_p_lk_tip        = b1->p_lk_tip_r;
	  b1->p_lk_tip_r       = b2->p_lk_tip_r;
	  b2->p_lk_tip_r       = buff_p_lk_tip;

	  buff_scale           = b1->sum_scale_f_rght;
	  b1->sum_scale_f_rght = b2->sum_scale_f_rght;
	  b2->sum_scale_f_rght = buff_scale;

	  buff_pars            = b1->pars_r;
	  b1->pars_r           = b2->pars_r;
	  b2->pars_r           = buff_pars;

	  buff_ui              = b1->ui_r;
	  b1->ui_r             = b2->ui_r;
	  b2->ui_r             = buff_ui;

	  buff_p_pars          = b1->p_pars_r;
	  b1->p_pars_r         = b2->p_pars_r;
	  b2->p_pars_r         = buff_p_pars;
	}
    }
  else
    {
      b1->left = v2;

      if(v2 == b2->left)
	{
	  buff_p_lk            = b1->p_lk_left;
	  b1->p_lk_left        = b2->p_lk_left;
	  b2->p_lk_left        = buff_p_lk;

	  buff_p_lk_tip        = b1->p_lk_tip_l;
	  b1->p_lk_tip_l       = b2->p_lk_tip_l;
	  b2->p_lk_tip_l       = buff_p_lk_tip;

	  buff_scale           = b1->sum_scale_f_left;
	  b1->sum_scale_f_left = b2->sum_scale_f_left;
	  b2->sum_scale_f_left = buff_scale;

	  buff_pars            = b1->pars_l;
	  b1->pars_l           = b2->pars_l;
	  b2->pars_l           = buff_pars;

	  buff_ui              = b1->ui_l;
	  b1->ui_l             = b2->ui_l;
	  b2->ui_l             = buff_ui;

	  buff_p_pars          = b1->p_pars_l;
	  b1->p_pars_l         = b2->p_pars_l;
	  b2->p_pars_l         = buff_p_pars;
	}
      else
	{
	  buff_p_lk            = b1->p_lk_left;
	  b1->p_lk_left        = b2->p_lk_rght; /* b2->p_lk_rght = NULL if b2->rght->tax */
	  b2->p_lk_rght        = buff_p_lk;

	  buff_p_lk_tip        = b1->p_lk_tip_l;
	  b1->p_lk_tip_l       = b2->p_lk_tip_r;
	  b2->p_lk_tip_r       = buff_p_lk_tip;

	  buff_scale           = b1->sum_scale_f_left;
	  b1->sum_scale_f_left = b2->sum_scale_f_rght;
	  b2->sum_scale_f_rght = buff_scale;

	  buff_pars            = b1->pars_l;
	  b1->pars_l           = b2->pars_r;
	  b2->pars_r           = buff_pars;

	  buff_ui              = b1->ui_l;
	  b1->ui_l             = b2->ui_r;
	  b2->ui_r             = buff_ui;

	  buff_p_pars          = b1->p_pars_l;
	  b1->p_pars_l         = b2->p_pars_r;
	  b2->p_pars_r         = buff_p_pars;
	}
    }

  For(i,3)
    if(v2->v[i] == a)
      {
	v2->v[i] = v1;
	v2->b[i] = b1;
	break;
      }

#ifdef DEBUG
  if(i == 3)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif

  For(i,3)
    if(v1->v[i] == a)
      {
	v1->v[i] = v2;
	break;
      }

#ifdef DEBUG
  if(i == 3)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif

  b1->l += b2->l;


  (v1 == b1->left)?
    (Make_Edge_Dirs(b1,v1,v2)):
    (Make_Edge_Dirs(b1,v2,v1));

  if(target)   (*target)   = b1;
  if(residual) (*residual) = b2;
}


/*********************************************************/

void Graft_Subtree(edge *target, node *link, edge *residual, arbre *tree)
{
  node *v1, *v2;
  int i, dir_v1, dir_v2;
  phydbl ***buff_p_lk, *buff_scale;
  int **buff_p_pars, *buff_pars; 
  short int **buff_p_lk_tip;
  unsigned int *buff_ui;
  edge *b_up;

  dir_v1 = dir_v2 = -1;
  b_up = NULL;
  For(i,3)
    {
      if(!link->v[i])
	{
	  if(dir_v1 < 0) dir_v1 = i;
	  else           dir_v2 = i;
	}
      else b_up = link->b[i];
    }

  if(target->left->num < target->rght->num)
    {
      v1                           = target->left;
      v2                           = target->rght;

      buff_p_lk                    = residual->p_lk_rght;
      residual->p_lk_rght          = target->p_lk_rght;
      target->p_lk_rght            = buff_p_lk;

      buff_p_lk_tip                = residual->p_lk_tip_r;
      residual->p_lk_tip_r         = target->p_lk_tip_r;
      target->p_lk_tip_r           = buff_p_lk_tip;

      buff_scale                   = residual->sum_scale_f_rght;
      residual->sum_scale_f_rght   = target->sum_scale_f_rght;
      target->sum_scale_f_rght     = buff_scale;

      buff_pars                    = residual->pars_r;
      residual->pars_r             = target->pars_r;
      target->pars_r               = buff_pars;

      buff_ui                      = residual->ui_r;
      residual->ui_r               = target->ui_r;
      target->ui_r                 = buff_ui;

      buff_p_pars                  = residual->p_pars_r;
      residual->p_pars_r           = target->p_pars_r;
      target->p_pars_r             = buff_p_pars;
    }
  else
    {
      v1                           = target->rght;
      v2                           = target->left;

      buff_p_lk                    = residual->p_lk_rght;
      residual->p_lk_rght          = target->p_lk_left;
      target->p_lk_left            = buff_p_lk;

      buff_p_lk_tip                = residual->p_lk_tip_r;
      residual->p_lk_tip_r         = target->p_lk_tip_l;
      target->p_lk_tip_l           = buff_p_lk_tip;

      buff_scale                   = residual->sum_scale_f_rght;
      residual->sum_scale_f_rght   = target->sum_scale_f_left;
      target->sum_scale_f_left     = buff_scale;

      buff_pars                    = residual->pars_r;
      residual->pars_r             = target->pars_l;
      target->pars_l               = buff_pars;

      buff_ui                      = residual->ui_r;
      residual->ui_r               = target->ui_l;
      target->ui_l                 = buff_ui;

      buff_p_pars                  = residual->p_pars_r;
      residual->p_pars_r           = target->p_pars_l;
      target->p_pars_l             = buff_p_pars;
    }

  For(i,3)
    if(v2->b[i] == target)
      {
	v2->v[i] = link;
	v2->b[i] = residual;
	break;
      }

  link->v[dir_v2] = v2;
  link->b[dir_v2] = residual;

  residual->left  = link;
  residual->rght  = v2;

  (v1 == target->left)?(target->rght = link):(target->left = link);

  link->v[dir_v1] = v1;
  link->b[dir_v1] = target;

  For(i,3)
    if(v1->v[i] == v2)
      {
	v1->v[i] = link;
	break;
      }

  target->l /= 2.;
  residual->l = target->l;

  Make_Edge_Dirs(target,target->left,target->rght);
  Make_Edge_Dirs(residual,residual->left,residual->rght);
  Make_Edge_Dirs(b_up,b_up->left,b_up->rght);
}

/*********************************************************/

void Pull_Subtree_From_Dead_Objects(node *a, node *d, arbre *tree)
{
  int i;

  For(i,3)
    if(a->v[i] == d)
      {
	tree->n_dead_nodes--;
	tree->n_dead_edges--;
#ifdef DEBUG
	  if((tree->n_dead_edges < 0) || (tree->n_dead_nodes < 0))
	    {
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
#endif
	break;
      }

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Pull_Subtree_From_Dead_Objects(d,d->v[i],tree);
	}
    }
}

/*********************************************************/

void Put_Subtree_In_Dead_Objects(node *a, node *d, arbre *tree)
{
  int i;


  For(i,3)
    {
      if(a->v[i] == d)
	{
#ifdef DEBUG
	  if((tree->n_dead_edges < 0) || (tree->n_dead_nodes < 0))
	    {
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
#endif
	  tree->t_dead_edges[tree->n_dead_edges] = a->b[i];
	  tree->t_dead_nodes[tree->n_dead_nodes] = d;
	  tree->n_dead_nodes++;
	  tree->n_dead_edges++;
	  break;
	}
    }

  if(d->tax)
    return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Put_Subtree_In_Dead_Objects(d,d->v[i],tree);
	}
    }
}

/*********************************************************/

void Reassign_Node_Nums(node *a, node *d, int *curr_ext_node, int *curr_int_node, arbre *tree)
{
  node *buff;
  int i;

  if(a->tax)
    {
      buff = tree->noeud[*curr_ext_node];
      tree->noeud[*curr_ext_node] = a;
      tree->noeud[a->num] = buff;
      buff->num = a->num;
      a->num = *curr_ext_node;
      (*curr_ext_node)++;
    }

  if(d->tax)
    {
      buff = tree->noeud[*curr_ext_node];
      tree->noeud[*curr_ext_node] = d;
      tree->noeud[d->num] = buff;
      buff->num = d->num;
      d->num = *curr_ext_node;
      (*curr_ext_node)++;
      return;
    }
  else
    {
      buff = tree->noeud[*curr_int_node];
      tree->noeud[*curr_int_node] = d;
      tree->noeud[d->num] = buff;
      buff->num = d->num;
      d->num = *curr_int_node;
      (*curr_int_node)++;
    }

  For(i,3)
    {
      if(d->v[i] != a)
	Reassign_Node_Nums(d,d->v[i],curr_ext_node,curr_int_node,tree);
    }
}

/*********************************************************/

void Reassign_Edge_Nums(node *a, node *d, int *curr_br, arbre *tree)
{
  edge *buff;
  int i,j;

  For(i,3)
    if(a->v[i] == d)
      {
	buff = tree->t_edges[*curr_br];
	For(j,2*N_MAX_OTU-3) if(tree->t_edges[j] == a->b[i]) break;
	if(j == 2*N_MAX_OTU-3)
	  {
	    printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }
	tree->t_edges[*curr_br] = a->b[i];
	tree->t_edges[j] = buff;
	a->b[i]->num = *curr_br;
	(*curr_br)++;
	break;
      }

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a)
	  Reassign_Edge_Nums(d,d->v[i],curr_br,tree);
    }
}

/*********************************************************/

void Make_List_Of_Reachable_Tips(arbre *tree)
{
  int i,j;

  For(i,2*tree->n_otu-2)
    {
      tree->noeud[i]->list_of_reachable_tips = (node ***)mCalloc(3,sizeof(node **));
      tree->noeud[i]->n_of_reachable_tips    = (int *)mCalloc(3,sizeof(int));
      For(j,3)
	tree->noeud[i]->list_of_reachable_tips[j] = (node **)mCalloc(tree->n_otu,sizeof(node *));
    }
}

/*********************************************************/

void Get_List_Of_Reachable_Tips(arbre *tree)
{
  int i,j;
  
  For(i,2*tree->n_otu-2)
    {
      tree->noeud[i]->n_of_reachable_tips[0] = 0;
      tree->noeud[i]->n_of_reachable_tips[1] = 0;
      tree->noeud[i]->n_of_reachable_tips[2] = 0;
      For(j,tree->n_otu)
	{
	  tree->noeud[i]->list_of_reachable_tips[0][j] = NULL;
	  tree->noeud[i]->list_of_reachable_tips[1][j] = NULL;
	  tree->noeud[i]->list_of_reachable_tips[2][j] = NULL;
	}
    }
  
  Get_List_Of_Reachable_Tips_Post(tree->noeud[0],
				  tree->noeud[0]->v[0],
				  tree);
  Get_List_Of_Reachable_Tips_Pre(tree->noeud[0],
				 tree->noeud[0]->v[0],
				 tree);
}

/*********************************************************/

void Get_List_Of_Reachable_Tips_Post(node *a, node *d, arbre *tree)
{
  int i,j,k,cpt;

  if(d->tax)
    {
      For(i,3)
	if(a->v[i] == d)
	  {
	    a->list_of_reachable_tips[i][0] = d;
	    a->n_of_reachable_tips[i]       = 1;
	    break;
	  }
      return;
    }
  else
    {
      For(i,3)
	if(d->v[i] != a)
	  Get_List_Of_Reachable_Tips_Post(d,d->v[i],tree);

      For(i,3)
	{
	  if(a->v[i] == d)
	    {
	      a->n_of_reachable_tips[i] = 0;
	      cpt                       = 0;
	      For(j,3)
		{
		  if(d->v[j] != a)
		    {
		      For(k,d->n_of_reachable_tips[j])
			{
			  a->list_of_reachable_tips[i][cpt] = d->list_of_reachable_tips[j][k];
			  a->n_of_reachable_tips[i]++;
			  cpt++;
			}
		    }
		}
	      break;
	    }
	}
    }
}

/*********************************************************/

void Get_List_Of_Reachable_Tips_Pre(node *a, node *d, arbre *tree)
{
  int i,j,k,cpt;

  For(i,3)
    {
      if(d->v[i] == a)
	{
	  if(a->tax)
	    {
	      d->list_of_reachable_tips[i][0] = a;
	      d->n_of_reachable_tips[i]       = 1;
	    }
	  else
	    {
	      d->n_of_reachable_tips[i] = 0;
	      cpt = 0;
	      For(j,3)
		{
		  if(a->v[j] != d)
		    {
		      For(k,a->n_of_reachable_tips[j])
			{
			  d->list_of_reachable_tips[i][cpt] = a->list_of_reachable_tips[j][k];
			  d->n_of_reachable_tips[i]++;
			  cpt++;
			}
		    }
		}
	    }
	  break;
	}
    }

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a)
	  Get_List_Of_Reachable_Tips_Pre(d,d->v[i],tree);

    }
}

/*********************************************************/

int Compare_List_Of_Reachable_Tips(node **list1, int size_list1, node **list2, int size_list2)
{
  int i,j,n_matches;

  n_matches = 0;
  For(i,size_list1)
    {
      For(j,size_list2)
	{
	  if(list1[i] == list2[j])
	    {
	      n_matches++;
	    }
	}
    }
  return n_matches;
}

/*********************************************************/

void Find_Mutual_Direction(node *n1, node *n2, int *dir_n1_to_n2, int *dir_n2_to_n1)
{
  int scores[3][3];
  int n_zero_line, n_zero_col;
  int i,j;

  For(i,3) For(j,3) scores[i][j] = 0;

  For(i,3)
    {
      For(j,3)
	{
	  scores[i][j] = Compare_List_Of_Reachable_Tips(n1->list_of_reachable_tips[i],
							n1->n_of_reachable_tips[i],
							n2->list_of_reachable_tips[j],
							n2->n_of_reachable_tips[j]);
	}
    }

  For(i,3)
    {
      n_zero_line = 0;
      For(j,3)
	{
	  if(!scores[i][j]) n_zero_line++;
	}
      if(n_zero_line != 2) {*dir_n1_to_n2 = i; break;}
    }


  For(i,3)
    {
      n_zero_col = 0;
      For(j,3)
	{
	  if(!scores[j][i]) n_zero_col++;
	}
      if(n_zero_col != 2) {*dir_n2_to_n1 = i; break;}
    }

}

/*********************************************************/

void Fill_Dir_Table(arbre *tree)
{
  int i,j,k,l;
  int found;

  Get_List_Of_Reachable_Tips(tree);

  For(i,tree->n_otu) For(j,2*tree->n_otu-2) tree->t_dir[i][j] = 0;

  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    For(j,tree->n_otu)
    {
      found = 0;
      For(k,3)
	{
	  For(l,tree->noeud[i]->n_of_reachable_tips[k])
	    {
	      if(tree->noeud[i]->list_of_reachable_tips[k][l] == tree->noeud[j])
		{
		  found = 1;
		  tree->t_dir[i][j] = k;
		  break;
		}
	    }
	  if(found) break;
	}
    }

  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    for(j=i;j<2*tree->n_otu-2;j++)
      {
	Find_Mutual_Direction(tree->noeud[i],tree->noeud[j],
			      &(tree->t_dir[i][j]),
			      &(tree->t_dir[j][i]));


      }
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/

int Get_Subtree_Size(node *a, node *d)
{
  int size,i;

  if(d->tax) return 1;
  else
    {
      size = 0;
      For(i,3)
	if(d->v[i] != a)
	  size += Get_Subtree_Size(d,d->v[i]);
    }
  return size;
}

/*********************************************************/


void Triple_Dist_Recur(node *a, node *d, arbre *tree)
{

  if(d->tax) return;
  else
    {
      int i;
      Triple_Dist(d,tree);

      For(i,3) if(d->v[i] != a)
	{
	  Update_P_Lk(tree,d->b[i],d);
	  Triple_Dist_Recur(d,d->v[i],tree);
	}
      For(i,3) if((d->v[i] == a) && !(d->v[i]->tax)) Update_P_Lk(tree,d->b[i],d);
    }
}

/*********************************************************/

void Fast_Br_Len_Recur(node *a, node *d, edge *b, arbre *tree)
{
  int i;

  Fast_Br_Len(b,tree);

  if(d->tax) return;
  else For(i,3) if(d->v[i] != a)
    {
      Update_P_Lk(tree,d->b[i],d);
      Fast_Br_Len_Recur(d,d->v[i],d->b[i],tree);
    }
  For(i,3) if((d->v[i] == a) && !(d->v[i]->tax)) Update_P_Lk(tree,d->b[i],d);

}

/*********************************************************/

void Fast_Br_Len(edge *b, arbre *tree)
{
  phydbl sum;
  phydbl ***prob, ****core, *F;
  int i, j, k, site;
  phydbl *pi;
  phydbl v_rght;


/*   Br_Len_Brent(10.*b->l,b->l,BL_MIN, */
/* 	       tree->mod->s_opt->min_diff_lk_local, */
/* 	       b,tree, */
/* 	       tree->mod->s_opt->brent_it_max); */


  core  = tree->triplet_struct->core;
  prob  = tree->triplet_struct->p_one_site;
  F     = tree->triplet_struct->F_bc;
  pi    = tree->triplet_struct->pi_bc;

  For(i,tree->mod->ns) pi[i] = tree->mod->pi[i];

  Update_PMat_At_Given_Edge(b,tree);

  For(i,tree->mod->ns) For(j,tree->mod->ns) For(k,tree->mod->n_catg)
    core[k][0][i][j] = b->Pij_rr[k][i][j]*tree->mod->pi[i]*tree->mod->gamma_r_proba[k];

  For(i,tree->mod->ns) For(j,tree->mod->ns) F[tree->mod->ns*i+j] = .0;

  For(site,tree->n_pattern)
    {
      For(i,tree->mod->ns) For(j,tree->mod->ns) prob[0][i][j] = .0;

      /* Joint probabilities of the states at the two ends of the edge */
      v_rght = -1.;
      For(i,tree->mod->ns)
	{
	  For(j,tree->mod->ns)
	    {
	      For(k,tree->mod->n_catg)
		{
		  v_rght = (b->rght->tax)?((phydbl)(b->p_lk_tip_r[site][j])):(b->p_lk_rght[site][k][j]);

		  prob[0][i][j]             +=
		    core[k][0][i][j]         *
		    b->p_lk_left[site][k][i] *
		    v_rght;
		}
	    }
	}

      sum = .0;
      For(i,tree->mod->ns) For(j,tree->mod->ns) sum += prob[0][i][j];

      /* Scaling */
      For(i,tree->mod->ns) For(j,tree->mod->ns) prob[0][i][j] /= sum;

      /* Expected number of each pair of states */
      For(i,tree->mod->ns) For(j,tree->mod->ns)
	F[tree->mod->ns*i+j] += tree->data->wght[site] * prob[0][i][j];
    }

  sum = .0;
  For(i,tree->mod->ns)
    {
      tree->mod->pi[i] = .0;
      For(j,tree->mod->ns)
	{
	  tree->mod->pi[i] += (F[tree->mod->ns*i+j] + F[tree->mod->ns*j+i])/2.;
	}
      tree->mod->pi[i] /= (phydbl)tree->data->init_len;
      sum += tree->mod->pi[i];
    }

  For(i,tree->mod->ns) tree->mod->pi[i] /= sum;

#ifdef DEBUG
  sum = .0;
  For(i,tree->mod->ns) sum += tree->mod->pi[i];
  if((int)rint(sum) != 1)
    {
      printf("\n. site %d sum = %f ",site,sum);
      printf("pi = ");
      For(i,tree->mod->ns) printf("%f ",tree->mod->pi[i]);
      printf("\n");
    }
#endif

  Divide_Cells(&F,(phydbl)tree->data->init_len,tree);
  Make_Symmetric(&F,tree->mod->ns);

#ifdef DEBUG
  phydbl lk_a, lk_b;
  if(b->l < BL_MIN) b->l = BL_MIN;
  else if(b->l > BL_MAX) b->l = BL_MAX;
  lk_b = Lk_Dist(F,b->l,tree->mod);
#endif

  Opt_Dist_F(&(b->l),F,tree->mod);
  
  if(b->l < BL_MIN) b->l = BL_MIN;
  else if(b->l > BL_MAX) b->l = BL_MAX;

#ifdef DEBUG
  lk_a = Lk_Dist(F,b->l,tree->mod);

  if(lk_b > lk_a + tree->mod->s_opt->min_diff_lk_local)
    {
      printf("\n. b->l = %f %d %f lk_b=%f lk_a=%f",
	     b->l,
	     tree->mod->n_catg,tree->mod->alpha,
	     lk_b,lk_a);
    }
#endif

  For(i,tree->mod->ns) tree->mod->pi[i] = pi[i];
}


/*********************************************************/

eigen *Make_Eigen_Struct(model *mod)
{
  eigen *eig;

  eig              = (eigen *)mCalloc(1,sizeof(eigen));
  eig->size        = mod->ns;
  eig->space       = (double *)mCalloc(2*mod->ns,sizeof(double));
  eig->space_int   = (int *)mCalloc(2*mod->ns,sizeof(int));
  eig->e_val       = (double *)mCalloc(mod->ns,sizeof(double));
  eig->e_val_im    = (double *)mCalloc(mod->ns,sizeof(double));
  eig->r_e_vect    = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
  eig->r_e_vect_im = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
  eig->l_e_vect    = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
  eig->q           = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));

  return eig;
}

/*********************************************************/

ttriplet *Make_Triplet_Struct(model *mod)
{
  int i,j,k;
  ttriplet *triplet_struct;

  triplet_struct                  = (ttriplet *)mCalloc(1,sizeof(ttriplet));
  triplet_struct->size            = mod->ns;
  triplet_struct->pi_bc           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->pi_cd           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->pi_bd           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->F_bc            = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl));
  triplet_struct->F_cd            = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl));
  triplet_struct->F_bd            = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl));
  triplet_struct->core            = (phydbl ****)mCalloc(mod->n_catg,sizeof(phydbl ***));
  triplet_struct->p_one_site      = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
  triplet_struct->sum_p_one_site  = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
  triplet_struct->eigen_struct    = (eigen *)Make_Eigen_Struct(mod);
  triplet_struct->mod             = mod;

  For(k,mod->n_catg)
    {
      triplet_struct->core[k]                = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
      For(i,mod->ns)
	{
	  triplet_struct->core[k][i]         = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
	  For(j,mod->ns)
	    triplet_struct->core[k][i][j]    = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
	}
    }

  For(i,mod->ns)
    {
      triplet_struct->p_one_site[i]          = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
      For(j,mod->ns)
	triplet_struct->p_one_site[i][j]     = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
    }

  For(i,mod->ns)
    {
      triplet_struct->sum_p_one_site[i]      = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
      For(j,mod->ns)
	triplet_struct->sum_p_one_site[i][j] = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
    }
  return triplet_struct;

}

/*********************************************************/

void Triple_Dist(node *a, arbre *tree)
{

  if(a->tax) return;
  else
    {
      Update_PMat_At_Given_Edge(a->b[1],tree);
      Update_PMat_At_Given_Edge(a->b[2],tree);
      Update_P_Lk(tree,a->b[0],a);
      Fast_Br_Len(a->b[0],tree);
      Update_PMat_At_Given_Edge(a->b[0],tree);

      Update_P_Lk(tree,a->b[1],a);
      Fast_Br_Len(a->b[1],tree);
      Update_PMat_At_Given_Edge(a->b[1],tree);

      Update_P_Lk(tree,a->b[2],a);
      Fast_Br_Len(a->b[2],tree);
      Update_PMat_At_Given_Edge(a->b[2],tree);

      Update_P_Lk(tree,a->b[1],a);
      Update_P_Lk(tree,a->b[0],a);


/*       node *b,*c,*d; */
/*       int _a,_b,_c,_d; */
/*       int dir_b,dir_c,dir_d; */
/*       int i, site, gamma; */
/*       phydbl ***p_lk_b,***p_lk_c,***p_lk_d; */
/*       double ***P_ab,***P_ac,***P_ad; */
/*       phydbl ***prob,***sum_prob,****core,*F_bc,*F_cd,*F_bd,*pi_bc,*pi_cd,*pi_bd; */
/*       phydbl sum, len; */
/*       phydbl d_bc, d_cd, d_bd; */
/*       double *eigen_val_real, *eigen_val_im, *eigen_vect_real, *eigen_vect_im, *space; */

/*       pi_bc           = tree->triplet_struct->pi_bc; */
/*       pi_cd           = tree->triplet_struct->pi_cd; */
/*       pi_bd           = tree->triplet_struct->pi_bd; */
/*       core            = tree->triplet_struct->core; */
/*       prob            = tree->triplet_struct->p_one_site; */
/*       sum_prob        = tree->triplet_struct->sum_p_one_site; */
/*       F_bc            = tree->triplet_struct->F_bc; */
/*       F_cd            = tree->triplet_struct->F_cd; */
/*       F_bd            = tree->triplet_struct->F_bd; */
/*       eigen_val_real  = tree->triplet_struct->eigen_struct->eigen_val_real; */
/*       eigen_val_im    = tree->triplet_struct->eigen_struct->eigen_val_im; */
/*       eigen_vect_real = tree->triplet_struct->eigen_struct->eigen_vect_real; */
/*       eigen_vect_im   = tree->triplet_struct->eigen_struct->eigen_vect_im; */
/*       space           = tree->triplet_struct->eigen_struct->space; */

/*       p_lk_b = p_lk_c = p_lk_d = NULL; */
/*       P_ab = P_ac = P_ad = NULL; */
/*       b = c = d = NULL; */
/*       dir_b = dir_c = dir_d = -1; */
/*       For(i,3) */
/* 	{ */
/* 	  if(!b) */
/* 	    { */
/* 	      b      = a->v[i]; */
/* 	      p_lk_b = (a == a->b[i]->left)?(a->b[i]->p_lk_rght):(a->b[i]->p_lk_left); */
/* 	      P_ab   = a->b[i]->Pij_rr; */
/* 	      dir_b  = i; */
/* 	    } */
/* 	  else if(!c) */
/* 	    { */
/* 	      c      = a->v[i]; */
/* 	      p_lk_c = (a == a->b[i]->left)?(a->b[i]->p_lk_rght):(a->b[i]->p_lk_left); */
/* 	      P_ac   = a->b[i]->Pij_rr; */
/* 	      dir_c  = i; */
/* 	    } */
/* 	  else if(!d) */
/* 	    { */
/* 	      d      = a->v[i]; */
/* 	      p_lk_d = (a == a->b[i]->left)?(a->b[i]->p_lk_rght):(a->b[i]->p_lk_left); */
/* 	      P_ad   = a->b[i]->Pij_rr; */
/* 	      dir_d  = i; */
/* 	    } */
/* 	} */

/*       For(_b,tree->mod->ns) For(_c,tree->mod->ns) For(_d,tree->mod->ns) sum_prob[_b][_c][_d] = .0; */

/*       For(i,tree->mod->n_catg) */
/* 	{ */
/* 	  len = a->b[dir_b]->l * tree->mod->rr[i]; */
/* 	  if(len < BL_MIN) len = BL_MIN; */
/* 	  else if(len > BL_MAX) len = BL_MAX; */
/* 	  PMat(len,tree->mod,&(P_ab[i])); */
/* 	} */

/*       For(i,tree->mod->n_catg) */
/* 	{ */
/* 	  len = a->b[dir_c]->l * tree->mod->rr[i]; */
/* 	  if(len < BL_MIN) len = BL_MIN; */
/* 	  else if(len > BL_MAX) len = BL_MAX; */
/* 	  PMat(len,tree->mod,&(P_ac[i])); */
/* 	} */

/*       For(i,tree->mod->n_catg) */
/* 	{ */
/* 	  len = a->b[dir_d]->l * tree->mod->rr[i]; */
/* 	  if(len < BL_MIN) len = BL_MIN; */
/* 	  else if(len > BL_MAX) len = BL_MAX; */
/* 	  PMat(len,tree->mod,&(P_ad[i])); */
/* 	} */

/* /\*       PMat(a->b[dir_b]->l,tree->mod,&(P_ab[0])); *\/ */
/* /\*       PMat(a->b[dir_c]->l,tree->mod,&(P_ac[0])); *\/ */
/* /\*       PMat(a->b[dir_d]->l,tree->mod,&(P_ad[0])); *\/ */

/*       For(gamma,tree->mod->n_catg) */
/* 	{ */
/* 	  For(_b,tree->mod->ns) */
/* 	    { */
/* 	      For(_c,tree->mod->ns) */
/* 		{ */
/* 		  For(_d,tree->mod->ns) */
/* 		    { */
/* 		      core[gamma][_b][_c][_d] = .0; */
/* 		      For(_a,tree->mod->ns) */
/* 			{ */
/* 			  core[gamma][_b][_c][_d] += */
/* 			    tree->mod->r_proba[gamma] * */
/* 			    tree->mod->pi[_a]         * */
/* 			    P_ab[gamma][_a][_b]       * */
/* 			    P_ac[gamma][_a][_c]       * */
/* 			    P_ad[gamma][_a][_d]       ; */
/* 			} */
/* 		    } */
/* 		} */
/* 	    } */
/* 	} */

/*       For(site,tree->n_pattern) */
/* 	{ */
/* 	  For(_b,tree->mod->ns) For(_c,tree->mod->ns) For(_d,tree->mod->ns) prob[_b][_c][_d] = .0; */

/* 	  For(gamma,tree->mod->n_catg) */
/* 	    For(_b,tree->mod->ns) */
/* 	    For(_c,tree->mod->ns) */
/* 	    For(_d,tree->mod->ns) */
/* 	    prob[_b][_c][_d]          += */
/* 	    core[gamma][_b][_c][_d]   * */
/* 	    p_lk_b[site][gamma][_b]   * */
/* 	    p_lk_c[site][gamma][_c]   * */
/* 	    p_lk_d[site][gamma][_d]   ; */


/* 	  sum = .0; */
/* 	  For(_b,tree->mod->ns) For(_c,tree->mod->ns) For(_d,tree->mod->ns) sum += prob[_b][_c][_d]; */

/* 	  For(_b,tree->mod->ns) For(_c,tree->mod->ns) For(_d,tree->mod->ns) prob[_b][_c][_d] /= sum; */

/* 	  For(_b,tree->mod->ns) For(_c,tree->mod->ns) For(_d,tree->mod->ns) */
/* 	    sum_prob[_b][_c][_d] += */
/* 	    tree->data->wght[site]* */
/* 	    prob[_b][_c][_d]; */
/* 	} */


/*       For(_b,tree->mod->ns) */
/* 	{ */
/* 	  pi_bc[_b] = .0; */
/* 	  For(_c,tree->mod->ns) */
/* 	    { */
/* 	      F_bc[tree->mod->ns*_b+_c] = .0; */
/* 	      For(_d,tree->mod->ns) */
/* 		{ */
/* 		  F_bc[tree->mod->ns*_b+_c] += sum_prob[_b][_c][_d]; */
/* 		} */
/* 	      pi_bc[_b] += (F_bc[tree->mod->ns*_b+_c]+F_bc[tree->mod->ns*_c+_b])/2.; */
/* 	    } */
/* 	  pi_bc[_b] /= (phydbl)tree->data->init_len; */
/* 	} */

/*       For(_c,tree->mod->ns) */
/* 	{ */
/* 	  pi_cd[_c] = .0; */
/* 	  For(_d,tree->mod->ns) */
/* 	    { */
/* 	      F_cd[tree->mod->ns*_c+_d] = .0; */
/* 	      For(_b,tree->mod->ns) */
/* 		{ */
/* 		  F_cd[tree->mod->ns*_c+_d] += sum_prob[_b][_c][_d]; */
/* 		} */
/* 	      pi_cd[_c] += (F_cd[tree->mod->ns*_c+_d]+F_cd[tree->mod->ns*_d+_c])/2.; */
/* 	    } */
/* 	  pi_cd[_c] /= (phydbl)tree->data->init_len; */
/* 	} */

/*       For(_b,tree->mod->ns) */
/* 	{ */
/* 	  pi_bd[_b] = .0; */
/* 	  For(_d,tree->mod->ns) */
/* 	    { */
/* 	      F_bd[tree->mod->ns*_b+_d] = .0; */
/* 	      For(_c,tree->mod->ns) */
/* 		{ */
/* 		  F_bd[tree->mod->ns*_b+_d] += sum_prob[_b][_c][_d]; */
/* 		} */
/* 	      pi_bd[_b] += (F_bd[tree->mod->ns*_b+_d]+F_bd[tree->mod->ns*_d+_b])/2.; */
/* 	    } */
/* 	  pi_bd[_b] /= (phydbl)tree->data->init_len; */
/* 	} */

/*       Divide_Cells(&F_bc,(phydbl)tree->data->init_len,tree); */
/*       Divide_Cells(&F_cd,(phydbl)tree->data->init_len,tree); */
/*       Divide_Cells(&F_bd,(phydbl)tree->data->init_len,tree); */

/*       Make_Symmetric(&F_bc,tree->mod->ns); */
/*       Make_Symmetric(&F_cd,tree->mod->ns); */
/*       Make_Symmetric(&F_bd,tree->mod->ns); */

/*       d_bc = d_cd = d_bd = 0.1; */


/* /\*       Dist_F_Brent(&(d_bc),F_bc,BL_MIN,d_bc,BL_MAX,1.E-4,&(d_bc),tree->mod,100); *\/ */
/* /\*       Dist_F_Brent(&(d_cd),F_cd,BL_MIN,d_cd,BL_MAX,1.E-4,&(d_cd),tree->mod,100); *\/ */
/* /\*       Dist_F_Brent(&(d_bd),F_bd,BL_MIN,d_bd,BL_MAX,1.E-4,&(d_bd),tree->mod,100); *\/ */

/*       Opt_Dist_F(&(d_bc),F_bc,tree->mod); */
/*       Opt_Dist_F(&(d_cd),F_cd,tree->mod); */
/*       Opt_Dist_F(&(d_bd),F_bd,tree->mod); */
      


/* /\*       d_bc = GTR_Dist(F_bc,(tree->mod->n_catg > 1)?(tree->mod->alpha):(-1.),tree->triplet_struct->eigen_struct); *\/ */
/* /\*       d_bd = GTR_Dist(F_bd,(tree->mod->n_catg > 1)?(tree->mod->alpha):(-1.),tree->triplet_struct->eigen_struct); *\/ */
/* /\*       d_cd = GTR_Dist(F_cd,(tree->mod->n_catg > 1)?(tree->mod->alpha):(-1.),tree->triplet_struct->eigen_struct); *\/ */

/*       a->b[dir_b]->l = (d_bc-d_cd+d_bd)/2.; */
/*       a->b[dir_c]->l = (d_bc-d_bd+d_cd)/2.; */
/*       a->b[dir_d]->l = (d_bd-d_bc+d_cd)/2.; */

/*       if(a->b[dir_b]->l < BL_MIN) a->b[dir_b]->l = BL_MIN; */
/*       else if(a->b[dir_b]->l > BL_MAX) a->b[dir_b]->l = BL_MAX; */
/*       if(a->b[dir_c]->l < BL_MIN) a->b[dir_c]->l = BL_MIN; */
/*       else if(a->b[dir_c]->l > BL_MAX) a->b[dir_c]->l = BL_MAX; */
/*       if(a->b[dir_d]->l < BL_MIN) a->b[dir_d]->l = BL_MIN; */
/*       else if(a->b[dir_d]->l > BL_MAX) a->b[dir_d]->l = BL_MAX; */
    }
}


/*********************************************************/

void Make_Symmetric(phydbl **F, int size)
{
  int i,j;

  For(i,size)
    {
      for(j=i+1;j<size;j++)
	{
	  (*F)[size*i+j] = ((*F)[size*i+j] + (*F)[size*j+i])/2.;
	  (*F)[size*j+i] = (*F)[size*i+j];
	}
    }
}

/*********************************************************/

void Round_Down_Freq_Patt(phydbl **F, arbre *tree)
{
  int i,j;

  For(i,tree->mod->ns)
    {
      For(j,tree->mod->ns)
	{
	  (*F)[tree->mod->ns*i+j] = rint((*F)[tree->mod->ns*i+j]);
	}
    }
}

/*********************************************************/

phydbl Get_Sum_Of_Cells(phydbl *F, arbre *tree)
{
  int i,j;
  phydbl sum = .0;

  For(i,tree->mod->ns)
    For(j,tree->mod->ns)
    sum += F[tree->mod->ns*i+j];

  return sum;
}


/*********************************************************/

void Divide_Cells(phydbl **F, phydbl div, arbre *tree)
{
  int i,j;

  For(i,tree->mod->ns)
    For(j,tree->mod->ns)
    (*F)[tree->mod->ns*i+j] /= div;
}

/*********************************************************/

void Divide_Mat_By_Vect(phydbl **F, phydbl *vect, int size)
{
  int i,j;
  For(i,size)
    For(j,size)
    (*F)[size*i+j] = (*F)[size*i+j] / vect[j];
}

/*********************************************************/

void Multiply_Mat_By_Vect(phydbl **F, phydbl *vect, int size)
{
  int i,j;
  For(i,size)
    For(j,size)
    (*F)[size*i+j] = (*F)[size*i+j] * vect[j];
}

/*********************************************************/

int Check_Spr_Move_Validity(spr *this_spr_move, arbre *tree)
{
  int match;

  match = 0;
  Found_In_Subtree(this_spr_move->n_link,
		   this_spr_move->n_opp_to_link,
		   this_spr_move->b_target->left,
		   &match,
		   tree);

  if(match) return 0;
  else      return 1;
}

/*********************************************************/

void Found_In_Subtree(node *a, node *d, node *target, int *match, arbre *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      if(d == target) *match =  1;
      For(i,3)
	{
	  if(d->v[i] != a)
	    Found_In_Subtree(d,d->v[i],target,match,tree);
	}
    }
}

/*********************************************************/

void Get_List_Of_Target_Edges(node *a, node *d, edge **list, int *list_size, arbre *tree)
{
  int i;

  For(i,3)
    {
      if(a->v[i] && a->v[i] == d)
	{
	  list[*list_size] = a->b[i];
	  (*list_size)++;
	}
    }

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Get_List_Of_Target_Edges(d,d->v[i],list,list_size,tree);
	}
    }
}

/*********************************************************/

void Fix_All(arbre *tree)
{
  int i;

  tree->mod->pinvar_old = tree->mod->pinvar;
  tree->mod->alpha_old  = tree->mod->alpha;
  tree->mod->kappa_old  = tree->mod->kappa;
  tree->mod->lambda_old = tree->mod->lambda;

  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    {
      tree->noeud[i]->b[0]->l_old = tree->noeud[i]->b[0]->l;
      tree->noeud[i]->b[1]->l_old = tree->noeud[i]->b[1]->l;
      tree->noeud[i]->b[2]->l_old = tree->noeud[i]->b[2]->l;
    }
}

/*********************************************************/

void Record_Br_Len(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l_old = tree->t_edges[i]->l;
}

/*********************************************************/

void Restore_Br_Len(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l = tree->t_edges[i]->l_old;
}

/*********************************************************/

void Get_Dist_Btw_Edges(node *a, node *d, arbre *tree)
{
  int i;
  edge *b_fcus;

  b_fcus = NULL;
  For(i,3) if(a->v[i] == d) {b_fcus = a->b[i]; break;}

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a)
	  {
	    d->b[i]->topo_dist_btw_edges = b_fcus->topo_dist_btw_edges + 1;
	    d->b[i]->dist_btw_edges      = b_fcus->dist_btw_edges + d->b[i]->l / 2.;
	    Get_Dist_Btw_Edges(d,d->v[i],tree);
	  }
    }


}

/*********************************************************/

void Detect_Polytomies(edge *b, phydbl l_thresh, arbre *tree)
{
  if((b->l < l_thresh) && (!b->left->tax) && (!b->rght->tax))
    {
      b->l               = 0.0;
      b->has_zero_br_len = 1;
    }
  else b->has_zero_br_len = 0;
}

/*********************************************************/

void Get_List_Of_Nodes_In_Polytomy(node *a, node *d, node ***list, int *size_list)
{
  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      if(!d->b[i]->has_zero_br_len)
		{
		  (*list)[*size_list] = d->v[i];
		  (*size_list)++;
		}

	      if(d->b[i]->has_zero_br_len)
		Get_List_Of_Nodes_In_Polytomy(d,d->v[i],list,size_list);
	    }
	}
    }

}


/*********************************************************/

void Check_Path(node *a, node *d, node *target, arbre *tree)
{
  printf("path---------\n");
  if(d==target) return;
  else Check_Path(d,d->v[tree->t_dir[d->num][target->num]],target,tree);
}


/*********************************************************/

void Connect_Two_Nodes(node *a, node *d)
{
  a->v[0] = d;
  d->v[0] = a;
}

/*********************************************************/

void Get_List_Of_Adjacent_Targets(node *a, node *d, node ***node_list, edge ***edge_list, int *list_size)
{
  int i;

  For(i,3)
    if(a->v[i] == d)
      {
	(*node_list)[*list_size] = a;
	(*edge_list)[*list_size] = a->b[i];
	(*list_size)++;
      }
  if(d->tax) return;
  else
    For(i,3)
      if(d->v[i] != a) Get_List_Of_Adjacent_Targets(d,d->v[i],node_list,edge_list,list_size);
}

/*********************************************************/

void Sort_List_Of_Adjacent_Targets(edge ***list, int list_size)
{
  edge *buff_edge;
  int i,j;

  buff_edge = NULL;

  For(i,list_size-1)
    {
      for(j=i+1;j<list_size;j++)
	if((*list)[j]->topo_dist_btw_edges < (*list)[i]->topo_dist_btw_edges)
	  {
	    buff_edge = (*list)[j];
	    (*list)[j] = (*list)[i];
	    (*list)[i] = buff_edge;
	  }
    }
}

/*********************************************************/

void Make_Best_Spr(arbre *tree)
{
  tree->best_spr = Make_One_Spr(tree);
  Init_One_Spr(tree->best_spr);
}

/*********************************************************/

void Make_Spr_List(arbre *tree)
{
  int i;

  tree->size_spr_list = 2*tree->n_otu-3;
  tree->spr_list = (spr **)mCalloc(2*tree->n_otu-2,sizeof(spr *));

  For(i,2*tree->n_otu-2)
    {
      tree->spr_list[i] = Make_One_Spr(tree);
      Init_One_Spr(tree->spr_list[i]);
    }
  tree->perform_spr_right_away = 0;
}

/*********************************************************/

void Init_One_Spr(spr *a_spr)
{
  a_spr->lnL             = UNLIKELY;
  a_spr->pars            = 1E+5;
  a_spr->depth_path      = 0;
  a_spr->dist            = 0;
  a_spr->init_target_l   = -1.;
  a_spr->l0              = -1.;
  a_spr->l1              = -1.;
  a_spr->l2              = -1.;
  a_spr->n_link          = NULL;
  a_spr->n_opp_to_link   = NULL;
  a_spr->b_opp_to_link   = NULL;
  a_spr->b_target        = NULL;
  a_spr->b_init_target   = NULL;
}

/*********************************************************/

spr *Make_One_Spr(arbre *tree)
{
  spr *a_spr;
  a_spr       = (spr *)mCalloc(1,sizeof(spr));
  a_spr->path = (node **)mCalloc(tree->n_otu,sizeof(node *));
  return a_spr;
}

/*********************************************************/

node *Common_Nodes_Btw_Two_Edges(edge *a, edge *b)
{
  if(a->left == b->left)      return b->left;
  else if(a->left == b->rght) return b->rght;
  else if(a->rght == b->left) return b->left;
  else if(a->rght == b->rght) return b->rght;

  printf("\n. First edge = %d (%d %d); Second edge = %d (%d %d)\n",
	 a->num,a->left->num,a->rght->num,
	 b->num,b->left->num,b->rght->num);
  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("");

  return NULL;
}

/*********************************************************/

int KH_Test(phydbl *site_lk_M1, phydbl *site_lk_M2, arbre *tree)
{
  phydbl *delta,mean,sd,obs_stat,threshold;
  int i;

  
  delta = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));

  threshold = .0;
  mean = .0;
  obs_stat = .0;
  For(i,tree->n_pattern)
    {
      delta[i] = site_lk_M1[i] - site_lk_M2[i];
      mean += ((int)tree->data->wght[i])*delta[i];
    }

  obs_stat = mean;

  mean /= tree->data->init_len;

  For(i,tree->data->init_len) delta[i] -= mean;

  sd = .0;
  For(i,tree->data->init_len) sd += pow(delta[i],2);
  sd /= (phydbl)(tree->data->init_len-1.);

/*   threshold = tree->dnorm_thresh*sqrt(sd*tree->data->init_len); */


/*   printf("\nObs stat = %f Threshold = %f\n",obs_stat,threshold); */
  Free(delta);

  if(obs_stat > threshold) return 1;
  else                     return 0;
}

/*********************************************************/

void Randomize_Spr_List(arbre *tree)
{
  int i,j;
  spr *buff;

  For(i,tree->size_spr_list)
    {
      j = (int)floor(rand()/(RAND_MAX+1.)*tree->size_spr_list);
      buff              = tree->spr_list[i];
      tree->spr_list[i] = tree->spr_list[j];
      tree->spr_list[j] = buff;
    }
}

/*********************************************************/

void Select_Compatible_Spr_Moves(arbre *tree)
{
  int i,j;


  For(i,tree->size_spr_list-1)
    {
      for(j=i+1;j<tree->size_spr_list;j++)
	{
	  if((tree->spr_list[j]->b_target) && (tree->spr_list[i]->b_target))
	    {
	      if(tree->spr_list[i]->n_link == tree->spr_list[j]->n_link)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->n_link == tree->spr_list[j]->n_opp_to_link)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->n_link == tree->spr_list[j]->b_target->left)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->n_link == tree->spr_list[j]->b_target->rght)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->n_opp_to_link == tree->spr_list[j]->n_link)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->n_opp_to_link == tree->spr_list[j]->n_opp_to_link)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->n_opp_to_link == tree->spr_list[j]->b_target->left)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->n_opp_to_link == tree->spr_list[j]->b_target->rght)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->b_target->left == tree->spr_list[j]->n_link)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->b_target->left == tree->spr_list[j]->n_opp_to_link)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->b_target->left == tree->spr_list[j]->b_target->left)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->b_target->left == tree->spr_list[j]->b_target->rght)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->b_target->rght == tree->spr_list[j]->n_link)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->b_target->rght == tree->spr_list[j]->n_opp_to_link)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->b_target->rght == tree->spr_list[j]->b_target->left)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	      else if(tree->spr_list[i]->b_target->rght == tree->spr_list[j]->b_target->rght)
		{
		  tree->spr_list[j]->b_target = NULL;
		}
	    }
	}
    }
}

/*********************************************************/

int Spr(phydbl init_lnL, arbre *tree)
{
  int spr_moves, best_move, br;

  tree->both_sides = 1;
  spr_moves        = 0;

  Reset_Spr_List(tree);

  For(br,2*tree->n_otu-3)
    {
      /* Subtree rooted by the node on the left */
      tree->n_moves = 0;

      if(!tree->t_edges[br]->left->tax) Test_All_Spr_Targets(tree->t_edges[br],tree->t_edges[br]->left,tree);

      if(tree->perform_spr_right_away)
	{
	  if(tree->n_moves)
	    {
	      best_move = Evaluate_List_Of_Regraft_Pos_Triple(tree->spr_list,
							      (int)ceil(tree->mod->s_opt->p_moves_to_examine*tree->n_moves),
							      tree);

	      if(tree->spr_list[best_move]->lnL > init_lnL) Try_One_Spr_Move_Triple(tree->spr_list[best_move],tree);
	      else
		{
		  tree->both_sides = 1;
		  Lk(tree);
		  Pars(tree);
		}
	    }
	  Reset_Spr_List(tree);
	}

      /* Subtree rooted by the node on the right */
      tree->n_moves = 0;

      if(!tree->t_edges[br]->rght->tax) Test_All_Spr_Targets(tree->t_edges[br],tree->t_edges[br]->rght,tree);

      if(tree->perform_spr_right_away)
	{
	  if(tree->n_moves)
	    {
	      best_move = Evaluate_List_Of_Regraft_Pos_Triple(tree->spr_list,
							      (int)ceil(tree->mod->s_opt->p_moves_to_examine*tree->n_moves),
							      tree);
	      

	      if(tree->spr_list[best_move]->lnL > init_lnL) Try_One_Spr_Move_Triple(tree->spr_list[best_move],tree);
	      else
		{
		  tree->both_sides = 1;
		  Lk(tree);
		  Pars(tree);
		}
	    }
	  Reset_Spr_List(tree);
	}
    }
  return 1;
}

/*********************************************************/

int Test_All_Spr_Targets(edge *b_pulled, node *n_link, arbre *tree)
{
  node *n_opp_to_link,*n_v1,*n_v2,*n_up;
  edge *b_target,*b_residual;
  int i,dir1,dir2;
  phydbl init_len_v1, init_len_v2, init_len_pulled;

  n_up = NULL;
  b_target = b_residual = NULL;
  n_opp_to_link  = (n_link == b_pulled->rght)?(b_pulled->left):(b_pulled->rght);

  init_len_pulled = b_pulled->l;
  dir1 = dir2 = -1;
  For(i,3)
    if(n_link->v[i] != n_opp_to_link)
      {
	if(dir1<0) dir1 = i;
	else       dir2 = i;
      }

  if(n_link->v[dir1]->num < n_link->v[dir2]->num)
    {
      n_v1        = n_link->v[dir1];
      n_v2        = n_link->v[dir2];
      init_len_v1 = n_link->b[dir1]->l;
      init_len_v2 = n_link->b[dir2]->l;
    }
  else
    {
      n_v1        = n_link->v[dir2];
      n_v2        = n_link->v[dir1];
      init_len_v1 = n_link->b[dir2]->l;
      init_len_v2 = n_link->b[dir1]->l;
    }


  Prune_Subtree(n_link,n_opp_to_link,&b_target,&b_residual,tree);

  //
/*   Update_PMat_At_Given_Edge(b_target,tree); */
  //

  tree->depth_curr_path = 1; tree->curr_path[0] = b_target->left;
  Test_One_Spr_Target_Recur(b_target->rght,
			    b_target->left,
			    b_pulled,n_link,b_residual,tree);

  tree->depth_curr_path = 1; tree->curr_path[0] = b_target->rght;
  Test_One_Spr_Target_Recur(b_target->left,
			    b_target->rght,
			    b_pulled,n_link,b_residual,tree);

  Graft_Subtree(b_target,n_link,b_residual,tree);

  if((n_link->v[dir1] != n_v1) || (n_link->v[dir2] != n_v2))
    printf("\n. Warning : -- SWITCH NEEDED -- ! \n");

  n_link->b[dir1]->l = init_len_v1; Update_PMat_At_Given_Edge(n_link->b[dir1],tree);
  n_link->b[dir2]->l = init_len_v2; Update_PMat_At_Given_Edge(n_link->b[dir2],tree);
  b_pulled->l = init_len_pulled;
  Update_PMat_At_Given_Edge(b_pulled,tree);

  //
/*   Update_P_Lk(tree,b_pulled,  n_link); */
/*   Update_P_Lk(tree,b_target,  n_link); */
/*   Update_P_Lk(tree,b_residual,n_link); */
  //

  Update_P_Pars(tree,b_pulled,  n_link);
  Update_P_Pars(tree,b_target,  n_link);
  Update_P_Pars(tree,b_residual,n_link);

  if(!tree->perform_spr_right_away)
    /* if perform_spr_right_away != 0 --> a spr move
     * will be performed anyway. Thus it is not necessary 
     * to update the partial likelihoods below 
     */
    {
      For(i,3)
	if(n_link->v[i] != n_opp_to_link)
	  {
/* 	    Pre_Order_Lk(n_link,n_link->v[i],tree); */
	    Pre_Order_Pars(n_link,n_link->v[i],tree);
	  }
    }
  return 0;
}

/*********************************************************/

void Test_One_Spr_Target_Recur(node *a, node *d, edge *pulled, node *link, edge *residual, arbre *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a)
	  {
	    //
/* 	    Update_P_Lk(tree,d->b[i],d); */
	    //
	    Update_P_Pars(tree,d->b[i],d);
	    tree->curr_path[tree->depth_curr_path] = d->v[i];
	    tree->depth_curr_path++;
	    Test_One_Spr_Target(d->b[i],pulled,link,residual,tree);
	    Test_One_Spr_Target_Recur(d,d->v[i],pulled,link,residual,tree);
	    tree->depth_curr_path--;
	  }
    }
}

/*********************************************************/

phydbl Test_One_Spr_Target(edge *b_target, edge *b_arrow, node *n_link, edge *b_residual, arbre *tree)
{
  phydbl init_target_len, init_arrow_len, init_residual_len;
  int i,dir_v0,dir_v1,dir_v2;
  phydbl l0,l1,l2;
  node *v1, *v2;
  phydbl init_lnL, move_lnL;
  int init_pars,move_pars;


  tree->n_moves++;

  move_lnL  = UNLIKELY;
  init_lnL  = tree->c_lnL;
  init_pars = tree->c_pars;

  Graft_Subtree(b_target,n_link,b_residual,tree);

  init_target_len   = b_target->l;
  init_arrow_len    = b_arrow->l;
  init_residual_len = b_residual->l;

  //
  /*   Triple_Dist(n_link,tree); */
  /*   Update_PMat_At_Given_Edge(b_target,tree); */
  /*   Update_PMat_At_Given_Edge(b_arrow,tree); */
  /*   Update_P_Lk(tree,b_residual,n_link); */
  /*   move_lnL = Lk_At_Given_Edge(b_residual,tree); */
  //

  Update_P_Pars(tree,b_residual,n_link);
  move_pars = Pars_At_Given_Edge(b_residual,tree);

  v1 = (b_residual->left == n_link)?(b_residual->rght):(b_residual->left);
  v2 = (b_target->left   == n_link)?(b_target->rght):(b_target->left);
  dir_v1 = dir_v2 = dir_v0 = -1;
  For(i,3)
    {
      if(n_link->v[i]      == v1) dir_v1 = i;
      else if(n_link->v[i] == v2) dir_v2 = i;
      else                        dir_v0 = i;
    }
  l0 = n_link->b[dir_v0]->l;
  if(n_link->v[dir_v1]->num > n_link->v[dir_v2]->num)
    {
      l1 = n_link->b[dir_v2]->l;
      l2 = n_link->b[dir_v1]->l;
    }
  else
    {
      l1 = n_link->b[dir_v1]->l;
      l2 = n_link->b[dir_v2]->l;
    }

  For(i,tree->depth_curr_path) tree->spr_list[tree->size_spr_list]->path[i] = tree->curr_path[i];
  tree->spr_list[tree->size_spr_list]->depth_path    = tree->depth_curr_path;
  tree->spr_list[tree->size_spr_list]->pars          = tree->c_pars;
  tree->spr_list[tree->size_spr_list]->lnL           = tree->c_lnL;
  tree->spr_list[tree->size_spr_list]->b_target      = b_target;
  tree->spr_list[tree->size_spr_list]->n_link        = n_link;
  tree->spr_list[tree->size_spr_list]->n_opp_to_link = (n_link==b_arrow->left)?(b_arrow->rght):(b_arrow->left);
  tree->spr_list[tree->size_spr_list]->b_opp_to_link = b_arrow;
  tree->spr_list[tree->size_spr_list]->l0            = l0;
  tree->spr_list[tree->size_spr_list]->l1            = l1;
  tree->spr_list[tree->size_spr_list]->l2            = l2;
  tree->spr_list[tree->size_spr_list]->dist          = b_target->topo_dist_btw_edges;

  Include_One_Spr_To_List_Of_Spr(tree->spr_list[tree->size_spr_list],tree);

  b_target->l   = init_target_len;
  b_arrow->l    = init_arrow_len;
  b_residual->l = init_residual_len;

  Prune_Subtree(n_link,
		(n_link==b_arrow->left)?(b_arrow->rght):(b_arrow->left),
		&b_target,
		&b_residual,
		tree);

  //
/*   Update_PMat_At_Given_Edge(b_target,tree); */
  //

  tree->c_lnL   = init_lnL;
  tree->c_pars  = init_pars;

  return .0;
}

/*********************************************************/

void Speed_Spr_Loop(arbre *tree)
{
  phydbl lk_old;

  do
    {
      lk_old = tree->c_lnL;
      Speed_Spr(tree);
      Simu(tree,1000);
      Check_NNI_Five_Branches(tree);
    }
  while(tree->c_lnL > lk_old + tree->mod->s_opt->min_diff_lk_global);

}

/*********************************************************/

void Speed_Spr(arbre *tree)
{
  int step,old_pars;
  phydbl old_lnL;

  if(tree->lock_topo)
    {
      printf("\n. The tree topology is locked.");
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  /* Optimise parameters of the Markov model */
  Optimiz_All_Free_Param(tree,0);

  tree->both_sides = 1;
  Pars(tree);
  Lk(tree);
  Record_Br_Len(tree);

  tree->best_lnL = tree->c_lnL;
  old_lnL        = tree->c_lnL;
  old_pars       = tree->c_pars;
  step           = 0;

  do
    {
      ++step;

      if(tree->mod->s_opt->print) printf("\n\n. Starting a SPR cycle... \n");

      Init_One_Spr(tree->best_spr);
      old_lnL  = tree->c_lnL;
      old_pars = tree->c_pars;

      tree->n_improvements         = 0;
      tree->perform_spr_right_away = 1;
      Spr(UNLIKELY,tree);

      /* Optimise parameters of the Markov model */
      Optimiz_All_Free_Param(tree,tree->mod->s_opt->print);

      /* Optimise branch lengths */
      Optimize_Br_Len_Serie(tree->noeud[0],
			    tree->noeud[0]->v[0],
			    tree->noeud[0]->b[0],
			    tree,
			    tree->data);

      /* Update partial likelihoods & parsimony */
      tree->both_sides = 1;
      Lk(tree);
      Pars(tree);

      /* Print log-likelihood and parsimony scores */
      if(tree->mod->s_opt->print) Print_Lk(tree,"[Branch lengths     ]");

      /* Record the current best log-likleihood  */
      tree->best_lnL = tree->c_lnL;

      if(tree->c_lnL < old_lnL)
	{
	  printf("\n. old_lnL = %f c_lnL = %f",old_lnL,tree->c_lnL); 
	  printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      /* Record the current best branch lengths  */
      Record_Br_Len(tree);

      /* Exit if no improvements after complete optimization */      
      if((!tree->n_improvements) || 
	 (fabs(old_lnL-tree->c_lnL) < tree->mod->s_opt->min_diff_lk_global)) break;

    }while(1);
}

/*********************************************************/

void Best_Spr(arbre *tree)
{
  int best_move,n_moves;
  phydbl best_move_lnL,init_lnL;
  int i;

  if(tree->mod->s_opt->print) printf("\n\n. Starting a SPR cycle... \n");

  /* Optimise parameters of the Markov model */
  Optimiz_All_Free_Param(tree,0);

  tree->both_sides = 1;
  Pars(tree);
  Lk(tree);
  Record_Br_Len(tree);
  
  init_lnL                     = tree->c_lnL;
  tree->best_lnL               = tree->c_lnL;
  best_move                    = -1;
  best_move_lnL                = UNLIKELY;
  tree->perform_spr_right_away = 0;

  Spr(UNLIKELY,tree);
  tree->both_sides = 1;
  Pars(tree);
  Lk(tree);

  n_moves = MAX(tree->size_spr_list,
		(int)ceil(2*tree->mod->s_opt->p_moves_to_examine*tree->size_spr_list));
  For(i,n_moves)
    {
      Evaluate_One_Regraft_Pos_Triple(tree->spr_list[i],tree);
      tree->both_sides = 1;
      Lk(tree);
      if(tree->spr_list[i]->lnL > best_move_lnL)
	{
	  best_move_lnL = tree->spr_list[i]->lnL;
	  best_move     = i;
	}
    }

  if(best_move > -1) Try_One_Spr_Move_Triple(tree->spr_list[best_move],tree);
  
  /* Optimise parameters of the Markov model */
  Optimiz_All_Free_Param(tree,tree->mod->s_opt->print);
  
  /* Optimise branch lengths */
  Optimize_Br_Len_Serie(tree->noeud[0],
			tree->noeud[0]->v[0],
			tree->noeud[0]->b[0],
			tree,
			tree->data);
  
  /* Update partial likelihoods & parsimony */
  tree->both_sides = 1;
  Lk(tree);
  Pars(tree);
  
  /* Print log-likelihood and parsimony scores */
  if(tree->mod->s_opt->print) Print_Lk(tree,"[Topology           ]");
  
  /* Record the current best log-likleihood  */
  tree->best_lnL = tree->c_lnL;
  
  if(tree->c_lnL < init_lnL)
    {
      printf("\n. init_lnL = %f c_lnL = %f\n",init_lnL,tree->c_lnL); 
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  /* Record the current best branch lengths  */
  Record_Br_Len(tree);  
}

/*********************************************************/

int Evaluate_List_Of_Regraft_Pos_Triple(spr **spr_list, int list_size, arbre *tree)
{
  spr *move;
  edge *init_target, *b_residual;
  int i,j,best_move;
  int dir_v0, dir_v1, dir_v2;
  phydbl recorded_l;
  phydbl move_lnL, best_lnL;

  move_lnL = best_lnL = UNLIKELY;
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


  For(i,2*tree->n_otu-3)
    {
      tree->t_edges[i]->is_p_lk_r_u2d = 0;
      tree->t_edges[i]->is_p_lk_l_u2d = 0;
    }
  
  recorded_l = -1.0;
  For(i,list_size)
    {
      move = spr_list[i];

      if(move->b_target)
	{
	  /* Record edge lengths */
	  Record_Br_Len(tree);

	  /* Prune subtree */
	  Prune_Subtree(move->n_link,move->n_opp_to_link,&init_target,&b_residual,tree);

	  if(recorded_l < 0.0)
	    {
	      /* Rough optimisation of the branch length at prune site
	       * We only need to perform this optimisation for the first
	       * element of spr_list because the pruned subtree is the
	       * same across all the elements of spr_list. It would not
	       * be true in the general case
	       */
	      Fast_Br_Len(init_target,tree);	      
	      /* Record branch length at prune site */
	      move->init_target_l = init_target->l;
	      recorded_l          = init_target->l;
	    }
	  else
	    {
	      init_target->l      = recorded_l;
	      move->init_target_l = recorded_l;
	    }

	  /* Update the change proba matrix at prune position */
	  Update_PMat_At_Given_Edge(init_target,tree);
	  
	  /* Update conditional likelihoods along the path from the prune to
	     the regraft position */
	  Update_P_Lk_Along_A_Path(move->path,move->depth_path,tree);
	  
	  /* Regraft subtree */
	  Graft_Subtree(move->b_target,move->n_link,b_residual,tree);
	  
	  Update_PMat_At_Given_Edge(move->b_target,tree);
	  Update_PMat_At_Given_Edge(b_residual,tree);
	  Update_P_Lk(tree,move->b_opp_to_link,move->n_link);
	  move_lnL = Lk_At_Given_Edge(move->b_opp_to_link,tree);
	 
	  if(move_lnL > best_lnL)
	    {
	      best_lnL = move_lnL;
	      best_move = i;
	    }
	  else
	    {
	      /* Estimate the three edge lengths at the regraft site */
	      Triple_Dist(move->n_link,tree);

	      move_lnL = Lk_At_Given_Edge(move->b_opp_to_link,tree);
	      
	      if(move_lnL > best_lnL)
		{
		  best_lnL = move_lnL;
		  best_move = i;
		}
	    }

	  if(move_lnL > tree->best_spr->lnL)
	    {
	      tree->best_spr->lnL           = move_lnL;
	      tree->best_spr->pars          = tree->c_pars;
	      tree->best_spr->b_target      = move->b_target;
	      tree->best_spr->n_link        = move->n_link;
	      tree->best_spr->n_opp_to_link = move->n_opp_to_link;
	      tree->best_spr->b_opp_to_link = move->b_opp_to_link;
	    }

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

	  /* Record likelihood */
	  move->lnL = tree->c_lnL;

	  /* Regraft the subtree at its original position */
	  Prune_Subtree(move->n_link,
			move->n_opp_to_link,
			&move->b_target,
			&b_residual,
			tree);

	  Graft_Subtree(init_target,
			move->n_link,
			b_residual,
			tree);

	  /* Restore branch lengths */
	  Restore_Br_Len(tree);

	  /* Update relevant change proba matrices */
	  Update_PMat_At_Given_Edge(move->b_target,tree);
	}
    }

  For(i,list_size)
    {
      move = spr_list[i];      
      if(move->b_target)
	{
	  For(j,3) Update_PMat_At_Given_Edge(move->n_link->b[j],tree);	  
	  For(j,3) Update_P_Lk(tree,move->n_link->b[j],move->n_link);
	  break;
	}
    }

#ifdef DEBUG
  if(best_move < 0)
    {
      printf("\n\n. Best_move < 0 !");

      printf("\n. List size = %d",list_size);
      For(i,list_size)
	{
	  move = spr_list[i];
	  printf("\n. %p %p",move,move->b_target);
	}

      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif

  return best_move;
}

/*********************************************************/

void Evaluate_One_Regraft_Pos_Triple(spr *move, arbre *tree)
{
  int j;
  edge *init_target, *b_residual;
  phydbl move_lnL;
  int dir_v0, dir_v1, dir_v2;

  move_lnL = UNLIKELY;
  init_target = b_residual = NULL;

  For(j,2*tree->n_otu-3)
    {
      tree->t_edges[j]->is_p_lk_r_u2d = 0;
      tree->t_edges[j]->is_p_lk_l_u2d = 0;
    }

  if(move->b_target)
    {
      /* Record edge lengths */
      Record_Br_Len(tree);
      
      /* Prune subtree */
      Prune_Subtree(move->n_link,move->n_opp_to_link,&init_target,&b_residual,tree);
      
      /* Rough optimisation of the branch length at prune site */
      Fast_Br_Len(init_target,tree);
      
      /* Record edge length at prune site */
      move->init_target_l = init_target->l;

      /* Update the change proba matrix at prune position */
      Update_PMat_At_Given_Edge(init_target,tree);
      
      /* Update partial likelihood along the path from the prune to
	 the regraft position */
      Update_P_Lk_Along_A_Path(move->path,move->depth_path,tree);
      
      /* Regraft subtree */
      Graft_Subtree(move->b_target,move->n_link,b_residual,tree);
      
      /* Estimate the three edge lengths at the regraft site */
      Triple_Dist(move->n_link,tree);

      move_lnL = Lk_At_Given_Edge(move->b_opp_to_link,tree);
            
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
      
      /* Record likelihood */
      move->lnL = tree->c_lnL;
      
      /* Regraft the subtree at its original position */
      Prune_Subtree(move->n_link,
		    move->n_opp_to_link,
		    &move->b_target,
		    &b_residual,
		    tree);
      
      Graft_Subtree(init_target,
		    move->n_link,
		    b_residual,
		    tree);
      
      /* Restore branch lengths */
      Restore_Br_Len(tree);
      
      /* Update relevant change proba matrices */
      Update_PMat_At_Given_Edge(move->b_target,tree);
      For(j,3) Update_PMat_At_Given_Edge(move->n_link->b[j],tree);
      
      /* Update relevant partial likelihoods */
      For(j,3) Update_P_Lk(tree,move->n_link->b[j],move->n_link);
    }
}

/*********************************************************/

int Try_One_Spr_Move_Triple(spr *move, arbre *tree)
{
  edge *init_target, *b_residual;
  int j;
  int dir_v0, dir_v1, dir_v2;


  Record_Br_Len(tree);

  Prune_Subtree(move->n_link,
		move->n_opp_to_link,
		&init_target,
		&b_residual,
		tree);

  init_target->l = move->init_target_l;

  Graft_Subtree(move->b_target,move->n_link,b_residual,tree);

  dir_v1 = dir_v2 = dir_v0 = -1;
  For(j,3)
    {
      if(move->n_link->v[j] == move->n_opp_to_link) dir_v0 = j;
      else if(dir_v1 < 0)                           dir_v1 = j;
      else                                          dir_v2 = j;
    }

  move->n_link->b[dir_v0]->l = move->l0;

  if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num)
    {
      move->n_link->b[dir_v2]->l = move->l1;
      move->n_link->b[dir_v1]->l = move->l2;
    }
  else
    {
      move->n_link->b[dir_v1]->l = move->l1;
      move->n_link->b[dir_v2]->l = move->l2;
    }

  if(move->lnL > tree->best_lnL)
    {
      time(&(tree->t_current));
      tree->both_sides = 1;
      Lk(tree);

      if(fabs(tree->c_lnL - move->lnL) > tree->mod->s_opt->min_diff_lk_local)
	{
	  if(tree->mod->s_opt->print) printf("\n. c_lnL = %f move_lnL = %f",
					     tree->c_lnL,move->lnL);
	  printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      Pars(tree);
      if(tree->mod->s_opt->print) Print_Lk(tree,"[Topology           ]");


      tree->n_improvements++;
      tree->best_lnL = tree->c_lnL;
      Record_Br_Len(tree);
      return 1;
    }
/*   else */
/*     { */
/*       Fast_Br_Len_Recur(tree->noeud[0], */
/* 			tree->noeud[0]->v[0], */
/* 			tree->noeud[0]->b[0], */
/* 			tree); */

/*       tree->both_sides = 1; */
/*       Lk(tree); */

/*       if(tree->c_lnL > tree->best_lnL) */
/* 	{ */
/* 	  Pars(tree); */
/* 	  if(tree->mod->s_opt->print) Print_Lk_And_Pars(tree); */
/* 	  tree->n_improvements++; */
/* 	  tree->best_lnL = tree->c_lnL; */
/* 	  Record_Br_Len(tree); */
/* 	  return 1; */
/* 	} */
/*     } */

  Prune_Subtree(move->n_link,
		move->n_opp_to_link,
		&move->b_target,
		&b_residual,
		tree);

  Graft_Subtree(init_target,
		move->n_link,
		b_residual,
		tree);

  Restore_Br_Len(tree);
  tree->both_sides = 1;
  Lk(tree);
  Pars(tree);
  return 0;
}

/*********************************************************/

int Try_One_Spr_Move_Full(spr *move, arbre *tree)
{
  edge *init_target, *b_residual;

  Record_Br_Len(tree);

  Prune_Subtree(move->n_link,
		move->n_opp_to_link,
		&init_target,
		&b_residual,
		tree);

  Graft_Subtree(move->b_target,move->n_link,b_residual,tree);

  tree->both_sides = 1;
  Lk(tree);

  Optimize_Br_Len_Serie(tree->noeud[0],
			tree->noeud[0]->v[0],
			tree->noeud[0]->b[0],
			tree,
			tree->data);
/*   Fast_Br_Len_Recur(tree->noeud[0], */
/* 		    tree->noeud[0]->v[0], */
/* 		    tree->noeud[0]->b[0], */
/* 		    tree); */
  
  tree->both_sides = 1;
  Lk(tree);

  if(tree->c_lnL > tree->best_lnL)
    {
      Pars(tree);
      if(tree->mod->s_opt->print) Print_Lk(tree,"[Topology           ]");
      tree->n_improvements++;
      tree->best_lnL = tree->c_lnL;
      Record_Br_Len(tree);
      return 1;
    }
  else
    {
      Prune_Subtree(move->n_link,
		    move->n_opp_to_link,
		    &move->b_target,
		    &b_residual,
		    tree);
      
      Graft_Subtree(init_target,
		    move->n_link,
		    b_residual,
		    tree);
      
      Restore_Br_Len(tree);
      tree->both_sides = 1;
      Lk(tree);
      Pars(tree);
      return 0;
    }

  return -1;
}

/*********************************************************/

void Include_One_Spr_To_List_Of_Spr(spr *move, arbre *tree)
{
  int i;
  spr *buff_spr;

  if(move->pars <= tree->spr_list[tree->size_spr_list-1]->pars)
/*   if(move->lnL > tree->spr_list[tree->size_spr_list-1]->lnL) */
    {
      tree->spr_list[tree->size_spr_list-1]->depth_path    = move->depth_path;
      tree->spr_list[tree->size_spr_list-1]->pars          = move->pars;
      tree->spr_list[tree->size_spr_list-1]->lnL           = move->lnL;
      tree->spr_list[tree->size_spr_list-1]->b_target      = move->b_target;
      tree->spr_list[tree->size_spr_list-1]->n_link        = move->n_link;
      tree->spr_list[tree->size_spr_list-1]->n_opp_to_link = move->n_opp_to_link;
      tree->spr_list[tree->size_spr_list-1]->b_opp_to_link = move->b_opp_to_link;
      tree->spr_list[tree->size_spr_list-1]->l0            = move->l0;
      tree->spr_list[tree->size_spr_list-1]->l1            = move->l1;
      tree->spr_list[tree->size_spr_list-1]->l2            = move->l2;
      tree->spr_list[tree->size_spr_list-1]->dist          = move->dist;

      For(i,tree->spr_list[tree->size_spr_list-1]->depth_path)
	tree->spr_list[tree->size_spr_list-1]->path[i] = move->path[i];

      for(i=tree->size_spr_list-1;i>0;i--)
	{
	  if(tree->spr_list[i]->pars <= tree->spr_list[i-1]->pars)
/* 	  if(tree->spr_list[i]->lnL > tree->spr_list[i-1]->lnL) */
	    {
	      buff_spr            = tree->spr_list[i-1];
	      tree->spr_list[i-1] = tree->spr_list[i];
	      tree->spr_list[i]   = buff_spr;
	    }
	  else  break;
	}
    }
}

/*********************************************************/

void Random_Tree(arbre *tree)
{
  int *is_available,*list_of_nodes;
  int i,node_num,step,n_available;

  is_available  = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  list_of_nodes = (int *)mCalloc(tree->n_otu,    sizeof(int));

  For(i,tree->n_otu) is_available[i]  = 1;
  For(i,tree->n_otu) list_of_nodes[i] = i;

  step = 0;
  do
    {
      node_num = (int)rint(rand()/(phydbl)(RAND_MAX+1.0)*(tree->n_otu-1-step));
      node_num = list_of_nodes[node_num];
      is_available[node_num] = 0;
      For(i,tree->n_otu) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}

      tree->noeud[node_num]->v[0] = tree->noeud[tree->n_otu+step];
      tree->noeud[tree->n_otu+step]->v[1] = tree->noeud[node_num];

      node_num = (int)rint(rand()/(phydbl)(RAND_MAX+1.0)*(tree->n_otu-2-step));
      node_num = list_of_nodes[node_num];
      is_available[node_num] = 0;
      For(i,tree->n_otu) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}

      tree->noeud[node_num]->v[0] = tree->noeud[tree->n_otu+step];
      tree->noeud[tree->n_otu+step]->v[2] = tree->noeud[node_num];

      is_available[tree->n_otu+step] = 1;
      For(i,tree->n_otu) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) list_of_nodes[n_available++] = i;

      step++;
    }while(step < tree->n_otu-2);

  tree->noeud[list_of_nodes[0]]->v[0] = tree->noeud[list_of_nodes[1]];
  tree->noeud[list_of_nodes[1]]->v[0] = tree->noeud[list_of_nodes[0]];

  tree->num_curr_branch_available = 0;
  Connect_Edges_To_Nodes_Recur(tree->noeud[0],tree->noeud[0]->v[0],tree);
/*   Print_Node(tree->noeud[0],tree->noeud[0]->v[0],tree); */
  
  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  
  Free(is_available);
  Free(list_of_nodes);
}

/*********************************************************/

void Random_Spr(int n_moves, arbre *tree)
{
  int i;
  int br_pulled, br_target;
  spr *spr_struct;
  edge *target, *residual;

  spr_struct = Make_One_Spr(tree);
  Init_One_Spr(spr_struct);
  target = residual = NULL;
  For(i,n_moves)
    {
      br_pulled = (int)((phydbl)rand()/RAND_MAX * (2*tree->n_otu-3-1));
      do
	{
	  br_target = (int)((phydbl)rand()/RAND_MAX * (2*tree->n_otu-3-1));
	}while(br_target == br_pulled);

      spr_struct->n_link        = tree->t_edges[br_pulled]->left;
      spr_struct->n_opp_to_link = tree->t_edges[br_pulled]->rght;
      spr_struct->b_opp_to_link = tree->t_edges[br_pulled];
      spr_struct->b_target      = tree->t_edges[br_target];
      spr_struct->b_init_target = NULL;

      if(!Check_Spr_Move_Validity(spr_struct,tree))
	{
	  spr_struct->n_link        = tree->t_edges[br_pulled]->rght;
	  spr_struct->n_opp_to_link = tree->t_edges[br_pulled]->left;
	}

#ifdef DEBUG
      if(!Check_Spr_Move_Validity(spr_struct,tree))
	{
	  Warn_And_Exit("\n. Could not find a valid move...\n");
	}
#endif

      Prune_Subtree(spr_struct->n_link,
		    spr_struct->n_opp_to_link,
		    &target,
		    &residual,
		    tree);

      Graft_Subtree(spr_struct->b_target,
		    spr_struct->n_link,
		    residual,tree);
    }
  Free(spr_struct);
}

/*********************************************************/

void Random_NNI(int n_moves, arbre *tree)
{
  int i,j;
  edge *b;
  node *n1,*n2,*n_target;

  n1 = n2 = NULL;
  b = NULL;
  For(i,n_moves)
    {
      n_target  = tree->noeud[tree->n_otu + (int)((phydbl)rand()/RAND_MAX * (2*tree->n_otu-3-tree->n_otu))];
      For(j,3) if(!n_target->v[j]->tax) {b = n_target->b[j]; break;}


      For(j,3) if(b->left->v[j] != b->rght) {n1 = b->left->v[j]; break;}
      For(j,3) if(b->rght->v[j] != b->left) {n2 = b->rght->v[j]; break;}

      Swap(n1,b->left,b->rght,n2,tree);
    }
}

/*********************************************************/

void Reset_Spr_List(arbre *tree)
{
  int i;

  tree->spr_list[0]->depth_path     = 0;
  tree->spr_list[0]->pars           = MAX_PARS;
  tree->spr_list[0]->lnL            = UNLIKELY;
  tree->spr_list[0]->n_link         = NULL;
  tree->spr_list[0]->n_opp_to_link  = NULL;
  tree->spr_list[0]->b_target       = NULL;

  for(i=1;i<tree->size_spr_list;i++)
    {
      tree->spr_list[i]->depth_path     = 0;
      tree->spr_list[i]->pars           = tree->spr_list[0]->pars;
      tree->spr_list[i]->lnL            = tree->spr_list[0]->lnL;
      tree->spr_list[i]->n_link         = NULL;
      tree->spr_list[i]->n_opp_to_link  = NULL;
      tree->spr_list[i]->b_target       = NULL;
    }
}

/*********************************************************/

void Print_Settings(option *io)
{
  int answer;

  printf("\n\n\n");
  printf("\n\n");

  printf("                                 ..........................                                      \n");
  printf(" ooooooooooooooooooooooooooooo        CURRENT SETTINGS        ooooooooooooooooooooooooooooooooooo\n");
  printf("                                 ..........................                                      \n");

  printf("\n                . Sequence filename : \t\t\t\t %s", io->in_seq_file);
  printf("\n                . Data type :             \t\t\t %s", (io->mod->datatype ? "aa" : "dna"));
  printf("\n                . Sequence format : \t\t\t\t %s", io->interleaved ? "interleaved" : "sequential");
  printf("\n                . Number of data sets : \t\t\t %d", io->n_data_sets);

  printf("\n                . Nb of bootstrapped data sets : \t\t %d", io->mod->bootstrap);

  if (io->mod->bootstrap > 0)
    printf("\n                . Compute approximate likelihood ratio test : \t no");
  else
    {
      if(io->ratio_test == 1)
	printf("\n                . Compute approximate likelihood ratio test : \t yes (aLRT statistics)");
      else if(io->ratio_test == 2)
	printf("\n                . Compute approximate likelihood ratio test : \t yes (Chi2-based parametric branch supports)");
      else if(io->ratio_test == 3)
	printf("\n                . Compute approximate likelihood ratio test : \t yes (Minimum of SH-like and Chi2-based branch supports)");
      else if(io->ratio_test == 4)
	printf("\n                . Compute approximate likelihood ratio test : \t yes (SH-like branch supports)");
    }

  printf("\n                . Model name : \t\t\t\t\t %s", io->mod->modelname);

  if (io->mod->datatype == NT)
    {
      if ((io->mod->whichmodel == K80)  ||
	  (io->mod->whichmodel == HKY85)||
	  (io->mod->whichmodel == F84)  ||
	  (io->mod->whichmodel == TN93))
	{
	  if (io->mod->s_opt->opt_kappa)
	    printf("\n                . Ts/tv ratio : \t\t\t\t estimated");
	  else
	    printf("\n                . Ts/tv ratio : \t\t\t\t %f", io->mod->kappa);
	}
    }

  if (io->mod->s_opt->opt_pinvar)
    printf("\n                . Proportion of invariable sites :\t\t estimated");
  else
    printf("\n                . Proportion of invariable sites :\t\t %f", io->mod->pinvar);


  printf("\n                . Number of subst. rate categs : \t\t %d", io->mod->n_catg);
  if(io->mod->s_opt->opt_alpha)
    printf("\n                . Gamma distribution parameter : \t\t estimated");
  else
    printf("\n                . Gamma distribution parameter : \t\t %f", io->mod->alpha);
  
  if(io->mod->datatype == AA)
    printf("\n                . Amino acid equilibrium frequencies : \t\t %s", (io->mod->s_opt->opt_state_freq) ? ("empirical"):("model"));
  else if(io->mod->datatype == NT)
    if((io->mod->whichmodel != JC69) &&
       (io->mod->whichmodel != K80)  &&
       (io->mod->whichmodel != F81))
      printf("\n                . Nucleotide equilibrium frequencies : \t\t %s", (io->mod->s_opt->opt_state_freq) ? ("ML"):("empirical"));


  printf("\n                . Optimise tree topology : \t\t\t %s", (io->mod->s_opt->opt_topo) ? "yes" : "no");
  if(io->mod->s_opt->opt_topo)
    {
      printf("\n                . Tree topology search : \t\t\t %s%s", 
	     (io->mod->s_opt->topo_search == NNI_MOVE)?("NNIs"):("SPRs"),
	     (io->mod->s_opt->spr_step_after_nnis)?("+SPRs"):(""));


      printf("\n                . Random input tree : \t\t\t\t %s", (io->mod->s_opt->random_input_tree) ? "yes" : "no");
      if(!io->mod->s_opt->random_input_tree)
	printf("\n                . Starting tree : \t\t\t\t %s", (!io->in_tree) ? "BIONJ" : io->in_tree_file);
      else
	printf("\n                . Number of random starting trees : \t\t %d", io->mod->s_opt->n_rand_starts);	
    }
  printf("\n                . Optimise branch lengths : \t\t\t %s", (io->mod->s_opt->opt_bl) ? "yes" : "no");

  answer = 0;
  if(io->mod->s_opt->opt_alpha  ||
     io->mod->s_opt->opt_kappa  ||
     io->mod->s_opt->opt_lambda ||
     io->mod->s_opt->opt_pinvar ||
     io->mod->s_opt->opt_rr) answer = 1;
  
  printf("\n                . Optimise substitution model parameters : \t %s", (answer) ? "yes" : "no");



  printf("\n\n oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");

  printf("\n\n");
  fflush(NULL);
}


/*********************************************************/

void Fill_Missing_Dist(matrix *mat)
{
  int i,j;
  For(i,mat->n_otu)
    {
      for(j=i+1;j<mat->n_otu;j++)
	{
	  if(i != j)
	    {
	      if(mat->dist[i][j] < .0) 
		{
		  Fill_Missing_Dist_XY(i,j,mat);
		  mat->dist[j][i] = mat->dist[i][j];
		}
	    }
	}
    }
}

/*********************************************************/

void Fill_Missing_Dist_XY(int x, int y, matrix *mat)
{

  int i,j;
  phydbl *local_mins,**S1S2;
  int cpt;
  int pos_best_estimate;
  phydbl min_crit, curr_crit;

  local_mins = (phydbl *)mCalloc(mat->n_otu*mat->n_otu,sizeof(phydbl ));
  S1S2       = (phydbl **)mCalloc(mat->n_otu*mat->n_otu,sizeof(phydbl *));
  For(i,mat->n_otu*mat->n_otu) S1S2[i] = (phydbl *)mCalloc(2,sizeof(phydbl));

  cpt = 0;
  For(i,mat->n_otu)
    {
      if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
	{
	  For(j,mat->n_otu)
	    {
	      if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
		{
		  if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
		    {
		      S1S2[cpt][0] = MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
		      S1S2[cpt][1] = MAX(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
		      cpt++;
		    }
		}
	    }
	}
    }

  Qksort_Matrix(S1S2,0,0,cpt-1);

  local_mins[0] = S1S2[0][1];
  for(i=1;i<cpt;i++) local_mins[i] = (i*local_mins[i-1] + S1S2[i][1])/(phydbl)(i+1);
 
  pos_best_estimate = 0;
  min_crit = curr_crit = MDBL_MAX;
	
  For(i,cpt-1)
    {
      if((local_mins[i] < S1S2[i+1][0]) && (local_mins[i] > S1S2[i][0]))
	{
	  curr_crit = Least_Square_Missing_Dist_XY(x,y,local_mins[i],mat);
	  if(curr_crit < min_crit)
	    {
	      min_crit = curr_crit;
	      pos_best_estimate = i;
	    }
	}
    }

  mat->dist[x][y] = local_mins[pos_best_estimate];
  mat->dist[y][x] = mat->dist[x][y];

  For(i,mat->n_otu*mat->n_otu) Free(S1S2[i]);
  Free(S1S2);
  Free(local_mins);
}

/*********************************************************/

phydbl Least_Square_Missing_Dist_XY(int x, int y, phydbl dxy, matrix *mat)
{
  int i,j;
  phydbl fit;

  fit = .0;
  For(i,mat->n_otu)
    {
      if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
	{
	  For(j,mat->n_otu)
	    {
	      if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
		{
		  if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
		    {
		      if(dxy < MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]))
			{
			  fit += pow((mat->dist[i][x] + mat->dist[j][y]) - (mat->dist[i][y] + mat->dist[j][x]),2);
			}
		      else if((mat->dist[i][x] + mat->dist[j][y]) < (mat->dist[i][y] + mat->dist[j][x]))
			{
			  fit += pow(dxy - (mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]),2);
			}
		      else
			{
			  fit += pow(dxy - (mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j]),2);
			}
		    }
		}
	    }
	}
    }
  return fit;
}

/*********************************************************/

void Print_Banner(FILE *fp)
{
  fprintf(fp,"\n");
  fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
  fprintf(fp,"                                                                                                  \n");
  fprintf(fp,"                                  ---  PhyML %s  ---                                             \n",VVERSION);
  fprintf(fp,"                                                                                                  \n");
  fprintf(fp,"    A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood    \n");
  fprintf(fp,"                            Stephane Guindon & Olivier Gascuel                                      \n");
  fprintf(fp,"                                                                                                  \n");
  fprintf(fp,"                                http://atgc.lirmm.fr/phyml                                          \n");
  fprintf(fp,"                                                                                                  \n");
  fprintf(fp,"                         Copyright CNRS - Universite Montpellier II                                 \n");
  fprintf(fp,"                                                                                                  \n");
  fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

/*********************************************************/


void Print_Data_Set_Number(option *io, FILE *fp)
{
  fprintf(fp,"\n");
  fprintf(fp,"                                                                                                  \n");
  fprintf(fp,"                                 [ Data set number %3d ]                                           \n",io->curr_gt+1);
  fprintf(fp,"                                                                                                  \n");
}
/*********************************************************/


void Check_Memory_Amount(arbre *tree)
{
  /* Rough estimate of the amount of memory that has to be used */

  int nbytes;
  model *mod;

  mod = tree->mod;

  nbytes = 0;


  /* Pmat */
  nbytes += (2*mod->n_otu-3) * mod->n_catg * mod->ns * mod->ns * sizeof(phydbl);
  nbytes += (2*mod->n_otu-3) * mod->n_catg * mod->ns * sizeof(phydbl *);
  nbytes += (2*mod->n_otu-3) * mod->n_catg * sizeof(phydbl **);
  nbytes += (2*mod->n_otu-3) * mod->n_catg * sizeof(phydbl ***);
  
  /* Partial Lk */
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * mod->n_catg * mod->ns * sizeof(phydbl);
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * mod->n_catg * sizeof(phydbl *);
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * sizeof(phydbl **);
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * sizeof(phydbl ***);
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * sizeof(phydbl);

  /* Partial Pars */
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * mod->ns * sizeof(int);
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * sizeof(int *);
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * sizeof(int );

  if(((phydbl)nbytes/(1.E+06)) > 256.)
    {
      char answer;
      printf("\n. WARNING: this analysis requires at least %.0fMo of memory space.\n",(phydbl)nbytes/(1.E+06));
      printf("\n. Do you really want to continue ? [Y/n] ");
      scanf("%c", &answer);
      if(answer == '\n') answer = 'Y';
      else if(answer == 'n' || answer == 'N') Warn_And_Exit("\n\n");
      else getchar();
    }
  else if(((phydbl)nbytes/(1.E+06)) > 100.)
    printf("\n. WARNING: this analysis will use at least %.0fMo of memory space...",(phydbl)nbytes/(1.E+06));
  else
    printf("\n. This analysis requires at least %.0fMo of memory space.",(phydbl)nbytes/(1.E+06));
  
}

/*********************************************************/

int Get_State_From_P_Lk(phydbl *p_lk, arbre *tree)
{
  int i;
  For(i,tree->mod->ns) if(p_lk[i] > .0) return i;
  return -1;
}

/*********************************************************/

int Get_State_From_P_Pars(short int *p_pars, arbre *tree)
{
  int i;
  For(i,tree->mod->ns) if(p_pars[i] > .0) return i;
  return -1;
}

/*********************************************************/

void Print_Lk(arbre *tree, char *string)
{
  time(&(tree->t_current));
  printf("\n. (%5d sec)  [%10.6f] %s",(int)(tree->t_current-tree->t_beg),tree->c_lnL,string);
#ifndef QUIET 
  fflush(NULL);
#endif
}

/*********************************************************/

void Print_Pars(arbre *tree)
{
  time(&(tree->t_current));
  printf("\n. (%5d sec)  [%5d]",(int)(tree->t_current-tree->t_beg),tree->c_pars);
#ifndef QUIET
  fflush(NULL);
#endif
}

/*********************************************************/

void Print_Lk_And_Pars(arbre *tree)
{	
  time(&(tree->t_current));

  printf("\n. (%5d sec) [%10.4f] [%5d]",
	 (int)(tree->t_current-tree->t_beg),
	 tree->c_lnL,tree->c_pars);
#ifndef QUIET
  fflush(NULL);
#endif
}

/*********************************************************/

void Check_Dirs(arbre *tree)
{
  int i;

  For(i,2*tree->n_otu-3)
    {
      if(!tree->t_edges[i]->left->tax)
	{
	  if(tree->t_edges[i]->left->v[tree->t_edges[i]->l_v1]->num <
	     tree->t_edges[i]->left->v[tree->t_edges[i]->l_v2]->num)
	    {
	      printf("\n. Edge %d ; v1=%d v2=%d",
		     tree->t_edges[i]->num,
		     tree->t_edges[i]->left->v[tree->t_edges[i]->l_v1]->num,
		     tree->t_edges[i]->left->v[tree->t_edges[i]->l_v2]->num);
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	}

      if(!tree->t_edges[i]->rght->tax)
	{
	  if(tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v1]->num <
	     tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v2]->num)
	    {
	      printf("\n. Edge %d ; v3=%d v4=%d",
		     tree->t_edges[i]->num,
		     tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v1]->num,
		     tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v2]->num);
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	}
    }
}

/*********************************************************/

void Warn_And_Exit(char *s)
{
  fprintf(stdout,"%s",s);
  fflush(NULL);
  #ifndef BATCH
  char c;
  fprintf(stdout,"\n. Type any key to exit.\n");
  fscanf(stdin,"%c",&c);
  #endif
  Exit("\n");
}

/*********************************************************/

void Read_Qmat(double *daa, phydbl *pi, FILE *fp)
{
  int i,j;
  phydbl sum;

  for(i=1;i<20;i++)
    {
      For(j,19)
	{
	  fscanf(fp,"%lf",&(daa[i*20+j]));
	  daa[j*20+i] = daa[i*20+j];
	  if(j == i-1) break; 
	}
    }

  For(i,20) fscanf(fp,"%lf",pi+i);
  sum = .0;
  For(i,20) sum += pi[i];
  if(fabs(sum - 1.) > 1.E-06)
    {
      printf("\n. Error : the rate matrix format is incorrect\n");
      Warn_And_Exit("\n");
    }
}

/*********************************************************/

void Print_Qmat_AA(double *daa, phydbl *pi)
{
  int i,j,cpt;

  cpt = 0;
  For(i,20)
    {
      for(j=0;j<i;j++)
	{
	  printf("daa[%2d*20+%2d] = %10f;  ",i,j,daa[i*20+j]);
	  cpt++;
	  if(!(cpt%4)) printf("\n");
	}
    }

  printf("\n\n");
  printf("for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];\n\n");
  For(i,20) printf("pi[%d] = %f; ",i,pi[i]);
  printf("\n");
  printf("Ala\tArg\tAsn\tAsp\tCys\tGln\tGlu\tGly\tHis\tIle\tLeu\tLys\tMet\tPhe\tPro\tSer\tThr\tTrp\tTyr\tVal\n");
}


/*********************************************************/

void Randomize_Sequence_Order(allseq *data)
{
  int i,exchange_with;
  phydbl buff_dbl;
  char *buff_name,*buff_seq;
  short int *buff_ambigu;
  
  exchange_with = -1;
  For(i,data->n_otu)
    {
      buff_dbl  = rand();
      buff_dbl /= (RAND_MAX+1.);
      buff_dbl *= data->n_otu;
      exchange_with = (int)floor(buff_dbl);
      
      buff_name = data->c_seq[i]->name;
      data->c_seq[i]->name = data->c_seq[exchange_with]->name;
      data->c_seq[exchange_with]->name = buff_name;

      buff_seq = data->c_seq[i]->state;
      data->c_seq[i]->state = data->c_seq[exchange_with]->state;
      data->c_seq[exchange_with]->state = buff_seq;

      buff_ambigu = data->c_seq[i]->is_ambigu;
      data->c_seq[i]->is_ambigu = data->c_seq[exchange_with]->is_ambigu;
      data->c_seq[exchange_with]->is_ambigu = buff_ambigu;
    }
}

/*********************************************************/

void Update_Root_Pos(arbre *tree)
{
  if(tree->n_root_pos > -1.0)
    {
      tree->n_root->l[0] = tree->e_root->l * tree->n_root_pos;
      tree->n_root->l[1] = tree->e_root->l * (1.-tree->n_root_pos);
    }
  else
    {
      tree->n_root->l[0] = tree->e_root->l / 2.;
      tree->n_root->l[1] = tree->e_root->l / 2.;
    }
}

/*********************************************************/

void Add_Root(edge *target, arbre *tree)
{
  printf("\n. Add root on edge %d left = %d right = %d",target->num,target->left->num,target->rght->num); fflush(NULL);
  tree->e_root = target;

  /* Create the root node if it does not exist yet */
  if((!tree->n_root) || (tree->n_root->num != 2*tree->n_otu-2))
    {      
      tree->n_root = (node *)Make_Node_Light(2*tree->n_otu-2);
    }

  /* Set the position of the root */
  tree->n_root->v[0] = tree->e_root->left;
  tree->n_root->v[1] = tree->e_root->rght;

  tree->n_root->b[0] = tree->e_root;
  tree->n_root->b[1] = tree->e_root;

  if(tree->n_root_pos > -1.0)
    {
/*       tree->n_root->l[0] = tree->e_root->l * (tree->n_root_pos/(1.+tree->n_root_pos)); */
/*       tree->n_root->l[1] = tree->e_root->l - tree->n_root->l[0]; */

      tree->n_root->l[0] = tree->e_root->l * tree->n_root_pos;
      tree->n_root->l[1] = tree->e_root->l * (1. - tree->n_root_pos);
    }
  else
    {
      tree->n_root->l[0] = tree->e_root->l / 2.;
      tree->n_root->l[1] = tree->e_root->l / 2.;
      tree->n_root_pos = 0.5;
    }
}

/*********************************************************/
/* Generate a random unrooted tree with 'n_otu' OTUs */
arbre *Generate_Random_Tree_From_Scratch(int n_otu)
{
  arbre *tree;
  int *is_available,*list_of_nodes;
  int i,node_num,step,n_available;
  phydbl curr_tree_height;
  phydbl rand_len;
  phydbl wanted_tree_length,curr_tree_length;  
  phydbl *rates;
  phydbl stick_prob;
  phydbl lim_inf, lim_sup;

  tree = (arbre *)Make_Tree(n_otu);
  Init_Tree(tree,tree->n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);
  Make_Tree_Path(tree);
  Make_List_Of_Reachable_Tips(tree);
  
  rates = (phydbl *)mCalloc(100,sizeof(phydbl));
  is_available  = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  list_of_nodes = (int *)mCalloc(tree->n_otu,    sizeof(int));

  For(i,tree->n_otu) is_available[i]  = 1;
  For(i,tree->n_otu) list_of_nodes[i] = i;
  For(i,tree->n_otu) 
    {
      tree->noeud[i]->t = .0;
      strcpy(tree->noeud[i]->name,"s");
      sprintf(tree->noeud[i]->name+1,"%d",i);
    }

  curr_tree_height = .0;
  step = 0;
  rand_len = .0;
  do
    {
      node_num = (int)rint(rand()/(phydbl)(RAND_MAX+1.0)*(tree->n_otu-1-step));
      node_num = list_of_nodes[node_num];
      is_available[node_num] = 0;
      For(i,tree->n_otu) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}

      tree->noeud[node_num]->v[0] = tree->noeud[tree->n_otu+step];
      tree->noeud[tree->n_otu+step]->v[1] = tree->noeud[node_num];

      node_num = (int)rint(rand()/(phydbl)(RAND_MAX+1.0)*(tree->n_otu-2-step));
      node_num = list_of_nodes[node_num];
      is_available[node_num] = 0;
      For(i,tree->n_otu) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}

      tree->noeud[node_num]->v[0] = tree->noeud[tree->n_otu+step];
      tree->noeud[tree->n_otu+step]->v[2] = tree->noeud[node_num];

      rand_len = rand();
      rand_len /= RAND_MAX;
      rand_len = -log(rand_len)/(tree->n_otu-step);
      tree->noeud[tree->n_otu+step]->t = curr_tree_height + rand_len;
      curr_tree_height += rand_len;

      is_available[tree->n_otu+step] = 1;
      For(i,tree->n_otu) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) list_of_nodes[n_available++] = i;

      step++;
    }while(step < tree->n_otu-2);

  tree->noeud[list_of_nodes[0]]->v[0] = tree->noeud[list_of_nodes[1]];
  tree->noeud[list_of_nodes[1]]->v[0] = tree->noeud[list_of_nodes[0]];

  tree->num_curr_branch_available = 0;
  Connect_Edges_To_Nodes_Recur(tree->noeud[0],tree->noeud[0]->v[0],tree);
  Add_Root(tree->noeud[list_of_nodes[0]]->b[0],tree);
  rand_len = rand();
  rand_len /= RAND_MAX;
  rand_len = -log(rand_len)/2;
  tree->n_root->t = curr_tree_height + rand_len;
  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  MC_Bl_From_T(tree);

  /* Scale branch lengths such that the sum of branch lengths is uniformly
     distributed in the [lim_inf,lim_sup] range */
  lim_inf =  2.0;
  lim_sup = 30.0;
  wanted_tree_length = rand();
  wanted_tree_length /= RAND_MAX;
  wanted_tree_length *= (lim_sup - lim_inf);
  wanted_tree_length += lim_inf;

  curr_tree_length = .0;
  For(i,2*tree->n_otu-3) curr_tree_length += tree->t_edges[i]->l;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l *= wanted_tree_length/curr_tree_length;

  rates[0] =  0.20;
  rates[1] =  0.50;
  rates[2] =  1.00;
  rates[3] =  2.00;
  rates[4] =  5.00;

  if(!(tree->e_root->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(tree->e_root);
  strcpy(tree->e_root->labels[tree->e_root->n_labels],"ROOT");
  tree->e_root->n_labels++;
  if(!(tree->e_root->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(tree->e_root);
  strcpy(tree->e_root->labels[tree->e_root->n_labels],"MEDIUM");
  tree->e_root->n_labels++;

  stick_prob = .8;
  Random_Lineage_Rates(tree->n_root,tree->n_root->v[0],NULL,stick_prob,rates,2,5,tree);
  Random_Lineage_Rates(tree->n_root,tree->n_root->v[1],NULL,stick_prob,rates,2,5,tree);

  tree->n_root = NULL;
  tree->e_root = NULL;
  tree->print_labels = 1;
  printf("\n\n%s\n",Write_Tree(tree));
  tree->print_labels = 0;
  printf("\n\n%s\n",Write_Tree(tree));
  Exit("\n");

  Free(is_available);
  Free(list_of_nodes);

  return tree;
}

/*********************************************************/

void Random_Lineage_Rates(node *a, node *d, edge *b, phydbl stick_prob, phydbl *rates, int curr_rate, int n_rates, arbre *tree)
{
  phydbl uni;
  int new_rate;
  int i;


  if(b)
    {
      uni  = rand();
      uni /= RAND_MAX;
      
      if(uni > stick_prob) /* Randomly pick a new rate */
	{
	  uni  = rand();
	  uni /= RAND_MAX;
	  uni = (phydbl)(uni * (n_rates-1));	  
	  if(uni-(int)(uni) > 0.5-MDBL_MAX) new_rate = (int)(uni)+1;
	  else new_rate = (int)(uni);	  
	}
      else
	{
	  new_rate = curr_rate;
	}

      For(i,3) 
	if(a->v[i] == d) 
	  {
	    a->b[i]->l *= rates[new_rate];
	    break;
	  }

      For(i,3)
	if(a->v[i] == d)
	  {
	    if(!(a->b[i]->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(a->b[i]);
	    if(rates[new_rate] > 1.0)      strcpy(a->b[i]->labels[a->b[i]->n_labels],"FAST");
	    else if(rates[new_rate] < 1.0) strcpy(a->b[i]->labels[a->b[i]->n_labels],"SLOW");
	    else                           strcpy(a->b[i]->labels[a->b[i]->n_labels],"MEDIUM");
	    a->b[i]->n_labels++;
	    break;
	  }
      curr_rate = new_rate;
    }
  
  if(d->tax) return;
  else
    {
      For(i,3) 
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  Random_Lineage_Rates(d,d->v[i],d->b[i],stick_prob,rates,curr_rate,n_rates,tree);
    }
}

/*********************************************************/

edge *Find_Edge_With_Label(char *label, arbre *tree)
{
  int i,j;

  For(i,2*tree->n_otu-3)
    {
      For(j,tree->t_edges[i]->n_labels)
	{
	  if(!strcmp(tree->t_edges[i]->labels[j],label)) return tree->t_edges[i];
	}
    }
  return NULL;
}

/*********************************************************/

void Print_Square_Matrix_Generic(int n, phydbl *mat)
{
  int i,j;

  printf("\n");
  For(i,n)
    {
      For(j,n)
	{
	  printf("%.3f ",mat[i*n+j]);
	}
      printf("\n");
    }
  printf("\n");
}

/*********************************************************/

void Evolve(allseq *data, model *mod, arbre *tree)
{
  int root_state, root_rate_class;
  int site,i;

  if(mod->use_m4mod) tree->print_labels = 1;
  
  /* Get the change probability matrices */
  Set_Model_Parameters(mod);
  For(i,2*tree->n_otu-3) Update_PMat_At_Given_Edge(tree->t_edges[i],tree);

  For(site,data->init_len)
    {
      root_state = root_rate_class = -1;

      /* Pick the root nucleotide/aa */
      root_state = Pick_State(mod->ns,mod->pi);
      data->c_seq[0]->state[site] = Reciproc_Assign_State(root_state,mod->datatype);

      /* Pick the rate class */
      root_rate_class = Pick_State(mod->n_catg,mod->gamma_r_proba);

      /* tree->noeud[0] is considered as the root node */
      Evolve_Recur(tree->noeud[0],
		   tree->noeud[0]->v[0],
		   tree->noeud[0]->b[0],
		   root_state,
		   root_rate_class,
		   site,
		   data,
		   mod,
		   tree);
      
/*       printf("%s\n",Write_Tree(tree)); */
      
      data->wght[site] = 1;
    }
  data->crunch_len = data->init_len;
}

/*********************************************************/

int Pick_State(int n, phydbl *prob)
{
  int pos;
  phydbl uni;

  do
    {
      pos  = rand();
      pos  = (pos % n);
      uni  = (phydbl)rand();
      uni /= (phydbl)RAND_MAX;
      if(uni < prob[pos]) break;
    }
  while(1);
  
  return (int)pos;
}

/*********************************************************/

void Evolve_Recur(node *a, node *d, edge *b, int a_state, int r_class, int site_num, allseq *gen_data, model *mod, arbre *tree)
{
  int d_state;

  d_state = Pick_State(mod->ns,b->Pij_rr[r_class][a_state]);
  
/*   printf("\n>> %c (%d,%d)",Reciproc_Assign_State(d_state,mod->datatype),d_state,(int)d_state/mod->m4mod->n_o); */

  if(mod->use_m4mod) 
    {
      phydbl rrate; /* relative rate of substitutions */
      
      rrate = mod->m4mod->multipl[(int)d_state/mod->m4mod->n_o];
      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
      if(rrate > 1.0) strcpy(b->labels[b->n_labels],"FASTER");
      else strcpy(b->labels[b->n_labels],"SLOWER");
      b->n_labels++;
    }

  if(d->tax) 
    {
      gen_data->c_seq[d->num]->state[site_num] = Reciproc_Assign_State(d_state,mod->datatype);
      return;
    }
  else
    {
      int i;
      For(i,3)
	if(d->v[i] != a)
	  Evolve_Recur(d,d->v[i],d->b[i],
		       d_state,r_class,site_num,gen_data,
		       mod,tree);
    }
}

/*********************************************************/

void Site_Diversity(arbre *tree)
{
  int i,j,k,ns;
  int *div,sum;

  ns = (tree->mod->datatype == NT)?(4):(20);

  div = (int *)mCalloc(ns,sizeof(int));

  Site_Diversity_Post(tree->noeud[0],tree->noeud[0]->v[0],tree->noeud[0]->b[0],tree);
  Site_Diversity_Pre (tree->noeud[0],tree->noeud[0]->v[0],tree->noeud[0]->b[0],tree);

  For(i,2*tree->n_otu-3)
    {
      For(j,ns)
	{
	  tree->t_edges[i]->div_post_pred_left[j] = 0;
	  tree->t_edges[i]->div_post_pred_rght[j] = 0;
	}
    }

  For(i,tree->n_pattern)
    {
      For(j,2*tree->n_otu-3)
	{
	  Binary_Decomposition(tree->t_edges[j]->ui_l[i],div,ns);
	  sum = 0;
	  For(k,ns) sum += div[k];
	  tree->t_edges[j]->div_post_pred_left[sum-1] += tree->data->wght[i];

	  Binary_Decomposition(tree->t_edges[j]->ui_r[i],div,ns);
	  sum = 0;
	  For(k,ns) sum += div[k];
	  tree->t_edges[j]->div_post_pred_rght[sum-1] += tree->data->wght[i];
	}
    }

/*   For(j,2*tree->n_otu-3) */
/*     { */
/*       printf("\n. Edge %4d   div_left = %4d %4d %4d %4d -- div_rght = %4d %4d %4d %4d", */
/* 	     j, */
/* 	     tree->t_edges[j]->div_post_pred_left[0], */
/* 	     tree->t_edges[j]->div_post_pred_left[1], */
/* 	     tree->t_edges[j]->div_post_pred_left[2], */
/* 	     tree->t_edges[j]->div_post_pred_left[3], */
/* 	     tree->t_edges[j]->div_post_pred_rght[0], */
/* 	     tree->t_edges[j]->div_post_pred_rght[1], */
/* 	     tree->t_edges[j]->div_post_pred_rght[2], */
/* 	     tree->t_edges[j]->div_post_pred_rght[3]); */
/*     } */
  
  Free(div);
}

/*********************************************************/

void Site_Diversity_Post(node *a, node *d, edge *b, arbre *tree)
{
  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	if(d->v[i] != a)
	  Site_Diversity_Post(d,d->v[i],d->b[i],tree);

      Subtree_Union(d,b,tree);
    }
}

/*********************************************************/

void Site_Diversity_Pre(node *a, node *d, edge *b, arbre *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      
      For(i,3)
	if(d->v[i] != a)
	  {
	    Subtree_Union(d,d->b[i],tree);
	    Site_Diversity_Pre(d,d->v[i],d->b[i],tree);
	  }
    }
}

/*********************************************************/

void Subtree_Union(node *n, edge *b_fcus, arbre *tree)
{
/*  
           |
	   |<- b_cus
	   |
	   n
          / \
       	 /   \
       	/     \
*/

  int site;
  unsigned int *ui, *ui_v1, *ui_v2;
  
  ui = ui_v1 = ui_v2 = NULL;

  if(n == b_fcus->left)
    {	     
      ui = b_fcus->ui_l;

      ui_v1 = 
      (n == n->b[b_fcus->l_v1]->left)?
      (n->b[b_fcus->l_v1]->ui_r):
      (n->b[b_fcus->l_v1]->ui_l);

      ui_v2 = 
      (n == n->b[b_fcus->l_v2]->left)?
      (n->b[b_fcus->l_v2]->ui_r):
      (n->b[b_fcus->l_v2]->ui_l);
    }
  else
    {
      ui = b_fcus->ui_r;
      
      ui_v1 = 
      (n == n->b[b_fcus->r_v1]->left)?
      (n->b[b_fcus->r_v1]->ui_r):
      (n->b[b_fcus->r_v1]->ui_l);

      ui_v2 = 
      (n == n->b[b_fcus->r_v2]->left)?
      (n->b[b_fcus->r_v2]->ui_r):
      (n->b[b_fcus->r_v2]->ui_l);
    }

  For(site,tree->n_pattern) ui[site] = ui_v1[site] | ui_v2[site];

}

/*********************************************************/

void Binary_Decomposition(int value, int *bit_vect, int size)
{
  int i,cumul;

  For(i,size) bit_vect[i] = 0;
  
  cumul = 0;
  for(i=size-1;i>=0;i--)
    {
      if(value - cumul < (int)pow(2,i))
	{
	  bit_vect[i] = 0;
	}
      else
	{
	  bit_vect[i] = 1;
	  cumul += (int)pow(2,i);
	}
    }
}

/*********************************************************/

void Print_Diversity_Header(FILE *fp, arbre *tree)
{
/*   fprintf(fp,"edge side mean\n");  */
  fprintf(fp,"edge side diversity count\n"); 
}

/*********************************************************/

void Print_Diversity(FILE *fp, arbre *tree)
{
  int ns;
  
  ns = (tree->mod->datatype == NT)?(4):(20);

  Print_Diversity_Pre(tree->noeud[0],
		      tree->noeud[0]->v[0],
		      tree->noeud[0]->b[0],
		      fp,
		      tree);

/*       mean_div_left = .0; */
/*       For(k,ns)  */
/* 	{ */
/* 	  mean_div_left += (k+1) * tree->t_edges[j]->div_post_pred_left[k]; */
/* 	} */
/*       mean_div_rght = .0; */
/*       For(k,ns) mean_div_rght += (k+1) * tree->t_edges[j]->div_post_pred_rght[k]; */

/*       mean_div_left /= (phydbl)tree->data->init_len; */
/*       mean_div_rght /= (phydbl)tree->data->init_len; */

/*       fprintf(fp,"%4d 0 %f\n",j,mean_div_left); */
/*       fprintf(fp,"%4d 1 %f\n",j,mean_div_rght); */


/*       mean_div_left = .0; */
/*       For(k,ns) mean_div_left += tree->t_edges[j]->div_post_pred_left[k]; */

/*       mean_div_rght = .0; */
/*       For(k,ns)  */
/* 	{ */
/* 	  mean_div_rght += tree->t_edges[j]->div_post_pred_rght[k]; */
/* 	} */

/*       if((mean_div_left != tree->data->init_len) || (mean_div_rght != tree->data->init_len)) */
/* 	{ */
/* 	  printf("\n. mean_div_left = %f mean_div_rght = %f init_len = %d", */
/* 		 mean_div_left,mean_div_rght,tree->data->init_len); */
/* 	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/* 	  Warn_And_Exit(""); */
/* 	} */
}

/*********************************************************/

void Print_Diversity_Pre(node *a, node *d, edge *b, FILE *fp, arbre *tree)
{
  int k,ns;


  if(d->tax) return;
  else
    {
      ns = (tree->mod->datatype == NT)?(4):(20);
      if(d == b->left) For(k,ns) fprintf(fp,"%4d 0 %2d %4d\n",b->num,k,b->div_post_pred_left[k]);
      else             For(k,ns) fprintf(fp,"%4d 1 %2d %4d\n",b->num,k,b->div_post_pred_rght[k]);

      For(k,3) if(d->v[k] != a) Print_Diversity_Pre(d,d->v[k],d->b[k],fp,tree);
    }

}

/*********************************************************/
/* Estimation of density using kernel smoothing. 
- where : point where I want to estimate the density,
- x : data vector, 
- sample_size :  number of data points in x
*/
phydbl Univariate_Kernel_Density_Estimate(phydbl where, phydbl *x, int sample_size)
{
  phydbl sd,h;
  phydbl density,sqrt2pi,cons;
  int i;

  sqrt2pi = 2.506628;

  sd = sqrt(Var(x,sample_size));
  h = 1.06 * sd * pow(sample_size,-1./5.); /* Quick and dirty way to set the bandwidth */
  
  cons = (1./sample_size) * (1./h) * (1./sqrt2pi);

  density = .0;
  For(i,sample_size) density += exp(-0.5 * pow((x[i] - where)/h,2));
  density *= cons;

  return density;
}

/*********************************************************/

/* Estimation of a multivariate density using kernel smoothing. 

- where : vector where I want to estimate the density,
- x : data matrix, i.e., sample of vectors, 
- sample_size : number of vectors,
- vect_size : vector length. 

See "Multivariate Density Estimation" by David Scott. pp 150.
*/
phydbl Multivariate_Kernel_Density_Estimate(phydbl *where, phydbl **x, int sample_size, int vect_size)
{
  phydbl sd,*h,cons,density,tmp;
  phydbl _2pi;
  int i,j;

  h = (phydbl *)mCalloc(vect_size,sizeof(phydbl));

  _2pi = 6.283185;
  
  For(i,vect_size)
    {
      sd = sqrt(Var(x[i],sample_size));
      h[i] = sd * pow(4./((vect_size+2.)*sample_size),1./(vect_size+4));
    }

  cons = sample_size;
  For(i,vect_size) cons *= h[i];
  cons *= pow(_2pi,vect_size/2.);
  cons = 1./cons;

  density = .0;
  For(i,sample_size)
    {
      tmp = 1.0;
      For(j,vect_size) tmp *= exp(-0.5 * pow((x[j][i] - where[j])/h[j],2));
      density += tmp;
    }
  
  density *= cons;

  Free(h);

  return density;
}

/*********************************************************/

phydbl Var(phydbl *x, int n)
{
  phydbl mean, sum2;
  int i;

  mean = Mean(x,n);

  sum2 = .0;
  For(i,n) sum2 += pow(x[i],2);
  
  return (1./n) * (sum2 - n * pow(mean,2));
}

/*********************************************************/

phydbl Mean(phydbl *x, int n)
{
  int i;
  phydbl sum;

  sum = .0;

  For(i,n) sum += x[i];

  return sum / n;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/

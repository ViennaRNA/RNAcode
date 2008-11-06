/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "utilities.h"
#include "optimiz.h"
#include "lk.h"
#include "free.h"
#include "models.h"
#include "mc.h"

/*********************************************************/

void Optimize_Single_Param_Generic(arbre *tree, phydbl *param, phydbl lim_inf, phydbl lim_sup, phydbl tol, int n_max_iter)
{
  phydbl ax,bx,cx;
  phydbl lk_init;
  
  lk_init = tree->c_lnL;

  ax =  lim_inf;
  bx = (*param);
  cx =  lim_sup;
  
  Generic_Brent(ax,bx,cx,tol,param,tree,n_max_iter);

  if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_global) 
    {
      printf("\n. %.10f < %.10f --> diff=%.10f param value = %f initial value = %f\n",
	     tree->c_lnL,lk_init,
	     tree->c_lnL-lk_init,
	     *param,bx);
      Exit("\n. Optimisation failed !\n");
    }
}

/*********************************************************/

int Generic_Brak(phydbl *param,
		 phydbl *ax, phydbl *bx, phydbl *cx, 
		 phydbl *fa, phydbl *fb, phydbl *fc,
		 phydbl lim_inf, phydbl lim_sup,
		 arbre *tree)
{
   phydbl ulim,u,r,q,fu,dum;

   u = 0.0;
   *param = *ax;

   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fa=-Return_Lk(tree);
   *param = *bx;
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fb=-Return_Lk(tree);
   if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
   }
   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   *param = fabs(*cx);
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fc=-Return_Lk(tree); 
   while (*fb > *fc) 
     {
        
       if(*ax > lim_sup) *ax = lim_sup;
       if(*ax < lim_inf) *ax = lim_inf;
       if(*bx > lim_sup) *bx = lim_sup;
       if(*bx < lim_inf) *bx = lim_inf;
       if(*cx > lim_sup) *cx = lim_sup;
       if(*cx < lim_inf) *cx = lim_inf;
       if(u   > lim_sup) u   = lim_sup;
       if(u   < lim_inf) u   = lim_inf;

       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(fabs(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > lim_inf) 
	 {
	   *param = fabs(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Return_Lk(tree);
	   if (fu < *fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       *fa=(*fb);
	       *fb=fu;
	       (*ax)=fabs(*ax);
	       (*bx)=fabs(*bx);
	       (*cx)=fabs(*cx);
	       return(0);
	     } 
	   else if (fu > *fb) 
	     {
	       *cx=u;
	       *fc=fu;	
	       (*ax)=fabs(*ax);
	       (*bx)=fabs(*bx);
	       (*cx)=fabs(*cx);
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = fabs(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Return_Lk(tree);
	 } 
       else if ((*cx-u)*(u-ulim) > lim_inf) 
	 {
	   *param = fabs(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Return_Lk(tree);
	   if (fu < *fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
	       *param = fabs(u); 
	       SHFT(*fb,*fc,fu,-Return_Lk(tree))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= lim_inf) 
	 {
	   u=ulim;
	   *param = fabs(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Return_Lk(tree);
	 } 
       else 
	 {
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = fabs(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Return_Lk(tree);
	 }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)


     }
   (*ax)=fabs(*ax);
   (*bx)=fabs(*bx);
   (*cx)=fabs(*cx);
   return(0);
}

/*********************************************************/

phydbl Generic_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		     phydbl *xmin, arbre *tree, int n_iter_max)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL;
  phydbl best_fu, best_x;

  
  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  (*xmin) = fabs(bx);
  fw=fv=fx=-Return_Lk(tree);
  init_lnL = -fw;
  best_fu = fw;
  best_x  = fabs(bx);;

/*   printf("init_lnL = %f\n",init_lnL); */

  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+BRENT_ZEPS);
      if(
	 ((fabs(tree->c_lnL-old_lnL) < tol) && (tree->c_lnL > init_lnL - tol)) || (iter > n_iter_max - 1))
	{
	  (*xmin) = best_x;
	  Lk(tree);	  
/* 	  printf("\n> iter=%3d max=%3d v=%f lnL=%f init_lnL=%f",iter,n_iter_max,(*xmin),tree->c_lnL,init_lnL); */
	  return tree->c_lnL;
	}
      
      if(fabs(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=fabs(q);
	  etemp=e;
	  e=d;
	  if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    {
	      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	      /*                   printf("Golden section step\n"); */
	    }
	  else
	    {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x);
	      /*                   printf("Parabolic step\n"); */
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  /*               printf("Golden section step (default)\n"); */
	}
      
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*xmin) = fabs(u);
      old_lnL = tree->c_lnL;
      fu = -Return_Lk(tree);
      
      if(fu < best_fu)
	{
	  best_fu = fu;
	  best_x = fabs(u);
	}

/*       printf("\n. iter=%d/%d param=%f loglk=%f (best=%f %f)",iter,BRENT_ITMAX,*xmin,tree->c_lnL,best_fu,best_x); */

      if(fu <= fx) 
	{
	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) 
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
	  else if (fu <= fv || v == x || v == w) 
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }

  Exit("\n. Too many iterations in BRENT !");
  return(-1);
  /* Not Reached ??  *xmin=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/

phydbl RRparam_GTR_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			  phydbl *xmin, arbre *tree, allseq *alldata, phydbl *param, int n_iter_max)
{
   phydbl f1,f2,x0,x1,x2,x3;
   int n_iter;


   x0=ax;
   x3=cx;
   if (fabs(cx-bx) > fabs(bx-ax)) 
     {
       x1=bx;
       x2=bx+GOLDEN_C*(cx-bx);
     } 
   else 
     {
       x2=bx;
       x1=bx-GOLDEN_C*(bx-ax);
     }
   (*param)=x1;

   Lk(tree);
   f1=-tree->c_lnL;
   (*param)=x2;

   Lk(tree);
   f2=-tree->c_lnL;

   n_iter = 0;
   while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) 
     {

       if (f2 < f1) 
	 {
	   SHFT3(x0,x1,x2,GOLDEN_R*x1+GOLDEN_C*x3)
	   (*param)=x2;
	   Lk(tree);
	   SHFT2(f1,f2,-tree->c_lnL)
	 } 
       else 
	 {
	   SHFT3(x3,x2,x1,GOLDEN_R*x2+GOLDEN_C*x0)
	   (*param)=x1;
	   Lk(tree);
	   SHFT2(f2,f1,-tree->c_lnL)
	 }
       
       if(n_iter++ > n_iter_max) break;
       
/*        printf("p=%E %f\n",(*param),tree->c_lnL); */
     }
   if (f1 < f2) 
    {
       *xmin=x1;
       return f1;
     } 
   else 
     {
       *xmin=x2;
       return f2;
     }
}

/*********************************************************/

phydbl Br_Len_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		     phydbl *xmin, edge *b_fcus, arbre *tree)
{
   phydbl f1,f2,x0,x1,x2,x3;

   x0=ax;
   x3=cx;
   if (fabs(cx-bx) > fabs(bx-ax)) 
     {
       x1=bx;
       x2=bx+GOLDEN_C*(cx-bx);
     } 
   else 
     {
       x2=bx;
       x1=bx-GOLDEN_C*(bx-ax);
     }
   
   b_fcus->l=x1;
   f1 = -Lk_At_Given_Edge(b_fcus,tree);
   b_fcus->l=x2;
   f2 = -Lk_At_Given_Edge(b_fcus,tree);
   while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) 
     {
       if (f2 < f1) 
	 {
	   SHFT3(x0,x1,x2,GOLDEN_R*x1+GOLDEN_C*x3)
	   b_fcus->l=x2;
	   SHFT2(f1,f2,-Lk_At_Given_Edge(b_fcus,tree))
	 } 
       else 
	 {
	   SHFT3(x3,x2,x1,GOLDEN_R*x2+GOLDEN_C*x0)
	   b_fcus->l=x1;
	   SHFT2(f2,f1,-Lk_At_Given_Edge(b_fcus,tree))
	 }
     }
   if (f1 < f2) 
     {
       *xmin=fabs(x1);
       return -f1;
     } 
   else 
     {
       *xmin=fabs(x2);
       return -f2;
     }
}

/*********************************************************/

int Br_Len_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		phydbl *fa, phydbl *fb, phydbl *fc, 
		edge *b_fcus, arbre *tree)
{
   phydbl ulim,u,r,q,fu,dum;

   b_fcus->l = *ax;
   *fa=-Lk_At_Given_Edge(b_fcus,tree);
   b_fcus->l = *bx;
   *fb=-Lk_At_Given_Edge(b_fcus,tree);
   if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
   }
   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   b_fcus->l = *cx;
   *fc=-Lk_At_Given_Edge(b_fcus,tree);
   while (*fb > *fc + tree->mod->s_opt->min_diff_lk_local) 
     {
       printf("fb=%f fc=%f\n",*fb,*fc);
       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(fabs(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > 0.0) 
	 {
	   b_fcus->l = u;
	   fu=-Lk_At_Given_Edge(b_fcus,tree);
	   if (fu < *fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       *fa=(*fb);
	       *fb=fu;
/* 	       (*ax)=fabs(*ax); */
/* 	       (*bx)=fabs(*bx); */
/* 	       (*cx)=fabs(*cx); */
	       return(0);
	     } 
	   else if (fu > *fb) 
	     {
	       *cx=u;
	       *fc=fu;	
/* 	       (*ax)=fabs(*ax); */
/* 	       (*bx)=fabs(*bx); */
/* 	       (*cx)=fabs(*cx); */
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   b_fcus->l = u;
	   fu=-Lk_At_Given_Edge(b_fcus,tree);
	 } 
       else if ((*cx-u)*(u-ulim) > 0.0) 
	 {
	   b_fcus->l = fabs(u);
	   fu=-Lk_At_Given_Edge(b_fcus,tree);
	   if (fu < *fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
	       b_fcus->l = u; 
	       SHFT(*fb,*fc,fu,-Lk_At_Given_Edge(b_fcus,tree))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= 0.0) 
	 {
	   u=ulim;
	   b_fcus->l = u;
	   fu=-Lk_At_Given_Edge(b_fcus,tree);
	 } 
       else 
	 {
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   b_fcus->l = u;
	   fu=-Lk_At_Given_Edge(b_fcus,tree);
	 }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)
      }
   (*ax)=fabs(*ax);
   (*bx)=fabs(*bx);
   (*cx)=fabs(*cx);
   return(0);
}

/*********************************************************/

phydbl Br_Len_Brent_Default(edge *b_fcus, arbre *tree)
{
  return Br_Len_Brent(10.*b_fcus->l,b_fcus->l,.10*b_fcus->l,tree->mod->s_opt->min_diff_lk_local,b_fcus,tree,1000);
}

/*********************************************************/

phydbl Br_Len_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
		    edge *b_fcus, arbre *tree, int n_iter_max)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL, init_lnL;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  b_fcus->l = fabs(bx);
  fw=fv=fx=fu=-Lk_At_Given_Edge(b_fcus,tree);
  init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+BRENT_ZEPS);
      if(
	 ((fabs(tree->c_lnL-old_lnL) < tol) && 
	  (tree->c_lnL > init_lnL - tol)) ||	 
	 (iter > n_iter_max - 1))
	{
	  b_fcus->l=x;
	  Lk_At_Given_Edge(b_fcus,tree);
/* 	  printf("\n. iter=%3d max=%3d l=%f lnL=%f init_lnL=%f",iter,n_iter_max,b_fcus->l,tree->c_lnL,init_lnL); */
	  return tree->c_lnL;
	}
      
      if(fabs(e) > tol1)
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=fabs(q);
	  etemp=e;
	  e=d;
	  if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  else{
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	  }
	}
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	}
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      if(u<BL_MIN) u = BL_MIN;
      b_fcus->l=fabs(u);
      old_lnL = tree->c_lnL;
      fu=-Lk_At_Given_Edge(b_fcus,tree);

/*       printf("\n. BRENT edge %3d l=%f lnL=%20f iter=%3d",b_fcus->num,b_fcus->l,fu,iter); */

      if(fu <= fx)
	{
	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	}
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    }
	  else if (fu <= fv || v == x || v == w) 
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }
  if(iter > BRENT_ITMAX) printf("\n. Too many iterations in BRENT (%d) (%f)",iter,b_fcus->l);
  return(-1);
  /* Not Reached ??  *xmin=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/

phydbl Node_Time_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
		       node *anc, node *des, arbre *tree, int n_iter_max)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL, init_lnL;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  des->t = fabs(bx);
  fw=fv=fx=fu=-Lk_Triplet(anc,des,tree);
  init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+BRENT_ZEPS);
      if(
	 ((fabs(tree->c_lnL-old_lnL) < tol) && 
	  (tree->c_lnL > init_lnL - tol)) ||	 
	 (iter > n_iter_max - 1))
	{
	  des->t=x;
	  Lk_Triplet(anc,des,tree);
/* 	  printf("\n. iter=%3d max=%3d l=%f lnL=%f init_lnL=%f",iter,n_iter_max,des->t,tree->c_lnL,init_lnL); */
	  return tree->c_lnL;
	}
      
      if(fabs(e) > tol1)
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=fabs(q);
	  etemp=e;
	  e=d;
	  if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  else{
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	  }
	}
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	}
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      if(u<BL_MIN) u = BL_MIN;
      des->t=fabs(u);
      old_lnL = tree->c_lnL;
      fu=-Lk_Triplet(anc,des,tree);

/*       printf("\n. BRENT node %3d t=%f lnL=%20f iter=%3d",des->num,des->t,fu,iter); */

      if(fu <= fx)
	{
	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	}
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    }
	  else if (fu <= fv || v == x || v == w) 
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }
  if(iter > BRENT_ITMAX) printf("\n. Too many iterations in BRENT (%d) (%f)",iter,des->t);
  return(-1);
  /* Not Reached ??  *xmin=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/

phydbl Time_Stamps_Mult_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
			      arbre *tree, int n_iter_max)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL, init_lnL;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  tree->mod->s_opt->tree_size_mult = fabs(bx);
  MC_Mult_Time_Stamps(tree);
  tree->both_sides = 1;
  fw=fv=fx=fu=-Return_Lk(tree);
  MC_Div_Time_Stamps(tree);
  init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+BRENT_ZEPS);
      if(
	 ((fabs(tree->c_lnL-old_lnL) < tol) && 
	  (tree->c_lnL > init_lnL - tol)) ||	 
	 (iter > n_iter_max - 1))
	{
	  tree->mod->s_opt->tree_size_mult=fabs(x);
	  MC_Mult_Time_Stamps(tree);
	  tree->both_sides = 1;
	  Lk(tree);
	  tree->mod->s_opt->tree_size_mult = 1.0;
	  return tree->c_lnL;
	}
      
      if(fabs(e) > tol1)
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=fabs(q);
	  etemp=e;
	  e=d;
	  if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  else{
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	  }
	}
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	}
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      old_lnL = tree->c_lnL;

      tree->mod->s_opt->tree_size_mult = fabs(u);
      MC_Mult_Time_Stamps(tree);
      tree->both_sides = 1;
      fu=-Return_Lk(tree);
      MC_Div_Time_Stamps(tree);

/*       printf("\n. BRENT mult=%f lnL=%20f", */
/* 	     tree->mod->s_opt->tree_size_mult, */
/* 	     tree->c_lnL); */

      if(fu <= fx)
	{
	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	}
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    }
	  else if (fu <= fv || v == x || v == w) 
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }
  tree->mod->s_opt->tree_size_mult = 1.0;
  if(iter > BRENT_ITMAX) printf("\n. Too many iterations in BRENT (%d)",iter);
  return(-1);
  /* Not Reached ??  *xmin=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/

phydbl Branch_Rate_Shape_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			       phydbl *xmin, arbre *tree, int n_iter_max)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_fu,init_fu;
  phydbl best_fu, best_x;

  
  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_fu = UNLIKELY;
  (*xmin) = fabs(bx);
  fw=fv=fx=-Lk_With_MAP_Branch_Rates(tree);
  init_fu = fw;
  best_fu = fw;
  fu      = fw;
  best_x  = fabs(bx);;

  printf("init_lnL = %f\n",fw);

  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+BRENT_ZEPS);
      if(
	 ((fabs(fu-old_fu) < tol) && (fu < init_fu-tol)) || (iter > n_iter_max - 1))
	{
	  (*xmin) = best_x;
	  printf("\n> iter=%d/%d param val=%f fu=%f (best_fu=%f best_x=%f old_fu=%f)",iter,BRENT_ITMAX,*xmin,fu,best_fu,best_x,old_fu);
	  return Lk_With_MAP_Branch_Rates(tree);
	}
      
      if(fabs(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=fabs(q);
	  etemp=e;
	  e=d;
	  if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    {
	      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	      /*                   printf("Golden section step\n"); */
	    }
	  else
	    {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x);
	      /*                   printf("Parabolic step\n"); */
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  /*               printf("Golden section step (default)\n"); */
	}
      
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*xmin) = fabs(u);
      old_fu = fu;
      fu = -Lk_With_MAP_Branch_Rates(tree);
      
      if(fu < best_fu)
	{
	  best_fu = fu;
	  best_x = fabs(u);
	}

      printf("\n. iter=%d/%d param val=%f fu=%f (best_fu=%f best_x=%f old_fu=%f)",iter,BRENT_ITMAX,*xmin,fu,best_fu,best_x,old_fu);

      if(fu <= fx) 
	{
	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) 
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
	  else if (fu <= fv || v == x || v == w) 
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }

  Exit("\n. Too many iterations in BRENT !");
  return(-1);
  /* Not Reached ??  *xmin=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/

int Dist_Seq_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		  phydbl *fa, phydbl *fb, phydbl *fc, 
		  allseq *data, int numseq1, int numseq2, 
		  model *mod)
{
   phydbl ulim,u,r,q,fu,dum;
   phydbl dist;
   phydbl lk;

   dist = *ax;
   *fa=-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk);
   dist = *bx;
   *fb=-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk);
   if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
   }
   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   dist = fabs(*cx);
   *fc=-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk);
   while (*fb > *fc) 
     {
       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(fabs(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > 0.0) 
	 {
	   dist = fabs(u);
	   fu=-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk);
	   if (fu < *fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       *fa=(*fb);
	       *fb=fu;
	       return(0);
	     } 
	   else if (fu > *fb) 
	     {
	       *cx=u;
	       *fc=fu;
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   dist = fabs(u);
	   fu=-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk);
	 } 
       else if ((*cx-u)*(u-ulim) > 0.0) 
	 {
	   dist = fabs(u);
	   fu=-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk);
	   if (fu < *fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
	       dist = fabs(u); 
	       SHFT(*fb,*fc,fu,-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= 0.0) 
	 {
	   u=ulim;
	   dist = fabs(u);
	   fu=-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk);
	 } 
       else 
	 {
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   dist = fabs(u);
	   fu=-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk);
	 }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)
      }
   return(0);
}

/*********************************************************/

phydbl Dist_Seq_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		      phydbl *xmin, allseq *data, 
		      int numseq1, int numseq2, model *mod)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl dist;
  phydbl lk;
  
  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  dist = fabs(bx);
  fw=fv=fx=-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk);
  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+BRENT_ZEPS);
      if(fabs(x-xm) <= (tol2-0.5*(b-a))) 
	{
	  *xmin=x;
	  return -fx;
	}
      
      if(fabs(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=fabs(q);
	  etemp=e;
	  e=d;
	  if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  else{
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	  }
	} else
	  {
	    d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  }
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      dist=fabs(u);
      fu=-Lk_Given_Two_Seq(data,numseq1,numseq2,dist,mod,&lk);
      if(fu <= fx) {
	if(u >= x) a=x; else b=x;
	SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	  } 
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) 
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
    } 
	  else if (fu <= fv || v == x || v == w) {
            v=u;
            fv=fu;
	  }
	}
    }
  printf("\n . BRENT method failed, trying Newton-Raphson");
  return(+1.0);
  /* Not Reached ??  *xmin=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/

phydbl Optimize_Dist(model *mod, phydbl init, allseq *twoseqs)
{
  phydbl d_infa,d_max,d_infb;
  phydbl lk_infa, lk_max, lk_infb, lk;

  d_infa = 100.*BL_MIN;
  d_max  = init;
  d_infb = 3.*init;
  if(init <= BL_MIN) {d_infa = -BL_START; d_max = .0; d_infb = BL_START;}
  lk_infa = lk_max = lk_infb = .0;

  Dist_Seq_Brak(&d_infa, &d_max, &d_infb,
		&lk_infa,&lk_max,&lk_infb,
		twoseqs,0,1,mod);
  
  lk = (phydbl) Dist_Seq_Brent(d_infa,d_max,d_infb,
			       1.e-5,&d_max,twoseqs,0,1,mod);
  if(lk > .0) return -1.0;
  else        return d_max;

}

/*********************************************************/

void Round_Optimize(arbre *tree, allseq *data)
{
  int n_round,each;
  phydbl lk_old, lk_new, tol;
  node *root;

  lk_new = tree->c_lnL;
  lk_old = UNLIKELY;
  n_round = 0;
  each = 1;
  tol = 1.e-2;
  root = tree->noeud[0];
  
  tree->both_sides = 1;
  Lk(tree);

  while(n_round < ROUND_MAX)
    {

      (!((n_round+2)%2))?(root=tree->noeud[0]):(root=tree->noeud[tree->n_otu-1]);
      
      if(tree->mod->s_opt->opt_bl) Optimize_Br_Len_Serie(root,root->v[0],root->b[0],tree,data);
            
      tree->both_sides = 1;
      Lk(tree);
      
      if(tree->mod->s_opt->print) Print_Lk(tree,"[Branch lengths     ]");

      if(!each)
	{
	  each = 1;
	  Optimiz_All_Free_Param(tree,tree->mod->s_opt->print);
	  tree->both_sides = 1;
	  Lk(tree);
	}
      
      lk_new = tree->c_lnL;      
      if(lk_new < lk_old - tree->mod->s_opt->min_diff_lk_global) Exit("\n. Optimisation failed ! (Round_Optimize)\n");
      if(fabs(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_global)  break;
/*       if(fabs(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_local)  break; */
      else lk_old  = lk_new;
      n_round++;
      each--;
    }
  
  Optimiz_All_Free_Param(tree,tree->mod->s_opt->print);
}

/*********************************************************/

void Optimize_Br_Len_Serie(node *a, node *d, edge *b_fcus, arbre *tree, allseq *alldata)
{
  int i;
  phydbl l_infa,l_max,l_infb;
  phydbl lk_init;
  
  
  lk_init = tree->c_lnL;
  
  l_infa = l_max  = l_infb = BL_MIN;
 
  l_infa = 10.*b_fcus->l;
  l_max  = b_fcus->l;
  l_infb = BL_MIN;
  
/*   Br_Len_Brent(l_infa,l_max,l_infb, */
/* 	       1.e-3, */
/* 	       &(b_fcus->l), */
/* 	       b_fcus,tree,500); */

  Br_Len_Brent(l_infa,l_max,l_infb,
	       tree->mod->s_opt->min_diff_lk_local,
	       b_fcus,tree,
	       tree->mod->s_opt->brent_it_max);

  if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local)
    {
      printf("\n. %f %f %f %f",l_infa,l_max,l_infb,b_fcus->l);
      printf("\n. %f -- %f",lk_init,tree->c_lnL);
      Warn_And_Exit("\n. Err. in Optimize_Br_Len_Serie\n");
    }
    
  if(d->tax) return;
  else For(i,3) if(d->v[i] != a)
    {
      Update_P_Lk(tree,d->b[i],d);
      Optimize_Br_Len_Serie(d,d->v[i],d->b[i],tree,alldata);
    }
  For(i,3) if((d->v[i] == a) && !(d->v[i]->tax)) Update_P_Lk(tree,d->b[i],d);
}

/*********************************************************/

void Optimiz_Ext_Br(arbre *tree)
{
  int i;
  edge *b;
  phydbl l_infa,l_max,l_infb;
  phydbl lk, lk_init,l_init;
  
  lk_init = tree->c_lnL;
  

  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      if((b->left->tax) || (b->rght->tax))
	{

	  l_init = b->l;

/* 	  Fast_Br_Len(b,tree); */
/* 	  lk = Lk_At_Given_Edge(tree,b); */

	  l_infa = 10.*b->l;
	  l_max  = b->l;
	  l_infb = BL_MIN;

	  lk = Br_Len_Brent(l_infa,l_max,l_infb,
			    tree->mod->s_opt->min_diff_lk_local,
			    b,tree,
			    tree->mod->s_opt->brent_it_max);

	  b->nni->best_l    = b->l;
	  b->nni->l0        = b->l;
	  b->nni->best_conf = 0;
	  b->l              = l_init;

	}
    }
  tree->c_lnL = lk_init; 
}

/*********************************************************/

void Optimiz_All_Free_Param(arbre *tree, int verbose)
{
  int  init_both_sides;

  init_both_sides  = tree->both_sides;
  tree->both_sides = 0;

  if((tree->mod->whichmodel == GTR) ||
     ((tree->mod->whichmodel == CUSTOM) && (tree->mod->s_opt->opt_rr) && (tree->mod->n_diff_rr > 1)))
    {
      int failed;
      
      failed = 0;
      
      tree->mod->update_eigen = 1;
      BFGS(tree,tree->mod->rr_val,tree->mod->n_diff_rr,1.e-5,1.e-7,
	   &Return_Abs_Lk,
	   &Num_Derivative_Several_Param,
	   &Lnsrch_RR_Param,&failed);
      
      if(failed)
	{
	  int i;
	  
/* 	  if(verbose) printf("\n. Optimising one-by-one...\n");	   */
	  For(i,tree->mod->n_diff_rr) 
	    if(i != 5)
	      {
		Optimize_Single_Param_Generic(tree,&(tree->mod->rr_val[i]),
					      1.E-20,1.E+10,
					      tree->mod->s_opt->min_diff_lk_global,
					      tree->mod->s_opt->brent_it_max);
	      }        
	}
      if(verbose) Print_Lk(tree,"[GTR parameters     ]");
      tree->mod->update_eigen = 0;
    }
  
  if(tree->mod->s_opt->opt_kappa)
    {
      tree->mod->update_eigen = 1;
      Optimize_Single_Param_Generic(tree,&(tree->mod->kappa),
				    .1,100.,
				    tree->mod->s_opt->min_diff_lk_global,
				    tree->mod->s_opt->brent_it_max);
      if(verbose) 
	{
	  Print_Lk(tree,"[Ts/ts ratio        ]");
	  printf("[%10f]",tree->mod->kappa);
	}
      tree->mod->update_eigen = 0;
    }
  
  if(tree->mod->s_opt->opt_lambda) 
    {
      Optimize_Single_Param_Generic(tree,&(tree->mod->lambda),.001,100.,
				    tree->mod->s_opt->min_diff_lk_global,
				    tree->mod->s_opt->brent_it_max);
      if(verbose) 
	{
	  Print_Lk(tree,"[Lambda             ]");
	  printf("[%10f]",tree->mod->lambda);
	}
    }
  
  if((tree->mod->s_opt->opt_pinvar) && (tree->mod->s_opt->opt_alpha))
    {
      Optimiz_Alpha_And_Pinv(tree);
      if(verbose) 
	{
	  Print_Lk(tree,"[Alpha              ]"); 
	  printf("[%10f]",tree->mod->alpha);
	  Print_Lk(tree,"[P-inv              ]"); 
	  printf("[%10f]",tree->mod->pinvar);
	}
    }
  else
    {
      if(tree->mod->s_opt->opt_pinvar)
	{
	  Optimize_Single_Param_Generic(tree,&(tree->mod->pinvar),
					.0001,0.9999,
					tree->mod->s_opt->min_diff_lk_global,
					tree->mod->s_opt->brent_it_max);
	  if(verbose) 
	    {
	      Print_Lk(tree,"[P-inv              ]");
	      printf("[%10f]",tree->mod->pinvar);
	    }
	}
      
      if(tree->mod->s_opt->opt_alpha)
	{
	  Optimize_Single_Param_Generic(tree,&(tree->mod->alpha),
					.01,100.,
					tree->mod->s_opt->min_diff_lk_global,
					tree->mod->s_opt->brent_it_max);
	  if(verbose) 
	    {
	      Print_Lk(tree,"[Alpha              ]");
	      printf("[%10f]",tree->mod->alpha);
	    }
	}
    }

  if((tree->mod->s_opt->opt_state_freq) && (tree->mod->datatype == NT))
    {
        int failed,i;
        
        failed = 0;
        tree->mod->update_eigen = 1;
        BFGS(tree,tree->mod->pi,4,1.e-5,1.e-7,&Return_Abs_Lk,&Num_Derivative_Several_Param,&Lnsrch_Nucleotide_Frequencies,&failed);

        if(failed)
	  {
	    For(i,4) 
	      Optimize_Single_Param_Generic(tree,&(tree->mod->pi_unscaled[i]),
					    -1000.,1000.,
					    tree->mod->s_opt->min_diff_lk_global,
					    tree->mod->s_opt->brent_it_max);
	  }
	if(verbose) Print_Lk(tree,"[Nucleotide freqs.  ]");
        tree->mod->update_eigen = 0;
    }

  if(tree->mod->use_m4mod)
    {
      int failed,i;

      if(tree->mod->s_opt->opt_cov_delta) 
	{
	  tree->mod->update_eigen = 1;
	  Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->delta),
					0.01,10.,
					tree->mod->s_opt->min_diff_lk_global,
					tree->mod->s_opt->brent_it_max);
	  if(verbose) 
	    {
	      Print_Lk(tree,"[Switching param.   ]");
	      printf("[%10f]",tree->mod->m4mod->delta);
	    }
	  tree->mod->update_eigen = 0;
	}
      
      if(tree->mod->s_opt->opt_cov_free_rates) 
	{
	  int rcat;
	  
	  tree->mod->update_eigen = 1;
	  For(rcat,tree->mod->m4mod->n_h)
	    {
	      Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->multipl_unscaled[rcat]),
					    .01,10.,
					    tree->mod->s_opt->min_diff_lk_global,
					    tree->mod->s_opt->brent_it_max);
	      
	      if(verbose) 
		{
		  Print_Lk(tree,"[Rel. subst. rate   ]");
		  printf("[%10f]",tree->mod->m4mod->multipl[rcat]);
		}
	    }
	  
	  For(rcat,tree->mod->m4mod->n_h)
	    {
	      Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->h_fq_unscaled[rcat]),
					    .01,100.,
					    tree->mod->s_opt->min_diff_lk_global,
					    tree->mod->s_opt->brent_it_max);
	      
	      if(verbose)
		{
		  Print_Lk(tree,"[Subst. class freq  ]");
		  printf("[%10f]",tree->mod->m4mod->h_fq[rcat]);
		}
	    }
	  tree->mod->update_eigen = 0;      
	}
      
      if(tree->mod->s_opt->opt_cov_alpha) 
	{
	  tree->mod->update_eigen = 1;
	  Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->alpha),
					.01,10.,
					tree->mod->s_opt->min_diff_lk_global,
					tree->mod->s_opt->brent_it_max);
	  if(verbose) 
	    {
	      Print_Lk(tree,"[Alpha (covarion)   ]");
	      printf("[%10f]",tree->mod->m4mod->alpha);
	    }
	  tree->mod->update_eigen = 0;      	  
	}
      
	  
      /* Substitutions between nucleotides are considered to follow a 
	 GTR model */
      
      failed = 0;
      tree->mod->update_eigen = 1;
      BFGS(tree,tree->mod->m4mod->o_rr,5,1.e-5,1.e-7,
	   &Return_Abs_Lk,
	   &Num_Derivative_Several_Param,
	   &Lnsrch_RR_Cov_Param,&failed);
      
      if(failed)
	{
	  For(i,5) 
	    {
	      Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->o_rr[i]),
					    1.E-20,1.E+10,
					    tree->mod->s_opt->min_diff_lk_global,
					    tree->mod->s_opt->brent_it_max);
	    }
	}
      if(verbose) Print_Lk(tree,"[GTR parameters     ]");
      tree->mod->update_eigen = 0;
    }

  tree->both_sides = init_both_sides;
}



#define ITMAX 200
#define EPS   3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

void BFGS(arbre *tree, phydbl *p, int n, phydbl gtol, phydbl step_size,
	  phydbl(*func)(), void (*dfunc)(), void (*lnsrch)(),int *failed)
{

  int check,i,its,j;
  phydbl den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fret;
  phydbl *dg,*g,*hdg,**hessin,*pnew,*xi;
  
  hessin = (phydbl **)mCalloc(n,sizeof(phydbl *));
  For(i,n) hessin[i] = (phydbl *)mCalloc(n,sizeof(phydbl));
  dg   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  g    = (phydbl *)mCalloc(n,sizeof(phydbl ));
  pnew = (phydbl *)mCalloc(n,sizeof(phydbl ));
  hdg  = (phydbl *)mCalloc(n,sizeof(phydbl ));
  xi   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  

  fp=(*func)(tree);
  (*dfunc)(tree,p,n,step_size,func,g);


  for (i=0;i<n;i++) 
    {
      for (j=0;j<n;j++) hessin[i][j]=0.0;
      hessin[i][i]=1.0;
      xi[i] = -g[i];
      sum += p[i]*p[i];
    }

  stpmax=STPMX*MAX(sqrt(sum),(phydbl)n);

  for(its=1;its<=ITMAX;its++) 
    {
      lnsrch(tree,n,p,fp,g,xi,pnew,&fret,stpmax,&check);

/*       printf("BFGS -> %f\n",tree->c_lnL); */

      fp = fret;
      
      for (i=0;i<n;i++) 
	{
	  xi[i]=pnew[i]-p[i];
	  p[i]=pnew[i];
	}

      test=0.0;
      for (i=0;i<n;i++) 
	{
	  temp=fabs(xi[i])/MAX(fabs(p[i]),1.0);
	  if (temp > test) test=temp;
	}
      if (test < TOLX) 
	{
	  (*func)(tree);
	  For(i,n) Free(hessin[i]);
	  free(hessin);
	  free(xi);
	  free(pnew);
	  free(hdg);
	  free(g);
	  free(dg);   

	  if(its == 1) 
	    {
/* 	      printf("\n. WARNING : BFGS failed ! \n"); */
	      *failed = 1;
	    }
	  return;
	}

      for (i=0;i<n;i++) dg[i]=g[i];

      (*dfunc)(tree,p,n,step_size,func,g);

      test=0.0;
      den=MAX(fret,1.0);
      for (i=0;i<n;i++) 
	{
	  temp=fabs(g[i])*MAX(fabs(p[i]),1.0)/den;
	  if (temp > test) test=temp;
	}
      if (test < gtol) 
	{
	  (*func)(tree);
	  For(i,n) Free(hessin[i]);
	  free(hessin);
	  free(xi);
	  free(pnew);
	  free(hdg);
	  free(g);
	  free(dg);   
	  return;
	}

    for (i=0;i<n;i++) dg[i]=g[i]-dg[i];

    for (i=0;i<n;i++) 
      {
	hdg[i]=0.0;
	for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
      }

    fac=fae=sumdg=sumxi=0.0;
    for (i=0;i<n;i++) 
      {
	fac += dg[i]*xi[i];
	fae += dg[i]*hdg[i];
	sumdg += SQR(dg[i]);
	sumxi += SQR(xi[i]);
      }
    
    if(fac*fac > EPS*sumdg*sumxi) 
      {
	fac=1.0/fac;
	fad=1.0/fae;
	for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
	for (i=0;i<n;i++) 
	  {
	    for (j=0;j<n;j++) 
	      {
		hessin[i][j] += fac*xi[i]*xi[j]
		  -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
	      }
	  }
      }
    for (i=0;i<n;i++) 
      {
	xi[i]=0.0;
	for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
      }
    }
  Exit("\n. Too many iterations in BFGS...\n");
  For(i,n) Free(hessin[i]);
  free(hessin);
  free(xi);
  free(pnew);
  free(hdg);
  free(g);
  free(dg);   
}

#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX

/*********************************************************/


#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_RR_Param(arbre *tree, int n, phydbl *xold, phydbl fold, 
		     phydbl *g, phydbl *p, phydbl *x,
		     phydbl *f, phydbl stpmax, int *check)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=fabs(p[i])/MAX(fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) 
	{
	  x[i]=fabs(local_xold[i]+alam*p[i]);
	}

      /**/
      for(i=0;i<n;i++)
	{
	  tree->mod->rr_val[i] = fabs(local_xold[i]+alam*p[i]);
	}
      /**/

      if(i==n) 
	{
	  *f=Return_Abs_Lk(tree);
/* 	  printf("loglk = %f\n",*f); */
	}
      else *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) 
	    {
	      x[i]=fabs(local_xold[i]);
	    }
	  /**/      
	  for(i=0;i<n;i++)
	    {
	      tree->mod->rr_val[i] = fabs(local_xold[i]);
	    }
	  /**/

	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a == 0.0) tmplam = -slope/(2.0*b);
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) tmplam = 0.5*alam;
		  else if(b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
		  else tmplam = -slope/(b+sqrt(disc));
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_RR_Cov_Param(arbre *tree, int n, phydbl *xold, phydbl fold, 
			 phydbl *g, phydbl *p, phydbl *x,
			 phydbl *f, phydbl stpmax, int *check)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=fabs(p[i])/MAX(fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) 
	{
	  x[i]=fabs(local_xold[i]+alam*p[i]);
	}

      /**/
      for(i=0;i<n;i++)
	{
	  tree->mod->m4mod->o_rr[i] = fabs(local_xold[i]+alam*p[i]);
	}
      /**/

      if(i==n) 
	{
	  *f=Return_Abs_Lk(tree);
/* 	  printf("loglk = %f\n",*f); */
	}
      else *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) 
	    {
	      x[i]=fabs(local_xold[i]);
	    }
	  /**/      
	  for(i=0;i<n;i++)
	    {
	      tree->mod->m4mod->o_rr[i] = fabs(local_xold[i]);
	    }
	  /**/

	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a == 0.0) tmplam = -slope/(2.0*b);
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) tmplam = 0.5*alam;
		  else if(b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
		  else tmplam = -slope/(b+sqrt(disc));
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_Nucleotide_Frequencies(arbre *tree, int n, phydbl *xold, phydbl fold, phydbl *g, phydbl *p, phydbl *x,
				   phydbl *f, phydbl stpmax, int *check)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=fabs(p[i])/MAX(fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) x[i]=fabs(local_xold[i]+alam*p[i]);
      /**/      
      for(i=0;i<n;i++) 
	{
	  tree->mod->pi[i]=fabs(local_xold[i]+alam*p[i]);
/* 	  if( */
/* 	     (tree->mod->pi[i] < 0.001) || */
/* 	     (tree->mod->pi[i] > 0.999) */
/* 	     ) */
/* 	    break; */
	}
      /**/
      if(i==n) 
	{
	  *f=Return_Abs_Lk(tree);
	}
      else     *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) x[i]=local_xold[i];
	  for (i=0;i<n;i++) tree->mod->pi[i]=local_xold[i];
	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a == 0.0) tmplam = -slope/(2.0*b);
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) 
		    {
		      disc=b*b-3.0*a*slope;
		      if (disc<0.0) tmplam = 0.5*alam;
		      else if(b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
		      else tmplam = -slope/(b+sqrt(disc));
		    }
		  else tmplam=(-b+sqrt(disc))/(3.0*a);
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

/*********************************************************/

/* void Optimize_Global_Rate(arbre *tree) */
/* { */
/*     printf("\n. Global rate (%f->)",tree->c_lnL); */
/*     Optimize_Single_Param_Generic(tree,&(tree->tbl),tree->tbl,BL_MIN,1.E+4,100); */
/*     printf("%f)\n",tree->c_lnL); */
/* } */


#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

void EM_Dist(model *mod, allseq *data)
{
  phydbl sum;
  phydbl **prob, *F;
  int i,j,site;
  phydbl d;
  phydbl ***Pij;
  int n_iter;
  phydbl p_diff;
  phydbl **p_lk_left,**p_lk_rght;
  phydbl *pi;


  p_lk_left = (phydbl **)mCalloc(data->c_seq[0]->len,sizeof(phydbl *));
  p_lk_rght = (phydbl **)mCalloc(data->c_seq[0]->len,sizeof(phydbl *));
  F = (phydbl *)mCalloc(mod->ns*mod->ns, sizeof(phydbl));
  pi = (phydbl *)mCalloc(mod->ns, sizeof(phydbl));
  prob = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
  For(i,mod->ns) prob[i] = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  Pij = mod->Pij_rr;

  For(i,data->c_seq[0]->len)
    {
      p_lk_left[i] = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
      p_lk_rght[i] = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
    }

  
  For(i,mod->ns) pi[i] = 0.25;
  d = 0.1;
  n_iter = 0;  
  do
    {
      PMat(d,mod,Pij);

      For(i,mod->ns) For(j,mod->ns) F[mod->ns*i+j] = .0;      

      For(site,data->c_seq[0]->len)
	{
	  For(i,mod->ns) For(j,mod->ns)
	    prob[i][j] =
	    pi[i] *
	    Pij[0][i][j]  *
	    p_lk_left[site][i] *
	    p_lk_rght[site][j] ;
	  
	  sum = .0;
	  For(i,mod->ns) For(j,mod->ns)
	    sum += prob[i][j];
	  For(i,mod->ns) For(j,mod->ns)
	    prob[i][j] /= sum;
	  
	  
	  For(i,mod->ns) For(j,mod->ns)
	    F[mod->ns*i+j] += 
	    data->wght[site] * 
	    prob[i][j];
	}
      
      For(i,mod->ns) For(j,mod->ns) 
	F[mod->ns*i+j] /= (phydbl)data->init_len;
      
      For(i,mod->ns) For(j,mod->ns) 
	F[mod->ns*i+j] = (F[mod->ns*i+j] +F[mod->ns*j+i])/2.;
      
      p_diff = .0;
      For(i,mod->ns) p_diff += F[mod->ns*i+i];
      p_diff = 1. - p_diff;
      d = -(3./4.)*log(1.-(4./3.)*p_diff);
      printf("\n. d = %f\n",d);
      n_iter++;
    }while(n_iter < 10);


  For(i,data->c_seq[0]->len)
    {
      Free(p_lk_left[i]);
      Free(p_lk_rght[i]);
    }
  Free(p_lk_left); Free(p_lk_rght);
  Free(pi);
}

/*********************************************************/

int Dist_F_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *F, phydbl *param, model *mod)
{
   phydbl ulim,u,r,q,dum;
   phydbl fa, fb, fc, fu;

   fa = -Lk_Dist(F,fabs(*ax),mod);
   fb = -Lk_Dist(F,fabs(*bx),mod);

   if(fb > fa) 
     {
       SHFT(dum,*ax,*bx,dum)
       SHFT(dum,fb,fa,dum)
     }

   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   fc = -Lk_Dist(F,fabs(*cx),mod);

   while (fb > fc) 
     {
       r=(*bx-*ax)*(fb-fc);
       q=(*bx-*cx)*(fb-fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(fabs(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > 0.0) 
	 {
	   fu = -Lk_Dist(F,fabs(u),mod);
	   if (fu < fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       fa=fb;
	       fb=fu;
	       return(0);
	     } 
	   else if (fu > fb) 
	     {
	       *cx=u;
	       fc=fu;
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   fu = -Lk_Dist(F,fabs(u),mod);
	 } 
       else if ((*cx-u)*(u-ulim) > 0.0) 
	 {
	   fu = -Lk_Dist(F,fabs(u),mod);
	   if (fu < fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
	       SHFT(fb,fc,fu,-Lk_Dist(F,fabs(u),mod))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= 0.0) 
	 {
	   u  = ulim;
	   fu = -Lk_Dist(F,fabs(u),mod);
	 } 
       else 
	 {
	   u  =(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   fu = -Lk_Dist(F,fabs(u),mod);
	 }

       SHFT(*ax,*bx,*cx,u)
       SHFT(fa,fb,fc,fu)
      }
   return(0);
}

/*********************************************************/

phydbl Dist_F_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
		    phydbl *param, phydbl *F, model *mod)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL, curr_lnL;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x = w = v = bx;
  old_lnL = UNLIKELY;
  fw = fv = fx = -Lk_Dist(F,fabs(bx),mod);
  curr_lnL = init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);

      tol2=2.0*(tol1=tol*fabs(x)+BRENT_ZEPS);

      if(
	 ((fabs(curr_lnL-old_lnL) < mod->s_opt->min_diff_lk_local) && 
	  (curr_lnL > init_lnL - mod->s_opt->min_diff_lk_local)) ||	 
	  (iter > n_iter_max - 1)
	 )	 
	{
	  *param = x;
	  curr_lnL = Lk_Dist(F,*param,mod);
	  return -curr_lnL;
	}
      
      if(fabs(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=fabs(q);
	  etemp=e;
	  e=d;
	  if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    {
	      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	      /*                   printf("Golden section step\n"); */
	    }
	  else
	    {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x);
	      /*                   printf("Parabolic step\n"); */
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  /*               printf("Golden section step (default)\n"); */
	}
      
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      if(u<BL_MIN) u = BL_MIN;
      (*param) = fabs(u);
      old_lnL = curr_lnL;
      fu = -Lk_Dist(F,fabs(u),mod);
      curr_lnL = -fu;      
/*       printf("param=%f loglk=%f\n",*param,fu); */
      
      if(fu <= fx) 
	{
	  if(iter > n_iter_max) return -fu;

	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) 
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
	  else if (fu <= fv || v == x || v == w) {
            v=u;
            fv=fu;
	  }
	}
    }
  Exit("\n. Too many iterations in BRENT !");
  return(-1);
}

/*********************************************************/

void Opt_Dist_F(phydbl *dist, phydbl *F, model *mod)
{
  phydbl ax,bx,cx;

  if(*dist < BL_MIN) *dist = BL_MIN;

  ax = *dist;
  bx = 1.5*(*dist);

  ax = 10.*(*dist);
  bx =     (*dist);
  cx = .10*(*dist);

/*   Dist_F_Brak(&ax,&bx,&cx,F,dist,mod); */
  Dist_F_Brent(ax,bx,cx,1.E-10,1000,dist,F,mod);
}

/*********************************************************/

int Missing_Dist_Brak(phydbl *ax, phydbl *bx, phydbl *cx, int x, int y, matrix *mat)
{
   phydbl ulim,u,r,q,dum;
   phydbl fa, fb, fc, fu;

   fa = Least_Square_Missing_Dist_XY(x,y,fabs(*ax),mat);
   fb = Least_Square_Missing_Dist_XY(x,y,fabs(*bx),mat);

   if(fb > fa) 
     {
       SHFT(dum,*ax,*bx,dum)
       SHFT(dum,fb,fa,dum)
     }

   *cx=(*bx)+MNBRAK_GOLD*((*bx)-(*ax));
   fc = Least_Square_Missing_Dist_XY(x,y,fabs(*cx),mat);

   while (fb > fc) 
     {
       r=((*bx)-(*ax))*(fb-fc);
       q=((*bx)-(*cx))*(fb-fa);
       u=(*bx)-(((*bx)-(*cx))*q-((*bx)-(*ax))*r)/
               (2.0*SIGN(MAX(fabs(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > 0.0) 
	 {
	   fu = Least_Square_Missing_Dist_XY(x,y,fabs(u),mat);
	   if (fu < fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       fa=fb;
	       fb=fu;
	       return(0);
	     } 
	   else if (fu > fb) 
	     {
	       *cx=u;
	       fc=fu;
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   fu = Least_Square_Missing_Dist_XY(x,y,fabs(u),mat);
	 } 
       else if ((*cx-u)*(u-ulim) > 0.0) 
	 {
	   fu = Least_Square_Missing_Dist_XY(x,y,fabs(u),mat);
	   if (fu < fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
		 SHFT(fb,fc,fu,Least_Square_Missing_Dist_XY(x,y,fabs(u),mat))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= 0.0) 
	 {
	   u  = ulim;
	   fu = Least_Square_Missing_Dist_XY(x,y,fabs(u),mat);
	 } 
       else 
	 {
	   u  =(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   fu = Least_Square_Missing_Dist_XY(x,y,fabs(u),mat);
	 }

       SHFT(*ax,*bx,*cx,u)
       SHFT(fa,fb,fc,fu)
      }
   return(0);
}

/*********************************************************/

phydbl Missing_Dist_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
			  int x, int y, matrix *mat)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,xx,xm;
  phydbl e=0.0;
  phydbl init_loglk, max_loglk;
  phydbl bestx;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  xx=w=v=bx;
  fx=Least_Square_Missing_Dist_XY(x,y,fabs(bx),mat);
  fw=fv=-fx;
  init_loglk = fw;
  max_loglk = UNLIKELY;
  bestx = bx;

  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(xx)+BRENT_ZEPS);

      if(fabs(xx-xm) <= (tol2-0.5*(b-a))) 
	{
	  mat->dist[x][y] = xx;
	  Least_Square_Missing_Dist_XY(x,y,mat->dist[x][y],mat);
	  return -fx;
	}
      
      if(fabs(e) > tol1) 
	{
	  r=(xx-w)*(fx-fv);
	  q=(xx-v)*(fx-fw);
	  p=(xx-v)*q-(xx-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=fabs(q);
	  etemp=e;
	  e=d;
	  if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-xx) || p >= q*(b-xx))
	    {
	      d=BRENT_CGOLD*(e=(xx >= xm ? a-xx : b-xx));
	      /*                   printf("Golden section step\n"); */
	    }
	  else
	    {
	      d=p/q;
	      u=xx+d;
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-xx);
	      /*                   printf("Parabolic step\n"); */
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(xx >= xm ? a-xx : b-xx));
	  /*               printf("Golden section step (default)\n"); */
	}
      
      u=(fabs(d) >= tol1 ? xx+d : xx+SIGN(tol1,d));
      fu = Least_Square_Missing_Dist_XY(x,y,fabs(u),mat);
            
/*       printf("param=%f loglk=%f\n",u,fu); */
      
      if(fu <= fx) 
	{
	  if(iter > n_iter_max) return -fu;

	  if(u >= xx) a=xx; else b=xx;
	  SHFT(v,w,xx,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < xx) a=u; else b=u;
	  if (fu <= fw || w == xx) 
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
	  else if (fu <= fv || v == xx || v == w) {
            v=u;
            fv=fu;
	  }
	}
    }
  Exit("\n. Too many iterations in BRENT !");
  return(-1);
}

/*********************************************************/

void Opt_Missing_Dist(int x, int y, matrix *mat)
{
  phydbl ax,bx,cx;

  ax = DIST_MAX;
  bx = DIST_MAX/4.;

  Missing_Dist_Brak(&ax,&bx,&cx,x,y,mat);
  printf("ax=%f bx=%f cx=%f\n",fabs(ax),fabs(bx),fabs(cx));
  Missing_Dist_Brent(fabs(ax),fabs(bx),fabs(cx),1.E-5,100,x,y,mat);
}

/*********************************************************/

int Optimiz_Alpha_And_Pinv(arbre *tree)
{

  int    iter;
  phydbl best_alpha, best_pinv;
  phydbl slope, intercept;
  phydbl lk_b, lk_a;
  phydbl lk_init, lk_final;
  phydbl f0,f1,f2,f3,x0,x1,x2,x3;
  phydbl pinv0, pinv1;
  phydbl a, b, c;
  phydbl fa, fb, fc;
  phydbl K;
  phydbl alpha0, alpha1;
  phydbl best_lnL;


  lk_final   = UNLIKELY;
  lk_b       = UNLIKELY;
  lk_a       = UNLIKELY;


/*   printf("\n\n. Init lnL = %f alpha=%f pinv=%f", */
/* 	 tree->c_lnL, */
/* 	 tree->mod->alpha, */
/* 	 tree->mod->pinvar); */
    
  /* Two (full) steps to compute  pinv_alpha_slope & pinv_alpha_intercept */

  tree->both_sides = 1;
  Lk(tree);
  lk_b = tree->c_lnL;

  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
			tree->noeud[0]->b[0],tree,tree->data);

  tree->both_sides = 0;
  Optimize_Single_Param_Generic(tree,&(tree->mod->pinvar),.0001,0.9999,
				tree->mod->s_opt->min_diff_lk_global,
				tree->mod->s_opt->brent_it_max);
  Optimize_Single_Param_Generic(tree,&(tree->mod->alpha),.01,100.,
				tree->mod->s_opt->min_diff_lk_global,
				tree->mod->s_opt->brent_it_max);
  
  pinv0  = tree->mod->pinvar;
  alpha0 = tree->mod->alpha;
  f0 = tree->c_lnL;


  tree->both_sides = 1;
  Lk(tree);

  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
			tree->noeud[0]->b[0],tree,tree->data);

  tree->both_sides = 0;
  Optimize_Single_Param_Generic(tree,&(tree->mod->pinvar),.0001,0.9999,
				tree->mod->s_opt->min_diff_lk_global,
				tree->mod->s_opt->brent_it_max);
  Optimize_Single_Param_Generic(tree,&(tree->mod->alpha),.01,100.,
				tree->mod->s_opt->min_diff_lk_global,
				tree->mod->s_opt->brent_it_max);

  lk_a = tree->c_lnL;

  pinv1  = tree->mod->pinvar;
  alpha1 = tree->mod->alpha;
  f1 = tree->c_lnL;
  best_lnL = f1;


  if(lk_a < lk_b - tree->mod->s_opt->min_diff_lk_local)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else if(fabs(lk_a - lk_b) < tree->mod->s_opt->min_diff_lk_local)
    {
      return 1;
    }
    
  Record_Br_Len(tree);
  best_alpha = tree->mod->alpha;
  best_pinv  = tree->mod->pinvar;
  lk_init    = tree->c_lnL;

/*   printf("\n\n. Init lnL after std opt = %f",tree->c_lnL); */

  slope     = (pinv1 - pinv0)/(alpha1 - alpha0);
  intercept = pinv1 - slope * alpha1;
  
  if((slope > 0.001) && (slope < 1./0.001))
    {
/*       printf("\n. pinv0 = %f, pinv1 = %f, alpha0 = %f, alpha1 = %f",pinv0,pinv1,alpha0,alpha1); */
/*       printf("\n. slope = %f intercept = %f",slope,intercept); */
      
      K = 0.381966;
      
      if(alpha1 < alpha0) 
	{
	  c  = alpha0;
	  b  = alpha1;
	  fc = f0;
	  fb = f1;
	  
	  a = (0.1 < alpha1)?(0.1):(0.5*alpha1);
	  tree->mod->alpha = a;
	  tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	  if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	  if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	  tree->both_sides = 1;
	  Lk(tree);
	  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				tree->noeud[0]->b[0],tree,tree->data);
	  fa = tree->c_lnL;
	  
	  iter = 0;	

/* 	  printf("\n. a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f (alpha=%f pinv=%f)",a,b,c,fa,fb,fc,tree->mod->alpha,tree->mod->pinvar); */

	  while(fa > fb)
	    {
	      a = a/5.;
	      tree->mod->alpha = a;
	      tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	      if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	      if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	      tree->both_sides = 1;
	      Lk(tree);
	      Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				    tree->noeud[0]->b[0],tree,tree->data);
	      fa = tree->c_lnL;
/* 	      printf("\n+ a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f",a,b,c,fa,fb,fc); */
	      if(iter++ > 10) return 0;
	    }
	}
      else
	{
	  a  = alpha0;
	  b  = alpha1;
	  fa = f0;
	  fb = f1;
	  
	  c = (alpha1 < 2.)?(2.0):(2.*alpha1);
	  tree->mod->alpha = c;
	  tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	  if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	  if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	  tree->both_sides = 1;
	  Lk(tree);
	  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				tree->noeud[0]->b[0],tree,tree->data);
	  fc = tree->c_lnL;

/* 	  printf("\n. a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f (alpha=%f pinv=%f)",a,b,c,fa,fb,fc,tree->mod->alpha,tree->mod->pinvar); */

	  iter = 0;
	  while(fc > fb)
	    {
	      c = c*2.;
	      tree->mod->alpha = c;
	      tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	      if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	      if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	      tree->both_sides = 1;
	      Lk(tree);
	      Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				    tree->noeud[0]->b[0],tree,tree->data);
	      fc = tree->c_lnL;
/* 	      printf("\n+ a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f",a,b,c,fa,fb,fc); */
	      if(iter++ > 10) return 0;
	    }
	}
      
      
      if(fabs(b - c) > fabs(a - b))
	{
	  x0 = a; x1 = b; x3 = c;
	  x2 = b + K * fabs(b - c);
	  
	  f0 = fa;
	  f1 = fb;
	  f3 = fc;
	  tree->mod->alpha = x2;
	  tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	  if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	  if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	  tree->both_sides = 1;
	  Lk(tree);
	  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				tree->noeud[0]->b[0],tree,tree->data);
	  f2 = tree->c_lnL;
	}
      else /* |b -c| < |a - b| */
	{
	  x0 = a; x2 = b; x3 = c;
	  x1 = b - K * fabs(b - a);
	  
	  f0 = fa;
	  f2 = fb;
	  f3 = fc;
	  tree->mod->alpha = x1;
	  tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	  if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	  if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	  tree->both_sides = 1;
	  Lk(tree);
	  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				tree->noeud[0]->b[0],tree,tree->data);
	  f1 = tree->c_lnL;
	}
      
      iter = 0;
      do
	{
/* 	  printf("\n. x0=%f, x1=%f, x2=%f, x3=%f, f0=%f, f1=%f, f2=%f, f3=%f", */
/* 		 x0,x1,x2,x3,f0,f1,f2,f3); */
	  
	  if(f1 > f2)
	    {
	      x3 = x2;
	      x2 = x1;
	      x1 = x2 - K * fabs(x2 - x0);
	      
	      f3 = f2;
	      f2 = f1;
	      
	      tree->mod->alpha = x1;
	      tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	      if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	      if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	      tree->both_sides = 1;
	      Lk(tree);
	      Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				    tree->noeud[0]->b[0],tree,tree->data);
	      f1 = tree->c_lnL;
	      if(f1 > best_lnL) 
		{
		  Record_Br_Len(tree);
		  best_alpha = tree->mod->alpha;
		  best_pinv = tree->mod->pinvar;
		}
/* 	      printf("\n> f1=%f",f1); */
	    }
	  else /* f1 < f2 */
	    {
	      x0 = x1;
	      x1 = x2;
	      x2 = x2 + K * fabs(x3 - x2);
	      
	      f0 = f1;
	      f1 = f2;
	      
	      tree->mod->alpha = x2;
	      tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	      if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	      if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	      tree->both_sides = 1;
	      Lk(tree);
	      Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				    tree->noeud[0]->b[0],tree,tree->data);
	      f2 = tree->c_lnL;
	      if(f2 > best_lnL) 
		{
		  Record_Br_Len(tree);
		  best_alpha = tree->mod->alpha;
		  best_pinv = tree->mod->pinvar;
		}
/* 	      printf("\n> f2=%f",f2); */
	    }
	  
	  if(fabs(f1 - f2) < 0.01) break;
	  
	  iter++;
	  
	}while(iter < 100);
    }
  
  tree->mod->alpha = best_alpha;
  tree->mod->pinvar = best_pinv;
  Restore_Br_Len(tree);      
  tree->both_sides = 1;
  Lk(tree);
/*   printf("\n\n. Init lnL after golden opt = %f",tree->c_lnL); */
  return 1;
}

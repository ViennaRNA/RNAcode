/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */

#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "options.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "mc.h"
#include "m4.h"
#include "draw.h"
#include "rates.h"




/*********************************************************/

void RATES_Monte_Carlo_Mean_Rates(arbre *tree)
{
  RATES_Monte_Carlo_Mean_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,1.0,tree);
  RATES_Monte_Carlo_Mean_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,1.0,tree);
}

/*********************************************************/

void RATES_Monte_Carlo_Mean_Rates_Pre(node *a, node *d, edge *b, phydbl curr_rate, arbre *tree)
{
  if(b)
    {
      phydbl curr_t, next_t, mean_rate;
      phydbl shape,exp;
      
      /**/
      shape = 1.0;
      /**/

      mean_rate = curr_rate;
      curr_t    = a->t - Rexp(tree->rates->lexp[b->num]);;
      

      while(curr_t > d->t)
	{
	  curr_rate = Rgamma(shape,tree->rates->br_r[b->num] / shape);
	  exp = Rexp(tree->rates->lexp[b->num]);
	  
/* 	  printf("\n. curr_rate=%f mean_r=%f exp=%f curr_t=%f lim=%f", */
/* 		 curr_rate,mean_rate,exp,curr_t,d->t); */

	  next_t = curr_t - exp;
	  
	  if(next_t > d->t)
	    {
	      mean_rate = (1./(a->t - next_t)) * (mean_rate * (a->t - curr_t) + curr_rate * (curr_t - next_t));
	    }
	  else
	    {
	      mean_rate = (1./(a->t - d->t)) * (mean_rate * (a->t - curr_t) + curr_rate * (curr_t - d->t));
	    }

	  curr_t = next_t;
	}

      tree->rates->mc_mr[d->num][tree->rates->curr_mc_run] = mean_rate;
/*       printf("\n. %3d %f %f",d->num,mean_rate,tree->rates->br_r[b->num]); */
    }

  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Monte_Carlo_Mean_Rates_Pre(d,d->v[i],d->b[i],curr_rate,tree);	      
	    }
	}
    }  
  return;
}

/*********************************************************/

void RATES_Print_Rates(arbre *tree)
{
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
}

/*********************************************************/

void RATES_Print_Rates_Pre(node *a, node *d, edge *b, arbre *tree)
{

  if(b) 
    {
      printf("\n. Edge %3d m_rate = %10f r_rate = %10f t_a = %15f t_d = %15f",
	     b->num,
	     tree->rates->mean_r,
	     tree->rates->br_r[b->num],
	     a->t,d->t);
    }
  if(d->tax) return;
  else
    {
      int i;

      For(i,3) 
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Print_Rates_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

/*********************************************************/

trate *RATES_Make_Rate_Struct(arbre *tree)
{
  trate *rates;

  rates = (trate *)mCalloc(1,sizeof(trate));
  rates->br_r = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
  rates->lexp = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
  rates->mc_mr = (phydbl **)mCalloc(2*tree->n_otu-1,sizeof(phydbl *));

  return rates;
}

/*********************************************************/

void RATES_Init_Rate_Struct(trate *rates, arbre *tree)
{
  int i;

  rates->n_mc_runs   = 500; 
  rates->mean_r      = 0.00001;
  rates->curr_mc_run = 0;

  For(i,2*tree->n_otu-3) 
    {
      rates->br_r[i] = 1.0;
      rates->lexp[i] = 0.001;
    }

  For(i,2*tree->n_otu-1)
    rates->mc_mr[i] = (phydbl *)mCalloc(rates->n_mc_runs,sizeof(phydbl));

}

/*********************************************************/

/*********************************************************/
/*********************************************************/

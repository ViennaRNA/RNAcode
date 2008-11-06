/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "m4.h"
#include "mc.h"

/* int    LIM_SCALE; */
/* phydbl LIM_SCALE_VAL; */
/* phydbl MDBL_MAX; */
/* phydbl MDBL_MIN; */

/*********************************************************/

void Init_Tips_At_One_Site_Nucleotides_Float(char state, phydbl *p_lk)
{
  switch(state)
    {
    case 'A' : p_lk[0]=1.; p_lk[1]=p_lk[2]=p_lk[3]=.0;
      break;
    case 'C' : p_lk[1]=1.; p_lk[0]=p_lk[2]=p_lk[3]=.0;
      break;
    case 'G' : p_lk[2]=1.; p_lk[1]=p_lk[0]=p_lk[3]=.0;
      break;
    case 'T' : p_lk[3]=1.; p_lk[1]=p_lk[2]=p_lk[0]=.0;
      break;
    case 'U' : p_lk[3]=1.; p_lk[1]=p_lk[2]=p_lk[0]=.0;
      break;
    case 'M' : p_lk[0]=p_lk[1]=1.; p_lk[2]=p_lk[3]=.0;
      break;
    case 'R' : p_lk[0]=p_lk[2]=1.; p_lk[1]=p_lk[3]=.0;
      break;
    case 'W' : p_lk[0]=p_lk[3]=1.; p_lk[1]=p_lk[2]=.0;
      break;
    case 'S' : p_lk[1]=p_lk[2]=1.; p_lk[0]=p_lk[3]=.0;
      break;
    case 'Y' : p_lk[1]=p_lk[3]=1.; p_lk[0]=p_lk[2]=.0;
      break;
    case 'K' : p_lk[2]=p_lk[3]=1.; p_lk[0]=p_lk[1]=.0;
      break;
    case 'B' : p_lk[1]=p_lk[2]=p_lk[3]=1.; p_lk[0]=.0;
      break;
    case 'D' : p_lk[0]=p_lk[2]=p_lk[3]=1.; p_lk[1]=.0;
      break;
    case 'H' : p_lk[0]=p_lk[1]=p_lk[3]=1.; p_lk[2]=.0;
      break;
    case 'V' : p_lk[0]=p_lk[1]=p_lk[2]=1.; p_lk[3]=.0;
      break;
    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
      p_lk[0]=p_lk[1]=p_lk[2]=p_lk[3]=1.;break;
    default :
      {
	printf("\n. Unknown character state : %c\n",state);
	Exit("\n. Init failed (check the data type)\n");
	break;
      }
    }
}

/*********************************************************/

void Init_Tips_At_One_Site_Nucleotides_Int(char state, short int *p_pars)
{
  switch(state)
    {
    case 'A' : p_pars[0]=1; p_pars[1]=p_pars[2]=p_pars[3]=0;
      break;
    case 'C' : p_pars[1]=1; p_pars[0]=p_pars[2]=p_pars[3]=0;
      break;
    case 'G' : p_pars[2]=1; p_pars[1]=p_pars[0]=p_pars[3]=0;
      break;
    case 'T' : p_pars[3]=1; p_pars[1]=p_pars[2]=p_pars[0]=0;
      break;
    case 'U' : p_pars[3]=1; p_pars[1]=p_pars[2]=p_pars[0]=0;
      break;
    case 'M' : p_pars[0]=p_pars[1]=1; p_pars[2]=p_pars[3]=0;
      break;
    case 'R' : p_pars[0]=p_pars[2]=1; p_pars[1]=p_pars[3]=0;
      break;
    case 'W' : p_pars[0]=p_pars[3]=1; p_pars[1]=p_pars[2]=0;
      break;
    case 'S' : p_pars[1]=p_pars[2]=1; p_pars[0]=p_pars[3]=0;
      break;
    case 'Y' : p_pars[1]=p_pars[3]=1; p_pars[0]=p_pars[2]=0;
      break;
    case 'K' : p_pars[2]=p_pars[3]=1; p_pars[0]=p_pars[1]=0;
      break;
    case 'B' : p_pars[1]=p_pars[2]=p_pars[3]=1; p_pars[0]=0;
      break;
    case 'D' : p_pars[0]=p_pars[2]=p_pars[3]=1; p_pars[1]=0;
      break;
    case 'H' : p_pars[0]=p_pars[1]=p_pars[3]=1; p_pars[2]=0;
      break;
    case 'V' : p_pars[0]=p_pars[1]=p_pars[2]=1; p_pars[3]=0;
      break;
    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
      p_pars[0]=p_pars[1]=p_pars[2]=p_pars[3]=1;break;
    default :
      {
	printf("\n. Unknown character state : %c\n",state);
	Exit("\n. Init failed (check the data type)\n");
	break;
      }
    }
}

/*********************************************************/

void Init_Tips_At_One_Site_AA_Float(char aa, phydbl *p_lk)
{
  int i;

  For(i,20) p_lk[i] = .0;

  switch(aa){
  case 'A' : p_lk[0]= 1.; break;/* Alanine */
  case 'R' : p_lk[1]= 1.; break;/* Arginine */
  case 'N' : p_lk[2]= 1.; break;/* Asparagine */
  case 'D' : p_lk[3]= 1.; break;/* Aspartic acid */
  case 'C' : p_lk[4]= 1.; break;/* Cysteine */
  case 'Q' : p_lk[5]= 1.; break;/* Glutamine */
  case 'E' : p_lk[6]= 1.; break;/* Glutamic acid */
  case 'G' : p_lk[7]= 1.; break;/* Glycine */
  case 'H' : p_lk[8]= 1.; break;/* Histidine */
  case 'I' : p_lk[9]= 1.; break;/* Isoleucine */
  case 'L' : p_lk[10]=1.; break;/* Leucine */
  case 'K' : p_lk[11]=1.; break;/* Lysine */
  case 'M' : p_lk[12]=1.; break;/* Methionine */
  case 'F' : p_lk[13]=1.; break;/* Phenylalanin */
  case 'P' : p_lk[14]=1.; break;/* Proline */
  case 'S' : p_lk[15]=1.; break;/* Serine */
  case 'T' : p_lk[16]=1.; break;/* Threonine */
  case 'W' : p_lk[17]=1.; break;/* Tryptophan */
  case 'Y' : p_lk[18]=1.; break;/* Tyrosine */
  case 'V' : p_lk[19]=1.; break;/* Valine */

  case 'B' : p_lk[2]= 1.; break;/* Asparagine */
  case 'Z' : p_lk[5]= 1.; break;/* Glutamine */

  case 'X' : case '?' : case '-' : For(i,20) p_lk[i] = 1.; break;
  default :
    {
      printf("\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (check the data type)\n");
      break;
    }
  }
}

/*********************************************************/

void Init_Tips_At_One_Site_AA_Int(char aa, short int *p_pars)
{
  int i;

  For(i,20) p_pars[i] = .0;

  switch(aa){
  case 'A' : p_pars[0]  = 1; break;/* Alanine */
  case 'R' : p_pars[1]  = 1; break;/* Arginine */
  case 'N' : p_pars[2]  = 1; break;/* Asparagine */
  case 'D' : p_pars[3]  = 1; break;/* Aspartic acid */
  case 'C' : p_pars[4]  = 1; break;/* Cysteine */
  case 'Q' : p_pars[5]  = 1; break;/* Glutamine */
  case 'E' : p_pars[6]  = 1; break;/* Glutamic acid */
  case 'G' : p_pars[7]  = 1; break;/* Glycine */
  case 'H' : p_pars[8]  = 1; break;/* Histidine */
  case 'I' : p_pars[9]  = 1; break;/* Isoleucine */
  case 'L' : p_pars[10] = 1; break;/* Leucine */
  case 'K' : p_pars[11] = 1; break;/* Lysine */
  case 'M' : p_pars[12] = 1; break;/* Methionine */
  case 'F' : p_pars[13] = 1; break;/* Phenylalanin */
  case 'P' : p_pars[14] = 1; break;/* Proline */
  case 'S' : p_pars[15] = 1; break;/* Serine */
  case 'T' : p_pars[16] = 1; break;/* Threonine */
  case 'W' : p_pars[17] = 1; break;/* Tryptophan */
  case 'Y' : p_pars[18] = 1; break;/* Tyrosine */
  case 'V' : p_pars[19] = 1; break;/* Valine */

  case 'B' : p_pars[2]  = 1; break;/* Asparagine */
  case 'Z' : p_pars[5]  = 1; break;/* Glutamine */

  case 'X' : case '?' : case '-' : For(i,20) p_pars[i] = 1; break;
  default :
    {
      printf("\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (check the data type)\n");
      break;
    }
  }
}

/*********************************************************/

void Get_All_Partial_Lk_Scale(arbre *tree, edge *b_fcus, node *a, node *d)
{
  if(d->tax) return;
  else
    {
      Update_P_Lk(tree,b_fcus,d);
    }
}

/*********************************************************/

void Post_Order_Lk(node *a, node *d, arbre *tree)
{
  int i,dir;

  dir = -1;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Post_Order_Lk(d,d->v[i],tree);
	  else dir = i;
	}      
      Get_All_Partial_Lk_Scale(tree,d->b[dir],a,d);
    }
}

/*********************************************************/

void Pre_Order_Lk(node *a, node *d, arbre *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Get_All_Partial_Lk_Scale(tree,d->b[i],d->v[i],d);
	      Pre_Order_Lk(d,d->v[i],tree);
	    }
	}
    }
}

/*********************************************************/

void Lk(arbre *tree)
{
  int br,site;
  int n_patterns;

  n_patterns = tree->n_pattern;

  tree->number_of_lk_calls++;

  Set_Model_Parameters(tree->mod);

#ifndef PHYML
  if(tree->bl_from_node_stamps) MC_Bl_From_T(tree);
#endif

  For(br,2*tree->n_otu-3)
    {
      if(!tree->t_edges[br]->rght->tax)
	For(site,n_patterns) tree->t_edges[br]->sum_scale_f_rght[site] = .0;

      if(!tree->t_edges[br]->left->tax)
	For(site,n_patterns) tree->t_edges[br]->sum_scale_f_left[site] = .0;
      
      Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
    }

  Post_Order_Lk(tree->noeud[0],tree->noeud[0]->v[0],tree);
  if(tree->both_sides)
    Pre_Order_Lk(tree->noeud[0],
		 tree->noeud[0]->v[0],
		 tree);

  tree->c_lnL     = .0;
  tree->curr_catg =  0;
  tree->curr_site =  0;
  For(site,n_patterns)
    {
      tree->c_lnL_sorted[site] = .0;
      tree->site_lk[site]      = .0;
      tree->curr_site          = site;
      Site_Lk(tree);
    }

/*   Qksort(tree->c_lnL_sorted,0,n_patterns-1); */

  tree->c_lnL = .0;
  For(site,n_patterns)
    {
      if(tree->c_lnL_sorted[site] < .0) /* WARNING : change cautiously */
	tree->c_lnL += tree->c_lnL_sorted[site];
    }
}

/*********************************************************/

void Site_Lk(arbre *tree)
{
  edge *eroot;

  eroot = tree->noeud[0]->b[0];

  if(!eroot->rght->tax)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(tree->data->wght[tree->curr_site] > MDBL_MIN) Lk_Core(eroot,tree);
  else tree->c_lnL_sorted[tree->curr_site] = 1.; /* WARNING : change cautiously */
}

/*********************************************************/

phydbl Lk_At_Given_Edge(edge *b_fcus, arbre *tree)
{
  int n_patterns;

  tree->number_of_branch_lk_calls++;

  n_patterns = tree->n_pattern;

#ifndef PHYML
  if(tree->bl_from_node_stamps) MC_Bl_From_T(tree);
#endif

  Update_PMat_At_Given_Edge(b_fcus,tree);

  if(b_fcus->left->tax)
    {
      printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  tree->c_lnL = .0;
  For(tree->curr_site,n_patterns)
    {
      if(tree->data->wght[tree->curr_site] > MDBL_MIN) Lk_Core(b_fcus,tree);
      else tree->c_lnL_sorted[tree->curr_site] = 1.; /* WARNING : change cautiously */
    }

/*   Qksort(tree->c_lnL_sorted,0,n_patterns-1); */

  tree->c_lnL = .0;
  For(tree->curr_site,n_patterns)
    if(tree->c_lnL_sorted[tree->curr_site] < .0) /* WARNING : change cautiously */
      tree->c_lnL += tree->c_lnL_sorted[tree->curr_site];

  return tree->c_lnL;
}

/*********************************************************/

phydbl Lk_Core(edge *b, arbre *tree)
{
  phydbl log_site_lk, site_lk, site_lk_cat;
  phydbl scale_left, scale_rght;
  phydbl sum;
  int ambiguity_check,state;
  int catg,ns,k,l,site;

  log_site_lk = site_lk = site_lk_cat = .0;
  ambiguity_check = state = -1;
  site = tree->curr_site;
  ns = tree->mod->ns;

  scale_left = 
    (b->sum_scale_f_left)?
    (b->sum_scale_f_left[site]):
    (0.0);
  
  scale_rght = 
    (b->sum_scale_f_rght)?
    (b->sum_scale_f_rght[site]):
    (0.0);

  if((b->rght->tax) && (!tree->mod->s_opt->greedy))
    {
      ambiguity_check = tree->data->c_seq[b->rght->num]->is_ambigu[site];      
      if(!ambiguity_check) state = Get_State_From_P_Pars(b->p_lk_tip_r[site],tree);
    }

  if(tree->mod->use_m4mod) ambiguity_check = 1;

  For(catg,tree->mod->n_catg)
    {
      site_lk_cat = .0;

      if((b->rght->tax) && (!tree->mod->s_opt->greedy))
	{
	  if(!ambiguity_check)
	    {
	      sum = .0;
	      For(l,ns)
		{
		  sum +=
		    (phydbl)(b->Pij_rr[catg][state][l]) *
		    b->p_lk_left[site][catg][l];
		}
	      site_lk_cat += sum * tree->mod->pi[state];
	    }
	  else
	    {
	      For(k,ns)
		{
		  sum = .0;
		  if(b->p_lk_tip_r[site][k] > .0)
		    {
		      For(l,ns)
			{
			  sum +=
			    (phydbl)(b->Pij_rr[catg][k][l]) *
			    b->p_lk_left[site][catg][l];
			}
		      site_lk_cat +=
			sum *
			tree->mod->pi[k] *
			(phydbl)(b->p_lk_tip_r[site][k]);
		    }
		}
	    }
	}
      else
	{
	  For(k,ns) 
	    {
	      sum = .0;
	      if(b->p_lk_rght[site][catg][k] > .0)
		{
		  For(l,ns)
		    {
		      sum +=
			(phydbl)(b->Pij_rr[catg][k][l]) *
			b->p_lk_left[site][catg][l];
		    }
		  site_lk_cat +=
		    sum *
		    tree->mod->pi[k] *
		    b->p_lk_rght[site][catg][k];
		}
	    }
	}

      tree->log_site_lk_cat[catg][site] = site_lk_cat;
      site_lk += site_lk_cat * tree->mod->gamma_r_proba[catg];

    }

  /* site_lk may be too small ? */
  if(site_lk < 1.E-300) site_lk = 1.E-300;

  if(!tree->mod->invar)
    {      
      log_site_lk = (phydbl)log(site_lk) + scale_left + scale_rght;
    }
  else
    {
      if((phydbl)tree->data->invar[site] > -0.5)
	{
	  if((scale_left + scale_rght > 0.0) || (scale_left + scale_rght < 0.0))
	    site_lk *= (phydbl)exp(scale_left + scale_rght);
	  
	  log_site_lk = (phydbl)log(site_lk*(1.0-tree->mod->pinvar) + tree->mod->pinvar*tree->mod->pi[tree->data->invar[site]]);
	}
      else
	{
	  log_site_lk = (phydbl)log(site_lk*(1.0-tree->mod->pinvar)) + scale_left + scale_rght;
	}
    }
  
  if(log_site_lk < -MDBL_MAX) Warn_And_Exit("\nlog_site_lk < -MDBL_MAX\n");

  For(catg,tree->mod->n_catg)
    tree->log_site_lk_cat[catg][site] = 
    (phydbl)log(tree->log_site_lk_cat[catg][site]) +
    scale_left + 
    scale_rght;
  
  tree->site_lk[site]      = log_site_lk;
  tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;
  return log_site_lk;
}

/*********************************************************/

phydbl Return_Lk(arbre *tree)
{
  Lk(tree);
  return tree->c_lnL;
}

/*********************************************************/

phydbl Return_Abs_Lk(arbre *tree)
{
  Lk(tree);
  return fabs(tree->c_lnL);
}

/*********************************************************/

matrix *ML_Dist(allseq *data, model *mod)
{
  int i,j,k,l;
  phydbl init;
  int n_catg;
  phydbl d_max,sum;
  matrix *mat;
  allseq *twodata,*tmpdata;
  int state0, state1,len;
  phydbl *F;
  eigen *eigen_struct;

  tmpdata             = (allseq *)mCalloc(1,sizeof(allseq));
  tmpdata->c_seq      = (seq **)mCalloc(2,sizeof(seq *));
  tmpdata->b_frq      = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  tmpdata->ambigu     = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  F                   = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl ));
  eigen_struct        = (eigen *)Make_Eigen_Struct(mod);

  tmpdata->n_otu      = 2;

  tmpdata->crunch_len = 0;
  tmpdata->init_len   = 0;
  For(i,data->crunch_len)
    {
      if(data->wght[i] > .0)
	{
	  tmpdata->crunch_len++;
	  tmpdata->init_len+=(int)data->wght[i];
	}
    }

  mat =
    (mod->datatype == NT) ?
    ((mod->whichmodel < 10)?(K80_dist(data,2000)):(JC69_Dist(data,mod))):
    (JC69_Dist(data,mod));

  For(i,mod->n_catg) /* Don't use gamma distribution */
    {
      mod->gamma_rr[i]      = 1.0;
      mod->gamma_r_proba[i] = 1.0;
    }

  n_catg = mod->n_catg;
  mod->n_catg = 1;

  For(j,data->n_otu-1)
    {
      tmpdata->c_seq[0]       = data->c_seq[j];
      tmpdata->c_seq[0]->name = data->c_seq[j]->name;
      tmpdata->wght           = data->wght;

      for(k=j+1;k<data->n_otu;k++)
	{
	  tmpdata->c_seq[1]       = data->c_seq[k];
	  tmpdata->c_seq[1]->name = data->c_seq[k]->name;

	  twodata = Compact_CSeq(tmpdata,mod);
	  For(l,mod->ns) twodata->b_frq[l] = data->b_frq[l];
	  Check_Ambiguities(twodata,mod->datatype,1);
	  Hide_Ambiguities(twodata);
	  
	  init = mat->dist[j][k];
	  if((init == DIST_MAX) || (init < .0)) init = 0.1;
	  	  
	  d_max = init;
	  
	  For(i,mod->ns*mod->ns) F[i]=.0;
	  len = 0;
	  For(l,twodata->c_seq[0]->len)
	    {
	      state0 = Assign_State(twodata->c_seq[0]->state+l,mod->datatype,mod->stepsize);
	      state1 = Assign_State(twodata->c_seq[1]->state+l,mod->datatype,mod->stepsize);
	      if((state0 > -1) && (state1 > -1))
		{
		  F[mod->ns*state0+state1] += twodata->wght[l];
		  len += (int)twodata->wght[l];
		}
	    }
	  if(len > .0) {For(i,mod->ns*mod->ns) F[i] /= (phydbl)len;}
	  	  
	  sum = 0.;
	  For(i,mod->ns*mod->ns) sum += F[i];
	  if(sum < .001) d_max = -1.;
	  else if((sum > 1. - .001) && (sum < 1. + .001)) Opt_Dist_F(&(d_max),F,mod);
	  else
	    {
	      printf("\n. sum = %f\n",sum);
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("");
	    }

	  /* BRENT */
	  /* d_max = Optimize_Dist(mod,init,twodata); */
	  
/* 	  printf("\n. Warning : not using the ML pairwise distances..."); */
/* 	  d_max = init; */
	  	  
	  if(d_max >= DIST_MAX)
	    {
/* 	      printf("\n. Large distance encountered between %s and %s sequences.", */
/* 		     tmpdata->c_seq[1]->name, */
/* 		     tmpdata->c_seq[0]->name); */
	      d_max = DIST_MAX;
	    }
	  
	  /* Do not correct for dist < BL_MIN, otherwise Fill_Missing_Dist 
	   *  will not be called 
	   */
	  
	  mat->dist[j][k] = d_max;
	  mat->dist[k][j] = mat->dist[j][k];
	  Free_Cseq(twodata);
	}
    }

  mod->n_catg = n_catg;
    
  Free(tmpdata->ambigu);
  Free(tmpdata->b_frq);
  Free(tmpdata->c_seq);
  free(tmpdata);
  Free_Eigen(eigen_struct);
  Free(F);

  return mat;
}

/*********************************************************/

phydbl Lk_Given_Two_Seq(allseq *data, int numseq1, int numseq2, phydbl dist, model *mod, phydbl *loglk)
{
  seq *seq1,*seq2;
  phydbl site_lk,log_site_lk;
  int i,j,k,l;
  phydbl **p_lk_l,**p_lk_r;
  phydbl len;

  DiscreteGamma(mod->gamma_r_proba, mod->gamma_rr, mod->alpha,
		mod->alpha,mod->n_catg,1);

  seq1 = data->c_seq[numseq1];
  seq2 = data->c_seq[numseq2];

  p_lk_l = (phydbl **)mCalloc(data->c_seq[0]->len,sizeof(phydbl *));
  p_lk_r = (phydbl **)mCalloc(data->c_seq[0]->len,sizeof(phydbl *));

  For(i,data->c_seq[0]->len)
    {
      p_lk_l[i] = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
      p_lk_r[i] = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
    }

  if(dist < BL_MIN) dist = BL_START;
  else if(dist > BL_MAX) dist = BL_START;

  For(i,mod->n_catg)
    {
      len = dist*mod->gamma_rr[i];
      if(len < BL_MIN) len = BL_MIN;
      else if(len > BL_MAX) len = BL_MAX;
      PMat(len,mod,&(mod->Pij_rr[i]));
    }

  if(mod->datatype == NT)
    {
      For(i,data->c_seq[0]->len)
	{
	  Init_Tips_At_One_Site_Nucleotides_Float(seq1->state[i],p_lk_l[i]);
	  Init_Tips_At_One_Site_Nucleotides_Float(seq2->state[i],p_lk_r[i]);
	}
    }
  else
    {
      For(i,data->c_seq[0]->len)
	{
	  Init_Tips_At_One_Site_AA_Float(seq1->state[i],p_lk_l[i]);
	  Init_Tips_At_One_Site_AA_Float(seq2->state[i],p_lk_r[i]);
	}
    }


  site_lk = .0;
  *loglk = 0;

  For(i,data->c_seq[0]->len)
    {
      if(data->wght[i])
	{
	  site_lk = log_site_lk = .0;
	  if(!data->ambigu[i])
	    {
	      For(k,mod->ns) {if(p_lk_l[i][k] > .0001) break;}
	      For(l,mod->ns) {if(p_lk_r[i][l] > .0001) break;}
	      For(j,mod->n_catg)
		{
		  site_lk +=
		    mod->gamma_r_proba[j] *
		    mod->pi[k] *
		    p_lk_l[i][k] *
		    (phydbl)mod->Pij_rr[j][k][l] *
		    p_lk_r[i][l];
		}
	    }
	  else
	    {
	      For(j,mod->n_catg)
		{
		  For(k,mod->ns) /*sort sum terms ? No global effect*/
		    {
		      For(l,mod->ns)
			{
			  site_lk +=
			    mod->gamma_r_proba[j] *
			    mod->pi[k] *
			    p_lk_l[i][k] *
			    (phydbl)mod->Pij_rr[j][k][l] *
			    p_lk_r[i][l];
			}
		    }
		}
	    }

/* 	  printf("'%c' '%c' -> %f\n",seq1->state[i],seq2->state[i],site_lk); */

	  if(site_lk <= .0)
	    {
	      printf("'%c' '%c'\n",seq1->state[i],seq2->state[i]);
	      Exit("\n. Err: site lk <= 0\n");
	    }

	  log_site_lk += (phydbl)log(site_lk);

	  *loglk += data->wght[i] * log_site_lk;/* sort sum terms ? No global effect*/
	}
    }

  For(i,data->c_seq[0]->len)
    {
      Free(p_lk_l[i]);
      Free(p_lk_r[i]);
    }
  Free(p_lk_l); Free(p_lk_r);
  return *loglk;
}

/*********************************************************/

void Unconstraint_Lk(arbre *tree)
{
  int i;

  tree->unconstraint_lk = .0;

  For(i,tree->data->crunch_len)
    {
      tree->unconstraint_lk +=
	tree->data->wght[i]*(phydbl)log(tree->data->wght[i]);
    }
  tree->unconstraint_lk -=
    tree->data->init_len*(phydbl)log(tree->data->init_len);
}

/*********************************************************/

void Update_P_Lk(arbre *tree, edge *b, node *d)
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
  node *n_v1, *n_v2;
  phydbl p1_lk1,p2_lk2;
  phydbl ***p_lk,***p_lk_v1,***p_lk_v2;
  double ***Pij1,***Pij2;
  phydbl max_p_lk;
  phydbl *sum_scale, *sum_scale_v1, *sum_scale_v2;
  phydbl scale_v1, scale_v2;
  int i,j;
  int catg,site;
  int dir1,dir2;
  int n_patterns;
  int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  
  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = -1;
  scale_v1 = scale_v2 = 0.0;
  p1_lk1 = p2_lk2 = .0;


  if(d->tax)
    {
      printf("\n. node %d is a leaf...",d->num);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

/*   n_patterns = (int)floor(tree->n_pattern*tree->prop_of_sites_to_consider); */
  n_patterns = tree->n_pattern;
  
  dir1=dir2=-1;
  For(i,3) if(d->b[i] != b) (dir1<0)?(dir1=i):(dir2=i);
  
  n_v1 = d->v[dir1];
  n_v2 = d->v[dir2];

  if(d == b->left)
    {
      p_lk = b->p_lk_left;
      sum_scale = b->sum_scale_f_left;
    }
  else
    {
      p_lk = b->p_lk_rght;
      sum_scale = b->sum_scale_f_rght;
    }
      
  if(d == d->b[dir1]->left)
    {
      p_lk_v1 = d->b[dir1]->p_lk_rght;
      sum_scale_v1 = d->b[dir1]->sum_scale_f_rght;
    }
  else
    {
      p_lk_v1 = d->b[dir1]->p_lk_left;
      sum_scale_v1 = d->b[dir1]->sum_scale_f_left;
    }
  
  if(d == d->b[dir2]->left)
    {
      p_lk_v2 = d->b[dir2]->p_lk_rght;
      sum_scale_v2 = d->b[dir2]->sum_scale_f_rght;
    }
  else
    {
      p_lk_v2 = d->b[dir2]->p_lk_left;
      sum_scale_v2 = d->b[dir2]->sum_scale_f_left;
    }
  
  Pij1 = d->b[dir1]->Pij_rr;
  Pij2 = d->b[dir2]->Pij_rr;
      
  For(site,n_patterns)
    {
      scale_v1 = (sum_scale_v1)?(sum_scale_v1[site]):(0.0);
      scale_v2 = (sum_scale_v2)?(sum_scale_v2[site]):(0.0); 
      sum_scale[site] = scale_v1 + scale_v2;

      max_p_lk = -MDBL_MAX;
      state_v1 = state_v2 = -1;
      ambiguity_check_v1 = ambiguity_check_v2 = -1;
      
      if(!tree->mod->s_opt->greedy)
	{
	  if(n_v1->tax)
	    {
	      ambiguity_check_v1 = tree->data->c_seq[n_v1->num]->is_ambigu[site];
	      if(!ambiguity_check_v1) state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r[site],tree);
	    }
	      
	  if(n_v2->tax)
	    {
	      ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site];
	      if(!ambiguity_check_v2) state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r[site],tree);
	    }
	}
      
      if(tree->mod->use_m4mod)
	{
	  ambiguity_check_v1 = 1;
	  ambiguity_check_v2 = 1;
	}

      For(catg,tree->mod->n_catg)
	{
	  For(i,tree->mod->ns)
	    {
	      p1_lk1 = .0;
	      
	      if((n_v1->tax) && (!tree->mod->s_opt->greedy))
		{
		  if(!ambiguity_check_v1)
		    {
		      p1_lk1 = (phydbl)Pij1[catg][i][state_v1];
		    }
		  else
		    {
		      For(j,tree->mod->ns)
			{
			  p1_lk1 += (phydbl)(Pij1[catg][i][j]) * (phydbl)(n_v1->b[0]->p_lk_tip_r[site][j]);
			}
		    }
		}
	      else
		{
		  For(j,tree->mod->ns)
		    {
		      p1_lk1 += (phydbl)(Pij1[catg][i][j]) * (phydbl)(p_lk_v1[site][catg][j]);
		    }
		}
	      
	      p2_lk2 = .0;
	      
	      if((n_v2->tax) && (!tree->mod->s_opt->greedy))
		{
		  if(!ambiguity_check_v2)
		    {
		      p2_lk2 = (phydbl)Pij2[catg][i][state_v2];
		    }
		  else
		    {
		      For(j,tree->mod->ns)
			{
			  p2_lk2 += (phydbl)(Pij2[catg][i][j]) * (phydbl)(n_v2->b[0]->p_lk_tip_r[site][j]);
			}
		    }
		}
	      else
		{
		  For(j,tree->mod->ns)
		    {
		      p2_lk2 += (phydbl)(Pij2[catg][i][j]) * (phydbl)(p_lk_v2[site][catg][j]);
		    }
		}
	      
	      p_lk[site][catg][i] = p1_lk1 * p2_lk2;
	      

	      if(p_lk[site][catg][i] > max_p_lk) max_p_lk = p_lk[site][catg][i];
	    }
	}
      
      if((max_p_lk < LIM_SCALE_VAL) || (max_p_lk > (1./LIM_SCALE_VAL)))
	{
	  For(catg,tree->mod->n_catg)
	    {
	      For(i,tree->mod->ns)
		{
		  p_lk[site][catg][i] /= max_p_lk;
		  
/* 		  if((p_lk[site][catg][i] > MDBL_MAX) || (p_lk[site][catg][i] < MDBL_MIN)) */
/* 		    { */
/* 		      printf("\n. Err in file %s at line %d",__FILE__,__LINE__); */
/* 		      printf("\n. p_lk[%3d][%2d][%3d] = %G max_p_lk = %G",site,catg,i,p_lk[site][catg][i],max_p_lk); */
/* 		      printf("\n. alpha=%f pinv=%f",tree->mod->alpha,tree->mod->pinvar); */
/* 		      For(i,tree->mod->n_catg) printf("\n. rr[%2d] = %G",i,tree->mod->rr[i]); */
/* 		      printf("\n. d->b[dir1]->l = %f, d->b[dir2]->l = %f",d->b[dir1]->l,d->b[dir2]->l); */
/* 		      printf("\n. d->v[dir1]->num = %d, d->v[dir2]->num = %d",d->v[dir1]->num,d->v[dir2]->num); */
/* 		      if(d->v[dir1]->tax) */
/* 			{ */
/* 			  printf("\n. Character observed at d->v[dir1] = %d",state_v1); */
/* 			} */
/* 		      if(d->v[dir2]->tax) */
/* 			{ */
/* 			  printf("\n. Character observed at d->v[dir2] = %d",state_v2); */
/* 			} */
/* 		      Warn_And_Exit("\n. Numerical precision problem ! (send me an e-mail : s.guindon@auckland.ac.nz)\n"); */
/* 		    } */
		}
	    }
	  sum_scale[site] += (phydbl)log(max_p_lk);
	}
    } 
}

/*********************************************************/


void Make_Tree_4_Lk(arbre *tree, allseq *alldata, int n_site)
{
  int i;

  tree->c_lnL_sorted = (phydbl *)mCalloc(tree->n_pattern, sizeof(phydbl));
  tree->site_lk      = (phydbl *)mCalloc(alldata->crunch_len,sizeof(phydbl));

  tree->log_site_lk_cat      = (phydbl **)mCalloc(tree->mod->n_catg,sizeof(phydbl *));
  For(i,tree->mod->n_catg)
    tree->log_site_lk_cat[i] = (phydbl *)mCalloc(alldata->crunch_len,sizeof(phydbl));

  tree->log_lks_aLRT = (phydbl **)mCalloc(3,sizeof(phydbl *));
  For(i,3) tree->log_lks_aLRT[i] = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));

  For(i,2*tree->n_otu-3)
    {
      Make_Edge_Lk(tree->t_edges[i],tree);
      Make_Edge_NNI(tree->t_edges[i]);
    }

  For(i,2*tree->n_otu-2) Make_Node_Lk(tree->noeud[i]);

  if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
  else Init_P_Lk_Tips_Int(tree);
}

/*********************************************************/

void Init_P_Lk_Tips_Double(arbre *tree)
{
  int curr_site,i,j,k;

  Fors(curr_site,tree->data->crunch_len,tree->mod->stepsize)
    {
      For(i,tree->n_otu)
	{
	  if (tree->mod->datatype == NT)
	    {
	      if(tree->noeud[i]->b[0]->rght->tax != 1)
		{
		  printf("\n. tree->noeud[i]->b[0]->rght->num = %d ; %f %f %f %f\n",
			 tree->noeud[i]->b[0]->rght->num,
			 tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][0],
			 tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][1],
			 tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][2],
			 tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][3]);
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}

	      Init_Tips_At_One_Site_Nucleotides_Float(tree->data->c_seq[i]->state[curr_site],
						      tree->noeud[i]->b[0]->p_lk_rght[curr_site][0]);
	    }
	  else
	    {
	      Init_Tips_At_One_Site_AA_Float(tree->data->c_seq[i]->state[curr_site],
					     tree->noeud[i]->b[0]->p_lk_rght[curr_site][0]);
	    }

	  for(j=1;j<tree->mod->n_catg;j++)
	    {
	      For(k,tree->mod->ns)
		{
		  tree->noeud[i]->b[0]->p_lk_rght[curr_site][j][k]=
		    tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][k];
		}
	    }
	}
    }

  #ifndef PHYML
  if(tree->mod->m4mod) M4_Init_P_Lk_Tips_Double(tree);
  #endif
}

/*********************************************************/

void Init_P_Lk_Tips_Int(arbre *tree)
{
  int curr_site,i;


  Fors(curr_site,tree->data->crunch_len,tree->mod->stepsize)
    {
      For(i,tree->n_otu)
	{
	  if(tree->mod->datatype == NT)
	    {
	      if(tree->noeud[i]->b[0]->rght->tax != 1)
		{
		  printf("\n. tree->noeud[i]->b[0]->rght->num = %d ; %f %f %f %f\n",
			 tree->noeud[i]->b[0]->rght->num,
			 tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][0],
			 tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][1],
			 tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][2],
			 tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][3]);
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}

	      Init_Tips_At_One_Site_Nucleotides_Int(tree->data->c_seq[i]->state[curr_site],
						    tree->noeud[i]->b[0]->p_lk_tip_r[curr_site]);
	    }
	  else
	    {
	      Init_Tips_At_One_Site_AA_Int(tree->data->c_seq[i]->state[curr_site],
					   tree->noeud[i]->b[0]->p_lk_tip_r[curr_site]);
	    }
	}
    }
 
 #ifndef PHYML
  if(tree->mod->m4mod) M4_Init_P_Lk_Tips_Int(tree);
  #endif

}

/*********************************************************/

void Init_P_Lk_At_One_Node(node *a, arbre *tree)
{
  int curr_site,j,k;

  if(!a->tax) Exit("\n. Node 'a' must be a tip (err. in Init_P_Lk_At_One_Node)\n");
  if(a->b[0]->rght != a) Exit("\n. Node 'a' must be on the right handside of the terminal branch (err. in Init_P_Lk_At_One_Node)\n");

  Fors(curr_site,tree->data->crunch_len,tree->mod->stepsize)
    {
      if(tree->mod->datatype == NT)
	Init_Tips_At_One_Site_Nucleotides_Float(tree->data->c_seq[a->num]->state[curr_site],
						tree->noeud[a->num]->b[0]->p_lk_rght[curr_site][0]);
      else
	Init_Tips_At_One_Site_AA_Float(tree->data->c_seq[a->num]->state[curr_site],
				       tree->noeud[a->num]->b[0]->p_lk_rght[curr_site][0]);

      for(j=1;j<tree->mod->n_catg;j++)
	{
	  For(k,tree->mod->ns)
	    {
	      a->b[0]->p_lk_rght[curr_site][j][k]=
		a->b[0]->p_lk_rght[curr_site][0][k];
	    }
	}
    }
}

/*********************************************************/

void Update_PMat_At_Given_Edge(edge *b_fcus, arbre *tree)
{
  int i;
  phydbl len;

  len = -1.0;

  if(b_fcus->l < BL_MIN)      b_fcus->l = BL_MIN;
  else if(b_fcus->l > BL_MAX) b_fcus->l = BL_MAX;

  For(i,tree->mod->n_catg)
    {
      if(b_fcus->has_zero_br_len) len = -1.0;
      else
	{
	  len = b_fcus->l*tree->mod->gamma_rr[i];
	  if(len < BL_MIN)      len = BL_MIN;
	  else if(len > BL_MAX) len = BL_MAX;
	}
      PMat(len,tree->mod,&b_fcus->Pij_rr[i]);
    }
}

/*********************************************************/

/* void Update_P_Lk_On_A_Path(node *a, node *d, edge *b, node *target_one_side, node *target_other_side, arbre *tree) */
/* { */


/*   /\* */
/*                 \               / */
/* 	   target\___________b_/ */
/* 		 /  \	\  \   \ */
/* 		/    \	 \  \	\ */

/*     target is the root of the subtree at which we want */
/*     the likelihood to be updated */
/*   *\/ */



/* /\*   printf("Update_p_lk on (%d %d) at %d (target=%d %d)\n", *\/ */
/* /\* 	 b->left->num, *\/ */
/* /\* 	 b->rght->num, *\/ */
/* /\* 	 a->num, *\/ */
/* /\* 	 target_one_side->num, *\/ */
/* /\* 	 target_other_side->num); *\/ */

/*   Update_P_Lk(tree,b,a); */
/*   if((a == target_one_side) && (d == target_other_side))  */
/*     return; */
/*   else */
/*     { */
/*       Update_P_Lk_On_A_Path(d, */
/* 			    d->v[tree->t_dir[d->num][target_other_side->num]], */
/* 			    d->b[tree->t_dir[d->num][target_other_side->num]], */
/* 			    target_one_side, */
/* 			    target_other_side, */
/* 			    tree);  */
/*     } */
/* } */

void Update_P_Lk_Along_A_Path(node **path, int path_length, arbre *tree)
{
  int i,j;

  For(i,path_length-1)
    {
      For(j,3)
	if(path[i]->v[j] == path[i+1])
	  {
	    if(path[i] == path[i]->b[j]->left)
	      {
		if(!path[i]->b[j]->is_p_lk_l_u2d) 
		  {
		    Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->left);		    
		  }
		path[i]->b[j]->is_p_lk_l_u2d = 1;
	      }

	    else if(path[i] == path[i]->b[j]->rght)
	      {
		if(!path[i]->b[j]->is_p_lk_r_u2d) 
		  {
		    Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->rght);
		  }
		path[i]->b[j]->is_p_lk_r_u2d = 1;
	      }
	    break;
	  }

#ifdef DEBUG
      if(j == 3)
	{
	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");
	}
#endif

    }
}

/*********************************************************/

void Sort_Sites_Based_On_Lk(arbre *tree)
{
  int i,j,k,l,m;
  phydbl buff_dbl;
  char buff_char;
  int buff_int;

  For(i,tree->data->crunch_len)
    {
      for(j=i+1;j<tree->data->crunch_len;j++)
	{
	  if(tree->site_lk[j]*tree->data->wght[j] < tree->site_lk[i]*tree->data->wght[i])
	    {
	      buff_dbl              = tree->site_lk[j];
	      tree->site_lk[j]      = tree->site_lk[i];
	      tree->site_lk[i]      = buff_dbl;

	      buff_dbl              = tree->data->wght[j];
	      tree->data->wght[j]   = tree->data->wght[i];
	      tree->data->wght[i]   = buff_dbl;

	      buff_int              = tree->data->invar[j];
	      tree->data->invar[j]  = tree->data->invar[i];
	      tree->data->invar[i]  = buff_int;

	      buff_int              = tree->data->ambigu[j];
	      tree->data->ambigu[j] = tree->data->ambigu[i];
	      tree->data->ambigu[i] = buff_int;

	      if(tree->data->c_seq[0]->state)
		{
		  For(k,tree->n_otu)
		    {
		      buff_char                      = tree->data->c_seq[k]->state[j];
		      tree->data->c_seq[k]->state[j] = tree->data->c_seq[k]->state[i];
		      tree->data->c_seq[k]->state[i] = buff_char;
		    }
		}

	      For(k,2*tree->n_otu-3)
		{
		  For(l,tree->mod->n_catg)
		    {
		      For(m,tree->mod->ns)
			{
			  buff_dbl                             = tree->t_edges[k]->p_lk_rght[j][l][m];
			  tree->t_edges[k]->p_lk_rght[j][l][m] = tree->t_edges[k]->p_lk_rght[i][l][m];
			  tree->t_edges[k]->p_lk_rght[i][l][m] = buff_dbl;

			  buff_dbl                             = tree->t_edges[k]->p_lk_left[j][l][m];
			  tree->t_edges[k]->p_lk_left[j][l][m] = tree->t_edges[k]->p_lk_left[i][l][m];
			  tree->t_edges[k]->p_lk_left[i][l][m] = buff_dbl;
			}
		    }

		  buff_dbl                              = tree->t_edges[k]->sum_scale_f_rght[j];
		  tree->t_edges[k]->sum_scale_f_rght[j] = tree->t_edges[k]->sum_scale_f_rght[i];
		  tree->t_edges[k]->sum_scale_f_rght[i] = buff_dbl;

		  buff_dbl                              = tree->t_edges[k]->sum_scale_f_left[j];
		  tree->t_edges[k]->sum_scale_f_left[j] = tree->t_edges[k]->sum_scale_f_left[i];
		  tree->t_edges[k]->sum_scale_f_left[i] = buff_dbl;
		}
	    }
	}
    }
}

/*********************************************************/

phydbl Lk_Dist(phydbl *F, phydbl dist, model *mod)
{
  int i,j,k;
  phydbl len,lnL,tmp;

  len = -1.;
  For(k,mod->n_catg)
    {
      len = dist*mod->gamma_rr[k];
      if(len < BL_MIN) len = BL_MIN;
      else if(len > BL_MAX) len = BL_MAX;
      PMat(len,mod,&(mod->Pij_rr[k]));
    }

  lnL = .0;
  For(i,mod->ns)
    {
      For(j,mod->ns)
	{
	  tmp = .0;
	  For(k,mod->n_catg)
	    {
	      tmp +=
		mod->gamma_r_proba[k] *
		mod->pi[i] *
		(phydbl)(mod->Pij_rr[k][i][j]);
	    }
	  lnL += F[mod->ns*i+j] * (phydbl)log(tmp);
	}
    }
  return lnL;
}

/*********************************************************/

phydbl Update_Lk_At_Given_Edge(edge *b_fcus, arbre *tree)
{
  if(!b_fcus->left->tax) Update_P_Lk(tree,b_fcus,b_fcus->left);
  if(!b_fcus->rght->tax) Update_P_Lk(tree,b_fcus,b_fcus->rght);
  tree->c_lnL = Lk_At_Given_Edge(b_fcus,tree);
  return tree->c_lnL;
}

/*********************************************************/

/*       root
           \
           /
          a
	  |
	  |
	  d
	 / \
        /   \
       w     x    	

       d->t has changed and we need to compute
       the likelihood.
       (1) update the three branch lengths l(ad), l(dw) and l(dx)
       (2) update the change proba matrices along these branches
       (3) update the likelihood of subtree (w,x) (WARNING: (a,x) and (a,w) are not updated)
*/
phydbl Lk_Triplet(node *a, node *d, arbre *tree)
{
  int i;
  phydbl max_height;
  phydbl up_bound, low_bound;

  if(d->tax)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  up_bound = low_bound = -1.0;
  max_height = -1.0;
  For(i,3)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  if(d->v[i]->t > max_height)
	    {
	      max_height = d->v[i]->t;
	    }
	}
      else
	{
	  up_bound = (a == tree->n_root)?(a->t):(d->v[i]->t);
	}
    }

  low_bound = max_height;

  if(up_bound < low_bound - 1.E-10)
    {
      printf("\n. a->num=%d d->num=%d",a->num,d->num);
      printf("\n. up_bound = %f, low_bound = %f",up_bound,low_bound);
      Warn_And_Exit("\n");
    }

  if(d->t < low_bound)     d->t = low_bound;
  else if(d->t > up_bound) d->t = up_bound;

  /* Step (1) */
  For(i,3)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  d->b[i]->l = 
	    (d->t - d->v[i]->t) * 
	    tree->rates->mean_r * 
	    tree->rates->br_r[d->b[i]->num];
	}
      else
	{
	  if(a == tree->n_root)
	    {
	      d->b[i]->l = 
		(tree->n_root->t - tree->n_root->v[0]->t + 
		 tree->n_root->t - tree->n_root->v[1]->t) * tree->rates->mean_r;
	    }
	  else
	    {
	      d->b[i]->l = (a->t - d->t) * tree->rates->mean_r * tree->rates->br_r[d->b[i]->num];	    
	    }
	}
    }
  
  /* Step (2) */
  For(i,3) Update_PMat_At_Given_Edge(d->b[i],tree);
  
  For(i,3) 
    if((d->v[i] == a) || (d->b[i] == tree->e_root))
      {
	Update_P_Lk(tree,d->b[i],d); 
	Lk_At_Given_Edge(d->b[i],tree);
	break;
      }

  return tree->c_lnL;
}

/*********************************************************/

void Print_Lk_Given_Edge_Recurr(node *a, node *d, edge *b, arbre *tree)
{
  printf("\n___ Edge %3d (left=%3d rght=%3d) lnL=%f",
	 b->num,
	 b->left->num,
	 b->rght->num,
	 Lk_At_Given_Edge(b,tree));

  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	if(d->v[i] != a)
	  Print_Lk_Given_Edge_Recurr(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/

/* Returns a vector containing the posterior probabilities of 
   the different branch rate classes 
*/
phydbl *Post_Prob_Rates_At_Given_Edge(edge *b, phydbl *post_prob, arbre *tree)
{
  phydbl norm_factor;
  int rcat, scale_int;
  phydbl sum,log2,lnL,worst_lnL,best_lnL,mid_lnL;


  For(rcat,tree->mod->n_rr_branch) post_prob[rcat] = .0;

  log2 = 0.6931472;
  
  best_lnL  = UNLIKELY;
  worst_lnL = .0;
  For(rcat,tree->mod->n_rr_branch)
    {
      tree->rates->br_r[b->num] = tree->mod->rr_branch[rcat];
      lnL = Lk_At_Given_Edge(b,tree);

      if(lnL < worst_lnL) worst_lnL = lnL;
      if(lnL > best_lnL)  best_lnL  = lnL;
      post_prob[rcat] = lnL;

      tree->rates->br_r[b->num] = 1.0;
    }

  mid_lnL = worst_lnL + fabs(worst_lnL - best_lnL)/2.;

  /* exp(log(P(D|r))) is hard to compute. Try exp(log(K*P(D|r)))
     instead. The value of K is chosen such that the most accurate
     estimates of the posterior probabilities are obtained for the
     most likely rate classes.
  */
  scale_int = 0;
  do scale_int++;
  while(best_lnL + scale_int * log2 < 20.0);

  norm_factor = .0;
  For(rcat,tree->mod->n_rr_branch)
    {
/*       printf("\n. best_lnL=%f curr_lnL=%f %f %E", */
/*  	     best_lnL, */
/* 	     post_prob[rcat] , */
/* 	     post_prob[rcat] + scale_int * log2, */
/* 	     exp(post_prob[rcat] + scale_int * log2)); */

      post_prob[rcat] = exp(post_prob[rcat] + scale_int * log2);


      post_prob[rcat] *= tree->mod->p_rr_branch[rcat];
      norm_factor += post_prob[rcat];
    }

  sum = .0;
  For(rcat,tree->mod->n_rr_branch)
    {
      post_prob[rcat] /= norm_factor;
/*       printf("%f ",post_prob[rcat]); */
      sum += post_prob[rcat];
    }
  
  if(sum < 0.999 || sum > 1.001)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  return post_prob;  
}

/*********************************************************/

phydbl Lk_With_MAP_Branch_Rates(arbre *tree)
{
  int br,rcat,best_post_prob_cat;
  phydbl *post_prob;
  phydbl best_post_prob;
  edge *b;


  post_prob = (phydbl *)mCalloc(tree->mod->n_rr_branch,sizeof(phydbl));
 
  Lk(tree);
  Record_Br_Len(tree);

  For(br,2*tree->n_otu-3)
    {
      b = tree->t_edges[br];

      /* Compute the posterior probability of each rate class on edge b */
      post_prob = (phydbl *)Post_Prob_Rates_At_Given_Edge(b,post_prob,tree);

      /* Find the most probable class */
      best_post_prob = UNLIKELY;
      best_post_prob_cat = -1;
      For(rcat,tree->mod->n_rr_branch)
	{
	  if(post_prob[rcat] > best_post_prob)
	    {
	      best_post_prob = post_prob[rcat];
	      best_post_prob_cat = rcat;
	    }
	}

      /* The relative rate on this branch corresponds to the most probable rate class */
      tree->rates->br_r[br] = tree->mod->rr_branch[best_post_prob_cat];
    }

  Lk(tree);
  
  For(br,2*tree->n_otu-3) tree->rates->br_r[br] = 1.0;

  Free(post_prob);

  return tree->c_lnL;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/

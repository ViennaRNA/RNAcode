/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "utilities.h"
#include "draw.h"


/*********************************************************/

void DR_Get_Tree_Coord(arbre *tree)
{
  DR_Init_Tdraw_Struct(tree->ps_tree);
  DR_Get_Tree_Box_Width(tree->ps_tree,tree);
  if(!tree->n_root) 
    {
      printf("\n. Adding root before rendering the tree.");
      Add_Root(tree->t_edges[0],tree);
    }
  else Update_Root_Pos(tree);
  DR_Dist_To_Root(tree->n_root,tree);
  tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree);
  DR_Get_X_Coord(tree->ps_tree,tree);
  DR_Get_Y_Coord(tree->ps_tree,tree);
}

/*********************************************************/

void DR_Print_Postscript_Header(int n_pages, FILE *fp)
{
  if(!fp)
    {
      printf("\n. Failed to open the postscript file.");
      printf("\n. Did you forget the '--ps' option ?.");
      Warn_And_Exit("\n");
    }

  fprintf(fp,"%%!PS-Adobe-3.0\n");
  fprintf(fp,"%%%%DocumentFonts: Times-Roman Times-Roman\n");
  fprintf(fp,"%%%%Creator: Stephane Guindon\n");
  fprintf(fp,"%%%%Title: tree\n");
  fprintf(fp,"%%%%BeginFeature: *PageSize\n"); 
  fprintf(fp,"a4\n");
  fprintf(fp,"%%%%EndFeature\n");
  fprintf(fp,"%%%%EndComments\n");
  fprintf(fp,"%%%%Pages: %d\n",n_pages);

  fprintf(fp,"/lt {lineto} bind def\n");
  fprintf(fp,"/mt {moveto} bind def\n");
  fprintf(fp,"/sc {setrgbcolor} bind def\n");

  fprintf(fp,"/clipbox\n");
  fprintf(fp,"{\n");
  fprintf(fp,"newpath\n");
  fprintf(fp,"40 40 moveto\n");
  fprintf(fp,"560 40 lineto\n");
  fprintf(fp,"560 820 lineto\n");
  fprintf(fp,"40 820 lineto\n");
  fprintf(fp,"40 40 lineto\n");
  fprintf(fp,"closepath\n");
  fprintf(fp,"clip\n");
  fprintf(fp,"} bind def\n");
  
  fprintf(fp,"/Times-Roman findfont\n");
  fprintf(fp,"12 scalefont\n");
  fprintf(fp,"setfont\n");


}

/*********************************************************/

void DR_Print_Postscript_EOF(FILE *fp)
{
  fprintf(fp,"%%%%Trailer\n");
  fprintf(fp,"%%%%EOF\n");
}

/*********************************************************/

void DR_Print_Tree_Postscript(int page_num, FILE *fp, arbre *tree)
{
  int i;
  tdraw *draw;
  node *n_root;
  
  draw = tree->ps_tree;
  DR_Get_Tree_Coord(tree);
  n_root = tree->n_root;

  fprintf(fp,"%%%%Page: %d %d\n",page_num,page_num); 
  fprintf(fp,"clipbox\n");
  fprintf(fp,"stroke\n");
  fprintf(fp,"50 50 translate\n");
  fprintf(fp,"newpath\n");

/*   if(b_root->prob_sel_regime <= 0.1) */
/*     fprintf(fp,".0 .0 1. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.1 && b_root->prob_sel_regime <= 0.2) */
/*     fprintf(fp,".0 .5 1. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.2 && b_root->prob_sel_regime <= 0.3) */
/*     fprintf(fp,".0 1. 1. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.3 && b_root->prob_sel_regime <= 0.4) */
/*     fprintf(fp,".0 1. .5 sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.4 && b_root->prob_sel_regime <= 0.5) */
/*     fprintf(fp,".0 1. .0 sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.5 && b_root->prob_sel_regime <= 0.6) */
/*     fprintf(fp,".5 1. .0 sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.6 && b_root->prob_sel_regime <= 0.7) */
/*     fprintf(fp,"1. 1. 0. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.7 && b_root->prob_sel_regime <= 0.8) */
/*     fprintf(fp,"1. .5 0. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.8 && b_root->prob_sel_regime <= 0.9) */
/*     fprintf(fp,"1. 0. 0. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.9) */
/*     fprintf(fp,"1. .0 .0 sc\n"); */

  fprintf(fp,"%d %d mt\n",draw->xcoord[n_root->v[0]->num],draw->ycoord[n_root->v[0]->num]);
  fprintf(fp,"%d %d lt\n",0,draw->ycoord[n_root->v[0]->num]);
  fprintf(fp,"%d %d lt\n",0,draw->ycoord[n_root->v[1]->num]);
  fprintf(fp,"%d %d lt\n",draw->xcoord[n_root->v[1]->num],draw->ycoord[n_root->v[1]->num]);
  fprintf(fp,"stroke\n");


  fprintf(fp,"%d %d mt\n",draw->xcoord[n_root->v[0]->num],draw->ycoord[n_root->v[0]->num]);
  if(n_root->v[0]->tax) fprintf(fp,"(%s) show\n",n_root->v[0]->name);
  else
    {
      For(i,3)
	if((n_root->v[0]->v[i]) && (n_root->v[0]->v[i] != n_root->v[1]))
	  DR_Print_Tree_Postscript_Pre(n_root->v[0],n_root->v[0]->v[i],fp,draw,tree);
    }

  fprintf(fp,"%d %d mt\n",draw->xcoord[n_root->v[1]->num],draw->ycoord[n_root->v[1]->num]);

  if(n_root->v[1]->tax) fprintf(fp,"(%s) show\n",n_root->v[1]->name);
  else
    {
      For(i,3)
	if((n_root->v[1]->v[i]) && (n_root->v[1]->v[i] != n_root->v[0]))
	  DR_Print_Tree_Postscript_Pre(n_root->v[1],n_root->v[1]->v[i],fp,draw,tree);
    }

  fprintf(fp,"closepath\n");
  fprintf(fp,"stroke\n");
  fprintf(fp,"showpage\n");
}

/*********************************************************/

void DR_Print_Tree_Postscript_Pre(node *a, node *d, FILE *fp, tdraw *w, arbre *tree)
{
  int i;

  fprintf(fp,"gsave\n");
  
  For(i,3)
    if(a->v[i] == d)
      {
/* 	if(a->b[i]->prob_sel_regime <= 0.1) */
/* 	  fprintf(fp,".0 .0 1. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.1 && a->b[i]->prob_sel_regime <= 0.2) */
/* 	  fprintf(fp,".0 .5 1. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.2 && a->b[i]->prob_sel_regime <= 0.3) */
/* 	  fprintf(fp,".0 1. 1. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.3 && a->b[i]->prob_sel_regime <= 0.4) */
/* 	  fprintf(fp,".0 1. .5 sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.4 && a->b[i]->prob_sel_regime <= 0.5) */
/* 	  fprintf(fp,".0 1. .0 sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.5 && a->b[i]->prob_sel_regime <= 0.6) */
/* 	  fprintf(fp,".5 1. .0 sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.6 && a->b[i]->prob_sel_regime <= 0.7) */
/* 	  fprintf(fp,"1. 1. 0. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.7 && a->b[i]->prob_sel_regime <= 0.8) */
/* 	  fprintf(fp,"1. .5 0. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.8 && a->b[i]->prob_sel_regime <= 0.9) */
/* 	  fprintf(fp,"1. 0. 0. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.9) */
/* 	  fprintf(fp,"1. .0 .0 sc\n"); */
	break;
      }

  fprintf(fp,"%d %d mt\n",w->xcoord[a->num],w->ycoord[a->num]);
  fprintf(fp,"%d %d lt\n",w->xcoord[a->num],w->ycoord[d->num]);
  fprintf(fp,"%d %d lt\n",w->xcoord[d->num],w->ycoord[d->num]);

  if(d->tax) 
    {
      fprintf(fp,"(%s) show \n",d->name);
      fprintf(fp,"stroke\n");
      fprintf(fp,"grestore\n");
      return;
    }
  else
    {
      fprintf(fp,"stroke\n");
      fprintf(fp,"grestore\n");
      For(i,3)
	if(d->v[i] != a) DR_Print_Tree_Postscript_Pre(d,d->v[i],fp,w,tree);
    }


  return;
}

/*********************************************************/

void DR_Dist_To_Root_Pre(node *a, node *d, edge *b, arbre *tree)
{
  int i;

  if(b) d->dist_to_root = a->dist_to_root + b->l;

  if(d->tax) return;
  else
    {
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  DR_Dist_To_Root_Pre(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/

void DR_Dist_To_Root(node *n_root, arbre *tree)
{  
  n_root->v[0]->dist_to_root = tree->e_root->l * tree->n_root_pos;
  n_root->v[1]->dist_to_root = tree->e_root->l * (1.-tree->n_root_pos);
  DR_Dist_To_Root_Pre(n_root,n_root->v[0],NULL,tree);
  DR_Dist_To_Root_Pre(n_root,n_root->v[1],NULL,tree);
}

/*********************************************************/

void DR_Get_X_Coord_Pre(node *a, node *d, edge *b, tdraw *w, arbre *tree)
{
  int i;

  if(b) w->xcoord[d->num] =  d->dist_to_root * (double)w->tree_box_width/w->max_dist_to_root;

  if(d->tax) return;
  else
    {
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  DR_Get_X_Coord_Pre(d,d->v[i],d->b[i],w,tree);
    }
}

/*********************************************************/

void DR_Get_X_Coord(tdraw *w, arbre *tree)
{
  w->xcoord[tree->n_root->v[0]->num] = tree->n_root->v[0]->dist_to_root * (double)w->tree_box_width/w->max_dist_to_root;
  w->xcoord[tree->n_root->v[1]->num] = tree->n_root->v[1]->dist_to_root * (double)w->tree_box_width/w->max_dist_to_root;
  DR_Get_X_Coord_Pre(tree->n_root,tree->n_root->v[0],NULL,w,tree);
  DR_Get_X_Coord_Pre(tree->n_root,tree->n_root->v[1],NULL,w,tree);
}


/*********************************************************/

void DR_Get_Y_Coord(tdraw *w, arbre *tree)
{
  int next_y_slot;
  next_y_slot = 0;
  DR_Get_Y_Coord_Post(tree->e_root->left,tree->e_root->rght,tree->e_root,&next_y_slot,w,tree);
  DR_Get_Y_Coord_Post(tree->e_root->rght,tree->e_root->left,tree->e_root,&next_y_slot,w,tree);
}

/*********************************************************/

void DR_Get_Y_Coord_Post(node *a, node *d, edge *b, int *next_y_slot, tdraw *w, arbre *tree)
{
  int i;

  if(d->tax) 
    {
      w->ycoord[d->num] = *next_y_slot + (int)(w->page_height / (2.*tree->n_otu));
      (*next_y_slot) += (int)(w->page_height / (tree->n_otu));
    }
  else
    {
      int d1, d2;

      d1 = d2 = -1;
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      DR_Get_Y_Coord_Post(d,d->v[i],d->b[i],next_y_slot,w,tree);
	      if(d1<0) d1 = i;
	      else     d2 = i;
	    }
	}
      w->ycoord[d->num] = (w->ycoord[d->v[d1]->num] + w->ycoord[d->v[d2]->num])/2.; 
    }
}

/*********************************************************/

tdraw *DR_Make_Tdraw_Struct(arbre *tree)
{
  tdraw *w;

  w = (tdraw *)mCalloc(1,sizeof(tdraw));
  w->xcoord = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  w->ycoord = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));

  return w;
}

/*********************************************************/

void DR_Init_Tdraw_Struct(tdraw *w)
{
  w->page_width  = 510;
  w->page_height = 770;
}

/*********************************************************/

void DR_Get_Tree_Box_Width(tdraw *w, arbre *tree)
{
  int i;
  int max_name_len, curr_len;

  max_name_len = curr_len = 0;
  For(i,tree->n_otu)
    {
      curr_len = (int)strlen(tree->noeud[i]->name);
      if(curr_len > max_name_len) max_name_len = curr_len;
    }

  w->tree_box_width = w->page_width - max_name_len * 8.66667;
}

/*********************************************************/

double DR_Get_Max_Dist_To_Root(arbre *tree)
{
  double mx;
  int i;

  mx = .0;
  For(i,tree->n_otu)
    {
      if(tree->noeud[i]->dist_to_root > mx)
	{
	  mx = tree->noeud[i]->dist_to_root;
	}
    }

  return mx;
}

/*********************************************************/

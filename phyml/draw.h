#ifndef DRAW_H
#define DRAW_H

void DR_Dist_To_Root_Pre(node *a, node *d, edge *b, arbre *tree);
void DR_Dist_To_Root(node *n_root, arbre *tree);
void DR_Get_X_Coord_Pre(node *a, node *d, edge *b, tdraw *w, arbre *tree);
void DR_Get_X_Coord(tdraw *w, arbre *tree);
tdraw *DR_Make_Tdraw_Struct(arbre *tree);
void DR_Init_Tdraw_Struct(tdraw *d);
void DR_Get_Tree_Box_Width(tdraw *w, arbre *tree);
void DR_Get_Y_Coord_Post(node *a, node *d, edge *b, int *next_y_slot, tdraw *w, arbre *tree);
void DR_Get_Y_Coord(tdraw *w, arbre *tree);
void DR_Get_Tree_Coord(arbre *tree);
double DR_Get_Max_Dist_To_Root(arbre *tree);
void DR_Print_Tree_Postscript(int tree_num, FILE *fp, arbre *tree);
void DR_Print_Tree_Postscript_Pre(node *a, node *d, FILE *fp, tdraw *w, arbre *tree);
void DR_Print_Postscript_EOF(FILE *fp);
void DR_Print_Postscript_Header(int n_pages, FILE *fp);


#endif

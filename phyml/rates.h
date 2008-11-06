/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef RATES_H
#define RATES_H

void RATES_Monte_Carlo_Mean_Rates(arbre *tree);
void RATES_Monte_Carlo_Mean_Rates_Pre(node *a, node *d, edge *b, phydbl curr_rate, arbre *tree);
void RATES_Print_Rates(arbre *tree);
void RATES_Print_Rates_Pre(node *a, node *d, edge *b, arbre *tree);
trate *RATES_Make_Rate_Struct(arbre *tree);
void RATES_Init_Rate_Struct(trate *rates, arbre *tree);

#endif

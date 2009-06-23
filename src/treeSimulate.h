/*  Copyright 2009, Stefan Washietl

    This file is part of RNAcode.

    RNAcode is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RNAcode is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RNAcode.  If not, see <http://www.gnu.org/licenses/>. */

#include "tree.h"
#include "rnaz_utils.h"

TNode* getLCA(TTree* tree, TNode* nodeA, TNode* nodeB);
void tree2aln(TTree* tree, struct aln *alignment[]);
TTree* string2tree(char* treeString);

void simulateTree(TTree* tree, float freqs[], float kap, int L);

float** getDistanceMatrix(TTree* tree, struct aln *alignment[]);

void freeSeqgenTree(TTree* tree);

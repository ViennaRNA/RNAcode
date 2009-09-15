#ifndef _POSTSCRIPT_H_
#define _POSTSCRIPT_H_

#include "rnaz_utils.h"
#include "code.h"

int colorAln(const char *filename, const struct aln *alignment[],segmentStats* region);

void colorHSS(const char *filename, const struct aln *alignment[],segmentStats *region, int b, int i);

void setParameters();


#endif

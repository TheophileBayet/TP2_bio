#ifndef _SWOUT_H_
#define _SWOUT_H_

#include "swCalc.h"
#include "swCost.h"
#include <string.h>
/* Find highest scoring local alignment(s) in mat, and print to stdout
   the corresponding best alignments.
   mat must have been filled with scores and prevs.
   cost is provided so mismatches with negative scores can be lowercased.
   */
void printAllBests(struct matrix *mat, struct cost *cost, char *s1, char *s2) ;
void printAllBestsGotoh(struct matrix *D, struct matrix *V,
        struct matrix *H, /*struct cost *cost,*/char *s1, char *s2);

void printAllPaths(int i, int j, double score, struct matrix *mat, struct cost *cost, char *s1, char *s2);
void printAllPathsGotoh(int i, int j, double score, struct matrix *D, struct matrix *V, struct matrix *H, /*struct cost *cost,*/char *s1, char *s2);
#endif

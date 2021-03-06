#ifndef _SWCALC_H_
#define _SWCALC_H_

#include <stdint.h>
#include "swCost.h"

/* defines the SW matrix datatypes */


/* in struct cell, the three last bits of prevs are used to
   say whether best paths come from Top, Left, or Diag (in that order).
   So eg if (prevs&4) then a best path comes from the
   cell above, and if (prevs&1) a best path comes from diag.
   The three bits are non-exclusive (can have multiple best paths).
   If score==0 then prevs must be 0.
*/
struct cell {
	double score;
	uint16_t prevs;
};

struct matrix {
	unsigned int w;
	unsigned int h;
	struct cell *cells; /* pointer to array of w*h cells
			       cells[w*i+j] contains cell (i,j) */
};

/* allocate and initialize (first row and col) a matrix for SW
   alignment of strings s1 and s2 */
struct matrix *swInitMat(char *s1, char *s2);

/* Fill the mat matrix, using Smith-Waterman with a linear indel model
   using cost->indelOpen, or Gotoh with an affine indel model using
   cost->indelOpen and cost->indelExtend.
   Preconditions:
   - mat is correctly allocated and initialized (by swInitMat)
   - cost->subst is defined for each pair of letters in s1 and s2
*/
void swFillMat(struct matrix *mat, struct cost *cost, char *s1, char *s2) ;


void swFillDVH(struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) ;


/* Fill the D matrix, using Gotoh with an affine indel model using
   cost->indelOpen and cost->indelExtend.
   Preconditions:
   - D is correctly allocated and initialized (by swInitMat)
   - cost->subst is defined for each pair of letters in s1 and s2
*/
void swFillDMat(int i, int j, struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) ;
//void swFillDMat(struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) ;

/* Fill the V matrix, using Gotoh with an affine indel model using
   cost->indelOpen and cost->indelExtend.
   Preconditions:
   - V is correctly allocated and initialized (by swInitMat)
   - cost->subst is defined for each pair of letters in s1 and s2
*/

void swFillVMat(int i, int j, struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) ;
//void swFillVMat(struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) ;

/* Fill the H matrix, using Gotoh with an affine indel model using
   cost->indelOpen and cost->indelExtend.
   Preconditions:
   - H is correctly allocated and initialized (by swInitMat)
   - cost->subst is defined for each pair of letters in s1 and s2
*/
void swFillHMat(int i, int j, struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) ;
//void swFillHMat(struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) ;

/* free all allocated memory in mat */
void swFreeMat(struct matrix *mat);

/* print contents of matrix, for debugging */
void swPrintMat(struct matrix *mat);

#endif

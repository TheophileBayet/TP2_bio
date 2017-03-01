#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "mem.h"
#include "swCalc.h"



/*******************************
  define extern functions
 *******************************/

/* allocate and initialize (first row and col) a matrix for SW
   alignment of strings s1 and s2 */
struct matrix *swInitMat(char *s1, char *s2){
  struct matrix *mat = mallocOrDie(sizeof(struct matrix), "E: failed to create matrix\n") ;
  int w = strlen(s1);
  int h = strlen(s2);
  mat -> w =  w;
  mat -> h = h ;
  struct cell *cells= mallocOrDie(w*h*sizeof(struct cell),"E: failed to create the cells tab\n");
  /* pointer to array of w*h cells
			       cells[w*i+j] contains cell (i,j) */
  for(int i = 0;i<w;i++){
    cells[i].score =0.0;
    cells[i].prevs=0;
    cells[w*i].score=0.0;
    cells[w*i].prevs=0;
  }
  return mat;
}


/* Fill the mat matrix, using Smith-Waterman with a linear indel model
   using cost->indelOpen, or Gotoh with an affine indel model using
   cost->indelOpen and cost->indelExtend.
   Preconditions:
   - mat is correctly allocated and initialized (by swInitMat)
   - cost->subst is defined for each pair of letters in s1 and s2
*/
void swFillMat(struct matrix *mat, struct cost *cost, char *s1, char *s2) {

}

/* free all allocated memory in mat */
void swFreeMat(struct matrix *mat){
    if (mat != NULL && mat->cells != NULL) {
        free(mat->cells);
    }
    if (mat != NULL) {
        free(mat);
    }
}

/* print contents of matrix, for debugging */
void swPrintMat(struct matrix *mat){
    for (unsigned int i = 0 ; i < mat->h ; i++) {
        for (unsigned int j = 0 ; j < mat->w ; j++) {
            fprintf(stdout, "%f\t", mat->cells[mat->w*i+j].score); 
        }
        fprintf(stdout, "\n"); 
    }
}

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
  /* récupérer les longueurs et initialiser la matrice */
  /* premier jet : taille 4x4*/
  struct matrix *mat = mallocOrDie(sizeof(struct matrix), "E: failed to create matrix\n") ;
  mat -> w = 4 ;
  mat -> h = 4 ;/*
  for (int i = 0;i<w*h-1;i++){
    struct cell cell = mallocOrDie(sizeof(struct cell),"E: failed to create cell")
    mat -> cells[w*i+h]=cell;
  }*/
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
}

/* print contents of matrix, for debugging */
void swPrintMat(struct matrix *mat){

}

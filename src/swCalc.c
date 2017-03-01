#include <stdio.h>
#include <string.h>
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
  int w = strlen(s1)+1;
  int h = strlen(s2)+1;
  mat -> w =  w;
  mat -> h = h ;
  struct cell *cells= mallocOrDie(w*h*sizeof(struct cell),"E: failed to create the cells tab\n");
  /* pointer to array of w*h cells
			       cells[w*i+j] contains cell (i,j) */
  for(int i = 0;i<w;i++){
    cells[i].score =0.0;
    cells[i].prevs=0;
  }
for(int i = 0;i<h;i++){
    cells[w*i].score=0.0;
    cells[w*i].prevs=0;
  }
  mat->cells =cells;
  return mat;
}

double* MAX(struct cell c1, struct cell c2, struct cell c3){
  double *ret = mallocOrDie(4*sizeof(double),"E:failed to create tab in MAX");
  //TODO : implémenter la fonction MAX pour qu'elle retourne bien ce qui est demandé.
  ret[0]=c1.score+c2.score+c3.score;
  ret[1]=0;
  ret[2]=0;
  ret[3]=0;
  return ret;
}


/* Fill the mat matrix, using Smith-Waterman with a linear indel model
   using cost->indelOpen, or Gotoh with an affine indel model using
   cost->indelOpen and cost->indelExtend.
   Preconditions:
   - mat is correctly allocated and initialized (by swInitMat)
   - cost->subst is defined for each pair of letters in s1 and s2
*/
void swFillMat(struct matrix *mat, struct cost *cost, char *s1, char *s2) {
  int w = strlen(s1)+1;
  int h = strlen(s2)+1;
  for(int i =1;i<w;i++){
    for(int j=1;j<h;j++){
      // Fonction MAX sur les trois précédents, conserver la valeur dans la cell
      // et modifier les prevs !
      double* max_tab;
      max_tab=MAX(mat->cells[w*(i-1)+j-1],mat->cells[w*(i-1)+j],mat->cells[w*i+j-1]);
      mat->cells[w*i+j].score=max_tab[0];
      mat->cells[w*i+j].prevs=0;
      if(max_tab[1]){
        mat->cells[w*i+j].prevs+=1;
      }
      if(max_tab[2]){
        mat->cells[w*i+j].prevs+=2;
      }
      if(max_tab[3]){
        mat->cells[w*i+j].prevs+=4;
      }
    }
  }
}



/* free all allocated memory in mat */
void swFreeMat(struct matrix *mat){
}

/* print contents of matrix, for debugging */
void swPrintMat(struct matrix *mat){

}

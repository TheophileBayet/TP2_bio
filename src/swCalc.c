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
  for(int i = 0 ; i<w ; i++){
      cells[i].score =0.0;
      cells[i].prevs=0;
  }
  for(int i = 0; i<h ;i++){
      cells[w*i].score=0.0;
      cells[w*i].prevs=0;
  }
  mat->cells =cells;
  return mat;
}

double* MAX(struct cell c1, struct cell c2, struct cell c3, struct cost *cost,char s1, char s2){
  double *ret = mallocOrDie(4*sizeof(double),"E:failed to create tab in MAX");
  ret[1]=0;
  ret[2]=0;
  ret[3]=0;
  int sc1 ;
  int sc2;
  int sc3;
  sc1 = c1.score + cost->subst(s1,s2);
  sc2 = c2.score + cost->indelOpen ;
  sc3 = c3.score + cost->indelOpen;
  if(sc1<sc2){
    if (sc2<sc3){ret[0]=sc3; ret[3]=1;}
    else{ret[0]=sc2; ret[2]=1;}
  }else{
    if(sc1<sc3){ret[0]=sc3; ret[3]=1;}
    else{ret[0]=sc1;ret[1]=1;}
  }
  if(ret[0]<0){
    ret[0]=0;
    ret[1]=0;
    ret[2]=0;
    ret[3]=0;
  }
  return ret;
}



double* MAX_D(struct cell c1, struct cell c2, struct cell c3, struct cost *cost,char s1, char s2){
  double *ret = mallocOrDie(4*sizeof(double),"E:failed to create tab in MAX");
  ret[1]=0;
  ret[2]=0;
  ret[3]=0;
  int sc1=0;
  int sc2=0;
  int sc3=0;
  sc1 = c1.score + cost->subst(s1,s2);
  sc2 = c2.score + cost->subst(s1,s2);
  sc3 = c3.score + cost->subst(s1,s2);
  if(sc1<sc2){
    if (sc2<sc3){ret[0]=sc3; ret[3]=1;}
    else{ret[0]=sc2; ret[2]=1;}
  }else{
    if(sc1<sc3){ret[0]=sc3; ret[3]=1;}
    else{ret[0]=sc1;ret[1]=1;}
  }
  if(ret[0]<0){
    ret[0]=0;
    ret[1]=0;
    ret[2]=0;
    ret[3]=0;
  }
  return ret;
}

double* MAX_V(struct cell c1, struct cell c2, struct cell c3, struct cost *cost){
  double *ret = mallocOrDie(4*sizeof(double),"E:failed to create tab in MAX");
  ret[1]=0;
  ret[2]=0;
  ret[3]=0;
  int sc1 ;
  int sc2;
  int sc3;
  sc1 = c1.score + cost->indelOpen;
  sc2 = c2.score + cost->indelExtend;
  sc3 = c3.score + cost->indelOpen;
  if(sc1<sc2){
    if (sc2<sc3){ret[0]=sc3; ret[3]=1;}
    else{ret[0]=sc2; ret[2]=1;}
  }else{
    if(sc1<sc3){ret[0]=sc3; ret[3]=1;}
    else{ret[0]=sc1;ret[1]=1;}
  }
  if(ret[0]<-1){
    ret[0]=-1;
    ret[1]=0;
    ret[2]=0;
    ret[3]=0;
  }
  return ret;
}

double* MAX_H(struct cell c1, struct cell c2, struct cell c3, struct cost *cost){
  double *ret = mallocOrDie(4*sizeof(double),"E:failed to create tab in MAX");
  ret[1]=0;
  ret[2]=0;
  ret[3]=0;
  int sc1 ;
  int sc2;
  int sc3;
  sc1 = c1.score + cost->indelOpen;
  sc2 = c2.score + cost->indelOpen;
  sc3 = c3.score + cost->indelExtend;
  if(sc1<sc2){
    if (sc2<sc3){ret[0]=sc3; ret[3]=1;}
    else{ret[0]=sc2; ret[2]=1;}
  }else{
    if(sc1<sc3){ret[0]=sc3; ret[3]=1;}
    else{ret[0]=sc1;ret[1]=1;}
  }
  if(ret[0]<-1){
    ret[0]=-1;
    ret[1]=0;
    ret[2]=0;
    ret[3]=0;
  }
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
  for(int i =1;i<h;i++){
    for(int j=1;j<w;j++){
      // Fonction MAX sur les trois précédents, conserver la valeur dans la cell
      // et modifier les prevs !
      double* max_tab;
      max_tab=MAX(mat->cells[w*(i-1)+j-1],mat->cells[w*(i-1)+j],mat->cells[w*i+j-1],cost,s1[j-1],s2[i-1]);
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
      free(max_tab);
    }
  }
}


void swFillDMat(struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) {
  int w = strlen(s1)+1;
  int h = strlen(s2)+1;
  // TODO : modifier les accès aux cellules et impélmenter le MAX_D
  for(int i =1;i<h;i++){
    for(int j=1;j<w;j++){
      // Fonction MAX sur les trois précédents, conserver la valeur dans la cell
      // et modifier les prevs !
      double* max_tab;
      max_tab=MAX_D(D->cells[w*(i-1)+j-1],V->cells[w*(i-1)+j-1],H->cells[w*(i-1)+j-1],cost,s1[j-1],s2[i-1]);
      D->cells[w*i+j].score=max_tab[0];
      D->cells[w*i+j].prevs=0;
      if(max_tab[1]){
        D->cells[w*i+j].prevs+=1;
      }
      if(max_tab[2]){
        D->cells[w*i+j].prevs+=2;
      }
      if(max_tab[3]){
        D->cells[w*i+j].prevs+=4;
      }
      free(max_tab);
    }
  }
}

void swFillVMat(struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) {
  int w = strlen(s1)+1;
  int h = strlen(s2)+1;
  for(int i =1;i<h;i++){
    for(int j=1;j<w;j++){
      // Fonction MAX sur les trois précédents, conserver la valeur dans la cell
      // et modifier les prevs !
      double* max_tab;
      max_tab=MAX_V(D->cells[w*(i-1)+j],V->cells[w*(i-1)+j],H->cells[w*(i-1)+j],cost);
      V->cells[w*i+j].score=max_tab[0];
      V->cells[w*i+j].prevs=0;
      if(max_tab[1]){
        V->cells[w*i+j].prevs+=1;
      }
      if(max_tab[2]){
        V->cells[w*i+j].prevs+=2;
      }
      if(max_tab[3]){
        V->cells[w*i+j].prevs+=4;
      }
      free(max_tab);
    }
  }
}

void swFillHMat(struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) {
  int w = strlen(s1)+1;
  int h = strlen(s2)+1;
  for(int i =1;i<h;i++){
    for(int j=1;j<w;j++){
      // Fonction MAX sur les trois précédents, conserver la valeur dans la cell
      // et modifier les prevs !
      double* max_tab;
      max_tab=MAX_H(D->cells[w*i+j-1],V->cells[w*i+j-1],H->cells[w*i+j-1],cost);
      H->cells[w*i+j].score=max_tab[0];
      H->cells[w*i+j].prevs=0;
      if(max_tab[1]){
        H->cells[w*i+j].prevs+=1;
      }
      if(max_tab[2]){
        H->cells[w*i+j].prevs+=2;
      }
      if(max_tab[3]){
        H->cells[w*i+j].prevs+=4;
      }
      free(max_tab);
    }
  }
}

/* free all allocated memory in mat */
void swFreeMat(struct matrix *mat){
    fprintf(stdout, "\t\t\tswFreeMat\n");
    if (mat != NULL && mat->cells != NULL) {
        free(mat->cells);
    }
    if (mat != NULL) {
        free(mat);
    }
}

/* print contents of matrix, for debugging */
void swPrintMat(struct matrix *mat){
    if (mat != NULL && mat->cells != NULL) {
        fprintf(stdout, "Matrice des scores\n");
        for (unsigned int i = 0 ; i < mat->h ; i++) {
            for (unsigned int j = 0 ; j < mat->w ; j++) {
                fprintf(stdout, "\t%f", mat->cells[mat->w*i+j].score);
            }
            fprintf(stdout, "\n");
        }
    } else {
        fprintf(stdout, "La matrice n'est pas correctement initialisée");
    }
    fprintf(stdout, "\n");
}

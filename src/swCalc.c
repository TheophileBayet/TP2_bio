#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

#include "mem.h"
#include "swCalc.h"




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

void setRet(double* ret, int sc1 , int sc2, int sc3) {
    ret[1]=0;
    ret[2]=0;
    ret[3]=0;
    if(sc1<sc2){
        if (sc2<sc3){ret[0]=sc3; ret[3]=1;}
        else{ret[0]=sc2; ret[2]=1;}
    }else{
        if(sc1<sc3){ret[0]=sc3; ret[3]=1;}
        else{ret[0]=sc1;ret[1]=1;}
    }
    if (ret[0] == sc1) {ret[1] = 1;}
    if (ret[0] == sc2) {ret[2] = 1;}
    if (ret[0] == sc3) {ret[3] = 1;}
    if(ret[0]<0){
        ret[0]=0;
        ret[1]=0;
        ret[2]=0;
        ret[3]=0;
    }
}
double* max(struct cell c1, struct cell c2, struct cell c3, struct cost *cost,char s1, char s2){
    double *ret = mallocOrDie(4*sizeof(double),"E:failed to create tab in max");
    int sc1 = c1.score + cost->subst(s1,s2);
    int sc2 = c2.score + cost->indelOpen;
    int sc3 = c3.score + cost->indelOpen;
    setRet(ret, sc1, sc2, sc3);
    return ret;
}

double* max_D(struct cell c1, struct cell c2, struct cell c3, struct cost *cost,char s1, char s2){
    double *ret = mallocOrDie(4*sizeof(double),"E:failed to create tab in max");
    int sc1 = c1.score + cost->subst(s1,s2);
    int sc2 = c2.score + cost->subst(s1,s2);
    int sc3 = c3.score + cost->subst(s1,s2);
    setRet(ret, sc1, sc2, sc3);
    return ret;
}

double* max_V(struct cell c1, struct cell c2, struct cell c3, struct cost *cost){
    double *ret = mallocOrDie(4*sizeof(double),"E:failed to create tab in max");
    int sc1 = c1.score + cost->indelOpen;
    int sc2 = c2.score + cost->indelExtend;
    int sc3 = c3.score + cost->indelOpen;
    setRet(ret, sc1, sc2, sc3);
    return ret;
}

double* max_H(struct cell c1, struct cell c2, struct cell c3, struct cost *cost){
    double *ret = mallocOrDie(4*sizeof(double),"E:failed to create tab in max");
    int sc1 = c1.score + cost->indelOpen;
    int sc2 = c2.score + cost->indelOpen;
    int sc3 = c3.score + cost->indelExtend;
    setRet(ret, sc1, sc2, sc3);
    return ret;
}

void swFillCell(struct matrix *mat, double *max_tab, int w, int i, int j) {
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
            // Fonction max sur les trois précédents, conserver la valeur dans la cell
            // et modifier les prevs !
            double* max_tab = max(mat->cells[w*(i-1)+j-1],mat->cells[w*(i-1)+j],mat->cells[w*i+j-1],cost,s1[j-1],s2[i-1]);
            swFillCell(mat, max_tab, w, i, j);
            free(max_tab);
        }
    }
}

void swFillDVH(struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) {
    int w = strlen(s1)+1;
    int h = strlen(s2)+1;

    for(int i =1;i<h;i++){
        for(int j=1;j<w;j++){
            swFillDMat(i,j, D, V, H, cost, s1, s2);
            swFillVMat(i,j, D, V, H, cost, s1, s2);
            swFillHMat(i,j, D, V, H, cost, s1, s2);
        }
    }
}

void swFillDMat(int i, int j, struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) {
    int w = strlen(s1)+1;
    // Fonction max sur les trois précédents, conserver la valeur dans la cell
    // et modifier les prevs !
    double* max_tab = max_D(D->cells[w*(i-1)+j-1],V->cells[w*(i-1)+j-1],H->cells[w*(i-1)+j-1],cost,s1[j-1],s2[i-1]);
    swFillCell(D, max_tab, w, i, j);
    free(max_tab);
}

void swFillVMat(int i, int j, struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) {
    int w = strlen(s1)+1;
    // Fonction max sur les trois précédents, conserver la valeur dans la cell
    // et modifier les prevs !
    double* max_tab = max_V(D->cells[w*(i-1)+j],V->cells[w*(i-1)+j],H->cells[w*(i-1)+j],cost);
    if (max_tab[0] == 0) {
        max_tab[0] = -1;
    }
    swFillCell(V, max_tab, w, i, j);
    free(max_tab);
}

void swFillHMat(int i, int j, struct matrix *D,struct matrix *V,struct matrix *H, struct cost *cost, char *s1, char *s2) {
    int w = strlen(s1)+1;
    // Fonction max sur les trois précédents, conserver la valeur dans la cell
    // et modifier les prevs !
    double* max_tab = max_H(D->cells[w*i+j-1],V->cells[w*i+j-1],H->cells[w*i+j-1],cost);
    if (max_tab[0] == 0) {
        max_tab[0] = -1;
    }
    swFillCell(H, max_tab, w, i, j);
    free(max_tab);
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

#include <stdio.h>

#include "swOut.h"

void printBestAlis(struct matrix *mat, struct cost *cost, char *s1, char *s2)
{
    fprintf(stdout, "printBestAlis\n");

    // On commence tout en bas
    unsigned int i = mat->h-1; // -1 ?
    unsigned int j = mat->w-1; // same

    char s1_res[strlen(s1)+strlen(s2)]; // résultat avec letters et '-' pour s1
    char s2_res[strlen(s1)+strlen(s2)]; // résultat avec letters et '-' pour s2

    strcpy(s1_res, s1);
    strcat(s1_res, s2);
    strcpy(s2_res, s1);
    strcat(s2_res, s2);

    unsigned int count = strlen(s1_res) - 1;

    // on commence à la case max
    //

    int max = 0;
    for (uint32_t k1 = 0; k1<mat->h; k1++){
        for (uint32_t k2 = 0; k2<mat->w; k2++){
            if (mat->cells[mat->w*k1+k2].score>max){
                i=k1;
                j=k2;
                max=mat->cells[mat->w*k1+k2].score;
            }
        }
    }
    printf("max = %d en %d, %d \n",max,i,j);
    int ideb=i-1;
    int jdeb=j-1;
    int ifin=0;
    int jfin=0;
    // en considérant la case à laquelle on commence, on peut remplir en partie les chaines de résultat
    // remplir à partir d'où jusqu'à où ? OU NE PAS REMPLIR LOL parce que ça ne compte pas


    while (i!=0 && j!=0){

        //un seul chemin, on privilégie la diag (prevs&1) (arbitraire)
        // à changer pour tous les chemins

        struct cell cur = mat->cells[i*mat->w+j];

        //printf("Résultat intermédiaire : %s\n", s1_res + count + 1);
        //printf("Résultat intermédiaire : %s\n\n", s1_res);
        // màj des chemins selon déplacement optimal trouvé

        switch(cur.prevs){

            //diag
            case 1 :
            case 3 :
            case 5 :
            case 7 :
                if (s1[j-1] == s2[i-1]) {
                    s1_res[count] = s1[j-1];
                    s2_res[count] = s2[i-1];
                } else {
                    s1_res[count] = s1[j-1] + 32;
                    s2_res[count] = s2[i-1] + 32;
                }
                i--; j--;
                break;

            case 2 :
            case 6 :
                s1_res[count] = '-';
                s2_res[count] = s2[i-1] + 32;
                i--;
                break;

            case 4 :
                s1_res[count] = s1[j-1] + 32;
                s2_res[count] = '-';
                j--;
                break;

            default :
                //fprintf(stderr, "%s\nScore %u\n", "Pas de cas correspondant pour prevs", cur.prevs);
                //exit(1);
                ifin = i ; jfin = j;
                i = 0; j = 0;
                count++;
                break;

        }
        count--;
    }
    //printf("Résultat final : \n\t%s\n\t%s\n", s1_res + count + 1, s2_res + count + 1);
    printf("Final result : \n");
    printf("s1_res : %s\ns2_res : %s\n", s1_res + count + 1, s2_res + count + 1);
    printf("Best match is at s1[%d:%d] and s2[%d:%d]\n", jfin, jdeb, ifin, ideb);

}


void printBestAlisGotoh(struct matrix *D, struct matrix *V, struct matrix *H, /*struct cost *cost,*/char *s1, char *s2)
{
    fprintf(stdout, "printBestAlisGotoh\n");

    unsigned int i = D->h-1;
    unsigned int j = D->w-1;

    char s1_res[strlen(s1)+strlen(s2)]; // résultat avec letters et '-' pour s1
    char s2_res[strlen(s1)+strlen(s2)]; // résultat avec letters et '-' pour s2

    strcpy(s1_res, s1);
    strcat(s1_res, s2);
    strcpy(s2_res, s1);
    strcat(s2_res, s2);

    unsigned int count = strlen(s1_res) - 1;

        int max = 0;
    for (uint32_t k1 = 0; k1<D->h; k1++){
        for (uint32_t k2 = 0; k2<D->w; k2++){
            if (D->cells[D->w*k1+k2].score>max){
                i=k1;
                j=k2;
                max=D->cells[D->w*k1+k2].score;
            }
        }
    }

        printf("max = %d en %d, %d \n",max,i,j);
    int ideb=i-1;
    int jdeb=j-1;
    struct cell cur = D->cells[i*D->w+j];
    while (i!=0 && j!=0){
        /*
        printf("Résultat intermédiaire : %s\n", s1_res + count + 1);
        printf("Résultat intermédiaire : %s\n\n", s1_res);
        // màj des chemins selon déplacement optimal trouvé        struct cell cur = mat->cells[i*mat->w+j];

        printf("Résultat intermédiaire : %s\n", s1_res + count + 1);
        printf("Résultat intermédiaire : %s\n\n", s1_res);
        // màj des chemins selon déplacement optimal trouvé
        */
        switch(cur.prevs){

            case 1 :
            case 3 :
            case 5 :
            case 7 :
                if (s1[j-1] == s2[i-1]) {
                    s1_res[count] = s1[j-1];
                    s2_res[count] = s2[i-1];
                } else {
                    s1_res[count] = s1[j-1] + 32;
                    s2_res[count] = s2[i-1] + 32;
                }
                i--; j--;
                cur = D->cells[i*D->w+j];
                break;

            case 2 :
            case 6 :
                s1_res[count] = '-';
                s2_res[count] = s2[i-1] + 32;
                i--;
                cur = V->cells[i*V->w+j];
                break;

            case 4 :
                s1_res[count] = s1[j-1] + 32;
                s2_res[count] = '-';
                j--;
                cur = H->cells[i*H->w+j];
                break;

            default :
                //fprintf(stderr, "%s\nScore %u\n", "Pas de cas correspondant pour prevs", cur.prevs);
                //exit(1);
                i = 0; j = 0;
                count++;
                break;
        }
    count--;
    }
    printf("Final result : \n");
    printf("s1_res : %s\ns2_res : %s\n", s1_res + count + 1, s2_res + count + 1);
    printf("Best match is at s1[:] and s2[:]\n");

}

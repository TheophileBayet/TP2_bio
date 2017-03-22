#include <stdio.h>

#include "swOut.h"

void printBestAlis(struct matrix *mat, struct cost *cost, char *s1, char *s2) 
{
    fprintf(stdout, "printBestAlis\n");

    // on commence tout en bas (for now)
    unsigned int i = mat->h-1; // -1 ?
    unsigned int j = mat->w-1; // same

    char s1_res[sizeof(s1)+sizeof(s2)]; // résultat avec letters et '-' pour s1
    char s2_res[sizeof(s1)+sizeof(s2)]; // résultat avec letters et '-' pour s1

    while (i!=0 || j!=0){

        //un seul chemin, on privilégie la diag (prevs&1) (arbitraire)
        // à changer pour tous les chemins

        struct cell cur = mat->cells[i*mat->w+j];

        printf("Résultat intermédiaire : %.*s\n", mat->h-1 + mat->w-1 - j+i, s1_res + j+i);
        // màj des chemins selon déplacement optimal trouvé

        switch(cur.prevs){

            //diag
            case 1 :
                printf("Case 1\n");
                s1_res[i+j] = s1[j-1];
                s2_res[i+j] = s2[i-1];
                //s1_res = strcat(s1[j-1], s1_res);
                //s2_res = strcat(s2[i-1], s2_res);
                i--; j--;
                break;

            case 2 :
                printf("Case 2\n");
                s1_res[i+j] = s1[j-1];
                s2_res[i+j] = '-';
                //s1_res = strcat(s1[j-1],s1_res);
                //s2_res = strcat("-",s2_res);
                j--;
                break;

            case 4 :
                printf("Case 4\n");
                s1_res[i+j] = '-';
                s2_res[i+j] = s2[i-1];
                //s1_res = strcat("-",s1_res);
                //s2_res = strcat(s2[i-1],s2_res);
                i--;
                break;

            default :
                fprintf(stderr, "%s\n", "Pas de cas correspondant pour prevs");
                exit(1);

        }
    }

}

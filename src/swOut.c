#include <stdio.h>

#include "swOut.h"

void printBestAlis(struct matrix *mat, struct cost *cost, char *s1, char *s2) 
{
    fprintf(stdout, "printBestAlis\n");

    // on commence tout en bas (for now)
    unsigned int i = mat->h-1; // -1 ?
    unsigned int j = mat->w-1; // same

    char s1_res[strlen(s1)+strlen(s2)]; // résultat avec letters et '-' pour s1
    char s2_res[strlen(s1)+strlen(s2)]; // résultat avec letters et '-' pour s1
    
    for (uint32_t k = 0 ; k < strlen(s1) + strlen(s2) ; k ++) {
        s1_res[k] = '.';
        s2_res[k] = '.';
    }
    
    while (i!=0 || j!=0){

        //un seul chemin, on privilégie la diag (prevs&1) (arbitraire)
        // à changer pour tous les chemins

        struct cell cur = mat->cells[i*mat->w+j];

        printf("Résultat intermédiaire : %s\n", s1_res);
        // màj des chemins selon déplacement optimal trouvé

        switch(cur.prevs){

            //diag
            case 1 :
            case 3 :
            case 5 :
            case 7 :
                printf("Case 1\n");
                s1_res[i+j] = s1[j-1];
                s2_res[i+j] = s2[i-1];
                //s1_res = strcat(s1[j-1], s1_res);
                //s2_res = strcat(s2[i-1], s2_res);
                i--; j--;
                break;

            case 2 :
            case 6 :
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
                fprintf(stderr, "%s\nScore %u\n", "Pas de cas correspondant pour prevs", cur.prevs);
                exit(1);

        }
    }

}

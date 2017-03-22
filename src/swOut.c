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
    
    /*
    for (uint32_t k = 0 ; k < strlen(s1) + strlen(s2) ; k ++) {
        s1_res[k] = '.';
        s2_res[k] = '.';
    }
    */
    
    strcpy(s1_res, s1);
    strcat(s1_res, s2);
    strcpy(s2_res, s1);
    strcat(s2_res, s2);
    
    unsigned int count = strlen(s1_res) - 1; 

    while (i!=0 || j!=0){

        //un seul chemin, on privilégie la diag (prevs&1) (arbitraire)
        // à changer pour tous les chemins

        struct cell cur = mat->cells[i*mat->w+j];

        printf("Résultat intermédiaire : %s\n", s1_res + count + 1);
        printf("Résultat intermédiaire : %s\n\n", s1_res);
        // màj des chemins selon déplacement optimal trouvé

        switch(cur.prevs){

            //diag
            case 1 :
            case 3 :
            case 5 :
            case 7 :
                s1_res[count] = s1[j-1];
                s2_res[count] = s2[i-1];
                i--; j--;
                break;

            case 2 :
            case 6 :
                s1_res[count] = s1[j-1];
                s2_res[count] = '-';
                j--;
                break;

            case 4 :
                s1_res[count] = '-';
                s2_res[count] = s2[i-1];
                i--;
                break;

            default :
                //fprintf(stderr, "%s\nScore %u\n", "Pas de cas correspondant pour prevs", cur.prevs);
                //exit(1);
                i = 0; j = 0;
                break;

        }
        count--;
    }
    //printf("Résultat final : \n\t%s\n\t%s\n", s1_res + count + 1, s2_res + count + 1);
    printf("Résultat final : \n\t%s\n\t%s\n", s1_res + count + 1, s2_res + count + 1);

}

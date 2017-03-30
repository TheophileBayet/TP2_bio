#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>           /* For O_* constants */
#include <sys/stat.h>        /* For mode constants */
#include <semaphore.h>

#include "swOut.h"

void update_res(char* s1, char* s2, char c1, char c2, int count) {
    s1[count] = c1;
    s2[count] = c2;
}

void match(int* i, int* j, char* s1, char* s2, char* s1_res, char* s2_res, int count) {
    if (s1[*j-1] == s2[*i-1]) {
        update_res(s1_res, s2_res, s1[*j-1], s2[*i-1], count);
    } else {
        update_res(s1_res, s2_res, s1[*j-1] + 32, s2[*i-1] + 32, count);
    }
    *i = *i - 1;
    *j = *j - 1;
}

void indel2(int* i, int* j, char* s1, char* s2, char* s1_res, char* s2_res, int count) {
    update_res(s1_res, s2_res, '-', s2[*i-1] + 32, count);
    *i = *i - 1;
}

void indel4(int* i, int* j, char* s1, char* s2, char* s1_res, char* s2_res, int count) {
    update_res(s1_res, s2_res, s1[*j-1] + 32, '-', count);
    *j = *j - 1;
}

void printAllBests(struct matrix *mat,  struct cost *cost, char *s1, char *s2) {

    fprintf(stdout, "printBestAlis\n");

    // search the max score
    int max = 0;
    for (uint32_t k1 = 0; k1<mat->h; k1++){
        for (uint32_t k2 = 0; k2<mat->w; k2++){
            if (mat->cells[mat->w*k1+k2].score>max){
                max=mat->cells[mat->w*k1+k2].score;
            }
        }
    }
    // print the path for all the cells of max score
    for (uint32_t k1 = 0; k1<mat->h; k1++){
        for (uint32_t k2 = 0; k2<mat->w; k2++){
            if (mat->cells[mat->w*k1+k2].score == max){
                printAllPaths(k1, k2, mat, cost, s1, s2);
            }
        }
    }
}

void printAllPaths(int i, int j, struct matrix *mat, struct cost *cost, char *s1, char *s2) {
    char s1_res[strlen(s1)+strlen(s2)]; // résultat avec letters et '-' pour s1
    char s2_res[strlen(s1)+strlen(s2)]; // résultat avec letters et '-' pour s2

    strcpy(s1_res, s1); strcat(s1_res, s2);
    strcpy(s2_res, s1); strcat(s2_res, s2);

    unsigned int count = strlen(s1_res) - 1;

    int ideb=i-1;
    int jdeb=j-1;
    int ifin=0;
    int jfin=0;

    while (i!=0 && j!=0){
        struct cell cur = mat->cells[i*mat->w+j];

        switch(cur.prevs){
            case 1 :
                match(&i, &j, s1, s2, s1_res, s2_res, count);
                break;
            case 2 :
                indel2(&i, &j, s1, s2, s1_res, s2_res, count);
                break;
            case 3 :
                switch(fork()) {
                    case 0 :
                        match(&i, &j, s1, s2, s1_res, s2_res, count);
                        break;
                    default:
                        indel2(&i, &j, s1, s2, s1_res, s2_res, count);
                        break;
                }
                break;
            case 4 :
                indel4(&i, &j, s1, s2, s1_res, s2_res, count);
                break;
            case 5 :
                switch(fork()) {
                    case 0 :
                        match(&i, &j, s1, s2, s1_res, s2_res, count);
                        break;
                    default:
                        indel4(&i, &j, s1, s2, s1_res, s2_res, count);
                        break;
                }
                break;
            case 6 :
                switch(fork()) {
                    case 0 :
                        indel2(&i, &j, s1, s2, s1_res, s2_res, count);
                        break;
                    default:
                        indel4(&i, &j, s1, s2, s1_res, s2_res, count);
                        break;
                }
                break;
            case 7 :
                switch(fork()) {
                    case 0 :
                        match(&i, &j, s1, s2, s1_res, s2_res, count);
                        break;
                    default:
                        switch(fork()) {
                            case 0 :
                                indel2(&i, &j, s1, s2, s1_res, s2_res, count);
                                break;
                            default:
                                indel4(&i, &j, s1, s2, s1_res, s2_res, count);
                                break;
                        }
                        break;
                }
                break;
            default :
                ifin = i ; jfin = j;
                i = 0; j = 0;
                // loop without updating res
                count++;
                break;
        }
        count--;
    }

    sem_t *sem = sem_open("sem", O_CREAT, S_IRWXU, 1);
    sem_wait(sem);

    printf("Best match at s1[%d:%d] and s2[%d:%d]\n", jfin, jdeb, ifin, ideb);
    printf("s1_res : %s\ns2_res : %s\n\n", s1_res + count + 1, s2_res + count + 1);

    sem_post(sem);

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
        switch(cur.prevs){

            case 1 :
            case 3 :
            case 5 :
            case 7 :
                match(&i, &j, s1, s2, s1_res, s2_res, count);
                cur = D->cells[i*D->w+j];
                break;

            case 2 :
            case 6 :
                indel2(&i, &j, s1, s2, s1_res, s2_res, count);
                cur = V->cells[i*V->w+j];
                break;

            case 4 :
                indel4(&i, &j, s1, s2, s1_res, s2_res, count);
                cur = H->cells[i*H->w+j];
                break;

            default :
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

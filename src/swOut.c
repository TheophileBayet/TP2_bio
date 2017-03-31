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
                printAllPaths(k1, k2, max, mat, cost, s1, s2);
            }
        }
    }
}

void printAllPaths(int i, int j, int score, struct matrix *mat, struct cost *cost, char *s1, char *s2) {
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

    printf("Best match at s1[%d:%d] and s2[%d:%d], score : %i\n", jfin, jdeb, ifin, ideb, score);

    for(int i = count + 1 ; i < strlen(s1_res) ; i += 80) {
        printf("s1_res: %.100s\ns2_res: %.100s\n\n", s1_res + i, s2_res + i);
    }
    sem_post(sem);

}


void printAllBestsGotoh(struct matrix *D, struct matrix *V, struct matrix *H,  /*struct cost *cost,*/ char *s1, char *s2) {

    fprintf(stdout, "printBestAlis\n");

    // search the max score
    int max = 0;
    for (uint32_t k1 = 0 ; k1 < D->h; k1++){
        for (uint32_t k2 = 0 ; k2 < D->w ; k2++){
            if (D->cells[D->w*k1+k2].score>max){
                max=D->cells[D->w*k1+k2].score;
            }
        }
    }
    // print the path for all the cells of max score
    for (uint32_t k1 = 0 ; k1 < D->h ; k1++){
        for (uint32_t k2 = 0 ; k2 < D->w ; k2++){
            if (D->cells[D->w*k1+k2].score == max){
                printAllPathsGotoh(k1, k2, max, D, V, H, /*cost,*/ s1, s2);
            }
        }
    }
}


void printAllPathsGotoh(int i, int j, int score, struct matrix *D, struct matrix *V, struct matrix *H, /*struct cost *cost,*/char *s1, char *s2)
{
    char s1_res[strlen(s1)+strlen(s2)]; // résultat avec letters et '-' pour s1
    char s2_res[strlen(s1)+strlen(s2)]; // résultat avec letters et '-' pour s2

    strcpy(s1_res, s1); strcat(s1_res, s2);
    strcpy(s2_res, s1); strcat(s2_res, s2);

    unsigned int count = strlen(s1_res) - 1;

    int ideb=i-1;
    int jdeb=j-1;
    int ifin=0;
    int jfin=0;

    struct cell cur = D->cells[i*D->w+j];

    while (i!=0 && j!=0){

        printf("cur.prevs : %i\n", cur.prevs);

        switch(cur.prevs){
            case 1 :
                match(&i, &j, s1, s2, s1_res, s2_res, count);
                cur = D->cells[i*D->w+j];
                break;
            case 2 :
                indel2(&i, &j, s1, s2, s1_res, s2_res, count);
                cur = V->cells[i*V->w+j];
                break;
            case 3 :
                switch(fork()) {
                    case 0 :
                        match(&i, &j, s1, s2, s1_res, s2_res, count);
                        cur = D->cells[i*D->w+j];
                        break;
                    default:
                        indel2(&i, &j, s1, s2, s1_res, s2_res, count);
                        cur = V->cells[i*V->w+j];
                        break;
                }
                break;
            case 4 :
                indel4(&i, &j, s1, s2, s1_res, s2_res, count);
                cur = H->cells[i*H->w+j];
                break;
            case 5 :
                switch(fork()) {
                    case 0 :
                        match(&i, &j, s1, s2, s1_res, s2_res, count);
                        cur = D->cells[i*D->w+j];
                        break;
                    default:
                        indel4(&i, &j, s1, s2, s1_res, s2_res, count);
                        cur = H->cells[i*H->w+j];
                        break;
                }
                break;
            case 6 :
                switch(fork()) {
                    case 0 :
                        indel2(&i, &j, s1, s2, s1_res, s2_res, count);
                        cur = V->cells[i*V->w+j];
                        break;
                    default:
                        indel4(&i, &j, s1, s2, s1_res, s2_res, count);
                        cur = H->cells[i*H->w+j];
                        break;
                }
                break;
            case 7 :
                switch(fork()) {
                    case 0 :
                        match(&i, &j, s1, s2, s1_res, s2_res, count);
                        cur = D->cells[i*D->w+j];
                        break;
                    default:
                        switch(fork()) {
                            case 0 :
                                indel2(&i, &j, s1, s2, s1_res, s2_res, count);
                                cur = V->cells[i*V->w+j];
                                break;
                            default:
                                indel4(&i, &j, s1, s2, s1_res, s2_res, count);
                                cur = H->cells[i*H->w+j];
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

    printf("Best match at s1[%d:%d] and s2[%d:%d], score : %i\n", jfin, jdeb, ifin, ideb, score);

    for(int i = count + 1 ; i < strlen(s1_res) ; i += 100) {
        printf("s1_res: %.100s\ns2_res: %.100s\n\n", s1_res + i, s2_res + i);
    }
    sem_post(sem);


}

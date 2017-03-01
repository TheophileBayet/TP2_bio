#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mem.h"
#include "swGetSeq.h"


/* allocated memory grows as needed in chunks of BLOCKSIZE bytes */
#define BLOCKSIZE 1000


/* Ask user to input a DNA (if type==0) or protein (if type==1)
   sequence, check the input, if AOK allocate mem for a string
   and store the sequence, changing all lowercase letters to 
   uppercase (eg atgc to ATGC for DNA).
   Return the string, or NULL if problem with input (in which
   case a descriptive message is printed to stderr).
 */
char *getSeq(int type)
{
	if ((type != 0) && (type != 1)) {
		fprintf(stderr, "E: getSeq called with type not 0 or 1, illegal.\n") ;
		exit(1);
	}
	printf("Please enter a ");
	if (type == 0) 
		printf("DNA sequence (atgcATGC");
	else
		printf("protein sequence (arndcqeghilkmfpstwyv upper or lower case");
	printf(", multiple lines ok, ignore spaces, end with blank line)\n") ;
	int prevNL = 0;
	int size = 0;
	int sizeAlloc = BLOCKSIZE*sizeof(char) ;
	char *s = (char*)mallocOrDie(sizeAlloc, "E: OOM in getSeq") ;
	char c ;
	while(1) {
		c = getchar() ;
		if (c=='\n') {
			if (prevNL) {
				if (size==0) {
					fprintf(stderr, "No sequence characters typed! Try again\n") ;
					free(s);
					return(NULL);
				}
				else {
					s[size]='\0' ;
					s = realloc(s,(size+1)*sizeof(char)) ;
					return(s) ;
				}
			}
			else {
				prevNL=1;
				continue;
			}
		}
		else if ((c=='\t')||(c==' '))
			continue;
		else {
			prevNL = 0;
			c = toupper(c) ;
			if ( ((type==0) && ((c=='A')||(c=='T')||(c=='G')||(c=='C'))) ||
			     ((type==1) && (c>='A') && (c<='Z') && (c!='B') && (c!='J') &&
			      (c!='O') && (c!='U') && (c!='X') && (c!='Z')) ) {
				if (size+1 == sizeAlloc) {
					sizeAlloc += BLOCKSIZE*sizeof(char);
					s = realloc(s,sizeAlloc) ;
					if (s==NULL) {
						fprintf(stderr, "E: OOM in getSeq") ;
						exit(1) ;
					}
				}
				s[size++] = c ;
			}
			else {
				fprintf(stderr, "Illegal char typed (%c), must be ", c) ;
				if (type==0)
					fprintf(stderr, "ATGCatgc, try again\n") ;
				else
					fprintf(stderr, "arndcqeghilkmfpstwyv upper or lower case, try again\n") ;
				free(s) ;
				// flush remaining chars from stdin
				while (getchar() != '\n') {
					// NOP
				}

				return(NULL) ;
			}
		}
	}
}


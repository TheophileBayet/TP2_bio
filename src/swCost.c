#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "mem.h"
#include "swCost.h"


/*******************************
  declare static functions
 *******************************/

/* 
   simple DNA substitution cost: +5 for match, -4 for mismatch
*/
static double dnaSubst(char x,char y) ;


/*
  Return the score of the x<->y protein substitution according to 
  the BLOSUM62 matrix.
  Precondition: x and y are uppercase letters representing amino acids.
*/
static double protSubstBlosum62(char x,char y) ;

/* 
   set the cost for the x<->y protein substitution in mat. 
   The convention is that this cost is stored at 
   mat[x-'A' +  26*(y-'A')] , where
  x and y must be uppercase letters representing amino acids.
 */
static void fillCost(double *mat, char x, char y, double cost) ;

/*
  Fill mat (a pre-allocated 26*26 matrix) with the BLOSUM62 substitution
  costs.
*/
static void fillBlosum62(double *mat) ;


/*******************************
  define static functions
 *******************************/


static double dnaSubst(char x,char y)
{
	if (x==y)
		return(5);
	else
		return(-4);
}



static double protSubstBlosum62(char x,char y)
{
	/* use a static blosum62 matrix that we fill once, so
	   later we just need to lookup the subst cost in the matrix */
	static double blosum62[26*26] = {DBL_MAX} ;
	if (blosum62[0] == DBL_MAX) {
		/* this is the first call to this function, blosum62 has not 
		   been initialized yet, do it */
		fillBlosum62(blosum62) ;
	}
	if ((x<'A') || (x>'Z') || (y<'A') || (y>'Z')) {
		fprintf(stderr, "E: protSubstBlosum62 called with %c and %c but both args must be UPPERCASE letters.\n", x, y);
		exit(1);
	}
	return(blosum62[x-'A' + 26*(y-'A')]) ;
}


static void fillCost(double *mat, char x, char y, double cost) {
	mat[x-'A' +  26*(y-'A')] = cost;
	if (x != y) 
		mat[y-'A' +  26*(x-'A')] = cost;
}

static void fillBlosum62(double *mat)
{
	/* initialize everything to DBL_MAX */
	for (int i = 0; i < 26; i++)
		for (int j = 0; j < 26; j++)
			mat[i+26*j] = DBL_MAX;
	
	/* now fill costs for actual amino acids */
	fillCost(mat,'C','C',9);
	fillCost(mat,'C','S',-1);
	fillCost(mat,'C','T',-1);
	fillCost(mat,'C','P',-3);
	fillCost(mat,'C','A',0);
	fillCost(mat,'C','G',-3);
	fillCost(mat,'C','N',-3);
	fillCost(mat,'C','D',-3);
	fillCost(mat,'C','E',-4);
	fillCost(mat,'C','Q',-3);
	fillCost(mat,'C','H',-3);
	fillCost(mat,'C','R',-3);
	fillCost(mat,'C','K',-3);
	fillCost(mat,'C','M',-1);
	fillCost(mat,'C','I',-1);
	fillCost(mat,'C','L',-1);
	fillCost(mat,'C','V',-1);
	fillCost(mat,'C','F',-2);
	fillCost(mat,'C','Y',-2);
	fillCost(mat,'C','W',-2);

	fillCost(mat,'S','S',4);
	fillCost(mat,'S','T',1);
	fillCost(mat,'S','P',-1);
	fillCost(mat,'S','A',1);
	fillCost(mat,'S','G',0);
	fillCost(mat,'S','N',1);
	fillCost(mat,'S','D',0);
	fillCost(mat,'S','E',0);
	fillCost(mat,'S','Q',0);
	fillCost(mat,'S','H',-1);
	fillCost(mat,'S','R',-1);
	fillCost(mat,'S','K',0);
	fillCost(mat,'S','M',-1);
	fillCost(mat,'S','I',-2);
	fillCost(mat,'S','L',-2);
	fillCost(mat,'S','V',-2);
	fillCost(mat,'S','F',-2);
	fillCost(mat,'S','Y',-2);
	fillCost(mat,'S','W',-3);

	fillCost(mat,'T','T',5);
	fillCost(mat,'T','P',-1);
	fillCost(mat,'T','A',0);
	fillCost(mat,'T','G',-2);
	fillCost(mat,'T','N',0);
	fillCost(mat,'T','D',-1);
	fillCost(mat,'T','E',-1);
	fillCost(mat,'T','Q',-1);
	fillCost(mat,'T','H',-2);
	fillCost(mat,'T','R',-1);
	fillCost(mat,'T','K',-1);
	fillCost(mat,'T','M',-1);
	fillCost(mat,'T','I',-1);
	fillCost(mat,'T','L',-1);
	fillCost(mat,'T','V',0);
	fillCost(mat,'T','F',-2);
	fillCost(mat,'T','Y',-2);
	fillCost(mat,'T','W',-2);

	fillCost(mat,'P','P',7);
	fillCost(mat,'P','A',-1);
	fillCost(mat,'P','G',-2);
	fillCost(mat,'P','N',-2);
	fillCost(mat,'P','D',-1);
	fillCost(mat,'P','E',-1);
	fillCost(mat,'P','Q',-1);
	fillCost(mat,'P','H',-2);
	fillCost(mat,'P','R',-2);
	fillCost(mat,'P','K',-1);
	fillCost(mat,'P','M',-2);
	fillCost(mat,'P','I',-3);
	fillCost(mat,'P','L',-3);
	fillCost(mat,'P','V',-2);
	fillCost(mat,'P','F',-4);
	fillCost(mat,'P','Y',-3);
	fillCost(mat,'P','W',-4);

	fillCost(mat,'A','A',4);
	fillCost(mat,'A','G',0);
	fillCost(mat,'A','N',-2);
	fillCost(mat,'A','D',-2);
	fillCost(mat,'A','E',-1);
	fillCost(mat,'A','Q',-1);
	fillCost(mat,'A','H',-2);
	fillCost(mat,'A','R',-1);
	fillCost(mat,'A','K',-1);
	fillCost(mat,'A','M',-1);
	fillCost(mat,'A','I',-1);
	fillCost(mat,'A','L',-1);
	fillCost(mat,'A','V',0);
	fillCost(mat,'A','F',-2);
	fillCost(mat,'A','Y',-2);
	fillCost(mat,'A','W',-3);

	fillCost(mat,'G','G',6);
	fillCost(mat,'G','N',0);
	fillCost(mat,'G','D',-1);
	fillCost(mat,'G','E',-2);
	fillCost(mat,'G','Q',-2);
	fillCost(mat,'G','H',-2);
	fillCost(mat,'G','R',-2);
	fillCost(mat,'G','K',-2);
	fillCost(mat,'G','M',-3);
	fillCost(mat,'G','I',-4);
	fillCost(mat,'G','L',-4);
	fillCost(mat,'G','V',-3);
	fillCost(mat,'G','F',-3);
	fillCost(mat,'G','Y',-3);
	fillCost(mat,'G','W',-2);

	fillCost(mat,'N','N',6);
	fillCost(mat,'N','D',1);
	fillCost(mat,'N','E',0);
	fillCost(mat,'N','Q',0);
	fillCost(mat,'N','H',1);
	fillCost(mat,'N','R',0);
	fillCost(mat,'N','K',0);
	fillCost(mat,'N','M',-2);
	fillCost(mat,'N','I',-3);
	fillCost(mat,'N','L',-3);
	fillCost(mat,'N','V',-3);
	fillCost(mat,'N','F',-3);
	fillCost(mat,'N','Y',-2);
	fillCost(mat,'N','W',-4);

	fillCost(mat,'D','D',6);
	fillCost(mat,'D','E',2);
	fillCost(mat,'D','Q',0);
	fillCost(mat,'D','H',-1);
	fillCost(mat,'D','R',-2);
	fillCost(mat,'D','K',-1);
	fillCost(mat,'D','M',-3);
	fillCost(mat,'D','I',-3);
	fillCost(mat,'D','L',-4);
	fillCost(mat,'D','V',-3);
	fillCost(mat,'D','F',-3);
	fillCost(mat,'D','Y',-3);
	fillCost(mat,'D','W',-4);

	fillCost(mat,'E','E',5);
	fillCost(mat,'E','Q',2);
	fillCost(mat,'E','H',0);
	fillCost(mat,'E','R',0);
	fillCost(mat,'E','K',1);
	fillCost(mat,'E','M',-2);
	fillCost(mat,'E','I',-3);
	fillCost(mat,'E','L',-3);
	fillCost(mat,'E','V',-2);
	fillCost(mat,'E','F',-3);
	fillCost(mat,'E','Y',-2);
	fillCost(mat,'E','W',-3);

	fillCost(mat,'Q','Q',5);
	fillCost(mat,'Q','H',0);
	fillCost(mat,'Q','R',1);
	fillCost(mat,'Q','K',1);
	fillCost(mat,'Q','M',0);
	fillCost(mat,'Q','I',-3);
	fillCost(mat,'Q','L',-3);
	fillCost(mat,'Q','V',-2);
	fillCost(mat,'Q','F',-3);
	fillCost(mat,'Q','Y',-1);
	fillCost(mat,'Q','W',-2);

	fillCost(mat,'H','H',8);
	fillCost(mat,'H','R',0);
	fillCost(mat,'H','K',-1);
	fillCost(mat,'H','M',-2);
	fillCost(mat,'H','I',-3);
	fillCost(mat,'H','L',-3);
	fillCost(mat,'H','V',-3);
	fillCost(mat,'H','F',-1);
	fillCost(mat,'H','Y',2);
	fillCost(mat,'H','W',-2);

	fillCost(mat,'R','R',5);
	fillCost(mat,'R','K',2);
	fillCost(mat,'R','M',-1);
	fillCost(mat,'R','I',-3);
	fillCost(mat,'R','L',-2);
	fillCost(mat,'R','V',-3);
	fillCost(mat,'R','F',-3);
	fillCost(mat,'R','Y',-2);
	fillCost(mat,'R','W',-3);

	fillCost(mat,'K','K',5);
	fillCost(mat,'K','M',-1);
	fillCost(mat,'K','I',-3);
	fillCost(mat,'K','L',-2);
	fillCost(mat,'K','V',-2);
	fillCost(mat,'K','F',-3);
	fillCost(mat,'K','Y',-2);
	fillCost(mat,'K','W',-3);

	fillCost(mat,'M','M',5);
	fillCost(mat,'M','I',1);
	fillCost(mat,'M','L',2);
	fillCost(mat,'M','V',1);
	fillCost(mat,'M','F',0);
	fillCost(mat,'M','Y',-1);
	fillCost(mat,'M','W',-1);

	fillCost(mat,'I','I',4);
	fillCost(mat,'I','L',2);
	fillCost(mat,'I','V',3);
	fillCost(mat,'I','F',0);
	fillCost(mat,'I','Y',-1);
	fillCost(mat,'I','W',-3);

	fillCost(mat,'L','L',4);
	fillCost(mat,'L','V',1);
	fillCost(mat,'L','F',0);
	fillCost(mat,'L','Y',-1);
	fillCost(mat,'L','W',-2);

	fillCost(mat,'V','V',4);
	fillCost(mat,'V','F',-1);
	fillCost(mat,'V','Y',-1);
	fillCost(mat,'V','W',-3);

	fillCost(mat,'F','F',6);
	fillCost(mat,'F','Y',3);
	fillCost(mat,'F','W',1);

	fillCost(mat,'Y','Y',7);
	fillCost(mat,'Y','W',2);

	fillCost(mat,'W','W',11);
}



/*******************************
  define extern functions
 *******************************/

struct cost *costDna(double indelOpen, double indelExtend)
{
	if ((indelOpen > indelExtend) || (indelExtend >= 0)) {
		fprintf(stderr, "E: in costDna, must have indelOpen <= indelExtend < 0\n");
		exit(1);
	}
	struct cost *cost = mallocOrDie(sizeof(struct cost), "E: OOM for cost\n") ;
	cost->indelOpen = indelOpen;
	cost->indelExtend = indelExtend ;
	cost->subst = &dnaSubst ;
	return(cost) ;
}

struct cost *costProt(double indelOpen, double indelExtend)
{
	if ((indelOpen > indelExtend) || (indelExtend >= 0)) {
		fprintf(stderr, "E: in costProt, must have indelOpen <= indelExtend < 0\n");
		exit(1);
	}
	struct cost *cost = mallocOrDie(sizeof(struct cost), "E: OOM for cost\n") ;
	cost->indelOpen = indelOpen;
	cost->indelExtend = indelExtend ;
	cost->subst = &protSubstBlosum62 ;
	return(cost) ;
}

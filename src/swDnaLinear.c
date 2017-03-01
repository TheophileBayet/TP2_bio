#include <stdio.h>
#include <stdlib.h>

#include "swCost.h"
#include "swCalc.h"
#include "swGetSeq.h"
#include "swOut.h"

int main(void)
{
	char *s1 ;
	while((s1 = getSeq(0)) == NULL) {
		// nothing to do
	}
	char *s2 ;
	while((s2 = getSeq(0)) == NULL) {
	}
	printf("Sequences read:\ns1\t%s\ns2\t%s\n\n", s1, s2) ;

	struct cost *cost = costDna(-10,-10);
	struct matrix *mat = swInitMat(s1,s2);
	swFillMat(mat,cost,s1,s2);
	/* for debugging you can uncomment:
	   swPrintMat(mat); */
	printBestAlis(mat,cost,s1,s2);

	swFreeMat(mat);
	free(cost);
	free(s1);
	free(s2);
	return(0);
}

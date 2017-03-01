#ifndef _SWCOST_H_
#define _SWCOST_H_


struct cost {
	/* if your implementation only deals with the linear cost model,
	   use indelOpen as the fixed indel cost and ignore indelExtend */
	double indelOpen;
	double indelExtend;
	/* subst: pointer to function returning the substitution score of
	   its args. Of course the function will be different for DNA and proteins.
	   Precondition: args are both uppercase letters */
	double (*subst)(char,char);
};


/* return a cost structure for DNA alignment.
   This cost structure is designed for an affine gap cost model (Gotoh),
   but implementations of the linear model (Smith-Waterman) should also use
   it and simply ignore indelExtend.
   Currently the subst function uses match==+5 and mismatch==-4, this
   can be easily customized by editing the static dnaSubst function.
   Requires indelOpen <= indelExtend < 0.
*/
struct cost *costDna(double indelOpen, double indelExtend);


/* return a cost structure for protein alignment.
   Again the struct is designed for Gotoh but should also be used for
   Smith-Waterman (just ignore indelExtend).
   The returned structure currently implements a BLOSUM62 substitution 
   cost.
   Requires indelOpen <= indelExtend < 0.
*/
struct cost *costProt(double indelOpen, double indelExtend);

#endif

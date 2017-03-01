#ifndef _SWGETSEQ_H_
#define _SWGETSEQ_H_


/* Ask user to input a DNA (if type==0) or protein (if type==1)
   sequence, check the input, if AOK allocate mem for a string
   and store the sequence, changing all lowercase letters to 
   uppercase (eg atgc to ATGC for DNA).
   Return the string, or NULL if problem with input (in which
   case a descriptive message is printed to stderr).
 */
char *getSeq(int type) ;

#endif

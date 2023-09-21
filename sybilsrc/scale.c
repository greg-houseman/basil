/*--------------------------------------------------------------------
 *    Basil / Sybil:   scale.c   1.1  1 October 1998
 *
 *    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/

/*
 * based on code in Steppe, T. "Composing Well-tempered Linear Scales",
 * Computer Language, Sept.,49,1989
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

static double pdSet[] = {1.0,2.0,5.0,10.0};
#define SET_LEN  (sizeof(pdSet)/sizeof(double)-1)

/* 
 * function declarations
 */
void scLewart();
void scCalcExtLabel();
double scFirstNiceNum();
double scNextNiceNum();
static double scPower();


/*
 * Lewart's algorithm
 */
void scLewart( min,max,approx_intrvls,scalemin,scalemax,actual )
double min, max, *scalemin, *scalemax;
int approx_intrvls, *actual;
{
    double intrvl_size, power, nicenum;
    int index, lomult, himult;

    assert (min < max);
    assert (approx_intrvls >= 2);

    /* calculate smallest potential interval */
    intrvl_size = (max - min)/approx_intrvls;

    /* find nice number closest to the smallest potential
     * interval. Use the geometric means of adjacent
     * multiplier values as break points
     */
    for (nicenum=scFirstNiceNum(intrvl_size,&index,&power);
            index < (SET_LEN-1) &&
            sqrt(pdSet[index]*pdSet[index+1])*power < intrvl_size;
                nicenum=scNextNiceNum(pdSet,SET_LEN,&index,&power)) ;

    /* produce the scale using the nice number */
    scCalcExtLabel(min,max,nicenum,&lomult,&himult);

    /* calculate scale limits */
    *scalemin = lomult * nicenum;
    *scalemax = himult * nicenum;
    *actual = himult - lomult;
}

void scCalcExtLabel(min,max,nicenum,lomult,himult)
double min, max, nicenum;
int *lomult, *himult;
{
    *lomult = (int)floor(min/nicenum);
    if ((double)(*lomult+1)*nicenum <= min) (*lomult)++;

    *himult = (int)ceil(max/nicenum);
    if ((double)(*lomult-1)*nicenum >= max) (*himult)--;
}

double scFirstNiceNum(size,index,power)
double size,*power;
int *index;
{
    int exp;

    /* calc initial power of 10 */
    exp = (int)floor(log10(size));
    *power = scPower(10.0,exp);
    if (*power * 10.0 <= size) *power *= 10.0;

    /* initial index is always zero */
    *index = 0;
    return(*power);
}

double scNextNiceNum(pdSet,num,index,power)
double *pdSet, *power;
int num, *index;
{
    (*index)++;
    
    /* if max index has been exceeded, reset index to zero
     * and increase the power of 10
     */
    if (*index > num) {
        *index = 0;
        *power *= 10.0;
    }
    return(pdSet[*index] * *power);
}

static double scPower( root, exponent )
double root;
int exponent;
{
    double res;

    /* for negative exponents, invert root and use positive */
    if (exponent < 0) {
        root = 1.0/root;
        exponent = -exponent;
    }
    res = 1.0;
    while (exponent) {
        if (exponent & 1) res *= root;
        exponent >>= 1;
        if (exponent) root *= root;
    }
    return( res );
}

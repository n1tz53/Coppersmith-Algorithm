#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"

bi_poly * init_poly(int size)
{
    bi_poly * ret = (bi_poly *) malloc(sizeof(bi_poly));
    ret->deg = 0;
    ret->sz = size;
    ret->coeff = (ui *) calloc(size, sizeof(ui));

    if (ret->coeff == NULL)
    {
        printf("Polynomial can not be created !!!");
    }

    return ret;
}

bi_poly * create_poly(int size, int * ary, int len)
{
    int i;
    bi_poly * ret = init_poly(size);

    for (i = 0; i < len; i++)
        flip_coeff(ret, ary[i]);

    ret->deg = ary[len - 1];

    return ret;
}

bi_poly * copy_poly(bi_poly * bp)
{
    int i = 0;
    bi_poly * ret = init_poly(bp->sz);
    ret->deg = bp->deg;

    for (i = 0; i < bp->sz; i++)
    {
        //printf("%d %d\n", i, bp->coeff[i]);
        ret->coeff[i] = bp->coeff[i];
    }

    return ret;
}

void flip_coeff(bi_poly * bp, int n)
{
    bp->coeff[n >> 5] ^= bits[n & 31];
}

void update_degree(bi_poly * bp)
{
    int i, j;
    bp->deg = 0;

    for (i = bp->sz - 1; i >= 0; i--)
        if ((bp->coeff[i] & 4294967295u)) break;

    if(i >= 0)
    {
        for (j = 31; j >= 0; j--)
        {
            if ((bp->coeff[i] & bits[j]) > 0)
            {
                bp->deg = i * 32 + j;
                break;
            }
        }

        j = bp->deg + 1;
        j = (j & 31) ? (j >> 5) + 1 : (j >> 5);
    }
    else j = 1;

    if (j < bp->sz)
    {
        bp->sz = j;
        bp->coeff = (ui *) realloc(bp->coeff, bp->sz * sizeof(ui));
    }
}

void show_poly(bi_poly * bp)
{
    int put = 0, i, j;

    for (i = bp->sz - 1; i >= 0; i--)
    {

        for (j = 31; j >= 0; j--)
        {
            if ((bp->coeff[i] & bits[j]) > 0)
            {
                if (put)
                    printf(" + ");

                if (j != 0 || i != 0)
                {
                    printf("x^");
                    printf("%d", i * 32 + j);
                }
                else printf("1");
                put = 1;
            }
        }
    }

    printf(" %d %d ", bp->deg, bp->sz);
    printf("\n");
}

void * free_poly(bi_poly * bp)
{
    free(bp->coeff);
    free(bp);
}




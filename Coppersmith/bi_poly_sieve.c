#include <stdio.h>
#include <stdlib.h>
#include "common.h"

inline int max(int a, int b) { return (a >= b) ? a : b; }

int get_lowest_bit(int n)
{
    int i;
    for (i = 0; i < 32; i++)
    {
        if (n & bits[i])
            return i;
    }

    return 0;
}

void gen_irred()
{
    int i, j;
    FILE * fp = fopen("primes.txt", "w");
    bi_poly ** elt = (bi_poly **) malloc((1 << 13) * sizeof(bi_poly *));

    for (i = 0; i < (1 << 13); i++)
    {
        elt[i] = init_poly(1);

        for (j = 0; j < 13; j++)
        {
            if ((i & (1 << j)))
                flip_coeff(elt[i], j);
        }

        update_degree(elt[i]);
    }

    char * is_prime = (char *) malloc((1 << 13) * sizeof(char));
    is_prime[3] = is_prime[2] = 1;
    bi_poly * tmp;

    for (i = 4; i < (1 << 13); i++)
    {
        is_prime[i] = 1;

        for (j = i - 1; j >= 2; j--)
        {
            tmp = copy_poly(elt[i]);
            reduce2(tmp, elt[j]);

            if (tmp->deg || (tmp->coeff[0] & 1));
            else
            {
                is_prime[i] = 0;
                break;
            }

            free_poly(tmp);
        }
    }

    for (i = 2; i < (1 << 13); i++)
    {
        if (is_prime[i])
        {
            fprintf(fp, "%d %d\n", elt[i]->deg, elt[i]->coeff[0]);
        }
    }

    fclose(fp);
}

void do_sieve(bi_poly * u1x, char * sieve, bi_poly ** irred, int * index)
{
    int i, j, d, dim;
    bi_poly * u2x, * g, * tmp;

    for (d = 1; d < 13; d++)
    {
        dim = max(11 - d, 0);

        for (j = index[d - 1]; j < index[d]; j++)
        {
            g = irred[j];
            u2x = shift_left(u1x, 32);
            reduce2(u2x, g);

            if (u2x->deg < 11)
            {
                for (i = 1; i <= (1 << dim); i++)
                {
                    sieve[u2x->coeff[0]] += d;
                    tmp = copy_poly(g);
                    tmp->coeff[0] <<= get_lowest_bit(i);
                    u2x->coeff[0] ^= tmp->coeff[0];
                    free_poly(tmp);
                }
            }

            free_poly(u2x);
        }
    }
}



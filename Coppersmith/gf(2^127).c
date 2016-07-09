#include <stdio.h>
#include <time.h>
#include "common.h"

void write_to_file(FILE * fp, bi_poly * bp)
{
    int i;
    fprintf(fp, "%x %x", bp->deg, bp->sz);

    for (i = 0; i < bp->sz; i++)
        fprintf(fp, " %x", bp->coeff[i]);
}

void coppersmith_gf2_127(bfa * gf, bi_poly ** elt, int len)
{
    int i, j;
    bi_poly * w1_x, * w2_x, * d, * tmp, * sw1_x;
    FILE * fp = fopen("coppersmith_pair_127.txt", "w");

    clock_t begin = clock(), end;
    int ctr = 0;

    for (i = 1; i < len; i++)
    {
        sw1_x = shift_left(elt[i], 32);

        for (j = 1; j < len; j++)
        {
            if (i != 1 && (i == j))
                continue;

            d = (i <= j) ? gcd(elt[i], elt[j]) : gcd(elt[j], elt[i]);

            if (d->deg == 0 && (d->coeff[0] & 1))
            {
                w1_x = add(sw1_x, elt[j]);

                if (smooth(gf, w1_x))
                {
                    tmp = sqr(gf, w1_x);
                    w2_x = sqr(gf, tmp);
                    free_poly(tmp);
                    reduce(gf, w2_x);

                    if (smooth(gf, w2_x))
                    {
                        ++ctr;
                        write_to_file(fp, w1_x);
                        write_to_file(fp, w2_x);
                    }

                    free_poly(w2_x);
                }

                free_poly(w1_x);
            }

            free_poly(d);
        }

        free_poly(sw1_x);
    }

    fclose(fp);
    end = clock() - begin;
    int msec = end * 1000 / CLOCKS_PER_SEC;

    printf("number of smooth (w1x, w2x) pair = %d time take = %d\n", ctr, msec / 1000);
}

void perform_sieve_gf2_127(bfa * gf, bi_poly ** elt, int len)
{
    int i, j, deg = 1, ctr = 0, k;
    char buffer[10], * rem;
    bi_poly * w1x, * w2x, * sw1x, * tmp, * u2x, * d;
    bi_poly ** irred = (bi_poly **) malloc(747 * sizeof(bi_poly *));
    int * index = (int *) malloc(13 * sizeof(int));
    index[0] = 0;
    FILE * fp = fopen("primes.txt", "r");

    /* reading irreducible from file */

    for (i = 0; i < 747; i++)
    {
        irred[i] = init_poly(1);
        fgets(buffer, sizeof(buffer), fp);
        irred[i]->deg = (int) strtoul(buffer, &rem, 10);
        ++rem;
        irred[i]->coeff[0] = (u_int32) atoi(rem);

        if (irred[i]->deg != deg)
        {
            index[deg] = i;
            deg = irred[i]->deg;
        }
    }

    index[12] = 746;

    fclose(fp);

    /* allocating and initializing sieve array */

    char * sieve = (char *) calloc((1 << 11), sizeof(char));
    fp = fopen("sieve_pair_127.txt", "w");

    /* perform sieving for each u1x */

    for (i = 1; i < len; i++)
    {
        do_sieve(elt[i], sieve, irred, index);
        sw1x = shift_left(elt[i], 32);

        for (j = 1; j < len; j++)
        {
            if (sieve[j] > elt[i]->deg + 20)
            {
                u2x = init_poly(1);
                u2x->coeff[0] = j;
                update_degree(u2x);
                w1x = copy_poly(sw1x);
                w1x->coeff[0] ^= j;
                update_degree(w1x);

                d = (elt[i]->deg <= u2x->deg) ? gcd(elt[i], u2x) : gcd(u2x, elt[i]);

                if (d->deg) continue;

                if (smooth(gf, w1x))
                {
                    w2x = sqr(gf, w1x);
                    tmp = w2x;
                    w2x = sqr(gf, tmp);
                    free_poly(tmp);
                    reduce(gf, w2x);

                    if (smooth(gf, w2x))
                    {
                        ++ctr;
                        write_to_file(fp, w1x);
                        write_to_file(fp, w2x);
                    }

                    free_poly(w2x);
                }

                free_poly(w1x);
            }
        }

        free_poly(sw1x);
        memset(sieve, 0, (1u << 11) * sizeof(char));
    }

    fclose(fp);

    printf("number of smooth pair found = %d\n", ctr);
}

void collect_relation_gf2_127()
{
    int i, j;
    bi_poly ** elt;
  /* initialize field GF(2^127) */
    int a[] = {0, 1, 127};
    bfa * gf = init_bfa(a, 3);

  /* generate all the binary polynomials of degree atmost 10 */
    elt = (bi_poly **) malloc((1 << 11) * sizeof(bi_poly *));

    for (i = 0; i < (1u << 11); i++)
    {
        elt[i] = init_poly(1);

        for (j = 0; j < 11; j++)
        {
            if ((i & (1u << j))) flip_coeff(elt[i], j);
        }

        update_degree(elt[i]);
    }

   /* Choose smoothness bound and initialize */

    smooth_parameter(gf, 12);

    /* Collect relation using Coppersmith method */

    //coppersmith_gf2_127(gf, elt, (1 << 11));        /* uncomment for using coppersmith approach */

    /* Collect relations using Gordon and McCurley Sieve */

    perform_sieve_gf2_127(gf, elt, (1 << 11));        /* Gordon and McCurley approach */

    /* free the memory used */

    for (i = 0; i < (1 << 11); i++)
        free_poly(elt[i]);

    free(elt);
    free_bfa(gf);
}

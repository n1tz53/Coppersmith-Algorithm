#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"

void test_binary_polynomial()
{
    int i, j;
    bi_poly ** elt = (bi_poly **) malloc((1u << 11) * sizeof(bi_poly *));

    for (i = 0; i < (1u << 11); i++)
    {
        elt[i] = init_poly(1);

        for (j = 0; j < 11; j++)
        {
            if ((i & (1u << j)) > 0)
                flip_coeff(elt[i], j);
        }

        update_degree(elt[i]);
    }

    int a[3] = {0, 1, 127};
    ffa* gf = init_ffa(a, 3);

    bi_poly * w1_x, * w2_x, * d, * tmp, * sw1_x;
    int b[1] = {128};
    w1_x = create_poly(5, b, 1);

    clock_t begin = clock(), end;
    int ctr = 0;

    for (i = 1; i < (1 << 11); i++)
    {
        sw1_x = shift_left(elt[i], 32);

        for (j = 1; j < (1 << 11); j++)
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
                        printf("%d\n", ++ctr);
                    }

                    free_poly(w2_x);
                }

                free_poly(w1_x);
            }

            free_poly(d);
        }

        free_poly(sw1_x);
    }

    for (i = 0; i < (1 << 11); i++)
        free_poly(elt[i]);

    free(elt);

    free_ffa(gf);

    end = clock() - begin;
    int msec = end * 1000 / CLOCKS_PER_SEC;

    printf("total relatively prime pairs = %d time = %d\n", ctr, msec / 1000);

}

void test_sieve()
{
    int i, j, deg = 1;
    char buffer[10];
    char * rem;
    bi_poly ** irred = (bi_poly **) malloc(747 * sizeof(bi_poly *));
    int * index = (int *) malloc(13 * sizeof(int));
    index[0] = 0;
    FILE * fp = fopen("primes.txt", "r");

    for (i = 0; i < 747; i++)
    {
        irred[i] = init_poly(1);
        fgets(buffer, sizeof(buffer), fp);
        irred[i]->deg = (int) strtoul(buffer, &rem, 10);
        ++rem;
        irred[i]->coeff[0] = (ui) atoi(rem);

        if (irred[i]->deg != deg)
        {
            index[deg] = i;
            deg = irred[i]->deg;
        }
    }

    index[12] = 746;

    fclose(fp);

    char * sieve = (char *) calloc((1 << 11), sizeof(char));
    bi_poly ** elt = (bi_poly **) malloc((1u << 11) * sizeof(bi_poly *));

    for (i = 0; i < (1u << 11); i++)
    {
        elt[i] = init_poly(1);

        for (j = 0; j < 11; j++)
        {
            if ((i & (1u << j)) > 0)
                flip_coeff(elt[i], j);
        }

        update_degree(elt[i]);
    }

    fp = fopen("poly_to_factor.txt", "w");

    int a[3] = {0, 1, 127};
    ffa* gf = init_ffa(a, 3);

    int avg;
    int ctr = 0, k;
    bi_poly * w1x, * w2x, * sw1x, * tmp, * u2x, * d;

    for (i = 1; i < (1u << 11); i++)
    {
        do_sieve(elt[i], sieve, irred, index);
        sw1x = shift_left(elt[i], 32);

        for (j = 1; j < (1 << 11); j++)
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

                if (d->deg)
                    continue;

                if (smooth(gf, w1x))
                {
                    w2x = sqr(gf, w1x);
                    tmp = w2x;
                    w2x = sqr(gf, tmp);
                    free_poly(tmp);
                    reduce(gf, w2x);

                    if (smooth(gf, w2x))
                    {
                        fprintf(fp, "%d %d", w1x->deg, w1x->sz);
                        for (k = 0; k < w1x->sz; k++)
                            fprintf(fp, " %d", w1x->coeff[k]);
                        fprintf(fp, "\n");
                        fprintf(fp, "%d %d", w2x->deg, w2x->sz);
                        for (k = 0; k < w2x->sz; k++)
                            fprintf(fp, " %d", w2x->coeff[k]);
                        fprintf(fp, "\n");
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

    for (i = 0; i < 747; i++)
        free_poly(irred[i]);
    free(irred);
}

void test_factorization()
{
    int i, deg, sz;
    char buffer[1024];
    char * rem, * nxt;
    FILE * fp = fopen("poly_to_factor.txt", "r");
    bi_poly * px, * qx, * tmp, * tmp2, * rx;

    int a[3] = {0, 1, 127};
    int b[7] = {0, 1, 2, 3, 4, 5, 6};
    ffa* gf = init_ffa(a, 3);


    while (fgets(buffer, 1024, fp) != NULL)
    {
        //fputs(buffer, stdout);
        deg = (int) strtoul(buffer, &rem, 10);
        ++rem;
        sz = (int) strtoul(rem, &nxt,10);

        px = init_poly(sz);
        px->deg = deg;

        for (i = 0; i < px->sz; i++)
        {
            rem = nxt;
            ++rem;
            px->coeff[i] = (ui) strtoul(rem, &nxt, 10);
        }

        //px = create_poly(1, b, 7);
        show_poly(px);

        factor_list * ft = sff(px);
        factor * itr;
        qx = init_poly(1);
        flip_coeff(qx, 0);

        printf("Factors are : \n");

        for (itr = ft->head; itr != NULL; itr = itr->nxt)
        {
            if (itr->divisor->deg)
            {
                //show_poly(itr->divisor);
                //printf("%d\n", itr->idx);
                tmp = raise(gf, itr->divisor, itr->idx);
                tmp2 = qx;
                qx = multiply(tmp2, tmp);
                free_poly(tmp2);
                free_poly(tmp);

                factor_list * dft = ddf(gf, itr->divisor);
                rx = init_poly(1);
                flip_coeff(rx, 0);

                factor * itr2;

                for (itr2 = dft->head; itr2 != NULL; itr2 = itr2->nxt)
                {
                    tmp = rx;
                    rx = multiply(tmp, itr2->divisor);
                    free(tmp);
                }

                free_factor_list(dft);
                free(itr2);

                if (! same(rx, itr->divisor))
                {
                    printf("Distinct Degree Factorization Failed !!!\n");
                    show_poly(itr->divisor);
                    show_poly(rx);
                }

                free_poly(rx);
            }
        }


        if (!same(px, qx))
        {
            printf("Square Free Factorization Failed !!!\n");
            show_poly(px);
            show_poly(qx);
        }

        free_poly(px);
        free_poly(qx);
        free_factor_list(ft);

        //break;

    }

    fclose(fp);

    free_ffa(gf);
}

int main()
{
    printf("Checking Square Free Factorization : \n");
    test_factorization();

    return 0;
}

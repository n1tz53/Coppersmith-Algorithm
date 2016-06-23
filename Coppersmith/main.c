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

 /*       show_poly(elt[i]);

        if (elt[i]->deg > 10)
        {
            printf("made some mistake !!! %d", i);
        }
*/

    }

    int a[3] = {0, 1, 127};
    ffa* gf = init_ffa(a, 3);

 /*   show_poly(gf->fx);

    show_poly(gf->rx);

    for (i = 0; i < 32; i++)
        show_poly(gf->ux[i]);

    for (i = 1; i < 13; i++)
        show_poly(gf->factors[i]);
*/
    bi_poly * w1_x, * w2_x, * d, * tmp, * sw1_x;

    clock_t begin = clock(), end;
    int ctr = 0;

  /*  for (i = 0; i < (1 << 11); i++)
    {
        //show_poly(elt[i]);
        tmp = shift_left(elt[i], 32);
        //show_poly(tmp);
        w1_x = sqr(gf, tmp);
        //show_poly(w1_x);
        free_poly(tmp);
        tmp = sqr(gf, w1_x);
        //show_poly(tmp);
        free_poly(w1_x);
        d = copy_poly(tmp);
        reduce(gf, tmp);
        //show_poly(tmp);
        w1_x = reduce2(d, gf->fx);
        //show_poly(w1_x);

        if (! same(tmp, w1_x))
        {
            printf("Error !!!\n");
            //show_poly(elt[i]);
            //show_poly(tmp);
            //show_poly(w1_x);
            //int hey = 1 / 0;
        }

        free_poly(tmp);
        free_poly(d);
        free_poly(w1_x);
    } */

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
                    //printf("%d\n", ++ctr);
                    //show_poly(w1_x);
                    tmp = sqr(gf, w1_x);
                    w2_x = sqr(gf, tmp);
                    free_poly(tmp);
                    reduce(gf, w2_x);
                    //show_poly(w2_x);

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

int main()
{
    printf("Number of Smooth Pairs : \n");
    test_binary_polynomial();
    return 0;
}

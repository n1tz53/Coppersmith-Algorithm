#include <stdlib.h>
#include "common.h"

factor_list * sff(bi_poly * p)
{
    int i = 1;
    factor_list * ret = init_factor_list();
    bi_poly * u, * w, * y, * drtv, * tmp;
    drtv = formal_derivative(p);

    if (drtv->deg || (drtv->coeff[0] & 1))
    {
        u = gcd(drtv, p);
        w = quotient(p, u);
        free_poly(drtv);

        while (w->deg || (w->coeff[0] & 1) == 0)
        {
            y = (w->deg <= u->deg) ? gcd(w, u) : gcd(u, w);
            add_factor(ret, init_factor(quotient(w, y), i++));
            free_poly(w);
            w = y;
            tmp = u;
            u = quotient(tmp, y);
            free_poly(tmp);
        }

        free_poly(w);

        if (u->deg || (u->coeff[0] & 1) == 0)
        {
            tmp = u;
            u = sqrt(tmp);
            free_poly(tmp);
            factor_list * part = sff(u);
            free_poly(u);
            factor * itr;
            for (itr = part->head; itr != NULL; itr = itr->nxt)
                { itr->idx *= 2; }
            free(itr);
            append(ret, part);
        }
    }
    else
    {
        tmp = sqrt(p);
        factor_list * part = sff(tmp);
        free_poly(tmp);
        factor * itr;
        for (itr = part->head; itr != NULL; itr = itr->nxt)
            { itr->idx *= 2; }
        free(itr);
        append(ret, part);
    }

    return ret;
}

factor_list * ddf(ffa * gf, bi_poly * p)
{
    int i = 1;
    factor_list * ret = init_factor_list();
    bi_poly * fstr = copy_poly(p), * tmp;
    bi_poly * g;

    while (fstr->deg >= 2 * i)
    {
        if (fstr->deg <= gf->factors[i]->deg)
            g = gcd(fstr, gf->factors[i]);
        else
            g = gcd(gf->factors[i], fstr);

        if (g->deg)
        {
            add_factor(ret, init_factor(g, 1));
            tmp = fstr;
            fstr = quotient(tmp, g);
            free_poly(tmp);
        }

        i++;
    }

    if (fstr->deg)
        add_factor(ret, init_factor(fstr, fstr->deg));

    if (!ret->size)
    {
        add_factor(ret, init_factor(p, 1));
        return ret;
    }
    else return ret;
}

void show_matrix(int ** a, int size)
{
    int i, j;

    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            printf("%d ", a[i][j]);
        }

        printf("\n");
    }
}

int ** null_space(int ** A, int size)
{
    int n = size, r = 0, k, j, p, q;
    int * c = (int *) malloc(n * sizeof(int));
    int ** ret = (int **) malloc(n * sizeof(int *));
    bool exist;
    memset(c, -1, n * sizeof(int));

    for (k = 0; k < n; k++)
        ret[k] = (int *) calloc(n, sizeof(int));

    for (k = 0; k < n; k++)
    {
        exist = false;

        for (j = 0; j < n; j++)
        {
            if (A[k][j] && c[j] < 0)
            {
                exist = true;

                for (p = k + 1; p < n; p++)
                {
                    for (q = 0; q < n; q++)
                    {
                        if (q == j || A[k][q] == 0)
                           continue;

                        A[p][q] ^= A[k][q] * A[p][j];
                    }
                }

                for (p = 0; p < n; p++)
                    if (p != j)
                       A[k][p] = 0;
                c[j] = k;
                break;
            }
        }

        if (! exist )
        {
            ret[r][k] = 1;

            for (j = 0; j < n; j++)
                if (c[j])
                   ret[r][c[j]] = A[k][j];
            r++;
        }
    }

    return ret;
}

bi_poly ** berlekamp(bi_poly * p)
{
    int n = p->deg, t, ctr, i, j;
    int ** Q = (int **) malloc(n * sizeof(int *));
    int * a = (int *) calloc(n, sizeof(int));
    a[0] = 1;
    ctr = 0;

    for (i = 0; i < n; i++)
        Q[i] = (int *) malloc(n * sizeof(int));

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            Q[i][j] = 0;

    for (i = 1; i <= 2 * n; i++)
    {
        t = a[n - 1];

        if ((ctr & 1) == 0)
        {
            for (j = 0; j < n; j++)
            {
                Q[ctr >> 1][j] = a[j];
            }
        }

        for (j = n - 1; j > 0; j--)
        {
            a[j] = a[j - 1];

            if (t == 1 && (p->coeff[j >> 5] & bits[j & 31]))
               a[j] ^= 1;
        }

        a[0] = t * (int) (p->coeff[0] & 1);

        ctr++;
    }

    for (i = 0; i < n; i++)
        Q[i][i] ^= 1;

    Q = null_space(Q, n);

    bi_poly ** ret = (bi_poly **) malloc(n * sizeof(bi_poly *));

    for (i = 0; i < n; i++)
    {
        ret[i] = init_poly(get_size(n));

        for (j = 0; j < n; j++)
        {
            if (Q[i][j])
            {
                flip_coeff(ret[i], j);
            }
        }

        update_degree(ret[i]);
    }

    return ret;
}

#include <math.h>
#include <stdlib.h>
#include "common.h"



extern const ui bits[] =

            { 1u << 0, 1u << 1, 1u << 2, 1u << 3, 1u << 4,
              1u << 5, 1u << 6, 1u << 7, 1u << 8,
              1u << 9, 1u << 10, 1u << 11, 1u << 12, 1u << 13, 1u << 14, 1u << 15,
              1u << 16, 1u << 17, 1u << 18, 1u << 19, 1u << 20, 1u << 21, 1u << 22, 1u << 23,
              1u << 24, 1u << 25, 1u << 26, 1u << 27, 1u << 28, 1u << 29, 1u << 30, 1u << 31
            };



ffa * init_ffa(int * coeff, int len)
{
    int i, j;
    ffa * gf = (ffa *) malloc(sizeof(ffa));
    gf->WORD_SIZE = 32;
    gf->m = coeff[len - 1];
    i = gf-> m + 1;
    gf->t = (i & 31) ? (i >> 5) + 1 : (i >> 5);
    gf->fx = (bi_poly *) create_poly(gf->t, coeff, len);
    i = coeff[len - 2] + 1;
    i = (i & 31) ? (i >> 5) + 1 : (i >> 5);
    gf->rx = (bi_poly *) create_poly(i, coeff, len - 1);
    gf->td = (int *) malloc((1 << 8) * sizeof(int));
    gf->ux = (bi_poly **) malloc(32 * sizeof(bi_poly *));
    gf->factors = (bi_poly **) malloc(13 * sizeof(bi_poly *));

    for (i = 0; i < 32; i++)
        gf->ux[i] = shift_left(gf->rx, i);

    for (i = 0; i < (1 << 8); i++)
    {
        gf->td[i] = 0;

        for (j = 0; j < 8; j++)
        {
            if ((i & (1 << j)) > 0)
                gf->td[i] ^= (1 << (2 * j));
        }
    }

    for (i = 1; i < 13; i++)
    {
        j = (1 << i) + 1;
        j = (j & 31) ? (j >> 5) + 1 : (j >> 5);
        gf->factors[i] = init_poly(j);
        flip_coeff(gf->factors[i], 1);
        flip_coeff(gf->factors[i], 1 << i);
        gf->factors[i]->deg = (1 << i);
    }

    return gf;
}

void free_ffa(ffa * gf)
{
    int i;

    free(gf->td);
    free_poly(gf->fx);
    free_poly(gf->rx);

    for (i =  0; i < 32; i++)
        free_poly(gf->ux[i]);

    free(gf->ux);

    for (i = 1; i < 13; i++)
        free_poly(gf->factors[i]);

    free(gf->factors);

    free(gf);
}



bi_poly * shift_left(bi_poly * bp, int shift)
{
    int a = bp->sz, b, i, j;
    bi_poly * ret;
    b = (bp->deg + shift + 1);
    b = (b & 31) ? (b >> 5) + 1 : (b >> 5);

    if (b > a) { ret = init_poly(b); }
    else { ret = init_poly(a); }

    for (i = a - 1; i >= 0; i--)
    {
        for (j = 31; j >= 0; j--)
        {
            if ((bp->coeff[i] & bits[j]))
            {
                ret->coeff[(shift + j + (i << 5)) >> 5] ^= bits[(shift + j) & 31];
            }
        }
    }

    ret->deg = bp->deg + shift;

    return ret;
}

bi_poly * add(bi_poly * p, bi_poly * q)
{
    int i;
    bi_poly * ret;

    if (p->deg <= q->deg)
    {
        ret = copy_poly(q);

        for (i = 0; i < p->sz; i++)
            ret->coeff[i] ^= p->coeff[i];
    }
    else
    {
        ret = copy_poly(p);

        for (i = 0; i < q->sz; i++)
            ret->coeff[i] ^= q->coeff[i];
    }

    update_degree(ret);

    return ret;
}


bi_poly * multiply(bi_poly * p, bi_poly * q)
{
    int i, j, k;
    i = (p->deg + q->deg + 1);
    i = (i & 31) ? (i >> 5) + 1 : (i >> 5);
    bi_poly * ret = init_poly(i);
    bi_poly * b = copy_poly(q), * tmp;

    for (k = 0; k < 32; k++)
    {
        for (j = 0; j < p->sz; j++)
        {
            if ((p->coeff[j] & bits[k]))
            {
                for (i = 0; i < b->sz; i++)
                    ret->coeff[j + i] ^= b->coeff[i];
            }
        }

        if (k != 31)
        {
            tmp = b;
            b = shift_left(tmp, 1);
            free_poly(tmp);
        }

    }

    ret->deg = p->deg + q->deg;
    free_poly(b);

    return ret;
}

bi_poly * sqr(ffa * gf, bi_poly * p)
{
    int i;
    ui num;
    i = (p->deg << 1) + 1;
    i = (i & 31) ? (i >> 5) + 1 : (i >> 5);
    bi_poly * ret = init_poly(i);

    for (i = 0; i < p->sz; i++)
    {
        num = p->coeff[i];
        if (num == 0) continue;
        ret->coeff[2 * i] |= gf->td[(int)(num & 255)];
        num >>= 8;
        if (num == 0) continue;
        ret->coeff[2 * i] |= (gf->td[(int)(num & 255)] << 16);
        num >>= 8;
        if (num == 0) continue;
        ret->coeff[2 * i + 1] |= gf->td[(int)(num & 255)];
        num >>= 8;
        if (num == 0) continue;
        ret->coeff[2 * i + 1] |= (gf->td[(int)(num & 255)] << 16);
    }

    ret->deg = p->deg << 1;

    return ret;
}

bi_poly * gcd(bi_poly * a, bi_poly * b)
{
    int j;
    bi_poly * u = copy_poly(a);
    bi_poly * v = copy_poly(b);
    bi_poly * g1 = init_poly(1);
    bi_poly * g2 = init_poly(1);
    bi_poly * h1 = init_poly(1);
    bi_poly * h2 = init_poly(1);
    bi_poly * tmp;
    bi_poly *d, *g, *h;

    flip_coeff(g1, 0);
    flip_coeff(h2, 0);

    while (u->deg || (u->coeff[0] & 1))
    {
        j = u->deg - v->deg;

        if (j < 0)
        {
            tmp = u; u = v; v = tmp;
            tmp = g1; g1 = g2; g2 = tmp;
            tmp = h1; h1 = h2; h2 = tmp;
            j = -j;
        }

        d = shift_left(v, j);
        g = shift_left(g2, j);
        h = shift_left(h2, j);

        tmp = u;
        u = add(tmp, d);
        free_poly(tmp);
        tmp = g1;
        g1 = add(tmp, g);
        free_poly(tmp);
        tmp = h1;
        h1 = add(tmp, h);
        free_poly(tmp);

        free_poly(d);
        free_poly(g);
        free_poly(h);
    }

    free_poly(u);
    free_poly(g1);
    free_poly(g2);
    free_poly(h1);
    free_poly(h2);

    return v;
}

void reduce(ffa * gf, bi_poly * p)
{
    int i, j, k, l;

    for (i = p->deg; i >= gf->m; i--)
    {
        if ((p->coeff[i >> 5] & bits[i & 31]))
        {
            j = (i - gf->m) >> 5;
            k = (i - gf->m) - (j << 5);

            for (l = 0; l < gf->ux[k]->sz; l++)
                p->coeff[l + j] ^= gf->ux[k]->coeff[l];

            p->coeff[i >> 5] ^= bits[i & 31];
        }
    }

    update_degree(p);
}

bi_poly * reduce2(bi_poly * p, bi_poly * q)
{
    int i, j, k;
    bi_poly * r = copy_poly(p);

    if (p->deg < q->deg)
        return r;

    for (k = p->deg - q->deg; k >= 0; k--)
    {
        if ((r->coeff[(q->deg + k) >> 5] & bits[(q->deg + k) & 31]))
        {
            for (j = q->deg + k - 1; j >= k; j--)
            {
                if ((q->coeff[(j - k) >> 5] & bits[(j - k) & 31]))
                {
                    r->coeff[j >> 5] ^= bits[j & 31];
                }
            }
        }
    }

    r->deg = 0;

    for (i = q->deg - 1; i >= 0; i--)
    {
        if ((r->coeff[i >> 5] & bits[i & 31]))
        {
            r->deg = i;
            break;
        }
    }

    i = r->deg + 1;
    r->sz = (i & 31) ? (i >> 5) + 1 : (i >> 5);
    r->coeff = (ui *) realloc(r->coeff, r->sz * sizeof(ui));

    if ((r->deg & 31) != 31)
        r->coeff[r->sz - 1] &= (bits[(r->deg & 31) + 1] - 1u);

    return r;
}

bi_poly * formal_derivative(bi_poly * p)
{
    int i, j;
    bi_poly * r = copy_poly(p);

    for (i = 0; i < r->sz; i++)
    {
        for (j = 0; j < 32; j++)
        {
            if ((r->coeff[i] & bits[j]))
            {
                if ((j & 1))
                {
                    r->coeff[((i << 5) + j - 1) >> 5] ^= bits[(j - 1) & 31];
                }

                r->coeff[i] ^= bits[j];
            }
        }
    }

    update_degree(r);

    return r;
}

bool smooth(ffa * gf, bi_poly * p)
{
    int i;
    bi_poly * q = formal_derivative(p), * tmp;

    for (i = 6; i <= 12; i++)
    {
        tmp = q;
        q = multiply(tmp, gf->factors[i]);
        //show_poly(tmp);
        //show_poly(gf->factors[i]);
        free_poly(tmp);
        //show_poly(q);
        tmp = q;
        //show_poly(tmp);
        q = reduce2(tmp, p);
        free_poly(tmp);

        if (q->deg == 0 && (q->coeff[0] & 1) == 0)
        {
            free_poly(q);
            return true;
        }
    }

    free_poly(q);

    return false;
}

bool same(bi_poly * p, bi_poly * q)
{
    int i;

    if (p->deg != q->deg || p->sz != q->sz)
        return false;

    for (i = 0; i < p->sz; i++)
    {
        if ((p->coeff[i] ^ q->coeff[i]) != 0) return false;
    }

    return true;
}





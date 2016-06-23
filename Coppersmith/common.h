#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED
#define MARK  printf("here\n");

typedef unsigned int ui;

extern const ui bits[32];

typedef enum {true = 1, false = 0} bool;

typedef struct
{
    int deg;
    int sz;
    ui * coeff;
} bi_poly;

bi_poly * init_poly(int size);

bi_poly * create_poly(int size, int * ary, int len);

bi_poly * copy_poly(bi_poly * bp);

void flip_coeff(bi_poly * bp, int n);

void update_degree(bi_poly * bp);

void show_poly(bi_poly * bp);


typedef struct
{
    int WORD_SIZE;
    int m, t;
    int * td;

    bi_poly * fx;
    bi_poly * rx;

    bi_poly ** ux;
    bi_poly ** factors;
} ffa;

ffa * init_ffa(int * coeff, int len);

void free_ffa(ffa * gf);

bi_poly * shift_left(bi_poly * bp, int shift);

bi_poly * add(bi_poly * p, bi_poly * q);

bi_poly * multiply(bi_poly * p, bi_poly * q);

bi_poly * sqr(ffa * gf, bi_poly * p);

bi_poly * gcd(bi_poly * a, bi_poly * b);

void reduce(ffa * gf, bi_poly * p);

bi_poly * reduce2(bi_poly * p, bi_poly * q);

bi_poly * formal_derivative(bi_poly * p);

bool smooth(ffa * gf, bi_poly * p);

#endif // COMMON_H_INCLUDED

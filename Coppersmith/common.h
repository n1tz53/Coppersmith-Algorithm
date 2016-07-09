#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

typedef unsigned int u_int32;

extern const u_int32 bits[32];

extern const u_int32 b[5];

extern const u_int32 s[5];

typedef enum {true = 1, false = 0} bool;

/* stucture for binary polynomial */

typedef struct
{
    int deg;
    int sz;
    u_int32 * coeff;

} bi_poly;

bi_poly * init_poly(int size);

bi_poly * create_poly(int size, int * ary, int len);

bi_poly * copy_poly(bi_poly * bp);

inline void flip_coeff(bi_poly * bp, int n);

inline int get_size(int n);

inline void free_poly(bi_poly * bp);

void update_degree(bi_poly * bp);

void show_poly(bi_poly * bp);

/* stucture for parameter of field $GF(2^n)$ */

typedef struct
{
    int WORD_SIZE;
    int m, t;
    int * td;
    int bound;             /* smoothness bound */

    bi_poly * fx;
    bi_poly * rx;

    bi_poly ** ux;
    bi_poly ** factors;

} bfa;

bfa * init_bfa(int * coeff, int len);

void smooth_parameter(bfa * gf, int bd);

void free_bfa(bfa * gf);

bi_poly * shift_left(bi_poly * bp, int shift);

bi_poly * add(bi_poly * p, bi_poly * q);

bi_poly * multiply(bi_poly * p, bi_poly * q);

bi_poly * sqr(bfa * gf, bi_poly * p);

bi_poly * gcd(bi_poly * a, bi_poly * b);

void reduce(bfa * gf, bi_poly * p);

void reduce2(bi_poly * p, bi_poly * q);

bi_poly * formal_derivative(bi_poly * p);

bi_poly * quotient(bi_poly * p, bi_poly * q);

bi_poly * raise(bfa * gf, bi_poly * p, int pow);

bool smooth(bfa * gf, bi_poly * p);

/* stucture for factor of binary polynomial */

struct factor
{
    bi_poly * divisor;
    int idx;
    struct factor * nxt;

};

typedef struct factor factor;

/* list of factors */

typedef struct
{
    factor * head;
    factor * tail;
    int size;

} factor_list;

factor * init_factor(bi_poly * bp, int n);

void free_factor(factor * fact);

factor_list * init_factor_list();

void add_factor(factor_list * fact_list, factor * fact);

void append(factor_list * fact_list, factor_list * toappend);

void free_factor_list(factor_list * fact_list);

factor_list * sff(bi_poly * p);

factor_list * ddf(bfa * gf, bi_poly * p);

bi_poly ** berlekamp(bi_poly * p);

#endif // COMMON_H_INCLUDED

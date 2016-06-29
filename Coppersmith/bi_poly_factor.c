#include <stdlib.h>
#include "common.h"

factor * init_factor(bi_poly * bp, int n)
{
    factor * fact = (factor *) malloc(sizeof(factor));
    fact->idx = n;
    fact->divisor = bp;
    fact->nxt = NULL;

    return fact;
}

void free_factor(factor * fact)
{
    free_poly(fact->divisor);
    free(fact);
}

factor_list * init_factor_list()
{
    factor_list * fact_list = (factor_list *) malloc(sizeof(factor_list));
    fact_list->head = NULL;
    fact_list->tail = NULL;
    fact_list->size = 0;

    return fact_list;
}

void add_factor(factor_list * fact_list, factor * fact)
{
    if (fact_list->size)
    {
        fact_list->tail->nxt = fact;
        fact_list->tail = fact;
        fact_list->size += 1;
    }
    else { fact_list->head = fact_list->tail = fact; fact_list->size = 1; }
}

void append(factor_list * fact_list, factor_list * toappend)
{
    if (toappend->size)
    {
        if (fact_list->size)
        {
            fact_list->tail->nxt = toappend->head;
            fact_list->tail = toappend->tail;
            fact_list->size += toappend->size;
        }
        else
        {
            fact_list->head = toappend->head;
            fact_list->tail = toappend->tail;
            fact_list->size = toappend->size;
        }
    }

    toappend->head = toappend->tail = NULL;
    free(toappend);
}

void free_factor_list(factor_list * fact_list)
{
    if (fact_list->size)
    {
        factor * itr;

        for (itr = fact_list->head; itr != NULL;)
        {
            fact_list->head = itr->nxt;
            free_factor(itr);
            itr = fact_list->head;
        }

        free(itr);
    }

    free(fact_list);
}

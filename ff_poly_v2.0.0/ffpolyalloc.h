#ifndef _FFPOLYALLOC_INCLUDE
#define _FFPOLYALLOC_INCLUDE

/*
    Copyright 2017 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <assert.h>
#include "ff.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FF_POLY_DEFAULT_STACK_SIZE      (1<<20) // number of ff_t's, bytes will be 8x this
#define FF_POLY_DEFAULT_STACK_ENTRIES   8192
#define FF_POLY_STACK_ALLOC_LIMIT       1024    // number of ff_t's

struct ff_poly_stack_context {
    ff_t *next, *base, *end;
    struct ff_poly_stack_entry {
        ff_t *ptr;
        int size;
        int onstack;
    } *tab;
    int entries;
    int next_entry;
    int offstack_entries;
};

void ff_poly_stack_init (int size, int entries);
void ff_poly_stack_clear (void);

static inline ff_t *ff_poly_stack_alloc (int d)
{
    extern struct ff_poly_stack_context *_ff_sctx;
    struct ff_poly_stack_entry *e;
    
    d++;
    if ( _ff_sctx->next_entry == _ff_sctx->entries ) { fprintf (stderr, "ff_poly_stack realloc\n"); _ff_sctx->entries *= 2; _ff_sctx->tab = realloc (_ff_sctx->tab, _ff_sctx->entries); }
    e = _ff_sctx->tab + _ff_sctx->next_entry++;
    e->size = d;
    if ( d > FF_POLY_STACK_ALLOC_LIMIT || _ff_sctx->next+d > _ff_sctx->end ) {
        // fprintf (stderr, "allocating offstack entry of size %d\n", d);
        e->onstack = 0; e->ptr = malloc(d*sizeof(ff_t));
        _ff_sctx->offstack_entries++;
    } else {
        e->onstack = 1;
        e->ptr = _ff_sctx->next;
        _ff_sctx->next += d;
    }
    return e->ptr;
}

static inline void ff_poly_stack_pop (ff_t *ptr)
{
    extern struct ff_poly_stack_context *_ff_sctx;
    int i,j;
    
    for ( i = _ff_sctx->next_entry-1 ; i >= 0 ; i-- ) if ( _ff_sctx->tab[i].ptr == ptr ) break;
    assert ( i >= 0 );
    if ( ! _ff_sctx->offstack_entries ) { _ff_sctx->next_entry = i; _ff_sctx->next = ptr; return; }
    for ( j = _ff_sctx->next_entry-1 ; j >= i ; j-- ) {
        if ( _ff_sctx->tab[j].onstack ) {
            assert (_ff_sctx->next == _ff_sctx->tab[j].ptr + _ff_sctx->tab[j].size );
            _ff_sctx->next = _ff_sctx->tab[j].ptr;
        } else {
            free (_ff_sctx->tab[j].ptr);
            _ff_sctx->offstack_entries--;
        }
    }
    _ff_sctx->next_entry = i;
}


#ifdef __cplusplus
}
#endif

#endif

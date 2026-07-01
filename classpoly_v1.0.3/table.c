#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include "cstd.h"
#include "table.h"

htab_t table_common_T;                          // table_* functions share this, not thread safe

/*
    Copyright 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

/*
    We allocate three sections of memory:
    
    (1) the hash table - 64 bits per entry
    (2) the lookaside table - 128 bits per entry
    (3) overflow lists - 128 bits per entry
    
    In most cases, only (1) gets used, and this is what we want to keep in cache.   We don't really
    care about the size of the rest.
    
    Currently the number of overflow entries is equal to twice the size of the hash table.  If we exceed this we
    are getting a whole lot of collisions and should increase the table size anway.
*/

void htab_alloc (htab_t *T, int bits)
{
    if ( bits > 30 ) { err_printf ("table_alloc: bits=%d too large!\n", bits);  exit (0); }
    T->htab = malloc(sizeof(*T->htab)*(1UL<<bits));               // use malloc rather than mem_alloc - we don't need the memory initialized
    if ( ! T->htab ) { printf ("Unable to allocate %ld bytes of memory.\n", sizeof(*T->htab)*(1<<bits));  abort(); }
    T->htab_list = malloc(sizeof(*T->htab_list)*(1UL<<(bits-1))); // ditto
    if ( ! T->htab_list ) { printf ("Unable to allocate %ld bytes of memory.\n", sizeof(*T->htab_list)*(1UL<<(bits-1)));  abort(); }
    T->htab_next = 1;
    T->htab_end = 1UL<<(bits-1);
    T->htab_bits = bits;
//printf("table alloced %ld bytes\n", (1<<bits)*(sizeof(*htab)+sizeof(*htab_list)/2));
}

void htab_list_extend (htab_t *T)
{
    T->htab_end = (3*T->htab_end)/2;
    T->htab_list = realloc (T->htab_list, T->htab_end*sizeof(*T->htab_list));
    if ( ! T->htab_list ) { printf ("Unable to allocate %ld bytes of memory.\n", T->htab_end*sizeof(*T->htab_list));  abort(); }
}

void htab_free (htab_t *T)
{
    free (T->htab); T->htab = 0;
    free (T->htab_list); T->htab_list = 0;
    T->htab_bits = 0;
}


void htab_init (htab_t *T, int bits)
{
    if ( bits > T->htab_bits ) {
        if ( T->htab_bits ) htab_free(T);
        htab_alloc (T,bits);
    }
    T->htab_mask = (1UL<<bits)-1UL;
    memset (T->htab, 0, (1UL<<bits)*sizeof(*T->htab));
    T->htab_next = 1;
}

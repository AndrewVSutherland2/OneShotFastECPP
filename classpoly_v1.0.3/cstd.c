int dbg_level;
unsigned long _cstd_seed;

int ui_qsort_cmp (const void *a, const void *b)
{
    if ( *((unsigned long *)a) < *((unsigned long*)b) ) return -1;
    if ( *((unsigned long *)a) > *((unsigned long*)b) ) return 1;
    return 0;
}

int dbl_qsort_cmp (const void *a, const void *b)
{
    if ( *((double *)a) < *((double*)b) ) return -1;
    if ( *((double *)a) > *((double*)b) ) return 1;
    return 0;
}

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <ctype.h>
#include "polyparse.h"
#include "cstd.h"

/*
    Copyright 2011-2017 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#define UI_POLY_BUFSIZE             65536       // maximum length of a polynomial expression with coefficients that fit in an unsigned long (this needs to fit on the stack!)

static inline int isqdigit(char c) { return (isdigit(c) || c == '/' ); }
static inline int issign(char c) { return c == '+' || c == '-'; }
static inline int isop(char c) { return c == '+' || c == '-' || c == '*' || c == '/' || c == '^' || c =='='; }
static inline int isvecchar(char c) { return isdigit(c) || issign(c) || isspace(c) || c == ','; }
static inline int ispolychar(char c) { return isdigit(c) || isalpha(c) || isop(c) || isspace(c); }
static inline int isparen(char c) { return c == '(' || c == ')'; }
static inline int isoparen(char c) { return c == '(' || c == '['; }
static inline int iscparen(char c) { return c == ')' || c == ']'; }
static inline char closeparen(char c) { return c=='(' ? ')' : (c == '[' ? ']' : 0); }

// given a ptr s to a monomial in variables v[0],...,v[n-1], computes degreees d[0],...,d[n-1] wrt these variables and returns total degree
// *s must be equal to one of v[0],...,v[n-1] and v[n] must be null.
// returns ptr to next non-whitespace char after monomial, or null if an error occurs
static inline char *monomial_degrees (int d[], char v[], char *s)
{
    char *t;
    int e;

    for ( t = v ; *t ; t++ ) d[t-v] = 0;
    while ( isspace(*s) ) s++;
    t = strchr(v,*s);
    if ( ! t || ! *t ) return 0;
    do {
        for ( s++ ; isspace(*s) ; s++ );
        if ( *s == '^' ) { 
            for ( s++ ; isspace(*s) ; s++ );
            if ( ! isdigit(*s) ) return 0;
            e = atoi(s);
            if ( e < 0 ) return 0;
            while ( isdigit(*s) ) s++;
        } else {
            e = 1;
        }
        d[t-v] += e;
        for ( ; isspace(*s) ; s++ );
        if ( *s == '*' ) for ( s++ ; isspace(*s) ; s++ );
    } while ( (t = strchr(v,*s)) && *t );
    if ( isalpha(*s) ) return 0;   // monomial contains unexpected variables
    return s;
}

void _nop_setzero (void *f, int i) { return; }
int _nop_addto (void *f, int i, mpq_t c, void *arg) { return 1; }
int _nop_iszero (void *f, int i) { return 0; }

// returns an upper bound on the degree of a polynomial expression (must be univariate, bivariate, or homogenous in 3 variables)
// a negative value indicates an error of some sort
int poly_parse_max_degree (char *expr)
{
    char c,v[4],*s;
    int d[3], deg, n;

    n = poly_parse_variables (v, expr);
    if ( n < 0 || n > 3 ) return -1;
    for ( s = expr ; isspace(*s) ; s++ );
    c = closeparen(*s);
    if ( c ) s++;
    deg = 0;
    if ( !n ) {
        while ( *s != c && isvecchar(*s) ) if ( *s++ == ',' ) deg++;
        return deg;
    }
    while ( *s != c && (ispolychar(*s) || isparen(*s)) ) {
        if ( strchr(v,*s) ) {
            s = monomial_degrees(d, v, s);
            switch (n) {
            case 0: case 1: if ( d[0] > deg ) deg = d[0]; break;
            case 2: if ( d[0]+d[1] > deg ) deg = d[0]+d[1];
            case 3: if ( n == 3 ) return d[0]+d[1]+d[2];
            }            
        } else {
            s++;
        }
    }
    return deg;
}

// return unique alpha character in string starting at b and ending at e-1 if there is one
// returns 0 if string contains no alpha characters and -1 if it contains more than one
static inline char poly_variable (char *b, char *e)
{
    register char *t;
    register char x;
    
    for ( t = b ; t < e && *t != ',' && ! isalpha(*t) ; t++ );
    if ( t == e || *t == ',' ) return 0;
    x = *t;
    for ( t++ ; t < e && *t != ',' ; t++ ) if ( isalpha(*t) && *t != x ) return -1;
    return x;
}

// returns sorted list of alpha chars in string specified by s that are not equal to a
// if there are more than 3 such alpha chars, -1 is returned
static inline int poly_variables (char v[4], char *s, char a)
{
    register int i,n;
    
    v[n=0] = '\0'; 
    for ( ; *s && *s != ',' ; s++ ) {
        if ( isalpha(*s) && *s != a ) {
            if ( ! strchr(v,*s) ) {
                if ( n == 3 ) return -1;
                v[n++] = *s; v[n] = '\0';
            }
        }
    }
    if ( ! n ) return n;
    // bubble sort variables
    for ( a = v[1] ; a ; a = 0 ) for ( i = 0 ; i < n-1 ; i++ ) if ( v[i] > v[i+1] ) { a = v[i]; v[i] = v[i+1]; v[i+1] = a; }
    return n;
}

// poly_parse_minpoly returns a pointer to the opening paren of the minpoly in [...]/(...) or null if poly-expr does not specify a number field
// note that if expr is invalid all bets are off, results are only guaranteed for valid expressions
char *poly_parse_minpoly (char *expr)
{
    register char *s;
    
    for ( s = expr ; isspace(*s) ; s++ );
    if ( *s++ != '[' ) return 0;
    for ( ; ispolychar(*s) || *s == '(' || *s == ')' || *s == ',' ; s++ );
    if ( *s != ']' ) return 0;
    for ( s++ ; isspace(*s) ; s++ );
    if ( *s != '/' ) return 0;
    for ( s++ ; isspace(*s) ; s++ );
    if ( isparen(*s) && strchr(s,closeparen(*s)) ) return s;
    return 0;
}

char poly_parse_minpoly_variable (char *expr)
{
    register char *s = poly_parse_minpoly (expr);
    if ( ! s ) return 0;
    while ( *s && ! isalpha(*s) ) s++;
    return *s;
}

// sets vars to null-terminated list of variable names (single characters) appearing in polynomial expression, sorted alphabetically.
// returns number of variables, or -1 if expr is invalid (or has more than 3 variables)
// if vars is null, just returns number of variables
int poly_parse_variables (char vars[3], char *expr)
{
    register char *s,*t;
    char a, v[4];
    int n;
    
    s = poly_parse_minpoly(expr);
    if ( s ) {
        t = strchr(s,closeparen(*s));
        a = poly_variable(s+1,t);
    } else {
        a = 0;
    }
    n = poly_variables (v, expr, a);
    if ( vars && n > 0 ) strcpy(vars,v);
    return n;
}

// parses sign, assumes whitespace has already been stripped, handles multiple signs, sets sign=0 if no sign present
static inline char *eat_sign (int *sign, char *s)
{
    if ( ! issign(*s) ) { *sign = 0; return s; }
    for ( *sign = 1 ; issign(*s) ; s++ ) if ( *s == '-' ) *sign = -*sign;
    return s;
}

static inline char *eat_equal (int *equal, char *s)
{
    if ( *s == '=' ) { *equal = 1; return s+1; } else { *equal = 0; return s; }
}

// parses a rational coefficient, assumes whitespace has already been stripped, applies sign and eats trailing * if present, returns null on error
static inline char *eat_coeff (mpq_t c, char *s, int sign)
{
    char *t, x;
    
    if ( isdigit(*s) ) {
        for ( t = s+1 ; isqdigit(*t) ; t++ );
        x = *t; *t = '\0';
        if ( mpq_set_str (c, s, 0) < 0 ) s = 0;
        *t = x;
        if ( ! s ) return 0;
        mpq_canonicalize (c);
        s = t;
        if ( *s == '*' ) s++;
    } else if ( isalpha(*s) ) {
        mpq_set_ui (c, 1, 1);
    } else {
        s = 0;
    }
    if ( sign < 0 ) mpq_neg(c,c);
    return s;
}

// parses univariate monomial, setting e to degree, returns null on error, assumes whitespace has been stripped
static inline char *eat_monomial (int *e, char *s, char x)
{
    if ( ! x || *s != x ) { *e = 0; return s; }
    do {
        s++;
        if ( *s != '^' ) { *e = 1; return s; }
        s++;
        if ( ! isdigit(*s) ) return 0;
        *e = atoi (s);
        if ( *e < 0 ) return 0;
        while ( isdigit(*s) ) s++;
        if ( *s == '*' ) s++;
    } while ( *s == x );
    return s;
}

static inline char *eat_comma (char *s)
{
    return ( *s == ',' ? s+1 : 0 );
}

// parses homogeneous monomial in v[0],v[1],v[2] of degree d (v[4] must be null), setting e[i] to degree of v[i]
// if v[2] is null or monomial does not contain v[2] and has total degree less than d, v[2] will be set to homogenize it to degree d (but otherwise d is ignored)
// assumes whitespace has been stripped, returns null on error
static inline char *eat_homogeneous_monomial (int e[3], char *s, char v[4], int d)
{
    char *t;
    int n,v2;

    v2 = 0;
    e[0]=e[1]=e[2]=0;
    t = strchr(v,*s);
    if ( ! t || ! *t ) { e[2] = d; return s; }
    do {
        if ( t-v == 2 ) v2 = 1;
        s++;
        if ( *s == '^' ) { 
            s++;
            if ( ! isdigit(*s) ) return 0;
            n = atoi(s);
            if ( n < 0 ) return 0;
            e[t-v] += n;
            while ( isdigit(*s) ) s++;
        } else {
            e[t-v]++;
        }
        if ( *s == '*' ) s++;
    } while ( (t = strchr(v,*s)) && *t );
    if ( isalpha(*s) ) return 0;   // monomial contains unexpected variables
    if ( e[0]+e[1]+e[2] < d ) { if ( !v2 ) e[2] = d-e[0]-e[1]; } // homogenize if appropriate
    return s;
}

static void _mpq_coeff_clear (void *f, int i) { mpq_set_ui (((mpq_t*)f)[i], 0, 1); }
static int _mpq_coeff_addto (void *f, int i, mpq_t c, void *unused) { mpq_add (((mpq_t*)f)[i], ((mpq_t*)f)[i], c); return 1; }
static int _mpq_coeff_iszero (void *f, int i) { return ! mpq_sgn (((mpq_t*)f)[i]); }
static int _mpq_poly_parse (mpq_t f[], int maxd, char *expr) { return poly_parse (f, maxd, expr, _mpq_coeff_clear, _mpq_coeff_addto, _mpq_coeff_iszero, 0); }

// parses a number-field coefficient, n is the degree of the field Q(a), assumes whitespace stripped, applies sign, eats trailing * if present
static inline char *eat_nf_coeff (mpq_t c[], char *s, int n, char a, int sign)
{
    register char *t, v;
    int i, e;
    
    for ( i = 0 ; i < n ; i++ ) mpq_set_ui(c[i], 0, 1);
    if ( *s == '(' ) {
        for ( t = s + 1 ; *t && *t != ')' ; t++ );
        if ( *t != ')' ) return 0;
        v = poly_variable(s,t);
        if ( v && v != a ) return 0;
        e = _mpq_poly_parse (c, n-1, s);
        if ( e < -1 ) return 0;
        s = t+1;
        if ( *s == '*' ) s++;
        if ( sign < 0 ) for ( i = 0 ; i < n ; i++ ) mpq_neg(c[i],c[i]);
        return s;
    }
    s = eat_coeff (c[0], s, sign);
    if ( ! s ) return 0;
    s = eat_monomial (&e, s, a);
    if ( ! s || e >= n ) return 0;
    if ( *s == '*' ) s++;
    if ( e > 0 ) { mpq_set (c[e],c[0]); mpq_set_ui(c[0], 0, 1); }
    mpq_canonicalize (c[e]);
    return s;
}

// copies polynomial expression into dynamically allocated space (caller must free), removing whitespace in the process
// valid strings: noparens, [noparens], (noparens), [string-containing-instances-of-(noparens)-but-no-other-parens]
char *clean_poly (char *s)
{
    register char *buf, *t, o, c;
    register int inparen;
    
    while ( isspace(*s) ) s++;
    buf = malloc(strlen(s)+1);
    if ( isoparen(*s) ) { o = *s++; c = closeparen(o); } else { o=c=0; }
    for ( t = buf, inparen=0 ; ispolychar(*s) || isparen(*s) || *s == ',' ; s++ ) {
        if ( *s == c ) break;
        if ( *s == o ) { free(buf); return 0; }
        if ( *s == '(' ) { if ( inparen ) { free(buf); return 0; } inparen = 1; }
        if ( *s == ')' ) { if ( !inparen ) { free(buf); return 0; } inparen = 0; }
        if ( ! isspace(*s) ) *t++ = *s;
    }
    if ( inparen ) { free(buf); return 0; }
    *t = '\0';
    return buf;
}

int poly_parse_array (void *f, int off, int max, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
    mpq_t c;
    char *buf;
    register char *s;
    char cp;
    int d, e, sign;

    // Do an initial validation of parentheses (if present) and strip whitespace
    buf = clean_poly(expr);
    if ( ! buf ) { return -3; }
    if ( ! *buf ) { free(buf); return -1; }
    if ( strchr(buf+1,'(') ) { free(buf); return -3; }
    
    // verify no variables are present
    if ( poly_variable (buf, buf+strlen(buf)) != 0 ) { free(buf); return -3; }
    mpq_init (c);
    for ( e = 0 ; e < max ; e++ ) (*setzero) (f, off+e);
    s = buf;
    if ( isoparen(*s) ) { cp = *s++; } else { cp = 0; }
    d = 0;
    while ( *s ) {
        s = eat_sign (&sign,s);
        s = eat_coeff (c, s, sign);
        if ( ! s ) goto done;
        if ( d >= max ) { d = -2; goto done; }
        if ( ! (*addto) (f, off+d, c, arg) ) { d = -3; goto done; }
        s = eat_comma (s);
        if ( !s || *s == cp ) goto done;
        d++;
    }
done:
    free(buf);  mpq_clear (c);
    return d;
}


int poly_parse (void *f, int maxd, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
    mpq_t c;
    char *buf;
    char x;
    register char *s;
    int d, e, md, sign;

    // Do an initial validation of parentheses and strip whitespace
    buf = clean_poly(expr);
    if ( ! buf ) { return -3; }
    if ( ! *buf ) { free(buf); return -1; }
    if ( strchr(buf+1,'(') ) { free(buf); return -3; }
    
    // get poly variable if present
    x = poly_variable (buf, buf+strlen(buf));
    if ( x < 0 ) { free(buf); return -3; }
    mpq_init (c);
    for ( e = 0 ; e <= maxd ; e++ ) (*setzero) (f, e);
    if ( !x ) {
        d = poly_parse_array (f, 0, maxd+1, expr, setzero, addto, iszero, arg);
    } else {
        s = strchr(buf,',');
        if ( s ) *s = '\0';
        s = buf;
        md = -1; d = -3;
        while ( *s ) {
            s = eat_sign (&sign,s);
            if ( ! sign && s != buf ) goto done;
            s = eat_coeff (c, s, sign);
            if ( ! s ) goto done;
            s = eat_monomial (&e, s, x);
            if ( ! s ) goto done;
            if ( e > maxd ) { d = -2; goto done; }
            if ( ! (*addto) (f, e, c, arg) ) { d = -3; goto done; }
            if ( e > md ) md = e;
        }
        d = md;
    }
done:
    while ( d >= 0 && (*iszero)(f, d) ) d--;
    free(buf);  mpq_clear (c);
    return d;
}

// parse a polynomial of degree at most maxd with coefficients in a field of degree at most maxn
// expected format is [(coeff-poly)*x^i + ...]/(nf-poly)
// sets g to the polynomial defining the number field, of degree n (sets *dg = n)
// sets f to a list of (d+1)*n coefficients (sets *df = d), first n coefficients specify constant coeff as poly in root of g (zero padded to degree n-1)
// returns d or -2 if d > maxd or n > maxn or -3 if an error
int nf_poly_parse (void *f, int *df, void *g, int *dg, int maxd, int maxn, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
    mpq_t *c;
    char *buf;
    char a, x;
    register char *s, *t;
    int i, d, e, n, md, sign;

    *df = *dg = -1;
    
    // parse minpoly defining number field (if not present assume Q)
    s = poly_parse_minpoly (expr);
    n = poly_parse (g, maxn, s ? s : "(a)", setzero, addto, iszero, arg);
    if ( n < 1 ) return -3;
    *dg = n;
    
    // if poly is over Q, just revert to poly_parse
    if ( n == 1 ) { d = poly_parse (f, maxd, expr, setzero, addto, iszero, arg); if ( d >= 0 ) *df = d; return d; }
    
    // get minpoly variable a
    assert ( s && isoparen(*s) ); t = strchr(s,closeparen(*s)); assert (t && iscparen(*t));
    a = poly_variable(s,t);
    
    // Do an initial validation of parentheses and strip whitespace
    buf = clean_poly(expr);
    if ( ! buf ) return -3;
    if ( ! *buf) { free(buf); return -1; }
    
    // get poly variable x
    for ( s = buf ; *s ; s++ ) if ( isalpha(*s) && *s != a ) break;
    x = *s;
    
    // handle constant poly
    if ( ! x ) {
        d = poly_parse (f, n-1, buf, setzero, addto, iszero, arg);
        if ( d < -1 ) { free(buf); return -3; }
        *df = ( d < 0 ? -1 : 0 );
        return *df;
    }
    
    // allocate coefficient and initialize f
    c = malloc ((n+1)*sizeof(*c));
    for ( i = 0 ; i <= n ; i++ ) mpq_init (c[i]);
    for ( i = 0 ; i < (maxd+1)*n ; i++ ) (*setzero)(f,i);

    // parse the nf-poly
    s = buf;
    md = -1; d = -3;
    while ( *s ) {
        s = eat_sign (&sign,s);
        if ( ! sign && s != buf ) goto done;
        s = eat_nf_coeff (c, s, n, a, sign);
        if ( ! s ) goto done;
        s = eat_monomial (&e, s, x);
        if ( ! s ) goto done;
        if ( e > maxd ) { d = -2; goto done; }
        for ( i = 0 ; i < n ; i++ ) if ( ! (*addto) (f, e*n+i, c[i], arg) ) { d = -3; goto done; }
        if ( e > md ) md = e;
    }
    
    // set degree -- we need to do this to handle cancellation
    for ( d = md ; d >= 0 ; d-- ) {
        for ( i = 0 ; i < n ; i++ ) if ( ! (*iszero) (f, d*n+i) ) break;
        if ( i < n ) break;
    }
    if ( d >= 0 ) *df= d;
done:
    for ( i = 0 ; i <= n ; i++ ) mpq_clear (c[i]);
    free (c);  free (buf);
    return d;
}

int poly_parse_hyperelliptic_curve (void *f, int *df, void *h, int *dh, int maxd, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
    register char *s, *t, *u, *v;
    int d;

    char buf[strlen(expr)+1];
    for ( s = buf, t = expr ; *t ; t++ ) if (!isspace(*t)) *s++ = *t;
    *s++ = '\0';
    // handles y^2=f(x) but not y^2+h(x)y=f(x)
    if ( (s=strchr(buf,'=')) ) {
        t = buf[0] == '[' ? buf+1 : buf;
        if ( isalpha(*t) && t[1] == '^' && t[2] == '2' ) {
            *df = poly_parse (f, maxd, s+1, setzero, addto, iszero, arg);
            return *df >= 3 ? (*df-1)/2 : 0;
        }
    }

    d = 0;
    *dh = -1;
    if ( ! isoparen(buf[0]) ) { *df = poly_parse (f, maxd, buf, setzero, addto, iszero, arg);  return *df >= 3 ? (*df-1)/2 : 0; }
    s = buf+1;
    if ( isoparen(*s) ) {
        t = strchr(s,closeparen(*s));
        if ( ! t++ ) return 0;
        u = t;
        if ( *u++ != ',' ) return 0;
        if ( ! isoparen(*u) ) return 0;
        v = strchr(u,closeparen(*u));
        if ( ! v++ ) return 0;
    } else {
        v = strchr(buf,closeparen(*buf));
        if ( !v ) return 0;
        for ( t = s ; t<v && ! isalpha(*t) ; t++ );
        if ( t == v ) { *df = poly_parse (f, maxd, buf, setzero, addto, iszero, arg);  return *df >= 3 ? (*df-1)/2 : 0; }
        if ( ! (t=strchr(s,',')) )  { *df = poly_parse (f, maxd, buf, setzero, addto, iszero, arg); return *df >= 3 ? (*df-1)/2 : 0; }
        u = t+1;
    }
    char fbuf[t-s+1]; memcpy(fbuf,s,t-s); fbuf[t-s] ='\0';
    char hbuf[v-u+1]; memcpy(hbuf,u,v-u); hbuf[v-u] ='\0';
    *df = poly_parse (f, maxd, fbuf, setzero, addto, iszero, arg);
    if ( *df < -1 ) return 0;
    *dh = poly_parse (h, maxd/2, hbuf, setzero, addto, iszero, arg);
    if ( *dh < -1 ) return 0;
    d = *df >= 2*(*dh) ? *df : 2*(*dh);
    if ( d < 3 ) return 0;
    return (d-1)/2;
}

int poly_parse_superelliptic_curve (int *m, void *f, int *df, int maxd, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
    register char *s, *t;

    char buf[strlen(expr)+1];
    for ( s = buf, t = expr ; *t ; t++ ) if (!isspace(*t)) *s++ = *t;
    *s++ = '\0';
    // handles y^m=f(x)
    if ( (s=strchr(buf,'=')) ) {
        t = buf[0] == '[' ? buf+1 : buf;
        if ( isalpha(*t) && t[1] == '^' ) {
            *m = atoi(t+2);
            *df = poly_parse (f, maxd, s+1, setzero, addto, iszero, arg);
            return (*m >= 2 && *df >= 3) ? ((*df-2)*(*m-1) + *m - ui_gcd(*m,*df)) / 2 : 0;
        }
    }

    s = buf;
    if ( isoparen(*s) ) {
        s++; t = strchr(s,closeparen(*s));
        if ( !t ) return 0;
        *t = 0;
    }
    if ( ! isdigit(buf[0]) ) return 0;
    *m = atoi(buf);  if ( *m < 2 ) return 0;
    for ( s = buf+1 ; isdigit(*s) ; s++ );
    if ( *s != ',' ) return 0;
    *df = poly_parse (f, maxd, buf, setzero, addto, iszero, arg);
    return (*m >= 2 && *df >= 3) ? ((*df-2)*(*m-1) + *m - ui_gcd(*m,*df)) / 2 : 0;
}

int nf_poly_parse_hyperelliptic_curve (void *f, int *df, void *h, int *dh, void *g, int *dg, int maxd, int maxn, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
    fprintf (stderr, "nf_poly_parse_hyperelliptic_curve not yet implemented!\n");
    abort();
}

/*

    For plane curves we use the lexicographic monomial ordering, e.g. for conics, x^2, x*y, x*z, y^2, y*z, z^2.
    This means exponents are in reverse lexicographic order, e.g. <2,0,0> comes first
    Trivariate polynomials f(x,y,z) must be homogenous
    Bivariate polynomials f(x,y) will automatically be homogenized to f(x,y,z) of degree d, if possible
    expr must have at least 2 variables (possibly with zero coefficients, so you can always achieve this)
*/
int poly_parse_plane_curve (void *f, int deg, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
    mpq_t c;
    char *buf;
    int off0[deg+1];
    char v[4];
    register char *s;
    int e[3];
    int i, d, sign, equal, smult;
    
    i = poly_variables (v, expr, 0);
    if ( i == 0 ) { // just a list of rational coefficients
        d = poly_parse (f, ((deg+1)*(deg+2))/2-1, expr, setzero, addto, iszero, arg);
        if ( d == ((deg+1)*(deg+2))/2-1 ) return deg;
        for ( i = 1 ; ((i+1)*(i+2))/2 < d ; i++ ) if ( ((i+1)*(i+2))/2 == d+1 ) return -2;
        return -3;
    }

    // Do an initial validation of parentheses, remove outer parens if present, strip whitespace
    buf = clean_poly(expr);
    if ( ! buf ) return -3;
    if ( ! *buf ) { free(buf); return -1; }
    if ( strchr(buf+1,'(') ) { free(buf); return -3; }
    s = strchr (buf, ',');
    if ( s ) *s = '\0';
    
    // get variables (must be 0 (coefficient list) or at least 2)
    if ( poly_variables (v,buf,0) < 2 ) { free(buf); return -3; }
    if ( ! v[2] ) { v[2] = v[0] < v[1] ? v[1]+1 : v[0]+1; v[3] = '\0'; } // add a third variable if necessary
    for ( i = 0 ; i < ((deg+1)*(deg+2))/2 ; i++ ) (*setzero)(f,i);

    // set off0[i] to the index of the first term in which v[0] has degree i
    for ( off0[deg] = 0, i = deg-1 ; i >= 0 ; i-- ) off0[i] = off0[i+1]+(deg-i);

    mpq_init(c);
    s = buf;
    d = -3;
    equal = 0;
    smult = 1;
    while ( *s ) {
        s = eat_equal (&equal,s);
        if ( equal ) { if ( smult < 0 ) goto done; else smult = -1; }
        s = eat_sign (&sign,s);
        if ( ! sign ) { if ( ! equal && s != buf ) goto done; sign = 1; }
        s = eat_coeff (c, s, sign*smult);
        if ( ! s ) goto done;
        s = eat_homogeneous_monomial (e, s, v, deg);
        if ( ! s ) goto done;
        if ( e[0]+e[1]+e[2] != deg ) { d = -2; goto done; }
        if ( ! (*addto) (f, off0[e[0]]+e[2], c, arg) ) { d = -3; goto done; }
    }
    d = deg;
done:
    free (buf); mpq_clear(c);
    return d;
}


int nf_poly_parse_plane_curve (void *f, int deg, void *g, int *dg, int maxn, char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
    mpq_t *c;
    char *buf;
    int off0[deg+1];
    char v[4];
    char a;
    int e[3];
    register char *s, *t;
    int i, d, n, sign;

    *dg = -1;
    
    // parse minpoly defining number field (if not present assume Q)
    s = poly_parse_minpoly (expr);
    n = poly_parse (g, maxn, s ? s : "(a)", setzero, addto, iszero, arg);
    if ( n < 1 ) return -3;
    *dg = n;

    // if poly is over Q, revert to poly_parse_plane_curve
    if ( n == 1 ) return poly_parse_plane_curve (f, deg, expr, setzero, addto, iszero, arg);

    // get minpoly variable a
    assert ( s && isoparen(*s) ); t = strchr(s,closeparen(*s)); assert (t && iscparen(*t));
    a = poly_variable(s,t);

    // Do an initial validation of parentheses and strip whitespace
    buf = clean_poly(expr);
    if ( ! buf ) return -3;
    if ( ! *buf) { free(buf); return -1; }

    // get variables (must be at least 2)
    if ( poly_variables (v, buf, a) < 2 ) { free(buf); return -3; }
    if ( ! v[2] ) { v[2] = v[0] < v[1] ? v[1]+1 : v[0]+1; v[3] = '\0'; } // add a third variable if necessary

    // allocate coefficient and initialize f
    c = malloc ((maxn+1)*sizeof(*c));
    for ( i = 0 ; i <= n ; i++ ) mpq_init (c[i]);
    for ( i = 0 ; i < ((deg+2)*(deg+1)/2)*n ; i++ ) (*setzero)(f,i);

    // set off0[i] to the index of the first term in which v[0] has degree i
    for ( off0[deg] = 0, i = deg-1 ; i >= 0 ; i-- ) off0[i] = off0[i+1]+deg-i;

    // parse the nf-poly
    s = buf;
    d = -3;
    while ( *s ) {
        s = eat_sign (&sign,s);
        if ( ! sign && s != buf ) goto done;
        s = eat_nf_coeff (c, s, n, a, sign);
        if ( ! s ) goto done;
        s = eat_homogeneous_monomial (e, s, v, deg);
        if ( ! s ) goto done;
        if ( e[0]+e[1]+e[2] != deg ) { d = -2; goto done; }
        for ( i = 0 ; i < n ; i++ ) if ( ! (*addto) (f, (off0[e[0]]+e[2])*n+i, c[i], arg) ) { d = -3; goto done; }
    }
    d = deg;
done
:
    for ( i = 0 ; i <= n ; i++ ) mpq_clear (c[i]);
    free (c); free (buf);
    return d;
}

/*
    Conic covers are curves in P^3 of the form g(x,y,z)=0, w^2=f(x,y,z), where g is homogeneous deg 2 (so 6 coeffs) and f is homogeneous deg d (binom(d+2,2) coeffs)
    We support 3 input formats:
        (1) [g(x,y,z),f(x,y,z)]
        (2) [[g1,...,g6],[f1,...,fn]] where n=binom(d+2,2), coefficient order is lex monomial (as with plane curves)
        (3) [[s,t],[h0,....h(d-1)],[f0,...,fd]], specifies the curve x^2-sy^2-tz^2=0, w^2=f(y,z)+xh(y,z), where deg(f)=d and deg(h)=d-1
*/
int poly_parse_conic_cover (void *g, void *f, int d ,char *expr, void (*setzero)(void *f, int i), int (*addto)(void *f, int i, mpq_t c, void *arg), int (*iszero)(void *f, int i), void *arg)
{
    char *buf;
    mpq_t c;
    char *s,*t;
    int i, r, sgn;

    r = -3;
    buf = malloc(strlen(expr)+1);
    for ( s = buf, t = expr ; *t ; t++ ) if ( ! isspace(*t) ) *s++ = *t;
    *s = '\0';
    mpq_init(c);
    if ( buf[1] != '[' ) { // case (1)
        s = strchr(buf,',');
        if ( (r=poly_parse_plane_conic(g, buf+1, setzero, addto, iszero, arg)) != 2 ) goto done;
        if ( (r=poly_parse_plane_curve(f, d, s+1, setzero, addto, iszero, arg)) != d ) goto done;
    } else {
        s = strchr(buf+2,','); if ( !s ) goto done;
        s = strchr(s+1,','); if ( !s ) goto done;
        t = strchr(buf+2,']'); if ( !t ) goto done;
        if ( *++t != ',' ) goto done;
        if ( s < t ) {  // case (2)
            if ( (r=poly_parse_plane_conic (g, buf+1, setzero, addto, iszero, arg)) != 2 ) goto done;
            if ( (r=poly_parse_plane_curve (f, d, t+1, setzero, addto, iszero, arg)) != d ) goto done;
        } else { // case (3)
            if ( s != t ) goto done;
            for ( i = 0 ; i < 6 ; i++ ) (*setzero) (g, i);
            mpq_set_ui(c,1,1);  (*addto) (g, 0, c, arg);
            // negate sgn to get x^2-s*y^2-t*z^2
            s = eat_sign (&sgn, buf+2); s = eat_coeff (c, s, -sgn); if ( ! s || *s++ != ',' ) goto done;
            (*addto) (g, 3, c, arg);
            s = eat_sign (&sgn, s); s = eat_coeff (c, s, -sgn); if ( ! s || *s++ != ']' ) goto done;
            (*addto) (g, 5, c, arg);
            if ( *s++ != ',' || *s != '[' ) goto done;
            for ( i = 0 ; i < 2*d+1 ; i++ ) (*setzero) (f, i);
            if ( (r=poly_parse_array (f, ((d+1)*(d+2))/2-(d+1)-d, d, s, setzero, addto, iszero, arg)) != d-1 ) goto done;
            s = strchr(s,']'); if ( !s++ || *s++ != ',' || *s != '[' ) goto done;
            r = poly_parse_array (f, ((d+1)*(d+2))/2-(d+1), d+1, s, setzero, addto, iszero, arg);
        }
    }
done:
    mpq_clear(c);
    free(buf);
    return r;
}


void _ui_coeff_setzero (void *f, int i) { ((unsigned long*)f)[i] = 0; }
int _ui_mod_p_coeff_addto (void *f, int i, mpq_t c, void *arg)
{
    unsigned long p, x;
    
    p = *((unsigned long *)arg);
    if (  mpz_cmp_ui (mpq_denref(c),1) == 0 ) {
        ((unsigned long *)f)[i] = mpz_fdiv_ui (mpq_numref(c), p);
    } else {
        x = ui_inverse (mpz_fdiv_ui(mpq_denref(c),p), p);
        if ( ! x ) return 0;
        // use the numerater of c as an mpz_t that we can use for workspace to avoid worrying about overflow (we are allowed to trash c)
        mpz_set_ui (mpq_numref(c), mpz_fdiv_ui (mpq_numref(c), p));
        mpz_mul_ui (mpq_numref(c), mpq_numref(c), x);
        mpz_set_ui (mpq_numref(c), mpz_fdiv_ui (mpq_numref(c),p));
        mpz_add_ui (mpq_numref(c), mpq_numref(c), ((unsigned long *)f)[i]);
        ((unsigned long *)f)[i] = mpz_fdiv_ui (mpq_numref(c), p);
    }
    return 1;
}
int _ui_coeff_iszero (void *f, int i) { return ((unsigned long *)f)[i] ? 0 : 1; }

int ui_poly_parse_mod_p (unsigned long f[], int maxd, char *expr, unsigned long p) { return poly_parse (f, maxd, expr, _ui_coeff_setzero, _ui_mod_p_coeff_addto, _ui_coeff_iszero, &p); }
int ui_poly_parse_plane_curve_mod_p (unsigned long f[], int d, char *expr, unsigned long p) { return poly_parse_plane_curve (f, d, expr, _ui_coeff_setzero, _ui_mod_p_coeff_addto, _ui_coeff_iszero, &p); }

void _i_coeff_setzero (void *f, int i) { ((long*)f)[i] = 0; }
int _i_coeff_addto (void *f, int i, mpq_t c, void *unused)
{ 
    long *pc;
    
    if ( mpz_cmp_ui (mpq_denref(c),1) != 0 ) return 0;
    pc = (long*)f + i;
    if ( *pc > 0 ) mpz_add_ui (mpq_numref(c), mpq_numref(c), *pc);
    else if ( *pc < 0 ) mpz_sub_ui  (mpq_numref(c), mpq_numref(c), -(*pc));
    if ( ! mpz_fits_slong_p (mpq_numref(c)) ) return 0;
    *pc = mpz_get_si(mpq_numref(c));
    return 1;
}
int _i_coeff_iszero (void *f, int i) { return ((long *)f)[i] ? 0 : 1; }

int i_poly_parse (long f[], int maxd, char *expr)
    { return poly_parse (f, maxd, expr, _i_coeff_setzero, _i_coeff_addto, _i_coeff_iszero, 0); }
int i_nf_poly_parse (long f[], int *df, long g[], int *dg, int maxd, int maxn, char *expr)
    { return nf_poly_parse (f, df, g, dg, maxd, maxn, expr, _i_coeff_setzero, _i_coeff_addto, _i_coeff_iszero, 0); }
int i_poly_parse_plane_curve (long f[], int d, char *expr)
    { return poly_parse_plane_curve (f, d, expr, _i_coeff_setzero, _i_coeff_addto, _i_coeff_iszero, 0); }
int i_poly_parse_hyperelliptic_curve (long f[], int *df, long h[], int *dh, int maxd, char *expr)
    { return poly_parse_hyperelliptic_curve (f, df, h, dh, maxd, expr, _i_coeff_setzero, _i_coeff_addto, _i_coeff_iszero, 0); }
int i_nf_poly_parse_hyperelliptic_curve (long f[], int *df, long h[], int *dh, long g[], int *dg, int maxd, int maxn, char *expr)
    { return nf_poly_parse_hyperelliptic_curve (f, df, h, dh, g, dg, maxd, maxn, expr, _i_coeff_setzero, _i_coeff_addto, _i_coeff_iszero, 0); }
int i_poly_parse_superelliptic_curve (int *m, long f[], int *df, int maxd, char *expr)
    { return poly_parse_superelliptic_curve (m, f, df, maxd, expr, _i_coeff_setzero, _i_coeff_addto, _i_coeff_iszero, 0); }
int i_poly_parse_conic_cover (long g[6], long f[], int d, char *expr)
    { return poly_parse_conic_cover (g, f, d, expr, _i_coeff_setzero, _i_coeff_addto, _i_coeff_iszero, 0); }

void ui_poly_print (unsigned long f[], int d_f)
{
    char buf[UI_POLY_BUFSIZE];

    ui_poly_sprint (buf, f, d_f);
    puts (buf);
}


int ui_poly_sprint (char *s, unsigned long f[], int d_f)
{
    char *t;
    int i;
    
    if ( d_f < 0 ) { strcpy (s, "[zero polynomial]");  return strlen(s); }
    t = s;
    if ( d_f >= 2 ) {
        if ( f[d_f] != 1 ) t += sprintf (t, "[%lux^%d", f[d_f], d_f); else  t += sprintf (t, "[x^%d", d_f);
    } else if ( d_f == 1 ) {
        if ( f[d_f] != 1 ) t += sprintf (t, "[%lux", f[d_f]); else  t += sprintf (t, "[x");
    } else {
        t += sprintf (t, "[%lu", f[d_f]);
    }
    for ( i = d_f-1 ; i >= 0 ; i-- ) {
        if ( f[i] ) {
            if ( i >= 2 ) {
                t += sprintf (t, " + %lux^%d", f[i], i);
            } else if ( i == 1 ) {
                t += sprintf (t, " + %lux", f[i]);
            } else {
                t += sprintf (t, " + %lu", f[i]);
            }
        }
    }
    *t++ = ']';
    *t= '\0';
    return t-s;
}

void i_poly_print (long f[], int df)
{
    char buf[UI_POLY_BUFSIZE];

    i_poly_sprint (buf, f, df);
    puts (buf);
}

int i_poly_sprint (char *s, long f[], int df)
{
    register char *t;
    register int i;
    
    while ( df >= 0 && ! f[df] ) df--;
    if ( df < 0 ) { strcpy (s, "[0]");  return strlen(s); }
    t = s;
    if ( df >= 2 ) {
        if ( f[df] == 1 ) { t += sprintf (t, "[x^%d", df); } else if ( f[df] == -1 ) { t += sprintf (t, "[-x^%d", df); } else { t += sprintf (t, "[%ld*x^%d", f[df], df); }
    } else if ( df == 1 ) {
        if ( f[df] == 1 ) { t += sprintf (t, "[x"); } else if ( f[df] == -1 ) { t += sprintf (t, "[-x"); } else { t += sprintf (t, "[%ld*x", f[df]); }
    } else {
        t += sprintf (t, "[%ld", f[df]);
    }
    for ( i = df-1 ; i >= 0 ; i-- ) {
        if ( f[i] > 1) {
            t += sprintf (t, " + %ld", f[i]); if ( i ) *t++ = '*';
        } else if ( f[i] == 1 ) {
            t += sprintf (t, " + ");
        } else if ( f[i] == -1 ) {
            t += sprintf (t, " - ");
        } else if ( f[i] < -1 ) {
            t += sprintf (t, " - %ld", -f[i]); if ( i ) *t++ = '*';
        }
        if ( f[i] ) {
            if ( i >= 2 ) {
                t += sprintf (t, "x^%d", i);
            } else if ( i == 1 ) {
                t += sprintf (t, "x");
            } else {
                if ( f[i] == 1 || f[i] == -1 ) *t++ = '1';
            }
        }
    }
    *t++ = ']'; *t = '\0';
    return t-s;
}

int i_poly_sprint_coeffs (char *buf, long f[], int df)
{
    char *s = buf; *s++ = '[';
    for ( int i = 0 ; i <= df ; i++ ) { if ( i ) { *s++ = ','; } s += sprintf(s,"%ld",f[i]); }
    *s++ = ']'; *s++ = '\0';
    return s-buf;
}

void i_poly_print_coeffs (long f[], int df)
{
    char *buf;
    int m, n;
    
    n = 3;
    for ( int i = 0 ; i <= df ; i++ ) n += 20;
    assert ( n > 0 );
    buf = malloc (n);
    m = i_poly_sprint_coeffs (buf,f,df);
    assert (m<n);
    puts (buf);
    free (buf);
}


void i_poly_print_plane_curve (long f[], int d)
{
    char buf[UI_POLY_BUFSIZE];
    i_poly_sprint_plane_curve (buf, f, d);
    puts (buf);
}

int i_poly_sprint_plane_curve (char *s, long f[], int d)
{
    register char *t;
    register int i, x, y, z;
    int first = 1;
    
    t = s;
    i = 0;
    *t++ = '[';
    for ( x = d ; x >= 0 ; x-- ) {
        for ( y = d-x ; y >= 0 ; y--, i++ ) {
            z = d-x-y;
            if ( f[i] > 1) {
                t += sprintf (t, " + %ld", f[i]); *t++ = '*';
            } else if ( f[i] == 1 ) {
                t += sprintf (t, " + ");
            } else if ( f[i] == -1 ) {
                t += sprintf (t, " - ");
            } else if ( f[i] < -1 ) {
                t += sprintf (t, " - %ld", -f[i]); *t++ = '*';
            }
            if ( f[i] ) {
                if ( first ) { if ( f[i] > 0 ) { memmove (s+1,s+4,t-s-3); t-=3; } else { s[1] = '-'; memmove(s+2,s+4,t-s-2); t-=2; } first = 0; }
                if ( x >= 2 ) {
                    t += sprintf (t, "x^%d", x);
                } else if ( x == 1 ) {
                    t += sprintf (t, "x");
                }
                if ( x && x < d ) *t++ = '*';
                if ( y >= 2 ) {
                    t += sprintf (t, "y^%d", y);
                } else if ( y == 1 ) {
                    t += sprintf (t, "y");
                }
                if ( y && x+y < d ) *t++ = '*';
                if ( z >= 2 ) {
                    t += sprintf (t, "z^%d", z);
                } else if ( z == 1 ) {
                    t += sprintf (t, "z");
                }
            }
        }
    }
    *t++ = ']'; *t = '\0';
    return t-s;
}

void i_poly_print_conic_cover (long g[6], long f[], int d)
{
    char buf[UI_POLY_BUFSIZE];
    char *s;
    
    i_poly_sprint_plane_curve (buf, g, 2);
    s = buf+strlen(buf)-1;
    *s = ',';
    i_poly_sprint_plane_curve (s+1, f, d);
    *(s+1) = ' ';
    puts (buf);
}

int i_poly_sprint_conic_cover (char *buf, long g[6], long f[], int d)
{
    char *s;

    i_poly_sprint_plane_curve (buf, g, 2);
    s = buf+strlen(buf)-1;
    *s = ',';
    i_poly_sprint_plane_curve (s+1, f, d);
    *(s+1) = ' ';
    return strlen(buf);
}

void i_hyperelliptic_curve_print (long f[], int df, long h[], int dh)
{
    char buf[UI_POLY_BUFSIZE];
    i_hyperelliptic_curve_sprint (buf, f, df, h, dh);
    puts (buf);
}

void i_hyperelliptic_curve_print_coeffs (long f[], int df, long h[], int dh)
{
    char buf[UI_POLY_BUFSIZE];
    i_hyperelliptic_curve_sprint_coeffs (buf, f, df, h, dh);
    puts (buf);
}

void i_superelliptic_curve_print (int m, long f[], int df)
{
    char buf[UI_POLY_BUFSIZE];
    i_superelliptic_curve_sprint (buf, m, f, df);
    puts (buf);
}

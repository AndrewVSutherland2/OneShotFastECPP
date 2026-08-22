/* short2.gp -- find certificates in the revised (August 2026) short ECPP format,
 * adapted from ShortPrimalityProofs/short.gp (Fable 5); the ntop-parameterized
 * entry points support mid-chain repairs.
 *
 * Revised format (n = ceil(log2 p_0), FIXED for the chain):
 *   - the smoothness bound is  B = ceil(n^2/log2 n);
 *   - the filler satisfies the radical cap  log2 rad(m) <= ceil(n/log2 n);
 *   - recursion requires  p_{i+1} > B^2; a rough cofactor below B^2 (hence 1 or
 *     prime -- it certifies itself by size) ends the chain.
 *
 * Usage:
 *   echo 'printshort2(nextprime(10^30))' | gp -q short2.gp
 *   echo 'printshort2from(P, NTOP)' | gp -q short2.gp     \\ repair: chain at modulus P
 *                                                            with top-level bit length NTOP
 */

default(parisizemax, 2^32);

SC_curves = 0;
SC_seacalls = 0;
SC_tlim = 20;                                            \\ seconds per rough-part factorization
SC_factorflags = 0;
SC_branchcurves = 64;
SC_maxcurves = 0;

/* B-smooth part s of N and the rough cofactor r = N/s. */
smoothpart(N, B) = {
  my(s = 1, r = N);
  forprime(q = 2, B, while(r % q == 0, r /= q; s *= q));
  [s, r];
};

scbound(p) = sqrtint(p) + 1 + sqrtint(4 * sqrtint(p));

scnonsquare(p) = {
  my(d);
  if(p % 4 == 3, return(p-1));
  d = 2; while(kronecker(d, p) != -1, d++);
  d;
};

/* radical cap: log2 rad(m) <= radlim (exact: a squarefree 2-power is just 2) */
sc2radok(m, radlim) = #binary(vecprod(factor(m)[,1])) <= radlim;

scpoint(E, N, o, fo) = {
  my(Q, ok);
  for(t = 1, 64,
    Q = ellmul(E, random(E), N/o);
    if(#Q == 1, next);
    ok = 1;
    for(i = 1, #fo, if(#ellmul(E, Q, o/fo[i]) == 1, ok = 0; break));
    if(ok, return(Q))
  );
  0;
};

/* Try one computed curve order.  B = smoothness bound, B2 = B^2 = recursion floor
 * (a rough cofactor below B2 is terminal), radlim = ceil(n/log2 n) radical cap. */
sctryorder2(p, B, B2, radlim, A, xden, E, N, L, rt) = {
  my(sr, s, R, dv, F, m, o, q, Q);
    sr = smoothpart(N, B); s = sr[1]; R = sr[2];
    dv = divisors(s);

    \\ (a) terminal level: o = m is B-smooth by itself, radical-capped, in the window
    for(i = 1, #dv, m = dv[i];
      if(m > L && m < factor(m)[1,1] * L && sc2radok(m, radlim),
        Q = scpoint(E, N, m, factor(m)[,1]);
        if(Q != 0, return([A, lift(Q[1]/Mod(xden, p)), m, 1]))));

    \\ (a') terminal level with a self-certifying prime: the rough part is below B^2,
    \\      hence automatically prime (no primality test), and stays inside o
    if(R > 1 && R < B2,
      for(i = 1, #dv, m = dv[i]; o = m * R;
        if(m > 1 && o > L && o < factor(m)[1,1] * L && sc2radok(m, radlim),
          Q = scpoint(E, N, o, concat(factor(m)[,1], [R]~));
          if(Q != 0, return([A, lift(Q[1]/Mod(xden, p)), o, 1])))));

    \\ (b) cheap descent: the rough part is itself prime, above the B^2 recursion floor
    if(s > L && R > B2 && R < rt && ispseudoprime(R),
      q = R;
      for(i = 1, #dv, m = dv[i]; o = m * q;
        if(m > 1 && o > L && o < factor(m)[1,1] * L && sc2radok(m, radlim),
          Q = scpoint(E, N, o, concat(factor(m)[,1], [q]~));
          if(Q != 0, return([A, lift(Q[1]/Mod(xden, p)), o, q])))));

    \\ (c) general descent: o = m*q, q a prime factor of the rough part dug out under
    \\     a time budget, m | s radical-capped
    if(SC_tlim > 0 && R > 1 && s * rt > L,
      F = iferr(alarm(SC_tlim, factorint(R, SC_factorflags)[,1]), e, 0);
      if(type(F) == "t_COL",
        for(j = 1, #F, q = F[j];
          if(q < rt && (q < B2 || ispseudoprime(q)),
            for(i = 1, #dv, m = dv[i]; o = m * q;
              if(m > 1 && o > L && o < factor(m)[1,1] * L && sc2radok(m, radlim),
                Q = scpoint(E, N, o, concat(factor(m)[,1], [q]~));
                if(Q != 0, return([A, lift(Q[1]/Mod(xden, p)), o, if(q < B2, 1, q)]))))))));
  0;
};

sclevel2(p, B, B2, radlim, {stopcurve = 0}) = {
  my(L = scbound(p), rt = sqrtint(p), d = scnonsquare(p), A, E, Et, N, lev);
  while(1,
    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    if(stopcurve && SC_curves >= stopcurve, return(0));
    A = random(p); E = ellinit([0, A, 0, 1, 0], p);
    SC_curves++;
    if(#E == 0, next);
    SC_seacalls++;
    N = ellcard(E);
    lev = sctryorder2(p, B, B2, radlim, A, 1, E, N, L, rt);
    if(lev != 0, return(lev));

    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    if(stopcurve && SC_curves >= stopcurve, return(0));
    Et = ellinit([0, (d*A)%p, 0, (d^2)%p, 0], p);
    SC_curves++;
    lev = sctryorder2(p, B, B2, radlim, A, d, Et, 2*p+2-N, L, rt);
    if(lev != 0, return(lev))
  );
};

scchain2(p, B, B2, radlim, stopcurve) = {
  if(p == 1, return([]));
  my(lev, tail, childstop);
  while(1,
    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    if(stopcurve && SC_curves >= stopcurve, return(0));
    lev = sclevel2(p, B, B2, radlim, stopcurve);
    if(lev == 0, return(0));
    childstop = if(SC_branchcurves, SC_curves + SC_branchcurves, 0);
    if(stopcurve && (!childstop || stopcurve < childstop), childstop = stopcurve);
    tail = scchain2(lev[4], B, B2, radlim, childstop);
    if(type(tail) == "t_VEC", return(concat([lev[1], lev[2], lev[3]], tail)))
  );
};

/* Chain from modulus p with top-level bit length ntop (repairs: p is a mid-chain
 * modulus of a chain whose p_0 has ntop bits; fresh chains: ntop = #binary(p)). */
shortcert2from(p, ntop) = {
  if(!ispseudoprime(p), error("short2: p is composite"));
  if(p < 5, error("short2: need p >= 5"));
  my(lg = log(ntop)/log(2), B = ceil(ntop^2/lg), B2 = B^2, radlim = ceil(ntop/lg), tail);
  SC_curves = 0; SC_seacalls = 0;
  tail = scchain2(p, B, B2, radlim, 0);
  if(type(tail) != "t_VEC", return(0));
  concat([p], tail);
};

printshort2from(p, ntop) = {
  my(c = shortcert2from(p, ntop));
  for(i = 1, #c, printf("%d%s", c[i], if(i < #c, " ", "\n")));
};

printshort2(p) = printshort2from(p, #binary(p));

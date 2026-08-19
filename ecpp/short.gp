/* short.gp -- compute a short ECPP certificate for a probable prime p >= 5.
 * Written by Fable 5.  A toy implementation, in the style of oneshot.gp; the goal is
 * clarity, not speed.
 *
 * Format: github.com/AndrewVSutherland/ShortPrimalityProofs.  The certificate is the flat
 * sequence  (p_0, A_0, x_0, o_0, A_1, x_1, o_1, ..., A_k, x_k, o_k)  with o_i = m_i*p_{i+1},
 * p_{k+1} = 1 and n = ceil(log_2 p_0) FIXED for the whole chain: on E_{A_i} : y^2 = x^3 +
 * A_i x^2 + x over F_{p_i} the point with x-coordinate x_i has order exactly o_i, where m_i
 * is n^2-smooth, p_{i+1} is a prime in (n^2, sqrt(p_i)) proven prime by the next level, and
 * L_i < o_i < r_i L_i with L_i = (p_i^{1/4}+1)^2 and r_i the least prime divisor of m_i.
 *
 * Method.  For each level we search random curves E_A/F_{p_i}, compute N = #E_A by SEA
 * (ellcard), and look for a divisor o | N inside the narrow window (L_i, r_i L_i) -- note
 * o is only about sqrt(p_i), so most of N is discarded.  A usable o must be of the shape
 * (n^2-smooth) * (one prime), so we trial-divide N up to n^2 to get its smooth part s,
 * and then look for the single large prime among the factors of the rough part N/s.  That
 * last step is the expensive one: the rough part is factored under an `alarm` time budget
 * (SC_tlim seconds) and the curve is abandoned if the budget runs out -- an early abort
 * that keeps the search on curves whose order factors easily.  A level with a fully
 * n^2-smooth o (p_{i+1} = 1) ends the chain.
 *
 * Usage:
 *     echo 'printshort(nextprime(10^30))' | gp -q short.gp
 *     echo 'SC_tlim = 60; printshort(nextprime(10^100))' | gp -q short.gp
 */

default(parisizemax, 2^30);                              \\ allow the stack to grow for large SEA computations

SC_curves = 0;                                           \\ curve orders tried (all levels)
SC_seacalls = 0;                                         \\ expensive ellcard calls (one call gives two twists)
SC_tlim = 20;                                            \\ seconds allowed per rough-part factorization
                                                         \\ (SC_tlim = 0 disables step (c): then the only
                                                         \\  large prime considered is a prime rough part)
SC_factorflags = 0;                                      \\ factorint flags; 9 is a fast partial-factor mode
SC_branchcurves = 64;                                    \\ backtrack after this many curves below a descent
                                                         \\ (0 retains the original unlimited commitment)
SC_maxcurves = 0;                                        \\ 0 = unlimited; else give up after this many curves

/* n^2-smooth part s of N and the rough cofactor r = N/s, by trial division over primes <= B */
smoothpart(N, B) = {
  my(s = 1, r = N);
  forprime(q = 2, B, while(r % q == 0, r /= q; s *= q));
  [s, r];
};

scbound(p) = sqrtint(p) + 1 + sqrtint(4 * sqrtint(p));   \\ integer form of L = (p^{1/4}+1)^2

/* A small quadratic nonresidue modulo p. */
scnonsquare(p) = {
  my(d);
  if(p % 4 == 3, return(p-1));
  d = 2; while(kronecker(d, p) != -1, d++);
  d;
};

/* A point of order exactly o on E (N = #E, o | N, fo = the primes of o), or 0 if none found. */
scpoint(E, N, o, fo) = {
  my(Q, ok);
  for(t = 1, 64,
    Q = ellmul(E, random(E), N/o);                       \\ order(Q) divides o
    if(#Q == 1, next);                                   \\ Q = O, resample
    ok = 1;
    for(i = 1, #fo, if(#ellmul(E, Q, o/fo[i]) == 1, ok = 0; break));
    if(ok, return(Q))
  );
  0;
};

/* Try one curve order that has already been computed.  A is the certificate parameter and
 * xden maps the model's x-coordinate to the certificate coordinate x/xden. */
sctryorder(p, n2, A, xden, E, N, L, rt) = {
  my(sr, s, R, dv, F, m, o, q, Q);
    sr = smoothpart(N, n2); s = sr[1]; R = sr[2];
    dv = divisors(s);

    \\ (a) terminal level: o = m is n^2-smooth all by itself, and lands in the window
    for(i = 1, #dv, m = dv[i];
      if(m > L && m < factor(m)[1,1] * L,
        Q = scpoint(E, N, m, factor(m)[,1]);
        if(Q != 0, return([A, lift(Q[1]/Mod(xden, p)), m, 1]))));

    \\ (b) cheap descent: the rough part is itself prime, so q = R needs no factoring.  Then the
    \\     discarded cofactor N/o = s/m divides s, i.e. the whole ~sqrt(p) worth of N that we
    \\     throw away must be n^2-smooth.  That forces s > L, of density rho(ln(sqrt p)/ln(n^2)):
    \\     measured on 200 random curves it is 25% at p ~ 10^10 but already 0/200 by p ~ 10^40,
    \\     so this cheap case alone does not scale -- hence (c).  Costs one primality test.
    if(s > L && R > 1 && R < rt && ispseudoprime(R),
      q = R;
      for(i = 1, #dv, m = dv[i]; o = m * q;
        if(m > 1 && o > L && o < factor(m)[1,1] * L,
          Q = scpoint(E, N, o, concat(factor(m)[,1], [q]~));
          if(Q != 0, return([A, lift(Q[1]/Mod(xden, p)), o, q])))));

    \\ (c) general descent: o = m*q with q any prime factor of the rough part, m | s, m > 1.
    \\     Here the discarded cofactor N/o is unconstrained, which is far likelier -- but q must
    \\     be dug out of R, so we factor R under a time budget and abandon slow curves.
    if(SC_tlim > 0 && R > 1 && s * rt > L,
      F = iferr(alarm(SC_tlim, factorint(R, SC_factorflags)[,1]), e, 0);
                                                         \\ flag 9 may leave a composite cofactor, filtered below
      if(type(F) == "t_COL",
        for(j = 1, #F, q = F[j];
          if(q > n2 && q < rt && ispseudoprime(q),
            for(i = 1, #dv, m = dv[i]; o = m * q;
              if(m > 1 && o > L && o < factor(m)[1,1] * L,
                Q = scpoint(E, N, o, concat(factor(m)[,1], [q]~));
                if(Q != 0, return([A, lift(Q[1]/Mod(xden, p)), o, q]))))))));
  0;
};

/* One level: search random curves over F_p for [A, x0, o, p_next].  p_next = 1 ends the chain.
 * If d is a quadratic nonresidue, y^2 = x^3+d*A*x^2+d^2*x is the quadratic twist of E_A,
 * and its model coordinate X maps to certificate coordinate x=X/d on d*y^2=x^3+A*x^2+x.
 * The two orders sum to 2p+2, so test both for the cost of one SEA computation. */
sclevel(p, n2, {stopcurve = 0}) = {
  my(L = scbound(p), rt = sqrtint(p), d = scnonsquare(p), A, E, Et, N, lev);
  while(1,
    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    if(stopcurve && SC_curves >= stopcurve, return(0));
    A = random(p); E = ellinit([0, A, 0, 1, 0], p);
    SC_curves++;
    if(#E == 0, next);                                   \\ singular (A = +-2 mod p)
    SC_seacalls++;
    N = ellcard(E);                                      \\ SEA
    lev = sctryorder(p, n2, A, 1, E, N, L, rt);
    if(lev != 0, return(lev));

    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    if(stopcurve && SC_curves >= stopcurve, return(0));
    Et = ellinit([0, (d*A)%p, 0, (d^2)%p, 0], p);
    SC_curves++;
    lev = sctryorder(p, n2, A, d, Et, 2*p+2-N, L, rt);
    if(lev != 0, return(lev))
  );
};

/* Recursively complete a chain, abandoning a child subtree at the absolute `stopcurve`.
 * A zero stop is unlimited (used only for the root), while every descent receives at most
 * the SC_branchcurves budget.  This avoids committing a worker forever to an unlucky child. */
scchain(p, n2, stopcurve) = {
  if(p == 1, return([]));
  my(lev, tail, childstop);
  while(1,
    if(SC_maxcurves && SC_curves >= SC_maxcurves, return(0));
    if(stopcurve && SC_curves >= stopcurve, return(0));
    lev = sclevel(p, n2, stopcurve);
    if(lev == 0, return(0));
    childstop = if(SC_branchcurves, SC_curves + SC_branchcurves, 0);
    if(stopcurve && (!childstop || stopcurve < childstop), childstop = stopcurve);
    tail = scchain(lev[4], n2, childstop);
    if(type(tail) == "t_VEC", return(concat([lev[1], lev[2], lev[3]], tail)))
  );
};

/* The full chain: returns the flat sequence (p, A_0, x_0, o_0, ..., A_k, x_k, o_k). */
shortcert(p) = {
  if(!ispseudoprime(p), error("short: p is composite"));
  if(p < 5, error("short: need p >= 5"));
  my(n = #binary(p), n2 = n^2, tail);
  SC_curves = 0; SC_seacalls = 0;
  tail = scchain(p, n2, 0);
  if(type(tail) != "t_VEC", return(0));
  concat([p], tail);
};

printshort(p) = {                                        \\ "p A0 x0 o0 A1 x1 o1 ..."
  my(c = shortcert(p));
  for(i = 1, #c, printf("%d%s", c[i], if(i < #c, " ", "\n")));
};

/* restricted.gp -- restricted (n^2-smooth) one-shot ECPP certificates.
 *
 * A restricted one-shot certificate is a one-shot ECPP (p, A, x0, m) whose point
 * order m is entirely n^2-smooth, n = ceil(log_2 p) -- the k = 0 case of the
 * OneShotPrimalityProofs format (no q_i list; voneshot.py verifies it with an
 * empty prime list).  Search is SEA sampling: random Montgomery A, ellcard, keep
 * a curve whose n^2-smooth part of #E exceeds L = (p^{1/4}+1)^2, then reduce a
 * point to minimal smooth order in (L, L*r).  Heuristic cost: rho(u/2)^{-1}
 * candidates with u = ln p / ln n^2, i.e. p^{1/4+o(1)} -- the square of the
 * unrestricted n^4 search's candidate count.  See reports/ecpp-varieties/.
 *
 * The search routine is oneshot.gp from
 * github.com/AndrewVSutherland/OneShotPrimalityProofs (Opus 4.8 Max), inlined
 * unchanged; restrictedcert() simply calls it with B = n^2 instead of n^4.
 *
 * Usage:
 *     echo 'printrestricted(nextprime(10^30))' | gp -q restricted.gp
 */

default(parisizemax, 2^30);

SC_curves = 0;

/* B-smooth part s of N with the rough cofactor r = N/s, by trial division */
smoothpart(N, B) = {
  my(s = 1, r = N);
  forprime(q = 2, B, while(r % q == 0, r /= q; s *= q));
  [s, r];
};

/* Try one random curve E_A over F_p; B is the smoothness bound (= n2 here).
 * Return [A, x0, m, qs] on success (qs = primes of m in (n2, B), empty when
 * B = n2), else 0.  Verbatim from oneshot.gp. */
sc_try(p, B, n2, bound) = {
  my(A = random(p), E = ellinit([0, A, 0, 1, 0], p));
  if(#E == 0, return(0));                                 \\ singular (A == +-2 mod p)
  my(N = ellcard(E), sr = smoothpart(N, B), s = sr[1], r = sr[2]);
  if(s <= bound, return(0));                              \\ smooth factor too small
  my(fs = factor(s)[, 1], P, Q, ord, q, d, fo, fd, lp, qs, Qm);
  for(t = 1, 64,
    P = random(E); Q = ellmul(E, P, r);                   \\ order(Q) divides s
    if(#Q == 1, next);                                    \\ Q = O, resample P
    ord = s;
    for(i = 1, #fs, q = fs[i]; while(ord % q == 0 && #ellmul(E, Q, ord/q) == 1, ord /= q));
    if(ord > bound,
      d = ord; fo = factor(ord)[, 1];                     \\ reduce to minimal smooth order > bound
      forstep(jj = #fo, 1, -1, q = fo[jj]; while(d % q == 0 && d/q > bound, d = d/q));
      fd = factor(d)[, 1]; lp = fd[1];
      if(d >= bound * lp, next);                          \\ enforce strict m < L*r
      qs = select(qq -> qq > n2, fd);
      Qm = ellmul(E, Q, ord/d);                           \\ ord(Qm) = d
      return([A, lift(Qm[1]), d, qs]))
  );
  0;
};

scbound(p) = sqrtint(p) + 1 + sqrtint(4 * sqrtint(p));    \\ integer form of L = (p^{1/4}+1)^2

/* n^2-smooth (k = 0) certificate: identical to oneshot.gp's smoothcert(), but
 * with smoothness bound n^2 instead of n^4, so the q_i list is empty. */
restrictedcert(p) = {
  if(!ispseudoprime(p), error("restricted: p is composite"));
  if(p <= 3, error("restricted: need p > 3"));
  my(n = #binary(p), n2 = n^2, bound = scbound(p), res);
  SC_curves = 0;
  while(1,
    SC_curves++;
    res = sc_try(p, n2, n2, bound);
    if(type(res) == "t_VEC", return([p, res[1], res[2], res[3], res[4]]));
    if(SC_curves % 200 == 0, printf("PROG %d\n", SC_curves)));
};

printrestricted(p) = {                                    \\ "p A x0 m" (+ curve count)
  my(c = restrictedcert(p));
  if(#c[5] != 0, error("nonempty q list!?"));
  printf("CERT %d %d %d %d\n", c[1], c[2], c[3], c[4]);
  printf("CURVES %d\n", SC_curves);
};

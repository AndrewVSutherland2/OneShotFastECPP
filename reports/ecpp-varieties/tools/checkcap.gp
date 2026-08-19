default(parisizemax, 2^30);
smoothpart(N, B) = { my(s=1, r=N); forprime(q=2, B, while(r%q==0, r/=q; s*=q)); [s, r]; };
check(v) = {
  my(p = v[1], ntop = #binary(v[1]), n2top = ntop^2, i = 4, lev = 0);
  my(maxcpl = 0, maxctop = 0, viol = 0, violt = 0);
  printf("p0=%.3g  n=%d\n", 1.0*p, ntop);
  while(i <= #v,
    my(o = v[i], ni = #binary(p), sr = smoothpart(o, n2top), m = sr[1], pn = sr[2]);
    if(pn == 0 || (pn > 1 && !ispseudoprime(pn)), printf("  BAD split\n"));
    my(F = factor(m), Pmax = F[#F[,1], 1]);
    my(cpl = log(m)/log(ni), ctop = log(m)/log(ntop));
    my(okpl = (Pmax <= ni^2), oktop = (Pmax <= n2top));
    printf("  lev %d: n_i=%3d  m: %2d bits, maxprime=%7d  n_i^2=%7d %s   c_perlevel=%.2f c_top=%.2f\n",
           lev, ni, #binary(m), Pmax, ni^2, if(okpl, "ok ", "VIOLATES"), cpl, ctop);
    maxcpl = max(maxcpl, cpl); maxctop = max(maxctop, ctop);
    if(!okpl, viol++); if(!oktop, violt++);
    if(pn == 1, break);
    p = pn; lev++; i += 3);
  [maxcpl, maxctop, viol];
};
{ my(gcpl = 0, gctop = 0, gviol = 0, lines = readstr("ShortPrimalityProofs/certs.csv"));
  for(j = 1, #lines,
    my(v = apply(eval, strsplit(lines[j], ",")));
    my(r = check(v));
    gcpl = max(gcpl, r[1]); gctop = max(gctop, r[2]); gviol += r[3]);
  printf("\nGLOBAL: max c (per-level n_i) = %.2f   max c (top n) = %.2f   levels violating per-level n_i^2-smoothness: %d\n", gcpl, gctop, gviol); }

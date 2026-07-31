# ShortPrimalityProofs
For the purpose of this repository, a **short ECPP** is a sequence of integers $(p_0,A_0,x_0,m_0p_1,A_1,x_1,\ldots,m_{k-1}p_k,A_k,x_k,m_kp_{k+1})$ in which
- we have $k \ge 0$, $p_0 \ge 5$, and put $n:=\lceil \log_2 p_0\rceil$,
- the $p_i$ are odd integers that satisfy $n^2 < p_{i+1} < \sqrt{p_i}$ for $0\le i < k$ and $p_{k+1}=1$,
- the $m_i$ are $n^2$-smooth integers satisfying $L_i < m_ip_{i+1} < r_iL_i$, where $L_i = q_i+1+\lfloor 2\sqrt{q_i}\rfloor$ with $q_i=\lfloor\sqrt{p_i}\rfloor$ and $r_i$ is the least prime divisor of $m_i$,
- each $A_i$ is a nonnegative integer less than $p_i$ with $\gcd(A_i^2-4,p_i)=1$,
- each $x_i$ is a nonnegative integer less than $p_i$,

such that for each $0\le i \le k$ there exist integers $B_i,y_i\in [0,p_i-1]$ with $\gcd(B_i,p_i)=1$ for which $(x_i,y_i)$ is a point of order $m_ip_{i+1}$ on the [Montgomery curve](https://en.wikipedia.org/wiki/Montgomery_curve) $B_iy^2 = x^3 + A_ix^2 +x$ modulo every prime divisor of $p_i$.

A short ECPP proves $p_0$ prime: if $\ell$ is a prime divisor of $p_i$, the reduction of $(x_i,y_i)$ has
order $m_ip_{i+1} > L_i$ in a group of order at most $\ell+1+\lfloor 2\sqrt{\ell}\rfloor$, which is at
most $L_i$ once $\ell\le \sqrt{p_i}$; so every prime divisor of $p_i$ exceeds $\sqrt{p_i}$ and $p_i$ is
prime, by downward induction from $p_{k+1}=1$.  Taking the order modulo every prime divisor of $p_i$ is
essential: the order of a point of $E(\mathbf{Z}/p_0\mathbf{Z})$ is only the least common multiple of its
orders modulo the primes dividing $p_0$, so for a composite $p_0=\ell\ell'$ the required order can be
split between the two factors with every other condition satisfied, as in
```
4410667997551 1365834658413 107710304518 4200232 199129 175565 880
```
where $p_0=2098153\cdot 2102167$ and $(x_0,1)$ has order exactly $8\cdot 525029 = m_0p_1$ in
$E(\mathbf{Z}/p_0\mathbf{Z})$ --- order $8$ modulo the first factor and $525029$ modulo the second ---
together with an honest level for the prime $p_1 = 525029$.  vsmallECPP.py rejects it, and the gcds it
takes along the way expose both factors of $p_0$.

Note that $m_i\ge 2$ always (otherwise $m_ip_{i+1} = p_{i+1} < \sqrt{p_i} < L_i$), so $r_i$ is well
defined, and that a valid $p_{i+1}$ is a prime exceeding $n^2$, hence exactly the $n^2$-rough part of the
product $m_ip_{i+1}$ carried in the certificate: a verifier recovers the $m_i$ and $p_{i+1}$ from the
sequence while trial-dividing only up to $n^2$.  The conditions $p_{i+1}<\sqrt{p_i}$ and
$m_ip_{i+1}<r_iL_i$ keep the certificate to $O(\log p_0)$ bits in $k=O(\log\log p_0)$ levels, and as with
a [one-shot ECPP](https://github.com/AndrewVSutherland/OneShotPrimalityProofs) the whole certificate can
be verified in quasi-quadratic time $O((\log p_0)^{2+o(1)})$.

This repository contains the following resources:
- vsmallECPP.py is a Python program that verifies a short ECPP in quasi-quadratic time.
- short8all.txt contains the 201,072 short ECPPs with $p\le 2^8$; every prime $5\le p\le 2^8$ admits at least one ($p=2$ and $p=3$ admit none, which is why the definition takes $p_0\ge 5$).
- short.gp is a GP script that uses SEA on random curves to search for short ECPPs.
- certs.csv is a list of short ECPPs for the primes listed in the table below.

**Challenge**

Below is a list of short ECPPs for the least prime $p>10^c$ for $c=10,20,\ldots,100$, each found by
<a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> on a single core.  Can you extend this list?

<details>
<summary>$p=10^{10}+19$,&nbsp; <a href="https://math.mit.edu/~drew/">AVS</a> and Claude Code (Fable 5) via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (&lt;1 CPU second, 2 levels).</summary>

```
10000000019 9322349340 1921958667 116108 7217 235 607
```
</details>
<details>
<summary>$p=10^{20}+39$,&nbsp; <a href="https://math.mit.edu/~drew/">AVS</a> and Claude Code (Fable 5) via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~1 CPU seconds, 2 levels).</summary>

```
100000000000000000039 89951393186720294033 26135327929659638076 11876954936 509599 419481 3373
```
</details>
<details>
<summary>$p=10^{30}+57$,&nbsp; <a href="https://math.mit.edu/~drew/">AVS</a> and Claude Code (Fable 5) via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~1 CPU seconds, 2 levels).</summary>

```
1000000000000000000000000000057 106991342299430347297585638871 188675017969977395328028624483 9586166844055967 14077590431530 22429269288900 4934867
```
</details>
<details>
<summary>$p=10^{40}+121$,&nbsp; <a href="https://math.mit.edu/~drew/">AVS</a> and Claude Code (Fable 5) via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~2 CPU seconds, 3 levels).</summary>

```
10000000000000000000000000000000000000121 4902559719197972567355860705693483269960 1108182596968252581482904098615171436619 136850847522421485837 4858820576325 5128710645172 5627651 10966 13232 267
```
</details>
<details>
<summary>$p=10^{50}+151$,&nbsp; <a href="https://math.mit.edu/~drew/">AVS</a> and Claude Code (Fable 5) via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~3 CPU seconds, 3 levels).</summary>

```
100000000000000000000000000000000000000000000000151 4168493225324236537663200519121316886619997624556 45590400574778487393338352639345746341371438526095 13601869893526828282016090 333961040995 170482403318 935683 40410 98277 1193
```
</details>
<details>
<summary>$p=10^{60}+7$,&nbsp; <a href="https://math.mit.edu/~drew/">AVS</a> and Claude Code (Fable 5) via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~3 CPU seconds, 3 levels).</summary>

```
1000000000000000000000000000000000000000000000000000000000007 582594560733647942864167554683101932559817757617443321218616 485107000290978194104100671498974843211121245518759155383025 1591141197950235962428006613308 13059006973177308093777996339 13059586929974475337549368522 1179275706066947 10409174197 23529769559 248836
```
</details>
<details>
<summary>$p=10^{70}+33$,&nbsp; <a href="https://math.mit.edu/~drew/">AVS</a> and Claude Code (Fable 5) via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~45 CPU seconds, 4 levels).</summary>

```
10000000000000000000000000000000000000000000000000000000000000000000033 9193638238367761016751675951961267177328031093306760076142420392370661 5154446894383162465732893340634268054473714420328923353004710402276202 177664059059812136202711786749417049 8649705917060546335363234070525383 42605680596938635666326245539183459 303885057908347530 9987604800409038 5847492647777907 1449221531 3881343 7598118 3432
```
</details>
<details>
<summary>$p=10^{80}+129$,&nbsp; <a href="https://math.mit.edu/~drew/">AVS</a> and Claude Code (Fable 5) via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~180 CPU seconds, 4 levels).</summary>

```
100000000000000000000000000000000000000000000000000000000000000000000000000000129 21738970501572142787334250995459457313706731047873015853200788703257646728746387 41991271598373420878627810537757671546773541549398474815901013810480424900203390 10475053568008371237180368561737331009607 6251483180616232563867118006 2989838890071853567090370520 109497485294994 64373892028 116329932035 967431233 82783 99666 1889
```
</details>
<details>
<summary>$p=10^{90}+289$,&nbsp; <a href="https://math.mit.edu/~drew/">AVS</a> and Claude Code (Fable 5) via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~350 CPU seconds, 3 levels).</summary>

```
1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000289 292390063605290292583797761883066555291812969778434353799388301927489872556173779693570186 201936337656195409865372057596506482368540069004779260915505630735608561247615266092698548 1495231769958646021642577303700692519195356608 1370055486155935701230818919713 631674234039720505250594205971 2824893046457853 65743 243188 50291
```
</details>
<details>
<summary>$p=10^{100}+267$,&nbsp; <a href="https://math.mit.edu/~drew/">AVS</a> and Claude Code (Fable 5) via <a href="https://github.com/AndrewVSutherland/ShortPrimalityProofs/blob/main/short.gp">short.gp</a> (~800 CPU seconds, 3 levels).</summary>

```
10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000267 9922776908302697864551893522633781539036740348871476332653810465626471659307742084010027213820226625 2836186773449440978191705565856816079355043173337102496065404967247933485306358946184015108636857285 159990589984239373593745787852818347617199852606512 989772819515561949937579106033692378969575949199 160869327784807806952128639954291424038666088760 2979396535159070660937591 221054722104 68038539734 625635
```
</details>

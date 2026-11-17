"""Does yafu's ECM find planted factors as often in numbers above 1040 bits as below?
For each size of N = f * G (G a product of two primes out of reach) and each (B1, bits of f): 96 primes f,
one batch of 8 curves each; hits = lanes that found f (yafu -v -v -v), by stage."""
import concurrent.futures as cf, random, re, subprocess, sys, time
Y = sys.argv[1]
FOUND = re.compile(r"found factor (\d+) in stage (\d) in thread \d+, vec position (\d+)")
rng = random.Random(2026)
def is_prime(n):
    if n % 2 == 0: return n == 2
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for _ in range(12):
        a = rng.randrange(2, n - 1); x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True
def prime(bits):
    while True:
        n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(n): return n
def run(args):
    b1, n, f = args
    r = subprocess.run([Y, "-B1ecm", str(b1), "-threads", "1", "-avxecm_stg2", "1", "-v", "-v", "-v", "-aprcl_p", "1", "-logfile", "/dev/null", "-session", "/dev/null"],
                       input="ecm(%d,8)\n" % n, capture_output=True, text=True, cwd="/tmp")
    l1, l2 = set(), set()
    for m in FOUND.finditer(r.stdout.replace("\r", "\n")):
        if int(m.group(1)) % f == 0:
            (l1 if m.group(2) == "1" else l2).add(int(m.group(3)))
    return len(l1), len(l1 | l2)
cells = [(50000, 50), (50000, 62), (1000000, 78)]
fs = {c: [prime(c[1]) for _ in range(96)] for c in cells}
for nbits in (780, 1030, 1130, 1330):
    gb = nbits - 60
    G = prime(gb // 2) * prime(gb - gb // 2)
    for c in cells:
        t = time.time()
        with cf.ThreadPoolExecutor(max_workers=24) as ex:
            res = list(ex.map(run, [(c[0], f * G, f) for f in fs[c]]))
        h1 = sum(r[0] for r in res); h = sum(r[1] for r in res)
        print("N ~%4d bits, B1=%-8d f %2d bits: %3d of 768 curves hit (%3d in stage 1)  rate %.3f   %.0f s" % (nbits, c[0], c[1], h, h1, h / 768.0, time.time() - t), flush=True)

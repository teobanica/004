"""The Schmidt alternative for the circulant graphs (Conjecture 8.4).  Every circulant graph X = Cay(Z_N, S) gets a label:
  Q : G(X) contains 2 non-trivial automorphisms with disjoint supports, so X has quantum symmetry (criterion of Schmidt);
  R : all the relations of the coherent algebra are commutative, by Lemma 6.3, so X has no quantum symmetry;
  P : X is of type E, with N = N1*N2, and Proposition 3.3 reduces the problem to the configurations of the groups
      K1, K2 on Z_N1, Z_N2, which have no quantum symmetry by Lemma 6.3 (applied recursively; not needed for N <= 43);
  M : as R or P, but using at some point the fact that a non-trivial circulant graph with p <= 250 vertices,
      p prime, has no quantum symmetry (Marashdeh);
  ? : undecided.
The conjecture states that the graphs which are not of type Q have no quantum symmetry.
Usage: python3 schmidt.py Nmin Nmax (default 2 43).   Output: schmidt.log (one line per N, appended when Nmin > 2)"""
import sys
from schur import *
from rigid import cayley_relations, all_commutative, orbit_coloring
from bichon import bichon_pair

SETUP = {}
def data(N):
    if N not in SETUP: SETUP[N] = setup((N,))
    return SETUP[N]

CONF = {}
def decide_conf(N, E):
    """The configuration of the orbits of E (a subgroup of Z_N^*, not necessarily containing -1) on Z_N."""
    E = tuple(sorted(E)); key = (N, E)
    if key in CONF: return CONF[key]
    _, neg, sub, cl = data(N)
    col = orbit_coloring(N, E)
    R = cayley_relations(N, sub, col)
    comm = all_commutative(R); pair = bichon_pair(N, sub, col, set()) is not None
    assert not (comm and pair), key                 # consistency of the lemmas with the criterion of Schmidt
    if comm: res = "R"
    elif pair: res = "Q"
    else:
        res = "?"
        for N1 in range(2, N):
            N2 = N // N1
            if N % N1 or N1 > N2 or math.gcd(N1, N2) > 1: continue
            K1 = {a % N1 for a in E if a % N2 == 1}
            K2 = {a % N2 for a in E if a % N1 == 1}
            r = {decide_conf(N1, K1), decide_conf(N2, K2)}
            if "Q" in r: res = "Q"; break               # Proposition 3.3: quantum symmetry of a factor
            if r <= {"R", "P"}: res = "P"; break
            if r <= {"R", "P", "M"}: res = "M"
        if res == "?" and len(E) < N - 1 and (N - 1) in E and N <= 250 and all(N % d for d in range(2, int(N ** 0.5) + 1)):
            res = "M"
    CONF[key] = res
    return res

PART = {}
def decide(N, S):
    _, neg, sub, cl = data(N)
    col = closure(N, sub, coloring(N, S))
    key = (N, tuple(col))
    if key in PART: return PART[key]
    pair = bichon_pair(N, sub, col, set()) is not None
    comm = all_commutative(cayley_relations(N, sub, col))
    assert not (comm and pair), (N, sorted(S))      # consistency of the lemmas with the criterion of Schmidt
    if pair: res = "Q"
    elif comm: res = "R"
    else:
        E = [a for a in range(1, N) if math.gcd(a, N) == 1 and all(col[a * x % N] == col[x] for x in range(N))]
        res = "?"
        if max(orbit_coloring(N, E)) == max(col):       # the coherent algebra is cyclotomic: X is of type E
            res = decide_conf(N, E)
            assert res != "Q"
    PART[key] = res
    return res

def count(N):
    """Number of sets S = -S with each label.  The sets are taken up to multiplication by units, and up to complement."""
    _, neg, sub, cl = data(N)
    n = len(cl); idx = {x: j for j, c in enumerate(cl) for x in c}
    units = [a for a in range(1, N) if math.gcd(a, N) == 1]
    perms = {tuple(idx[a * c[0] % N] for c in cl) for a in units}
    full = (1 << n) - 1
    seen = bytearray(1 << n); tot = {}; und = []
    for mask in range(1 << n):
        if seen[mask]: continue
        orbit = set()
        for p in perms:
            im = 0
            for j in range(n):
                if mask >> j & 1: im |= 1 << p[j]
            orbit.add(im); orbit.add(full ^ im)
        for o in orbit: seen[o] = 1
        S = {x for j, c in enumerate(cl) if mask >> j & 1 for x in c}
        res = decide(N, S)
        tot[res] = tot.get(res, 0) + len(orbit)
        if res == "?": und.append(sorted(S))
    return 1 << n, tot, und

if __name__ == "__main__":
    lo = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    hi = int(sys.argv[2]) if len(sys.argv) > 2 else 43
    out = open("schmidt.log", "w" if lo <= 2 else "a")
    def say(*a):
        print(*a); print(*a, file=out); out.flush()
    if lo <= 2:
        say("CIRCULANT GRAPHS: number of sets S = -S with each label")
        say("   N    sets       Q       R       P       M       ?")
    for N in range(lo, hi + 1):
        sets, tot, und = count(N)
        say("  %2d %7d %7d %7d %7d %7d %7d  %s" % (N, sets, tot.get("Q", 0), tot.get("R", 0), tot.get("P", 0),
            tot.get("M", 0), tot.get("?", 0), ("undecided, up to units and complement: %s" % und) if und else ""))
    out.close()

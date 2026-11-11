"""Proportion of Cayley graphs Cay(A,S), S = -S, which have saturated spectrum (eigenvalues distinct up to
chi <-> chi^-1) and which are of dihedral type (coherent algebra = orbital algebra of A x| {+-1}).
Section 5.   Output: dihedral.log"""
import math, random, sys
from schur import *
random.seed(2026)

def stats(dims, samples):
    N, neg, sub, cl = setup(dims)
    m = len(cl)
    els = list(itertools.product(*[range(n) for n in dims]))
    def val(x, s): return math.cos(2 * math.pi * sum(a * b / n for a, b, n in zip(x, s, dims)))
    tab = [[(val(x, els[c[0]]) + val(x, els[c[1]])) if c[0] != c[1] else val(x, els[c[0]]) for x in els] for c in cl]
    reps = [0] + [c[0] for c in cl]
    exact = 2 ** m - 2 <= samples
    masks = range(1, 2 ** m - 1) if exact else [random.getrandbits(m) for _ in range(samples)]
    sat = dih = tot = 0
    for mask in masks:
        if mask in (0, 2 ** m - 1): continue
        tot += 1
        f = sorted(sum(tab[c][j] for c in range(m) if mask >> c & 1) for j in reps)
        sat += all(f[i + 1] - f[i] > 1e-7 for i in range(len(f) - 1))
        dih += rank(N, sub, sets_from_mask(mask, cl)) == m + 1
    return N, tot, sat, dih, exact

if __name__ == "__main__":
    out = open("dihedral.log", "w")
    def say(*a):
        print(*a); print(*a, file=out); out.flush()
    groups = [((12,), 300), ((16,), 300), ((24,), 300), ((36,), 300), ((60,), 150), ((81,), 60), ((100,), 60), ((128,), 40),
              ((3, 3), 300), ((3, 3, 3), 300), ((3, 3, 3, 3), 60), ((5, 5), 300), ((7, 7), 150), ((9, 9), 60), ((11, 11), 40), ((5, 5, 5), 40),
              ((2,) * 4, 300), ((2,) * 5, 300), ((2,) * 6, 120), ((2,) * 7, 40),
              ((4, 2, 2), 300), ((4, 2, 2, 2), 300), ((4, 2, 2, 2, 2), 60), ((4, 4, 4), 120)]
    for g, n in groups:
        N, tot, sat, dih, exact = stats(g, n)
        say("A=%-16s N=%3d  %-11s saturated spectrum: %5.1f%%   dihedral type: %5.1f%%" % (
            "x".join("Z%d" % k for k in g), N, "(all S)" if exact else "(%d S)" % tot, 100 * sat / tot, 100 * dih / tot))
    out.close()

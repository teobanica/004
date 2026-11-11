"""Circulant graphs with p^n vertices, n >= 2 (Section 4).  For every S = -S we test
  (Q) S meets the units Z_N^* in a union of cosets of the subgroup H of order p  -> quantum symmetry
  (T) Cay(Z_N,S) is of type E, with 1 + p^(n-1) not in E                          -> no quantum symmetry
and check that exactly one of them holds.   Output: primepower.log"""
import random, sys
from schur import *
random.seed(11)

def run(p, n, say, samples=None):
    N = p ** n
    _, neg, sub, cl = setup((N,))
    m = len(cl)
    units = [a for a in range(1, N) if a % p]
    h = p ** (n - 1)
    masks = range(0, 2 ** m) if samples is None else [random.getrandbits(m) for _ in range(samples)]
    cnt = {"Q": 0, "T": 0, "both": 0, "neither": 0}; types = {}
    for mask in masks:
        S = sets_from_mask(mask, cl)
        Su = {s for s in S if s % p}
        Q = all((s + h) % N in Su for s in Su)
        E = [a for a in units if {a * s % N for s in S} == S]
        T = False
        if (1 + h) % N not in E:
            orbits = len({frozenset(a * x % N for a in E) for x in range(1, N)})
            T = rank(N, sub, S) == orbits + 1
            if T: types[len(E)] = types.get(len(E), 0) + 1
        cnt["both" if Q and T else "Q" if Q else "T" if T else "neither"] += 1
    say("N=%3d (%s): quantum symmetry (Q): %d, no quantum symmetry (T): %d, both: %d, neither: %d;  types |E|: %s" % (
        N, "all S" if samples is None else "%d random S" % samples, cnt["Q"], cnt["T"], cnt["both"], cnt["neither"], dict(sorted(types.items()))))

if __name__ == "__main__":
    out = open("primepower.log", "w")
    def say(*a):
        print(*a); print(*a, file=out); out.flush()
    run(2, 2, say); run(2, 3, say); run(3, 2, say); run(2, 4, say); run(5, 2, say); run(3, 3, say)
    run(2, 5, say, 400); run(7, 2, say, 400); run(5, 3, say, 60)
    out.close()

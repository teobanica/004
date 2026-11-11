"""Circulant graphs: G(X) = D_n (order 2n) if and only if X is of dihedral type, for small n (Section 5).
Output: normal.log"""
from schur import *

if __name__ == "__main__":
    out = open("normal.log", "w")
    def say(*a):
        print(*a); print(*a, file=out); out.flush()
    for n in range(5, 15):
        N, neg, sub, cl = setup((n,)); m = len(cl)
        both = autonly = dihonly = neither = 0
        for mask in range(1, 2 ** m - 1):
            S = sets_from_mask(mask, cl)
            dih = rank(N, sub, S) == m + 1
            adj = [[(b - a) % N in S for b in range(N)] for a in range(N)]
            mn = aut_order(N, adj) == 2 * n
            both += dih and mn; autonly += mn and not dih; dihonly += dih and not mn; neither += not (dih or mn)
        say("n=%2d: |G|=2n and dihedral type: %3d;  |G|=2n only: %d;  dihedral type only: %d;  neither: %3d" % (n, both, autonly, dihonly, neither))
    out.close()

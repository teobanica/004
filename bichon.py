"""The quantum group of Bichon G*(X), between G(X) and G+(X) (Section 7).
For a circulant graph X = Cay(Z_N, S) the script tries to decide whether G*(X) = G(X), as follows:
  C (classical): X complete, or X empty with N <= 3, or reduction to a smaller graph when the complement
     is disconnected (Proposition 7.3), or all the relations are commutative, by starting with the relations
     contained in S, and by applying Lemma 6.3.
  Q (quantum): there are 2 non-trivial automorphisms with disjoint supports, and no edges between the supports
     (Proposition 7.1), found by backtracking.
  ? : undecided.  Conjecture 7.4 predicts C for these graphs.
Also checks Theorem 7.2 (dihedral type) on the groups Z_4 x Z_2^m.   Output: bichon.log"""
import sys, random
from schur import *
from rigid import cayley_relations, commutative_relations

def all_commutative(N, sub, col, S):
    """Are all the relations commutative for G*(X)?  Starts with the relations contained in S (Bichon),
    and applies Lemma 6.3."""
    R = cayley_relations(N, sub, col)
    return len(commutative_relations(R, {col[s] for s in S})) == len(set(col))

def bichon_pair(N, sub, col, S):
    """Searches for automorphisms sigma, tau with disjoint supports and no edges between the supports.
    By translating we can assume that sigma fixes 0.  Returns the support of sigma, as a bitmask, when such a pair
    exists, or None.  The existence of tau, whose support must avoid the support of sigma and its neighbors, is
    checked by tau_inside, and tau itself is not returned."""
    R = [[col[sub[y][x]] for y in range(N)] for x in range(N)]
    full = (1 << N) - 1
    nb = [0] * N
    for x in range(N):
        for y in range(N):
            if sub[y][x] in S: nb[x] |= 1 << y
    cache = {}
    def tau_inside(F):
        """Is there a non-trivial automorphism whose support is contained in the set F?"""
        if F in cache: return cache[F]
        vs = [v for v in range(N) if F >> v & 1]
        img = list(range(N)); out = [v for v in range(N) if not F >> v & 1]
        def rec(i, used, moved):
            if i == len(vs): return moved
            v = vs[i]
            for w in vs:
                if used >> w & 1: continue
                if any(R[v][o] != R[w][o] for o in out): continue
                if any(R[v][vs[j]] != R[w][img[vs[j]]] for j in range(i)): continue
                img[v] = w
                if rec(i + 1, used | 1 << w, moved or w != v): return True
            return False
        cache[F] = rec(0, 0, False)
        return cache[F]
    img = [0] * N
    def rec(i, used, M):
        if M:
            zone = M
            for v in range(N):
                if M >> v & 1: zone |= nb[v]
            free = full & ~zone
            if bin(free).count("1") < 2: return None
        if i == N:
            if M and tau_inside(free): return M
            return None
        for w in [i] + [w for w in range(1, N) if w != i]:
            if used >> w & 1: continue
            if any(R[i][j] != R[w][img[j]] for j in range(i)): continue
            img[i] = w
            r = rec(i + 1, used | 1 << w, M | ((1 << i) | (1 << w) if w != i else 0))
            if r: return r
        return None
    return rec(1, 1, 0)

def decide(N, S):
    """Returns (status, reason) for the circulant graph Cay(Z_N, S)."""
    S = set(S)
    if len(S) == N - 1: return "C", "complete"
    if not S: return ("C" if N <= 3 else "Q"), "empty"
    T = set(range(1, N)) - S
    d = math.gcd(N, *T)
    if d > 1:                                   # complement disconnected: Proposition 7.3
        st, why = decide(N // d, {s // d for s in S if s % d == 0})
        return st, "join"
    _, neg, sub, cl = setup((N,))
    col = closure(N, sub, coloring(N, S))
    c = all_commutative(N, sub, col, S)
    q = bichon_pair(N, sub, col, S) is not None
    assert not (c and q), (N, sorted(S))
    return ("C", "relations") if c else ("Q", "pair") if q else ("?", "")

def symmetric_sets(cl):
    for mask in range(1 << len(cl)):
        yield {x for j, c in enumerate(cl) if mask >> j & 1 for x in c}

if __name__ == "__main__":
    P = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    out = open("bichon.log", "w")
    def say(*a):
        print(*a); print(*a, file=out); out.flush()
    say("CIRCULANT GRAPHS with N <= %d vertices: is G*(X) = G(X)?" % P)
    say("   N  sets     C     Q     ?")
    und = []
    for N in range(2, P + 1):
        _, neg, sub, cl = setup((N,))
        cnt = {"C": 0, "Q": 0, "?": 0}
        for S in symmetric_sets(cl):
            st, why = decide(N, S)
            cnt[st] += 1
            if st == "?": und.append((N, sorted(S)))
        say("  %2d %5d %5d %5d %5d" % (N, 1 << len(cl), cnt["C"], cnt["Q"], cnt["?"]))
    say("  undecided graphs:", und if und else "none")
    say("THEOREM 7.5, graphs with p^2 vertices having quantum symmetry (S meets the units in a union of cosets of pZ):")
    for p in (2, 3):
        N = p * p
        _, neg, sub, cl = setup((N,))
        for S in symmetric_sets(cl):
            U = {s for s in S if s % p}
            if any((s + p) % N not in U for s in U): continue
            say("  N=%d S=%-22s units in S: %-5s  %s" % (N, sorted(S), len(U) == N - p, decide(N, S)))
    say("THEOREM 7.2, dihedral type on Z_4 x Z_2^m: all the relations commutative for G*(X)?")
    random.seed(1)
    for dims, samples in [((4,), 0), ((4, 2), 0), ((4, 2, 2), 0), ((4, 2, 2, 2), 400)]:
        N, neg, sub, cl = setup(dims)
        want = 1 + len(cl); tot = good = 0; exc = []
        sets = symmetric_sets(cl) if not samples else \
            ({x for c in cl if random.random() < 0.5 for x in c} for _ in range(samples))
        for S in sets:
            if not S: continue
            col = closure(N, sub, coloring(N, S))
            if max(col) + 1 != want: continue
            tot += 1
            if all_commutative(N, sub, col, S): good += 1
            else: exc.append(sorted(S))
        say("  A = Z_%s: graphs of dihedral type %d (%s), with G* = G proved %d, exceptions %s" %
            ("x".join(map(str, dims)), tot, "sampled" if samples else "all", good, exc))
    out.close()

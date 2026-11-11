"""Transitive permutation groups (Theorem 8.1).  For a group G acting on the cosets of a subgroup H, computes
N = |G/H|, k = |H|, the maximal number m of fixed points of the non-trivial elements, the orbital configuration,
and checks that all its relations are commutative (Lemma 6.3) when N > 2(k-1)(km-1) + 2.
Also searches for a graph (union of orbitals) whose coherent algebra is the commutant of G.   Output: transitive.log"""
import random
from rigid import commutative_relations

def mul(p, q): return tuple(p[i] for i in q)                      # (p*q)(x) = p(q(x))

def generate(gens):
    e = tuple(range(len(gens[0]))); els = {e}; todo = [e]
    while todo:
        g = todo.pop()
        for s in gens:
            h = mul(s, g)
            if h not in els: els.add(h); todo.append(h)
    return els

def coset_action(gens, hgens):
    G = generate(gens); H = generate(hgens)
    rep = {}; cosets = []
    for g in sorted(G):
        if g in rep: continue
        for h in H: rep[mul(g, h)] = len(cosets)
        cosets.append(g)
    N = len(cosets)
    act = {g: tuple(rep[mul(g, c)] for c in cosets) for g in G}
    return G, H, N, act

def orbitals(N, act, gens):
    R = [[-1] * N for _ in range(N)]; c = 0
    for k in range(N):
        if R[0][k] >= 0: continue
        todo = [(0, k)]; R[0][k] = c
        while todo:
            i, j = todo.pop()
            for s in gens:
                a, b = act[s][i], act[s][j]
                if R[a][b] < 0: R[a][b] = c; todo.append((a, b))
        c += 1
    return R, c

def wl(N, C):
    """Weisfeiler-Leman on a colored complete directed graph."""
    while True:
        r = max(max(row) for row in C) + 1
        sig = {}
        new = [[0] * N for _ in range(N)]
        for i in range(N):
            Ci = C[i]
            for j in range(N):
                cnt = {}
                for v in range(N):
                    key = Ci[v] * r + C[v][j]
                    cnt[key] = cnt.get(key, 0) + 1
                s = (Ci[j], tuple(sorted(cnt.items())))
                new[i][j] = sig.setdefault(s, len(sig))
        if len(sig) == r: return r
        C = new

def study(name, gens, hgens, say):
    G, H, N, act = coset_action(gens, hgens)
    k = len(H)
    e = tuple(range(len(gens[0])))
    m = max(sum(1 for i in range(N) if act[g][i] == i) for g in G if g != e)
    R, rank = orbitals(N, act, gens)
    bound = 2 * (k - 1) * (k * m - 1) + 2
    comm = len(commutative_relations(R)) == rank
    assert comm or N <= bound
    # a graph with this coherent algebra: union of orbitals, closed under transposition
    tr = {R[i][j]: R[j][i] for i in range(N) for j in range(N)}
    random.seed(1); found = None
    for attempt in range(12):
        pick = {r for r in range(1, rank) if random.random() < 0.5}
        pick |= {tr[r] for r in pick}
        if not pick or len(pick) == rank - 1: continue
        C = [[0 if i == j else (1 if R[i][j] in pick else 2) for j in range(N)] for i in range(N)]
        if wl(N, C) == rank:
            found = sum(1 for j in range(N) if R[0][j] in pick); break
    say("%-28s |G|=%4d N=%3d k=%d m=%d rank=%2d bound=%3d applies=%-5s all commutative=%-5s graph with this algebra: %s" % (
        name, len(G), N, k, m, rank, bound, N > bound, comm, "valency %d" % found if found else "not found"))

def perm(n, *cycles):
    p = list(range(n))
    for c in cycles:
        for a, b in zip(c, c[1:] + c[:1]): p[a] = b
    return tuple(p)

if __name__ == "__main__":
    out = open("transitive.log", "w")
    def say(*a):
        print(*a); print(*a, file=out); out.flush()
    # PSL(2,7) on the projective line 0..6, infinity = 7
    t = tuple([(x + 1) % 7 for x in range(7)] + [7])
    inv = {0: 7, 7: 0}
    for x in range(1, 7): inv[x] = (-pow(x, 5, 7)) % 7
    s = tuple(inv[x] for x in range(8))
    study("PSL(2,7) / Z_2", [t, s], [s], say)
    study("PSL(2,7) / Z_3", [t, s], [mul(s, t)], say)
    a5 = [perm(5, (0, 1, 2, 3, 4)), perm(5, (0, 1, 2))]
    study("A_5 / Z_2", a5, [perm(5, (0, 1), (2, 3))], say)
    study("A_5 / Z_3", a5, [perm(5, (0, 1, 2))], say)
    s5 = [perm(5, (0, 1, 2, 3, 4)), perm(5, (0, 1))]
    study("S_5 / <(12)>", s5, [perm(5, (0, 1))], say)
    study("S_5 / <(12)(34)>", s5, [perm(5, (0, 1), (2, 3))], say)
    study("S_5 / <(123)>", s5, [perm(5, (0, 1, 2))], say)
    study("S_4 / Z_2 (transposition)", [perm(4, (0, 1, 2, 3)), perm(4, (0, 1))], [perm(4, (0, 1))], say)
    out.close()

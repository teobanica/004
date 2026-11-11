"""Examples quoted in the paper: type and automorphism group of some explicit graphs.   Output: examples.log"""
from schur import *

def check(N, S, E, say, name):
    _, neg, sub, cl = setup((N,))
    orbits = len({frozenset(a * x % N for a in E) for x in range(1, N)})
    r = rank(N, sub, S)
    if N > 70:                                  # too big for the naive automorphism count
        say("%s: rank of coherent algebra = %d, number of orbitals of Z_N x| E = %d" % (name, r, orbits + 1)); return
    adj = [[(b - a) % N in S for b in range(N)] for a in range(N)]
    say("%s: rank of coherent algebra = %d, number of orbitals of Z_N x| E = %d, |Aut| = %d, N|E| = %d" % (
        name, r, orbits + 1, aut_order(N, adj), N * len(E)))

if __name__ == "__main__":
    out = open("examples.log", "w")
    def say(*a):
        print(*a); print(*a, file=out); out.flush()
    check(25, {1, 7, 18, 24}, [1, 7, 18, 24], say, "Cay(Z_25, {+-1,+-7})")
    check(49, {1, 18, 19, 30, 31, 48}, [1, 18, 19, 30, 31, 48], say, "Cay(Z_49, {+-1,+-18,+-19})")
    check(65, {1, 8, 57, 64}, [1, 8, 57, 64], say, "Cay(Z_65, {+-1,+-8})")
    check(93, {1, 25, 26, 67, 68, 92}, [1, 25, 26, 67, 68, 92], say, "Cay(Z_93, {+-1,+-25,+-26})")
    check(6, {1, 3, 5}, [1, 5], say, "K_33 = Cay(Z_6, {1,3,5})")
    check(4, {1, 3}, [1, 3], say, "C_4")
    out.close()

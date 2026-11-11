# Cayley graphs of abelian groups and their quantum symmetries

**Teo Banica and Claude Fable**, November 2026

A graph $X$ has no quantum symmetry when its quantum automorphism group $G^+(X)$ coincides with the usual automorphism group $G(X)$. Building on the recent work of Marashdeh on the circulant graphs of prime order, this paper studies the Cayley graphs of the finite abelian groups $A$, by using the Fourier transform over $A$ and the coherent algebra of the graph, and then more general graphs, by a method of rigid triangles. Among the results:

- if the coherent algebra is the orbital algebra of $A\rtimes E$, with $E$ of order $k$, whose non-trivial elements have at most $m$ fixed points, the graph has no quantum symmetry when $|A|\ge k^2m$;
- the circulant graphs with $p^n$ vertices having quantum symmetry are classified;
- for $E=\pm1$ the graph has quantum symmetry exactly when $A=\mathbb Z_4\times\mathbb Z_2^n$, and almost all circulant graphs have no quantum symmetry;
- a similar bound, $N\ge 2k^2m$, holds for the graphs with $N$ vertices coming from transitive permutation groups;
- conjecturally, a circulant graph has quantum symmetry exactly when it has 2 non-trivial automorphisms with disjoint supports, and this is verified up to 43 vertices.

## Contents

- `004.pdf` — the paper (30 pages)
- `004.tex` — its LaTeX source (self-contained; compile with `pdflatex`, twice)
- `code/` — every script behind the computational statements, with their output logs. `code/README.md` maps the scripts to the sections and statements of the paper.

## Running the code

Python 3, standard library only. Each script explains what it does in its first lines, and writes a log file next to itself. The computations behind the statements of the paper are exact, in integer arithmetic. Most scripts run in seconds to minutes, the longest job being the census of `schmidt.py` up to 43 vertices.

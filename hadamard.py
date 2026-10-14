"""Check and build Hadamard matrices in pure Python (no dependencies).

A Hadamard matrix of order n is an n x n matrix with entries +1/-1 whose rows are
pairwise orthogonal, i.e. H H^T = n I.

Commands:
    python hadamard.py check FILE   verify a matrix (rows of +/- characters, or 1/-1 tokens)
    python hadamard.py gs FILE      build and verify a Goethals-Seidel matrix from four
                                    +/- sequences of equal length n (order 4n)
    python hadamard.py gaps N       multiples of 4 up to N that Sylvester, Paley I, Paley II
                                    and Kronecker products do NOT reach
    python hadamard.py build N      print a Hadamard matrix of order N built from those
                                    elementary constructions, if they reach N

By MV2 of Munim, Inc. Made with AI agents working in MV2's name. Apache-2.0 licence.
"""

import argparse
import sys


# ---------------------------------------------------------------- verification

def _row_mask(row):
    """Encode a +1/-1 row as an int whose set bits mark the -1 entries."""
    mask = 0
    for j, x in enumerate(row):
        if x == -1:
            mask |= 1 << j
    return mask


def check(rows):
    """Return (ok, message) for a candidate Hadamard matrix given as a list of rows."""
    n = len(rows)
    if n == 0:
        return False, "empty matrix"
    for i, row in enumerate(rows):
        if len(row) != n:
            return False, "row %d has length %d, expected %d (matrix must be square)" % (i, len(row), n)
        for j, x in enumerate(row):
            if x != 1 and x != -1:
                return False, "entry (%d, %d) is %r, expected +1 or -1" % (i, j, x)
    masks = [_row_mask(r) for r in rows]
    for i in range(n):
        mi = masks[i]
        for j in range(i + 1, n):
            # inner product = (#agreements) - (#disagreements) = n - 2 * popcount(xor)
            dot = n - 2 * bin(mi ^ masks[j]).count("1")
            if dot != 0:
                return False, "rows %d and %d have inner product %d, expected 0" % (i, j, dot)
    return True, "Hadamard matrix of order %d (H H^T = %d I)" % (n, n)


def is_hadamard(rows):
    return check(rows)[0]


# ---------------------------------------------------------------- basic constructions

def kron(a, b):
    """Kronecker product of two matrices given as lists of rows."""
    m = len(b)
    out = []
    for arow in a:
        for k in range(m):
            brow = b[k]
            out.append([x * y for x in arow for y in brow])
    return out


H2 = [[1, 1], [1, -1]]


def sylvester(k):
    """Sylvester's Hadamard matrix of order 2**k."""
    h = [[1]]
    for _ in range(k):
        h = kron(h, H2)
    return h


def prime_power(q):
    """Return (p, k) with q == p**k and p prime, or None."""
    if q < 2:
        return None
    p = 2
    while p * p <= q:
        if q % p == 0:
            break
        p += 1
    else:
        return (q, 1)
    k = 0
    while q % p == 0:
        q //= p
        k += 1
    return (p, k) if q == 1 else None


class _Field:
    """GF(p**k); elements are ints 0..q-1 whose base-p digits are polynomial coefficients."""

    def __init__(self, q):
        pk = prime_power(q)
        if pk is None:
            raise ValueError("%d is not a prime power" % q)
        self.q = q
        self.p, self.k = pk
        powers = self._primitive_powers()
        self.squares = set(powers[0::2])

    def _digits(self, v):
        p = self.p
        return [(v // p ** i) % p for i in range(self.k)]

    def _code(self, digits):
        return sum(c * self.p ** i for i, c in enumerate(digits))

    def _primitive_powers(self):
        """Powers 1, x, x^2, ..., x^(q-2) of x modulo a monic f of degree k in which x
        has order q - 1. Such an f is irreducible and x is a primitive element."""
        p, k, q = self.p, self.k, self.q
        one = [1] + [0] * (k - 1)
        for low in range(q):
            f_low = self._digits(low)  # f = x^k + f_low[k-1] x^(k-1) + ... + f_low[0]
            if f_low[0] == 0:
                continue
            powers = []
            cur = one
            ok = True
            for i in range(q - 1):
                powers.append(self._code(cur))
                top = cur[-1]
                shifted = [0] + cur[:-1]
                cur = [(shifted[t] - top * f_low[t]) % p for t in range(k)]
                if cur == one and i != q - 2:
                    ok = False
                    break
            if ok and cur == one:
                return powers
        raise ArithmeticError("no primitive polynomial found for q = %d" % q)

    def sub(self, a, b):
        if self.k == 1:
            return (a - b) % self.p
        da, db = self._digits(a), self._digits(b)
        return self._code([(x - y) % self.p for x, y in zip(da, db)])

    def chi(self, a):
        """Quadratic character: 0 at 0, +1 on nonzero squares, -1 otherwise."""
        if a == 0:
            return 0
        return 1 if a in self.squares else -1


def jacobsthal(q):
    """Jacobsthal matrix Q[a][b] = chi(a - b) over GF(q), q an odd prime power."""
    if q % 2 == 0:
        raise ValueError("q must be odd")
    f = _Field(q)
    return [[f.chi(f.sub(a, b)) for b in range(q)] for a in range(q)]


def paley1(q):
    """Paley I: Hadamard matrix of order q + 1 for a prime power q = 3 (mod 4)."""
    if q % 4 != 3 or prime_power(q) is None:
        raise ValueError("Paley I needs a prime power q = 3 (mod 4), got %d" % q)
    Q = jacobsthal(q)
    n = q + 1
    h = [[1] * n]
    for i in range(q):
        row = [-1] + [Q[i][j] + (1 if i == j else 0) for j in range(q)]
        h.append(row)
    return h


def paley2(q):
    """Paley II: Hadamard matrix of order 2(q + 1) for a prime power q = 1 (mod 4)."""
    if q % 4 != 1 or prime_power(q) is None:
        raise ValueError("Paley II needs a prime power q = 1 (mod 4), got %d" % q)
    Q = jacobsthal(q)
    m = q + 1
    # symmetric conference matrix C = [[0, 1^T], [1, Q]]
    c = [[0] + [1] * q] + [[1] + Q[i] for i in range(q)]
    zero_block = [[1, -1], [-1, -1]]
    h = []
    for i in range(m):
        for r in range(2):
            row = []
            for j in range(m):
                if c[i][j] == 0:
                    row.extend(zero_block[r])
                else:
                    row.extend([c[i][j] * x for x in H2[r]])
            h.append(row)
    return h


# ---------------------------------------------------------------- Goethals-Seidel

def periodic_autocorrelation(seq, s):
    n = len(seq)
    return sum(seq[i] * seq[(i + s) % n] for i in range(n))


def goethals_seidel(a, b, c, d):
    """Goethals-Seidel array of order 4n from four +/-1 sequences of length n whose
    periodic autocorrelations sum to 0 at every nonzero shift."""
    n = len(a)
    if not (len(b) == len(c) == len(d) == n) or n == 0:
        raise ValueError("the four sequences must have the same positive length")
    for s in range(1, n):
        total = sum(periodic_autocorrelation(x, s) for x in (a, b, c, d))
        if total != 0:
            raise ValueError("periodic autocorrelations sum to %d at shift %d, need 0" % (total, s))

    def circ(x):
        return [[x[(j - i) % n] for j in range(n)] for i in range(n)]

    def times_r(m):  # M R, with R the back-diagonal permutation
        return [[row[n - 1 - j] for j in range(n)] for row in m]

    def transpose(m):
        return [list(col) for col in zip(*m)]

    def neg(m):
        return [[-v for v in row] for row in m]

    A, B, C, D = circ(a), circ(b), circ(c), circ(d)
    BR, CR, DR = times_r(B), times_r(C), times_r(D)
    BtR, CtR, DtR = times_r(transpose(B)), times_r(transpose(C)), times_r(transpose(D))
    blocks = [
        [A, BR, CR, DR],
        [neg(BR), A, DtR, neg(CtR)],
        [neg(CR), neg(DtR), A, BtR],
        [neg(DR), CtR, neg(BtR), A],
    ]
    h = []
    for brow in blocks:
        for i in range(n):
            h.append([v for blk in brow for v in blk[i]])
    return h


# ---------------------------------------------------------------- elementary orders

def elementary_recipes(limit):
    """Orders <= limit reached by Sylvester, Paley I, Paley II and Kronecker products,
    each mapped to one recipe that builds it."""
    recipes = {1: ("one",), 2: ("h2",)}
    for q in range(3, limit + 1, 2):
        if prime_power(q) is None:
            continue
        if q % 4 == 3 and q + 1 <= limit:
            recipes.setdefault(q + 1, ("paley1", q))
        if q % 4 == 1 and 2 * (q + 1) <= limit:
            recipes.setdefault(2 * (q + 1), ("paley2", q))
    changed = True
    while changed:
        changed = False
        known = sorted(recipes)
        for x in known:
            if x == 1:
                continue
            for y in known:
                if y == 1 or x * y > limit:
                    continue
                if x * y not in recipes:
                    recipes[x * y] = ("kron", x, y)
                    changed = True
    return recipes


def elementary_gaps(limit):
    recipes = elementary_recipes(limit)
    return [n for n in range(4, limit + 1, 4) if n not in recipes]


def build(order, _recipes=None):
    """Hadamard matrix of the given order from the elementary constructions."""
    recipes = _recipes if _recipes is not None else elementary_recipes(order)
    if order not in recipes:
        raise ValueError("order %d is not reached by Sylvester, Paley or Kronecker products" % order)
    r = recipes[order]
    if r[0] == "one":
        return [[1]]
    if r[0] == "h2":
        return [row[:] for row in H2]
    if r[0] == "paley1":
        return paley1(r[1])
    if r[0] == "paley2":
        return paley2(r[1])
    return kron(build(r[1], recipes), build(r[2], recipes))


# ---------------------------------------------------------------- input / output

def parse_rows(text):
    """Parse rows of '+'/'-' characters, or whitespace/comma separated 1/-1 tokens.
    Blank lines and lines starting with '#' are skipped. A single line whose length is
    a perfect square is split into a square matrix."""
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if any(ch in line for ch in " \t,;"):
            tokens = line.replace(",", " ").replace(";", " ").split()
            row = []
            for t in tokens:
                if t in ("1", "+1", "+"):
                    row.append(1)
                elif t in ("-1", "-"):
                    row.append(-1)
                else:
                    raise ValueError("unexpected token %r" % t)
        else:
            row = []
            for ch in line:
                if ch == "+":
                    row.append(1)
                elif ch == "-":
                    row.append(-1)
                else:
                    raise ValueError("unexpected character %r (use + and -)" % ch)
        rows.append(row)
    if len(rows) == 1:
        flat = rows[0]
        n = int(round(len(flat) ** 0.5))
        if n > 1 and n * n == len(flat):
            rows = [flat[i * n:(i + 1) * n] for i in range(n)]
    return rows


def format_rows(rows):
    return "\n".join("".join("+" if x == 1 else "-" for x in row) for row in rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd")
    p_check = sub.add_parser("check", help="verify a matrix file")
    p_check.add_argument("file")
    p_gs = sub.add_parser("gs", help="build and verify a Goethals-Seidel matrix from 4 sequences")
    p_gs.add_argument("file")
    p_gs.add_argument("--print", action="store_true", help="print the matrix")
    p_gaps = sub.add_parser("gaps", help="orders not reached by elementary constructions")
    p_gaps.add_argument("limit", type=int)
    p_build = sub.add_parser("build", help="print an elementary Hadamard matrix")
    p_build.add_argument("order", type=int)
    args = parser.parse_args(argv)

    if args.cmd == "check":
        with open(args.file, encoding="utf-8") as fh:
            rows = parse_rows(fh.read())
        ok, msg = check(rows)
        print(("OK: " if ok else "FAIL: ") + msg)
        return 0 if ok else 1
    if args.cmd == "gs":
        with open(args.file, encoding="utf-8") as fh:
            seqs = parse_rows(fh.read())
        if len(seqs) != 4:
            print("FAIL: expected 4 sequences, got %d" % len(seqs))
            return 1
        try:
            h = goethals_seidel(*seqs)
        except ValueError as e:
            print("FAIL: " + str(e))
            return 1
        ok, msg = check(h)
        print(("OK: " if ok else "FAIL: ") + msg)
        if args.print:
            print(format_rows(h))
        return 0 if ok else 1
    if args.cmd == "gaps":
        gaps = elementary_gaps(args.limit)
        print("%d multiples of 4 up to %d not reached by Sylvester/Paley/Kronecker:" % (len(gaps), args.limit))
        print(" ".join(str(g) for g in gaps))
        return 0
    if args.cmd == "build":
        try:
            h = build(args.order)
        except ValueError as e:
            print("FAIL: " + str(e))
            return 1
        print(format_rows(h))
        return 0
    parser.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())

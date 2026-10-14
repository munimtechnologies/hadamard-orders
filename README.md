# The Hadamard conjecture in October 2026

**What is proved, what the August 2026 constructions settled, and what is still open.** This is an
exposition with complete proofs of the classical results and a small, tested Python companion
(`hadamard.py`) that checks any claimed Hadamard matrix and rebuilds the classical constructions.

By MV2 of Munim, Inc. This article and the code were researched and written with AI agents working in
MV2's name. Every statement below is either proved here or cited to a source; section 7 separates
what is proved from what is open, and section 8 lists what is new here. Corrections are welcome by
issue or pull request.

Version 1.0, 2026-10-14. Text CC BY 4.0, code MIT.

Related problem records: Justin Sun Prize JSP-000013 (MathDB #315867), Epoch AI's FrontierMath
open-problem page for order 668, OEIS A007299.

---

## 1. The question

A **Hadamard matrix** of order *n* is an *n* × *n* matrix *H* with entries ±1 whose rows are pairwise
orthogonal:

    H Hᵀ = n I.

Examples: (1), and

    H₂ = [ 1   1 ]
         [ 1  -1 ]

**Hadamard conjecture.** For every positive integer *k* there is a Hadamard matrix of order 4*k*.

Hadamard (1893) met these matrices as the extremal case of his determinant inequality (section 2) and
showed that their order must be 1, 2 or a multiple of 4 (section 3). The conjecture that the converse
holds carries his name, although it is often attributed to Paley (1933). It is still open: nobody knows a proof that every multiple of 4 is an order, and infinitely
many multiples of 4 still have no known construction (section 6).

## 2. Why these matrices: Hadamard's determinant bound

**Theorem (Hadamard 1893).** If *M* is a real *n* × *n* matrix with all |*m*ᵢⱼ| ≤ 1, then
|det *M*| ≤ *n*^{n/2}. A ±1 matrix attains the bound exactly when it is a Hadamard matrix.

*Proof.* |det *M*| is the volume of the parallelepiped spanned by the rows *r*₁, ..., *r*ₙ. Gram–Schmidt
writes that volume as ∏ ‖*r*ᵢ*\**‖, where *r*ᵢ*\** is the component of *r*ᵢ orthogonal to the earlier rows,
so |det *M*| ≤ ∏ ‖*r*ᵢ‖ ≤ (√*n*)ⁿ. For a ±1 matrix ‖*r*ᵢ‖ = √*n* exactly, and equality in the first
step holds exactly when every row is orthogonal to the earlier ones. ∎

So the conjecture asks: in which dimensions can the maximal determinant *n*^{n/2} actually be
reached? The same matrices are the best possible two-level experimental designs and the basis of
the code Mariner 9 used to send images from Mars, which is why the question has been studied far beyond pure
combinatorics.

## 3. The order must be 1, 2 or a multiple of 4

Multiplying rows or columns by −1 and permuting them preserves the Hadamard property, so we may assume
the first row is all +1.

**Theorem.** If a Hadamard matrix of order *n* ≥ 3 exists, then 4 divides *n*.

*Proof.* Take the first three rows (first row all +1). Sort the columns by the signs of rows 2 and 3:
let *a*, *b*, *c*, *d* count the columns with signs (+,+), (+,−), (−,+), (−,−). Orthogonality gives

- row 1 · row 2 = 0: *a* + *b* = *c* + *d*
- row 1 · row 3 = 0: *a* + *c* = *b* + *d*
- row 2 · row 3 = 0: *a* + *d* = *b* + *c*

Subtracting the first two equations gives *b* = *c*, then *a* = *d*, and the third gives *a* = *b*. So
*a* = *b* = *c* = *d* and *n* = 4*a*. ∎

## 4. Constructions that reach most orders, with proofs

### 4.1 Kronecker products (Sylvester 1867)

If *A* and *B* are Hadamard of orders *m* and *n*, then (*A* ⊗ *B*)(*A* ⊗ *B*)ᵀ = *A Aᵀ* ⊗ *B Bᵀ* =
*mI* ⊗ *nI* = *mn I*, so *A* ⊗ *B* is Hadamard of order *mn*. With *H*₂ this gives every power of 2
(Sylvester's matrices), and in general the set of Hadamard orders is closed under multiplication.

### 4.2 Paley's constructions (Paley 1933)

Let *q* be an odd prime power and χ the quadratic character of the field GF(*q*): χ(0) = 0, χ(*x*) = 1
if *x* is a nonzero square, χ(*x*) = −1 otherwise. χ is multiplicative, and Σₓ χ(*x*) = 0 because half
of the nonzero elements are squares.

**Lemma.** For *a* ≠ *b* in GF(*q*), Σₓ χ(*x* − *a*) χ(*x* − *b*) = −1.

*Proof.* The term *x* = *b* is 0. For *x* ≠ *b*, χ(*x* − *a*)χ(*x* − *b*) = χ((*x* − *a*)/(*x* − *b*)),
because χ(*y*) = χ(*y*⁻¹). The map *x* ↦ (*x* − *a*)/(*x* − *b*) is a bijection from GF(*q*) ∖ {*b*} onto
GF(*q*) ∖ {1}, so the sum is Σ_{y≠1} χ(*y*) = 0 − χ(1) = −1. ∎

Let *Q* be the *q* × *q* **Jacobsthal matrix** *Q*ₐᵦ = χ(*a* − *b*). The lemma and Σχ = 0 give

    Q Qᵀ = q I − J,    Q J = J Q = 0,

where *J* is the all-ones matrix. Since χ(−1) = (−1)^{(q−1)/2}, *Q* is antisymmetric when *q* ≡ 3 (mod 4)
and symmetric when *q* ≡ 1 (mod 4).

**Paley I (*q* ≡ 3 mod 4): order *q* + 1.** Put *S* = [[0, **1**ᵀ], [−**1**, *Q*]] (bordered by a row of
+1 and a column of −1). Row by row, *S Sᵀ* = *q I*: the border row has *q* ones; it is orthogonal to the
other rows because the rows of *Q* sum to 0; and two other rows *i*, *j* give 1 + (*q* δᵢⱼ − 1). *S* is
antisymmetric, so with *H* = *I* + *S*,

    H Hᵀ = (I + S)(I − S) = I − S² = I + S Sᵀ = (q + 1) I.

The diagonal of *Q* is χ(0) = 0, so *H* has entries ±1: a Hadamard matrix of order *q* + 1.

**Paley II (*q* ≡ 1 mod 4): order 2(*q* + 1).** Now *C* = [[0, **1**ᵀ], [**1**, *Q*]] is symmetric with
*C Cᵀ* = *q I* (same row computation) and zero diagonal. Let *A* = [[1, 1], [1, −1]] and
*B* = [[1, −1], [−1, −1]]. Then *A Aᵀ* = *B Bᵀ* = 2*I* and *A Bᵀ* + *B Aᵀ* = 0, so

    H = C ⊗ A + I ⊗ B,    H Hᵀ = 2q I + 2 I + C ⊗ (A Bᵀ + B Aᵀ) = 2(q + 1) I,

and *H* has entries ±1 because *C* is ±1 off the diagonal and 0 on it.

### 4.3 Where the easy constructions stop

Close {1, 2} ∪ {*q* + 1 : *q* ≡ 3 mod 4} ∪ {2(*q* + 1) : *q* ≡ 1 mod 4} (prime powers *q*) under
multiplication. Every multiple of 4 up to 88 is reached, for instance 28 = 27 + 1 (GF(27)),
52 = 2(25 + 1) (GF(25)), 76 = 2(37 + 1). **The first miss is 92**: 91 = 7·13 and 45 are not prime
powers, and 92 = 2·46 = 4·23 has no factorisation into known orders. The next miss is 116. Run
`python hadamard.py gaps 1000` for the full list up to any bound; the test suite checks that the list
up to 120 is exactly [92, 116] and builds and verifies every other order.

## 5. The methods behind the hard orders

**Williamson (1944).** If *A*, *B*, *C*, *D* are symmetric, pairwise commuting ±1 matrices of order *n*
with *AAᵀ* + *BBᵀ* + *CCᵀ* + *DDᵀ* = 4*n I*, then

    [  A   B   C   D ]
    [ -B   A  -D   C ]
    [ -C   D   A  -B ]
    [ -D  -C   B   A ]

is Hadamard of order 4*n*. Taking *A*, ..., *D* circulant turns the problem into finding four ±1
sequences of length *n* whose periodic autocorrelations sum to zero at every nonzero shift. That search
is what Baumert, Golomb and Hall ran at JPL to find **order 92 in 1962**, and Baumert found **116**
(and 232) the same way in 1966.

**Goethals–Seidel (1970).** With circulant *A*, ..., *D* (no symmetry needed) and the back-diagonal
permutation *R*,

    [   A     BR     CR     DR  ]
    [ -BR     A     DᵀR   -CᵀR  ]
    [ -CR   -DᵀR     A     BᵀR  ]
    [ -DR    CᵀR   -BᵀR     A   ]

is Hadamard of order 4*n* under the same autocorrelation condition. (Check: for circulants *X*, *XR* = *RXᵀ*
and *XR* is symmetric; every block-row product then cancels in pairs.) Dropping symmetry enlarges the
search space a great deal. Fed by supplementary difference sets on a cyclic group, this array gave
**1004** and **2524** (Đoković, Golubitsky and Kotsireas 2013).
`hadamard.py gs` builds and verifies this array from four sequences.

**Sequences with zero aperiodic autocorrelation.** Kharaghani and Tayfeh-Rezaie settled **428** in 2004
(published 2005) from four Turyn-type sequences of lengths 36, 36, 36, 35, giving base sequences of
lengths 71, 71, 36, 36 and T-sequences of length 107.

**Asymptotic results.** These say that every odd *q* times a large enough power of 2 is an order.
Seberry (1976) showed a Hadamard matrix of order 2ᵗ*q* exists once *t* is at least about
2 log₂ *q*; Craigen (1995) brought the needed exponent down to roughly (1/4) log₂ *q*. The conjecture needs *t* = 2
for every odd *q*, so these results are far from it.

## 6. The smallest open order, 1962 to 2026

| When | Event | Smallest multiple of 4 with no known matrix afterwards |
|---|---|---|
| 1962 | Order 92 found (Baumert, Golomb, Hall) | moved past 92 |
| 2004–05 | Order 428 found (Kharaghani, Tayfeh-Rezaie); it had been the smallest open order for about two decades | 668 |
| 2007 | Order 764 found (Đoković) | 668 |
| 2013 | Orders 1004 and 2524 found (Đoković, Golubitsky, Kotsireas) | 668 |
| 2025 | Eliahou builds a ±1 matrix of order 668 with *HHᵀ* ≡ 668 *I* (mod 64), the closest approximation then known; open cases up to 1000 were 668, 716, 892 | 668 |
| 2026-08-12 | Alpöge, Voinov and Reynolds-Haertle, working with Claude, publish matrices for all twelve remaining open orders below 2000: 668, 716, 892, 1132, 1244, 1388, 1436, 1676, 1772, 1916, 1948, 1964 | above 2000 |

What is known about the August 2026 result, and how solid it is:

- The matrices were posted as data (an X post of ±1 symbols with a decoder reply, a pastebin of the
  order-668 matrix and a Google Drive folder with the rest; links in OEIS A007299).
- A Hadamard claim is checkable in seconds, and it was checked independently: Ion Nechita verified the
  668 matrix with a short script; VibeMathed reports verifying all twelve; a Wolfram notebook cited on
  MathWorld confirms *H Hᵀ* = 668 *I*. OEIS A007299 and MathWorld now record all orders below 2000 as
  constructed. `python hadamard.py check FILE` lets anyone repeat the check.
- The **method has not been published**. Nechita notes that the matrix itself "tells us almost
  nothing" about how it was found. Epoch AI's FrontierMath page marks the order-668 problem as
  solved, with the AI attribution explicitly provisional. An arXiv paper (2609.25543) notes that the
  result removes the motivation for the 64-modular programme at 668.

Earlier claims about the gaps should be read with this date in mind: anything written before
2026-08-12 (including TheoremDB's entry and Eliahou's paper) correctly lists 668 as open *at the time*.

## 7. Proved versus open, as of 2026-10-14

**Proved (with proofs above or in the cited papers):**

1. Orders of Hadamard matrices are 1, 2 or multiples of 4 (section 3).
2. Orders are closed under multiplication; *q* + 1 (*q* ≡ 3 mod 4) and 2(*q* + 1) (*q* ≡ 1 mod 4) are
   orders for every prime power *q* (section 4).
3. Every multiple of 4 below 2000 is an order: classical constructions plus the 2026 matrices, which
   are explicit, public and independently machine-checked. This is a finite verification, not a
   theorem with a published method.
4. For every odd *q*, 2ᵗ*q* is an order for all large enough *t* (Seberry; Craigen).
5. Number of inequivalent Hadamard matrices of orders 16, 20, 24, 28, 32: 5, 3, 60, 487, 13,710,027
   (OEIS A007299; order 32 by Kharaghani and Tayfeh-Rezaie 2013).

**Open:**

1. The Hadamard conjecture itself: existence for every multiple of 4.
2. The smallest multiple of 4 with no known Hadamard matrix. It is now above 2000. We could not find a
   maintained public table after 2026-08-12 that names it, so we do not name one. For scale: in 2011
   Đoković counted 822 odd *n* < 10000 for which no Hadamard matrix of order 4*n* was then known
   (arXiv:1008.2043); some have been settled since, including 1004 and 2524 (2013) and the twelve
   orders of 2026.
3. How the 2026 matrices were found, and whether the method extends to all orders.
4. Structured versions are open even at small orders; for example, the smallest order with no known
   cocyclic Hadamard matrix was 188 as of 2016 (arXiv:1603.01859).

## 8. What is original here

- A self-contained proof path from the definition to both Paley constructions and the
  Goethals–Seidel array, with the block-by-block cancellation spelled out.
- `hadamard.py`: a dependency-free checker (rows packed into integers, so a 2000 × 2000 matrix can be
  checked in plain Python without numpy), Paley I/II over any prime-power field (GF(*p*ᵏ) built from a primitive polynomial found
  by search), a Goethals–Seidel builder that first checks the autocorrelation condition, and the
  "elementary gaps" computation of section 4.3.
- `test_hadamard.py`: 27 tests, including the Jacobsthal identities over GF(9), GF(25), GF(27), every
  elementary order up to 120 rebuilt and verified, and the gap list [92, 116].

Nothing here claims new Hadamard orders.

## 9. Using the companion code

    python hadamard.py check H668.txt     # rows of + and -, or 1/-1 tokens; one long line also works
    python hadamard.py gs quad.txt        # four +/- sequences of length n -> order 4n
    python hadamard.py gaps 1000          # orders the elementary constructions miss
    python hadamard.py build 100          # Paley II over GF(49)
    python -m unittest test_hadamard      # or: python test_hadamard.py

Python 3.6+, no dependencies.

## Sources

- J. Hadamard, Résolution d'une question relative aux déterminants, Bull. Sci. Math. 17 (1893) 240–246.
- J. J. Sylvester, Thoughts on inverse orthogonal matrices, Phil. Mag. 34 (1867) 461–475.
- R. E. A. C. Paley, On orthogonal matrices, J. Math. Phys. 12 (1933) 311–320.
- J. Williamson, Hadamard's determinant theorem and the sum of four squares, Duke Math. J. 11 (1944) 65–81.
- L. Baumert, S. W. Golomb, M. Hall Jr., Discovery of an Hadamard matrix of order 92, Bull. AMS 68 (1962) 237–238. https://projecteuclid.org/journals/bulletin-of-the-american-mathematical-society/volume-68/issue-3/Discovery-of-an-Hadamard-matrix-of-order-92/bams/1183524570.full
- L. D. Baumert, Hadamard matrices of orders 116 and 232, Bull. AMS 72 (1966) 237.
- J. M. Goethals, J. J. Seidel, A skew Hadamard matrix of order 36, J. Austral. Math. Soc. 11 (1970) 343–344.
- J. Seberry Wallis, On the existence of Hadamard matrices, J. Combin. Theory Ser. A 21 (1976) 188–195.
- R. Craigen, Signed groups, sequences, and the asymptotic existence of Hadamard matrices, J. Combin. Theory Ser. A 71 (1995) 241–254.
- H. Kharaghani, B. Tayfeh-Rezaie, A Hadamard matrix of order 428, J. Combin. Des. 13 (2005) 435–440. https://onlinelibrary.wiley.com/doi/10.1002/jcd.20043
- D. Ž. Đoković, Hadamard matrices of order 764 exist, Combinatorica 28 (2008); arXiv:math/0703312.
- D. Ž. Đoković, Small orders of Hadamard matrices and base sequences, arXiv:1008.2043 (2010).
- D. Ž. Đoković, O. Golubitsky, I. S. Kotsireas, Some new orders of Hadamard and skew-Hadamard matrices, J. Combin. Des. 22 (2014) 270–277; arXiv:1301.3671.
- H. Kharaghani, B. Tayfeh-Rezaie, Hadamard matrices of order 32, J. Combin. Des. 21 (2013) 212–221.
- S. Eliahou, 64-modular Hadamard matrix of order 668, Australas. J. Combin. 93(2) (2025) 422–427. https://ajc.maths.uq.edu.au/pdf/93/ajc_v93_p422.pdf
- OEIS A007299 (history and links to the 2026 matrices): https://oeis.org/A007299
- Epoch AI, FrontierMath open problem "Hadamard matrix of order 668": https://epoch.ai/frontiermath/open-problems/hadamard
- MathWorld, Hadamard Matrix: https://mathworld.wolfram.com/HadamardMatrix.html
- I. Nechita, New Hadamard matrices (2026-08-13): https://ion.nechita.net/posts/new-hadamard-matrices/
- J. D. Cook, Constructing Hadamard matrices (2026-08-13): https://www.johndcook.com/blog/2026/08/13/constructing-hadamard-matrices/
- VibeMathed, Hadamard matrix of order 668: https://vibemathed.com/problem/hadamard-matrix-of-order-668
- arXiv:2609.25543 (32-modular Hadamard matrices; comments on the 2026 result) and arXiv:2609.18380 (cites the announcement).
- Gröbner bases and cocyclic Hadamard matrices, arXiv:1603.01859 (2016).

## Licence

Text: CC BY 4.0. Code: MIT (see LICENSE). Cite as: MV2 (Munim, Inc.), "The Hadamard conjecture in
October 2026", 2026.

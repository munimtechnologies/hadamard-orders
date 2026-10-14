import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import hadamard as h  # noqa: E402


def naive_is_hadamard(rows):
    n = len(rows)
    if any(len(r) != n for r in rows):
        return False
    if any(x not in (1, -1) for r in rows for x in r):
        return False
    for i in range(n):
        for j in range(n):
            dot = sum(rows[i][k] * rows[j][k] for k in range(n))
            if dot != (n if i == j else 0):
                return False
    return True


class TestCheck(unittest.TestCase):
    def test_small_known(self):
        self.assertTrue(h.is_hadamard([[1]]))
        self.assertTrue(h.is_hadamard([[1, 1], [1, -1]]))

    def test_rejects_bad_entry(self):
        ok, msg = h.check([[1, 0], [1, -1]])
        self.assertFalse(ok)
        self.assertIn("entry", msg)

    def test_rejects_non_square(self):
        self.assertFalse(h.is_hadamard([[1, 1], [1, -1], [1, 1]]))

    def test_rejects_non_orthogonal(self):
        m = h.sylvester(3)
        m[5][2] = -m[5][2]
        ok, msg = h.check(m)
        self.assertFalse(ok)
        self.assertIn("inner product", msg)

    def test_fast_check_matches_naive(self):
        for m in (h.sylvester(3), h.paley1(11), h.paley2(5)):
            self.assertEqual(h.is_hadamard(m), naive_is_hadamard(m))
            bad = [row[:] for row in m]
            bad[1][1] = -bad[1][1]
            self.assertEqual(h.is_hadamard(bad), naive_is_hadamard(bad))
            self.assertFalse(h.is_hadamard(bad))


class TestConstructions(unittest.TestCase):
    def test_sylvester(self):
        for k in range(6):
            m = h.sylvester(k)
            self.assertEqual(len(m), 2 ** k)
            self.assertTrue(h.is_hadamard(m))

    def test_prime_power(self):
        self.assertEqual(h.prime_power(2), (2, 1))
        self.assertEqual(h.prime_power(27), (3, 3))
        self.assertEqual(h.prime_power(49), (7, 2))
        self.assertEqual(h.prime_power(43), (43, 1))
        self.assertIsNone(h.prime_power(1))
        self.assertIsNone(h.prime_power(45))
        self.assertIsNone(h.prime_power(91))

    def test_jacobsthal_properties(self):
        for q in (5, 7, 9, 11, 13, 25, 27):
            Q = h.jacobsthal(q)
            for a in range(q):
                self.assertEqual(sum(Q[a]), 0)
                for b in range(q):
                    dot = sum(Q[a][x] * Q[b][x] for x in range(q))
                    self.assertEqual(dot, q - 1 if a == b else -1)
                    sign = 1 if q % 4 == 1 else -1
                    self.assertEqual(Q[b][a], sign * Q[a][b])

    def test_paley1(self):
        for q in (3, 7, 11, 19, 23, 27, 31, 43):
            m = h.paley1(q)
            self.assertEqual(len(m), q + 1)
            self.assertTrue(h.is_hadamard(m), q)

    def test_paley2(self):
        for q in (5, 9, 13, 17, 25):
            m = h.paley2(q)
            self.assertEqual(len(m), 2 * (q + 1))
            self.assertTrue(h.is_hadamard(m), q)

    def test_paley_rejects_wrong_q(self):
        with self.assertRaises(ValueError):
            h.paley1(5)
        with self.assertRaises(ValueError):
            h.paley2(7)
        with self.assertRaises(ValueError):
            h.paley1(15)

    def test_kron(self):
        m = h.kron(h.paley1(11), h.H2)
        self.assertEqual(len(m), 24)
        self.assertTrue(h.is_hadamard(m))


class TestGoethalsSeidel(unittest.TestCase):
    def test_order_4(self):
        m = h.goethals_seidel([1], [1], [1], [1])
        self.assertTrue(h.is_hadamard(m))

    def test_williamson_order_12(self):
        m = h.goethals_seidel([1, 1, 1], [1, -1, -1], [1, -1, -1], [1, -1, -1])
        self.assertEqual(len(m), 12)
        self.assertTrue(h.is_hadamard(m))

    def test_order_20(self):
        a = [1, -1, -1, -1, -1]
        b = [1, 1, -1, -1, 1]
        c = [1, -1, 1, 1, -1]
        d = [1, -1, -1, -1, -1]
        m = h.goethals_seidel(a, b, c, d)
        self.assertEqual(len(m), 20)
        self.assertTrue(h.is_hadamard(m))

    def test_non_symmetric_sequences(self):
        # cyclic shifts keep periodic autocorrelations, so the array must still work
        a = [-1, -1, -1, -1, 1]
        b = [1, -1, -1, 1, 1]
        c = [-1, 1, 1, -1, 1]
        d = [-1, -1, 1, -1, -1]
        self.assertTrue(h.is_hadamard(h.goethals_seidel(a, b, c, d)))

    def test_rejects_bad_quadruple(self):
        with self.assertRaises(ValueError):
            h.goethals_seidel([1, 1, 1], [1, 1, 1], [1, -1, -1], [1, -1, -1])
        with self.assertRaises(ValueError):
            h.goethals_seidel([1, 1], [1], [1], [1])


class TestElementaryOrders(unittest.TestCase):
    def test_first_gap_is_92(self):
        self.assertEqual(h.elementary_gaps(100), [92])

    def test_gaps_up_to_120(self):
        # 92 and 116 need a method beyond Sylvester, Paley and Kronecker products
        self.assertEqual(h.elementary_gaps(120), [92, 116])

    def test_build_everything_reached_up_to_120(self):
        recipes = h.elementary_recipes(120)
        for n in range(4, 121, 4):
            if n in (92, 116):
                self.assertNotIn(n, recipes)
                continue
            m = h.build(n)
            self.assertEqual(len(m), n)
            self.assertTrue(h.is_hadamard(m), n)

    def test_build_rejects_gap(self):
        with self.assertRaises(ValueError):
            h.build(92)


class TestParse(unittest.TestCase):
    def test_compact(self):
        rows = h.parse_rows("# comment\n++\n+-\n")
        self.assertEqual(rows, [[1, 1], [1, -1]])

    def test_tokens(self):
        rows = h.parse_rows("1 1\n1, -1\n")
        self.assertEqual(rows, [[1, 1], [1, -1]])

    def test_single_line_square(self):
        rows = h.parse_rows("+++-")
        self.assertEqual(rows, [[1, 1], [1, -1]])

    def test_round_trip(self):
        m = h.paley1(7)
        self.assertEqual(h.parse_rows(h.format_rows(m)), m)

    def test_rejects_garbage(self):
        with self.assertRaises(ValueError):
            h.parse_rows("+x\n")


class TestCli(unittest.TestCase):
    def test_gaps_and_build(self):
        self.assertEqual(h.main(["gaps", "100"]), 0)
        self.assertEqual(h.main(["build", "92"]), 1)
        self.assertEqual(h.main(["build", "28"]), 0)


if __name__ == "__main__":
    unittest.main()

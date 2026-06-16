class QuantumCombinatorics:
    def __init__(self):
        # Cache to store previously computed q-binomials: keys are (n, k)
        self.binomial_cache = {}
        self.multinomial_cache = {}

    def _add_polynomials(self, p1, p2):
        """Element-wise addition of two polynomials."""
        length = max(len(p1), len(p2))
        result = [0] * length
        for i in range(len(p1)):
            result[i] += p1[i]
        for i in range(len(p2)):
            result[i] += p2[i]
        return result

    def _multiply_polynomials(self, p1, p2):
        """Multiplies two polynomials using discrete convolution."""
        if not p1 or not p2:
            return []
        result = [0] * (len(p1) + len(p2) - 1)
        for i, c1 in enumerate(p1):
            for j, c2 in enumerate(p2):
                result[i + j] += c1 * c2
        return result

    def q_binomial(self, n, k):
        """
        Computes the q-binomial coefficient (n choose k)_q.
        Returns a list of coefficients.
        """
        # Base cases
        if k < 0 or k > n:
            return [0]
        if k == 0 or k == n:
            return [1]  # Equals 1 (q^0)
        
        # Check cache
        if (n, k) in self.binomial_cache:
            return self.binomial_cache[(n, k)]

        # Recurrence: binom(n, k) = binom(n-1, k-1) + q^k * binom(n-1, k)
        left_term = self.q_binomial(n - 1, k - 1)
        right_term = self.q_binomial(n - 1, k)
        
        # Multiplying by q^k means shifting the polynomial array by k zeros
        shifted_right = [0] * k + right_term

        # Add them together and cache the result
        result = self._add_polynomials(left_term, shifted_right)
        self.binomial_cache[(n, k)] = result
        
        return result

    def q_multinomial(self, n, k_list):
        """
        Computes the q-multinomial coefficient.
        k_list is a list of the bottom parameters [k_1, k_2, ..., k_m].
        """
        if sum(k_list) != n:
            raise ValueError("The sum of elements in k_list must equal n.")
            
        result = [1]
        current_sum = k_list[0]
        
        # Factor into a product of q-binomials:
        # binom(k1+k2, k2) * binom(k1+k2+k3, k3) * ... * binom(n, k_m)
        for k_i in k_list[1:]:
            current_sum += k_i
            binom = self.q_binomial(current_sum, k_i)
            result = self._multiply_polynomials(result, binom)
            
        return result
    
    def get_multinomial(self, d):
        d_tuple = tuple(sorted(d))
        if d_tuple not in self.multinomial_cache:
            self.multinomial_cache[d_tuple] = self.q_multinomial(sum(d), list(d))
        return self.multinomial_cache[d_tuple]

    def format_polynomial(self, poly):
        """Helper to print the array as a readable mathematical string."""
        terms = []
        for power, coeff in enumerate(poly):
            if coeff == 0:
                continue
            if power == 0:
                terms.append(str(coeff))
            elif power == 1:
                terms.append(f"{coeff if coeff > 1 else ''}q")
            else:
                terms.append(f"{coeff if coeff > 1 else ''}q^{power}")
        return " + ".join(terms) if terms else "0"

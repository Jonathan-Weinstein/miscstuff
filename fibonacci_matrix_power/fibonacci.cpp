/*
g++ -std=c++14 -Wall -Wextra -Wshadow -Werror -Os fibonacci.cpp -o fib.out && ./fib.out

   fib(n) {
      fib(0): 0
      fib(1): 1
      fib(n): fib(n-2) + fib(n-1) | for n >= 2
   } -> 0,1,1,2,3,6,9,...

   fib_exp(e) = fib(n + 1) | for e >= -1

   constexpr mat2x2 M = { {1, 1}, {1, 0} }

   fib_exp(e) = (M^e).v0.x

   -----------------------------

   Modeling this via "complex" numbers might be interesting.
   Instead of the "imag" component representing the coef sqrt(-1),
   it'd represent the coef of sqrt(5).
*/
#include <stdint.h>
#include <assert.h>
#include <stdio.h>

#define bsr(v) (__builtin_clz(v) ^ 31)
#define bsf(v) __builtin_ctz(v)

struct diagonal_matrix_2v2u64 {
   uint64_t a; // v0.x, fib(i)
   uint64_t k; // v0.y == v1.x, fib(i - 1)
   uint64_t d; // v1.y
};

struct vec2u64 {
   uint64_t x, y;
};

vec2u64 mul(diagonal_matrix_2v2u64 m, vec2u64 v)
{
   return {
      m.a*v.x + m.k*v.y,
      m.k*v.x + m.d*v.y
   };
}

// returns mul(m, m)
diagonal_matrix_2v2u64 square(diagonal_matrix_2v2u64 m)
{
   uint64_t const t = m.k * m.k;
   return {
      m.a*m.a + t,
      m.k*(m.a + m.d),
      m.d*m.d + t
   };
}

/*
 * Should be approx the 5th root of UINT64_MAX:
 */
#define FIB_N_MAX_U64 93   // fib(n=93) last that doesn't overflow uint64_t,
                           // fib(n=94) first to overflow uint64_t

// Method A:
// O(popcount(n))
// as fast as you can get (I think) without doing `static const uint64_t FinalResult[94]`
// n must be < 94
uint64_t fib_bititerate_lut(unsigned n)
{
   // diagonal_matrix_2v2u64 m = { 1, 1, 0 }; // M^1 == M^2^0
   // for (int i = 0; i <= bsr(FIB_N_MAX_U64); ++i) { print(m); m = square(m); }
   static const diagonal_matrix_2v2u64 Table[7] = {
      { 0x1ull, 0x1ull, 0x0ull }, // M^2^0
      { 0x2ull, 0x1ull, 0x1ull }, // M^2^1
      { 0x5ull, 0x3ull, 0x2ull }, // M^2^2
      { 0x22ull, 0x15ull, 0xDull }, // M^2^3
      { 0x63Dull, 0x3DBull, 0x262ull }, // M^2^4
      { 0x35C7E2ull, 0x213D05ull, 0x148ADDull }, // M^2^5
      { 0xF9D297A859Dull, 0x9A661CA20BBull, 0x5F6C7B064E2ull } // M^2^6
   };

   assert(n <= FIB_N_MAX_U64);
   if (n <= 1) { return n; }
   unsigned e = n - 1; // e >= 1

   vec2u64 r = { 1, 0 };
   do {
      unsigned const shift = bsf(e);
      r = mul(Table[shift], r);
      e ^= 1u<<shift;
   } while (e);
   return r.x;
}

// Method B:
// O(log2(n))
// no look-up-table
// n can be >= 94, but the returned value is the actual answer mod 2**64
uint64_t fib_mod_2p64(unsigned n)
{
   if (n <= 1) { return n; }
   unsigned e = n - 1; // e >= 1

   diagonal_matrix_2v2u64 m = { 1, 1, 0 };
   vec2u64 r = { 1, 0 };
   do {
      if (e & 1u) { r = mul(m, r); }
      m = square(m);
      e >>= 1u;
   } while (e);
   return r.x;
}

/*
 * This is not related to fibonacci, but this is an alternative power method with
 * a different order that might be more suitable for floating-point:
 */
double pow_fp64_uint(double b, unsigned e)
{
   if (e == 0) { return 1; }
   double r0 = pow_fp64_uint(b, e/2u);
   double r1 = e & 1u ? r0*b : r0;
   return r0 * r1;
}

int main()
{
   uint64_t fib_ns2 = 0, fib_ns1 = 1, fib_n;
   unsigned n = 2;
   while ((fib_n = fib_ns2 + fib_ns1) >= fib_ns1) {
      assert(fib_n == fib_bititerate_lut(n));
      assert(fib_n == fib_mod_2p64(n));
      fib_ns2 = fib_ns1;
      fib_ns1 = fib_n;
      n++;
   }
   assert(n-1 == FIB_N_MAX_U64);
   puts("Hooray");

   // other random thing:
   const double b = 3.0;
   double r = 1.0;
   for (int e = 0; e <= 17; ++e) {
      assert(r == pow_fp64_uint(b, e));
      r *= b;
   }
   puts("Huzzah");
   return 0;
}

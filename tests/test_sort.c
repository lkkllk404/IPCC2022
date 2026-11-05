#define main pivot_program_main
#include "../pivot.c"
#undef main
#include <assert.h>

int main(void)
{
    unsigned state = 2022;
    const int sizes[] = {0, 1, 2, 15, 16, 17, 31, 63, 64, 500, 1024, 2048, 2049, 16384};
    const double extremes[] = {-DBL_MAX, -1e300, -DBL_MIN, -0.0, 0.0, DBL_MIN, 1e300, DBL_MAX};
    int checks = 0;
    for (int s = 0; s < 14; ++s)
        for (int mode = 0; mode < 8; ++mode)
        {
            int n = sizes[s], depth = 0;
            double *a = malloc(sizeof(double) * (n + 1));
            double *b = malloc(sizeof(double) * (n + 1));
            for (int i = 0; i < n; ++i)
            {
                state = state * 1664525u + 1013904223u;
                a[i] = mode == 0 ? (double)(int)state / 13.0
                     : mode == 1 ? (double)i : mode == 2 ? (double)(n-i)
                     : mode == 3 ? 0.0 : mode == 4 ? (double)(i % 5 - 2)
                     : mode == 5 ? (double)(i <= n/2 ? i : n-i)
                     : mode == 6 ? extremes[i % 8]
                     : ldexp((double)((int)(state >> 16)-32768), (i % 1800)-900);
                b[i] = a[i];
            }
            for (int size = n; size > 1; size >>= 1) depth += 2;
            SortProjection(a, n, depth);
            qsort(b, (size_t)n, sizeof(*b), CompareDouble);
            for (int i = 0; i < n; ++i) assert(a[i] == b[i]);
            free(b); free(a); ++checks;
        }
    puts("PASS 112 sorting cases: random/sorted/reversed/equal/repeated/organ-pipe/extremes/mixed-exponents, n=0..16384.");
    assert(checks == 112);
    return 0;
}

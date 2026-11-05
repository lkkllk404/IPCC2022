/* Test the optimized score, not just the final ranking. */
#define main pivot_program_main
#include "../pivot.c"
#undef main
#include <assert.h>

int main(void)
{
    unsigned state = 2022;
    const int sizes[] = {1, 2, 17, 64, 500};
    const double scales[] = {1e-12, 1.0, 1e12};
    int checked = 0;
    for (int s = 0; s < 5; ++s)
        for (int scale = 0; scale < 3; ++scale)
            for (int mode = 0; mode < 3; ++mode)
            {
                int n = sizes[s];
                double *xy = malloc(sizeof(double) * 2 * n);
                double *d = malloc(sizeof(double) * (size_t)n * n);
                double *scratch = malloc(sizeof(double) * 4 * n);
                for (int i = 0; i < 2 * n; ++i)
                {
                    state = state * 1664525u + 1013904223u;
                    xy[i] = mode == 0 ? ((double)(state >> 8) / 16777216.0 - 0.5) * scales[scale]
                          : mode == 1 ? (double)(i % 7) * scales[scale] : 0.0;
                }
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        double x = xy[2*i] - xy[2*j], y = xy[2*i+1] - xy[2*j+1];
                        d[(size_t)i*n+j] = sqrt(x*x+y*y);
                    }
                for (int k = 1; k <= 4 && k <= n; ++k)
                    for (int trial = 0; trial < 7; ++trial)
                    {
                        int pivots[4];
                        for (int p = 0; p < k; ++p) pivots[p] = (trial + p) % n;
                        long double expected = 0;
                        for (int i = 0; i < n; ++i)
                            for (int j = 0; j < n; ++j)
                            {
                                double maximum = 0;
                                for (int p = 0; p < k; ++p)
                                {
                                    double delta = fabs(d[(size_t)pivots[p]*n+i] - d[(size_t)pivots[p]*n+j]);
                                    if (delta > maximum) maximum = delta;
                                }
                                expected += maximum;
                            }
                        double actual = SumDistance(k, n, d, pivots, scratch);
                        if (fabsl(actual-expected) > 1e-11L * fabsl(expected) + 1e-280L)
                        {
                            fprintf(stderr, "score mismatch: n=%d k=%d scale=%g mode=%d actual=%.17g expected=%.21Lg\n",
                                    n, k, scales[scale], mode, actual, expected);
                            return 1;
                        }
                        ++checked;
                    }
                free(scratch); free(d); free(xy);
            }
    printf("PASS %d direct score comparisons (long-double reference, scales 1e-12..1e12).\n", checked);
    return 0;
}

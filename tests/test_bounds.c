#define main pivot_program_main
#include "../pivot.c"
#undef main
#include <assert.h>

int main(void)
{
    unsigned state = 2022;
    int checked = 0, skipped = 0;
    const int sizes[] = {1, 2, 17, 64, 500, 2048};
    const double scales[] = {1e-12, 1, 1e12};
    for (int s = 0; s < 6; ++s)
        for (int mode = 0; mode < 8; ++mode)
            for (int t = 0; t < 3; ++t)
            {
                int n = sizes[s];
                double *a = malloc(sizeof(double) * n), *copy = malloc(sizeof(double) * n);
                for (int i = 0; i < n; ++i)
                {
                    state = state * 1664525u + 1013904223u;
                    double x = (double)(state >> 8) / 16777216.0;
                    a[i] = mode == 0 ? x * scales[t] : mode == 1 ? (double)(i % 7) * scales[t]
                         : mode == 2 ? 0 : mode == 3 ? 1e12 + x * 1e-3
                         : mode == 4 ? ((i % 2) ? -x : x) * scales[t]
                         : mode == 5 ? x * DBL_MIN
                         : mode == 6 ? nextafter((double)(i % 17) * scales[t] / 16.0,
                                                 (i % 2) ? INFINITY : -INFINITY)
                         : (i == 0 ? scales[t] : x * DBL_EPSILON * scales[t]);
                }
                double lower, upper;
                ProjectionBounds(a, n, &lower, &upper);
                long double exact = 0;
                for (int i = 0; i < n; ++i)
                    for (int j = i + 1; j < n; ++j) exact += fabsl((long double)a[i] - a[j]);
                memcpy(copy, a, sizeof(double) * n);
                double algorithm = ProjectionSum(copy, n);
                assert(lower <= exact && exact <= upper);
                assert(lower <= algorithm && algorithm <= upper);
                free(copy); free(a); ++checked;
            }
    /* Nonfinite range and overflow must conservatively decline to prune. */
    double extremes[] = {-DBL_MAX, DBL_MAX}, lower, upper;
    ProjectionBounds(extremes, 2, &lower, &upper);
    assert(lower == 0 && isinf(upper));
    assert(PairMayRank(DBL_MAX, DBL_MAX, 64, DBL_MAX, 0));

    const int rankingSizes[] = {1, 8, 32};
    for (int mode = 0; mode < 3; ++mode)
        for (int scale = 0; scale < 3; ++scale)
            for (int rank = 0; rank < 3; ++rank)
            {
                int m = rankingSizes[rank];
                const int n = 64;
                double xy[128], d[4096], scratch[128], singles[64];
                for (int i = 0; i < 128; ++i)
                {
                    state = state * 1664525u + 1013904223u;
                    xy[i] = scales[scale] * (mode == 0 ? (double)(state >> 8) / 16777216.0
                                           : mode == 1 ? (double)(i % 5) : (double)(i % 2));
                }
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        double x = xy[2*i]-xy[2*j], y = xy[2*i+1]-xy[2*j+1];
                        d[i*n+j] = sqrt(x*x+y*y);
                    }
                    memcpy(scratch, d + i*n, sizeof(double) * n);
                    singles[i] = 2 * ProjectionSum(scratch, n);
                }
                double *expected = malloc(sizeof(double) * 2*m), *actual = malloc(sizeof(double) * 2*m);
                int *oldp = calloc((size_t)4*m, sizeof(int)), *newp = calloc((size_t)4*m, sizeof(int));
                for (int i = 0; i < m; ++i)
                {
                    expected[i] = actual[i] = 0;
                    expected[m+i] = actual[m+i] = DBL_MAX;
                }
                for (int i = 0; i < n; ++i)
                    for (int j = i + 1; j < n; ++j)
                    {
                        int pivots[] = {i, j};
                        double exact = SumDistance(2, n, d, pivots, scratch), score = 0;
                        int use = PairMayRank(singles[i], singles[j], n, actual[m-1], actual[2*m-1]) &&
                                  PairScoreIfUseful(n, d, pivots, scratch, actual[m-1], actual[2*m-1], &score);
                        if (use)
                        {
                            assert(score == exact);
                            InsertRanking(actual, newp, m, 2, score, pivots, 1);
                            InsertRanking(actual+m, newp+2*m, m, 2, score, pivots, 0);
                        }
                        else
                        {
                            assert(exact <= actual[m-1] && exact >= actual[2*m-1]);
                            ++skipped;
                        }
                        InsertRanking(expected, oldp, m, 2, exact, pivots, 1);
                        InsertRanking(expected+m, oldp+2*m, m, 2, exact, pivots, 0);
                        assert(memcmp(expected, actual, sizeof(double) * 2*m) == 0);
                        assert(memcmp(oldp, newp, sizeof(int) * 4*m) == 0);
                        ++checked;
                    }
                free(newp); free(oldp); free(actual); free(expected);
            }
    printf("PASS %d interval/ranking checks; %d exact evaluations pruned without changing any score/ID.\n", checked, skipped);
    return 0;
}

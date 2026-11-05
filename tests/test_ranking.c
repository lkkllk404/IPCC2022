#define main pivot_program_main
#include "../pivot.c"
#undef main
#include <assert.h>

int main(void)
{
    unsigned state = 2022;
    const int sizes[] = {1, 8, 1000};
    int checked = 0;
    for (int s = 0; s < 3; ++s)
        for (int k = 1; k <= 4; ++k)
            for (int descending = 0; descending <= 1; ++descending)
            {
                int m = sizes[s];
                double *legacy = malloc(sizeof(double) * (m + 1));
                double *actual = malloc(sizeof(double) * m);
                int *oldp = calloc((size_t)(m + 1) * k, sizeof(int));
                int *newp = calloc((size_t)m * k, sizeof(int));
                for (int i = 0; i < m; ++i) legacy[i] = actual[i] = descending ? 0 : DBL_MAX;
                for (int t = 0; t < 5000; ++t)
                {
                    state = state * 1664525u + 1013904223u;
                    double score = t < 1000 ? (double)(t % 32) : (double)((state >> 16) % 32);
                    int pivots[4];
                    for (int p = 0; p < k; ++p) pivots[p] = t * k + p;
                    legacy[m] = score;
                    memcpy(oldp + (size_t)m * k, pivots, (size_t)k * sizeof(int));
                    for (int a = m; a > 0; --a)
                        if (descending ? legacy[a] > legacy[a-1] : legacy[a] < legacy[a-1])
                        {
                            double tmp = legacy[a]; legacy[a] = legacy[a-1]; legacy[a-1] = tmp;
                            for (int p = 0; p < k; ++p)
                            {
                                int tmpi = oldp[a*k+p]; oldp[a*k+p] = oldp[(a-1)*k+p]; oldp[(a-1)*k+p] = tmpi;
                            }
                        }
                    InsertRanking(actual, newp, m, k, score, pivots, descending);
                    assert(memcmp(actual, legacy, sizeof(double) * m) == 0);
                    assert(memcmp(newp, oldp, sizeof(int) * (size_t)m * k) == 0);
                    ++checked;
                }
                free(newp); free(oldp); free(actual); free(legacy);
            }
    printf("PASS %d ranking updates, exact legacy scores and stable-tie IDs (M=1/8/1000, k=1..4).\n", checked);
    return 0;
}

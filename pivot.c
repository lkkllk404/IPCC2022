#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#include <float.h>

static int CompareDouble(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

// Typed comparisons avoid qsort's indirect calls in the common path.
// Recurse on the smaller partition, with a depth-limited qsort fallback.
static void ComparisonSort(double *a, int n, int depth)
{
    while (n > 16)
    {
        if (depth-- == 0)
        {
            qsort(a, (size_t)n, sizeof(*a), CompareDouble);
            return;
        }
        double x = a[0], y = a[n / 2], z = a[n - 1];
        double pivot = x < y ? (y < z ? y : (x < z ? z : x))
                             : (x < z ? x : (y < z ? z : y));
        int i = 0, j = n - 1;
        while (i <= j)
        {
            while (a[i] < pivot) ++i;
            while (a[j] > pivot) --j;
            if (i <= j)
            {
                double tmp = a[i]; a[i] = a[j]; a[j] = tmp;
                ++i; --j;
            }
        }
        int left = j + 1, right = n - i;
        if (left < right)
        {
            ComparisonSort(a, left, depth);
            a += i; n = right;
        }
        else
        {
            ComparisonSort(a + i, right, depth);
            n = left;
        }
    }
    for (int i = 1; i < n; ++i)
    {
        double value = a[i];
        int j = i;
        while (j > 0 && a[j - 1] > value) { a[j] = a[j - 1]; --j; }
        a[j] = value;
    }
}

// Bounded stack workspace; use comparison sort outside this size range.
static void SortProjection(double *a, int n, int depth)
{
#if DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 && UINT64_MAX == 18446744073709551615ULL
    if (sizeof(double) == sizeof(uint64_t) && n >= 64 && n <= 2048)
    {
        uint64_t keys[2048], temporary[2048];
        uint64_t varying = 0;
        for (int i = 0; i < n; ++i)
        {
            uint64_t bits;
            memcpy(&bits, a + i, sizeof(bits));
            keys[i] = bits ^ ((bits >> 63) ? UINT64_MAX : (UINT64_C(1) << 63));
            varying |= keys[i] ^ keys[0];
        }
        uint64_t *source = keys, *target = temporary;
        // Sort the high 32 bits first; exact comparisons refine collisions.
        for (int shift = 32; shift < 64; shift += 8)
        {
            if (((varying >> shift) & 255) == 0) continue;
            unsigned counts[256] = {0};
            for (int i = 0; i < n; ++i) ++counts[(source[i] >> shift) & 255];
            unsigned offset = 0;
            for (int b = 0; b < 256; ++b)
            {
                unsigned count = counts[b]; counts[b] = offset; offset += count;
            }
            for (int i = 0; i < n; ++i)
                target[counts[(source[i] >> shift) & 255]++] = source[i];
            uint64_t *swap = source; source = target; target = swap;
        }
        for (int i = 0; i < n; ++i)
        {
            uint64_t bits = source[i] ^ ((source[i] >> 63) ? (UINT64_C(1) << 63) : UINT64_MAX);
            memcpy(a + i, &bits, sizeof(bits));
        }
        for (int first = 0; first < n; )
        {
            int end = first + 1;
            while (end < n && (source[first] >> 32) == (source[end] >> 32)) ++end;
            if (end - first > 1) ComparisonSort(a + first, end - first, depth);
            first = end;
        }
        return;
    }
#endif
    ComparisonSort(a, n, depth);
}

// Sum over unordered absolute differences of a sorted scalar projection.
static double ProjectionSum(double *values, int n)
{
    int depth = 0;
    for (int size = n; size > 1; size >>= 1) depth += 2;
    SortProjection(values, n, depth);
    double sum = 0;
    for (int i = 1; i < n; ++i)
        sum += (double)i * (n - i) * (values[i] - values[i - 1]);
    return sum;
}

// O(n log n) for k=2; O(n^2 k) for the generic fallback.
double SumDistance(const int k, const int n, const double *distances,
                   const int *pivots, double *rebuiltCoord)
{
    // Pivot-major rows: consecutive j values are consecutive memory accesses.
    if (k == 2)
    {
        const double *a = distances + (size_t)pivots[0] * n;
        const double *b = distances + (size_t)pivots[1] * n;
        // max(|dx|,|dy|) = (|dx+dy| + |dx-dy|)/2.
        // The ordered-pair factor 2 cancels the identity's factor 1/2.
        double *plus = rebuiltCoord, *minus = rebuiltCoord + n;
        for (int i = 0; i < n; ++i)
        {
            plus[i] = a[i] + b[i];
            minus[i] = a[i] - b[i];
        }
        return ProjectionSum(plus, n) + ProjectionSum(minus, n);
    }

    int i;
    for (i = 0; i < n; ++i)
        for (int p = 0; p < k; ++p)
            rebuiltCoord[i * k + p] = distances[(size_t)pivots[p] * n + i];

    // Calculate the sum of Chebyshev distance with rebuilt coordinates between every points
    double chebyshevSum = 0;
    for (i = 0; i < n; i++)
    {
        int j;
        // Unordered pairs suffice: d(i,j) == d(j,i), d(i,i) == 0.
        for (j = i + 1; j < n; j++)
        {
            double chebyshev = 0;
            int ki;
            for (ki = 0; ki < k; ki++)
            {
                double dis = fabs(rebuiltCoord[i * k + ki] - rebuiltCoord[j * k + ki]);
                chebyshev = dis > chebyshev ? dis : chebyshev;
            }
            chebyshevSum += chebyshev;
        }
    }

    return 2.0 * chebyshevSum;
}

// Strict admission and upper-bound insertion preserve the legacy stable ties.
static void InsertRanking(double *scores, int *stored, const int M, const int k,
                          double score, const int *pivots, int descending)
{
    if (descending ? score <= scores[M - 1] : score >= scores[M - 1]) return;
    int lo = 0, hi = M;
    while (lo < hi)
    {
        int mid = lo + (hi - lo) / 2;
        if (descending ? score > scores[mid] : score < scores[mid]) hi = mid;
        else lo = mid + 1;
    }
    size_t count = (size_t)(M - lo - 1);
    memmove(scores + lo + 1, scores + lo, count * sizeof(*scores));
    memmove(stored + (size_t)(lo + 1) * k, stored + (size_t)lo * k,
            count * k * sizeof(*stored));
    scores[lo] = score;
    memcpy(stored + (size_t)lo * k, pivots, (size_t)k * sizeof(*stored));
}

// Recursive function Combination() : combine pivots and calculate the sum of distance while combining different pivots.
// ki  : current depth of the recursion
// k   : number of pivots
// n   : number of points
// dim : dimension of metric space
// M   : number of combinations to store
// coord  : coordinates of points
// pivots : indexes of pivots
// maxDistanceSum  : the largest M distance sum
// maxDisSumPivots : the top M pivots combinations
// minDistanceSum  : the smallest M distance sum
// minDisSumPivots : the bottom M pivots combinations
void Combination(int ki, const int k, const int n, const int M,
                 const double *distances, double *scratch, int *pivots,
                 double *maxDistanceSum, int *maxDisSumPivots, double *minDistanceSum, int *minDisSumPivots)
{
    if (ki == k - 1)
    {
        int i;
        for (i = pivots[ki - 1] + 1; i < n; i++)
        {
            pivots[ki] = i;

            // Calculate sum of distance while combining different pivots.
            double distanceSum = SumDistance(k, n, distances, pivots, scratch);

            InsertRanking(maxDistanceSum, maxDisSumPivots, M, k, distanceSum, pivots, 1);
            InsertRanking(minDistanceSum, minDisSumPivots, M, k, distanceSum, pivots, 0);
        }
        return;
    }

    // Recursively call Combination() to combine pivots
    int i;
    for (i = pivots[ki - 1] + 1; i < n; i++)
    {
        pivots[ki] = i;
        Combination(ki + 1, k, n, M, distances, scratch, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);

        /** Iteration Log : pivots computed, best pivots, max distance sum, min distance sum pivots, min distance sum
        *** You can delete the logging code. **/
        if (ki == k - 2)
        {
            int kj;
            for (kj = 0; kj < k; kj++)
            {
                printf("%d ", pivots[kj]);
            }
            putchar('\t');
            for (kj = 0; kj < k; kj++)
            {
                printf("%d ", maxDisSumPivots[kj]);
            }
            printf("%lf\t", maxDistanceSum[0]);
            for (kj = 0; kj < k; kj++)
            {
                printf("%d ", minDisSumPivots[kj]);
            }
            printf("%lf\n", minDistanceSum[0]);
        }
    }
}

int main(int argc, char *argv[])
{
    // filename : input file namespace
    char *filename = (char *)"uniformvector-2dim-5h.txt";
    if (argc == 2)
    {
        filename = argv[1];
    }
    else if (argc != 1)
    {
        printf("Usage: ./pivot <filename>\n");
        return -1;
    }
    // M : number of combinations to store
    const int M = 1000;
    // dim : dimension of metric space
    int dim;
    // n : number of points
    int n;
    // k : number of pivots
    int k;

    // Read parameter
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("%s file not found.\n", filename);
        return -1;
    }
    if (fscanf(file, "%d%d%d", &dim, &n, &k) != 3 || dim < 1 || n < 1 ||
        k < 1 || k > n || n > INT_MAX / dim || n > INT_MAX / k ||
        (size_t)n > SIZE_MAX / sizeof(double) / (size_t)n)
    {
        fprintf(stderr, "Invalid input dimensions\n");
        fclose(file);
        return 1;
    }
    printf("dim = %d, n = %d, k = %d\n", dim, n, k);

    // Read Data
    double *coord = (double *)malloc(sizeof(double) * dim * n);
    if (!coord) { fprintf(stderr, "Out of memory\n"); fclose(file); return 1; }
    int i;
    for (i = 0; i < n; i++)
    {
        int j;
        for (j = 0; j < dim; j++)
        {
            if (fscanf(file, "%lf", &coord[i * dim + j]) != 1 || !isfinite(coord[i * dim + j]))
            {
                fprintf(stderr, "Invalid coordinate\n");
                fclose(file); free(coord); return 1;
            }
        }
    }
    fclose(file);

    // Start timing
    struct timeval start;
    struct timeval Time_point1, Time_point2;
    gettimeofday(&start, NULL);

    gettimeofday(&Time_point1, NULL);
    // maxDistanceSum : the largest M distance sum
    double *maxDistanceSum = (double *)malloc(sizeof(double) * (M + 1));
    if (!maxDistanceSum) { fprintf(stderr, "Out of memory\n"); return 1; }
    for (i = 0; i < M; i++)
    {
        maxDistanceSum[i] = 0;
    }
    // maxDisSumPivots : the top M pivots combinations
    int *maxDisSumPivots = (int *)malloc(sizeof(int) * (size_t)k * (M + 1));
    if (!maxDisSumPivots) { fprintf(stderr, "Out of memory\n"); return 1; }
    for (i = 0; i < M; i++)
    {
        int ki;
        for (ki = 0; ki < k; ki++)
        {
            maxDisSumPivots[i * k + ki] = 0;
        }
    }
    // minDistanceSum : the smallest M distance sum
    double *minDistanceSum = (double *)malloc(sizeof(double) * (M + 1));
    if (!minDistanceSum) { fprintf(stderr, "Out of memory\n"); return 1; }
    for (i = 0; i < M; i++)
    {
        minDistanceSum[i] = __DBL_MAX__;
    }
    // minDisSumPivots : the bottom M pivots combinations
    int *minDisSumPivots = (int *)malloc(sizeof(int) * (size_t)k * (M + 1));
    if (!minDisSumPivots) { fprintf(stderr, "Out of memory\n"); return 1; }
    for (i = 0; i < M; i++)
    {
        int ki;
        for (ki = 0; ki < k; ki++)
        {
            minDisSumPivots[i * k + ki] = 0;
        }
    }
    gettimeofday(&Time_point2, NULL);
    printf("Using time : %f ms\n", (Time_point2.tv_sec - Time_point1.tv_sec) * 1000.0 + (Time_point2.tv_usec - Time_point1.tv_usec) / 1000.0);

    // temp : indexes of pivots with dummy array head
    int *temp = (int *)malloc(sizeof(int) * (k + 1));
    if (!temp) { fprintf(stderr, "Out of memory\n"); return 1; }
    temp[0] = -1;

    // Build reusable pivot-major distance rows once, not per combination.
    double *distances = malloc(sizeof(double) * (size_t)n * n);
    double *scratch = malloc(sizeof(double) * (size_t)n * k);
    if (!distances || !scratch) { fprintf(stderr, "Out of memory\n"); return 1; }
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q)
        {
            double squared = 0;
            for (int d = 0; d < dim; ++d)
            {
                double delta = coord[p * dim + d] - coord[q * dim + d];
                squared += delta * delta;
            }
            if (!isfinite(squared) || squared > __DBL_MAX__ / 4)
            { fprintf(stderr, "Coordinate range too large\n"); return 1; }
            distances[(size_t)p * n + q] = sqrt(squared);
        }

    Combination(0, k, n, M, distances, scratch, &temp[1], maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);

    // End timing
    struct timeval end;
    gettimeofday(&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0);

    // Store the result
    FILE *out = fopen("result.txt", "w");
    if (!out) { perror("result.txt"); return 1; }
    for (i = 0; i < M; i++)
    {
        int ki;
        for (ki = 0; ki < k - 1; ki++)
        {
            fprintf(out, "%d ", maxDisSumPivots[i * k + ki]);
        }
        fprintf(out, "%d\n", maxDisSumPivots[i * k + k - 1]);
    }
    for (i = 0; i < M; i++)
    {
        int ki;
        for (ki = 0; ki < k - 1; ki++)
        {
            fprintf(out, "%d ", minDisSumPivots[i * k + ki]);
        }
        fprintf(out, "%d\n", minDisSumPivots[i * k + k - 1]);
    }
    fclose(out);

    // Log
    int ki;
    printf("max : ");
    for (ki = 0; ki < k; ki++)
    {
        printf("%d ", maxDisSumPivots[ki]);
    }
    printf("%lf\n", maxDistanceSum[0]);
    printf("min : ");
    for (ki = 0; ki < k; ki++)
    {
        printf("%d ", minDisSumPivots[ki]);
    }
    printf("%lf\n", minDistanceSum[0]);
    // for(i=0; i<M; i++){
    // int ki;
    // for(ki=0; ki<k; ki++){
    // printf("%d\t", maxDisSumPivots[i*k + ki]);
    // }
    // printf("%lf\n", maxDistanceSum[i]);
    // }
    // for(i=0; i<M; i++){
    // int ki;
    // for(ki=0; ki<k; ki++){
    // printf("%d\t", minDisSumPivots[i*k + ki]);
    // }
    // printf("%lf\n", minDistanceSum[i]);
    // }

    free(distances);
    free(scratch);
    free(temp);
    free(coord);
    free(maxDistanceSum);
    free(minDistanceSum);
    free(maxDisSumPivots);
    free(minDisSumPivots);
    return 0;
}

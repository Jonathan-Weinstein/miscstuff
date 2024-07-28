/*
## The problem, which is supposedly over 10 years old:

https://math.stackexchange.com/questions/2209/expected-number-of-loops
"There are [N] ropes in a bag. In each step, two rope ends are picked at random,
tied together and put back into a bag. The process is repeated until there are no free ends.
What is the expected number of loops at the end of the process?"

Descriptions on some other sites:
https://www.quora.com/You-have-a-bag-with-x-number-of-pieces-of-strings-in-it-You-randomly-grab-two-ends-and-tie-them-together-until-there-are-no-more-loose-ends-When-you-finish-doing-this-what-is-the-expected-number-of-loops-strings
https://www.reddit.com/r/riddles/comments/2byu4g/you_have_a_bag_of_with_n_number_of_strings_at/

Solution F(N) for N >= 2 (F(1) = 1) strings (ropes):
   = (1/(2N - 1)) * (1 + F(N - 1)) + (1 - 1/(2N-1)) * F(N - 1)
   = (probability of forming a loop on first step) * (1 + F(N - 1)) +
     (1 - that probability) * F(N-1). When tieing two strings together that dont
     form a loop, we can pretend they result in just 1 longer string, since it
     still has exactly 2 open ends.

   = F(N - 1) + 1/(2N - 1)
   = 1 + 1/3 + 1/5 + 1/7 + 1/9 + ... + 1/N
   = Harmonic(2N) - Harmonic(N)/2

Approximate solution, using H(n) ~= ln(n) + Gamma + 1/(2n) - 1/(12n^2)
and Gamma ~= 0.57721566490 for N = 1,000,000 is ~= 7.889510,
which is quite close to this programs Monte Carlo simulation result of
7.850000 for 1000 iterations with the default random seed.

A rougher approximation for the expected number of loops is ~= ln(N)/2 + 1.

## This program's design and view of the problem:

Each of the N strings has 2 "nodes" connected to each end of the string,
so there are 2N nodes in total. Each node is represented by a unique integer
in [0, 2N). Each even numbered node is implicitly connected to the odd node
above it, and each odd node to the even node below it (NodeID i
is implicitly connected to NodeID i xor 1).

Simulation: while there exists a node that hasn't been "welded" to another node,
pick 2 random nodes and weld them. A node cannot be welded to itself.
This is done by reading from a sequence of unique, non-repeating random
numbers in [0, 2N) so the picking process time is O(N). The sequence
of random numbers is generated by shuffling the sequence 0,1,2,...,2N-1.

Each node will have exactly 2 connections, one of each distinct connection type:
    - An implicit/string connection.
    - An explicit/weld connection.

Once all welds have been made, the number of loops (AKA connected-components)
are computed in O(N) time and O(1) extra space from the resulting graph
represented only by `NodeID welds[NumStrings*2]`.

See SimpleTest() for an example.

Skip to the @end_util marker in this document for the interesting parts.

---

Some possible ways to compile this:

Debug:
    gcc -Wall -Wextra -Wshadow -std=c99 -D_DEBUG -g string_weld_loop_count_simulation.c -o prog_debug.exe
Release:
    gcc -Wall -Wextra -Wshadow -std=c99 -DNDEBUG -O2 -s string_weld_loop_count_simulation.c -o prog_rel.exe

Linux release, notice -lm:
   gcc -Wall -Wextra -Wshadow -std=c99 -DNDEBUG -O2 -s string_weld_loop_count_simulation.c -lm -o prog_rel.out

This should compile as c++ too, e.g:
   g++ -x c++ -std=c++14 -Wall -Wextra -Wshadow -D_DEBUG -g string_weld_loop_count_simulation.c -o prog_debug.exe
*/


// Example usage/output:
// Processor: Intel(R) Core(TM) i5-7200U CPU @ 2.50GHz   2.70 GHz
// prog_rel.exe iters=1000 strings=1000000

/*
SimpleTest: result=5
Succeeded to change process priority class from NORMAL -> ABOVE_NORMAL.
Average milliseconds each iter (denom=997, first 3 iters discarded):
RandGen: 28.188206, CountConnectedComponents: 56.485805, Everything: 104.909090

Results for strings=1000000 iters=1000 seed=[[TODO]]:
  1:     0
  2:     6 ******
  3:    15 ***************
  4:    52 ****************************************************
  5:   102 ******************************************************************************************************
  6:   148 ****************************************************************************************************************************************************
  7:   152 ********************************************************************************************************************************************************
  8:   156 ************************************************************************************************************************************************************
  9:   124 ****************************************************************************************************************************
 10:    97 *************************************************************************************************
 11:    66 ******************************************************************
 12:    36 ************************************
 13:    26 **************************
 14:    11 ***********
 15:     5 *****
 16:     4 ****
Average answer = 7.850000
NDEBUG defined
*/


#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#if !defined(_POSIX_C_SOURCE) || _POSIX_C_SOURCE < 199309L
#define _POSIX_C_SOURCE 199309L
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#define ASSERT assert


// OS util: -------------------------------------------------------------------
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>
int64_t TimerTicksPerSecond(void) { LARGE_INTEGER li; QueryPerformanceFrequency(&li); return li.QuadPart; }
int64_t TimerGetTicks(void) { LARGE_INTEGER li; QueryPerformanceCounter(&li); return li.QuadPart; }
#else
#include <time.h>
#define TimerTicksPerSecond() 1000000000ll // ((long long)1e9)
int64_t TimerGetTicks(void)
{
   struct timespec ts;
   clock_gettime(CLOCK_MONOTONIC, &ts);
   return TimerTicksPerSecond() * (uint64_t)ts.tv_sec + ts.tv_nsec;
}
#endif
// ----------------------------------------------------------------------------


#if defined __GNUC__
#define NO_RETURN __attribute__((noreturn))
#define bsr(v) (__builtin_clz(v) ^ 31)
#elif defined _MSC_VER
#define NO_RETURN __declspec(noreturn)
#include <intrin.h>
#pragma intrinsic(_BitScanReverse)
static __forceinline int bsr(uint32_t v)
{
    unsigned long i;
    _BitScanReverse(&i, v);
    return i;
}
#else
#error "unkown compiler"
#endif


// misc util ------------------------------------------------------------------
void *xmalloc(size_t nbytes)
{
    void *p = malloc(nbytes);
    if (!p) {
        perror("malloc returned NULL");
        exit(EXIT_FAILURE);
    }
    return p;
}
#define XMALLOC_TYPED(T, n) ((T *)xmalloc((n) * sizeof(T)))

static NO_RETURN void FailedVerify(const char *expr, int line)
{
    printf("VERIFY(%s) failed on line %d\n", expr, line);
    exit(EXIT_FAILURE);
}
#define VERIFY(x) ((x) ? (void)0 : FailedVerify(#x, __LINE__))
// ----------------------------------------------------------------------------


// ============================================================================
// PCG by Melissa E. O'Neill
// Slight modifications where made here from original: https://github.com/imneme/pcg-c-basic

struct pcg_state_setseq_64 {    // Internals are *Private*.
    uint64_t state;             // RNG state.  All values are possible.
    uint64_t inc;               // Controls which RNG sequence (stream) is
                                // selected. Must *always* be odd.
};
typedef struct pcg_state_setseq_64 pcg32_random_t;
// If you *must* statically initialize it, here's one.
#define PCG32_INITIALIZER   { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL }

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq);
uint32_t pcg32_random_r(pcg32_random_t* rng);
uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound);

// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)
void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}

// pcg32_random_r(rng)
//     Generate a uniformly distributed 32-bit random number
uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    uint32_t xorshifted = (uint32_t) (((oldstate >> 18u) ^ oldstate) >> 27u);
    int32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// pcg32_boundedrand_r(rng, bound):
//     Generate a uniformly distributed number, r, where 0 <= r < bound
uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound)
{
    // To avoid bias, we need to make the range of the RNG a multiple of
    // bound, which we do by dropping output less than a threshold.
    // A naive scheme to calculate the threshold would be to do
    //
    //     uint32_t threshold = 0x100000000ull % bound;
    //
    // but 64-bit div/mod is slower than 32-bit div/mod (especially on
    // 32-bit platforms).  In essence, we do
    //
    //     uint32_t threshold = (0x100000000ull-bound) % bound;
    //
    // because this version will calculate the same modulus, but the LHS
    // value is less than 2^32.

    uint32_t threshold = -(int32_t)bound % bound;

    // Uniformity guarantees that this loop will terminate.  In practice, it
    // should usually terminate quickly; on average (assuming all bounds are
    // equally likely), 82.25% of the time, we can expect it to require just
    // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
    // (i.e., 2147483649), which invalidates almost 50% of the range.  In
    // practice, bounds are typically small and only a tiny amount of the range
    // is eliminated.
    for (;;) {
        uint32_t r = pcg32_random_r(rng);
        if (r >= threshold)
            return r % bound;
    }
}
// end PCG
// ============================================================================


// @end_util:


// Fisher-Yates/Knuth shuffle:
void Shuffle(uint32_t *bag, int32_t n, pcg32_random_t *rng)
{
    while (n >= 2) {
        uint32_t j = pcg32_boundedrand_r(rng, n); // result in [0, n)
        uint32_t i = --n;
        /* swap(a[i], a[j]): */
        uint32_t x = bag[i];
        uint32_t y = bag[j];
        bag[i] = y;
        bag[j] = x;
    }
}


#define MAX_STRINGS_LG2 27
#define MAX_NODES_LG2 (MAX_STRINGS_LG2 + 1)
typedef uint32_t NodeId;
struct ThreadParam {
    uint32_t nNodes; // nStrings * 2

    /*
     * welds[x] == y means x is welded to y, and
     * if that is the case then welds[y] == x must be true.
     * Each value is a unique integer in [0, nNodes),
     * and welds[i] != i
     */
    NodeId *welds; // length == nNodes;

    uint32_t *randomBag; // length == nNodes;
};
typedef struct ThreadParam ThreadParam;

#define StringConnection(x) ((x) ^ 1) // starting string connections, implicit

#ifdef _DEBUG
#define WELD_CHECK
#endif

void Construct(ThreadParam *g, uint32_t nNodes)
{
    VERIFY(!(nNodes & 1)); // must be even
    VERIFY(nNodes >= 2 && nNodes <= (1 << MAX_NODES_LG2));
    /* welds (NodeId[]) and randomBag (uint32_t[]) in same alloc: */
    char *base = XMALLOC_TYPED(char, nNodes * (sizeof(NodeId) + sizeof(uint32_t)));
    g->nNodes = nNodes;
    g->welds = (NodeId *)base;
    g->randomBag = (uint32_t *) ((NodeId *)base + nNodes);
}

void Destroy(ThreadParam *g)
{
    free(g->welds);
}

void Weld(NodeId a, NodeId b, NodeId *welds, uint32_t nNodes)
{
    ASSERT(a < nNodes && b < nNodes);
    ASSERT(a != b); // cant weld node to itself
#ifdef WELD_CHECK
    VERIFY(welds[a] == (NodeId)-1); // should not already be welded
    VERIFY(welds[b] == (NodeId)-1); // should not already be welded
#endif
   (void)nNodes; // only used for assert

    welds[a] = b;
    welds[b] = a;
}

/*
    The "Destructive" part means that the welds[] is not preserved.
    This allows not needing an extra visited[], though this isn't a big win.

    The neat part is that usually for flood-fill or graph traversal
    (DFS, BFS, etc) algorithms, a data structure is needed to store
    nodes to visit later (this data structure might be the program stack
    for recursive methods). However, since we know here that all the connected
    components form a circle, going in the same direction (e.g "counter-clockwise")
    each time is enough to traverse every node in the loop. Thus, we don't need a
    data structure to store nodes to visit.

    This algorithm uses O(1) extra space. A non-destructive version
    is doeable with O(NumNodes/2) space for a visited[] bool array, or 1/8 that
    for a bitset.

    The runtime complexity of this is O(NumNodes).

    For smaller node counts, this appears faster than Shuffle(), but is slower
    for larger node counts.
*/
uint32_t CountConnectedComponentsDestructive(uint32_t *welds, uint32_t nNodes)
{
    uint32_t result = 0;
    uint32_t nStringsAllComponents = 0;
    for (uint32_t outerIter = 0; outerIter < nNodes; outerIter += 2) { // NOTE: step by 2
        if ((int32_t)welds[outerIter] >= 0) { // unlikely
            int32_t currentNodeId = outerIter;
            uint32_t nStringsThisComponent = 0;
            uint32_t weldConn;
            do {
                weldConn = welds[currentNodeId];
                ASSERT((int32_t)weldConn >= 0);
                welds[currentNodeId & -2] = (NodeId)-1; // mark even node of currentNodeId's string
                nStringsThisComponent++;
                /* Traverse 1 weld then 1 string in the same direction (e.g counter-clockwise): */
                currentNodeId = StringConnection(weldConn);
            } while ((weldConn & -2) != outerIter); // (x & -2) clears lowest bit in x
            ASSERT((int32_t)welds[currentNodeId & -2] < 0);
            result++;
            nStringsAllComponents += nStringsThisComponent;
        }
    }
    ASSERT(nStringsAllComponents * 2 == nNodes);
    for (uint32_t i = 0; i < nNodes; i += 2) { // NOTE: step of 2
        ASSERT((int32_t)welds[i] < 0);
        ASSERT((int32_t)welds[i + 1] >= 0);
    }
    return result;
}

/*
 * Array "a" is of length nBalls and is pre-allocated.
 *
 * In the implementation, a[i] represents the current free ball of the _composite_ string
 * the ball of a[i] is attached to. After connecting verts "x" and "y", a[x] and a[y] do not matter.
 */
uint32_t SimulateAndCountLoopsMethod2(uint32_t nVerts, const uint32_t *__restrict shuffledDeck, uint32_t * __restrict a)
{
    // Starting state
    for (uint32_t i = 0; i < nVerts; ++i) {
        a[i] = i ^ 1;
    }

    uint32_t nLoops = 0;

    for (uint32_t i2 = 0; i2 < nVerts; i2 += 2) {
        /*
            Connect v1 to v2:

            *---*     *---*     ->    *------*
            v0  v1    v2  v3          v0     v3
        */
        uint32_t const v1 = shuffledDeck[i2];
        uint32_t const v2 = shuffledDeck[i2 + 1];

        uint32_t const v0 = a[v1];
        uint32_t const v3 = a[v2];

        a[v0] = v3;
        a[v3] = v0;

        ASSERT((v0 == v2) == (v3 == v1));

        if (v0 == v2) {
            nLoops++;
        }
    }

    return nLoops;
}

static void SimpleTest(void)
{
    static const struct {
        uint8_t a, b;
    } WeldsToMake[] = {
        // 1-string loop:
        { 0, 1 },
        // 3-string loop:
        { 7, 2 },
        { 3, 4 },
        { 5, 6 },
        // 4-string loop:
        { 15, 8 },
        { 9, 10 },
        { 11, 12 },
        { 13, 14 },
        // 1-string loop:
        { 16, 17 },
        // 2-string loop:
        { 21, 18 },
        { 19, 20 },
    };
    enum {
       NumWeldsToMake = sizeof(WeldsToMake) / sizeof(WeldsToMake[0]),
       NumNodes = NumWeldsToMake * 2
    };
    NodeId welds[NumNodes];
#ifdef WELD_CHECK
    memset(welds, 0xff, NumNodes * sizeof(NodeId)); // set everything to -1
#endif
    for (int i = 0; i < NumWeldsToMake; ++i) {
        Weld(WeldsToMake[i].a, WeldsToMake[i].b, welds, NumNodes);
    }
    int const result = CountConnectedComponentsDestructive(welds, NumNodes);
    if (result != 5) {
        puts(" *** TEST FAILED ***");
        exit(EXIT_FAILURE);
    }
}


//https://www.johndcook.com/blog/2017/04/18/computing-harmonic-numbers/
//http://fredrik-j.blogspot.com/2009/02/how-not-to-compute-harmonic-numbers.html
#include <math.h> // leeds to link with -lm on Linux
#define EulerGammaConstant 0.5772156649015328606065120
static double Harmonic(int32_t n_int)
{
   ASSERT(n_int >= 1);
   if (n_int <= 16) {
      double sum = 0.0;
      do {
         sum += 1.0 / (n_int);
      } while (--n_int);
      return sum;
   } else {
      // H(n) ~= ln(n) + Gamma + 1/(2n) - 1/(12n^2) + ...
      double const n = n_int;
      double r = log(n);
      r += EulerGammaConstant + 1.0/(n + n) - 1.0/(12*n*n);
      return r;
   }
}

static void Simulate(uint32_t const numStrings,
                     uint32_t const numIters,
                     double const MillisecondsPerTickF64,
                     bool const bTest)
{
    VERIFY(numStrings && numStrings <= (1u << MAX_STRINGS_LG2));
    VERIFY(numIters && numIters <= (1u << 15));

    enum { NumIterTimingDiscard = 3 };

    uint32_t const numNodes = numStrings * 2;

    double randgenMsAccum = 0.0;
    double countConnCompsMsAccum = 0.0;
    double totalMsAccum = 0.0;
    double method2MsAccum = 0.0;

    ThreadParam g;
    Construct(&g, numNodes);
    uint32_t *const bag = g.randomBag;
    pcg32_random_t rng = PCG32_INITIALIZER;

    unsigned highestTrackedAnswer = 0;
    enum { HistoCap = 128 };
    uint16_t histo[HistoCap] = { 0 };

#define TIME_IF(accum, expr, cond) { \
    int64_t _ticksBegin = TimerGetTicks(); \
    (expr); \
    if (cond) (accum) += (double) (TimerGetTicks() - _ticksBegin) * MillisecondsPerTickF64; \
}

    for (uint32_t currentIter = 0; currentIter < numIters; ++currentIter) {
        int64_t const totalTicksBegin = TimerGetTicks();
    #ifdef WELD_CHECK
        memset(g.welds, 0xff, g.nNodes * sizeof(NodeId)); // set everything to -1
    #endif
        for (uint32_t i = 0; i < numNodes; ++i) {
            bag[i] = i;
        }
        TIME_IF(randgenMsAccum, Shuffle(bag, numNodes, &rng), currentIter >= NumIterTimingDiscard);
        /* do numNodes/2 welds, each closes 2 nodes: */
        for (uint32_t i = 0; i < numNodes; i += 2) { // NOTE: step = 2
            Weld(bag[i], bag[i + 1], g.welds, g.nNodes);
        }
        unsigned answer;
        TIME_IF(countConnCompsMsAccum,
                answer = CountConnectedComponentsDestructive(g.welds, g.nNodes),
                currentIter >= NumIterTimingDiscard);
        if (bTest && currentIter >= NumIterTimingDiscard) {
            int64_t const t0 = TimerGetTicks();
            unsigned const answerMethod2 = SimulateAndCountLoopsMethod2(numNodes, bag, g.welds); // not "welds here"
            int64_t const t1 = TimerGetTicks();
            method2MsAccum += (double)(t1 - t0) * MillisecondsPerTickF64;
            VERIFY(answerMethod2 == answer);
        }
        if (answer < HistoCap) {
            histo[answer] += 1;
            highestTrackedAnswer = highestTrackedAnswer >= answer ? highestTrackedAnswer : answer;
        }
        else {
            printf("Very large answer of %d, not tracked in histo!\n", answer);
        }
        if (currentIter >= NumIterTimingDiscard) {
            totalMsAccum += (double) (TimerGetTicks() - totalTicksBegin) * MillisecondsPerTickF64;
        }
    }

    if (numIters > NumIterTimingDiscard) {
        int denom = numIters - NumIterTimingDiscard;
        double const r = 1.0 / denom;
        printf("Average milliseconds each iter (denom=%d, first %d iters discarded):\n"
            "RandGen: %f, CountConnectedComponents: %f, Method2: %f, Everything: %f\n\n",
            denom, NumIterTimingDiscard,
            randgenMsAccum * r,
            countConnCompsMsAccum *r,
            method2MsAccum * r,
            totalMsAccum * r);
    }

    printf("Results for strings=%d iters=%d seed=[[TODO]]:\n", numStrings, numIters);
    VERIFY(histo[0] == 0);
    uint32_t sumOfCounts = 0;
    int64_t sumOfAllAnswers = 0;
    for (unsigned answer = 1; answer <= highestTrackedAnswer; ++answer) {
        int count = histo[answer];
        sumOfCounts += count;
        sumOfAllAnswers += answer * count;
        printf("%3d: %5d ", answer, count);
        if (count > 160) {
            puts(" ~ a lot ~");
        }
        else {
            while (count--) {
                putchar('*');
            }
            puts("");
        }
    }
    VERIFY(sumOfCounts == numIters);
    printf("Average answer = %f, approx expected answer = %f\n",
            (double)sumOfAllAnswers / (double)numIters,
            Harmonic(numStrings * 2) - Harmonic(numStrings)/2.0);

    Destroy(&g);
}

int main(int argc, char **argv)
{
    SimpleTest();

    // uint64_t seed = ...;
    uint32_t numIters = 30;
    uint32_t numStrings = 1u << 15;
    bool bDebugTest = false;

    for (int i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        int val;
        if (sscanf(arg, "iters=%d", &val) == 1) {
            if (val <= 0 || val > (1 << 15)) {
                puts("bad iter count");
                return EXIT_FAILURE;
            }
            numIters = val;
        }
        else if (sscanf(arg, "strings=%d", &val) == 1) {
            if (val <= 0 || val > (1 << MAX_STRINGS_LG2)) {
                puts("bad string count");
                return EXIT_FAILURE;
            }
            numStrings = val;
        }
        else if (strcmp(arg, "test") == 0) {
            bDebugTest = true;
        }
        else {
            printf("unknown argv[%d]=%s\n", i, arg);
            return EXIT_FAILURE;
        }
    }

#if !defined(_DEBUG) && defined(_WIN32)
    /* Increase priority to try and make timers more consistent due to less context switches(?): */
    if (numIters >= 8) {
        HANDLE const hThisProcess = GetCurrentProcess();
        DWORD const startPriorityClass = GetPriorityClass(hThisProcess);
        if (startPriorityClass == NORMAL_PRIORITY_CLASS) {
            BOOL ok = SetPriorityClass(hThisProcess, ABOVE_NORMAL_PRIORITY_CLASS);
            printf("%s to change process priority class from NORMAL -> ABOVE_NORMAL.\n", ok ? "Succeeded" : "FAILED");
        }
        else {
            printf("Start priority class == 0x%X.\n", (unsigned int) startPriorityClass);
        }
    }
#endif

    const double MillisecondsPerTickF64 = 1000.0 / (double) TimerTicksPerSecond();

    if (!bDebugTest) {
        Simulate(numStrings, numIters, MillisecondsPerTickF64, false);
    }
    else {
        puts("doing debug test");
        for (int32_t x = -17; x <= 17; ++x) {
            int32_t newIterCount = numIters + x;
            if (newIterCount >= 1) {
                Simulate(numStrings, newIterCount, MillisecondsPerTickF64, true);
            }
        }
    }

    if (sizeof(void *) != 8) {
        puts("not 64 bit");
    }
#ifdef _DEBUG
    puts("_DEBUG defined");
#endif
#ifdef NDEBUG
    puts("NDEBUG defined");
#endif
    return 0;
}

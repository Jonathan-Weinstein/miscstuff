/*

Setup: N strings, each having 2 "nodes" conected but the string.

Then while there exists a node that hasn't been helded to another node
pick 2 random nodes and weld them. A node cannot be welded to itself.

After this is done, how many loops are formed?
See SimpleTest() for an example.

A possible way to compile this:
    g++ -x c -std=c99 -Wall -Wextra -Wshadow -DNDEBUG -O2 -s loop_simulation.cpp -o prog.exe
*/

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#define ASSERT assert


// OS util: -------------------------------------------------------------------
int64_t TimerTicksPerSecond() { LARGE_INTEGER li; QueryPerformanceFrequency(&li); return li.QuadPart; }
int64_t TimerGetTicks() { LARGE_INTEGER li; QueryPerformanceCounter(&li); return li.QuadPart; }
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


// Fisher–Yates/Knuth shuffle:
void Shuffle(uint32_t *bag, int32_t n, pcg32_random_t *rng)
{
    /*
        for (int i = n-1; i > 0; i--) {
            int j = UniformRandomExclusive(i + 1);
            Swap(ref bag[i], ref bag[j]);
        }
    */
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
struct StringGraph {
    uint32_t nNodes; // nStrings * 2
    NodeId *explicitConnection; // made by welds
};
typedef struct StringGraph StringGraph;
#define ImplicitConnection(x) ((x) ^ 1) // starting string connections

void StringGraph_Contruct(StringGraph *sg, uint32_t nNodes)
{
    ASSERT(!(nNodes & 1)); // must be even
    ASSERT(nNodes >= 2 && nNodes <= (1 << MAX_NODES_LG2));
    NodeId *p = XMALLOC_TYPED(NodeId, nNodes);
    sg->nNodes = nNodes;
    sg->explicitConnection = p;
    for (uint32_t i = 0; i < nNodes; ++i) {
        p[i] = -1;
    }
}

void StringGraph_Destroy(StringGraph *sg)
{
    free(sg->explicitConnection);
}

void StringGraph_Weld(StringGraph *sg, NodeId a, NodeId b)
{
    ASSERT(a < sg->nNodes && b < sg->nNodes);
    ASSERT(a != b); // cant weld node to itself
    ASSERT(sg->explicitConnection[a] == (NodeId)-1); // should not already be welded
    ASSERT(sg->explicitConnection[b] == (NodeId)-1); // should not already be welded

    sg->explicitConnection[a] = b;
    sg->explicitConnection[b] = a;
}

uint32_t StringGraph_CountConnectedComponentsDestructive(uint32_t *excon, uint32_t nNodes)
{
    uint32_t result = 0;
    uint32_t nStringsAllComponents = 0;
    for (uint32_t outerIter = 0; outerIter < nNodes; outerIter += 2) { // NOTE: step by 2
        if ((int32_t)excon[outerIter] >= 0) {
            int32_t currentNodeId = outerIter;
            uint32_t nStringsThisComponent = 0;
            uint32_t weldConn;
            do {
                weldConn = excon[currentNodeId];
                ASSERT((int32_t)weldConn >= 0);
                excon[currentNodeId & -2] = (NodeId)-1; // mark even node of currentNodeId's string
                nStringsThisComponent++;
                /* Traverse 1 weld then 1 string in the same direction (e.g counter-clockwise): */
                currentNodeId = ImplicitConnection(weldConn);
            } while ((weldConn & -2) != outerIter); // (x & -2) clears lowest bit in x
            ASSERT((int32_t)excon[currentNodeId & -2] < 0);
            result++;
            nStringsAllComponents += nStringsThisComponent;
        }
    }
    ASSERT(nStringsAllComponents * 2 == nNodes);
    for (uint32_t i = 0; i < nNodes; i += 2) { // NOTE: step of 2
        ASSERT((int32_t)excon[i] < 0);
        ASSERT((int32_t)excon[i + 1] >= 0);
    }
    return result;
}


uint32_t StringGraph_CountConnectedComponents(const StringGraph *g)
{
    uint32_t const nNodes = g->nNodes;
    NodeId const *const excon = g->explicitConnection;
    uint32_t const toVisitCapacity = nNodes * 2 + 2; // hmm
    NodeId *const toVisit = XMALLOC_TYPED(NodeId, toVisitCapacity);
    uint8_t *const visited = (uint8_t *)calloc(nNodes, sizeof(uint8_t)); // destructive idea
    if (visited == NULL) {
        exit(EXIT_FAILURE);
    }
    uint32_t result = 0;
    uint32_t nNodesAllComponents = 0;
    for (uint32_t outerIter = 0; outerIter < nNodes; outerIter += 2) { // NOTE: can step by 2
        if (!visited[outerIter]) {
            /* flood-fill */
            uint32_t currentNodeId = outerIter;
            uint32_t nToVisit = 0;
            uint32_t nNodesThisComponent = 0;
            for (;;) {
                if (!visited[currentNodeId]) {
                    uint32_t const a = ImplicitConnection(currentNodeId);
                    uint32_t const b = excon[currentNodeId];
                    /* NOTE: a == b is possible here. */
                    assert(b != currentNodeId);
                    assert(a < nNodes && b < nNodes);
                    visited[currentNodeId] = true;
                    nNodesThisComponent++;
                    if (!visited[a]) {
                    #if 0
                        if (nToVisit == toVisitCapacity) {
                            puts("oops a");
                            exit(EXIT_FAILURE);
                        }
                    #endif
                        toVisit[nToVisit++] = a;
                    }
                    if (!visited[b] && a != b) {
                    #if 0
                        if (nToVisit == toVisitCapacity) {
                            puts("oops b");
                            exit(EXIT_FAILURE);
                        }
                    #endif
                        toVisit[nToVisit++] = b;
                    }
                }

                if (nToVisit) {
                    currentNodeId = toVisit[--nToVisit];
                }
                else {
                    break;
                }
            }
            result++;
            nNodesAllComponents += nNodesThisComponent;
        }
    }

    ASSERT(nNodes == nNodesAllComponents);
#ifdef _DEBUG
    for (uint32_t i = 0; i < nNodes; ++i) {
        ASSERT(visited[i]);
    }
#endif
    free(visited);
    free(toVisit);
    ASSERT(StringGraph_CountConnectedComponentsDestructive(g->explicitConnection, g->nNodes) == result);
    return result;
}

void SimpleTest()
{
    static const struct {
        uint8_t a, b;
    } welds[] = {
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
    enum { NumWelds = sizeof(welds) / sizeof(welds[0]) };

    StringGraph g;
    StringGraph_Contruct(&g, 22);
    for (int i = 0; i < NumWelds; ++i) {
        StringGraph_Weld(&g, welds[i].a, welds[i].b);
    }
    int result = StringGraph_CountConnectedComponents(&g);
    printf("%s: result=%d\n", __FUNCTION__, result);
    if (result != 5) {
        puts(" *** TEST FAILED ***");
        exit(EXIT_FAILURE);
    }
    StringGraph_Destroy(&g);
}


static void Simulate(uint32_t const numStrings, uint32_t const numIters, double const MillisecondsPerTickF64)
{
    ASSERT(numStrings && numStrings <= (1u << MAX_STRINGS_LG2));
    ASSERT(numIters && numIters <= (1u << 15));

    enum { NumIterTimingDiscard = 3 };

    uint32_t const numNodes = numStrings * 2;

    double randgenMsAccum = 0.0;
    double countConnCompsMsAccum = 0.0;
    double totalMsAccum = 0.0;

    StringGraph g;
    StringGraph_Contruct(&g, numNodes);
    uint32_t *const bag = XMALLOC_TYPED(uint32_t, numNodes);
    pcg32_random_t rng = PCG32_INITIALIZER;

    unsigned highestTrackedAnswer = 0;
    enum { HistoCap = 128 };
    uint16_t histo[HistoCap] = { };

#define TIME_IF(accum, expr, cond) { \
    int64_t _ticksBegin = TimerGetTicks(); \
    (expr); \
    if (cond) (accum) += (double) (TimerGetTicks() - _ticksBegin) * MillisecondsPerTickF64; \
}

    for (uint32_t currentIter = 0; currentIter < numIters; ++currentIter) {
        int64_t const totalTicksBegin = TimerGetTicks();
        memset(g.explicitConnection, 0xff, numNodes * sizeof(NodeId));
        for (uint32_t i = 0; i < numNodes; ++i) {
            bag[i] = i;
        }
        TIME_IF(randgenMsAccum, Shuffle(bag, numNodes, &rng), currentIter >= NumIterTimingDiscard);
        /* do numNodes/2 welds, each closes 2 nodes: */
        for (uint32_t i = 0; i < numNodes; i += 2) { // NOTE: step = 2
            StringGraph_Weld(&g, bag[i], bag[i + 1]);
        }
        unsigned answer;
        //TIME_IF(countConnCompsMsAccum, answer = StringGraph_CountConnectedComponents(&g), currentIter >= NumIterTimingDiscard);
        TIME_IF(countConnCompsMsAccum,
                answer = StringGraph_CountConnectedComponentsDestructive(g.explicitConnection, g.nNodes),
                currentIter >= NumIterTimingDiscard);
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
            "RandGen: %f, CountConnectedComponents: %f, Everything: %f\n\n",
            denom, NumIterTimingDiscard,
            randgenMsAccum * r,
            countConnCompsMsAccum *r,
            totalMsAccum * r);
    }

    printf("Results for strings=%d iters=%d seed=???:\n", numStrings, numIters);
    ASSERT(histo[0] == 0);
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
    if (sumOfCounts != numIters) {
        puts("\noops!\n");
        exit(EXIT_FAILURE);
    }
    printf("Average answer = %f\n", (double)sumOfAllAnswers / (double)numIters);

    free(bag);
    StringGraph_Destroy(&g);
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

#if !defined(_DEBUG)
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
        Simulate(numStrings, numIters, MillisecondsPerTickF64);
    }
    else {
        puts("doing debug test");
        for (int32_t x = -17; x <= 17; ++x) {
            int32_t newIterCount = numIters + x;
            if (newIterCount >= 1) {
                Simulate(numStrings, newIterCount, MillisecondsPerTickF64);
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
// ^^^ revert that back.
// do the amllocs once.
// flood fill loop idea, no need to visit thing.
// WELD check


/*

Results for strings=1000000, iters=1000, seed=???:
  0:     0
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

*/

/*
Succeeded to change process priority class from NORMAL -> ABOVE_NORMAL.
Average milliseconds each iter (denom=128, first 3 iters discarded):
RandGen: 0.787859, CountConnectedComponents: 1.083199, Everything: 2.260526

Results for strings=64123 iters=131 seed=???:
  1:     0
  2:     3 ***
  3:     9 *********
  4:     9 *********
  5:    26 **************************
  6:    16 ****************
  7:    28 ****************************
  8:    22 **********************
  9:     8 ********
 10:     7 *******
 11:     3 ***
*/
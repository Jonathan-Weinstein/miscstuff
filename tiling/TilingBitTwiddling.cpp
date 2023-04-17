// utility begin ---------------------------------------------------------------
#include <stdint.h>
#include <assert.h>
#define ASSERT assert
typedef unsigned uint;

#ifdef _MSC_VER
extern "C" unsigned char _BitScanForward(unsigned long* Index, unsigned long Mask);
extern "C" unsigned char _BitScanReverse(unsigned long* Index, unsigned long Mask);
#pragma intrinsic(_BitScanForward, _BitScanReverse)
__forceinline unsigned long bsf(unsigned long v) { unsigned long i; _BitScanForward(&i, v); return i; }
__forceinline unsigned long bsr(unsigned long v) { unsigned long i; _BitScanReverse(&i, v); return i; }
#else
#define bsf(v) __builtin_ctz(v)
#define bsr(v) (__builtin_clz(v) ^ 31)
#endif

template<class T> T min(T a, T b) { return b < a ? b : a; }

inline uint32_t mask(uint n)
{
	ASSERT(32u >= n);
	return n ? ~uint32_t(0) >> (32 - n) : 0;
}

inline uint32_t ubfe(uint32_t base, uint offset, uint width)
{
	return (base >> offset) & mask(width);
}

// Overwrite bits selected by src0 into src2 with values taken from src1:
// dst = (src0 & src1) | (~src0 & src2)
//
// selector mask : 10111
//        source : 1x011
//          base : 0y101
//                 -----
//        result : 1y011
inline uint32_t v_bfi_b32(uint32_t selector, uint32_t source, uint32_t base)
{
	return (source & selector) | (base & ~selector);
}

struct uint2 {
	uint32_t x, y;
	uint32_t operator[](unsigned i) const { return i ? y : x; };
};
inline uint2 operator| (uint2 a, uint2 b) { return { a.x |  b.x, a.y |  b.y }; }
inline uint2 operator+ (uint2 a, uint2 b) { return { a.x +  b.x, a.y +  b.y }; }
inline uint2 operator- (uint2 a, uint2 b) { return { a.x -  b.x, a.y -  b.y }; }
inline uint2 operator* (uint2 a, uint2 b) { return { a.x *  b.x, a.y *  b.y }; }
inline uint2 operator<<(uint2 a, uint2 b) { return { a.x << b.x, a.y << b.y }; }
inline uint2 operator<<(uint2 a, uint sc) { return { a.x <<  sc, a.y <<  sc }; }
// utility end -----------------------------------------------------------------

#include <stdio.h>

/*
 * Simple 1-level of tiling when the blocks (AKA tiles) perfectly fit the rectangle of dimensions
 * { numItemsX, numItemsY }. An "item" is just a name given to the base-unit, it could be a thread, or something else.
 * 
 * itemId is a row-major input and would be the only "run-time" value;
 * all other values are known at "compile-time".
 * 
 * The output blocks are row-major, which means itemId.y is useful.
 * The output is also row-major within a block.
 * 
 * No divisions are used.
 * The second path for non-POT numItemsX only has a single multiply that could be done via u24.
 * 
 * See @bfi_method for explanation for first path.
**/
static uint2
Tile(uint2    itemId,
	 uint     lx, // = exact log2(ItemsPerTileX)
	 uint     ly, // = exact log2(ItemsPerTileY)
	 uint32_t numItemsX,
	 uint32_t numItemsY = 1u << 16) // optional, might allow for a smaller selectorMask immediate in the first path
{
	assert(numItemsX >= 2);
	assert((numItemsX & mask(lx)) == 0);

	// First path: used with numItemsX is known at compile time and is a power-of-2:
	if ((numItemsX & (numItemsX - 1)) == 0)
	{
		assert(numItemsY >= 2);
		const uint c = bsr(numItemsX);         // exact log2
		const uint d = bsr(numItemsY - 1) + 1; // ceil log2

		// RDNA shader ISA needs another trailing DWORD for immediates outside the range [-16, 64].
		// Prefer whatever would not need a 32-bit literal, or the zext version if both would:
		const int32_t  sext = (1 - numItemsX) | mask(lx); // 1 - numItemsX == ~mask(c) == -1<<c;
		const uint32_t zext = sext & mask(d + c - ly);
		const uint32_t selectorMask = (zext <= 64u || sext < -16) ? zext : sext;

		uint2 result;

		uint32_t flat        = itemId.y << c | itemId.x; // shl_or, could also be shl_add since carry-less
		uint32_t flat_shr_lx = flat >> lx; // ushr
		uint32_t flat_shr_ly = flat >> ly; // when lx == ly, 4 instructions instead of 5 total
		result.y             = v_bfi_b32(mask(ly), flat_shr_lx, itemId.y);
		result.x             = v_bfi_b32(selectorMask, itemId.x, flat_shr_ly);

		return result;

		// Thinking about this way to much maybe...
		/**
		 * if (lx == ly) // can be done in 4 instructions instead of 5:
		 * {
		 *     if (c >= (lx + ly)) // (numItems.x % (itemsPerBlock.x * itemsPerBlock.y)) == 0:
		 *     {
		 *         // Longer distances between dependent instructions,
		 *         // but uses 1 more register; probably not worth caring about(?), since another condition to check:
		 *         x_shr_lx    = x >> lx;
		 *         flat_shr_lk = y << (c - lx) | x_shr_lx; // == (flat >> lx) == (flat == ly) // waits[-1]
		 *         result.y    = v_bfi(mask(l), x_shr_lx, itemId.y);                          // waits[-2], shorter wait
		 *         result.x    = v_bfi(selmask, itemId.x, flat_shr_lk);                       // waits[-2]
		 *     }
		 *     else
		 *     {
		 *         // Shorter distances between dependent instructions,
		 *         // but uses 1 less register (flat not needed after flat_shr_lk calculated):
		 *         flat        = x << c | y;
		 *         flat_shr_lk = flat >> lx;               // == (flat >> lx) == (flat == ly) // waits[-1]
		 *         result.y    = v_bfi(mask(l), flat_shr_lk, itemId.y);                       // waits[-1], longer wait
		 *         result.x    = v_bfi(selmask, itemId.x, flat_shr_lk);                       // waits[-2]
		 *     }
		 * }
		**/
	}
	else // numItemsX is a runtime value or not a power-of-2:
	{
		uint2 idInBlock;
		uint2 blockId;

		uint itemIdFlatInThisRowOfTiles = (itemId.y & mask(ly)) * numItemsX + itemId.x; // and, umad24

		blockId.y = itemId.y >> ly;                                                     // ushr
		blockId.x = itemIdFlatInThisRowOfTiles >> (lx + ly);                            // ushr

		idInBlock.x = itemIdFlatInThisRowOfTiles & mask(lx);                            // and
		idInBlock.y = ubfe(itemIdFlatInThisRowOfTiles, lx, ly);                         // ubfe

		return (blockId << uint2{ lx, ly }) | idInBlock; // shl_or (or shl_add, carry-less) for each component
	}
}

static uint2
TileReferenceVersion(
	uint2 itemId,
	uint  itemsPerBlockX,
	uint  itemsPerBlockY,
	uint  numItemsX,
	uint  numItemsY)
{
	const uint blocksPerRow = numItemsX / itemsPerBlockX;
	const uint itemsPerBlockFlat = itemsPerBlockX * itemsPerBlockY;

	uint2 idInBlock;
	uint2 blockId;

	const uint itemIdFlat = itemId.y * numItemsX + itemId.x;
	const uint blockIdFlat = itemIdFlat / itemsPerBlockFlat;

	idInBlock.x =  itemIdFlat % itemsPerBlockX;
	idInBlock.y = (itemIdFlat / itemsPerBlockX) % itemsPerBlockY;

	blockId.x = blockIdFlat % blocksPerRow;
	blockId.y = blockIdFlat / blocksPerRow;

	return blockId * uint2{ itemsPerBlockX, itemsPerBlockY } + idInBlock;
}

#include <string.h>
#include <stdlib.h>

static void TestTile()
{
	{
		static const uint2 LogItemsPerBlock[] = {
			{ 1, 1 }, // 2x2 (quads)
			{ 2, 2 }, // 4x4
			{ 3, 3 }, // 8x8
			{ 3, 2 }, // 8x4
			{ 2, 3 }, // 8x3
		};
		static const uint2 NumBlocks[] = {
			{ 1, 1 },
			{ 1, 3 },
			{ 1, 4 },
			{ 3, 1 },
			{ 4, 1 },
			{ 4, 4 },
			{ 6, 4 },
			{ 4, 6 },
			{ 4, 2 },
			{ 2, 4 },
		};

		enum { W = 64, H = 64 };
		uint16_t cameFrom[H][W];

		for (const uint2 logItemsPerBlock : LogItemsPerBlock) {
			for (const uint2 numBlocks : NumBlocks) {
				const uint2 numItems = numBlocks << logItemsPerBlock;
				memset(cameFrom, 0xff, sizeof cameFrom);
				ASSERT(W >= numItems.x && H >= numItems.y);

				printf("\nLogItemsPerBlock={%d,%d}, ItemsPerBlock={%d,%d}, NumItems={%d,%d}\n",
					logItemsPerBlock.x, logItemsPerBlock.y,
					1u << logItemsPerBlock.x, 1u << logItemsPerBlock.y,
					numItems.x, numItems.y);
				fflush(stdout);

				for (uint y = 0; y < numItems.y; ++y) {
					for (uint x = 0; x < numItems.x; ++x) {
						const uint2 ref =
							TileReferenceVersion({ x, y }, 1u << logItemsPerBlock.x, 1u << logItemsPerBlock.y, numItems.x, numItems.y);
						const uint2 output = Tile({ x, y }, logItemsPerBlock.x, logItemsPerBlock.y, numItems.x, numItems.y);
						ASSERT(output.x < numItems.x && output.y < numItems.y);
						ASSERT(output.x == ref.x && output.y == ref.y);
						cameFrom[output.y][output.x] = y * numItems.x + x;
					}
				}

				for (uint y = 0; y < numItems.y; ++y) {
					for (uint x = 0; x < numItems.x; ++x) {
						const uint32_t flatId = cameFrom[y][x];
						ASSERT(flatId < W * H);
						printf("%4x", flatId);
					}
					putchar('\n');
				}
			}
		}
	}

#if 0
	printf("\n\nDoing stress test...");
	fflush(stdout);
	for (uint ly = 1; ly < 7; ++ly) {
		for (uint lx = 1; lx < 7; ++lx) {
			for (uint numBlocksY = 1; numBlocksY < 17; ++numBlocksY) {
				for (uint numBlocksX = 1; numBlocksX < 17; ++numBlocksX) {
					const uint2 numItems = uint2{ numBlocksX, numBlocksY } << uint2{ lx, ly };
					for (uint y = 0; y < numItems.y; ++y) {
						for (uint x = 0; x < numItems.x; ++x) {
							const uint2 ref = TileReferenceVersion({ x, y }, 1u << lx, 1u << ly, numItems.x, numItems.y);
							const uint2 output = Tile({ x, y },                    lx,       ly, numItems.x, numItems.y);
							if (ref.x != output.x || ref.y != output.y) {
								puts("failed.");
								exit(3);
							}
						}
					}
				}
			}
		}
	}
	puts("done.\n");
#endif
}

// z-order curve:
// input = dcba
// 
// output.x = ca 
// output.y = db 
//
// If would be doing more bits, could use s_bitreplicate?
static uint2 Interleave2x2BitsFromFlatId(uint dcba)
{
	ASSERT(dcba <= 0xf);
	// _ = 0
	const uint _dcb = dcba >> 1;
	const uint d__a = dcba & 0b1001;

	return {
		(_dcb  & 2) | (d__a & 1),
		(d__a >> 2) | (_dcb & 1)
	};
}

static uint2
SwizzleWithPotentialRemainderNoDivision(uint2 input, uint2 numItems)
{
	constexpr uint lx = 2; // log ItemsPerBlockX
	constexpr uint ly = 2; // log ItemsPerBlockY

	uint2 output = input;
	uint2 tileId;

	uint numFullTilesX = numItems.x >> lx;
	uint numTilesY     = numItems.y >> ly;
	tileId.y           = input.y    >> ly;

	if (tileId.y < numTilesY) {
		uint2 idInTile;
		uint flatIdInThisRowOfTiles = (input.y & mask(ly)) * numItems.x + input.x;
		tileId.x = flatIdInThisRowOfTiles >> (lx + ly);
		uint idInTileFlat = flatIdInThisRowOfTiles & mask(lx + ly);
		ASSERT(tileId.x <= numFullTilesX);

		if (tileId.x == numFullTilesX) {
			// The one partial tile on the right in this tile-row
			// arrange items in column-major order to avoid div/mod.
			idInTile.x = idInTileFlat >> ly;      // idInTileFlat / ItemsPerBlockY
			idInTile.y = idInTileFlat & mask(ly); // idInTileFlat % ItemsPerBlockY
		}
		else {
			// main section of full-tiles:
			static_assert(lx == ly && lx == 2, "");
			idInTile = Interleave2x2BitsFromFlatId(idInTileFlat);
		}
		output = (tileId << uint2{ lx, ly }) + idInTile;
	} // else: bottom section, do nothing.

	return output;
}

static uint2
SwizzleWithPotentialRemainderWithDivs(uint2 input, uint2 numItems)
{
	constexpr uint lx = 2; // log ItemsPerBlockX
	constexpr uint ly = 2; // log ItemsPerBlockY

	uint2 output;

	uint bottomSectionBeginY = (numItems.y & ~mask(ly));
	uint itemIdFlat = input.y * numItems.x + input.x;
	if (input.y < bottomSectionBeginY) {
		uint numTilesX = numItems.x >> lx;
		uint numTilesY = numItems.y >> ly;

		uint sideSectionFlatBegin = (numTilesX * numTilesY) << (lx + ly);

		if (itemIdFlat >= sideSectionFlatBegin) {
			uint localFlatId = itemIdFlat - sideSectionFlatBegin;
			uint sideWidth = numItems.x & mask(lx);
			// Arrange individual items in row-major order:
			output.y =  localFlatId / sideWidth;
			output.x = (localFlatId % sideWidth) + (numItems.x & ~mask(lx));
		}
		else {
			// main section of full-tiles:
			uint tileIdFlat   = itemIdFlat >> (lx + ly);
			uint idInTileFlat = itemIdFlat & mask(lx + ly);
			uint2 tileId = {
				tileIdFlat % numTilesX,
				tileIdFlat / numTilesX
			};
			static_assert(lx == ly && lx == 2, "");
			uint2 idInTile = Interleave2x2BitsFromFlatId(idInTileFlat);
			output = (tileId << uint2{ lx, ly }) + idInTile;
		}
	}
	else {
		// bottom section, arrange individual items in column-major order:
		uint localFlatId = itemIdFlat - bottomSectionBeginY * numItems.x;
		uint bottomHeight = numItems.y & mask(ly);
		output.x =  localFlatId / bottomHeight;
		output.y = (localFlatId % bottomHeight) + bottomSectionBeginY;
	}

	return output;
}

// TODO(?): Hybrid method like SwizzleWithPotentialRemainderNoDivision,
// but row-major in the potentially partial tile in each row of tiles,
// and column-major in the bottom section. Gets ride of the div in the main section,
// where most of the items should be.
//
// Another idea: could also pass in the magic values to do a divide via mul_hi, shifts, adds, etc.

static void
TestSwizzleWithPotentialRemainder()
{
	enum { W = 64, H = 64 };
	uint16_t cameFrom[H][W];

	const uint2 numItems = { 11, 10 };
	ASSERT(W >= numItems.x && H >= numItems.y);


	for (int method = 0; method < 2; ++method) {
		printf("\nNumItems={%d,%d}, method=\"%s\"\n", numItems.x, numItems.y, method == 0 ? "with divs" : "no divs");
		fflush(stdout);

		memset(cameFrom, 0xff, sizeof cameFrom);

		for (uint y = 0; y < numItems.y; ++y) {
			for (uint x = 0; x < numItems.x; ++x) {
				const uint2 output = method == 0 ? SwizzleWithPotentialRemainderWithDivs(  { x,y }, numItems)
					                             : SwizzleWithPotentialRemainderNoDivision({ x,y }, numItems);
				ASSERT(output.x < numItems.x && output.y < numItems.y);
				cameFrom[output.y][output.x] = y * numItems.x + x;
			}
		}

		for (uint y = 0; y < numItems.y; ++y) {
			for (uint x = 0; x < numItems.x; ++x) {
				const uint32_t flatId = cameFrom[y][x];
				ASSERT(flatId < W * H);
				printf("%4x", flatId);
			}
			putchar('\n');
		}
	}
}
/*
NumItems={11,10}, method="with divs"
   0   1   4   5  10  11  14  15  40  41  42
   2   3   6   7  12  13  16  17  43  44  45
   8   9   c   d  18  19  1c  1d  46  47  48
   a   b   e   f  1a  1b  1e  1f  49  4a  4b
  20  21  24  25  30  31  34  35  4c  4d  4e
  22  23  26  27  32  33  36  37  4f  50  51
  28  29  2c  2d  38  39  3c  3d  52  53  54
  2a  2b  2e  2f  3a  3b  3e  3f  55  56  57
  58  5a  5c  5e  60  62  64  66  68  6a  6c
  59  5b  5d  5f  61  63  65  67  69  6b  6d

NumItems={11,10}, method="no divs"
   0   1   4   5  10  11  14  15  20  24  28
   2   3   6   7  12  13  16  17  21  25  29
   8   9   c   d  18  19  1c  1d  22  26  2a
   a   b   e   f  1a  1b  1e  1f  23  27  2b
  2c  2d  30  31  3c  3d  40  41  4c  50  54
  2e  2f  32  33  3e  3f  42  43  4d  51  55
  34  35  38  39  44  45  48  49  4e  52  56
  36  37  3a  3b  46  47  4a  4b  4f  53  57
  58  59  5a  5b  5c  5d  5e  5f  60  61  62
  63  64  65  66  67  68  69  6a  6b  6c  6d
*/

int main(int argc, char **argv)
{
	TestTile();
	puts("\n\n--- potentially non-divisible by tile size: ---\n");
	TestSwizzleWithPotentialRemainder();
	return 0;
}


/*
@bfi_method:

General formula for {lx, ly} == log2(ItemsPerBlock.xy), as long as NumItems.x is a POT.

NOTE: output blocks in row-major order

const c = exact log2(NumItemsX), c >= lx, though c == lx is a nop.

itemIdFlat                  =  itemId.y             << c | itemId.x;
itemIdFlatInThisRowOfBlocks = (itemId.y & mask(ly)) << c | itemId.x; // not used

// Y axis: -------------------------------------------------------------------------------------------------------------

blockId.y   = itemId.y >> ly;           // itemId.y / ItemsPerBlockY

idInBlock.y = ubfe(itemIdFlat, lx, ly); // (itemIdFlat / ItemsPerBlockX) % ItemsPerBlockY
		if (c >= (lx + ly))             // (numItems.x % (itemsPerBlock.x * itemsPerBlock.y)) == 0
			= ubfe(itemId.x,   lx, ly);
		else
			= ubfe(itemId.y << c      | itemId.x,    lx,        ly);
			=     (itemId.y << c      | itemId.x) >> lx  & mask(ly);
			=     (itemId.y << (c-lx) | itemId.x  >> lx) & mask(ly);
			    
result.y    = blockId.y << ly        | idInBlock.y; // blockId.y*ItemsPerBlockY + idInBlock.y;
            = (itemId.y >> ly) << ly | (itemIdFlat >> lx & mask(ly));
			= (itemId.y & ~mask(ly)) | ...
			= (itemId.y & ~mask(ly)) | ...
		if (c >= (lx + ly))             // (numItems.x % (itemsPerBlock.x * itemsPerBlock.y)) == 0
			= v_bfi(mask(ly), itemId.x   >> lx, itemId.y);
		else
		    = v_bfi(mask(ly), itemIdFlat >> lx, itemId.y);


// X axis: -------------------------------------------------------------------------------------------------------------

const LogBlocksPerRow = (c - lx); // NumItemsX / ItemsPerBlockX

blockId.x   = itemIdFlat >> (lx + ly) & mask(LogBlocksPerRow);
            = ubfe(itemIdFlat, lx + ly, c - lx);

idInBlock.x = itemId.x & mask(lx);

result.x = blockId.x                                                                   << lx  | idInBlock.x;
		 = ubfe( itemIdFlat,                     lx + ly,        c - lx )              << lx  | (itemId.x & mask(lx));
		 =     ( itemIdFlat                  >> (lx + ly) & mask(c - lx))              << lx  | (itemId.x & mask(lx));
		 // change bfe via shift-then-mask to mask-then-shift:
		 =     ( (itemIdFlat       & (mask(c - lx) << (lx + ly)) ) >> (lx + ly) )      << lx  | (itemId.x & mask(lx));
		 =     ( (itemIdFlat       & (mask(c - lx) << (lx + ly)) ) >>       ly) ) & ~mask(lx) | (itemId.x & mask(lx));
		 =       itemIdFlat >> ly  & (mask(c - lx) <<  lx)                        & ~mask(lx) | (itemId.x & mask(lx));
		 =       itemIdFlat >> ly  &  mask(c) &   mask(lx)                        & ~mask(lx) | (itemId.x & mask(lx));
		 
		 = v_bfi(        mask(lx), itemId.x, itemIdFlat >> ly & mask(c) & mask(lx));
		 = v_bfi(        mask(lx), itemId.x, itemIdFlat >> ly & mask(c)); // since the v_bfi will write over bits [lx-1 : 0].

		 = v_bfi(-1<<c | mask(lx), itemId.x, itemIdFlat >> ly);           // since bits at pos >= c in itemId.x are 0, Note[2]


Note[1]: If lx == ly, (itemId.x >> lx) used for result.y and (itemId.x >> ly) used for result.x are the same
and can be reused, making 4 instructions total to compute both result.x and result.y.

Note[2]: The v_bfi's sign-extended selectorMask of (-1<<c | mask(lx)) could also be AND'd with mask(d + c - ly)
to yield a zero-extended immediate if that would result in a smaller immediate encoding.


Some examples of where the bits of x and y go:

v[2] means the bit-value at bit 2 of v

"::" means concat bits:

	NumItems = { 8, 8 } and { lx, ly } == { 1, 1 }: // ItemsPerBlock = { 2, 2 } (quads)
	{
		out.y = y[2] :: y[1] :: x[1];
		out.x = y[0] :: x[2] :: x[0];
	}
*/

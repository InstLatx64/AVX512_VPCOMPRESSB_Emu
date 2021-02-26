// AVX512BW_VPCOMPRESSB_Emu
//

#include "stdafx.h"

using namespace std;

#define RETRY			1000

CPU_Props cpu_props;

extern "C" size_t __fastcall Test_VPCOMPRESSD_Asm(const char* src, char* dst, size_t len);
extern "C" size_t __fastcall Test_VPCOMPRESSB_Asm(const char* src, char* dst, size_t len);
extern "C" size_t __fastcall Test_VPCOMPRESSB_ymm_Asm(const char* src, char* dst, size_t len);
extern "C" size_t __fastcall Test_VPCOMPRESSB2_Asm(const char* src, char* dst, size_t len);
extern "C" size_t __fastcall Test_VPCOMPRESSB4_Asm(const char* src, char* dst, size_t len);
extern "C" size_t __fastcall Test_VPCOMPRESSB2_ymm_Asm(const char* src, char* dst, size_t len);

typedef size_t(__fastcall* TEST_PTR)(const char*, char*, size_t);

size_t remove_spaces_scalar(const char* src, char* dst, size_t n);
size_t remove_spaces_avx512vbmi(const char* src, char* dst, size_t n);
size_t remove_spaces_avx512bw(const char* src, char* dst, size_t n);
size_t remove_spaces_avx512vbmi_zach(const char* src, char* dst, size_t n);
size_t despace_avx2_vpermd(const char* src_void, char* dst_void, size_t length);
size_t despace_sse41_lut_216(const char* src_void, char* dst_void, size_t length);
size_t despace_avx2_vpermd2(const char* src_void, char* dst_void, size_t length);
size_t sse4_despace_branchless(const char * src_void, char * dst_void, size_t howmany);
size_t sse4_despace_branchless_u4(const char* src_void, char* dst_void, size_t length);
size_t sse4_despace_branchless_u2(const char* src_void, char* dst_void, size_t length);
size_t avx2_despace_branchless(const char* src_void, char* dst_void, size_t length);
size_t avx2_despace_branchless_u2(const char* src_void, char* dst_void, size_t length);
size_t despace_branchless(unsigned char* dst_void, unsigned char * src_void, size_t length);

typedef struct {
	const char 	name[32];
	const char 	isaName[16];
	int			length;
	int			retry;
	TEST_PTR	timed;
	ISAs		isa;
	bool		test;
} methods;

methods m[] = {
	{"Scalar              ",	"X64        ",	 1, 500,	remove_spaces_scalar,				ISA_RDTSC,			true},
	{"SSE4_DESPACE_BRANCHLESS        ",	"SSSE3      ",	16,	RETRY, sse4_despace_branchless,				ISA_SSE41,					true},
	{"SSE4_DESPACE_BRANCHLESS_U4     ",	"SSSE3      ",	16,	RETRY, sse4_despace_branchless_u4,			ISA_SSE41,					true},
	{"SSE4_DESPACE_BRANCHLESS_U2     ",	"SSSE3      ",	16,	RETRY, sse4_despace_branchless_u2,			ISA_SSE41,					true},
	{"SSE41                          ",	"SSE41      ",	16,	RETRY,	despace_sse41_lut_216,				ISA_SSE41,			true},
	{"AVX2_DESPACE_BRANCHLESS        ",	"AVX2       ",	32,	RETRY, avx2_despace_branchless,				ISA_AVX2,					true},
	{"AVX2_DESPACE_BRANCHLESS_U2     ",	"AVX2       ",	32,	RETRY, avx2_despace_branchless_u2,			ISA_AVX2,					true},
	{"AVX2_VPERMD         ",	"AVX2       ",	32, RETRY,	despace_avx2_vpermd,				ISA_AVX2,			true},
	{"AVX2_VPERMD2                   ",	"AVX2       ",	32,	RETRY,	despace_avx2_vpermd2,				ISA_AVX2,			true},
	{"AVX512BW_VPCOMPRESSD           ",	"AVX512BW   ",	64,	RETRY,	remove_spaces_avx512bw,				ISA_AVX512BW,		true},
	{"AVX512BW_VPCOMPRESSD_Asm       ",	"AVX512BW   ",	64,	RETRY,	Test_VPCOMPRESSD_Asm,				ISA_AVX512BW,		true},
	{"AVX512VBMI_Intrin              ",	"AVX512VBMI ",	64, RETRY,	remove_spaces_avx512vbmi,			ISA_AVX512VBMI,		true},
	{"AVX512VBMI_Zach                ",	"AVX512VBMI ",	64, RETRY,	remove_spaces_avx512vbmi_zach,		ISA_AVX512VBMI,		true},
	{"AVX512VBMI2_VPCOMPRESSB_Asm    ",	"AVX512VBMI2",	64,	RETRY,	Test_VPCOMPRESSB_Asm,				ISA_AVX512_VBMI2,	true},
	{"AVX512VBMI2_VPCOMPRESSB2_Asm   ",	"AVX512VBMI2",	128,RETRY,	Test_VPCOMPRESSB2_Asm,				ISA_AVX512_VBMI2,	true},
	{"AVX512VBMI2_VPCOMPRESSB4_Asm   ",	"AVX512VBMI2",	256,RETRY,	Test_VPCOMPRESSB4_Asm,				ISA_AVX512_VBMI2,	true},
	{"AVX512VBMI2_VPCOMPRESSB_ymm_Asm",	"AVX512VBMI2",	32,	RETRY,	Test_VPCOMPRESSB_ymm_Asm,			ISA_AVX512_VBMI2,	true},
	{"AVX512VBMI2_VPCOMPRESSB2ymm_Asm",	"AVX512VBMI2",	64,	RETRY,	Test_VPCOMPRESSB2_ymm_Asm,			ISA_AVX512_VBMI2,	true},
};

//credit: Wojciech Mu³a http://0x80.pl/notesen/2019-01-05-avx512vbmi-remove-spaces.html
size_t remove_spaces_scalar(const char* src, char* dst, size_t n) {
	char* startdst = dst;
	for (size_t i = 0; i < n; i++) {
		if (src[i] == ' ' || src[i] == '\r' || src[i] == '\n')
			continue;

		*dst++ = src[i];
	}

	return dst - startdst;
}

size_t remove_spaces_avx512vbmi(const char* src_void, char* dst_void, size_t length) {
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	//assert(n % 64 == 0);

	// values 0..63
	const __m512i no_gaps_indices = _mm512_setr_epi32(
		0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c,
		0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c,
		0x23222120, 0x27262524, 0x2b2a2928, 0x2f2e2d2c,
		0x33323130, 0x37363534, 0x3b3a3938, 0x3f3e3d3c
	);

	const __m512i ones = _mm512_set1_epi8(1);
	const __m512i NL = _mm512_set1_epi8('\n');
	const __m512i CR = _mm512_set1_epi8('\r');
	const __m512i spaces = _mm512_set1_epi8(' ');

	size_t len;
	for (size_t i = 0; i < (length & ~0x3f); i += 64) {

		const __m512i input = _mm512_loadu_si512((const __m512i*)(src));
		__m512i output;

		uint64_t mask = _mm512_cmpeq_epi8_mask(input, spaces)
			| _mm512_cmpeq_epi8_mask(input, NL)
			| _mm512_cmpeq_epi8_mask(input, CR);

		if (mask) {
			len = 64 - _mm_popcnt_u64(mask);
			__m512i indices = no_gaps_indices;
			__m512i increment = ones;

			uint64_t first;
			uint64_t prev;

			first = (mask & (0 - mask));
			prev = first;
			mask ^= first;
			mask >>= 1;

			while (mask) {
				const uint64_t curr = (mask & (0 - mask));
				mask ^= curr;
				mask >>= 1;

				if (prev == curr) {
					increment = _mm512_add_epi8(increment, ones);
					prev = curr;
				}
				else {
					indices = _mm512_mask_add_epi8(indices, ~(first - 1), indices, increment);

					first = curr;
					prev = curr;
					increment = ones;
				}
			}

			indices = _mm512_mask_add_epi8(indices, ~(first - 1), indices, increment);

			output = _mm512_permutexvar_epi8(indices, input);
		}
		else {
			output = input;
			len = 64;
		}

		_mm512_storeu_si512((__m512i*)(dst), output);
		dst += len;
		src += 64;
	}

	dst += despace_branchless(dst, src, length & 63);
	return (size_t)(dst - ((uint8_t*)dst_void));
}

//idea: https://branchfree.org/2018/05/22/bits-to-indexes-in-bmi2-and-avx-512/
//geofflangdale
//vpcmpeqb -> pshufb credit aqrit https://gist.github.com/aqrit/6e73ca6ff52f72a2b121d584745f89f3
#define VPCOMPRESSD_Core(addr) \
	_mm512_mask_cvtepi32_storeu_epi8(dst, full, _mm512_mask_compress_epi32(zero, __mmask16(mask), _mm512_cvtepu8_epi32(*(const __m128i*)(addr))));	\
	_mm512_mask_cvtepi32_storeu_epi8(dst + _mm_popcnt_u64((WORD)mask), full, _mm512_mask_compress_epi32(zero, __mmask16(mask >> 16), _mm512_cvtepu8_epi32(*(const __m128i*)(addr + 16)))); \
	_mm512_mask_cvtepi32_storeu_epi8(dst + _mm_popcnt_u64((DWORD)mask), full, _mm512_mask_compress_epi32(zero, __mmask16(mask >> 32), _mm512_cvtepu8_epi32(*(const __m128i*)(addr + 32)))); \
	_mm512_mask_cvtepi32_storeu_epi8(dst + _mm_popcnt_u64(mask << 16), full, _mm512_mask_compress_epi32(zero, __mmask16(mask >> 48), _mm512_cvtepu8_epi32(*(const __m128i*)(addr + 48))));

size_t remove_spaces_avx512bw(const char* src, char* dst, size_t n) {
	char* startdst = dst;

	const __m512i zero = _mm512_setzero_si512();
	const __m512i lut_cntrl2	= _mm512_broadcast_i32x4(_mm_setr_epi8(' ', 0, 0, 0, 0, 0, 0, 0, 0, '\t', '\n', 0, 0, '\r', 0, 0));
	const __mmask16 full = 0xffff;
	size_t i = 0;
	uint64_t mask;
	for (size_t i = 0; i < (n & ~0x3f); i += 64) {
		const __m512i input = _mm512_loadu_si512((const __m512i*)(src + i));
		mask = _mm512_cmpneq_epi8_mask(_mm512_shuffle_epi8(lut_cntrl2, input), input);

		VPCOMPRESSD_Core(src + i)
		dst += _mm_popcnt_u64(mask);
	}
	if ((n & 0x3f) != 0) {
		mask = _bzhi_u64(~0ULL, (int)(n & 0x3f));
		const char * addr = src + (n & ~0x3f);
		const __m512i input = _mm512_maskz_loadu_epi8(mask, (const __m512i*)(addr));
		mask = _mm512_mask_cmpneq_epi8_mask(mask, _mm512_shuffle_epi8(lut_cntrl2, input), input);

		VPCOMPRESSD_Core(addr)
		dst += _mm_popcnt_u64(mask);
	}
	return dst - startdst;
}

size_t despace_branchless(unsigned char* dst_void, unsigned char * src_void, size_t length)
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	for (; length != 0; length--) {
		uint8_t c = *src++;
		*dst = c;
		dst += !!((c != 0x20) && (c != 0x0A) && (c != 0x0D) && (c != 0x09));
	}
	return (size_t)(dst - ((uint8_t*)dst_void));
}


//copy from http://0x80.pl/notesen/2019-01-05-avx512vbmi-remove-spaces.html
//credit https://github.com/zwegner/
size_t remove_spaces_avx512vbmi_zach(const char* src_void, char* dst_void, size_t length) {
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
	//assert(n % 64 == 0);

	const __m512i NL = _mm512_set1_epi8('\n');
	const __m512i CR = _mm512_set1_epi8('\r');
	const __m512i spaces = _mm512_set1_epi8(' ');

	uint64_t index_masks[6] = {
		0xaaaaaaaaaaaaaaaa,
		0xcccccccccccccccc,
		0xf0f0f0f0f0f0f0f0,
		0xff00ff00ff00ff00,
		0xffff0000ffff0000,
		0xffffffff00000000,
	};

	const __m512i index_bits[6] = {
		_mm512_set1_epi8(1),
		_mm512_set1_epi8(2),
		_mm512_set1_epi8(4),
		_mm512_set1_epi8(8),
		_mm512_set1_epi8(16),
		_mm512_set1_epi8(32),
	};

	size_t len;
	for (size_t i = 0; i < (length & ~0x3f); i += 64) {
		const __m512i input = _mm512_loadu_si512((const __m512i*)(src));
		__m512i output;

		uint64_t mask = _mm512_cmpeq_epi8_mask(input, spaces)
			| _mm512_cmpeq_epi8_mask(input, NL)
			| _mm512_cmpeq_epi8_mask(input, CR);

		if (mask) {
			//len = 64 - __builtin_popcountll(mask);
			len = 64 - _mm_popcnt_u64(mask);
			mask = ~mask;
			__m512i indices = _mm512_set1_epi8(0);
			for (size_t index = 0; index < 6; index++) {
				uint64_t m = _pext_u64(index_masks[index], mask);
				indices = _mm512_mask_add_epi8(indices, m, indices, index_bits[index]);
			}

			output = _mm512_permutexvar_epi8(indices, input);
		}
		else {
			output = input;
			len = 64;
		}

		_mm512_storeu_si512((__m512i*)(dst), output);
		dst += len;
		src += 64;
	}

	dst += despace_branchless(dst, src, length & 63);
	return (size_t)(dst - ((uint8_t*)dst_void));
}

//credit : Michael Howard
size_t despace_avx2_vpermd(const char* src_void, char* dst_void, size_t length)
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
	const __m256i lut_cntrl = _mm256_setr_epi8(
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, -1, -1, 0x00, 0x00, -1, 0x00, 0x00,
		//
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, -1, -1, 0x00, 0x00, -1, 0x00, 0x00
	);
	// hi and lo dwords are combined into a qword using xor
	// qword = dword(lo[0b000] ^ hi[id]), dword(lo[id] ^ hi[0b000])
	// byte_0 of lo always gets shifted off unless the id is 0b111
	// so if lo[0b111] ^ hi[0b000] then byte_0 must be zero
	const __m256i lut_lo = _mm256_set_epi32(
		0x03020107, // 111
		0x03020000, // 110
		0x03010000, // 101
		0x03000000, // 100
		0x02010000, // 011
		0x02000000, // 010
		0x01000000, // 001
		0x00000000  // 000
	);
	const __m256i lut_hi = _mm256_set_epi32(
		0x07060504, // 111
		0x00070605, // 110
		0x00070604, // 101
		0x00000706, // 100
		0x00070504, // 011
		0x00000705, // 010
		0x00000704, // 001
		0x00000007  // 000
	);
	const __m256i mask_20 = _mm256_set1_epi8(0x20);
	const __m256i mask_70 = _mm256_set1_epi8(0x70);
	const __m256i mask_sad = _mm256_set1_epi64x(0x80A0908884828180);
	__m256i temp_08 = _mm256_srli_epi64(mask_20, 2);
	const __m256i mask_range = _mm256_blend_epi32(temp_08, _mm256_setzero_si256(), 0x33); // 0x08080808080808080000000000000000
	const __m256i mask_shift = _mm256_blend_epi32(temp_08, _mm256_setzero_si256(), 0xAA); // 0x00000000080808080000000008080808
	for (; length >= 32; length -= 32) {
		__m256i r0, r1, r2, r3, r4;
		r0 = _mm256_loadu_si256((__m256i*)src); // asrc
		src += 32;
		//
		r1 = _mm256_adds_epu8(mask_70, r0);
		r2 = _mm256_cmpeq_epi8(mask_20, r0);
		r1 = _mm256_shuffle_epi8(lut_cntrl, r1);
		r1 = _mm256_or_si256(r1, r2); // bytemask of spaces
		//
		r2 = _mm256_andnot_si256(r1, mask_sad);
		r1 = _mm256_and_si256(r1, mask_shift);
		r2 = _mm256_sad_epu8(r2, _mm256_setzero_si256()); // non-space -> bitmap[5:0],population count[15:7]
		r1 = _mm256_sad_epu8(r1, _mm256_setzero_si256()); // bit shift amount to drop space bytes on low end of qword
		r3 = _mm256_permutevar8x32_epi32(lut_lo, r2); // lo
		r2 = _mm256_slli_epi64(r2, 29); // move hi index to 2nd dword
		r4 = _mm256_permutevar8x32_epi32(lut_hi, r2); // hi
		r2 = _mm256_srli_epi64(r2, 36); // popcnt
		r3 = _mm256_or_si256(r3, mask_range); // change range of hi_qw from 7:0 to 15:8
		r3 = _mm256_xor_si256(r3, r4); // merge lo & hi dwords
		r3 = _mm256_srlv_epi64(r3, r1);
		//
		r0 = _mm256_shuffle_epi8(r0, r3);
		*((uint64_t*)dst) = _mm256_extract_epi64(r0, 0);
		dst += _mm256_extract_epi64(r2, 0);
		*((uint64_t*)dst) = _mm256_extract_epi64(r0, 1);
		dst += _mm256_extract_epi64(r2, 1);
		*((uint64_t*)dst) = _mm256_extract_epi64(r0, 2);
		dst += _mm256_extract_epi64(r2, 2);
		*((uint64_t*)dst) = _mm256_extract_epi64(r0, 3);
		dst += _mm256_extract_epi64(r2, 3);
	}
	dst += despace_branchless(dst, src, length & 31);
	return (size_t)(dst - ((uint8_t*)dst_void));
}

//credit aqrit https://gist.github.com/aqrit/6e73ca6ff52f72a2b121d584745f89f3
size_t despace_sse41_lut_216(const char* src_void, char* dst_void, size_t length)
{
    uint8_t* src = (uint8_t*)src_void;
    uint8_t* dst = (uint8_t*)dst_void;

    if (length >= 16) {
        // table of control characters (space, tab, newline, carriage return)
        const __m128i lut_cntrl = _mm_setr_epi8(' ', 0, 0, 0, 0, 0, 0, 0, 0, '\t', '\n', 0, 0, '\r', 0, 0);

        // bits[4:0] = index -> ((trit_d * 0) + (trit_c * 9) + (trit_b * 3) + (trit_a * 1))
        // bits[15:7] = popcnt
        const __m128i sadmask = _mm_set1_epi64x(0x8080898983838181);

        // adding 8 to each shuffle index is cheaper than extracting the high qword
        const __m128i offset = _mm_cvtsi64_si128(0x0808080808080808);

        // shuffle control indices
        static const uint64_t table[27] = {
            0x0000000000000706, 0x0000000000070600, 0x0000000007060100, 
            0x0000000000070602, 0x0000000007060200, 0x0000000706020100, 
            0x0000000007060302, 0x0000000706030200, 0x0000070603020100, 
            0x0000000000070604, 0x0000000007060400, 0x0000000706040100,
            0x0000000007060402, 0x0000000706040200, 0x0000070604020100, 
            0x0000000706040302, 0x0000070604030200, 0x0007060403020100, 
            0x0000000007060504, 0x0000000706050400, 0x0000070605040100, 
            0x0000000706050402, 0x0000070605040200, 0x0007060504020100,
            0x0000070605040302, 0x0007060504030200, 0x0706050403020100
		};

        const uint8_t* end = &src[length & ~15];
        do {
            __m128i v = _mm_loadu_si128((__m128i*)src);
            src += 16;

            // detect spaces
            __m128i mask = _mm_cmpeq_epi8(_mm_shuffle_epi8(lut_cntrl, v), v);

            // shift w/blend: each word now only has 3 states instead of 4
            // which takes us from 128 possiblities, per qword, down to 27
            v = _mm_blendv_epi8(v, _mm_srli_epi16(v, 8), mask);

            // extract bitfields describing each qword: index, popcnt
            __m128i desc = _mm_sad_epu8(_mm_and_si128(mask, sadmask), sadmask);
            size_t lo_desc = (size_t)_mm_cvtsi128_si32(desc);
            size_t hi_desc = (size_t)_mm_extract_epi16(desc, 4);

            // load shuffle control indices from pre-computed table
            __m128i lo_shuf = _mm_loadl_epi64((__m128i*)&table[lo_desc & 0x1F]);
            __m128i hi_shuf = _mm_or_si128(_mm_loadl_epi64((__m128i*)&table[hi_desc & 0x1F]), offset);

            // store an entire qword then advance the pointer by how ever
            // many of those bytes are actually wanted. Any trailing
            // garbage will be overwritten by the next store.
            // note: little endian byte memory order
            _mm_storel_epi64((__m128i*)dst, _mm_shuffle_epi8(v, lo_shuf));
            dst += (lo_desc >> 7);
            _mm_storel_epi64((__m128i*)dst, _mm_shuffle_epi8(v, hi_shuf));
            dst += (hi_desc >> 7);
        } while (src != end);
    }

    // tail loop
    length &= 15;
    if (length != 0) {
        const uint64_t bitmap = 0xFFFFFFFEFFFFC1FF;
        do {
            uint64_t c = *src++;
            *dst = (uint8_t)c;
            dst += ((bitmap >> c) & 1) | ((c + 0xC0) >> 8);
        } while (--length);
    }

    // return pointer to the location after the last element in dst
    return (size_t)(dst - ((uint8_t*)dst_void));
}

//credit aqrit https://gist.github.com/aqrit/6e73ca6ff52f72a2b121d584745f89f3
// no lookup tables
// probably needs improvment...
size_t despace_avx2_vpermd2(const char* src_void, char* dst_void, size_t length)
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;

	const __m256i lut_cntrl2	= _mm256_broadcastsi128_si256(_mm_setr_epi8(' ', 0, 0, 0, 0, 0, 0, 0, 0, '\t', '\n', 0, 0, '\r', 0, 0));
	const __m256i permutation_mask = _mm256_set1_epi64x( 0x0020100884828180 );
	const __m256i invert_mask = _mm256_set1_epi64x( 0x0020100880808080 ); 
	const __m256i zero = _mm256_setzero_si256();
	const __m256i fixup = _mm256_set_epi32(
		0x08080808, 0x0F0F0F0F, 0x00000000, 0x07070707,
		0x08080808, 0x0F0F0F0F, 0x00000000, 0x07070707
	);
	const __m256i lut = _mm256_set_epi32(
		0x04050607, // 0x03020100', 0x000000'07
		0x04050704, // 0x030200'00, 0x0000'0704
		0x04060705, // 0x030100'00, 0x0000'0705
		0x04070504, // 0x0300'0000, 0x00'070504
		0x05060706, // 0x020100'00, 0x0000'0706
		0x05070604, // 0x0200'0000, 0x00'070604
		0x06070605, // 0x0100'0000, 0x00'070605
		0x07060504  // 0x00'000000, 0x'07060504
	);

	// hi bits are ignored by pshufb, used to reject movement of low qword bytes
	const __m256i shuffle_a = _mm256_set_epi8(
		0x7F, 0x7E, 0x7D, 0x7C, 0x7B, 0x7A, 0x79, 0x78, 0x07, 0x16, 0x25, 0x34, 0x43, 0x52, 0x61, 0x70,
		0x7F, 0x7E, 0x7D, 0x7C, 0x7B, 0x7A, 0x79, 0x78, 0x07, 0x16, 0x25, 0x34, 0x43, 0x52, 0x61, 0x70
	);

	// broadcast 0x08 then blendd...
	const __m256i shuffle_b = _mm256_set_epi32(
		0x08080808, 0x08080808, 0x00000000, 0x00000000,
		0x08080808, 0x08080808, 0x00000000, 0x00000000
	);

	for( uint8_t* end = &src[(length & ~31)]; src != end; src += 32){
		__m256i r0,r1,r2,r3,r4;
		unsigned int s0,s1;

		r0 = _mm256_loadu_si256((__m256i *)src); // asrc

		// detect spaces
		r1 = _mm256_cmpeq_epi8(_mm256_shuffle_epi8(lut_cntrl2, r0), r0);

		r2 = _mm256_sad_epu8(zero, r1);
		s0 = _mm256_movemask_epi8(r1);
		r1 = _mm256_andnot_si256(r1, permutation_mask);

		r1 = _mm256_sad_epu8(r1, invert_mask); // index_bitmap[0:5], low32_spaces_count[7:15]

		r2 = _mm256_shuffle_epi8(r2, zero);

		r2 = _mm256_sub_epi8(shuffle_a, r2); // add space cnt of low qword
		s0 = ~s0;

		r3 = _mm256_slli_epi64(r1, 29); // move top part of index_bitmap to high dword
		r4 = _mm256_srli_epi64(r1, 7); // number of spaces in low dword 

		r4 = _mm256_shuffle_epi8(r4, shuffle_b);
		r1 = _mm256_or_si256(r1, r3);

		r1 = _mm256_permutevar8x32_epi32(lut, r1);
		s1 = _mm_popcnt_u32(s0);
		r4 = _mm256_add_epi8(r4, shuffle_a);
		s0 = s0 & 0xFFFF; // isolate low oword

		r2 = _mm256_shuffle_epi8(r4, r2);
		s0 = _mm_popcnt_u32(s0);

		r2 = _mm256_max_epu8(r2, r4); // pin low qword bytes

		r1 = _mm256_xor_si256(r1, fixup);

		r1 = _mm256_shuffle_epi8(r1, r2); // complete shuffle mask

		r0 = _mm256_shuffle_epi8(r0, r1); // despace!

		_mm_storeu_si128((__m128i*)dst, _mm256_castsi256_si128(r0));
		_mm_storeu_si128((__m128i*)&dst[s0], _mm256_extracti128_si256(r0,1));
		dst += s1;
	}
	dst += despace_branchless(dst, src, length & 31);
	return (size_t)(dst - ((uint8_t*)dst_void));
}

size_t sse4_despace_branchless(const char* src_void, char* dst_void, size_t length)
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
	size_t pos = 0;
	__m128i spaces = _mm_set1_epi8(' ');
	__m128i newline = _mm_set1_epi8('\n');
	__m128i carriage = _mm_set1_epi8('\r');
	size_t i = 0;
	for (; i + 15 < length; i += 16) {
		__m128i x = _mm_loadu_si128((const __m128i *)(src + i));
		__m128i xspaces = _mm_cmpeq_epi8(x, spaces);
		__m128i xnewline = _mm_cmpeq_epi8(x, newline);
		__m128i xcarriage = _mm_cmpeq_epi8(x, carriage);
		__m128i anywhite = _mm_or_si128(_mm_or_si128(xspaces, xnewline), xcarriage);
		uint64_t mask16 = _mm_movemask_epi8(anywhite);
		x = _mm_shuffle_epi8(x, *((__m128i *)despace_mask16 + (mask16 & 0x7fff)));
		_mm_storeu_si128((__m128i *)(dst + pos), x);
		pos += 16 - _mm_popcnt_u64(mask16);
	}
	for (; i < length; i++) {
		char c = src[i];
		if (c == '\r' || c == '\n' || c == ' ') {
		  continue;
		}
		dst[pos++] = c;
	}
	return pos;
}

static inline __m128i cleanm128(__m128i x, __m128i spaces, __m128i newline,
                                __m128i carriage, int *mask16) {
  __m128i xspaces = _mm_cmpeq_epi8(x, spaces);
  __m128i xnewline = _mm_cmpeq_epi8(x, newline);
  __m128i xcarriage = _mm_cmpeq_epi8(x, carriage);
  __m128i anywhite = _mm_or_si128(_mm_or_si128(xspaces, xnewline), xcarriage);
  *mask16 = _mm_movemask_epi8(anywhite);
  return _mm_shuffle_epi8(
      x, _mm_loadu_si128((const __m128i *)despace_mask16 + (*mask16 & 0x7fff)));
}

//static inline size_t sse4_despace_branchless_u4(char *bytes, size_t howmany) {
size_t sse4_despace_branchless_u4(const char* src_void, char* dst_void, size_t length)
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
  size_t pos = 0;
  __m128i spaces = _mm_set1_epi8(' ');
  __m128i newline = _mm_set1_epi8('\n');
  __m128i carriage = _mm_set1_epi8('\r');
  size_t i = 0;
  for (; i + 64 - 1 < length; i += 64) {
    __m128i x1 = _mm_loadu_si128((const __m128i *)(src + i));
    __m128i x2 = _mm_loadu_si128((const __m128i *)(src + i + 16));
    __m128i x3 = _mm_loadu_si128((const __m128i *)(src + i + 32));
    __m128i x4 = _mm_loadu_si128((const __m128i *)(src + i + 48));

    int mask16;
    x1 = cleanm128(x1, spaces, newline, carriage, &mask16);
    _mm_storeu_si128((__m128i *)(dst + pos), x1);
    pos += 16 - _mm_popcnt_u32(mask16);

    x2 = cleanm128(x2, spaces, newline, carriage, &mask16);
    _mm_storeu_si128((__m128i *)(dst + pos), x2);
    pos += 16 - _mm_popcnt_u32(mask16);

    x3 = cleanm128(x3, spaces, newline, carriage, &mask16);
    _mm_storeu_si128((__m128i *)(dst + pos), x3);
    pos += 16 - _mm_popcnt_u32(mask16);

    x4 = cleanm128(x4, spaces, newline, carriage, &mask16);
    _mm_storeu_si128((__m128i *)(dst + pos), x4);
    pos += 16 - _mm_popcnt_u32(mask16);
  }
  for (; i + 16 - 1 < length; i += 16) {
    __m128i x = _mm_loadu_si128((const __m128i *)(src + i));
    int mask16;
    x = cleanm128(x, spaces, newline, carriage, &mask16);
    _mm_storeu_si128((__m128i *)(dst + pos), x);
    pos += 16 - _mm_popcnt_u32(mask16);
  }
  for (; i < length; i++) {
    char c = src[i];
    if (c == '\r' || c == '\n' || c == ' ') {
      continue;
    }
    dst[pos++] = c;
  }
  return pos;
}

//static inline size_t sse4_despace_branchless_u2(char *bytes, size_t howmany) {
size_t sse4_despace_branchless_u2(const char* src_void, char* dst_void, size_t length)
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
  size_t pos = 0;
  __m128i spaces = _mm_set1_epi8(' ');
  __m128i newline = _mm_set1_epi8('\n');
  __m128i carriage = _mm_set1_epi8('\r');
  size_t i = 0;
  for (; i + 32 - 1 < length; i += 32) {
    __m128i x1 = _mm_loadu_si128((const __m128i *)(src + i));
    __m128i x2 = _mm_loadu_si128((const __m128i *)(src + i + 16));
    int mask16;
    x1 = cleanm128(x1, spaces, newline, carriage, &mask16);
    _mm_storeu_si128((__m128i *)(dst + pos), x1);
    pos += 16 - _mm_popcnt_u32(mask16);

    x2 = cleanm128(x2, spaces, newline, carriage, &mask16);
    _mm_storeu_si128((__m128i *)(dst + pos), x2);
    pos += 16 - _mm_popcnt_u32(mask16);
  }
  for (; i + 16 - 1 < length; i += 16) {
    __m128i x = _mm_loadu_si128((const __m128i *)(src + i));
    int mask16;
    x = cleanm128(x, spaces, newline, carriage, &mask16);
    _mm_storeu_si128((__m128i *)(dst + pos), x);
    pos += 16 - _mm_popcnt_u32(mask16);
  }
  for (; i < length; i++) {
    char c = src[i];
    if (c == '\r' || c == '\n' || c == ' ') {
      continue;
    }
    dst[pos++] = c;
  }
  return pos;
}

static inline __m256i cleanm256(__m256i x, __m256i spaces, __m256i newline,
                                __m256i carriage, unsigned int *mask1,
                                unsigned int *mask2) {
  __m256i xspaces = _mm256_cmpeq_epi8(x, spaces);
  __m256i xnewline = _mm256_cmpeq_epi8(x, newline);
  __m256i xcarriage = _mm256_cmpeq_epi8(x, carriage);
  __m256i anywhite =
      _mm256_or_si256(_mm256_or_si256(xspaces, xnewline), xcarriage);
  unsigned int mask32 = _mm256_movemask_epi8(anywhite);
  unsigned int maskhigh = (mask32) >> 16;
  unsigned int masklow = (mask32)&0xFFFF;
  assert(maskhigh < (1 << 16));
  assert(masklow < (1 << 16));
  *mask1 = masklow;
  *mask2 = maskhigh;
  __m256i mask = _mm256_loadu2_m128i((const __m128i *)despace_mask16 + (maskhigh & 0x7fff),
                                     (const __m128i *)despace_mask16 + (masklow & 0x7fff));
  return _mm256_shuffle_epi8(x, mask);
}

size_t avx2_despace_branchless(const char* src_void, char* dst_void, size_t length)
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
  size_t pos = 0;
  __m128i spaces = _mm_set1_epi8(' ');
  __m128i newline = _mm_set1_epi8('\n');
  __m128i carriage = _mm_set1_epi8('\r');

  __m256i spaces256 = _mm256_set1_epi8(' ');
  __m256i newline256 = _mm256_set1_epi8('\n');
  __m256i carriage256 = _mm256_set1_epi8('\r');

  size_t i = 0;
  for (; i + 32 - 1 < length; i += 32) {
    __m256i x = _mm256_loadu_si256((const __m256i *)(src + i));
    unsigned int masklow, maskhigh;
    x = cleanm256(x, spaces256, newline256, carriage256, &masklow, &maskhigh);
    int offset1 = 16 - _mm_popcnt_u32(masklow);
    int offset2 = 16 - _mm_popcnt_u32(maskhigh);
    _mm256_storeu2_m128i((__m128i *)(dst + pos + offset1),
                         (__m128i *)(dst + pos), x);
    pos += offset1 + offset2;
  }
  for (; i + 16 - 1 < length; i += 16) {
    __m128i x = _mm_loadu_si128((const __m128i *)(src + i));
    int mask16;
    x = cleanm128(x, spaces, newline, carriage, &mask16);
    _mm_storeu_si128((__m128i *)(dst + pos), x);
    pos += 16 - _mm_popcnt_u32(mask16);
  }
  for (; i < length; i++) {
    char c = src[i];
    if (c == '\r' || c == '\n' || c == ' ') {
      continue;
    }
    dst[pos++] = c;
  }
  return pos;
}

size_t avx2_despace_branchless_u2(const char* src_void, char* dst_void, size_t length)
{
	uint8_t* src = (uint8_t*)src_void;
	uint8_t* dst = (uint8_t*)dst_void;
  size_t pos = 0;
  __m128i spaces = _mm_set1_epi8(' ');
  __m128i newline = _mm_set1_epi8('\n');
  __m128i carriage = _mm_set1_epi8('\r');

  __m256i spaces256 = _mm256_set1_epi8(' ');
  __m256i newline256 = _mm256_set1_epi8('\n');
  __m256i carriage256 = _mm256_set1_epi8('\r');

  size_t i = 0;
  for (; i + 64 - 1 < length; i += 64) {
    __m256i x1, x2;
    int offset11, offset12, offset21, offset22;
    unsigned int masklow1, maskhigh1, masklow2, maskhigh2;

    x1 = _mm256_loadu_si256((const __m256i *)(src + i));
    x2 = _mm256_loadu_si256((const __m256i *)(src + i + 32));

    x1 = cleanm256(x1, spaces256, newline256, carriage256, &masklow1,
                   &maskhigh1);
    offset11 = 16 - _mm_popcnt_u32(masklow1);
    offset12 = 16 - _mm_popcnt_u32(maskhigh1);
    x2 = cleanm256(x2, spaces256, newline256, carriage256, &masklow2,
                   &maskhigh2);
    offset21 = 16 - _mm_popcnt_u32(masklow2);
    offset22 = 16 - _mm_popcnt_u32(maskhigh2);

    _mm256_storeu2_m128i((__m128i *)(dst + pos + offset11),
                         (__m128i *)(dst + pos), x1);
    pos += offset11 + offset12;

    _mm256_storeu2_m128i((__m128i *)(dst + pos + offset21),
                         (__m128i *)(dst + pos), x2);
    pos += offset21 + offset22;
  }

  for (; i + 32 - 1 < length; i += 32) {
    unsigned int masklow, maskhigh;

    int offset1, offset2;
    __m256i x = _mm256_loadu_si256((const __m256i *)(src + i));
    x = cleanm256(x, spaces256, newline256, carriage256, &masklow, &maskhigh);
    offset1 = 16 - _mm_popcnt_u32(masklow);
    offset2 = 16 - _mm_popcnt_u32(maskhigh);
    _mm256_storeu2_m128i((__m128i *)(dst + pos + offset1),
                         (__m128i *)(dst + pos), x);
    pos += offset1 + offset2;
  }
  for (; i + 16 - 1 < length; i += 16) {
    __m128i x = _mm_loadu_si128((const __m128i *)(src + i));
    int mask16;
    x = cleanm128(x, spaces, newline, carriage, &mask16);
    _mm_storeu_si128((__m128i *)(dst + pos), x);
    pos += 16 - _mm_popcnt_u32(mask16);
  }
  for (; i < length; i++) {
    char c = src[i];
    if (c == '\r' || c == '\n' || c == ' ') {
      continue;
    }
    dst[pos++] = c;
  }
  return pos;
}

unsigned __int64 Test(int n, char* inbuf, char* outbuf, size_t len, char* refbuf, size_t* compressed, double refTime = 0.0) {
	unsigned __int64 startTime = 0, endTime = 0, diff = MAXULONG64;
	TEST_PTR t = m[n].timed;
	memset(outbuf, 0, len);
	*compressed = (t)(inbuf, outbuf, len);
#if !defined(_DEBUG)
	for (int retries = 0; retries < m[n].retry; retries++) {
		startTime = __rdtsc();
		(t)(inbuf, outbuf, len);
		endTime = __rdtsc();
		diff = min(diff, endTime - startTime);
	}
#endif

	cout << std::setw(2) << right << n << ' ';
	cout <<"TSC (" << std::setw(31) << left << m[n].name << "):" << dec << std::setw(9) << right << std::setprecision(3) << fixed << diff;
	if (refTime != 0) {
		cout << "  SpeedUp (Scalar/" << std::setw(20) << left << m[n].name << "):";
		cout << std::setw(8)<< right << std::setprecision(3) << fixed << refTime / (double)diff << ' ';
		cout << std::setw(8)<< right << std::setprecision(3) << fixed << (double)diff / ((double)len / (double)m[n].length) << ' ';
		cout << std::setw(8)<< right << std::setprecision(3) << fixed << (double)diff / (double)len << ' ';
		cout << m[n].length << endl;
	} else {
		cout << endl;
	}
	if ((m[n].test) && (refbuf != 0)) {
		int res = memcmp(outbuf, refbuf, *compressed);
		if (res != 0) {
			cout << endl << "Difference!" << endl;
			for (size_t i = 0; i < *compressed; i++)
				if (outbuf[i] != refbuf[i]) {
					cout << "index:" << dec << i << " outbuf:" << hex << outbuf[i] << " refbuf:" << hex << refbuf[i] << endl;
					break;
				}
		}
	}

	return diff;
}

int main(int argc, char* argv[])
{
	if (argc == 3) {
		char* inbuf = 0;
		char* outbuf0 = 0;
		char* outbuf1 = 0;
		unsigned __int64 diff = 0;
		ifstream textfile;
		ofstream textfile_nospace;
		streampos inlength;
		size_t compressed = 0;

		textfile.open(argv[1], ios::in | ios::binary | ios::ate);
		textfile_nospace.open(argv[2], ios::out | ios::binary);
		inlength = textfile.tellg();
		textfile.seekg(0, ios::beg);

		inbuf = new char[inlength];
		outbuf0 = new char[inlength];
		outbuf1 = new char[inlength];

		textfile.read(inbuf, inlength);
		if (cpu_props.IsFeat(ISA_AVX)) {
			_mm256_zeroall();
		}
		cout << "----------------------------------------" << endl;
		diff = Test(0, inbuf, outbuf0, inlength, 0, &compressed, 0.0);
		unsigned __int64 bestDiff = MAXULONG64, diff2 = 0;
		int best = 0;
		for (int b = 1; b < sizeof(m) / sizeof(methods); b++) {
			if (cpu_props.IsFeat(m[b].isa)) {
				diff2 = Test(b, inbuf, outbuf1, inlength, outbuf0, &compressed, (double)diff);
				if (diff2 < bestDiff) {
					bestDiff = diff2;
					best = b;
				}
			}
		}
		cout << "========================================" << endl;
		cout << std::setw(2) << left << best << ' ';
		cout <<"TSC (" << std::setw(31) << right << m[best].name << "):" << dec << std::setw(9) << right << std::setprecision(3) << fixed << bestDiff << endl;

		textfile_nospace.write(outbuf1, compressed);
		textfile_nospace.close();
		textfile.close();
	}
}

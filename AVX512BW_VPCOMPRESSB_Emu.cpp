// AVX512BW_VPCOMPRESSB_Emu
//

#include "stdafx.h"

using namespace std;

#define RETRY			1000

CPU_Props cpu_props;

extern "C" UINT32 __fastcall CheckISA();
extern "C" size_t __fastcall Test_VPCOMPRESSD_Asm(const char* src, char* dst, size_t len);
extern "C" size_t __fastcall Test_VPCOMPRESSB_Asm(const char* src, char* dst, size_t len);

typedef size_t(__fastcall* TEST_PTR)(const char*, char*, size_t);

size_t remove_spaces_scalar(const char* src, char* dst, size_t n);
size_t remove_spaces_avx512vbmi(const char* src, char* dst, size_t n);
size_t remove_spaces_avx512bw(const char* src, char* dst, size_t n);
size_t remove_spaces_avx512vbmi_zach(const char* src, char* dst, size_t n);
size_t despace_avx2_vpermd(const char* src_void, char* dst_void, size_t length);

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
	{"Scalar              ",	"X64        ",	32, 500,	remove_spaces_scalar,				ISA_RDTSC,			true},
	{"AVX2_VPERMD         ",	"AVX2       ",	32, RETRY,	despace_avx2_vpermd,				ISA_AVX2,			true},
	{"AVX512BW_Intrin     ",	"AVX512BW   ",	32, RETRY,	remove_spaces_avx512bw,				ISA_AVX512BW,		true},
	{"AVX512BW_Asm        ",	"AVX512BW   ",	32, RETRY,	Test_VPCOMPRESSD_Asm,				ISA_AVX512BW,		true},
	{"AVX512VBMI_Intrin   ",	"AVX512VBMI ",	32, RETRY,	remove_spaces_avx512vbmi,			ISA_AVX512VBMI,		true},
	{"AVX512VBMI_Zach     ",	"AVX512VBMI ",	32, RETRY,	remove_spaces_avx512vbmi_zach,		ISA_AVX512VBMI,		true},
	{"AVX512VBMI2_Asm     ",	"AVX512VBMI2",	32, RETRY,	Test_VPCOMPRESSB_Asm,				ISA_AVX512_VBMI2,	true},
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

size_t remove_spaces_avx512vbmi(const char* src, char* dst, size_t n) {
	char* startdst = dst;

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
	for (size_t i = 0; i < n; i += 64) {

		const __m512i input = _mm512_loadu_si512((const __m512i*)(src + i));
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
	}

	return dst - startdst;
}

#define VPCOMPRESSD_Core(addr) \
	_mm512_mask_cvtepi32_storeu_epi8(dst, full, _mm512_mask_compress_epi32(zero, __mmask16(mask), _mm512_cvtepu8_epi32(*(const __m128i*)(addr))));	\
	_mm512_mask_cvtepi32_storeu_epi8(dst + _mm_popcnt_u64((WORD)mask), full, _mm512_mask_compress_epi32(zero, __mmask16(mask >> 16), _mm512_cvtepu8_epi32(*(const __m128i*)(addr + 16)))); \
	_mm512_mask_cvtepi32_storeu_epi8(dst + _mm_popcnt_u64((DWORD)mask), full, _mm512_mask_compress_epi32(zero, __mmask16(mask >> 32), _mm512_cvtepu8_epi32(*(const __m128i*)(addr + 32)))); \
	_mm512_mask_cvtepi32_storeu_epi8(dst + _mm_popcnt_u64(mask << 16), full, _mm512_mask_compress_epi32(zero, __mmask16(mask >> 48), _mm512_cvtepu8_epi32(*(const __m128i*)(addr + 48))));

size_t remove_spaces_avx512bw(const char* src, char* dst, size_t n) {
	char* startdst = dst;

	const __m512i NL = _mm512_set1_epi8('\n');
	const __m512i CR = _mm512_set1_epi8('\r');
	const __m512i spaces = _mm512_set1_epi8(' ');
	const __m512i zero = _mm512_setzero_si512();
	const __mmask16 full = 0xffff;
	size_t i = 0;
	uint64_t mask;
	for (size_t i = 0; i < (n & ~0x3f); i += 64) {
		const __m512i input = _mm512_loadu_si512((const __m512i*)(src + i));
		mask = _mm512_cmpneq_epi8_mask(input, spaces)
			& _mm512_cmpneq_epi8_mask(input, NL)
			& _mm512_cmpneq_epi8_mask(input, CR);

		VPCOMPRESSD_Core(src + i)
		dst += _mm_popcnt_u64(mask);
	}
	if ((n & 0x3f) != 0) {
		mask = _bzhi_u64(~0ULL, (int)(n & 0x3f));
		const char * addr = src + (n & ~0x3f);
		const __m512i input = _mm512_maskz_loadu_epi8(mask, (const __m512i*)(addr));
		mask &= _mm512_cmpneq_epi8_mask(input, spaces)
			& _mm512_cmpneq_epi8_mask(input, NL)
			& _mm512_cmpneq_epi8_mask(input, CR);

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

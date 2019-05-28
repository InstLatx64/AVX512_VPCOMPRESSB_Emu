.data
.code

; based on Intel Architecture Instruction Set Extensions Programming Reference
; 319433_015.pdf, p.23
CheckISA proc
; in   - void
; out  - eax[0] - RDTSC
;	   - eax[1] - RDTSCP
;	   - eax[2] - POPCNT
;	   - eax[3] - RDRAND
;	   - eax[4] - AVX2 & OS
;	   - eax[5] - AVX512F & OS
;	   - eax[6] - AVX512BW & AVX512F & OS
;	   - eax[7] - AVX512VL & AVX512F & OS
;	   - eax[8] - AVX512VBMI & AVX512F & OS
;	   - eax[9] - VPOPCNT & AVX512F & OS
;	   - eax[10]- BITALG & AVX512F & OS
;	   - eax[11]- VBMI2 & AVX512F & OS
; used - 
;	   - rax, rcx, rdx, rbx, r8, r9, r10, r11
		push			rbx
		push			rdx
		mov				eax, 1
		cpuid
		mov				r8d, ecx
		mov				r9d, edx
		xor				eax, eax
		bt  			ecx, 27					; check OSXSAVE flag
		jnc				no_osxsave

		mov				eax, 80000001h
		xor				ecx, ecx
		cpuid									; edx[27]: RDTSCP
		mov				r11d, edx

		mov				eax, 7
		xor				ecx, ecx
		cpuid									; ebx[30]: AVX512BW ebx[16]: AVX512F ebx[5]: AVX2
		mov				r10d, ecx				; ecx[14]: VPOPCNTDQ ecx[12]: BITALG

		xor				ecx, ecx				; specify 0 for XFEATURE_ENABLED_MASK register
		xgetbv									; result in EDX:EAX
		; edx, ecx free
		xor				ecx, ecx
		xor				edx, edx
		and				eax, 0e6h
		cmp				eax, 0e6h				; verify that XCR0[7:5] = ‘111b’ (OPMASK state, upper 256-bit of ZMM0-ZMM15 and 
												; ZMM16-ZMM31 state are enabled by OS) and that XCR0[2:1] = ‘11b’ (XMM state and YMM state are enabled by OS)
		sete			dl
		and				eax, 06h				; XCR0[2:1] = ‘11b’ (XMM state and YMM state are enabled by OS)
		cmp				eax, 06h
		sete			cl
		xor				eax, eax
		bt				r10d, 6					; check VBMI2
		rcl				eax, 1
		bt				r10d, 12				; check BITALG
		rcl				eax, 1
		bt				r10d, 14				; check VPOPCNTDQ
		rcl				eax, 1
		bt				r10d, 1					; check AVX512VBMI
		rcl				eax, 1
		bt				ebx, 31					; check AVX512VL
		rcl				eax, 1
		bt				ebx, 30					; check AVX512BW
		rcl				eax, 1

		xor				r10, r10
		bt				ebx, 16					; check AVX512F
		setc			r10b
		rcl				eax, 1

		and				edx, r10d
		imul			eax, edx				; clear w/o OS support

		bt				ebx, 5					; check AVX2 flag
		rcl				eax, 1
		imul			eax, ecx				; clear w/o OS support
no_osxsave:
		bt				r8d, 30					; check RDRAND flag
		rcl				eax, 1
		bt				r8d, 23					; check HW POPCNT flag
		rcl				eax, 1
		bt				r11d, 27				; check RDTSCP flag
		rcl				eax, 1
		bt				r9d, 1					; check RDTSC flag
		rcl				eax, 1

		pop				rdx
		pop				rbx
		ret
CheckISA endp

;idea: @geofflangdale https://branchfree.org/2018/05/22/bits-to-indexes-in-bmi2-and-avx-512/
Test_VPCOMPRESSD_Asm	proc
; in
; rcx  - inbuf
; rdx  - outbuf
; r8   - length, remain
; r9   - popcnt(15-0 bits)
; r10  - popcnt(31-0 bits)
; r11  - popcnt(47-0 bits)
; r12  - popcnt(63-0 bits)
; r13 - const 16 for shrx
; out
; rax - compressed size
	push				r15
	push				r14
	push				r13
	push				r12
	mov					r9d, 0dh
	mov					r10d, 0ah
	mov					r11d, 020h
	mov					r13d, 16
	xor					eax, eax
	vpbroadcastb		zmm31, r9
	vpbroadcastb		zmm30, r10
	vpbroadcastb		zmm29, r11
align 16
mainloop:
	cmp					r8, 40h
	jge					full512
	test				r8, r8
	jl					ready;
	mov					r9, -1
	bzhi				r9, r9, r8
	kmovq				k4, r9
	vmovdqu8			zmm25 {k4}{z}, [rcx]
	vpcmpub				k1 {k4}, zmm31, zmm25, 4
	vpcmpub				k2 {k4}, zmm30, zmm25, 4
	vpcmpub				k3 {k4}, zmm29, zmm25, 4
	jmp					compress
align 16
full512:
	vmovdqu8			zmm25, [rcx]						;P23    ; 1- 8;
	vpcmpub				k1, zmm31, zmm25, 4					;P5     ; 9-11;
	vpcmpub				k2, zmm30, zmm25, 4					;P5     ;10-12;
	vpcmpub				k3, zmm29, zmm25, 4					;P5     ;11-13;
compress:
	vpmovzxbd			zmm25, xmm25						;P5     ;12-12;
	vpmovzxbd			zmm26, xmmword ptr [rcx + 10h]		;P237P5 ;13-13;
	kandq				k1, k1, k2							;P0     ;13-13;
	vpmovzxbd			zmm27, xmmword ptr [rcx + 20h]		;P237P5 ;14-14;
	kandq				k1, k1, k3							;P0     ;14-14;
	vpmovzxbd			zmm28, xmmword ptr [rcx + 30h]		;P237P5 ;15-15;
	kmovq				r12, k1								;P0     ;15-15;
	movzx				r9, r12w							;P0156  ;16-16;15_0 bits
	mov					r10d, r12d							;<>     ;16-16;31_0 bits
	shlx				r11, r12, r13						;P06    ;16-16;47_0 bits
	popcnt				r9, r9								;P1     ;17-19;popcnt(15-0 bits)
	kshiftrq			k2, k1, 16							;P5     ;17-20;
	popcnt				r10, r10							;P1     ;18-20;popcnt(31-0 bits)
	kshiftrq			k3, k1, 32							;P5     ;18-21;
	popcnt				r11, r11							;P1     ;19-21;popcnt(47-0 bits)
	kshiftrq			k4, k1, 48							;P5     ;19-22;
	popcnt				r12, r12							;P1     ;16-18;popcnt(63-0 bits)
	vpcompressd			zmm25 {k1}, zmm25					;P5     ;20-22;
	vpcompressd			zmm26 {k2}, zmm26					;P5     ;22-24;
	vpcompressd			zmm27 {k3}, zmm27					;P5     ;24-26;
	vpcompressd			zmm28 {k4}, zmm28					;P5     ;26-28;
	vpmovdb				xmmword ptr [rdx], zmm25			;P237P5 ;28-
	vpmovdb				xmmword ptr [rdx + r9], zmm26		;P237P5 ;30-
	vpmovdb				xmmword ptr [rdx + r10], zmm27		;P237P5 ;32-
	vpmovdb				xmmword ptr [rdx + r11], zmm28		;P237P5 ;34-
	lea					rcx, [rcx + 40h]					;P15
	lea					rdx, [rdx + r12]					;P15
	add					rax, r12							;P0156
	sub					r8, 40h								;P0156
	jmp					mainloop							;P6
ready:
	pop					r12
	pop					r13
	pop					r14
	pop					r15
	ret
Test_VPCOMPRESSD_Asm endp

; supported from Intel Sunny Cove core
Test_VPCOMPRESSB_Asm	proc
; in
; rcx - inbuf
; rdx - outbuf
; r8  - length, remain
; out
; rax - compressed size
	push				r15
	push				r14
	push				r13
	push				r12
	mov					r9d, 0dh
	mov					r10d, 0ah
	mov					r11d, 020h
	xor					eax, eax
	vpbroadcastb		zmm31, r9
	vpbroadcastb		zmm30, r10
	vpbroadcastb		zmm29, r11
align 16
mainloop:
	cmp					r8, 40h
	jge					full512
	test				r8, r8
	jl					ready;
	mov					r9, -1
	bzhi				r9, r9, r8
	kmovq				k4, r9
	vmovdqu8			zmm25 {k4}{z}, [rcx]
	vpcmpub				k1 {k4}, zmm31, zmm25, 4
	vpcmpub				k2 {k4}, zmm30, zmm25, 4
	vpcmpub				k3 {k4}, zmm29, zmm25, 4
	jmp					compress
align 16
full512:
	vmovdqu8			zmm28, [rcx]
	vpcmpub				k1, zmm31, zmm28, 4
	vpcmpub				k2, zmm30, zmm28, 4
	vpcmpub				k3, zmm29, zmm28, 4
compress:
	kandq				k1, k1, k2
	kandq				k1, k1, k3
	vpcompressb			zmm28 {k1}, zmm28
	kmovq				r15, k1
	popcnt				r15, r15
	vmovdqu8			zmmword ptr [rdx], zmm28
	lea					rcx, [rcx + 40h]
	lea					rdx, [rdx + r15]
	add					rax, r15
	sub					r8, 40h
	jmp					mainloop
ready:
	pop					r12
	pop					r13
	pop					r14
	pop					r15
	ret
Test_VPCOMPRESSB_Asm endp

end
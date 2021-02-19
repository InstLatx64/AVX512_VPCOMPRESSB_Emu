.data

const_lut	db 20h, 00h, 00h, 00h, 00h, 00h, 00h, 00h, 00h, 09h, 0ah, 00h, 00h, 0dh, 00h, 00h

.code

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
	vbroadcasti64x2		zmm31, xmmword ptr const_lut
align 16
mainloop:
	cmp					r8, 40h
	jge					full512
	test				r8, r8
	jle					ready;
	mov					r9, -1
	bzhi				r9, r9, r8
	kmovq				k4, r9
	vmovdqu8			zmm25 {k4}{z}, [rcx]
	vpshufb				zmm30 {k4}, zmm31, zmm25
	vpcmpub				k1 {k4}, zmm30, zmm25, 4
	jmp					compress
align 16
full512:
	vmovdqu8			zmm25, [rcx]						;P23    ; 1- 8;
	vpshufb				zmm30, zmm31, zmm25
	vpcmpub				k1, zmm30, zmm25, 4
compress:
	vpmovzxbd			zmm25, xmm25						;P5     ;12-12;
	vpmovzxbd			zmm26, xmmword ptr [rcx + 10h]		;P237P5 ;13-13;
	vpmovzxbd			zmm27, xmmword ptr [rcx + 20h]		;P237P5 ;14-14;
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
	add					rcx, 40h							;P0156
	sub					r8, 40h								;P0156
	vpmovdb				xmmword ptr [rdx], zmm25			;P237P5 ;28-
	vpmovdb				xmmword ptr [rdx + r9], zmm26		;P237P5 ;30-
	vpmovdb				xmmword ptr [rdx + r10], zmm27		;P237P5 ;32-
	vpmovdb				xmmword ptr [rdx + r11], zmm28		;P237P5 ;34-
	add					rdx, r12							;P0156
	add					rax, r12							;P0156
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
	xor					eax, eax
	kxnorq				k4, k4, k4
	vbroadcasti64x2		zmm31, xmmword ptr const_lut
align 16
mainloop:
	cmp					r8, 40h
	jnge				remain
	vmovdqu8			zmm28, [rcx]
	vpshufb				zmm30, zmm31, zmm28
	vpcmpub				k1{k4}, zmm30, zmm28, 4
compress:
	vpcompressb			zmm28 {k1}, zmm28
	kmovq				r10, k1
	popcnt				r10, r10
	vmovdqu8			zmmword ptr [rdx]{k4}, zmm28
	add					rcx, 40h
	sub					r8, 40h
	add					rax, r10
	add					rdx, r10
	jmp					mainloop
remain:
	test				r8, r8
	jle					ready;
	or					r9, -1
	bzhi				r9, r9, r8
	kmovq				k4, r9
	vmovdqu8			zmm28 {k4}{z}, [rcx]
	vpshufb				zmm30 {k4}, zmm31, zmm28
	vpcmpub				k1 {k4}, zmm30, zmm28, 4
	vpcompressb			zmm28 {k1}, zmm28
	kmovq				r10, k1
	popcnt				r10, r10
	add					rax, r10
	vmovdqu8			zmmword ptr [rdx]{k4}, zmm28
ready:
	ret
Test_VPCOMPRESSB_Asm endp

Test_VPCOMPRESSB2_Asm	proc
; in
; rcx - inbuf
; rdx - outbuf
; r8  - length, remain
; out
; rax - compressed size
	xor					eax, eax
	vbroadcasti64x2		zmm31, xmmword ptr const_lut
align 16
mainloop:
	cmp					r8, 80h
	jnge				remain
	vmovdqu8			zmm28, [rcx]
	vmovdqu8			zmm27, [rcx + 40h]
	vpshufb				zmm30, zmm31, zmm28
	vpshufb				zmm29, zmm31, zmm27
	vpcmpub				k1, zmm30, zmm28, 4
	vpcmpub				k2, zmm29, zmm27, 4
compress:
	vpcompressb			zmm28 {k1}, zmm28
	vpcompressb			zmm27 {k2}, zmm27
	kmovq				r10, k1
	kmovq				r11, k2
	popcnt				r10, r10
	popcnt				r11, r11
	vmovdqu8			zmmword ptr [rdx], zmm28
	vmovdqu8			zmmword ptr [rdx + r10], zmm27
	add					rcx, 80h
	sub					r8, 80h
	add					rax, r10
	add					rdx, r10
	add					rax, r11
	add					rdx, r11
	jmp					mainloop
remain:
	test				r8, r8
	jle					ready;
	mov					r10, -1
	xor					r11, r11
	bzhi				r9, r10, r8
	sub					r8, 40h
	cmovnge				r8, r11
	bzhi				r10, r10, r8
	kmovq				k4, r9
	kmovq				k5, r10
	vmovdqu8			zmm28 {k4}{z}, [rcx]
	vmovdqu8			zmm27 {k5}{z}, [rcx + 40h]
	vpshufb				zmm30, zmm31, zmm28
	vpshufb				zmm29, zmm31, zmm27
	vpcmpub				k1 {k4}, zmm30, zmm28, 4
	vpcmpub				k2 {k5}, zmm29, zmm27, 4
	vpcompressb			zmm28 {k1}{z}, zmm28
	vpcompressb			zmm27 {k2}{z}, zmm27
	kmovq				r10, k1
	kmovq				r11, k2
	popcnt				r10, r10
	popcnt				r11, r11
	add					rax, r10
	add					rax, r11
	vmovdqu8			zmmword ptr [rdx]{k4}, zmm28
	vmovdqu8			zmmword ptr [rdx + r10]{k5}, zmm27
ready:
	ret
Test_VPCOMPRESSB2_Asm endp

Test_VPCOMPRESSB4_Asm	proc
; in
; rcx - inbuf
; rdx - outbuf
; r8  - length, remain
; out
; rax - compressed size
	push				rbx
	xor					eax, eax
	vbroadcasti64x2		zmm31, xmmword ptr const_lut
align 16
mainloop:
	cmp					r8, 100h
	jnge				remain
	vmovdqu8			zmm28, [rcx]
	vmovdqu8			zmm27, [rcx + 040h]
	vmovdqu8			zmm26, [rcx + 080h]
	vmovdqu8			zmm25, [rcx + 0C0h]
	vpshufb				zmm30, zmm31, zmm28
	vpshufb				zmm29, zmm31, zmm27
	vpshufb				zmm24, zmm31, zmm26
	vpshufb				zmm23, zmm31, zmm25
	vpcmpub				k1, zmm30, zmm28, 4
	vpcmpub				k2, zmm29, zmm27, 4
	vpcmpub				k3, zmm24, zmm26, 4
	vpcmpub				k4, zmm23, zmm25, 4
compress:
	vpcompressb			zmm28 {k1}{z}, zmm28
	vpcompressb			zmm27 {k2}{z}, zmm27
	vpcompressb			zmm26 {k3}{z}, zmm26
	vpcompressb			zmm25 {k4}{z}, zmm25
	kmovq				r9, k1
	kmovq				r10, k2
	kmovq				r11, k3
	kmovq				rbx, k4
	popcnt				r9, r9
	popcnt				r10, r10
	popcnt				r11, r11
	popcnt				rbx, rbx

	add					r10, r9
	add					r11, r10
	add					rbx, r11

	vmovdqu8			zmmword ptr [rdx], zmm28
	vmovdqu8			zmmword ptr [rdx + r9], zmm27
	vmovdqu8			zmmword ptr [rdx + r10], zmm26
	vmovdqu8			zmmword ptr [rdx + r11], zmm25

	add					rcx, 100h
	sub					r8, 100h

	add					rax, rbx
	add					rdx, rbx

	jmp					mainloop
remain:
	test				r8, r8
	jle					ready;
	mov					r10, -1
	xor					r11, r11
	bzhi				r9, r10, r8
	kmovq				k4, r9
	vmovdqu8			zmm28 {k4}{z}, [rcx]
	sub					r8, 40h
	cmovnge				r8, r11
	bzhi				r9, r10, r8
	kmovq				k5, r9
	vmovdqu8			zmm27 {k5}{z}, [rcx + 040h]
	sub					r8, 40h
	cmovnge				r8, r11
	bzhi				r9, r10, r8
	kmovq				k6, r9
	vmovdqu8			zmm26 {k6}{z}, [rcx + 080h]
	sub					r8, 40h
	cmovnge				r8, r11
	bzhi				r9, r10, r8
	kmovq				k7, r9
	vmovdqu8			zmm25 {k7}{z}, [rcx + 0C0h]

	vpshufb				zmm30, zmm31, zmm28
	vpshufb				zmm29, zmm31, zmm27
	vpshufb				zmm24, zmm31, zmm26
	vpshufb				zmm23, zmm31, zmm25
	vpcmpub				k1 {k4}, zmm30, zmm28, 4
	vpcmpub				k2 {k5}, zmm29, zmm27, 4
	vpcmpub				k3 {k6}, zmm24, zmm26, 4
	vpcmpub				k0 {k7}, zmm23, zmm25, 4
	vpcompressb			zmm28 {k1}{z}, zmm28
	kmovq				r9, k1
	kmovq				k1, k0
	vpcompressb			zmm27 {k2}{z}, zmm27
	vpcompressb			zmm26 {k3}{z}, zmm26
	vpcompressb			zmm25 {k1}{z}, zmm25
	kmovq				r10, k2
	kmovq				r11, k3
	kmovq				rbx, k0

	popcnt				r9, r9
	popcnt				r10, r10
	popcnt				r11, r11
	popcnt				rbx, rbx

	add					r10, r9
	add					r11, r10
	add					rbx, r11
	add					rax, rbx
	
	vmovdqu8			zmmword ptr [rdx]{k4}, zmm28
	vmovdqu8			zmmword ptr [rdx + r9]{k5}, zmm27
	vmovdqu8			zmmword ptr [rdx + r10]{k6}, zmm26
	vmovdqu8			zmmword ptr [rdx + r11]{k7}, zmm25

ready:
	pop					rbx
	ret
Test_VPCOMPRESSB4_Asm endp

Test_VPCOMPRESSB_ymm_Asm	proc
; in
; rcx - inbuf
; rdx - outbuf
; r8  - length, remain
; out
; rax - compressed size
	xor					eax, eax
	kxnorq				k4, k4, k4
	vbroadcasti128		ymm2, xmmword ptr const_lut
align 16
mainloop:
	cmp					r8, 20h
	jnge				remain
	vmovdqu8			ymm0, [rcx]
	vpshufb				ymm1, ymm2, ymm0
	vpcmpub				k1, ymm1, ymm0, 4
	vpcompressb			ymm0 {k1}, ymm0
	kmovd				r10d, k1
	popcnt				r10d, r10d
	vmovdqu8			ymmword ptr [rdx], ymm0
	add					rcx, 20h
	sub					r8, 20h
	add					rax, r10
	add					rdx, r10
	jmp					mainloop
remain:
	test				r8, r8
	jle					ready;
	mov					r9d, -1
	bzhi				r9d, r9d, r8d
	kmovd				k4, r9d
	vmovdqu8			ymm0 {k4}{z}, [rcx]
	vpshufb				ymm1 {k4}, ymm2, ymm0
	vpcmpub				k1 {k4}, ymm1, ymm0, 4
	vpcompressb			ymm0 {k1}, ymm0
	kmovd				r10d, k1
	popcnt				r10d, r10d
	add					rax, r10
	vmovdqu8			ymmword ptr [rdx]{k4}, ymm0
ready:
	ret
Test_VPCOMPRESSB_ymm_Asm endp

Test_VPCOMPRESSB2_ymm_Asm	proc
; in
; rcx - inbuf
; rdx - outbuf
; r8  - length, remain
; out
; rax - compressed size
	xor					eax, eax
	vbroadcasti128		ymm2, xmmword ptr const_lut
align 16
mainloop:
	cmp					r8, 40h
	jnge				remain
	vmovdqu8			ymm0, [rcx]
	vmovdqu8			ymm3, [rcx + 20h]
	vpshufb				ymm1, ymm2, ymm0
	vpshufb				ymm4, ymm2, ymm3
	vpcmpub				k1, ymm1, ymm0, 4
	vpcmpub				k2, ymm4, ymm3, 4
	vpcompressb			ymm0 {k1}, ymm0
	vpcompressb			ymm3 {k2}, ymm3
	kmovd				r10d, k1
	kmovd				r11d, k2
	popcnt				r10d, r10d
	popcnt				r11d, r11d
	vmovdqu8			ymmword ptr [rdx], ymm0
	vmovdqu8			ymmword ptr [rdx + r10], ymm3
	add					rcx, 40h
	sub					r8, 40h
	add					rax, r10
	add					rdx, r10
	add					rax, r11
	add					rdx, r11
	jmp					mainloop
remain:
	test				r8, r8
	jle					ready;
	mov					r10d, -1
	xor					r11d, r11d
	bzhi				r9d, r10d, r8d
	sub					r8, 20h
	cmovnge				r8d, r11d
	bzhi				r10d, r10d, r8d
	kmovq				k4, r9
	kmovq				k5, r10
	vmovdqu8			ymm0 {k4}{z}, [rcx]
	vmovdqu8			ymm3 {k5}{z}, [rcx + 20h]
	vpshufb				ymm1 {k4}, ymm2, ymm0
	vpshufb				ymm4 {k5}, ymm2, ymm3
	vpcmpub				k1 {k4}, ymm1, ymm0, 4
	vpcmpub				k2 {k5}, ymm4, ymm3, 4
	vpcompressb			ymm0 {k1}, ymm0
	vpcompressb			ymm3 {k2}, ymm3
	kmovd				r10d, k1
	kmovd				r11d, k2
	popcnt				r10d, r10d
	popcnt				r11d, r11d
	vmovdqu8			ymmword ptr [rdx]{k4}, ymm0
	vmovdqu8			ymmword ptr [rdx + r10]{k5}, ymm3
	add					rax, r10
	jmp					ready
ready:
	ret
Test_VPCOMPRESSB2_ymm_Asm endp

end
/* profanity.cl
 * ============
 * Contains multi-precision arithmetic functions and iterative elliptical point
 * addition which is the heart of profanity.
 *
 * Terminology
 * ===========
 * 
 *
 * Cutting corners
 * ===============
 * In some instances this code will produce the incorrect results. The elliptical
 * point addition does for example not properly handle the case of two points
 * sharing the same X-coordinate. The reason the code doesn't handle it properly
 * is because it is very unlikely to ever occur and the performance penalty for
 * doing it right is too severe. In the future I'll introduce a periodic check
 * after N amount of cycles that verifies the integrity of all the points to
 * make sure that even very unlikely event are at some point rectified.
 * 
 * Currently, if any of the points in the kernels experiences the unlikely event
 * of an error then that point is forever garbage and your runtime-performance
 * will in practice be (i*I-N) / (i*I). i and I here refers to the values given
 * to the program via the -i and -I switches (default values of 255 and 16384
 * respectively) and N is the number of errornous points.
 *
 * So if a single error occurs you'll lose 1/(i*I) of your performance. That's
 * around 0.00002%. The program will still report the same hashrate of course,
 * only that some of that work is entirely wasted on this errornous point.
 *
 * Initialization of main structure
 * ================================
 *
 * Iteration
 * =========
 *
 *
 * TODO
 * ====
 *   * Update comments to reflect new optimizations and structure
 *
 */

/* ------------------------------------------------------------------------ */
/* Multiprecision functions                                                 */
/* ------------------------------------------------------------------------ */
#define MP_WORDS 8
#define MP_BITS 32
#define bswap32(n) (rotate(n & 0x00FF00FF, 24U)|(rotate(n, 8U) & 0x00FF00FF))

typedef uint mp_word;
typedef struct {
	mp_word d[MP_WORDS];
} mp_number;

// mod              = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
__constant const mp_number mod              = { {0xfffffc2f, 0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff} };

// tripleNegativeGx = 0x92c4cc831269ccfaff1ed83e946adeeaf82c096e76958573f2287becbb17b196
__constant const mp_number tripleNegativeGx = { {0xbb17b196, 0xf2287bec, 0x76958573, 0xf82c096e, 0x946adeea, 0xff1ed83e, 0x1269ccfa, 0x92c4cc83 } };

// doubleNegativeGy = 0x6f8a4b11b2b8773544b60807e3ddeeae05d0976eb2f557ccc7705edf09de52bf
__constant const mp_number doubleNegativeGy = { {0x09de52bf, 0xc7705edf, 0xb2f557cc, 0x05d0976e, 0xe3ddeeae, 0x44b60807, 0xb2b87735, 0x6f8a4b11} };

// negativeGy       = 0xb7c52588d95c3b9aa25b0403f1eef75702e84bb7597aabe663b82f6f04ef2777
__constant const mp_number negativeGy       = { {0x04ef2777, 0x63b82f6f, 0x597aabe6, 0x02e84bb7, 0xf1eef757, 0xa25b0403, 0xd95c3b9a, 0xb7c52588 } };


// Multiprecision subtraction. Underflow signalled via return value.
mp_word mp_sub(mp_number * const r, const mp_number * const a, const mp_number * const b) {
	mp_word t, c = 0;

	for (mp_word i = 0; i < MP_WORDS; ++i) {
		t = a->d[i] - b->d[i] - c;
		c = t > a->d[i] ? 1 : (t == a->d[i] ? c : 0);

		r->d[i] = t;
	}

	return c;
}

// Multiprecision subtraction of the modulus saved in mod. Underflow signalled via return value.
mp_word mp_sub_mod(mp_number * const r) {
	mp_number mod = { {0xfffffc2f, 0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff} };

	mp_word t, c = 0;

	for (mp_word i = 0; i < MP_WORDS; ++i) {
		t = r->d[i] - mod.d[i] - c;
		c = t > r->d[i] ? 1 : (t == r->d[i] ? c : 0);

		r->d[i] = t;
	}

	return c;
}

// Multiprecision subtraction modulo M, M = mod.
// This function is often also used for additions by subtracting a negative number. I've chosen
// to do this because:
//   1. It's easier to re-use an already existing function
//   2. A modular addition would have more overhead since it has to determine if the result of
//      the addition (r) is in the gap M <= r < 2^256. This overhead doesn't exist in a
//      subtraction. We immediately know at the end of a subtraction if we had underflow
//      or not by inspecting the carry value. M refers to the modulus saved in variable mod.
void mp_mod_sub(mp_number * const r, const mp_number * const a, const mp_number * const b) {
	mp_word i, t, c = 0;

	for (i = 0; i < MP_WORDS; ++i) {
		t = a->d[i] - b->d[i] - c;
		c = t < a->d[i] ? 0 : (t == a->d[i] ? c : 1);

		r->d[i] = t;
	}

	if (c) {
		c = 0;
		for (i = 0; i < MP_WORDS; ++i) {
			r->d[i] += mod.d[i] + c;
			c = r->d[i] < mod.d[i] ? 1 : (r->d[i] == mod.d[i] ? c : 0);
		}
	}
}

// Multiprecision subtraction modulo M from a constant number.
// I made this in the belief that using constant address space instead of private address space for any
// constant numbers would lead to increase in performance. Judges are still out on this one.
void mp_mod_sub_const(mp_number * const r, __constant const mp_number * const a, const mp_number * const b) {
	mp_word i, t, c = 0;

	for (i = 0; i < MP_WORDS; ++i) {
		t = a->d[i] - b->d[i] - c;
		c = t < a->d[i] ? 0 : (t == a->d[i] ? c : 1);

		r->d[i] = t;
	}

	if (c) {
		c = 0;
		for (i = 0; i < MP_WORDS; ++i) {
			r->d[i] += mod.d[i] + c;
			c = r->d[i] < mod.d[i] ? 1 : (r->d[i] == mod.d[i] ? c : 0);
		}
	}
}

// Multiprecision subtraction modulo M of G_x from a number.
// Specialization of mp_mod_sub in hope of performance gain.
void mp_mod_sub_gx(mp_number * const r, const mp_number * const a) {
	mp_word i, t, c = 0;

	t = a->d[0] - 0x16f81798; c = t < a->d[0] ? 0 : (t == a->d[0] ? c : 1); r->d[0] = t;
	t = a->d[1] - 0x59f2815b - c; c = t < a->d[1] ? 0 : (t == a->d[1] ? c : 1); r->d[1] = t;
	t = a->d[2] - 0x2dce28d9 - c; c = t < a->d[2] ? 0 : (t == a->d[2] ? c : 1); r->d[2] = t;
	t = a->d[3] - 0x029bfcdb - c; c = t < a->d[3] ? 0 : (t == a->d[3] ? c : 1); r->d[3] = t;
	t = a->d[4] - 0xce870b07 - c; c = t < a->d[4] ? 0 : (t == a->d[4] ? c : 1); r->d[4] = t;
	t = a->d[5] - 0x55a06295 - c; c = t < a->d[5] ? 0 : (t == a->d[5] ? c : 1); r->d[5] = t;
	t = a->d[6] - 0xf9dcbbac - c; c = t < a->d[6] ? 0 : (t == a->d[6] ? c : 1); r->d[6] = t;
	t = a->d[7] - 0x79be667e - c; c = t < a->d[7] ? 0 : (t == a->d[7] ? c : 1); r->d[7] = t;

	if (c) {
		c = 0;
		for (i = 0; i < MP_WORDS; ++i) {
			r->d[i] += mod.d[i] + c;
			c = r->d[i] < mod.d[i] ? 1 : (r->d[i] == mod.d[i] ? c : 0);
		}
	}
}

// Multiprecision subtraction modulo M of G_y from a number.
// Specialization of mp_mod_sub in hope of performance gain.
void mp_mod_sub_gy(mp_number * const r, const mp_number * const a) {
	mp_word i, t, c = 0;

	t = a->d[0] - 0xfb10d4b8; c = t < a->d[0] ? 0 : (t == a->d[0] ? c : 1); r->d[0] = t;
	t = a->d[1] - 0x9c47d08f - c; c = t < a->d[1] ? 0 : (t == a->d[1] ? c : 1); r->d[1] = t;
	t = a->d[2] - 0xa6855419 - c; c = t < a->d[2] ? 0 : (t == a->d[2] ? c : 1); r->d[2] = t;
	t = a->d[3] - 0xfd17b448 - c; c = t < a->d[3] ? 0 : (t == a->d[3] ? c : 1); r->d[3] = t;
	t = a->d[4] - 0x0e1108a8 - c; c = t < a->d[4] ? 0 : (t == a->d[4] ? c : 1); r->d[4] = t;
	t = a->d[5] - 0x5da4fbfc - c; c = t < a->d[5] ? 0 : (t == a->d[5] ? c : 1); r->d[5] = t;
	t = a->d[6] - 0x26a3c465 - c; c = t < a->d[6] ? 0 : (t == a->d[6] ? c : 1); r->d[6] = t;
	t = a->d[7] - 0x483ada77 - c; c = t < a->d[7] ? 0 : (t == a->d[7] ? c : 1); r->d[7] = t;

	if (c) {
		c = 0;
		for (i = 0; i < MP_WORDS; ++i) {
			r->d[i] += mod.d[i] + c;
			c = r->d[i] < mod.d[i] ? 1 : (r->d[i] == mod.d[i] ? c : 0);
		}
	}
}

// Multiprecision addition. Overflow signalled via return value.
mp_word mp_add(mp_number * const r, const mp_number * const a) {
	mp_word c = 0;

	for (mp_word i = 0; i < MP_WORDS; ++i) {
		r->d[i] += a->d[i] + c;
		c = r->d[i] < a->d[i] ? 1 : (r->d[i] == a->d[i] ? c : 0);
	}

	return c;
}

// Multiprecision addition of the modulus saved in mod. Overflow signalled via return value.
mp_word mp_add_mod(mp_number * const r) {
	mp_word c = 0;

	for (mp_word i = 0; i < MP_WORDS; ++i) {
		r->d[i] += mod.d[i] + c;
		c = r->d[i] < mod.d[i] ? 1 : (r->d[i] == mod.d[i] ? c : 0);
	}

	return c;
}

// Multiprecision addition of two numbers with one extra word each. Overflow signalled via return value.
mp_word mp_add_more(mp_number * const r, mp_word * const extraR, const mp_number * const a, const mp_word * const extraA) {
	const mp_word c = mp_add(r, a);
	*extraR += *extraA + c;
	return *extraR < *extraA ? 1 : (*extraR == *extraA ? c : 0);
}

// Multiprecision greater than or equal (>=) operator
mp_word mp_gte(const mp_number * const a, const mp_number * const b) {
	mp_word l = 0, g = 0;

	for (mp_word i = 0; i < MP_WORDS; ++i) {
		if (a->d[i] < b->d[i]) l |= (1 << i);
		if (a->d[i] > b->d[i]) g |= (1 << i);
	}

	return g >= l;
}

// Bit shifts a number with an extra word to the right one step
void mp_shr_extra(mp_number * const r, mp_word * const e) {
	r->d[0] = (r->d[1] << 31) | (r->d[0] >> 1);
	r->d[1] = (r->d[2] << 31) | (r->d[1] >> 1);
	r->d[2] = (r->d[3] << 31) | (r->d[2] >> 1);
	r->d[3] = (r->d[4] << 31) | (r->d[3] >> 1);
	r->d[4] = (r->d[5] << 31) | (r->d[4] >> 1);
	r->d[5] = (r->d[6] << 31) | (r->d[5] >> 1);
	r->d[6] = (r->d[7] << 31) | (r->d[6] >> 1);
	r->d[7] = (*e << 31) | (r->d[7] >> 1);
	*e >>= 1;
}

// Bit shifts a number to the right one step
void mp_shr(mp_number * const r) {
	r->d[0] = (r->d[1] << 31) | (r->d[0] >> 1);
	r->d[1] = (r->d[2] << 31) | (r->d[1] >> 1);
	r->d[2] = (r->d[3] << 31) | (r->d[2] >> 1);
	r->d[3] = (r->d[4] << 31) | (r->d[3] >> 1);
	r->d[4] = (r->d[5] << 31) | (r->d[4] >> 1);
	r->d[5] = (r->d[6] << 31) | (r->d[5] >> 1);
	r->d[6] = (r->d[7] << 31) | (r->d[6] >> 1);
	r->d[7] >>= 1;
}

// Multiplies a number with a word and adds it to an existing number with an extra word, overflow of the extra word is signalled in return value
// This is a special function only used for modular multiplication
mp_word mp_mul_word_add_extra(mp_number * const r, const mp_number * const a, const mp_word w, mp_word * const extra) {
	mp_word cM = 0; // Carry for multiplication
	mp_word cA = 0; // Carry for addition
	mp_word tM = 0; // Temporary storage for multiplication

	for (mp_word i = 0; i < MP_WORDS; ++i) {
		tM = (a->d[i] * w + cM);
		cM = mul_hi(a->d[i], w) + (tM < cM);

		r->d[i] += tM + cA;
		cA = r->d[i] < tM ? 1 : (r->d[i] == tM ? cA : 0);
	}

	*extra += cM + cA;
	return *extra < cM ? 1 : (*extra == cM ? cA : 0);
}

// Multiplies a number with a word, potentially adds modhigher to it, and then subtracts it from en existing number, no extra words, no overflow
// This is a special function only used for modular multiplication
void mp_mul_mod_word_sub(mp_number * const r, const mp_word w, const bool withModHigher) {
	// Having these numbers declared here instead of using the global values in __constant address space seems to lead
	// to better optimizations by the compiler on my GTX 1070.
	mp_number mod = { { 0xfffffc2f, 0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff} };
	mp_number modhigher = { {0x00000000, 0xfffffc2f, 0xfffffffe, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff} };

	mp_word cM = 0; // Carry for multiplication
	mp_word cS = 0; // Carry for subtraction
	mp_word tS = 0; // Temporary storage for subtraction
	mp_word tM = 0; // Temporary storage for multiplication
	mp_word cA = 0; // Carry for addition of modhigher

	for (mp_word i = 0; i < MP_WORDS; ++i) {
		tM = (mod.d[i] * w + cM);
		cM = mul_hi(mod.d[i], w) + (tM < cM);

		tM += (withModHigher ? modhigher.d[i] : 0) + cA;
		cA = tM < (withModHigher ? modhigher.d[i] : 0) ? 1 : (tM == (withModHigher ? modhigher.d[i] : 0) ? cA : 0);

		tS = r->d[i] - tM - cS;
		cS = tS > r->d[i] ? 1 : (tS == r->d[i] ? cS : 0);

		r->d[i] = tS;
	}
}

// Modular multiplication. Based on Algorithm 3 (and a series of hunches) from this article:
// https://www.esat.kuleuven.be/cosic/publications/article-1191.pdf
// When I first implemented it I never encountered a situation where the additional end steps
// of adding or subtracting the modulo was necessary. Maybe it's not for the particular modulo
// used in secp256k1, maybe the overflow bit can be skipped in to avoid 8 subtractions and
// trade it for the final steps? Maybe the final steps are necessary but seldom needed?
// I have no idea, for the time being I'll leave it like this, also see the comments at the
// beginning of this document under the title "Cutting corners".
void mp_mod_mul(mp_number * const r, const mp_number * const X, const mp_number * const Y) {
	mp_number Z = { {0} };
	mp_word extraWord;

	for (int i = MP_WORDS - 1; i >= 0; --i) {
		// Z = Z * 2^32
		extraWord = Z.d[7]; Z.d[7] = Z.d[6]; Z.d[6] = Z.d[5]; Z.d[5] = Z.d[4]; Z.d[4] = Z.d[3]; Z.d[3] = Z.d[2]; Z.d[2] = Z.d[1]; Z.d[1] = Z.d[0]; Z.d[0] = 0;

		// Z = Z + X * Y_i
		bool overflow = mp_mul_word_add_extra(&Z, X, Y->d[i], &extraWord);

		// Z = Z - qM
		mp_mul_mod_word_sub(&Z, extraWord, overflow);
	}

	*r = Z;
}

// Modular inversion of a number. 
void mp_mod_inverse(mp_number * const r) {
	mp_number A = { { 1 } };
	mp_number C = { { 0 } };
	mp_number v = mod;

	mp_word extraA = 0;
	mp_word extraC = 0;

	while (r->d[0] || r->d[1] || r->d[2] || r->d[3] || r->d[4] || r->d[5] || r->d[6] || r->d[7]) {
		while (!(r->d[0] & 1)) {
			mp_shr(r);
			if (A.d[0] & 1) {
				extraA += mp_add_mod(&A);
			}

			mp_shr_extra(&A, &extraA);
		}

		while (!(v.d[0] & 1)) {
			mp_shr(&v);
			if (C.d[0] & 1) {
				extraC += mp_add_mod(&C);
			}

			mp_shr_extra(&C, &extraC);
		}

		if (mp_gte(r, &v)) {
			mp_sub(r, r, &v);
			mp_add_more(&A, &extraA, &C, &extraC);
		}
		else {
			mp_sub(&v, &v, r);
			mp_add_more(&C, &extraC, &A, &extraA);
		}
	}

	while (extraC) {
		extraC -= mp_sub_mod(&C);
	}

	v = mod;
	mp_sub(r, &v, &C);
}

/* ------------------------------------------------------------------------ */
/* Elliptic point and addition (with caveats).                              */
/* ------------------------------------------------------------------------ */
typedef struct {
	mp_number x;
	mp_number y;
} point;

// Elliptical point addition
// Does not handle points sharing X coordinate, this is a deliberate design choice.
// For more information on this choice see the beginning of this file.
void point_add(point * const r, point * const p, point * const o) {
	mp_number tmp;
	mp_number newX;
	mp_number newY;

	mp_mod_sub(&tmp, &o->x, &p->x);

	mp_mod_inverse(&tmp);

	mp_mod_sub(&newX, &o->y, &p->y);
	mp_mod_mul(&tmp, &tmp, &newX);

	mp_mod_mul(&newX, &tmp, &tmp);
	mp_mod_sub(&newX, &newX, &p->x);
	mp_mod_sub(&newX, &newX, &o->x);

	mp_mod_sub(&newY, &p->x, &newX);
	mp_mod_mul(&newY, &newY, &tmp);
	mp_mod_sub(&newY, &newY, &p->y);

	r->x = newX;
	r->y = newY;
}

/* ------------------------------------------------------------------------ */
/* Profanity.                                                               */
/* ------------------------------------------------------------------------ */
typedef struct {
	uint found;
	uint foundId;
	uchar foundHash[20];
} result;

void profanity_init_seed(__global const point * const precomp, point * const p, bool * const pIsFirst, const size_t precompOffset, const ulong seed) {
	point o;

	for (uchar i = 0; i < 8; ++i) {
		const uchar shift = i * 8;
		const uchar byte = (seed >> shift) & 0xFF;

		if (byte) {
			o = precomp[precompOffset + i * 255 + byte - 1];
			if (*pIsFirst) {
				*p = o;
				*pIsFirst = false;
			}
			else {
				point_add(p, p, &o);
			}
		}
	}
}

__kernel void profanity_init(__global const point * const precomp, __global mp_number * const pDeltaX, __global mp_number * const pPrevLambda, __global result * const pResult, const ulong4 seed, const ulong4 seedX, const ulong4 seedY) {
	const size_t id = get_global_id(0);
	point p = {
		.x = {.d = {
			seedX.x & 0xFFFFFFFF, seedX.x >> 32,
			seedX.y & 0xFFFFFFFF, seedX.y >> 32,
			seedX.z & 0xFFFFFFFF, seedX.z >> 32,
			seedX.w & 0xFFFFFFFF, seedX.w >> 32,
		}},
		.y = {.d = {
			seedY.x & 0xFFFFFFFF, seedY.x >> 32,
			seedY.y & 0xFFFFFFFF, seedY.y >> 32,
			seedY.z & 0xFFFFFFFF, seedY.z >> 32,
			seedY.w & 0xFFFFFFFF, seedY.w >> 32,
		}},
	};
	point p_random;
	bool bIsFirst = true;

	mp_number tmp1, tmp2;
	point tmp3;

	// Calculate k*G where k = seed.wzyx (in other words, find the point indicated by the private key represented in seed)
	profanity_init_seed(precomp, &p_random, &bIsFirst, 8 * 255 * 0, seed.x);
	profanity_init_seed(precomp, &p_random, &bIsFirst, 8 * 255 * 1, seed.y);
	profanity_init_seed(precomp, &p_random, &bIsFirst, 8 * 255 * 2, seed.z);
	profanity_init_seed(precomp, &p_random, &bIsFirst, 8 * 255 * 3, seed.w + id);
	point_add(&p, &p, &p_random);

	// Calculate current lambda in this point
	mp_mod_sub_gx(&tmp1, &p.x);
	mp_mod_inverse(&tmp1);

	mp_mod_sub_gy(&tmp2, &p.y); 
	mp_mod_mul(&tmp1, &tmp1, &tmp2);

	// Jump to next point (precomp[0] is the generator point G)
	tmp3 = precomp[0];
	point_add(&p, &tmp3, &p);

	// pDeltaX should contain the delta (x - G_x)
	mp_mod_sub_gx(&p.x, &p.x);

	pDeltaX[id] = p.x;
	pPrevLambda[id] = tmp1;

	for (uchar i = 0; i < PROFANITY_MAX_SCORE + 1; ++i) {
		pResult[i].found = 0;
	}
}

// This kernel calculates several modular inversions at once with just one inverse.
// It's an implementation of Algorithm 2.11 from Modern Computer Arithmetic:
// https://members.loria.fr/PZimmermann/mca/pub226.html 
//
// My RX 480 is very sensitive to changes in the second loop and sometimes I have
// to make seemingly non-functional changes to the code to make the compiler
// generate the most optimized version.
__kernel void profanity_inverse(__global const mp_number * const pDeltaX, __global mp_number * const pInverse) {
	const size_t id = get_global_id(0) * PROFANITY_INVERSE_SIZE;

	// negativeDoubleGy = 0x6f8a4b11b2b8773544b60807e3ddeeae05d0976eb2f557ccc7705edf09de52bf
	mp_number negativeDoubleGy = { {0x09de52bf, 0xc7705edf, 0xb2f557cc, 0x05d0976e, 0xe3ddeeae, 0x44b60807, 0xb2b87735, 0x6f8a4b11 } };

	mp_number copy1, copy2;
	mp_number buffer[PROFANITY_INVERSE_SIZE];
	mp_number buffer2[PROFANITY_INVERSE_SIZE];

	// We initialize buffer and buffer2 such that:
	// buffer[i] = pDeltaX[id] * pDeltaX[id + 1] * pDeltaX[id + 2] * ... * pDeltaX[id + i]
	// buffer2[i] = pDeltaX[id + i]
	buffer[0] = pDeltaX[id];
	for (uint i = 1; i < PROFANITY_INVERSE_SIZE; ++i) {
		buffer2[i] = pDeltaX[id + i];
		mp_mod_mul(&buffer[i], &buffer2[i], &buffer[i - 1]);
	}

	// Take the inverse of all x-values combined
	copy1 = buffer[PROFANITY_INVERSE_SIZE - 1];
	mp_mod_inverse(&copy1);

	// We multiply in -2G_y together with the inverse so that we have:
	//            - 2 * G_y
	//  ----------------------------
	//  x_0 * x_1 * x_2 * x_3 * ...
	mp_mod_mul(&copy1, &copy1, &negativeDoubleGy);

	// Multiply out each individual inverse using the buffers
	for (uint i = PROFANITY_INVERSE_SIZE - 1; i > 0; --i) {
		mp_mod_mul(&copy2, &copy1, &buffer[i - 1]);
		mp_mod_mul(&copy1, &copy1, &buffer2[i]);
		pInverse[id + i] = copy2;
	}

	pInverse[id] = copy1;
}

// This kernel performs en elliptical curve point addition. See:
// https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication#Point_addition
// I've made one mathematical optimization by never calculating x_r,
// instead I directly calculate the delta (x_q - x_p). It's for this
// delta we calculate the inverse and that's already been done at this
// point. By calculating and storing the next delta we don't have to
// calculate the delta in profanity_inverse_multiple which saves us
// one call to mp_mod_sub per point, but inversely we have to introduce
// an addition (or addition by subtracting a negative number) in
// profanity_end to retrieve the actual x-coordinate instead of the
// delta as that's what used for calculating the public hash.
//
// One optimization is when calculating the next y-coordinate. As
// given in the wiki the next y-coordinate is given by:
//   y_r = λ²(x_p - x_r) - y_p
// In our case the other point P is the generator point so x_p = G_x,
// a constant value. x_r is the new point which we never calculate, we
// calculate the new delta (x_q - x_p) instead. Let's denote the delta
// with d and new delta as d' and remove notation for points P and Q and
// instead refeer to x_p as G_x, y_p as G_y and x_q as x, y_q as y.
// Furthermore let's denote new x by x' and new y with y'.
//
// Then we have:
//   d = x - G_x <=> x = d + G_x
//   x' = λ² - G_x - x <=> x_r = λ² - G_x - d - G_x = λ² - 2G_x - d
//   
//   d' = x' - G_x = λ² - 2G_x - d - G_x = λ² - 3G_x - d
//
// So we see that the new delta d' can be calculated with the same
// amount of steps as the new x'; 3G_x is still just a single constant.
//
// Now for the next y-coordinate in the new notation:
//   y' =  λ(G_x - x') - G_y
//
// If we expand the expression (G_x - x') we can see that this
// subtraction can be removed! Saving us one call to mp_mod_sub!
//   G_x - x' = -(x' - G_x) = -d'
// It has the same value as the new delta but negated! We can avoid
// having to perform the negation by:
//   y' = λ * -d' - G_y = -G_y - (λ * d')
//
// We can just precalculate the constant -G_y and we get rid of one
// subtraction. Woo!
//
// But we aren't done yet! Let's expand the expression for the next
// lambda, λ'. We have:
//   λ' = (y' - G_y) / d'
//      = (-λ * d' - G_y - G_y) / d' 
//      = (-λ * d' - 2*G_y) / d' 
//      = -λ - 2*G_y / d' 
//
// So the next lambda value can be calculated from the old one. This in
// and of itself is not so interesting but the fact that the term -2 * G_y
// is a constant is! Since it's constant it'll be the same value no matter
// which point we're currently working with. This means that this factor
// can be multiplied in during the inversion, and just with one call per
// inversion instead of one call per point! This is small enough to be
// negligible and thus we've reduced our point addition from three
// multi-precision multiplications to just two! Wow. Just wow.
//
// There is additional overhead introduced by storing the previous lambda
// but it's still a net gain. To additionally decrease memory access
// overhead I never any longer store the Y coordinate. Instead I
// calculate it at the end directly from the lambda and deltaX.
// 
// In addition to this some algebraic re-ordering has been done to move
// constants into the same argument to a new function mp_mod_sub_const
// in hopes that using constant storage instead of private storage
// will aid speeds.
//
// After the above point addition this kernel calculates the public address
// corresponding to the point and stores it in pInverse which is used only
// as interim storage as it won't otherwise be used again this cycle.
//
// One of the scoring kernels will run after this and fetch the address
// from pInverse.
__kernel void profanity_iterate(__global mp_number * const pDeltaX, __global mp_number * const pInverse, __global mp_number * const pPrevLambda) {
	const size_t id = get_global_id(0);

	// negativeGx = 0x8641998106234453aa5f9d6a3178f4f8fd640324d231d726a60d7ea3e907e497
	mp_number negativeGx = { {0xe907e497, 0xa60d7ea3, 0xd231d726, 0xfd640324, 0x3178f4f8, 0xaa5f9d6a, 0x06234453, 0x86419981 } };

	ethhash h = { { 0 } };

	mp_number dX = pDeltaX[id];
	mp_number tmp = pInverse[id];
	mp_number lambda = pPrevLambda[id];

	// λ' = - (2G_y) / d' - λ <=> lambda := pInversedNegativeDoubleGy[id] - pPrevLambda[id]
	mp_mod_sub(&lambda, &tmp, &lambda);

	// λ² = λ * λ <=> tmp := lambda * lambda = λ²
	mp_mod_mul(&tmp, &lambda, &lambda);

	// d' = λ² - d - 3g = (-3g) - (d - λ²) <=> x := tripleNegativeGx - (x - tmp)
	mp_mod_sub(&dX, &dX, &tmp);
	mp_mod_sub_const(&dX, &tripleNegativeGx, &dX);

	pDeltaX[id] = dX;
	pPrevLambda[id] = lambda;

	// Calculate y from dX and lambda
	// y' = (-G_Y) - λ * d' <=> p.y := negativeGy - (p.y * p.x)
	mp_mod_mul(&tmp, &lambda, &dX);
	mp_mod_sub_const(&tmp, &negativeGy, &tmp);

	// Restore X coordinate from delta value
	mp_mod_sub(&dX, &dX, &negativeGx);

	// Initialize Keccak structure with point coordinates in big endian
	h.d[0] = bswap32(dX.d[MP_WORDS - 1]);
	h.d[1] = bswap32(dX.d[MP_WORDS - 2]);
	h.d[2] = bswap32(dX.d[MP_WORDS - 3]);
	h.d[3] = bswap32(dX.d[MP_WORDS - 4]);
	h.d[4] = bswap32(dX.d[MP_WORDS - 5]);
	h.d[5] = bswap32(dX.d[MP_WORDS - 6]);
	h.d[6] = bswap32(dX.d[MP_WORDS - 7]);
	h.d[7] = bswap32(dX.d[MP_WORDS - 8]);
	h.d[8] = bswap32(tmp.d[MP_WORDS - 1]);
	h.d[9] = bswap32(tmp.d[MP_WORDS - 2]);
	h.d[10] = bswap32(tmp.d[MP_WORDS - 3]);
	h.d[11] = bswap32(tmp.d[MP_WORDS - 4]);
	h.d[12] = bswap32(tmp.d[MP_WORDS - 5]);
	h.d[13] = bswap32(tmp.d[MP_WORDS - 6]);
	h.d[14] = bswap32(tmp.d[MP_WORDS - 7]);
	h.d[15] = bswap32(tmp.d[MP_WORDS - 8]);
	h.d[16] ^= 0x01; // length 64

	sha3_keccakf(&h);

	// Save public address hash in pInverse, only used as interim storage until next cycle
	pInverse[id].d[0] = h.d[3];
	pInverse[id].d[1] = h.d[4];
	pInverse[id].d[2] = h.d[5];
	pInverse[id].d[3] = h.d[6];
	pInverse[id].d[4] = h.d[7];
}

void profanity_result_update(const size_t id, __global const uchar * const hash, __global result * const pResult, const uchar score, const uchar scoreMax) {
	if (score && score > scoreMax) {
		uchar hasResult = atomic_inc(&pResult[score].found); // NOTE: If "too many" results are found it'll wrap around to 0 again and overwrite last result. Only relevant if global worksize exceeds MAX(uint).

		// Save only one result for each score, the first.
		if (hasResult == 0) {
			pResult[score].foundId = id;

			for (int i = 0; i < 20; ++i) {
				pResult[score].foundHash[i] = hash[i];
			}
		}
	}
}


__kernel void profanity_score_benchmark(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;
	int score = 0;

	profanity_result_update(id, hash, pResult, score, scoreMax);
}


/* ------------------------------------------------------------------------ */
/* TRON Address Generation and Scoring                                       */
/* ------------------------------------------------------------------------ */

// Base58 alphabet
__constant const char BASE58_ALPHABET[58] = "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz";

// SHA256 constants
__constant const uint SHA256_K[64] = {
	0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
	0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
	0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
	0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
	0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
	0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
	0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
	0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

#define SHA256_ROTR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))
#define SHA256_CH(x, y, z) (((x) & (y)) ^ (~(x) & (z)))
#define SHA256_MAJ(x, y, z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define SHA256_EP0(x) (SHA256_ROTR(x, 2) ^ SHA256_ROTR(x, 13) ^ SHA256_ROTR(x, 22))
#define SHA256_EP1(x) (SHA256_ROTR(x, 6) ^ SHA256_ROTR(x, 11) ^ SHA256_ROTR(x, 25))
#define SHA256_SIG0(x) (SHA256_ROTR(x, 7) ^ SHA256_ROTR(x, 18) ^ ((x) >> 3))
#define SHA256_SIG1(x) (SHA256_ROTR(x, 17) ^ SHA256_ROTR(x, 19) ^ ((x) >> 10))

// SHA256 for small inputs (up to 55 bytes, single block after padding)
void sha256_small(const uchar * data, int len, uchar * hash) {
	uint state[8] = {
		0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
		0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
	};

	// Prepare block with padding
	uchar block[64];
	for (int i = 0; i < 64; ++i) block[i] = 0;
	for (int i = 0; i < len; ++i) block[i] = data[i];
	block[len] = 0x80;

	// Length in bits (big endian) at the end
	ulong bits = (ulong)len * 8;
	block[63] = bits & 0xff;
	block[62] = (bits >> 8) & 0xff;
	block[61] = (bits >> 16) & 0xff;
	block[60] = (bits >> 24) & 0xff;
	block[59] = (bits >> 32) & 0xff;
	block[58] = (bits >> 40) & 0xff;
	block[57] = (bits >> 48) & 0xff;
	block[56] = (bits >> 56) & 0xff;

	// Process block
	uint w[64];
	for (int i = 0; i < 16; ++i) {
		w[i] = ((uint)block[i * 4] << 24) | ((uint)block[i * 4 + 1] << 16) |
		       ((uint)block[i * 4 + 2] << 8) | ((uint)block[i * 4 + 3]);
	}
	for (int i = 16; i < 64; ++i) {
		w[i] = SHA256_SIG1(w[i - 2]) + w[i - 7] + SHA256_SIG0(w[i - 15]) + w[i - 16];
	}

	uint a = state[0], b = state[1], c = state[2], d = state[3];
	uint e = state[4], f = state[5], g = state[6], h = state[7];

	for (int i = 0; i < 64; ++i) {
		uint t1 = h + SHA256_EP1(e) + SHA256_CH(e, f, g) + SHA256_K[i] + w[i];
		uint t2 = SHA256_EP0(a) + SHA256_MAJ(a, b, c);
		h = g; g = f; f = e; e = d + t1;
		d = c; c = b; b = a; a = t1 + t2;
	}

	state[0] += a; state[1] += b; state[2] += c; state[3] += d;
	state[4] += e; state[5] += f; state[6] += g; state[7] += h;

	// Output hash
	for (int i = 0; i < 8; ++i) {
		hash[i * 4] = (state[i] >> 24) & 0xff;
		hash[i * 4 + 1] = (state[i] >> 16) & 0xff;
		hash[i * 4 + 2] = (state[i] >> 8) & 0xff;
		hash[i * 4 + 3] = state[i] & 0xff;
	}
}

// Convert 21-byte TRON address (0x41 + 20-byte hash) to Base58Check
// Returns the length of encoded string, stores result in 'out' (max 34 chars for TRON)
int tron_base58_encode(const uchar * const input21, char * out) {
	// Calculate checksum: double SHA256
	uchar hash1[32], hash2[32];
	sha256_small(input21, 21, hash1);
	sha256_small(hash1, 32, hash2);

	// Create 25-byte data: 21 bytes + 4 byte checksum
	uchar data25[25];
	for (int i = 0; i < 21; ++i) data25[i] = input21[i];
	for (int i = 0; i < 4; ++i) data25[21 + i] = hash2[i];

	// Work with a copy as a big integer (using 32-bit words)
	// 25 bytes = 200 bits, need 7 * 32 = 224 bits
	uint num[7];
	num[0] = ((uint)data25[21] << 24) | ((uint)data25[22] << 16) | ((uint)data25[23] << 8) | (uint)data25[24];
	num[1] = ((uint)data25[17] << 24) | ((uint)data25[18] << 16) | ((uint)data25[19] << 8) | (uint)data25[20];
	num[2] = ((uint)data25[13] << 24) | ((uint)data25[14] << 16) | ((uint)data25[15] << 8) | (uint)data25[16];
	num[3] = ((uint)data25[9] << 24) | ((uint)data25[10] << 16) | ((uint)data25[11] << 8) | (uint)data25[12];
	num[4] = ((uint)data25[5] << 24) | ((uint)data25[6] << 16) | ((uint)data25[7] << 8) | (uint)data25[8];
	num[5] = ((uint)data25[1] << 24) | ((uint)data25[2] << 16) | ((uint)data25[3] << 8) | (uint)data25[4];
	num[6] = (uint)data25[0]; // 0x41 prefix

	char temp[35];
	int len = 0;

	// Convert to base58 by repeatedly dividing by 58
	while (num[6] || num[5] || num[4] || num[3] || num[2] || num[1] || num[0]) {
		ulong remainder = 0;
		for (int i = 6; i >= 0; --i) {
			ulong value = (remainder << 32) | num[i];
			num[i] = (uint)(value / 58);
			remainder = value % 58;
		}
		temp[len++] = BASE58_ALPHABET[remainder];
	}

	// Add leading '1's for leading zero bytes
	int leadingZeros = 0;
	for (int i = 0; i < 25 && data25[i] == 0; ++i) {
		leadingZeros++;
	}
	for (int i = 0; i < leadingZeros; ++i) {
		temp[len++] = '1';
	}

	// Reverse the result
	for (int i = 0; i < len; ++i) {
		out[i] = temp[len - 1 - i];
	}
	out[len] = 0;

	return len;
}

// 豹子号 (Leopard number): Score on trailing repeated characters
// data1[0] = minimum required repeat count
__kernel void profanity_score_tron_repeat(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;

	// Build TRON address: 0x41 prefix + 20 bytes hash
	uchar tronAddr[21];
	tronAddr[0] = 0x41;
	for (int i = 0; i < 20; ++i) {
		tronAddr[i + 1] = hash[i];
	}

	// Encode to Base58
	char base58[35];
	int len = tron_base58_encode(tronAddr, base58);

	// Count trailing repeated characters
	int score = 1;
	if (len >= 2) {
		char lastChar = base58[len - 1];
		for (int i = len - 2; i >= 0; --i) {
			if (base58[i] == lastChar) {
				score++;
			} else {
				break;
			}
		}
	}

	profanity_result_update(id, hash, pResult, score, scoreMax);
}

// 顺子号 (Sequential number): Score on trailing sequential characters (ascending or descending)
// Supports sequences like: 123456, 654321, abcdef, fedcba
__kernel void profanity_score_tron_sequential(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;

	// Build TRON address
	uchar tronAddr[21];
	tronAddr[0] = 0x41;
	for (int i = 0; i < 20; ++i) {
		tronAddr[i + 1] = hash[i];
	}

	// Encode to Base58
	char base58[35];
	int len = tron_base58_encode(tronAddr, base58);

	if (len < 2) {
		profanity_result_update(id, hash, pResult, 0, scoreMax);
		return;
	}

	// Check for ascending sequence from the end
	int ascScore = 1;
	for (int i = len - 2; i >= 0; --i) {
		if (base58[i] + 1 == base58[i + 1]) {
			ascScore++;
		} else {
			break;
		}
	}

	// Check for descending sequence from the end
	int descScore = 1;
	for (int i = len - 2; i >= 0; --i) {
		if (base58[i] - 1 == base58[i + 1]) {
			descScore++;
		} else {
			break;
		}
	}

	int score = ascScore > descScore ? ascScore : descScore;

	profanity_result_update(id, hash, pResult, score, scoreMax);
}

// 谐音靓号/自定义后缀 (Custom suffix matching)
// Supports multiple patterns separated by null bytes
// data1 = pattern bytes (patterns separated by 0x00, max total 20 bytes)
// data2[0] = total length of all patterns including separators
// data2[1] = number of patterns
// Uses 'X' or 'x' as wildcard
__kernel void profanity_score_tron_suffix(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;

	// Build TRON address
	uchar tronAddr[21];
	tronAddr[0] = 0x41;
	for (int i = 0; i < 20; ++i) {
		tronAddr[i + 1] = hash[i];
	}

	// Encode to Base58
	char base58[35];
	int addrLen = tron_base58_encode(tronAddr, base58);

	int totalDataLen = data2[0];
	int numPatterns = data2[1];

	if (numPatterns == 0 || totalDataLen == 0) {
		profanity_result_update(id, hash, pResult, 0, scoreMax);
		return;
	}

	int bestScore = 0;
	int patternStart = 0;

	// Try each pattern
	for (int p = 0; p < numPatterns && patternStart < totalDataLen; ++p) {
		// Find pattern length (until next 0x00 or end)
		int patternLen = 0;
		for (int i = patternStart; i < totalDataLen && data1[i] != 0; ++i) {
			patternLen++;
		}

		if (patternLen > 0 && patternLen <= addrLen) {
			// Match pattern from the end of address
			int score = 0;
			int startPos = addrLen - patternLen;
			bool fullMatch = true;

			for (int i = 0; i < patternLen; ++i) {
				char patternChar = data1[patternStart + i];
				char addrChar = base58[startPos + i];

				// 'X' or 'x' is wildcard
				if (patternChar == 'X' || patternChar == 'x') {
					score++;
				} else if (patternChar == addrChar) {
					score++;
				} else {
					fullMatch = false;
					break;
				}
			}

			// Only count if all pattern characters match
			if (fullMatch && score > bestScore) {
				bestScore = score;
			}
		}

		// Move to next pattern (skip current pattern + null separator)
		patternStart += patternLen + 1;
	}

	profanity_result_update(id, hash, pResult, bestScore, scoreMax);
}

// 多种谐音靓号模式匹配 (Multiple pattern matching for Chinese lucky numbers)
// Matches patterns like: 5211314, 1314521, 168888, 888888, etc.
// data1[0] = number of patterns
// Patterns are pre-defined lucky number combinations
__kernel void profanity_score_tron_lucky(__global mp_number * const pInverse, __global result * const pResult, __constant const uchar * const data1, __constant const uchar * const data2, const uchar scoreMax) {
	const size_t id = get_global_id(0);
	__global const uchar * const hash = pInverse[id].d;

	// Build TRON address
	uchar tronAddr[21];
	tronAddr[0] = 0x41;
	for (int i = 0; i < 20; ++i) {
		tronAddr[i + 1] = hash[i];
	}

	// Encode to Base58
	char base58[35];
	int len = tron_base58_encode(tronAddr, base58);

	int score = 0;

	// Check for lucky patterns at the end
	// Pattern 1: 5211314 (我爱你一生一世) - 7 chars
	if (len >= 7) {
		if (base58[len-7] == '5' && base58[len-6] == '2' && base58[len-5] == '1' &&
		    base58[len-4] == '1' && base58[len-3] == '3' && base58[len-2] == '1' && base58[len-1] == '4') {
			score = 7;
		}
	}

	// Pattern 2: 1314521 (一生一世我爱你)
	if (len >= 7 && score < 7) {
		if (base58[len-7] == '1' && base58[len-6] == '3' && base58[len-5] == '1' &&
		    base58[len-4] == '4' && base58[len-3] == '5' && base58[len-2] == '2' && base58[len-1] == '1') {
			score = 7;
		}
	}

	// Pattern 3: 168888 (一路发发发发) - 6 chars
	if (len >= 6 && score < 6) {
		if (base58[len-6] == '1' && base58[len-5] == '6' &&
		    base58[len-4] == '8' && base58[len-3] == '8' && base58[len-2] == '8' && base58[len-1] == '8') {
			score = 6;
		}
	}

	// Pattern 4: 888888 - 6 chars
	if (len >= 6 && score < 6) {
		bool all8 = true;
		for (int i = 1; i <= 6; ++i) {
			if (base58[len-i] != '8') all8 = false;
		}
		if (all8) score = 6;
	}

	// Pattern 5: 666666 - 6 chars
	if (len >= 6 && score < 6) {
		bool all6 = true;
		for (int i = 1; i <= 6; ++i) {
			if (base58[len-i] != '6') all6 = false;
		}
		if (all6) score = 6;
	}

	// Pattern 6: 520 (我爱你) - 3 chars
	if (len >= 3 && score < 3) {
		if (base58[len-3] == '5' && base58[len-2] == '2' && base58[len-1] == '0') {
			score = 3;
		}
	}

	// Pattern 7: 1314 (一生一世) - 4 chars
	if (len >= 4 && score < 4) {
		if (base58[len-4] == '1' && base58[len-3] == '3' && base58[len-2] == '1' && base58[len-1] == '4') {
			score = 4;
		}
	}

	profanity_result_update(id, hash, pResult, score, scoreMax);
}


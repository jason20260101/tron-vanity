#ifndef HPP_KEYGENERATOR
#define HPP_KEYGENERATOR

#include <string>
#include <random>
#include <sstream>
#include <iomanip>
#include <cstdint>

// secp256k1 curve parameters
// G point (generator)
static const uint8_t SECP256K1_GX[32] = {
    0x79, 0xBE, 0x66, 0x7E, 0xF9, 0xDC, 0xBB, 0xAC,
    0x55, 0xA0, 0x62, 0x95, 0xCE, 0x87, 0x0B, 0x07,
    0x02, 0x9B, 0xFC, 0xDB, 0x2D, 0xCE, 0x28, 0xD9,
    0x59, 0xF2, 0x81, 0x5B, 0x16, 0xF8, 0x17, 0x98
};

static const uint8_t SECP256K1_GY[32] = {
    0x48, 0x3A, 0xDA, 0x77, 0x26, 0xA3, 0xC4, 0x65,
    0x5D, 0xA4, 0xFB, 0xFC, 0x0E, 0x11, 0x08, 0xA8,
    0xFD, 0x17, 0xB4, 0x48, 0xA6, 0x85, 0x54, 0x19,
    0x9C, 0x47, 0xD0, 0x8F, 0xFB, 0x10, 0xD4, 0xB8
};

// Simple 256-bit big integer operations for secp256k1
class BigInt256 {
public:
    uint64_t d[4]; // Little endian

    BigInt256() { d[0] = d[1] = d[2] = d[3] = 0; }

    BigInt256(const uint8_t* bytes32) {
        // Big endian bytes to little endian words
        for (int i = 0; i < 4; i++) {
            d[3-i] = ((uint64_t)bytes32[i*8+0] << 56) | ((uint64_t)bytes32[i*8+1] << 48) |
                     ((uint64_t)bytes32[i*8+2] << 40) | ((uint64_t)bytes32[i*8+3] << 32) |
                     ((uint64_t)bytes32[i*8+4] << 24) | ((uint64_t)bytes32[i*8+5] << 16) |
                     ((uint64_t)bytes32[i*8+6] << 8)  | ((uint64_t)bytes32[i*8+7]);
        }
    }

    void toBytes(uint8_t* bytes32) const {
        for (int i = 0; i < 4; i++) {
            bytes32[(3-i)*8+0] = (d[i] >> 56) & 0xFF;
            bytes32[(3-i)*8+1] = (d[i] >> 48) & 0xFF;
            bytes32[(3-i)*8+2] = (d[i] >> 40) & 0xFF;
            bytes32[(3-i)*8+3] = (d[i] >> 32) & 0xFF;
            bytes32[(3-i)*8+4] = (d[i] >> 24) & 0xFF;
            bytes32[(3-i)*8+5] = (d[i] >> 16) & 0xFF;
            bytes32[(3-i)*8+6] = (d[i] >> 8) & 0xFF;
            bytes32[(3-i)*8+7] = d[i] & 0xFF;
        }
    }

    std::string toHex() const {
        uint8_t bytes[32];
        toBytes(bytes);
        std::ostringstream ss;
        ss << std::hex << std::setfill('0');
        for (int i = 0; i < 32; i++) {
            ss << std::setw(2) << (int)bytes[i];
        }
        return ss.str();
    }

    bool isZero() const {
        return d[0] == 0 && d[1] == 0 && d[2] == 0 && d[3] == 0;
    }

    int compare(const BigInt256& other) const {
        for (int i = 3; i >= 0; i--) {
            if (d[i] > other.d[i]) return 1;
            if (d[i] < other.d[i]) return -1;
        }
        return 0;
    }

    // Add with carry
    BigInt256 add(const BigInt256& other) const {
        BigInt256 r;
        uint64_t carry = 0;
        for (int i = 0; i < 4; i++) {
            __uint128_t sum = (__uint128_t)d[i] + other.d[i] + carry;
            r.d[i] = (uint64_t)sum;
            carry = (uint64_t)(sum >> 64);
        }
        return r;
    }

    // Subtract with borrow
    BigInt256 sub(const BigInt256& other) const {
        BigInt256 r;
        uint64_t borrow = 0;
        for (int i = 0; i < 4; i++) {
            __uint128_t diff = (__uint128_t)d[i] - other.d[i] - borrow;
            r.d[i] = (uint64_t)diff;
            borrow = (diff >> 64) ? 1 : 0;
        }
        return r;
    }

    // Multiply and return lower 256 bits
    BigInt256 mulLow(const BigInt256& other) const {
        BigInt256 r;
        __uint128_t acc[4] = {0, 0, 0, 0};

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i + j < 4) {
                    __uint128_t prod = (__uint128_t)d[i] * other.d[j];
                    acc[i + j] += (uint64_t)prod;
                    if (i + j + 1 < 4) {
                        acc[i + j + 1] += (uint64_t)(prod >> 64);
                    }
                }
            }
        }

        // Propagate carries
        uint64_t carry = 0;
        for (int i = 0; i < 4; i++) {
            acc[i] += carry;
            r.d[i] = (uint64_t)acc[i];
            carry = (uint64_t)(acc[i] >> 64);
        }

        return r;
    }
};

// secp256k1 field modulus p = 2^256 - 2^32 - 977
static const BigInt256 SECP256K1_P = []() {
    BigInt256 p;
    p.d[0] = 0xFFFFFFFEFFFFFC2FULL;
    p.d[1] = 0xFFFFFFFFFFFFFFFFULL;
    p.d[2] = 0xFFFFFFFFFFFFFFFFULL;
    p.d[3] = 0xFFFFFFFFFFFFFFFFULL;
    return p;
}();

// secp256k1 curve order n
static const BigInt256 SECP256K1_N = []() {
    BigInt256 n;
    n.d[0] = 0xBFD25E8CD0364141ULL;
    n.d[1] = 0xBAAEDCE6AF48A03BULL;
    n.d[2] = 0xFFFFFFFFFFFFFFFEULL;
    n.d[3] = 0xFFFFFFFFFFFFFFFFULL;
    return n;
}();

// Field element modular reduction (simple, not optimized)
inline BigInt256 modP(const BigInt256& a) {
    BigInt256 r = a;
    while (r.compare(SECP256K1_P) >= 0) {
        r = r.sub(SECP256K1_P);
    }
    return r;
}

// Modular addition
inline BigInt256 addMod(const BigInt256& a, const BigInt256& b) {
    BigInt256 r = a.add(b);
    if (r.compare(SECP256K1_P) >= 0) {
        r = r.sub(SECP256K1_P);
    }
    return r;
}

// Modular subtraction
inline BigInt256 subMod(const BigInt256& a, const BigInt256& b) {
    if (a.compare(b) < 0) {
        return a.add(SECP256K1_P).sub(b);
    }
    return a.sub(b);
}

// Modular multiplication (simple)
inline BigInt256 mulMod(const BigInt256& a, const BigInt256& b) {
    // For simplicity, use double-and-add method
    BigInt256 result;
    BigInt256 temp = a;

    for (int i = 0; i < 4; i++) {
        uint64_t bits = b.d[i];
        for (int j = 0; j < 64; j++) {
            if (bits & 1) {
                result = addMod(result, temp);
            }
            temp = addMod(temp, temp);
            bits >>= 1;
        }
    }
    return result;
}

// Modular inverse using Fermat's little theorem: a^(-1) = a^(p-2) mod p
inline BigInt256 invMod(const BigInt256& a) {
    // p - 2
    BigInt256 exp = SECP256K1_P;
    exp.d[0] -= 2;

    BigInt256 result;
    result.d[0] = 1;
    BigInt256 base = a;

    for (int i = 0; i < 4; i++) {
        uint64_t bits = exp.d[i];
        for (int j = 0; j < 64; j++) {
            if (bits & 1) {
                result = mulMod(result, base);
            }
            base = mulMod(base, base);
            bits >>= 1;
        }
    }
    return result;
}

// Point on secp256k1 curve
struct ECPoint {
    BigInt256 x, y;
    bool infinity;

    ECPoint() : infinity(true) {}
    ECPoint(const BigInt256& _x, const BigInt256& _y) : x(_x), y(_y), infinity(false) {}
};

// Point doubling
inline ECPoint pointDouble(const ECPoint& p) {
    if (p.infinity) return p;

    // lambda = (3 * x^2) / (2 * y)
    BigInt256 x2 = mulMod(p.x, p.x);
    BigInt256 num = addMod(addMod(x2, x2), x2); // 3 * x^2
    BigInt256 denom = addMod(p.y, p.y); // 2 * y
    BigInt256 lambda = mulMod(num, invMod(denom));

    // x3 = lambda^2 - 2*x
    BigInt256 lambda2 = mulMod(lambda, lambda);
    BigInt256 x3 = subMod(subMod(lambda2, p.x), p.x);

    // y3 = lambda * (x - x3) - y
    BigInt256 y3 = subMod(mulMod(lambda, subMod(p.x, x3)), p.y);

    return ECPoint(x3, y3);
}

// Point addition
inline ECPoint pointAdd(const ECPoint& p, const ECPoint& q) {
    if (p.infinity) return q;
    if (q.infinity) return p;

    if (p.x.compare(q.x) == 0) {
        if (p.y.compare(q.y) == 0) {
            return pointDouble(p);
        } else {
            return ECPoint(); // Point at infinity
        }
    }

    // lambda = (y2 - y1) / (x2 - x1)
    BigInt256 num = subMod(q.y, p.y);
    BigInt256 denom = subMod(q.x, p.x);
    BigInt256 lambda = mulMod(num, invMod(denom));

    // x3 = lambda^2 - x1 - x2
    BigInt256 lambda2 = mulMod(lambda, lambda);
    BigInt256 x3 = subMod(subMod(lambda2, p.x), q.x);

    // y3 = lambda * (x1 - x3) - y1
    BigInt256 y3 = subMod(mulMod(lambda, subMod(p.x, x3)), p.y);

    return ECPoint(x3, y3);
}

// Scalar multiplication using double-and-add
inline ECPoint scalarMult(const ECPoint& p, const BigInt256& k) {
    ECPoint result;
    ECPoint temp = p;

    for (int i = 0; i < 4; i++) {
        uint64_t bits = k.d[i];
        for (int j = 0; j < 64; j++) {
            if (bits & 1) {
                result = pointAdd(result, temp);
            }
            temp = pointDouble(temp);
            bits >>= 1;
        }
    }
    return result;
}

class KeyGenerator {
public:
    std::string privateKey;
    std::string publicKey;

    void generate() {
        // Generate random 256-bit private key
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<uint64_t> dis;

        BigInt256 privKey;
        do {
            privKey.d[0] = dis(gen);
            privKey.d[1] = dis(gen);
            privKey.d[2] = dis(gen);
            privKey.d[3] = dis(gen);
        } while (privKey.isZero() || privKey.compare(SECP256K1_N) >= 0);

        privateKey = privKey.toHex();

        // Compute public key = privKey * G
        BigInt256 gx(SECP256K1_GX);
        BigInt256 gy(SECP256K1_GY);
        ECPoint G(gx, gy);

        ECPoint pubPoint = scalarMult(G, privKey);

        // Public key is x || y (uncompressed, without 04 prefix)
        publicKey = pubPoint.x.toHex() + pubPoint.y.toHex();
    }
};

#endif /* HPP_KEYGENERATOR */


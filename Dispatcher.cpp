#include "Dispatcher.hpp"

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <cstring>

#if defined(__APPLE__) || defined(__MACOSX)
#include <machine/endian.h>
#else
#include <arpa/inet.h>
#endif

#include "precomp.hpp"

#ifndef htonll
#define htonll(x) ((((uint64_t)htonl(x)) << 32) | htonl((x) >> 32))
#endif

// Base58 alphabet
static const char* BASE58_ALPHABET = "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz";

// SHA256 implementation for checksum
static void sha256_transform(uint32_t state[8], const uint8_t data[64]) {
	static const uint32_t k[64] = {
		0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
		0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
		0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
		0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
		0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
		0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
		0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
		0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
	};

	#define ROTR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))
	#define CH(x, y, z) (((x) & (y)) ^ (~(x) & (z)))
	#define MAJ(x, y, z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
	#define EP0(x) (ROTR(x, 2) ^ ROTR(x, 13) ^ ROTR(x, 22))
	#define EP1(x) (ROTR(x, 6) ^ ROTR(x, 11) ^ ROTR(x, 25))
	#define SIG0(x) (ROTR(x, 7) ^ ROTR(x, 18) ^ ((x) >> 3))
	#define SIG1(x) (ROTR(x, 17) ^ ROTR(x, 19) ^ ((x) >> 10))

	uint32_t w[64];
	for (int i = 0; i < 16; i++) {
		w[i] = ((uint32_t)data[i * 4] << 24) | ((uint32_t)data[i * 4 + 1] << 16) |
		       ((uint32_t)data[i * 4 + 2] << 8) | ((uint32_t)data[i * 4 + 3]);
	}
	for (int i = 16; i < 64; i++) {
		w[i] = SIG1(w[i - 2]) + w[i - 7] + SIG0(w[i - 15]) + w[i - 16];
	}

	uint32_t a = state[0], b = state[1], c = state[2], d = state[3];
	uint32_t e = state[4], f = state[5], g = state[6], h = state[7];

	for (int i = 0; i < 64; i++) {
		uint32_t t1 = h + EP1(e) + CH(e, f, g) + k[i] + w[i];
		uint32_t t2 = EP0(a) + MAJ(a, b, c);
		h = g; g = f; f = e; e = d + t1;
		d = c; c = b; b = a; a = t1 + t2;
	}

	state[0] += a; state[1] += b; state[2] += c; state[3] += d;
	state[4] += e; state[5] += f; state[6] += g; state[7] += h;

	#undef ROTR
	#undef CH
	#undef MAJ
	#undef EP0
	#undef EP1
	#undef SIG0
	#undef SIG1
}

static void sha256(const uint8_t* data, size_t len, uint8_t hash[32]) {
	uint32_t state[8] = {
		0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
		0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
	};

	uint8_t block[64];
	size_t i = 0;

	while (len >= 64) {
		sha256_transform(state, data + i);
		i += 64;
		len -= 64;
	}

	memset(block, 0, 64);
	memcpy(block, data + i, len);
	block[len] = 0x80;

	if (len >= 56) {
		sha256_transform(state, block);
		memset(block, 0, 64);
	}

	uint64_t bits = (i + len) * 8;
	block[63] = bits & 0xff;
	block[62] = (bits >> 8) & 0xff;
	block[61] = (bits >> 16) & 0xff;
	block[60] = (bits >> 24) & 0xff;
	block[59] = (bits >> 32) & 0xff;
	block[58] = (bits >> 40) & 0xff;
	block[57] = (bits >> 48) & 0xff;
	block[56] = (bits >> 56) & 0xff;

	sha256_transform(state, block);

	for (int j = 0; j < 8; j++) {
		hash[j * 4] = (state[j] >> 24) & 0xff;
		hash[j * 4 + 1] = (state[j] >> 16) & 0xff;
		hash[j * 4 + 2] = (state[j] >> 8) & 0xff;
		hash[j * 4 + 3] = state[j] & 0xff;
	}
}

// Base58Check encode for TRON address
static std::string toBase58Check(const uint8_t* data21) {
	// Double SHA256 for checksum
	uint8_t hash1[32], hash2[32];
	sha256(data21, 21, hash1);
	sha256(hash1, 32, hash2);

	// Append 4-byte checksum
	uint8_t data25[25];
	memcpy(data25, data21, 21);
	memcpy(data25 + 21, hash2, 4);

	// Count leading zeros
	int leadingZeros = 0;
	for (int i = 0; i < 25 && data25[i] == 0; i++) {
		leadingZeros++;
	}

	// Convert to base58
	std::string result;

	// Use big integer division
	uint8_t temp[25];
	memcpy(temp, data25, 25);

	while (true) {
		bool allZero = true;
		for (int i = 0; i < 25; i++) {
			if (temp[i] != 0) {
				allZero = false;
				break;
			}
		}
		if (allZero) break;

		int remainder = 0;
		for (int i = 0; i < 25; i++) {
			int value = remainder * 256 + temp[i];
			temp[i] = value / 58;
			remainder = value % 58;
		}
		result = BASE58_ALPHABET[remainder] + result;
	}

	// Add leading '1's
	for (int i = 0; i < leadingZeros; i++) {
		result = '1' + result;
	}

	return result;
}

static std::string::size_type fromHex(char c) {
	if (c >= 'A' && c <= 'F') {
		c += 'a' - 'A';
	}

	const std::string hex = "0123456789abcdef";
	const std::string::size_type ret = hex.find(c);
	return ret;
}

static cl_ulong4 fromHex(const std::string & strHex) {
	uint8_t data[32];
	std::fill(data, data + sizeof(data), cl_uchar(0));

	auto index = 0;
	for(size_t i = 0; i < strHex.size(); i += 2) {
		const auto indexHi = fromHex(strHex[i]);
		const auto indexLo = i + 1 < strHex.size() ? fromHex(strHex[i+1]) : std::string::npos;

		const auto valHi = (indexHi == std::string::npos) ? 0 : indexHi << 4;
		const auto valLo = (indexLo == std::string::npos) ? 0 : indexLo;

		data[index] = valHi | valLo;
		++index;
	}

	cl_ulong4 res = {
		.s = {
			htonll(*(uint64_t *)(data + 24)),
			htonll(*(uint64_t *)(data + 16)),
			htonll(*(uint64_t *)(data + 8)),
			htonll(*(uint64_t *)(data + 0)),
		}
	};
	return res;
}

// Add two hex private keys (without 0x prefix), return result mod N (secp256k1 curve order)
static std::string addPrivateKeys(const std::string& keyA, const std::string& keyB) {
	// secp256k1 curve order N
	static const char* N_HEX = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";

	// Convert hex string to bytes
	auto hexToBytes = [](const std::string& hex, uint8_t* bytes, size_t len) {
		for (size_t i = 0; i < len && i * 2 < hex.size(); i++) {
			char hi = hex[i * 2];
			char lo = (i * 2 + 1 < hex.size()) ? hex[i * 2 + 1] : '0';
			auto hexVal = [](char c) -> uint8_t {
				if (c >= '0' && c <= '9') return c - '0';
				if (c >= 'a' && c <= 'f') return c - 'a' + 10;
				if (c >= 'A' && c <= 'F') return c - 'A' + 10;
				return 0;
			};
			bytes[i] = (hexVal(hi) << 4) | hexVal(lo);
		}
	};

	uint8_t a[32] = {0}, b[32] = {0}, n[32] = {0};
	hexToBytes(keyA, a, 32);
	hexToBytes(keyB, b, 32);
	hexToBytes(N_HEX, n, 32);

	// Add a + b
	uint8_t sum[33] = {0};
	uint16_t carry = 0;
	for (int i = 31; i >= 0; i--) {
		uint16_t s = (uint16_t)a[i] + (uint16_t)b[i] + carry;
		sum[i + 1] = s & 0xFF;
		carry = s >> 8;
	}
	sum[0] = carry;

	// Compare sum with N
	auto compare = [](const uint8_t* x, size_t xLen, const uint8_t* y, size_t yLen) -> int {
		// Skip leading zeros
		while (xLen > 0 && x[0] == 0) { x++; xLen--; }
		while (yLen > 0 && y[0] == 0) { y++; yLen--; }
		if (xLen != yLen) return (xLen > yLen) ? 1 : -1;
		for (size_t i = 0; i < xLen; i++) {
			if (x[i] != y[i]) return (x[i] > y[i]) ? 1 : -1;
		}
		return 0;
	};

	// If sum >= N, subtract N
	uint8_t* result = sum + 1;  // Point to 32-byte result
	if (compare(sum, 33, n, 32) >= 0) {
		uint16_t borrow = 0;
		for (int i = 31; i >= 0; i--) {
			int16_t diff = (int16_t)sum[i + 1] - (int16_t)n[i] - borrow;
			if (diff < 0) {
				diff += 256;
				borrow = 1;
			} else {
				borrow = 0;
			}
			result[i] = diff & 0xFF;
		}
	}

	// Convert back to hex
	std::ostringstream ss;
	ss << std::hex << std::setfill('0');
	for (int i = 0; i < 32; i++) {
		ss << std::setw(2) << (int)result[i];
	}
	return ss.str();
}

// Save address and private key to file
static void saveToFile(const std::string& address, const std::string& privateKey) {
	std::string outputDir = "output";
	#if defined(__APPLE__) || defined(__MACOSX) || defined(__linux__)
	system(("mkdir -p " + outputDir).c_str());
	#else
	system(("mkdir " + outputDir + " 2>nul").c_str());
	#endif

	std::string filename = outputDir + "/" + address + ".txt";
	std::ofstream outFile(filename);
	if (outFile.is_open()) {
		outFile << privateKey << std::endl;
		outFile.close();
		std::cout << "  已保存: " << filename << std::endl;
	} else {
		std::cerr << "  错误: 无法保存文件: " << filename << std::endl;
	}
}

static void printResult(cl_ulong4 seed, cl_ulong round, result r, cl_uchar score, const std::chrono::time_point<std::chrono::steady_clock> & timeStart, const Mode & mode, const std::string & seedPrivateKey = "") {
	const auto seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - timeStart).count();

	// Format private key offset
	cl_ulong carry = 0;
	cl_ulong4 seedRes;

	seedRes.s[0] = seed.s[0] + round; carry = seedRes.s[0] < round;
	seedRes.s[1] = seed.s[1] + carry; carry = !seedRes.s[1];
	seedRes.s[2] = seed.s[2] + carry; carry = !seedRes.s[2];
	seedRes.s[3] = seed.s[3] + carry + r.foundId;

	std::ostringstream ss;
	ss << std::hex << std::setfill('0');
	ss << std::setw(16) << seedRes.s[3] << std::setw(16) << seedRes.s[2] << std::setw(16) << seedRes.s[1] << std::setw(16) << seedRes.s[0];
	const std::string strPrivate = ss.str();

	// Generate TRON address (Base58Check encoding)
	uint8_t tronAddr[21];
	tronAddr[0] = 0x41;
	for (int i = 0; i < 20; i++) {
		tronAddr[i + 1] = r.foundHash[i];
	}
	std::string strAddress = toBase58Check(tronAddr);

	// Print result
	const std::string strVT100ClearLine = "\33[2K\r";
	std::cout << strVT100ClearLine << "  时间: " << std::setw(5) << seconds << "s 分数: " << std::setw(2) << (int)score
	          << " 地址: " << strAddress << std::endl;

	// Calculate final private key and save to file
	if (!seedPrivateKey.empty()) {
		std::string finalPrivateKey = addPrivateKeys(seedPrivateKey, strPrivate);
		// 加密显示私钥：只显示前6位和后6位，中间用*号代替
		std::string maskedKey = finalPrivateKey.substr(0, 6) + std::string(52, '*') + finalPrivateKey.substr(58, 6);
		std::cout << "  私钥: 0x" << maskedKey << std::endl;
		saveToFile(strAddress, "0x" + finalPrivateKey);
	}
}

unsigned int getKernelExecutionTimeMicros(cl_event & e) {
	cl_ulong timeStart = 0, timeEnd = 0;
	clWaitForEvents(1, &e);
	clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof(timeStart), &timeStart, NULL);
	clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_END, sizeof(timeEnd), &timeEnd, NULL);
	return (timeEnd - timeStart) / 1000;
}

Dispatcher::OpenCLException::OpenCLException(const std::string s, const cl_int res) :
	std::runtime_error( s + " (res = " + toString(res) + ")"),
	m_res(res)
{

}

void Dispatcher::OpenCLException::OpenCLException::throwIfError(const std::string s, const cl_int res) {
	if (res != CL_SUCCESS) {
		throw OpenCLException(s, res);
	}
}

cl_command_queue Dispatcher::Device::createQueue(cl_context & clContext, cl_device_id & clDeviceId) {
#ifdef PROFANITY_DEBUG
	cl_command_queue_properties p = CL_QUEUE_PROFILING_ENABLE;
#else
	cl_command_queue_properties p = 0;
#endif

#ifdef CL_VERSION_2_0
	const cl_command_queue ret = clCreateCommandQueueWithProperties(clContext, clDeviceId, &p, NULL);
#else
	const cl_command_queue ret = clCreateCommandQueue(clContext, clDeviceId, p, NULL);
#endif
	return ret == NULL ? throw std::runtime_error("failed to create command queue") : ret;
}

cl_kernel Dispatcher::Device::createKernel(cl_program & clProgram, const std::string s) {
	cl_kernel ret  = clCreateKernel(clProgram, s.c_str(), NULL);
	return ret == NULL ? throw std::runtime_error("failed to create kernel \"" + s + "\"") : ret;
}

cl_ulong4 Dispatcher::Device::createSeed() {
#ifdef PROFANITY_DEBUG
	cl_ulong4 r;
	r.s[0] = 1;
	r.s[1] = 1;
	r.s[2] = 1;
	r.s[3] = 1;
	return r;
#else
	// We do not need really safe crypto random here, since we inherit safety
	// of the key from the user-provided seed public key.
	// We only need this random to not repeat same job among different devices
	std::random_device rd;

	cl_ulong4 diff;
	diff.s[0] = (((uint64_t)rd()) << 32) | rd();
	diff.s[1] = (((uint64_t)rd()) << 32) | rd();
	diff.s[2] = (((uint64_t)rd()) << 32) | rd();
	diff.s[3] = (((uint64_t)rd() & 0x0000ffff) << 32) | rd(); // zeroing 2 highest bytes to prevent overflowing sum private key after adding to seed private key
	return diff;
#endif
}

Dispatcher::Device::Device(Dispatcher & parent, cl_context & clContext, cl_program & clProgram, cl_device_id clDeviceId, const size_t worksizeLocal, const size_t size, const size_t index, const Mode & mode, cl_ulong4 clSeedX, cl_ulong4 clSeedY) :
	m_parent(parent),
	m_index(index),
	m_clDeviceId(clDeviceId),
	m_worksizeLocal(worksizeLocal),
	m_clScoreMax(0),
	m_clQueue(createQueue(clContext, clDeviceId)),
	m_kernelInit(createKernel(clProgram, "profanity_init")),
	m_kernelInverse(createKernel(clProgram, "profanity_inverse")),
	m_kernelIterate(createKernel(clProgram, "profanity_iterate")),
	m_kernelScore(createKernel(clProgram, mode.kernel)),
	m_memPrecomp(clContext, m_clQueue, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, sizeof(g_precomp), g_precomp),
	m_memPointsDeltaX(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, size, true),
	m_memInversedNegativeDoubleGy(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, size, true),
	m_memPrevLambda(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, size, true),
	m_memResult(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY, PROFANITY_MAX_SCORE + 1),
	m_memData1(clContext, m_clQueue, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, 20),
	m_memData2(clContext, m_clQueue, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, 20),
	m_clSeed(createSeed()),
	m_clSeedX(clSeedX),
	m_clSeedY(clSeedY),
	m_round(0),
	m_speed(PROFANITY_SPEEDSAMPLES),
	m_sizeInitialized(0),
	m_eventFinished(NULL)
{

}

Dispatcher::Device::~Device() {

}

Dispatcher::Dispatcher(cl_context & clContext, cl_program & clProgram, const Mode mode, const size_t worksizeMax, const size_t inverseSize, const size_t inverseMultiple, const cl_uchar clScoreQuit, const std::string & seedPublicKey, const std::string & seedPrivateKey)
	: m_clContext(clContext)
	  , m_clProgram(clProgram)
	  , m_mode(mode)
	  , m_worksizeMax(worksizeMax)
	  , m_inverseSize(inverseSize)
	  , m_size(inverseSize * inverseMultiple)
	  , m_clScoreMax(mode.score)
	  , m_clScoreQuit(clScoreQuit)
	  , m_eventFinished(nullptr)
	  , m_countPrint(0)
	  , m_countRunning(0), m_sizeInitTotal(0), m_sizeInitDone(0), m_quit(false),
	  m_publicKeyX(fromHex(seedPublicKey.substr(0, 64)))
	  , m_publicKeyY(fromHex(seedPublicKey.substr(64, 64)))
	  , m_seedPrivateKey(seedPrivateKey) {
}

Dispatcher::~Dispatcher() = default;

void Dispatcher::addDevice(cl_device_id clDeviceId, const size_t worksizeLocal, const size_t index) {
	auto * pDevice = new Device(*this, m_clContext, m_clProgram, clDeviceId, worksizeLocal, m_size, index, m_mode, m_publicKeyX, m_publicKeyY);
	m_vDevices.push_back(pDevice);
}

void Dispatcher::run() {
	m_eventFinished = clCreateUserEvent(m_clContext, NULL);
	timeStart = std::chrono::steady_clock::now();

	init();

	const auto timeInitialization = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - timeStart).count();
	std::cout << "初始化耗时: " << timeInitialization << " 秒" << std::endl;

	m_quit = false;
	m_countRunning = m_vDevices.size();

	std::cout << "开始搜索..." << std::endl;
	std::cout << std::endl;

	for (const auto & m_vDevice : m_vDevices) {
		dispatch(*m_vDevice);
	}

	clWaitForEvents(1, &m_eventFinished);
	clReleaseEvent(m_eventFinished);
	m_eventFinished = nullptr;
}

void Dispatcher::init() {
	std::cout << "初始化设备..." << std::endl;

	const auto deviceCount = m_vDevices.size();
	m_sizeInitTotal = m_size * deviceCount;
	m_sizeInitDone = 0;

	cl_event * const pInitEvents = new cl_event[deviceCount];

	for (size_t i = 0; i < deviceCount; ++i) {
		pInitEvents[i] = clCreateUserEvent(m_clContext, NULL);
		m_vDevices[i]->m_eventFinished = pInitEvents[i];
		initBegin(*m_vDevices[i]);
	}

	clWaitForEvents(deviceCount, pInitEvents);
	for (size_t i = 0; i < deviceCount; ++i) {
		m_vDevices[i]->m_eventFinished = NULL;
		clReleaseEvent(pInitEvents[i]);
	}

	delete[] pInitEvents;

	std::cout << std::endl;
}

void Dispatcher::initBegin(Device & d) {
	// Set mode data
	for (auto i = 0; i < 20; ++i) {
		d.m_memData1[i] = m_mode.data1[i];
		d.m_memData2[i] = m_mode.data2[i];
	}

	// Write precompute table and mode data
	d.m_memPrecomp.write(true);
	d.m_memData1.write(true);
	d.m_memData2.write(true);

	// Kernel arguments - profanity_begin
	d.m_memPrecomp.setKernelArg(d.m_kernelInit, 0);
	d.m_memPointsDeltaX.setKernelArg(d.m_kernelInit, 1);
	d.m_memPrevLambda.setKernelArg(d.m_kernelInit, 2);
	d.m_memResult.setKernelArg(d.m_kernelInit, 3);
	CLMemory<cl_ulong4>::setKernelArg(d.m_kernelInit, 4, d.m_clSeed);
	CLMemory<cl_ulong4>::setKernelArg(d.m_kernelInit, 5, d.m_clSeedX);
	CLMemory<cl_ulong4>::setKernelArg(d.m_kernelInit, 6, d.m_clSeedY);

	// Kernel arguments - profanity_inverse
	d.m_memPointsDeltaX.setKernelArg(d.m_kernelInverse, 0);
	d.m_memInversedNegativeDoubleGy.setKernelArg(d.m_kernelInverse, 1);

	// Kernel arguments - profanity_iterate
	d.m_memPointsDeltaX.setKernelArg(d.m_kernelIterate, 0);
	d.m_memInversedNegativeDoubleGy.setKernelArg(d.m_kernelIterate, 1);
	d.m_memPrevLambda.setKernelArg(d.m_kernelIterate, 2);


	// Kernel arguments - profanity_score_*
	d.m_memInversedNegativeDoubleGy.setKernelArg(d.m_kernelScore, 0);
	d.m_memResult.setKernelArg(d.m_kernelScore, 1);
	d.m_memData1.setKernelArg(d.m_kernelScore, 2);
	d.m_memData2.setKernelArg(d.m_kernelScore, 3);

	CLMemory<cl_uchar>::setKernelArg(d.m_kernelScore, 4, d.m_clScoreMax); // Updated in handleResult()

	// Seed device
	initContinue(d);
}

void Dispatcher::initContinue(Device & d) {
	size_t sizeLeft = m_size - d.m_sizeInitialized;
	const size_t sizeInitLimit = m_size / 20;

	// Print progress
	const size_t percentDone = m_sizeInitDone * 100 / m_sizeInitTotal;
	std::cout << "  " << percentDone << "%\r" << std::flush;

	if (sizeLeft) {
		cl_event event;
		const size_t sizeRun = std::min(sizeInitLimit, std::min(sizeLeft, m_worksizeMax));
		const auto resEnqueue = clEnqueueNDRangeKernel(d.m_clQueue, d.m_kernelInit, 1, &d.m_sizeInitialized, &sizeRun, NULL, 0, NULL, &event);
		OpenCLException::throwIfError("kernel queueing failed during initilization", resEnqueue);

		// See: https://www.khronos.org/registry/OpenCL/sdk/1.2/docs/man/xhtml/clSetEventCallback.html
		// If an application needs to wait for completion of a routine from the above list in a callback, please use the non-blocking form of the function, and
		// assign a completion callback to it to do the remainder of your work. Note that when a callback (or other code) enqueues commands to a command-queue,
		// the commands are not required to begin execution until the queue is flushed. In standard usage, blocking enqueue calls serve this role by implicitly
		// flushing the queue. Since blocking calls are not permitted in callbacks, those callbacks that enqueue commands on a command queue should either call
		// clFlush on the queue before returning or arrange for clFlush to be called later on another thread.
		clFlush(d.m_clQueue); 

		std::lock_guard<std::mutex> lock(m_mutex);
		d.m_sizeInitialized += sizeRun;
		m_sizeInitDone += sizeRun;

		const auto resCallback = clSetEventCallback(event, CL_COMPLETE, staticCallback, &d);
		OpenCLException::throwIfError("failed to set custom callback during initialization", resCallback);
	} else {
		const std::string strOutput = "  GPU" + toString(d.m_index) + " 初始化完成";
		std::cout << strOutput << std::endl;
		clSetUserEventStatus(d.m_eventFinished, CL_COMPLETE);
	}
}

void Dispatcher::enqueueKernel(cl_command_queue & clQueue, cl_kernel & clKernel, size_t worksizeGlobal, const size_t worksizeLocal, cl_event * pEvent = NULL) {
	const size_t worksizeMax = m_worksizeMax;
	size_t worksizeOffset = 0;
	while (worksizeGlobal) {
		const size_t worksizeRun = std::min(worksizeGlobal, worksizeMax);
		const size_t * const pWorksizeLocal = (worksizeLocal == 0 ? NULL : &worksizeLocal);
		const auto res = clEnqueueNDRangeKernel(clQueue, clKernel, 1, &worksizeOffset, &worksizeRun, pWorksizeLocal, 0, NULL, pEvent);
		OpenCLException::throwIfError("kernel queueing failed", res);

		worksizeGlobal -= worksizeRun;
		worksizeOffset += worksizeRun;
	}
}

void Dispatcher::enqueueKernelDevice(Device & d, cl_kernel & clKernel, size_t worksizeGlobal, cl_event * pEvent = NULL) {
	try {
		enqueueKernel(d.m_clQueue, clKernel, worksizeGlobal, d.m_worksizeLocal, pEvent);
	} catch ( OpenCLException & e ) {
		// If local work size is invalid, abandon it and let implementation decide
		if ((e.m_res == CL_INVALID_WORK_GROUP_SIZE || e.m_res == CL_INVALID_WORK_ITEM_SIZE) && d.m_worksizeLocal != 0) {
			std::cout << std::endl << "warning: local work size abandoned on GPU" << d.m_index << std::endl;
			d.m_worksizeLocal = 0;
			enqueueKernel(d.m_clQueue, clKernel, worksizeGlobal, d.m_worksizeLocal, pEvent);
		}
		else {
			throw;
		}
	}
}

void Dispatcher::dispatch(Device & d) {
	cl_event event;
	d.m_memResult.read(false, &event);

#ifdef PROFANITY_DEBUG
	cl_event eventInverse;
	cl_event eventIterate;

	enqueueKernelDevice(d, d.m_kernelInverse, m_size / m_inverseSize, &eventInverse);
	enqueueKernelDevice(d, d.m_kernelIterate, m_size, &eventIterate);
#else
	enqueueKernelDevice(d, d.m_kernelInverse, m_size / m_inverseSize);
	enqueueKernelDevice(d, d.m_kernelIterate, m_size);
#endif

	enqueueKernelDevice(d, d.m_kernelScore, m_size);
	clFlush(d.m_clQueue);

#ifdef PROFANITY_DEBUG
	clFinish(d.m_clQueue);
	std::cout << "Timing: profanity_inverse = " << getKernelExecutionTimeMicros(eventInverse) << "us, profanity_iterate = " << getKernelExecutionTimeMicros(eventIterate) << "us" << std::endl;
#endif

	const auto res = clSetEventCallback(event, CL_COMPLETE, staticCallback, &d);
	OpenCLException::throwIfError("failed to set custom callback", res);
}

void Dispatcher::handleResult(Device & d) {
	for (auto i = PROFANITY_MAX_SCORE; i > m_clScoreMax; --i) {
		result & r = d.m_memResult[i];

		if (r.found > 0 && i >= d.m_clScoreMax) {
			d.m_clScoreMax = i;
			CLMemory<cl_uchar>::setKernelArg(d.m_kernelScore, 4, d.m_clScoreMax);

			std::lock_guard<std::mutex> lock(m_mutex);
			if (i >= m_clScoreMax) {
				m_clScoreMax = i;

				if (m_clScoreQuit && i >= m_clScoreQuit) {
					m_quit = true;
				}

				++m_totalFound;
				printResult(d.m_clSeed, d.m_round, r, i, timeStart, m_mode, m_seedPrivateKey);
			}

			break;
		}
	}
}

void Dispatcher::onEvent(cl_event event, cl_int status, Device & d) {
	if (status != CL_COMPLETE) {
		std::cout << "Dispatcher::onEvent - Got bad status: " << status << std::endl;
	}
	else if (d.m_eventFinished != NULL) {
		initContinue(d);
	} else {
		++d.m_round;
		handleResult(d);

		bool bDispatch = true;
		{
			std::lock_guard<std::mutex> lock(m_mutex);
			m_totalSearched += m_size;  // 累加已搜索数量
			d.m_speed.sample(m_size);
			printSpeed();

			if( m_quit ) {
				bDispatch = false;
				if(--m_countRunning == 0) {
					clSetUserEventStatus(m_eventFinished, CL_COMPLETE);
				}
			}
		}

		if (bDispatch) {
			dispatch(d);
		}
	}
}

// This is run when m_mutex is held.
void Dispatcher::printSpeed() {
	++m_countPrint;
	if (m_countPrint > m_vDevices.size()) {
		// 计算运行时间
		const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - timeStart).count();
		const int hours = elapsed / 3600;
		const int minutes = (elapsed % 3600) / 60;
		const int seconds = elapsed % 60;

		// 计算总速度
		double speedTotal = 0;
		std::string strGPUs;
		for (auto & e : m_vDevices) {
			const auto curSpeed = e->m_speed.getSpeed();
			speedTotal += curSpeed;
			strGPUs += " GPU" + toString(e->m_index) + ":" + formatSpeed(curSpeed);
		}

		// 格式化已搜索数量
		std::string strSearched;
		if (m_totalSearched >= 1000000000000ULL) {
			strSearched = toString(m_totalSearched / 1000000000000ULL) + "T";
		} else if (m_totalSearched >= 1000000000ULL) {
			strSearched = toString(m_totalSearched / 1000000000ULL) + "G";
		} else if (m_totalSearched >= 1000000ULL) {
			strSearched = toString(m_totalSearched / 1000000ULL) + "M";
		} else if (m_totalSearched >= 1000ULL) {
			strSearched = toString(m_totalSearched / 1000ULL) + "K";
		} else {
			strSearched = toString(m_totalSearched);
		}

		const std::string strVT100ClearLine = "\33[2K\r";
		std::cerr << strVT100ClearLine
		          << "[" << std::setfill('0') << std::setw(2) << hours << ":"
		          << std::setw(2) << minutes << ":" << std::setw(2) << seconds << "] "
		          << "速度:" << formatSpeed(speedTotal)
		          << " | 已搜索:" << strSearched
		          << " | 已找到:" << m_totalFound
		          << " |" << strGPUs
		          << '\r' << std::flush;
		m_countPrint = 0;
	}
}

void CL_CALLBACK Dispatcher::staticCallback(cl_event event, cl_int event_command_exec_status, void * user_data) {
	Device * const pDevice = static_cast<Device *>(user_data);
	pDevice->m_parent.onEvent(event, event_command_exec_status, *pDevice);
	clReleaseEvent(event);
}

std::string Dispatcher::formatSpeed(double f) {
	const std::string S = " KMGT";

	unsigned int index = 0;
	while (f > 1000.0f && index < S.size()) {
		f /= 1000.0f;
		++index;
	}

	std::ostringstream ss;
	ss << std::fixed << std::setprecision(3) << (double)f << " " << S[index] << "H/s";
	return ss.str();
}

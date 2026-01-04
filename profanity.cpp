#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <thread>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#include <OpenCL/cl_ext.h>
#else
#include <CL/cl.h>
#include <CL/cl_ext.h>
#endif

#define CL_DEVICE_PCI_BUS_ID_NV  0x4008
#define CL_DEVICE_PCI_SLOT_ID_NV 0x4009

// AMD OpenCL extension definitions (may not be available on all platforms)
#ifndef CL_DEVICE_TOPOLOGY_AMD
#define CL_DEVICE_TOPOLOGY_AMD 0x4037
#endif

#ifndef CL_DEVICE_TOPOLOGY_TYPE_PCIE_AMD
#define CL_DEVICE_TOPOLOGY_TYPE_PCIE_AMD 1

typedef union {
    struct { cl_uint type; cl_uint data[5]; } raw;
    struct { cl_uint type; cl_char unused[17]; cl_char bus; cl_char device; cl_char function; } pcie;
} cl_device_topology_amd;

#endif

#include "Dispatcher.hpp"
#include "ArgParser.hpp"
#include "Mode.hpp"
#include "help.hpp"
#include "KeyGenerator.hpp"

static std::string readFile(const char * const szFilename) {
	std::ifstream in(szFilename, std::ios::in | std::ios::binary);
	std::ostringstream contents;
	contents << in.rdbuf();
	return contents.str();
}

static std::vector<cl_device_id> getAllDevices(cl_device_type deviceType) {
	std::vector<cl_device_id> vDevices;

	cl_uint platformIdCount = 0;
	clGetPlatformIDs(0, NULL, &platformIdCount);
	if (platformIdCount == 0) return vDevices;

	std::vector<cl_platform_id> platformIds(platformIdCount);
	clGetPlatformIDs(platformIdCount, platformIds.data(), NULL);

	for (const auto & platformId : platformIds) {
		cl_uint countDevice = 0;
		cl_int err = clGetDeviceIDs(platformId, deviceType, 0, NULL, &countDevice);
		if (err != CL_SUCCESS || countDevice == 0) continue;

		std::vector<cl_device_id> deviceIds(countDevice);
		err = clGetDeviceIDs(platformId, deviceType, countDevice, deviceIds.data(), &countDevice);
		if (err == CL_SUCCESS) {
			std::copy(deviceIds.begin(), deviceIds.end(), std::back_inserter(vDevices));
		}
	}

	return vDevices;
}

template <typename T, typename U, typename V, typename W>
static T clGetWrapper(U function, V param, W param2) {
	T t;
	function(param, param2, sizeof(t), &t, NULL);
	return t;
}

template <typename U, typename V, typename W>
static std::string clGetWrapperString(U function, V param, W param2) {
	size_t len;
	function(param, param2, 0, NULL, &len);
	char * const szString = new char[len];
	function(param, param2, len, szString, NULL);
	std::string r(szString);
	delete[] szString;
	return r;
}

template <typename T, typename U, typename V, typename W>
static std::vector<T> clGetWrapperVector(U function, V param, W param2) {
	size_t len;
	function(param, param2, 0, NULL, &len);
	len /= sizeof(T);
	std::vector<T> v;
	if (len > 0) {
		T * pArray = new T[len];
		function(param, param2, len * sizeof(T), pArray, NULL);
		for (size_t i = 0; i < len; ++i) {
			v.push_back(pArray[i]);
		}
		delete[] pArray;
	}
	return v;
}

static std::vector<std::string> getBinaries(cl_program & clProgram) {
	std::vector<std::string> vReturn;
	auto vSizes = clGetWrapperVector<size_t>(clGetProgramInfo, clProgram, CL_PROGRAM_BINARY_SIZES);
	if (!vSizes.empty()) {
		unsigned char ** pBuffers = new unsigned char *[vSizes.size()];
		for (size_t i = 0; i < vSizes.size(); ++i) {
			pBuffers[i] = new unsigned char[vSizes[i]];
		}

		clGetProgramInfo(clProgram, CL_PROGRAM_BINARIES, vSizes.size() * sizeof(unsigned char *), pBuffers, NULL);
		for (size_t i = 0; i < vSizes.size(); ++i) {
			std::string strData(reinterpret_cast<char *>(pBuffers[i]), vSizes[i]);
			vReturn.push_back(strData);
			delete[] pBuffers[i];
		}

		delete[] pBuffers;
	}

	return vReturn;
}

static unsigned int getUniqueDeviceIdentifier(const cl_device_id & deviceId) {
	auto topology = clGetWrapper<cl_device_topology_amd>(clGetDeviceInfo, deviceId, CL_DEVICE_TOPOLOGY_AMD);
	if (topology.raw.type == CL_DEVICE_TOPOLOGY_TYPE_PCIE_AMD) {
		return (topology.pcie.bus << 16) + (topology.pcie.device << 8) + topology.pcie.function;
	}
	cl_int bus_id = clGetWrapper<cl_int>(clGetDeviceInfo, deviceId, CL_DEVICE_PCI_BUS_ID_NV);
	cl_int slot_id = clGetWrapper<cl_int>(clGetDeviceInfo, deviceId, CL_DEVICE_PCI_SLOT_ID_NV);
	return (bus_id << 16) + slot_id;
}

template <typename T>
static bool printResult(const T & t, const cl_int & err) {
	std::cout << ((t == NULL) ? toString(err) : "OK") << std::endl;
	return t == NULL;
}

static bool printResult(const cl_int err) {
	std::cout << ((err != CL_SUCCESS) ? toString(err) : "OK") << std::endl;
	return err != CL_SUCCESS;
}

static std::string getDeviceCacheFilename(cl_device_id & d, const size_t & inverseSize, bool isCPU) {
	const auto uniqueId = getUniqueDeviceIdentifier(d);
	return "cache-opencl." + toString(inverseSize) + "." + toString(uniqueId) + (isCPU ? ".cpu" : "");
}

// 根据设备能力自动优化参数
struct DeviceOptimization {
	size_t worksizeLocal;
	size_t inverseMultiple;
	cl_ulong maxMemory;
	cl_uint computeUnits;
	bool isGPU;
};

static DeviceOptimization optimizeForDevice(cl_device_id deviceId, bool isGPU, size_t gpuMemPercent, size_t cpuCores) {
	DeviceOptimization opt;
	opt.isGPU = isGPU;

	// 获取设备信息
	opt.computeUnits = clGetWrapper<cl_uint>(clGetDeviceInfo, deviceId, CL_DEVICE_MAX_COMPUTE_UNITS);
	opt.maxMemory = clGetWrapper<cl_ulong>(clGetDeviceInfo, deviceId, CL_DEVICE_GLOBAL_MEM_SIZE);
	size_t maxWorkGroupSize = clGetWrapper<size_t>(clGetDeviceInfo, deviceId, CL_DEVICE_MAX_WORK_GROUP_SIZE);

	if (isGPU) {
		// GPU 优化：最大化并行度
		opt.worksizeLocal = std::min(maxWorkGroupSize, (size_t)256);
		// 根据显存限制调整
		cl_ulong effectiveMemory = opt.maxMemory * gpuMemPercent / 100;
		size_t maxItems = effectiveMemory / (1024 * 4);
		opt.inverseMultiple = std::min(maxItems / 255, (size_t)65536);
		opt.inverseMultiple = std::max(opt.inverseMultiple, (size_t)4096);
	} else {
		// CPU 优化：根据核心数调整
		opt.worksizeLocal = 1;
		// 如果指定了核心数限制
		cl_uint effectiveCores = (cpuCores > 0 && cpuCores < opt.computeUnits) ? cpuCores : opt.computeUnits;
		opt.inverseMultiple = effectiveCores * 256;
		opt.inverseMultiple = std::min(opt.inverseMultiple, (size_t)8192);
		opt.inverseMultiple = std::max(opt.inverseMultiple, (size_t)512);
		// 更新有效核心数
		opt.computeUnits = effectiveCores;
	}

	return opt;
}

int main(int argc, char ** argv) {
	try {
		ArgParser argp(argc, argv);

		// 基本参数
		bool bHelp = false;
		bool bBenchmark = false;
		std::string strPublicKey;

		// TRON 模式参数
		bool bModeTronRepeat = false;
		bool bModeTronSequential = false;
		std::string strModeTronSuffix;
		bool bModeTronLucky = false;
		size_t repeatCount = 4;     // 豹子号最少位数，默认4
		size_t sequentialCount = 4; // 顺子号最少位数，默认4

		// 设备和性能参数
		size_t worksizeLocal = 0;  // 0 表示自动检测
		size_t worksizeMax = 0;
		bool bNoCache = false;
		size_t inverseSize = 255;
		size_t inverseMultiple = 0;  // 0 表示自动检测

		// 设备控制参数
		std::string strDeviceType;  // gpu, cpu
		size_t cpuCores = 0;        // 0 表示使用全部核心
		size_t gpuMemPercent = 100; // GPU显存使用百分比

		// 注册参数
		argp.addSwitch('h', "help", bHelp);
		argp.addSwitch('0', "benchmark", bBenchmark);
		argp.addSwitch('z', "publicKey", strPublicKey);

		// TRON 模式
		argp.addSwitch('R', "tron-repeat", bModeTronRepeat);
		argp.addSwitch('S', "tron-sequential", bModeTronSequential);
		argp.addSwitch('T', "tron-suffix", strModeTronSuffix);
		argp.addSwitch('L', "tron-lucky", bModeTronLucky);
		argp.addSwitch('r', "repeat-count", repeatCount);
		argp.addSwitch('s', "sequential-count", sequentialCount);

		// 设备控制
		argp.addSwitch('n', "no-cache", bNoCache);
		argp.addSwitch('d', "device", strDeviceType);
		argp.addSwitch('c', "cpu-cores", cpuCores);
		argp.addSwitch('g', "gpu-mem", gpuMemPercent);

		// 性能调优
		argp.addSwitch('w', "work", worksizeLocal);
		argp.addSwitch('W', "work-max", worksizeMax);
		argp.addSwitch('i', "inverse-size", inverseSize);
		argp.addSwitch('I', "inverse-multiple", inverseMultiple);

		if (!argp.parse()) {
			std::cout << "错误: 参数解析失败" << std::endl;
			return 1;
		}

		if (bHelp) {
			std::cout << g_strHelp << std::endl;
			return 0;
		}

		// 验证参数
		if (gpuMemPercent < 1 || gpuMemPercent > 100) {
			std::cout << "错误: GPU显存百分比必须在 1-100 之间" << std::endl;
			return 1;
		}

		// 选择模式
		Mode mode = Mode::benchmark();
		if (bBenchmark) {
			mode = Mode::benchmark();
		} else if (bModeTronRepeat) {
			mode = Mode::tronRepeat(repeatCount);
		} else if (bModeTronSequential) {
			mode = Mode::tronSequential(sequentialCount);
		} else if (!strModeTronSuffix.empty()) {
			mode = Mode::tronSuffix(strModeTronSuffix);
		} else if (bModeTronLucky) {
			mode = Mode::tronLucky();
		} else {
			std::cout << g_strHelp << std::endl;
			return 0;
		}
		
		// 自动生成密钥对
		std::string generatedPrivateKey;
		if (strPublicKey.empty()) {
			std::cout << "自动生成密钥对..." << std::endl;
			KeyGenerator keyGen;
			keyGen.generate();
			strPublicKey = keyGen.publicKey;
			generatedPrivateKey = keyGen.privateKey;
			// 加密显示种子私钥
			std::string maskedKey = generatedPrivateKey.substr(0, 6) + std::string(52, '*') + generatedPrivateKey.substr(58, 6);
			std::cout << "种子私钥: 0x" << maskedKey << std::endl;
			std::cout << std::endl;
		}

		if (strPublicKey.length() != 128) {
			std::cout << "错误: 公钥必须是128位十六进制字符" << std::endl;
			return 1;
		}

		std::cout << "模式: " << mode.name << std::endl;

		// 根据 --device 参数确定设备类型
		bool forceGPU = (strDeviceType == "gpu");
		bool forceCPU = (strDeviceType == "cpu");

		std::vector<cl_device_id> vFoundDevices;

		if (forceCPU) {
			// 强制使用 CPU
			vFoundDevices = getAllDevices(CL_DEVICE_TYPE_CPU);
			if (vFoundDevices.empty()) {
				std::cout << "错误: 未找到可用的 CPU OpenCL 设备" << std::endl;
				return 1;
			}
		} else if (forceGPU) {
			// 强制使用 GPU
			vFoundDevices = getAllDevices(CL_DEVICE_TYPE_GPU);
			if (vFoundDevices.empty()) {
				std::cout << "错误: 未找到可用的 GPU 设备" << std::endl;
				return 1;
			}
		} else {
			// 自动检测：优先 GPU，没有则用 CPU
			vFoundDevices = getAllDevices(CL_DEVICE_TYPE_GPU);
			if (vFoundDevices.empty()) {
				std::cout << "未检测到 GPU，自动使用 CPU..." << std::endl;
				vFoundDevices = getAllDevices(CL_DEVICE_TYPE_CPU);
			}
			if (vFoundDevices.empty()) {
				std::cout << "错误: 未找到任何可用的 OpenCL 设备" << std::endl;
				return 1;
			}
		}

		std::vector<cl_device_id> vDevices;
		std::map<cl_device_id, size_t> mDeviceIndex;
		std::vector<std::string> vDeviceBinary;
		std::vector<size_t> vDeviceBinarySize;
		cl_int errorCode;
		bool bUsedCache = false;

		// 用于存储优化参数
		DeviceOptimization bestOpt = {64, 16384, 0, 0, true};

		// 显示资源限制信息
		if (cpuCores > 0) {
			std::cout << "CPU核心限制: " << cpuCores << std::endl;
		}
		if (gpuMemPercent < 100) {
			std::cout << "GPU显存限制: " << gpuMemPercent << "%" << std::endl;
		}

		std::cout << "设备:" << std::endl;
		for (size_t i = 0; i < vFoundDevices.size(); ++i) {
			cl_device_id & deviceId = vFoundDevices[i];
			cl_device_type deviceType = clGetWrapper<cl_device_type>(clGetDeviceInfo, deviceId, CL_DEVICE_TYPE);
			bool isGPU = (deviceType == CL_DEVICE_TYPE_GPU);

			const auto strName = clGetWrapperString(clGetDeviceInfo, deviceId, CL_DEVICE_NAME);
			const auto computeUnits = clGetWrapper<cl_uint>(clGetDeviceInfo, deviceId, CL_DEVICE_MAX_COMPUTE_UNITS);
			const auto globalMemSize = clGetWrapper<cl_ulong>(clGetDeviceInfo, deviceId, CL_DEVICE_GLOBAL_MEM_SIZE);
			const auto maxClockFreq = clGetWrapper<cl_uint>(clGetDeviceInfo, deviceId, CL_DEVICE_MAX_CLOCK_FREQUENCY);
			bool precompiled = false;

			// 自动优化参数（考虑资源限制）
			DeviceOptimization opt = optimizeForDevice(deviceId, isGPU, gpuMemPercent, cpuCores);
			if (opt.computeUnits > bestOpt.computeUnits || (isGPU && !bestOpt.isGPU)) {
				bestOpt = opt;
			}

			if (!bNoCache) {
				std::ifstream fileIn(getDeviceCacheFilename(deviceId, inverseSize, !isGPU), std::ios::binary);
				if (fileIn.is_open()) {
					vDeviceBinary.push_back(std::string((std::istreambuf_iterator<char>(fileIn)), std::istreambuf_iterator<char>()));
					vDeviceBinarySize.push_back(vDeviceBinary.back().size());
					precompiled = true;
				}
			}

			std::string devTypeName = isGPU ? "GPU" : "CPU";
			std::cout << "  " << devTypeName << i << ": " << strName << std::endl;

			// 显示有效资源（考虑限制）
			cl_uint effectiveCores = isGPU ? computeUnits : opt.computeUnits;
			cl_ulong effectiveMem = isGPU ? (globalMemSize * gpuMemPercent / 100) : globalMemSize;

			std::cout << "      内存: " << (effectiveMem / 1024 / 1024) << " MB";
			if (isGPU && gpuMemPercent < 100) {
				std::cout << " (限制" << gpuMemPercent << "%)";
			}
			std::cout << ", 计算单元: " << effectiveCores;
			if (!isGPU && cpuCores > 0 && cpuCores < computeUnits) {
				std::cout << "/" << computeUnits << " (限制)";
			}
			std::cout << ", 频率: " << maxClockFreq << " MHz"
			          << (precompiled ? " [cached]" : "") << std::endl;

			vDevices.push_back(vFoundDevices[i]);
			mDeviceIndex[vFoundDevices[i]] = i;
		}

		if (vDevices.empty()) {
			std::cout << "错误: 未找到可用的计算设备" << std::endl;
			return 1;
		}

		// 如果用户没有指定，使用自动优化的参数
		if (worksizeLocal == 0) {
			worksizeLocal = bestOpt.worksizeLocal;
		}
		if (inverseMultiple == 0) {
			inverseMultiple = bestOpt.inverseMultiple;
		}

		std::cout << std::endl;
		std::cout << "优化参数: 工作组大小=" << worksizeLocal
		          << ", 并行度=" << inverseMultiple
		          << ", 总工作项=" << (inverseSize * inverseMultiple) << std::endl;

		std::cout << std::endl;
		std::cout << "初始化 OpenCL..." << std::endl;
		std::cout << "  创建上下文..." << std::flush;
		auto clContext = clCreateContext(NULL, vDevices.size(), vDevices.data(), NULL, NULL, &errorCode);
		if (printResult(clContext, errorCode)) {
			return 1;
		}

		cl_program clProgram;
		if (vDeviceBinary.size() == vDevices.size()) {
			bUsedCache = true;
			std::cout << "  加载缓存内核..." << std::flush;
			const unsigned char ** pKernels = new const unsigned char *[vDevices.size()];
			for (size_t i = 0; i < vDeviceBinary.size(); ++i) {
				pKernels[i] = reinterpret_cast<const unsigned char *>(vDeviceBinary[i].data());
			}
			cl_int * pStatus = new cl_int[vDevices.size()];
			clProgram = clCreateProgramWithBinary(clContext, vDevices.size(), vDevices.data(), vDeviceBinarySize.data(), pKernels, pStatus, &errorCode);
			delete[] pKernels;
			delete[] pStatus;
			if (printResult(clProgram, errorCode)) {
				return 1;
			}
		} else {
			std::cout << "  编译内核..." << std::flush;
			const std::string strKeccak = readFile("keccak.cl");
			const std::string strVanity = readFile("profanity.cl");
			const char * szKernels[] = { strKeccak.c_str(), strVanity.c_str() };

			clProgram = clCreateProgramWithSource(clContext, sizeof(szKernels) / sizeof(char *), szKernels, NULL, &errorCode);
			if (printResult(clProgram, errorCode)) {
				return 1;
			}
		}

		std::cout << "  构建程序..." << std::flush;
		const std::string strBuildOptions = "-D PROFANITY_INVERSE_SIZE=" + toString(inverseSize) + " -D PROFANITY_MAX_SCORE=" + toString(PROFANITY_MAX_SCORE);
		if (printResult(clBuildProgram(clProgram, vDevices.size(), vDevices.data(), strBuildOptions.c_str(), NULL, NULL))) {
			return 1;
		}

		if (!bUsedCache && !bNoCache) {
			std::cout << "  保存缓存..." << std::flush;
			auto binaries = getBinaries(clProgram);
			for (size_t i = 0; i < binaries.size(); ++i) {
				cl_device_type devType = clGetWrapper<cl_device_type>(clGetDeviceInfo, vDevices[i], CL_DEVICE_TYPE);
				bool isCPU = (devType == CL_DEVICE_TYPE_CPU);
				std::ofstream fileOut(getDeviceCacheFilename(vDevices[i], inverseSize, isCPU), std::ios::binary);
				fileOut.write(binaries[i].data(), binaries[i].size());
			}
			std::cout << "OK" << std::endl;
		}

		std::cout << std::endl;

		Dispatcher d(clContext, clProgram, mode, worksizeMax == 0 ? inverseSize * inverseMultiple : worksizeMax, inverseSize, inverseMultiple, 0, strPublicKey, generatedPrivateKey);
		for (auto & i : vDevices) {
			d.addDevice(i, worksizeLocal, mDeviceIndex[i]);
		}

		d.run();
		clReleaseContext(clContext);
		return 0;

	} catch (std::runtime_error & e) {
		std::cout << "运行时错误: " << e.what() << std::endl;
	} catch (...) {
		std::cout << "未知错误" << std::endl;
	}

	return 1;
}


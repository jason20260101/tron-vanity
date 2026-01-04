#include "Mode.hpp"
#include <algorithm>

Mode::Mode() : score(0) {
	std::fill(data1, data1 + sizeof(data1), cl_uchar(0));
	std::fill(data2, data2 + sizeof(data2), cl_uchar(0));
}

Mode Mode::benchmark() {
	Mode r;
	r.name = "benchmark";
	r.kernel = "profanity_score_benchmark";
	return r;
}

// TRON vanity address modes

Mode Mode::tronRepeat() {
	Mode r;
	r.name = "tron-repeat (豹子号)";
	r.kernel = "profanity_score_tron_repeat";
	return r;
}

Mode Mode::tronSequential() {
	Mode r;
	r.name = "tron-sequential (顺子号)";
	r.kernel = "profanity_score_tron_sequential";
	return r;
}

Mode Mode::tronSuffix(const std::string suffix) {
	Mode r;
	r.name = "tron-suffix (自定义后缀)";
	r.kernel = "profanity_score_tron_suffix";

	// Parse comma-separated patterns
	// Store patterns in data1, separated by 0x00
	// data2[0] = total length including separators
	// data2[1] = number of patterns

	size_t dataPos = 0;
	size_t patternCount = 0;
	size_t i = 0;

	while (i < suffix.length() && dataPos < 19) {
		// Skip leading commas
		while (i < suffix.length() && suffix[i] == ',') {
			i++;
		}

		if (i >= suffix.length()) break;

		// Find pattern end
		size_t patternStart = i;
		while (i < suffix.length() && suffix[i] != ',') {
			i++;
		}
		size_t patternLen = i - patternStart;

		if (patternLen > 0) {
			// Check if pattern fits
			if (dataPos + patternLen + 1 > 20) {
				patternLen = 20 - dataPos - 1;
				if (patternLen <= 0) break;
			}

			// Copy pattern
			for (size_t j = 0; j < patternLen; ++j) {
				r.data1[dataPos++] = static_cast<cl_uchar>(suffix[patternStart + j]);
			}

			// Add null separator
			r.data1[dataPos++] = 0;
			patternCount++;
		}
	}

	r.data2[0] = static_cast<cl_uchar>(dataPos);
	r.data2[1] = static_cast<cl_uchar>(patternCount);

	return r;
}

Mode Mode::tronLucky() {
	Mode r;
	r.name = "tron-lucky (谐音靓号)";
	r.kernel = "profanity_score_tron_lucky";
	return r;
}

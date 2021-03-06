#include "CommonApprox.h"
#include <vector>
#include <algorithm>

const floattype dfYOffset = 10000.0;

CCommonApprox::CCommonApprox(IGlucoseLevels *levels) : mEnumeratedLevels(levels) {
	if (mEnumeratedLevels != NULL) {
		mEnumeratedLevels->AddRef();
		mEnumeratedLevels->GetLevels(&this->levels);
		mEnumeratedLevels->GetLevelsCount(&size);
	}
}

CCommonApprox::~CCommonApprox() {
	if (mEnumeratedLevels != NULL) mEnumeratedLevels->Release();
}

HRESULT get_time_interval(TGlucoseLevel *levels, size_t size, floattype time, int *index) {
	std::vector<TGlucoseLevel> vec(levels, levels + size);
	std::vector<TGlucoseLevel>::const_reverse_iterator cri = std::lower_bound(vec.rbegin(),
		vec.rend(), time, [&](const TGlucoseLevel &gl, const floattype &time) { return gl.datetime > time; });
	if (cri != vec.rend()) {
		*index = static_cast<int>(vec.rend() - cri - 1); // retrieve index of the found glucoselevel
		return S_OK;
	}
	return S_FALSE;
}

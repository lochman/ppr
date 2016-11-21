#pragma once

#include "CommonApprox.h"
#include<vector>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class CubicApprox : public CCommonApprox {
	std::vector<TGlucoseLevel *> masks[255];
	std::vector<floattype> a[255], b[255], c[255], d[255];
	size_t sizes[255];
public:
	CubicApprox(IGlucoseLevels *levels) : CCommonApprox(levels) { };
	virtual HRESULT Approximate(TApproximationParams * params);
	virtual HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count,
		floattype *levels, size_t *filled, size_t derivationorder);
};

#pragma warning( pop )
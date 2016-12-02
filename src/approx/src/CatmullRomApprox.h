#pragma once

#include "CommonApprox.h"
#include <vector>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class CatmullRomApprox : public CCommonApprox {
	std::vector<TGlucoseLevel> c, d;
	floattype iterate(TGlucoseLevel &p0, TGlucoseLevel &p1, TGlucoseLevel &p2, TGlucoseLevel &p3, floattype n, size_t *i);
public:
	CatmullRomApprox(IGlucoseLevels *levels) : CCommonApprox(levels) { };
	virtual HRESULT Approximate(TApproximationParams * params);
	virtual HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count,
		floattype *levels, size_t *filled, size_t derivationorder);
};

#pragma warning( pop )
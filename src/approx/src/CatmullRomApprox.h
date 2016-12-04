#pragma once

#include "CommonApprox.h"
#include <vector>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class CatmullRomApprox : public CCommonApprox {
	std::vector<floattype> a, b, c, d;
	void get_tangent(const floattype &p0, const floattype &p1, const floattype &p2,
		const floattype &p3, const floattype &t0, const floattype &t1, const floattype &t2);
	void get_coefficients(const floattype &p0, const floattype &p1, const floattype t0, const floattype t1);
	void extrapolation();
	void iterate(const TGlucoseLevel &p0, const TGlucoseLevel &p1, const TGlucoseLevel &p2, const TGlucoseLevel &p3);
	void approximate_gpu();
public:
	CatmullRomApprox(IGlucoseLevels *levels) : CCommonApprox(levels) { };
	virtual HRESULT Approximate(TApproximationParams * params);
	virtual HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count,
		floattype *levels, size_t *filled, size_t derivationorder);
};

#pragma warning( pop )
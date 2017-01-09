#pragma once

#include "CommonApprox.h"
#include <vector>
#include <deque>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class AkimaApprox : public CCommonApprox {
	std::vector<floattype> p1, p2, p3;
	std::deque<floattype> m;
	floattype get_m(size_t i);
	HRESULT iterate(floattype &m_next, size_t i, floattype *ti);
	void approximate_gpu(floattype *ti);
public:
	AkimaApprox(IGlucoseLevels *levels) : CCommonApprox(levels) { };
	virtual HRESULT Approximate(TApproximationParams * params);
	virtual HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count,
		floattype *levels, size_t *filled, size_t derivationorder);
	virtual std::string get_name();
};

#pragma warning( pop )
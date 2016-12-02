#pragma once

#include "CommonApprox.h"
#include <vector>
#include <deque>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class AkimaApprox : public CCommonApprox {
	std::vector<floattype> p1, p2, p3;
	std::deque<floattype> m;
	//HRESULT get_m(int i, floattype *m);
	floattype get_m(int i);
	//floattype get_t(floattype m0, floattype m1, floattype m2, floattype m3) restrict(cpu, amp);
	//HRESULT get_t(floattype *t);
	HRESULT iterate(floattype &m_next, int i, floattype *ti);
public:
	AkimaApprox(IGlucoseLevels *levels) : CCommonApprox(levels) { };
	virtual HRESULT Approximate(TApproximationParams * params);
	virtual HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count,
		floattype *levels, size_t *filled, size_t derivationorder);
};

#pragma warning( pop )
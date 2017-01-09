#pragma once

#include "common\rtl\referencedImpl.h"
#include "common/iface/ApproxIface.h"
#include "approx\src\GlucoseLevels.h"
#include "approx\src\CommonApprox.h"
#include "MaskService.h"
#include <vector>
#include <chrono>
#include <map>
#include <sstream>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class Statistics : public virtual CReferenced {
	MaskService *mask_service;
	TGlucoseLevel *ref_lvls;
	std::map<floattype, floattype> *ref_devs;
	std::stringstream output;
	boolean graph;
	size_t size;
public:
	Statistics(MaskService *mask_service, std::map<floattype, floattype> *ref_devs, boolean graph);
	HRESULT get_stats(const int mask, CCommonApprox *approx);
	HRESULT get_errors(TGlucoseLevel *levels, size_t size, const int mask, CCommonApprox *approx, boolean graph);
	HRESULT print_stats(std::vector<floattype> &errors);
	std::string get_output();
};

class Timer : public virtual CReferenced {
	std::chrono::steady_clock::time_point begin, end;
	std::string message;
public:
	Timer(std::string message) : message(message) { }
	HRESULT start();
	HRESULT stop();
};

#pragma warning( pop )
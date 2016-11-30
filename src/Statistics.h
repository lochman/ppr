#pragma once

#include "common\rtl\referencedImpl.h"
#include "common/iface/ApproxIface.h"
#include "approx\src\GlucoseLevels.h"
#include "approx\src\CommonApprox.h"
#include "MaskService.h"
#include <vector>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class Statistics : public virtual CReferenced {
	TGlucoseLevel *ref_values;
	size_t size;
public:
	Statistics(MaskService *mask_service, int &mask, CCommonApprox *approx);
	void get_errors(TGlucoseLevel *ref_values, int &mask, CCommonApprox *approx);
	void print_stats(std::vector<floattype> &errors);
};

#pragma warning( pop )
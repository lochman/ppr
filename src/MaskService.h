#pragma once

#include "common\rtl\referencedImpl.h"
#include "common/iface/ApproxIface.h"
#include <vector>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class MaskService : public virtual CReferenced {
	TGlucoseLevel *levels;
	size_t size;
public:
	MaskService(IGlucoseLevels *mEnumeratedLevels);
	void get_masked_values(std::vector<TGlucoseLevel *> *glucose_levels, uint8_t mask);
	void get_inverse_masked_values(std::vector<TGlucoseLevel *> *glucose_levels, uint8_t mask);
};

#pragma warning( pop )
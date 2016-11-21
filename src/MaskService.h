#pragma once

#include "common\rtl\referencedImpl.h"
#include "common/iface/ApproxIface.h"

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class MaskService : public virtual CReferenced {
public:
	void get_masked_values(IGlucoseLevels *masked, IGlucoseLevels &mEnumeratedLevels, uint8_t mask);
	void get_inverse_masked_values(IGlucoseLevels *masked, IGlucoseLevels &mEnumeratedLevels, uint8_t mask);
};

#pragma warning( pop )
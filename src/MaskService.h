#pragma once

#include "common\rtl\referencedImpl.h"
#include "common/iface/ApproxIface.h"
#include "approx\src\GlucoseLevels.h"
#include <vector>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class MaskService : public virtual CReferenced {
	std::vector<IGlucoseLevels *> levels;
	size_t segment_size;
	HRESULT get_masked_values(const TGlucoseLevel *levels, uint8_t mask);
public:
	MaskService(IGlucoseLevels *levels);
	~MaskService();
	HRESULT get_mask(IGlucoseLevels **levels, uint8_t mask);
	HRESULT get_inverse_mask(IGlucoseLevels **levels, uint8_t mask);
	HRESULT get_segment_size(size_t *size);
	HRESULT get_mask_count(size_t *size);
};

#pragma warning( pop )
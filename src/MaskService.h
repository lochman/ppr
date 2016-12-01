#pragma once

#include "common\rtl\referencedImpl.h"
#include "common/iface/ApproxIface.h"
#include "approx\src\GlucoseLevels.h"
#include <vector>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class MaskService : public virtual CReferenced {
	CGlucoseLevels masks[255];
	TGlucoseLevel *levels;
	size_t size;
public:
	MaskService(TGlucoseLevel *levels, size_t const &size);
	void get_masked_values(std::vector<TGlucoseLevel> &glucose_levels, uint8_t mask) const;
	void get_mask(IGlucoseLevels **levels, uint8_t mask);
	void get_inverse_mask(IGlucoseLevels **levels, uint8_t mask);
	void get_levels(TGlucoseLevel **levels) const;
	void get_levels_size(size_t *size) const;
};

#pragma warning( pop )
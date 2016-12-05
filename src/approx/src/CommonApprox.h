#pragma once

#include "../../common/iface/ApproxIface.h"
#include "..\..\common\rtl\hresult.h"
#include "..\..\common\rtl\referencedImpl.h"
#include <string>


extern const floattype dfYOffset; //some interpolation requires negative values
				//and it is impossible to compute real ln of a negative number
				//100.0 and lower do not work always

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class CCommonApprox : public IApproximatedGlucoseLevels, public virtual CReferenced {
protected:
	IGlucoseLevels *mEnumeratedLevels;
	TGlucoseLevel *levels;
	size_t size;
public:
	CCommonApprox(IGlucoseLevels *levels);
	virtual ~CCommonApprox();
	//dctor has to be virtual, even if it is empty, due to the inheritance by dominance	
};

HRESULT get_time_interval(TGlucoseLevel *levels, size_t size, floattype time, int *index);

#pragma warning( pop )
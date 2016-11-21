#pragma once

#include <vector>

#include "common\rtl\hresult.h"
#include"approx/src/GlucoseLevels.h"

class DataService {
public:
	virtual ~DataService() { };
	virtual HRESULT get_segments(std::vector<std::shared_ptr<CGlucoseLevels>> &segments) = 0;
};
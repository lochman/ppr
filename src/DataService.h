#pragma once

#include <vector>

#include "common\rtl\hresult.h"
#include"approx/src/GlucoseLevels.h"

class DataService {
protected:
	std::string path;
public:
	DataService(std::string path) : path(path) { };
	virtual ~DataService() { };
	virtual HRESULT get_segments(std::vector<IGlucoseLevels *> &segments, std::vector<std::string> &segment_ids) = 0;
};
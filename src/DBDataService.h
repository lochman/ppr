#pragma once

#include <vector>

#include "common\rtl\referencedImpl.h"
#include "approx/src/GlucoseLevels.h"
#include "sqlite3/sqlite3.h"
#include "DataService.h"

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class DBDataService : public DataService, public virtual CReferenced {
	sqlite3 *db;
	HRESULT get_glucose_levels(std::vector<TGlucoseLevel> &glucose_levels, std::string segmentid);
public:
	DBDataService(std::string path);
	~DBDataService();
	HRESULT get_segments(std::vector<IGlucoseLevels *> &segments, std::vector<std::string> &segment_ids);
};

#pragma warning( pop )
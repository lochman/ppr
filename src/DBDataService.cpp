#include "DBDataService.h"
#include <string>
#include <iostream>

double QDateTime2RatTime(const int64_t unixepochtime) {
	const int64_t diffFrom1970To1900 = 2208988800L; //2209161600000;
	const double SecsPerDay = 24.0 * 60.0 * 60.0; //*1000.0;
	const double InvSecsPerDay = 1.0 / SecsPerDay;

	int64_t diff = unixepochtime + diffFrom1970To1900;
	return ((double)diff) * InvSecsPerDay;
}

DBDataService::DBDataService(std::string path) : DataService(path) {
	if (sqlite3_open(path.c_str(), &db)) {
		std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
	} else {
		std::cerr << "Database " << path << " opened successfully" << std::endl;
	}
}

DBDataService::~DBDataService() {
	if (sqlite3_close(db)) {
		std::cerr << "Can't close database: " << sqlite3_errmsg(db) << std::endl;
	} else {
		std::cerr << "Database " << this->path << " closed successfully" << std::endl;
	}
}

HRESULT fill_lvl(IGlucoseLevels *level, std::vector<TGlucoseLevel> &gl_levels) {
	TGlucoseLevel *lvl;
	level->SetLevelsCount(gl_levels.size());
	level->GetLevels(&lvl);
	level->AddRef();
	memcpy(lvl, gl_levels.data(), gl_levels.size() * sizeof TGlucoseLevel);
	return S_OK;
}

HRESULT DBDataService::get_segments(std::vector<IGlucoseLevels *> &segments, std::vector<std::string> &segment_ids) {
	std::vector<TGlucoseLevel> glucose_levels;
	sqlite3_stmt *res;
	const char *tail;
	std::string query = "SELECT id FROM timesegment";
	IGlucoseLevels *lvls;

	if (sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.length()), &res, &tail) != SQLITE_OK) {
		std::cerr << sqlite3_errmsg(db) << std::endl;
		return S_FALSE;
	}

	while (sqlite3_step(res) == SQLITE_ROW) {
		//std::shared_ptr<CGlucoseLevels> segment = make_shared_reference(new CGlucoseLevels(), true);
		//std::vector<TGlucoseLevel> glucose_levels;
		std::string segmentid = std::string(reinterpret_cast<const char*>(sqlite3_column_text(res, 0)));
		if (get_glucose_levels(glucose_levels, segmentid) == S_OK && glucose_levels.size() > 0) {
			lvls = new CGlucoseLevels();
			fill_lvl(lvls, glucose_levels);
			segments.push_back(lvls);
			segment_ids.push_back(segmentid);
		}
		glucose_levels.clear();
	}

	sqlite3_finalize(res);
	return S_OK;
}

HRESULT DBDataService::get_glucose_levels(std::vector<TGlucoseLevel> &glucose_levels, std::string segmentid) {
	sqlite3_stmt *res;
	const char *tail;
	//std::string id = std::string(reinterpret_cast<const char*>(segmentid));
	std::string query = "SELECT strftime('%s', measuredat), ist FROM measuredvalue WHERE segmentid = " + segmentid + " AND ist IS NOT NULL";
	//std::vector<TGlucoseLevel> glucose_levels; //julianday(measuredat)

	if (sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.length()), &res, &tail) != SQLITE_OK) {
		std::cerr << sqlite3_errmsg(db) << std::endl;
		return S_FALSE;
	}

	while (sqlite3_step(res) == SQLITE_ROW) {
		TGlucoseLevel g_level;
		g_level.datetime = QDateTime2RatTime(std::stoll(reinterpret_cast<const char*>(sqlite3_column_text(res, 0))));
		//g_level.datetime = std::stod(reinterpret_cast<const char*>(sqlite3_column_text(res, 0)));
		g_level.level = std::stod(reinterpret_cast<const char*>(sqlite3_column_text(res, 1)));
		glucose_levels.push_back(g_level);
	}
	sqlite3_finalize(res);


	/*
	segment->SetLevelsCount(glucose_levels.size());
	TGlucoseLevel *lvl;
	segment->GetLevels(&lvl);
	memcpy(lvl, glucose_levels.data(), glucose_levels.size() * sizeof TGlucoseLevel);
	*/
	return S_OK;
}
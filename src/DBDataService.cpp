#include"DBDataService.h"
#include<string>

double QDateTime2RatTime(const int64_t unixepochtime) {
	const int64_t diffFrom1970To1900 = 2208988800L; //2209161600000;
	const double SecsPerDay = 24.0*60.0*60.0; //*1000.0;
	const double InvSecsPerDay = 1.0 / SecsPerDay;

	int64_t diff = unixepochtime + diffFrom1970To1900;
	return ((double)diff)*InvSecsPerDay;
}

DBDataService::DBDataService(const char *path) {
	if (sqlite3_open(path, &db)) {
		fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
	}
	else {
		fprintf(stderr, "Database opened successfully\n");
	}
}

DBDataService::~DBDataService() {
	if (sqlite3_close(db)) {
		fprintf(stderr, "Can't close database: %s\n", sqlite3_errmsg(db));
	}
	else {
		fprintf(stderr, "Database closed successfully\n");
	}
}

HRESULT DBDataService::get_segments(std::vector<std::vector<TGlucoseLevel>> &segments) {
	std::vector<TGlucoseLevel> glucose_levels;
	sqlite3_stmt *res;
	const char *tail;

	std::string query = "SELECT id FROM timesegment";

	if (sqlite3_prepare_v2(db, query.c_str(), query.length(), &res, &tail) != SQLITE_OK) {
		fprintf(stderr, "%s\n", sqlite3_errmsg(db));
		return S_FALSE;
	}

	while (sqlite3_step(res) == SQLITE_ROW) {
		//std::shared_ptr<CGlucoseLevels> segment = make_shared_reference(new CGlucoseLevels(), true);
		//std::vector<TGlucoseLevel> glucose_levels;
		if (get_glucose_levels(glucose_levels, sqlite3_column_text(res, 0)) == S_OK && glucose_levels.size() > 0) {
			segments.push_back(glucose_levels);
		}
		glucose_levels.clear();
	}

	sqlite3_finalize(res);
	return S_OK;
}

HRESULT DBDataService::get_glucose_levels(std::vector<TGlucoseLevel> &glucose_levels, const unsigned char *segmentid) {
	sqlite3_stmt *res;
	const char *tail;
	std::string id = std::string(reinterpret_cast<const char*>(segmentid));
	std::string query = "SELECT strftime('%s', measuredat), ist FROM measuredvalue WHERE segmentid = " + id + " AND ist IS NOT NULL";
	//std::vector<TGlucoseLevel> glucose_levels; //julianday(measuredat)

	if (sqlite3_prepare_v2(db, query.c_str(), query.length(), &res, &tail) != SQLITE_OK) {
		fprintf(stderr, "%s\n", sqlite3_errmsg(db));
		return S_FALSE;
	}

	while (sqlite3_step(res) == SQLITE_ROW) {
		TGlucoseLevel g_level;
		g_level.datetime = QDateTime2RatTime(std::stod(reinterpret_cast<const char*>(sqlite3_column_text(res, 0))));
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
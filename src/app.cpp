
#include "common/iface/ApproxIface.h"
#include "approx/src/GlucoseLevels.h"
#include <vector>
#include "DBDataService.h"
#include "approx\src\CubicApprox.h"
#include "approx\src\AkimaApprox.h"
#include "approx\src\CatmullRomApprox.h"
#include "MaskService.h"
#include "Statistics.h"
#include <tbb/tbb.h>
#include <chrono>
#include <iostream>
#include "ArgParser.h"
#include "defs.h"

HRESULT load_segments(const std::string &filename, std::vector<std::vector<TGlucoseLevel>> &segments) {
	DBDataService dbservice(filename.c_str());
	DataService *data_service = &dbservice;

	if (data_service->get_segments(segments) != S_OK) {
		printf("Failed to load segments from database.\n");
		return S_FALSE;
	}
	return S_OK;
}

HRESULT approx_mask(MaskService *mask_service, const std::string &method, int i) {
	CCommonApprox *approx;
	IGlucoseLevels *levels;
	mask_service->get_mask(&levels, i);
	if (str_compare(method, "akima") || str_compare(method, "a")) {
		approx = new AkimaApprox(levels);
	} else if (str_compare(method, "catmull") || str_compare(method, "cr")) {
		approx = new CatmullRomApprox(levels);
	} else if (str_compare(method, "cubic") || str_compare(method, "c")) {
		approx = new CubicApprox(levels);
	} else {
		approx = new AkimaApprox(levels);
	}
	approx->Approximate(nullptr);
	Statistics stats(mask_service, i, approx);
	delete approx;
	return S_OK;
}

HRESULT approx_all_masks(MaskService *mask_service, const std::string &method) {
#ifdef TBB
	tbb::parallel_for(1, MASK_COUNT + 1, 1, [&, mask_service](int mask) {
		approx_mask(mask_service, method, mask);
	});
#else
	for (int i = MASK_COUNT; i > 0; i--) { // i > 0
		approx_mask(mask_service, method, i);
	}
#endif
	return S_OK;
}

HRESULT handle_all_segments(const std::string &filename, const std::string &method) {
	std::vector<std::vector<TGlucoseLevel>> segments;
	if (load_segments(filename, segments) == S_FALSE) { return S_FALSE; }
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	
	for (size_t i = 0; i < segments.size(); i++) {
		printf("Getting all masks for segment %zd\n", i);
		MaskService mask_service(&segments[i][0], segments[i].size());
		printf("Got all masks for segment %zd\n", i);
		approx_all_masks(&mask_service, method);
		printf("Counted all masks for segment %zd\n", i);
	}
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Total time = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "ms" << std::endl;
	return S_OK;
}

void print_help(char *name) {
	std::cout << name <<" usage:\n" <<
		"\t-h show help\n" <<
		"\t-f FILE\n\t\tinput file location, by default: ./data/direcnet.sqlite\n" <<
		"\t-m METHOD\n\t\tapproximation method, available values:\n" <<
		"\t\t\t\'akima\' or 'a' for Akima spline (default)\n" <<
		"\t\t\t\'catmull\' or 'cr' for CatmullRom spline\n" <<
		"\t\t\t\'cubic\' or 'c' for Cubic spline.\n" <<
		"Example:\n" <<
		"\t" << name << " -f ../../data/input.sqlite -m a\n";
}

int main(int argc, char *argv[]) {
	//system("pause");
	std::string filename, method;
	ArgParser parser(argc, argv);
	if (parser.check_option("-h")) {
		print_help(argv[0]);
		return 0;
	}
	filename = parser.get_option("-f");
	if (filename.empty()) { filename = DB_FILE; }
	method = parser.get_option("-m");

	
	handle_all_segments(filename, method);
	system("pause");
	return 0;
}
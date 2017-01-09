#include "common/iface/ApproxIface.h"
#include "approx/src/GlucoseLevels.h"
#include <vector>
#include "DBDataService.h"
#include "MaskService.h"
#include "Statistics.h"
#include <tbb/tbb.h>
#include <iostream>
#include <algorithm>
#include <iterator>
#include "ArgParser.h"
#include "defs.h"

HRESULT load_segments(const std::string &filename, std::vector<IGlucoseLevels *> &segments, std::vector<std::string> &segment_ids) {
	DBDataService dbservice(filename);
	DataService *data_service = &dbservice;
	if (data_service->load_segments(segments, segment_ids) != S_OK) {
		std::cerr << "Failed to load segments from database." << std::endl;
		return S_FALSE;
	}
	return S_OK;
}

HRESULT approx_mask(MaskService *mask_service, const std::string &method, CCommonApprox **approx, int mask) {
	IGlucoseLevels *levels;
	mask_service->get_mask(&levels, mask);
	parse_method(method, levels, approx);
	(*approx)->Approximate(nullptr);
	return S_OK;
}

HRESULT get_ref_devs(MaskService *mask_service, CCommonApprox *approx, std::map<floattype, floattype> &ref_devs) {
	IGlucoseLevels *ref_lvls;
	TGlucoseLevel *lvls;
	size_t size, filled, mask_count;
	mask_service->get_mask_count(&mask_count);
	mask_service->get_segment_size(&size);
	mask_service->get_mask(&ref_lvls, static_cast<uint8_t>(mask_count));
	ref_lvls->GetLevels(&lvls);
	floattype dev;
	for (size_t i = 0; i < size; i++) {
		approx->GetLevels(lvls[i].datetime, 0, 1, &dev, &filled, 1);
		ref_devs.insert(std::make_pair(lvls[i].datetime, dev));
	}
	return S_OK;
}

HRESULT approx_all_masks(MaskService *mask_service, const std::string &method) {
	size_t mask_count;
	mask_service->get_mask_count(&mask_count);
	std::vector<CCommonApprox *> approxs(mask_count);
	std::map<floattype, floattype> ref_devs;
#ifdef TBB
	tbb::parallel_for(1, MASK_COUNT + 1, 1, [&, mask_service](int mask) {
		approx_mask(mask_service, method, &approxs[mask - 1], mask);
	});
#else
	for (int i = static_cast<int>(mask_count); i > 0; i--) {
		approx_mask(mask_service, method, &approxs[i - 1], i);
	}
#endif
	get_ref_devs(mask_service, approxs[mask_count - 1], ref_devs);
	std::cout << method << std::endl;
	Statistics stats(mask_service, &ref_devs, false);
	for (int i = static_cast<int>(mask_count); i > 0; i--) {
		stats.get_stats(i, approxs[i - 1]);
		approxs[i - 1]->Release();
	}
	return S_OK;
}

HRESULT get_single_result(IGlucoseLevels *lvls, std::string method, int mask) {
	std::map<floattype, floattype> ref_devs;
	CCommonApprox *approx = NULL;
	MaskService mask_service(lvls);
	approx_mask(&mask_service, method, &approx, mask);
	get_ref_devs(&mask_service, approx, ref_devs);
	Statistics stats(&mask_service, &ref_devs, true);
	std::cout << method << std::endl;
	stats.get_stats(mask, approx);
	approx->Release();
	return S_OK;
}

HRESULT compute(const std::string &filename, std::string method, std::string mask, std::string segment) {
	std::vector<IGlucoseLevels *> segments;
	std::vector<std::string> segment_ids;
	if (load_segments(filename, segments, segment_ids) == S_FALSE) { return S_FALSE; }
	Timer timer("Total time");
	timer.start();
	if (!segment.empty() || !mask.empty()) {
		if (mask.empty()) { mask = "255"; }
		int seg_id;
		auto it = std::find(segment_ids.begin(), segment_ids.end(), segment);
		if (it == segment_ids.end()) { seg_id = 0; }
		else { 
			seg_id = static_cast<int>(std::distance(segment_ids.begin(), it));
		}
		get_single_result(segments[seg_id], method, std::stol(mask, 0, 10));
	} else {
		for (size_t i = 0; i < segments.size(); i++) {
			MaskService mask_service(segments[i]);
			approx_all_masks(&mask_service, method);
		}
	}
	for (size_t i = 0; i < segments.size(); i++) {
		segments[i]->Release();
	}
	timer.stop();
	return S_OK;
}

void print_help(char *name) {
	std::cout << name << " usage:\n" <<
		"\t-h show help\n" <<
		"\t-f FILE\n\t\tinput file location, by default: ./data/direcnet.sqlite\n" <<
		"\t-m METHOD\n\t\tapproximation method, available values:\n" <<
		"\t\t\t\'akima\' or 'a' for Akima spline (default)\n" <<
		"\t\t\t\'catmull\' or 'cr' for CatmullRom spline\n" <<
		"\t\t\t\'cubic\' or 'c' for Cubic spline.\n" <<
		"\t-mask MASK [1-255]\n\t\twhen specified, prints csv graph and stats for the chosen mask and segment (default mask is 255)\n" <<
		"\t-s or -segment SEGMENTID [2-127]\n\t\twhen specified, prints csv graph and stats for the segment and mask (default segmentID is 2)\n" <<
		"Example:\n" <<
		"\t" << name << " -f ../../data/input.sqlite -m akima\n" <<
		"\t" << name << " -f ../../data/input.sqlite -m c -mask 128 -s 7\n" <<
		"\t(Prints out a csv graph for segmentID 7 and mask 128 using Cubic Spline)\n";
}

int main(int argc, char *argv[]) {
	std::string filename, method;
	ArgParser parser(argc, argv);
	if (parser.check_option("-h")) {
		print_help(argv[0]);
		return 0;
	}
	filename = parser.get_option("-f");
	if (filename.empty()) { filename = DB_FILE; }
	method = parser.get_option("-m");
	
	compute(filename, method, parser.get_option("-mask"), parser.get_segment());
	system("pause");
	return 0;
}
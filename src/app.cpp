#include "common/iface/ApproxIface.h"
#include "approx/src/GlucoseLevels.h"
#include <vector>
#include "DBDataService.h"
#include "MaskService.h"
#include "Statistics.h"
#include <tbb/tbb.h>
#include <iostream>
#include "ArgParser.h"
#include "defs.h"

HRESULT load_segments(const std::string &filename, std::vector<IGlucoseLevels *> &segments, std::vector<std::string> &segment_ids) {
	DBDataService dbservice(filename);
	DataService *data_service = &dbservice;

	if (data_service->get_segments(segments, segment_ids) != S_OK) {
		std::cerr << "Failed to load segments from database." << std::endl;
		return S_FALSE;
	}
	return S_OK;
}

HRESULT approx_mask(MaskService *mask_service, const std::string &method, CCommonApprox **approx, int mask) {
	IGlucoseLevels *levels;
	mask_service->get_mask(&levels, mask);
	//if (levels != NULL) levels->AddRef();
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
	mask_service->get_mask(&ref_lvls, mask_count);
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
		approx_mask(mask_service, method, &approxs[i - 1], mask);
	});
#else
	for (int i = mask_count; i > 0; i--) {
		approx_mask(mask_service, method, &approxs[i - 1], i);
	}
#endif
	get_ref_devs(mask_service, approxs[mask_count - 1], ref_devs);
	std::cout << method << std::endl;
	for (int i = mask_count; i > 0; i--) {
		Statistics stats(mask_service, &ref_devs, i, approxs[i - 1], false);
		std::cout << stats.get_output();
		delete approxs[i - 1];
	}
	return S_OK;
}

HRESULT one_result(IGlucoseLevels *lvls, const std::string method, int mask) {
	std::map<floattype, floattype> ref_devs;
	CCommonApprox *approx = NULL;
	MaskService mask_service(lvls);
	approx_mask(&mask_service, method, &approx, mask);
	get_ref_devs(&mask_service, approx, ref_devs);
	Statistics stats(&mask_service, &ref_devs, mask, approx, true);
	std::cout << method << std::endl;
	std::cout << stats.get_output();
	delete approx;
	return S_OK;
}

HRESULT handle_all_segments(const std::string &filename, std::string method, unsigned int m) {
	std::vector<IGlucoseLevels *> segments;
	std::vector<std::string> segment_ids;
	if (load_segments(filename, segments, segment_ids) == S_FALSE) { return S_FALSE; }
	Timer timer;
	timer.start();
	if (m != MASK_COUNT + 1) {
		one_result(segments[0], method, m);
	} else {
		for (size_t i = 0; i < 1; i++) {
			MaskService mask_service(segments[0]);
			approx_all_masks(&mask_service, method);
			std::cout << "Counted all masks for segment " << i << std::endl;
		}
	}
	timer.stop();
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
		"\t-mask [1-255] prints graph csv and stats for chosen mask on first segment" <<
		"Example:\n" <<
		"\t" << name << " -f ../../data/input.sqlite -m a\n";
}

int main(int argc, char *argv[]) {
	system("pause");
	std::string filename, method, mask;
	unsigned int m = MASK_COUNT + 1;
	ArgParser parser(argc, argv);
	if (parser.check_option("-h")) {
		print_help(argv[0]);
		return 0;
	}
	filename = parser.get_option("-f");
	if (filename.empty()) { filename = DB_FILE; }
	mask = parser.get_option("-mask");
	method = parser.get_option("-m");
	if (!mask.empty()) {
		m = std::stol(mask, 0, 10);
		std::cout << "mask " << m << std::endl;
		std::cout << "mask2 " << mask.data() << std::endl;
	}
	
	handle_all_segments(filename, method, m);
	system("pause");
	return 0;
}

#include "common/iface/ApproxIface.h"
#include "approx/src/GlucoseLevels.h"
#include <vector>
#include "DBDataService.h"
#include "approx\src\CubicApprox.h"
#include "approx\src\AkimaApprox.h"
#include "MaskService.h"
#include "Statistics.h"
#include <tbb/tbb.h>

void print_segments(const std::vector<std::shared_ptr<CGlucoseLevels>> &segments) {
	size_t size;
	TGlucoseLevel *levels;
	for (int i = 0; i < 1; i++) {
		segments[i]->GetLevelsCount(&size);
		segments[i]->GetLevels(&levels);
		for (int j = 0; j < size; j++) {
			fprintf(stdout, "%d %zd %f %f\n", i, size, levels[j].datetime, levels[j].level);
		}
	}
}

void load_segments(std::vector<std::vector<TGlucoseLevel>> &segments) {
	DBDataService dbservice("data\\direcnet.sqlite");
	DataService *data_service = &dbservice;

	if (data_service->get_segments(segments) != S_OK) {
		printf("Failed to load segments from database.\n");
		system("pause");
		exit(1);
	}
}

HRESULT approx_all_masks(MaskService *mask_service) {
	CCommonApprox *approx;
	IGlucoseLevels *levels;
	TApproximationParams params;
	/*tbb::parallel_for(0, 255, 1, [&, mask_service](int i) {
		mask_service->get_mask(&levels, i);
		if (levels) printf("null\n");
		CubicApprox cubic(levels);
		approx = &cubic;
		approx->Approximate(&params);
		Statistics stats(mask_service, i, approx);
		printf("Got mask %d\n", i);
	});*/

	for (int i = 1; i > 0; i--) { // 255
		mask_service->get_mask(&levels, i);
		CubicApprox cubic(levels);
		approx = &cubic;
		approx->Approximate(&params);
		Statistics stats(mask_service, i, approx);
	}

	return S_OK;
}

HRESULT handle_all_segments() {
	std::vector<std::vector<TGlucoseLevel>> segments;
	MaskService *mask_service;

	load_segments(segments);
	for (size_t i = 0; i < 1; i++) { //size
		printf("Getting all masks for segment %zd\n", i);
		mask_service = new MaskService(&segments[1][0], segments[1].size());
		printf("Got all masks for segment %zd\n", i);
		approx_all_masks(mask_service);
		printf("Counted all masks for segment %zd\n", i);
		delete mask_service;
	}
	printf("Done all segments.\n");
	return S_OK;
}

int main() {
	system("pause");
	handle_all_segments();
	system("pause");
	return 0;
}
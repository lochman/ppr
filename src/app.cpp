
#include "common/iface/ApproxIface.h"
#include "approx/src/GlucoseLevels.h"
#include <vector>
#include "DBDataService.h"
#include "approx\src\CubicApprox.h"
#include "MaskService.h"

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

void load_segments(std::vector<std::shared_ptr<CGlucoseLevels>> &segments) {
	
	std::shared_ptr<DataService> data_service = make_shared_reference(new DBDataService("data\\direcnet.sqlite"), true);

	if (data_service->get_segments(segments) != S_OK) {
		system("pause");
		exit(1);
	}

	//print_segments(segments);
}

int main() {
	std::vector<std::shared_ptr<CGlucoseLevels>> segments;
	IGlucoseLevels *masked = new CGlucoseLevels;
	MaskService mask_service;
	uint8_t mask = 250;
	load_segments(segments);

	fprintf(stdout, "Creating mask\n");
	mask_service.get_masked_values(masked, *segments[0], mask);
	fprintf(stdout, "Mask ready\n");
	CCommonApprox *ca = new CubicApprox(segments[0].get());
	TApproximationParams params;
	floattype *levels;
	size_t count = 10, filled = 0;
	ca->Approximate(&params);
	levels = (floattype *) malloc(count * sizeof floattype);
	if (levels == NULL) { printf("Alloc error\n"); return 1; }
	//ca->GetLevels(2451545.282639, 0.02, count, levels, &filled, filled);
	free(levels);
	/*
	mask_service.get_inverse_masked_values(masked, *segments[0], mask);
	CCommonApprox *ca = new CubicApprox(segments[0].get());
	TApproximationParams params;
	floattype *levels;
	size_t count = 10, filled = 0;
	ca->Approximate(&params);
	levels = (floattype *) malloc(count * sizeof floattype);
	if (levels == NULL) { printf("Alloc error\n"); return 1; }
	
	free(levels);
	*/
	system("pause");
	return 0;
}
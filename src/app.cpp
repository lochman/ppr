
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

//#define TBB
#define MASKS 255 // 255
#define SEGMENTS segments.size() // segments.size()
#define SEGMENT i // i

HRESULT approx_all_masks(MaskService *mask_service) {
#ifdef TBB
	tbb::parallel_for(1, MASKS + 1, 1, [&, mask_service](int mask) { // 256
		IGlucoseLevels *levels;
		//printf("getting mask %d\n", mask);
		mask_service->get_mask(&levels, mask);
		AkimaApprox approx(levels);
		approx.Approximate(nullptr);
		Statistics stats(mask_service, mask, &approx);
		//printf("got mask %d\n", mask);
	});
#else
	for (int i = MASKS; i > 0; i--) {
		CCommonApprox *approx;
		IGlucoseLevels *levels;
		mask_service->get_mask(&levels, i);
		CubicApprox cubic(levels);
		approx = &cubic;
		approx->Approximate(nullptr);
		Statistics stats(mask_service, i, approx);
	}
#endif
	return S_OK;
}

HRESULT handle_all_segments() {
	std::vector<std::vector<TGlucoseLevel>> segments;

	load_segments(segments);

	for (size_t i = 0; i < SEGMENTS; i++) {
		printf("Getting all masks for segment %zd\n", i);
		MaskService mask_service(&segments[SEGMENT][0], segments[SEGMENT].size());
		printf("Got all masks for segment %zd\n", i);
		approx_all_masks(&mask_service);
		printf("Counted all masks for segment %zd\n", i);
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
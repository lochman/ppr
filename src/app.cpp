
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
#define MASKS 1		  // 1 | 255
#define SEGMENTS 1		  // 1 | segments.size()
#define SEGMENT 1		  // 1 | i

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
	for (int i = MASKS; i == MASKS; i--) { // i > 0
		CCommonApprox *approx;
		IGlucoseLevels *levels;
		mask_service->get_mask(&levels, i);
		//CatmullRomApprox cubic(levels);
		AkimaApprox cubic(levels);
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
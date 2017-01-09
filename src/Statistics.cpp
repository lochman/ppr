#include "Statistics.h"
#include <algorithm>
#include <numeric>
#include <iostream>

#undef min
#undef max	

floattype min_error(std::vector<floattype> &errors) {
	return *std::min_element(errors.begin(), errors.end());
}

floattype max_error(std::vector<floattype> &errors) {
	return *std::max_element(errors.begin(), errors.end());
}

double quantil(std::vector<floattype> &errors, const floattype quantil) {
	size_t size = errors.size(), q_ceil, q_floor;
	floattype result;
	if (size == 1) { return errors[0]; }
	q_ceil = static_cast<size_t>(std::ceil(size * quantil)); // quantil may be counted from 2 values, it depends on the vector size
										// so it is better to get two nearby elements and count the average
	q_ceil = std::min(q_ceil, errors.size() - 2); // prevent overflow
	q_floor = std::max(q_ceil - 1, size_t(0));	// prevent underflow
	std::nth_element(errors.begin(), errors.begin() + q_ceil, errors.end());
	result = errors[q_ceil];

	std::nth_element(errors.begin(), errors.begin() + q_floor, errors.begin() + q_ceil);
	result = (errors[q_floor] + result) * 0.5;
	return result;
}

double mean(std::vector<floattype> &errors) {
	return std::accumulate(errors.begin(), errors.end(), 0.0) / (float) errors.size();
}

double std_deviation(std::vector<floattype> &errors, const floattype mean) {
	double square_sum, dif;
	if (errors.size() == 1) { return 0; }
	std::vector<double> diff(errors.size());
	std::transform(errors.begin(), errors.end(), diff.begin(), [mean](double x) { return x - mean; });
	square_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	dif = square_sum / (float) (errors.size() - 1);
	return dif > 0 ? std::sqrt(dif) : 0;
}

void print_graph(TGlucoseLevel *levels, const int mask, std::vector<floattype> &errors) {
	printf("time,reference,\"mask %d\"\n", mask);
	for (size_t i = 0; i < errors.size(); i++) {
		printf("%f,%f,%f\n", levels[i].datetime, levels[i].level, errors[i]);
	}
}

void print_graph_new(TGlucoseLevel *levels, size_t &size, size_t &steps, const int mask, CCommonApprox *approx) {
	std::vector<floattype> approx_lvls(size);
	floattype from = levels[0].datetime,
			  to = levels[size - 1].datetime,
			  stepsize = (to - from) / (floattype)steps;
	size_t filled, j = 0;
	approx->GetLevels(from, stepsize, steps, &approx_lvls[0], &filled, 0);
	
	printf("time,reference\n");
	for (size_t i = 0; i < size; i++) {
		printf("%f,%f\n", levels[i].datetime, levels[i].level);
	}
	printf("time,\"mask %d\"\n", mask);
	for (size_t i = 0; i < approx_lvls.size(); i++) {
		printf("%f,%f\n", from + i * stepsize, approx_lvls[i]);
	}
}

HRESULT Statistics::get_errors(TGlucoseLevel *levels, size_t size, const int mask, CCommonApprox *approx, boolean graph) {
	std::vector<floattype> approx_lvls(size), abs_errors(size), rel_errors(size),
			approx_lvls_dev(size), abs_errors_dev(size);
	size_t filled;

	for (size_t i = 0; i < size; i++) {
		approx->GetLevels(levels[i].datetime, 0, 1, &approx_lvls[i], &filled, 0);
		abs_errors[i] = (std::abs(ref_lvls[i].level - approx_lvls[i]));
		rel_errors[i] = (abs_errors[i] / ref_lvls[i].level);
		filled = 0;
		approx->GetLevels(levels[i].datetime, 0, 1, &approx_lvls_dev[i], &filled, 1);
		abs_errors_dev[i] = (std::abs((*ref_devs)[levels[i].datetime] - approx_lvls_dev[i]));
	}
	if (graph) { print_graph(levels, mask, approx_lvls); }
	//print_graph_new(levels, size, steps, mask, approx);
	output << "\tabs: ";
	print_stats(abs_errors);
	output << "\trel: ";
	print_stats(rel_errors);
	output << "\t1.der: ";
	print_stats(abs_errors_dev);
	return S_OK;
}

Statistics::Statistics(MaskService *mask_service, std::map<floattype, floattype> *ref_devs, const int mask, CCommonApprox *approx, boolean graph) : ref_devs(ref_devs) {
	IGlucoseLevels *glevels;
	TGlucoseLevel *lvls;
	size_t mask_count, lvl_count;
	mask_service->get_mask_count(&mask_count);
	//printf("Mask count is %zd", mask_count);
	mask_service->get_mask(&glevels, static_cast<unsigned int>(mask_count));
	glevels->GetLevels(&ref_lvls);
	glevels->GetLevelsCount(&size);
	//printf("Mask size is %zd", size);
	//mask_service->get_segment_size(&size);
	output << "  mask 0x" << std::hex << mask << ":" << std::endl;
	output << "    all:" << std::endl;
	get_errors(ref_lvls, size, mask, approx, graph);

	output << "    bit 0:" << std::endl;
	mask_service->get_inverse_mask(&glevels, mask);
	glevels->GetLevels(&lvls);
	glevels->GetLevelsCount(&lvl_count);
	get_errors(lvls, lvl_count, mask, approx, false);

	output << "    bit 1:" << std::endl;
	mask_service->get_mask(&glevels, mask);
	glevels->GetLevels(&lvls);
	glevels->GetLevelsCount(&lvl_count);
	get_errors(lvls, lvl_count, mask, approx, false);
}

HRESULT Statistics::print_stats(std::vector<floattype> &errors) {
	if (errors.size() == 0) { return S_FALSE; }
	floattype m = mean(errors);
	//printf("mean,min,Q1,median,Q3,max,std_dev\n");
	//printf("%f,%f,%f,%f,%f,%f,%f\n", m, min_error(errors), quantil(errors, 0.25), quantil(errors, 0.5),quantil(errors, 0.75), max_error(errors), std_deviation(errors, m));
	output << m << "," << min_error(errors) << "," << quantil(errors, 0.25) << "," << quantil(errors, 0.5)
		<< "," << quantil(errors, 0.75) << "," << max_error(errors) << "," << std_deviation(errors, m) << std::endl;
	return S_OK;
}

std::string Statistics::get_output() {
	return output.str();
}

HRESULT Timer::start() {
	begin = std::chrono::steady_clock::now();
	return S_OK;
}

HRESULT Timer::stop() {
	end = std::chrono::steady_clock::now();
	std::cout << message << " = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds" << std::endl;
	return S_OK;
}
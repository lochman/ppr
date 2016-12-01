#include "Statistics.h"
#include <algorithm>
#include <numeric>

#undef min
#undef max	

floattype min_error(std::vector<floattype> &errors) {
	return *std::min_element(errors.begin(), errors.end());
}

floattype max_error(std::vector<floattype> &errors) {
	return *std::max_element(errors.begin(), errors.end());
}

floattype quantil(std::vector<floattype> &errors, floattype const &quantil) {
	size_t size = errors.size(), q_ceil, q_floor;
	floattype result;
	if (size == 1) { return errors[0]; } //TODO fix ceil
	q_ceil = std::ceil(size * quantil); // quantil may be counted from 2 values, it depends on the vector size
										// so it is better to get two nearby elements and count the average
	q_ceil = std::min(q_ceil, errors.size() - 2); // prevent overflow
	q_floor = std::max(q_ceil - 1, size_t(0)); // prevent underflow
	std::nth_element(errors.begin(), errors.begin() + q_ceil, errors.end());
	//printf("q = %f, size = %zd, floor = %zd, ceil = %zd\n", quantil, errors.size(), q_floor, q_ceil);
	result = errors[q_ceil];

	std::nth_element(errors.begin(), errors.begin() + q_floor, errors.begin() + q_ceil);
	result = (errors[q_floor] + result) * 0.5;
	return result;
}

floattype mean(std::vector<floattype> &errors) {
	return std::accumulate(errors.begin(), errors.end(), 0.0) / errors.size();
}

floattype std_deviation(std::vector<floattype> &errors, floattype &mean) {
	double square_sum;
	std::vector<double> diff(errors.size());
	std::transform(errors.begin(), errors.end(), diff.begin(), [mean](double x) { return x - mean; });
	square_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	return std::sqrt(square_sum / errors.size() - 1);
}

void print_reference(TGlucoseLevel *ref_values, size_t &size) {
	for (size_t i = 0; i < size; i++) {
		printf("%f,%f\n", ref_values[i].datetime, ref_values[i].level);
	}
}

void print_graph(TGlucoseLevel *levels, int &mask, std::vector<floattype> &errors) {
	printf("time,reference,\"mask %d\"\n", mask);
	for (size_t i = 0; i < errors.size(); i++) {
		printf("%f,%f,%f\n", levels[i].datetime, levels[i].level, errors[i]);
	}
}

void Statistics::get_errors(TGlucoseLevel *levels, int &mask, CCommonApprox *approx) {
	std::vector<floattype> approx_lvls(size), abs_errors, rel_errors;
	size_t filled;

	for (size_t i = 0; i < size; i++) {
		approx->GetLevels(levels[i].datetime, 0, 1, &approx_lvls[i], &filled, 0);
		abs_errors.push_back(std::abs(levels[i].level - approx_lvls[i]));
		rel_errors.push_back(abs_errors[i] / levels[i].level);
	}
	//printf("errors:\n");
	//print_graph(levels, mask, approx_lvls);
	/*
	print_stats(abs_errors);
	print_stats(rel_errors);
	*/
}

Statistics::Statistics(const MaskService *mask_service, int &mask, CCommonApprox *approx) {
	IGlucoseLevels *glevels;
	TGlucoseLevel *levels;
	mask_service->get_levels(&ref_values);
	mask_service->get_levels_size(&size);
	//print_reference(ref_values, size);
	get_errors(ref_values, mask, approx);
	/*
	mask_service->get_mask(&glevels, mask);
	glevels->GetLevels(&levels);
	get_errors(levels, approx);
	mask_service->get_inverse_mask(&glevels, mask);
	glevels->GetLevels(&levels);
	get_errors(levels, approx);
	*/
}

void Statistics::print_stats(std::vector<floattype> &errors) {
	if (errors.size() == 0) { return; }
	floattype m = mean(errors);
	//printf("mean,min,Q1,median,Q3,max,std_dev\n");
	printf("%f,%f,%f,%f,%f,%f,%f\n", m, min_error(errors), quantil(errors, 0.25), quantil(errors, 0.5),
		quantil(errors, 0.75), max_error(errors), std_deviation(errors, m));
}
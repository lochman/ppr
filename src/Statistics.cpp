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

HRESULT Statistics::get_errors(TGlucoseLevel *levels, size_t size, const int mask, CCommonApprox *approx, boolean graph, std::stringstream &output) {
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
	output << "\t\t\tabs: ";
	output << print_stats(abs_errors);
	output << "\t\t\trel: ";
	output << print_stats(rel_errors);
	output << "\t\t\t1.der: ";
	output << print_stats(abs_errors_dev);
	return S_OK;
}

Statistics::Statistics(MaskService *mask_service, std::map<floattype, floattype> *ref_devs, boolean graph) :
			ref_devs(ref_devs), mask_service(mask_service), graph(graph) {
	IGlucoseLevels *glevels;
	size_t mask_count;
	mask_service->get_mask_count(&mask_count);
	mask_service->get_mask(&glevels, static_cast<unsigned int>(mask_count));
	glevels->GetLevels(&ref_lvls);
	glevels->GetLevelsCount(&size);
}

std::string Statistics::get_stats(const int mask, CCommonApprox *approx) {
	IGlucoseLevels *glevels;
	TGlucoseLevel *lvls;
	size_t lvl_count;
	std::stringstream output;
	//output.str(std::string());
	//output.clear();
	output << "\tmask 0x" << std::hex << mask << ":\n";
	output << "\t\tall:\n";
	get_errors(ref_lvls, size, mask, approx, graph, output);

	output << "\t\tbit 0:\n";
	mask_service->get_inverse_mask(&glevels, mask);
	glevels->GetLevels(&lvls);
	glevels->GetLevelsCount(&lvl_count);
	get_errors(lvls, lvl_count, mask, approx, false, output);

	output << "\t\tbit 1:\n";
	mask_service->get_mask(&glevels, mask);
	glevels->GetLevels(&lvls);
	glevels->GetLevelsCount(&lvl_count);
	get_errors(lvls, lvl_count, mask, approx, false, output);
	return output.str();
}

std::string Statistics::print_stats(std::vector<floattype> &errors) {
	std::stringstream output;
	if (errors.size() == 0) { return ""; }
	floattype m = mean(errors);
	output << m << "," << min_error(errors) << "," << quantil(errors, 0.25) << "," << quantil(errors, 0.5)
		<< "," << quantil(errors, 0.75) << "," << max_error(errors) << "," << std_deviation(errors, m) << "\n";
	return output.str();
}

HRESULT Timer::start() {
	begin = std::chrono::steady_clock::now();
	return S_OK;
}

long long Timer::stop() {
	end = std::chrono::steady_clock::now();
	return std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
}
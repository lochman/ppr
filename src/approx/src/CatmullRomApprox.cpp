#include "CatmullRomApprox.h"
#include <amp.h>
#include <amp_math.h>

const floattype alpha = 0.5f; // centripetal parametrization

floattype get_distance(const TGlucoseLevel &p0, const TGlucoseLevel &p1) {
	const floattype time_diff = std::pow(p1.datetime - p0.datetime, 2),
					lvl_diff = std::pow(p1.level - p0.level, 2);
	return std::pow(std::sqrt(time_diff + lvl_diff), alpha);
}

floattype get_distance(const floattype x0, const floattype y0, const floattype x1, const floattype y1) restrict(amp) {
	const floattype time_diff = concurrency::precise_math::pow(x1 - x0, 2),
					lvl_diff = concurrency::precise_math::pow(y1 - y0, 2);
	return concurrency::precise_math::pow(concurrency::precise_math::sqrt(time_diff + lvl_diff), 0.5); // aplha
}

HRESULT CatmullRomApprox::get_coefficients(const floattype p0, const floattype p1, const floattype t0, const floattype t1, const size_t i) {
	b[i] = t0;
	c[i] = -3 * p0 + 3 * p1 - 2 * t0 - t1;
	d[i] = 2 * p0 - 2 * p1 + t0 + t1;
	return S_OK;
}

HRESULT extrapolated_point(const TGlucoseLevel &p0, const TGlucoseLevel &p1, TGlucoseLevel *ext) {
	ext->datetime = p1.datetime + 0.5 * (p1.datetime - p0.datetime);
	ext->level = p1.level + 0.5 * (p1.level - p0.level);
	return S_OK;
}

HRESULT CatmullRomApprox::extrapolation() {
	TGlucoseLevel ext1, ext2, ext3;
	extrapolated_point(levels[size - 2], levels[size - 1], &ext1);
	extrapolated_point(levels[size - 1], ext1, &ext2);
	extrapolated_point(ext1, ext2, &ext3);

	iterate(levels[size - 3], levels[size - 2], levels[size - 1], ext1, size - 3);
	iterate(levels[size - 2], levels[size - 1], ext1, ext2, size - 2);
	iterate(levels[size - 1], ext1, ext2, ext3, size - 1);
	return S_OK;
}

floattype get_tangent(const floattype p0, const floattype p1, const floattype p2,
		const floattype t0, const floattype t1, const floattype t2) restrict(cpu, amp) {
	return ((p1 - p0) / t0 - (p2 - p0) / (t0 + t1) + (p2 - p1) / t1) *t2;
}

HRESULT CatmullRomApprox::iterate(const TGlucoseLevel &p0, const TGlucoseLevel &p1, const TGlucoseLevel &p2, const TGlucoseLevel &p3, const size_t i) {
	floattype t0, t1, t2, tan1, tan2;
	t0 = get_distance(p0, p1);
	t1 = get_distance(p1, p2);
	t2 = get_distance(p2, p3);
	tan1 = get_tangent(p0.level, p1.level, p2.level, t0, t1, t1);
	tan2 = get_tangent(p1.level, p2.level, p3.level, t1, t2, t1);

	get_coefficients(p0.level, p1.level, tan1, tan2, i);
	return S_OK;
}

HRESULT CatmullRomApprox::approximate_gpu() {
	std::vector<floattype> vec_times(size), vec_lvls(size);
	for (size_t i = 0; i < size; i++) {
		vec_times[i] = levels[i].datetime;
		vec_lvls[i] = levels[i].level;
	}

	concurrency::extent<1> ext(size - 3);
	const concurrency::array_view<const floattype, 1> times(size, vec_times), lvls(size, vec_lvls);
	concurrency::array_view<floattype, 1> b_view(size, b), c_view(size, c), d_view(size, d);
	//a_view.discard_data();
	b_view.discard_data();
	c_view.discard_data();
	d_view.discard_data();
	concurrency::parallel_for_each(ext, [=](concurrency::index<1> idx) restrict(amp) {
		const int i = idx[0];
		const floattype t0 = get_distance(lvls[i], lvls[i], times[i + 1], times[i + 1]),
						t1 = get_distance(lvls[i + 1], lvls[i + 1], times[i + 2], times[i + 2]),
						t2 = get_distance(lvls[i + 2], lvls[i + 2], times[i + 3], times[i + 3]),
						tan1 = get_tangent(lvls[i], lvls[i + 1], lvls[i + 2], t0, t1, t1),
						tan2 = get_tangent(lvls[i + 1], lvls[i + 2], lvls[i + 3], t1, t2, t1);
		//a_view[idx] = lvls[idx];
		b_view[i] = tan1;
		c_view[i] = -3 * lvls[i] + 3 * lvls[i + 1] - 2 * tan1 - tan2;
		d_view[i] = 2 * lvls[i] - 2 * lvls[i + 1] + tan1 + tan2;
	});
	//a_view.synchronize();
	b_view.synchronize();
	c_view.synchronize();
	d_view.synchronize();
	return S_OK;
}

#define GPU

HRESULT CatmullRomApprox::Approximate(TApproximationParams * params) {
	//a.reserve(size);
	b.resize(size);
	c.resize(size);
	d.resize(size);
#ifdef GPU	// GPU
	approximate_gpu();
#else
	for (size_t i = 0; i < size - 3; i++) {
		iterate(levels[i], levels[i + 1], levels[i + 2], levels[i + 3], i);
	}
#endif
	extrapolation();
	return S_OK;
}

HRESULT CatmullRomApprox::GetLevels(floattype desiredtime, floattype stepping, size_t count,
	floattype *levels, size_t *filled, size_t derivationorder) {
	TGlucoseLevel *gl = this->levels;
	floattype time = desiredtime, x;
	int i;
	for (size_t j = 0; j < count; j++) {
		if (get_time_interval(gl, size, time, &i) == S_FALSE) { levels[j] = 0.0; continue; }
		x = (time - gl[i].datetime) / (gl[i + 1].datetime - gl[i].datetime);
		switch (derivationorder) {
		case 1:
			levels[j] = b[i] + 2 * c[i] * x + 3 * d[i] * pow(x, 2);
			break;
		case 2:
			levels[j] = 2 * c[i] + 6 * d[i] * x;
			break;
		default:
			levels[j] = gl[i].level + b[i] * x + c[i] * pow(x, 2) + d[i] * pow(x, 3);
		}
		time += stepping;
		(*filled)++;
	}
	return S_OK;
}
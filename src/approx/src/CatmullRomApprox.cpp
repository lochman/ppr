#include "CatmullRomApprox.h"
#include <amp.h>
#include <amp_math.h>

const floattype alpha = 0.5f; // centripetal parametrization

floattype get_distance(const TGlucoseLevel &p0, const TGlucoseLevel &p1) {
	const floattype time_diff = std::pow(p1.datetime - p0.datetime, 2),
					lvl_diff = std::pow(p1.level - p0.level, 2);
	return std::pow(std::sqrt(time_diff + lvl_diff), alpha);
}

floattype get_distance(const floattype &x0, const floattype &x1, const floattype &y0, const floattype &y1) restrict(amp) {
	const floattype time_diff = concurrency::precise_math::pow(x1 - x0, 2),
					lvl_diff = concurrency::precise_math::pow(y1 - y0, 2);
	return concurrency::precise_math::pow(concurrency::precise_math::sqrt(time_diff + lvl_diff), 0.5); // aplha
}

void CatmullRomApprox::get_coefficients(const floattype &p0, const floattype &p1, const floattype t0, const floattype t1) {
	a.push_back(p0);
	b.push_back(t0);
	c.push_back(-3 * p0 + 3 * p1 - 2 * t0 - t1);
	d.push_back(2 * p0 - 2 * p1 + t0 + t1);
}

floattype extrapolated_point(const floattype &p0, const floattype &p1) {
	return p1 + 0.5 * (p1 - p0);
}

void CatmullRomApprox::extrapolation() {
	TGlucoseLevel ext1, ext2, ext3;
	ext1.datetime = extrapolated_point(levels[size - 2].datetime, levels[size - 1].datetime);
	ext1.level = extrapolated_point(levels[size - 2].level, levels[size - 1].level);
	
	ext2.datetime = extrapolated_point(levels[size - 1].datetime, ext1.datetime);
	ext2.level = extrapolated_point(levels[size - 1].level, ext1.level);

	ext3.datetime = extrapolated_point(ext1.datetime, ext2.datetime);
	ext3.level = extrapolated_point(ext1.level, ext2.level);

	iterate(levels[size - 3], levels[size - 2], levels[size - 1], ext1);
	iterate(levels[size - 2], levels[size - 1], ext1, ext2);
	iterate(levels[size - 1], ext1, ext2, ext3);
}

void CatmullRomApprox::get_tangent(const floattype &p0, const floattype &p1, const floattype &p2,
		const floattype &p3, const floattype &t0, const floattype &t1, const floattype &t2) {
	floattype tan1 = (p1 - p0) / t0 - (p2 - p0) / (t0 + t1) + (p2 - p1) / t1,
			  tan2 = (p2 - p1) / t1 - (p3 - p1) / (t1 + t2) + (p3 - p2) / t2;
	tan1 *= t1;
	tan2 *= t1;
	get_coefficients(p0, p1, tan1, tan2);
}

void CatmullRomApprox::iterate(const TGlucoseLevel &p0, const TGlucoseLevel &p1, const TGlucoseLevel &p2, const TGlucoseLevel &p3) {
	floattype t0, t1, t2;
	t0 = get_distance(p0, p1);
	t1 = get_distance(p1, p2);
	t2 = get_distance(p2, p3);
	get_tangent(p0.level, p1.level, p2.level, p3.level, t0, t1, t2);
}

void CatmullRomApprox::approximate_gpu() {
	std::vector<floattype> vec_times(size), vec_lvls(size);
	floattype m_next;
	for (size_t i = 0; i < size; i++) {
		vec_times[i] = levels[i].datetime;
		vec_lvls[i] = levels[i].level;
	}

	concurrency::extent<1> ext(size - 3);
	const concurrency::array_view<const floattype, 1> times(size, vec_times), lvls(size, vec_lvls);
	concurrency::array_view<floattype, 1> a_view(size, a), b_view(size, b), c_view(size, c), d_view(size, d);
	a_view.discard_data();
	b_view.discard_data();
	c_view.discard_data();
	d_view.discard_data();
	concurrency::parallel_for_each(ext, [=](concurrency::index<1> idx) restrict(amp) {
		floattype t0, t1, t2;
		t0 = get_distance(lvls[idx], lvls[idx], times[idx + 1], times[idx + 1]);
		t1 = get_distance(lvls[idx + 1], lvls[idx + 1], times[idx + 2], times[idx + 2]);
		t2 = get_distance(lvls[idx + 2], lvls[idx + 2], times[idx + 3], times[idx + 3]);

		a_view[idx];
		b_view[idx];
		c_view[idx];
		d_view[idx];
	});
	a_view.synchronize();
	b_view.synchronize();
	c_view.synchronize();
	d_view.synchronize();
}

#define GPU
HRESULT CatmullRomApprox::Approximate(TApproximationParams * params) {
	a.reserve(size);
	b.reserve(size);
	c.reserve(size);
	d.reserve(size);
#ifdef GPU	// GPU
	approximate_gpu();
#else
	for (size_t i = 0; i < size - 3; i++) {
		iterate(levels[i], levels[i + 1], levels[i + 2], levels[i + 3]);
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
		default:
			levels[j] = a[i] + b[i] * x + c[i] * pow(x, 2) + d[i] * pow(x, 3);
		}
		time += stepping;
		(*filled)++;
	}
	return S_OK;
}
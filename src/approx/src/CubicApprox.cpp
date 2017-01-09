#include "CubicApprox.h"
#include "../../defs.h"
#include <math.h>
#include <iterator>
#include <algorithm>
#include <amp.h>
#include <deque>

HRESULT CubicApprox::Approximate(TApproximationParams *params) {
	//printf("Starting Cubic with size %zd\n", size); n = size-1
	
	std::vector<floattype> h(size), a(size), l(size), u(size), z(size);
	
	b.resize(size);
	c.resize(size);
	d.resize(size);
	
	l[0] = 1.0;
	u[0] = 0.0;
	z[0] = 0.0;
	l[size - 1] = 1.0;
	z[size - 1] = 0.0;
	c[size - 1] = 0.0;
	h[0] = levels[1].datetime - levels[0].datetime;

#ifdef GPU
	h[1] = levels[2].datetime - levels[1].datetime;
	a[0] = (levels[1].level - levels[0].level) / h[0];
	a[1] = (levels[2].level - levels[1].level) / h[1];
	u[1] = 2.0 * (h[0] + h[1]);
	l[1] = 6.0 * (a[1] - a[0]);
	for (int i = 2; i < size - 1; i++) {
		h[i] = levels[i + 1].datetime - levels[i].datetime;
		a[i] = (levels[i + 1].level - levels[i].level) / h[i];
		u[i] = 2.0*(h[i] + h[i - 1]) - h[i - 1] * h[i - 1] / u[i - 1];
		l[i] = 6.0*(a[i] - a[i - 1]) - h[i - 1] * l[i - 1] / u[i - 1];
	}
	z[size - 1] = 0;
	for (int i = static_cast<int>(size - 2); i > 0; i--) {
		z[i] = (l[i] - h[i] * z[i + 1]) / u[i];
	}
	std::vector<floattype> vec_times(size), vec_lvls(size);

	for (size_t i = 0; i < size; i++) {
		vec_times[i] = levels[i].datetime;
		vec_lvls[i] = levels[i].level;
	}
	int s = static_cast<int>(size);
	concurrency::extent<1> ext(s - 1);
	concurrency::array_view<const floattype> times(s, vec_times), lvls(s, vec_lvls), zs(s, z);
	concurrency::array_view<floattype> bs(ext, b), cs(ext, c), ds(ext, d);
	bs.discard_data();
	cs.discard_data();
	ds.discard_data();
	concurrency::parallel_for_each(ext, [=](concurrency::index<1> i) restrict(amp) {
		const floattype h = times[i + 1] - times[i];
		bs[i] = (-h / 6) * zs[i + 1] - (h / 3) * zs[i] + (1 / h) * (lvls[i + 1] - lvls[i]);
		cs[i] = zs[i] / 2;
		ds[i] = (1 / (6 * h)) * (zs[i + 1] - zs[i]);
	});
	bs.synchronize();
	cs.synchronize();
	ds.synchronize();
#else // end GPU

	for (size_t i = 1; i < size - 1; i++) {
		h[i] = levels[i + 1].datetime - levels[i].datetime;
		a[i] = (3 / h[i]) * (levels[i + 1].level - levels[i].level) -
			   (3 / h[i - 1]) * (levels[i].level - levels[i - 1].level);
		l[i] = 2 * (levels[i + 1].datetime - levels[i - 1].datetime) - h[i - 1] * u[i - 1];
		u[i] = h[i] / l[i];
		z[i] = (a[i] - h[i - 1] * z[i - 1]) / l[i];
	}

	for (int i = static_cast<int>(size - 2); i > -1; i--) {
		c[i] = z[i] - u[i] * c[i + 1];
		b[i] = (levels[i + 1].level - levels[i].level) / h[i] - h[i] * ((c[i + 1] + 2 * c[i]) / 3);
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
	}

#endif
	return S_OK;
}

HRESULT CubicApprox::GetLevels(floattype desiredtime, floattype stepping, size_t count,
	floattype *levels, size_t *filled, size_t derivationorder) {
	TGlucoseLevel *gl = this->levels;
	floattype time = desiredtime, x;
	int i;
	for (size_t j = 0; j < count; j++) {
		if (get_time_interval(gl, size, time, &i) == S_FALSE) { levels[j] = 0.0; continue; }
		x = time - gl[i].datetime;
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

std::string CubicApprox::get_name() { return std::string("Cubic spline"); }
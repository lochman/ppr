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

	std::vector<floattype> vec_times(size), vec_lvls(size);
	std::deque<int> indices(size);

	for (size_t i = 0; i < size; i++) {
		vec_times[i] = levels[i].datetime;
		vec_lvls[i] = levels[i].level;
		indices[i] = i;
	}
	indices.pop_back(); // make sure to iterate on interval <1, size - 1)
	indices.pop_front();

	concurrency::extent<1> ext(&indices[0]);
	concurrency::array_view<const floattype> times(ext, vec_times), lvls(ext, vec_lvls);
	concurrency::array_view<floattype> hs(ext, h), as(ext, a), ls(ext, l), us(ext, u), zs(ext, z);
	as.discard_data();
	concurrency::parallel_for_each(ext, [=](concurrency::index<1> i) restrict(amp) {
		hs[i] = times[i + 1] - times[i];
		as[i] = (3 / hs[i]) * (lvls[i + 1] - lvls[i]) -
				(3 / hs[i - 1]) * (lvls[i] - lvls[i - 1]);
		ls[i] = 2 * (times[i + 1] - times[i - 1]) - hs[i - 1] * us[i - 1];
		us[i] = hs[i] / ls[i];
		zs[i] = (as[i] - hs[i - 1] * zs[i - 1]) / ls[i];
	});
	hs.synchronize();
	as.synchronize();
	ls.synchronize();
	us.synchronize();
	zs.synchronize();

	indices.push_front(0); // add index zero
	printf("%d, %d, %d\n", indices[0], indices[1], indices[2]);
	std::reverse(indices.begin(), indices.end()); // invert indices order

	concurrency::extent<1> ext1(&indices[0]);
	concurrency::array_view<floattype> bs(ext1, b), cs(ext1, c), ds(ext1, d);
	bs.discard_data();
	//cs.discard_data();
	ds.discard_data();
	concurrency::parallel_for_each(ext1, [=](concurrency::index<1> i) restrict(amp) {
		cs[i] = zs[i] - us[i] * cs[i + 1];
		bs[i] = (lvls[i + 1] - lvls[i]) / hs[i] - hs[i] * ((cs[i + 1] + 2 * cs[i]) / 3);
		ds[i] = (cs[i + 1] - cs[i]) / (3 * hs[i]);
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

	for (int i = size - 2; i > -1; i--) {
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
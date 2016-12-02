#include "CatmullRomApprox.h"
//#include <amp.h>
//#include <amp_math.h>

const floattype alpha = 0.5f; // centripetal parametrization
const floattype point_count = 50.0f;

floattype get_next_t(TGlucoseLevel &xi, TGlucoseLevel &xi1, floattype &t) {
	floattype time_diff = std::pow(xi1.datetime - xi.datetime, 2),
			  lvl_diff = std::pow(xi1.level - xi.level, 2);
	return std::pow(std::sqrt(time_diff + lvl_diff), alpha) + t;
}

void get_a(TGlucoseLevel *lvl, TGlucoseLevel &pi, TGlucoseLevel &pi1, floattype &ti, floattype &ti1, floattype &t) {
	lvl->datetime = (ti1 - t) / (ti1 - ti) * pi.datetime + (t - ti) / (ti1 - ti) * pi1.datetime;
	lvl->level = (ti1 - t) / (ti1 - ti) * pi.level + (t - ti) / (ti1 - ti) * pi1.level;
}

void get_a_derivative(TGlucoseLevel *lvl, TGlucoseLevel &pi, TGlucoseLevel &pi1, floattype &ti, floattype &ti1) {
	lvl->datetime = (pi1.datetime - pi.datetime) / (ti1 - ti);
	lvl->level = (pi1.level - pi.level) / (ti1 - ti);
}



floattype CatmullRomApprox::iterate(TGlucoseLevel &p0, TGlucoseLevel &p1, TGlucoseLevel &p2, TGlucoseLevel &p3, floattype n, size_t *i) {
	floattype t, t0, t1, t2, t3, step;
	TGlucoseLevel a1, a2, a3, a1d, a2d, a3d, b1, b2, b1d, b2d, cp, cd, tmp;

	t0 = 0.0f;
	t1 = get_next_t(p0, p1, t0);
	t2 = get_next_t(p1, p2, t1);
	t3 = get_next_t(p2, p3, t2);
	step = (t2 - t1) / n;
	for (floattype t = t1; t < t2; t += step) {
		get_a(&a1, p0, p1, t0, t1, t);
		get_a_derivative(&a1d, p0, p1, t0, t1);

		get_a(&a2, p1, p2, t1, t2, t);
		get_a_derivative(&a2d, p1, p2, t1, t2);

		get_a(&a3, p2, p3, t2, t3, t);
		get_a_derivative(&a3d, p2, p3, t2, t3);

		get_a(&b1, a1, a2, t0, t2, t);
		get_a_derivative(&tmp, a1, a2, t0, t2);
		get_a(&b1d, a1d, a2d, t0, t2, t);
		b1d.datetime += tmp.datetime;
		b1d.level += tmp.level;

		get_a(&b2, a2, a3, t1, t3, t);
		get_a_derivative(&tmp, a2, a3, t1, t3);
		get_a(&b2d, a2d, a3d, t1, t3, t);
		b2d.datetime += tmp.datetime;
		b2d.level += tmp.level;

		get_a(&cp, b1, b2, t1, t2, t);
		get_a_derivative(&tmp, b1, b2, t1, t2); // *t
		get_a(&cd, b1d, b2d, t1, t2, t);
		b1d.datetime += tmp.datetime;
		b1d.level += tmp.level;

		c.push_back(cp);
		d.push_back(cd);
		//printf("%f, %f\n", c[*i].datetime, c[*i].level);
		(*i)++;
	}
	return 0;
}

HRESULT CatmullRomApprox::Approximate(TApproximationParams * params) {
	c.resize(size * point_count * 2);
	d.resize(size * point_count * 2);
	size_t ind = 0;
	for (size_t i = 0; i < size - 3; i++) {
		iterate(levels[i], levels[i + 1], levels[i + 2], levels[i + 3], point_count, &ind);
	}
	//for (size_t i = 0; i < ind; i++) { printf("%f %f\n", c[i].datetime, c[i].level); }
	printf("c is %zd, total is %zd\n", c.size(), size);
	//printf("%f, %f, %f, %f\n", c[0].datetime, c[0].level, c[1].datetime, c[1].level);
	return S_OK;
}

HRESULT CatmullRomApprox::GetLevels(floattype desiredtime, floattype stepping, size_t count,
	floattype *levels, size_t *filled, size_t derivationorder) {
	TGlucoseLevel *gl = this->levels;
	floattype time = desiredtime, x;
	int i;
	for (size_t j = 0; j < count; j++) {
		if (get_time_interval(&c[0], c.size(), time, &i) == S_FALSE) { printf("FAIL\n");levels[j] = 0.0; continue; }
		//x = time - c[i].datetime;
		switch (derivationorder) {
		case 1:
			levels[j] = d[i].level;
			break;
		default:
			levels[j] = c[i].level;
		}
		time += stepping;
		(*filled)++;
	}
	return S_OK;
}
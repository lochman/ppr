#include "CubicApprox.h"
#include "../../MaskService.h"
#include <math.h>
#include <iterator>
#include <algorithm>

HRESULT CubicApprox::Approximate(TApproximationParams *params) {
	TGlucoseLevel *levels;
	size_t size;
	mEnumeratedLevels->GetLevels(&levels);
	mEnumeratedLevels->GetLevelsCount(&size);
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
	/*
	for (size_t i = 0; i < size - 1; i++) {
		h[i] = levels[i + 1].datetime - levels[i].datetime;
	}

	float first, second;
	for (size_t i = 1; i < size - 1; i++) {
		first = 3 * (levels[i + 1].level - levels[i].level) / h[i];
		second = 3 * (levels[i].level - levels[i - 1].level) / h[i - 1];
		a[i] = first - second;
	}

	l[0] = 1.0;
	u[0] = 0.0;
	z[0] = 0.0;

	for (size_t i = 1; i < size - 1; i++) {
		l[i] = 2 * (levels[i + 1].datetime - levels[i - 1].datetime) - h[i - 1] * u[i - 1];
		u[i] = h[i] / l[i];
		z[i] = (a[i] - h[i - 1] * z[i - 1]) / l[i];
	}

	l[size - 1] = 1.0;
	z[size - 1] = 0.0;
	c[size - 1] = 0.0;

	for (int i = size - 2; i >= 0; i--) {
		//fprintf(stdout, "%d\n", i);
		c[i] = z[i] - u[i] * c[i + 1];
		b[i] = (levels[i + 1].level - levels[i].level) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3);
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
	}
	*/
	return S_OK;
}

HRESULT CubicApprox::GetLevels(floattype desiredtime, floattype stepping, size_t count,
	floattype *levels, size_t *filled, size_t derivationorder) {
	TGlucoseLevel *gl;
	size_t size;
	mEnumeratedLevels->GetLevelsCount(&size);
	mEnumeratedLevels->GetLevels(&gl);
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
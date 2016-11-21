#include"CubicApprox.h"
#include<vector>
#include <math.h>

HRESULT CubicApprox::Approximate(TApproximationParams *params) {
	TGlucoseLevel *levels;
	size_t size;

	mEnumeratedLevels->GetLevels(&levels);
	mEnumeratedLevels->GetLevelsCount(&size);
	std::vector<floattype> h(size), l(size), u(size), z(size);
	a.resize(size);
	b.resize(size);
	c.resize(size);
	d.resize(size);

	for (size_t i = 0; i < size - 1; i++) {
		h[i] = levels[i + 1].datetime - levels[i].datetime;
	}
	//fprintf(stdout, "1. for\n");
	float first, second;
	for (size_t i = 1; i < size - 1; i++) {
		first = 3 * (levels[i + 1].level - levels[i].level) / h[i];
		second = 3 * (levels[i].level - levels[i - 1].level) / h[i - 1];
		a[i] = first - second;
	}
	//fprintf(stdout, "2. for\n");
	l[0] = 1.0;
	u[0] = 0.0;
	z[0] = 0.0;

	for (size_t i = 1; i < size - 1; i++) {
		l[i] = 2 * (levels[i + 1].datetime - levels[i - 1].datetime) - h[i - 1] * u[i - 1];
		u[i] = h[i] / l[i];
		z[i] = (a[i] - h[i - 1] * z[i - 1]) / l[i];
	}
	//fprintf(stdout, "3. for\n");
	l[size - 1] = 1.0;
	z[size - 1] = 0.0;
	c[size - 1] = 0.0;

	for (int i = size - 2; i >= 0; i--) {
		//fprintf(stdout, "%d\n", i);
		c[i] = z[i] - u[i] * c [i + 1];
		b[i] = (levels[i + 1].level - levels[i].level) / (h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3);
		d[i] = (c[i + 1] - c[i]) / 3 * h[i];
	}
	//fprintf(stdout, "4. for\n");
	/*
	floattype x;
	
	for (size_t i = 0; i < size - 1; i++) {
		x = 0.002;
		printf("%f %f %f\n", levels[i].level, levels[i].level + b[i] * x + c[i] * pow(x, 2) + d[i] * pow(x, 3), levels[i + 1].level);
	}
	*/
	return S_OK;
};

HRESULT get_time_interval(TGlucoseLevel *levels, size_t size, floattype time, int *index) {
	floattype diff, last = std::abs(levels[0].datetime - time);
	for (size_t i = 0; i < size - 1; i++) {
		diff = std::abs(levels[i].datetime - time);
		if (diff > last) {
			*index = i;
			return S_OK;
		}
		last = diff;
	}
	return S_FALSE;
}

HRESULT CubicApprox::GetLevels(floattype desiredtime, floattype stepping, size_t count,
	floattype *levels, size_t *filled, size_t derivationorder) {
	TGlucoseLevel *gl;
	size_t size;
	mEnumeratedLevels->GetLevelsCount(&size);
	mEnumeratedLevels->GetLevels(&gl);
	floattype time = desiredtime, x;
	int i;
	printf("%zd %zd %zd %zd %zd\n", size, a.size(), b.size(), c.size(), d.size());
	for (size_t j = 0; j < count; j++) {
		if (get_time_interval(gl, size, time, &i) == S_FALSE) { continue; }
		x = time - gl[i].datetime;
		levels[j] = gl[i].level + b[i] * x + c[i] * pow(x, 2) + d[i] * pow(x, 3);
		printf("%f %f\n", gl[i].level, levels[j]);
		time += stepping;
		(*filled)++;
	}
	
	return S_OK;
};
#include "AkimaApprox.h"

HRESULT AkimaApprox::get_m(int i, floattype *m) {
	if (i < 0 || i > size - 2) { printf("m out of bounds: %d %d\n", i, size);return S_FALSE; }
	*m = (levels[i + 1].level - levels[i].level) / (levels[i + 1].datetime - levels[i].datetime);
	return S_OK;
}

HRESULT AkimaApprox::get_t(floattype *t) {
	floattype numerator, denominator;
	if (m.size() != 4) { printf("sizenot4\n");return S_FALSE; }
	numerator = std::abs(m[3] - m[2]) * m[1] + std::abs(m[1] - m[0]) * m[2];
	denominator = std::abs(m[3] - m[2]) + std::abs(m[1] - m[0]);
	if (denominator == 0) printf("denomis 0\n");
	*t = (denominator == 0) ? 0.5 * (m[2] + m[1]) : numerator / denominator;
	return S_OK;
}

floattype get_p2(floattype &xx, floattype &yy, floattype &ti, floattype &ti1) {
	return (3 * yy / xx - 2 * ti - ti1) / xx;
}

floattype get_p3(floattype &xx, floattype &yy, floattype &ti, floattype &ti1) {
	return (ti + ti1 - 2 * yy / xx) / std::pow(xx, 2);
}

HRESULT AkimaApprox::iterate(floattype &m_next, int i, floattype *ti) {
	floattype yy, xx, ti1;
	m.push_back(m_next);
	m.pop_front();
	if (get_t(&ti1) == S_FALSE) { return S_FALSE; }

	yy = levels[i + 1].level - levels[i].level;
	xx = levels[i + 1].datetime - levels[i].datetime;

	p1[i] = *ti;
	p2[i] = get_p2(xx, yy, *ti, ti1);
	p3[i] = get_p3(xx, yy, *ti, ti1);
	*ti = ti1;
	return S_OK;
}

HRESULT AkimaApprox::Approximate(TApproximationParams * params) {
	floattype ti, ti1, m_next;

	mEnumeratedLevels->GetLevels(&levels);
	mEnumeratedLevels->GetLevelsCount(&size);
	if (size < 4) { return S_FALSE; }

	p1.resize(size);
	p2.resize(size);
	p3.resize(size);

	if (get_m(0, &m_next) == S_FALSE) { return S_FALSE; }
	m.push_back(m_next);
	if (get_m(1, &m_next) == S_FALSE) { return S_FALSE; }
	m.push_back(m_next);
	m.push_front(2 * m[0] - m[1]); // m-1
	m.push_front(2 * m[0] - m[1]); // m-2

	if (get_t(&ti) == S_FALSE) { return S_FALSE; }
	
	for (size_t i = 0; i < size - 3; i++) {
		if (get_m(i + 2, &m_next) == S_FALSE) { return S_FALSE; }
		this->iterate(m_next, i, &ti);
	}
	m_next = 2 * m[2] - m[3];
	this->iterate(m_next, size - 3, &ti);
	m_next = 2 * m[2] - m[3];
	this->iterate(m_next, size - 2, &ti);

	return S_OK;
}

HRESULT AkimaApprox::GetLevels(floattype desiredtime, floattype stepping, size_t count,
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
			levels[j] = p1[i] + 2 * p2[i] * x + 3 * p3[i] * pow(x, 2);
			break;
		case 2:
			levels[j] = 2 * p2[i] + 6 * p3[i] * x;
			break;
		default:
			levels[j] = gl[i].level + p1[i] * x + p2[i] * pow(x, 2) + p3[i] * pow(x, 3);
		}
		time += stepping;
		(*filled)++;
	}
	return S_OK;
}
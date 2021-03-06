#include "AkimaApprox.h"
#include "../../defs.h"
#include <amp.h>
#include <amp_math.h>

floattype AkimaApprox::get_m(size_t i) {
	return (levels[i + 1].level - levels[i].level) / (levels[i + 1].datetime - levels[i].datetime);
}

floattype get_t(floattype m0, floattype m1, floattype m2, floattype m3) restrict(cpu) {
	floattype numerator, denominator;
	numerator = std::abs(m3 - m2) * m1 + std::abs(m1 - m0) * m2;
	denominator = std::abs(m3 - m2) + std::abs(m1 - m0);
	return (denominator == 0) ? 0.5 * (m2 + m1) : numerator / denominator;
}

floattype get_t(floattype m0, floattype m1, floattype m2, floattype m3) restrict(amp) {
	floattype numerator, denominator;
	numerator = concurrency::precise_math::fabs(m3 - m2) * m1 + concurrency::precise_math::fabs(m1 - m0) * m2;
	denominator = concurrency::precise_math::fabs(m3 - m2) + concurrency::precise_math::fabs(m1 - m0);
	return (denominator == 0) ? 0.5 * (m2 + m1) : numerator / denominator;
}

floattype get_p2(floattype &xx, floattype &yy, floattype &ti, floattype &ti1) restrict(cpu, amp) {
	return (3 * yy / xx - 2 * ti - ti1) / xx;
}

floattype get_p3(floattype &xx, floattype &yy, floattype &ti, floattype &ti1) restrict(cpu) {
	return (ti + ti1 - 2 * yy / xx) / std::pow(xx, 2);
}

floattype get_p3(floattype &xx, floattype &yy, floattype &ti, floattype &ti1) restrict(amp) {
	return (ti + ti1 - 2 * yy / xx) / concurrency::precise_math::pow(xx, 2);
}

HRESULT AkimaApprox::iterate(floattype &m_next, size_t i, floattype *ti) {
	floattype yy, xx, ti1;
	m.push_back(m_next);
	m.pop_front();
	ti1 = get_t(m[0], m[1], m[2], m[3]);
	yy = levels[i + 1].level - levels[i].level;
	xx = levels[i + 1].datetime - levels[i].datetime;

	p1[i] = *ti;
	p2[i] = get_p2(xx, yy, *ti, ti1);
	p3[i] = get_p3(xx, yy, *ti, ti1);
	*ti = ti1;
	return S_OK;
}

void AkimaApprox::approximate_gpu(floattype *ti) {
	int size = static_cast<int>(this->size);
	std::vector<floattype> vec_times(size), vec_lvls(size);
	floattype m_next;
	for (size_t i = 0; i < size; i++) {
		vec_times[i] = levels[i].datetime;
		vec_lvls[i] = levels[i].level;
	}

	concurrency::extent<1> ext(size - 5);
	const concurrency::array_view<const floattype, 1> times(size, vec_times), lvls(size, vec_lvls);
	concurrency::array_view<floattype, 1> p1s(size, p1), p2s(size, p2), p3s(size, p3);
	p1s.discard_data();
	p2s.discard_data();
	p3s.discard_data();
	concurrency::parallel_for_each(ext, [=](concurrency::index<1> idx) restrict(amp) {
		floattype yy, xx, ti, ti1, m0, m1, m2, m3, m4;
		const int i = idx[0] + 2; // offset
		m0 = (lvls[i - 1] - lvls[i - 2]) / (times[i - 1] - times[i - 2]);
		m1 = (lvls[i] - lvls[i - 1]) / (times[i] - times[i - 1]);
		m2 = (lvls[i + 1] - lvls[i]) / (times[i + 1] - times[i]);
		m3 = (lvls[i + 2] - lvls[i + 1]) / (times[i + 2] - times[i + 1]);
		m4 = (lvls[i + 3] - lvls[i + 2]) / (times[i + 3] - times[i + 2]);
		ti = get_t(m0, m1, m2, m3);
		ti1 = get_t(m1, m2, m3, m4);

		yy = lvls[i + 1] - lvls[i];
		xx = times[i + 1] - times[i];

		p1s[i] = ti;
		p2s[i] = get_p2(xx, yy, ti, ti1);
		p3s[i] = get_p3(xx, yy, ti, ti1);
	});
	p1s.synchronize();
	p2s.synchronize();
	p3s.synchronize();

	// calc with m-1, m-2
	m_next = get_m(2);
	this->iterate(m_next, 0, ti);
	m_next = get_m(3);
	this->iterate(m_next, 1, ti);
	// prepare for mn+1, mn+2
	m.clear();
	m.push_back(get_m(size - 5));
	m.push_back(get_m(size - 4));
	m.push_back(get_m(size - 3));
	m.push_back(get_m(size - 2));
	*ti = get_t(m[0], m[1], m[2], m[3]);
}

HRESULT AkimaApprox::Approximate(TApproximationParams * params) {
	floattype ti, m_next;
	if (size < 4) { return S_FALSE; }

	p1.resize(size);
	p2.resize(size);
	p3.resize(size);

	m.push_back(get_m(0));
	m.push_back(get_m(1));
	m.push_front(2 * m[0] - m[1]); // m-1
	m.push_front(2 * m[0] - m[1]); // m-2

	ti = get_t(m[0], m[1], m[2], m[3]);

#ifdef GPU	// GPU
	approximate_gpu(&ti);
#else // end GPU
	for (size_t i = 0; i < size - 3; i++) {
		m_next = get_m(i + 2);
		this->iterate(m_next, i, &ti);
	}
#endif
	m_next = 2 * m[2] - m[3];
	this->iterate(m_next, size - 3, &ti);
	m_next = 2 * m[2] - m[3];
	this->iterate(m_next, size - 2, &ti);
	return S_OK;
}

HRESULT AkimaApprox::GetLevels(floattype desiredtime, floattype stepping, size_t count,
	floattype *levels, size_t *filled, size_t derivationorder) {
	TGlucoseLevel *gl = this->levels;
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

std::string AkimaApprox::get_name() { return std::string("Akima spline"); }
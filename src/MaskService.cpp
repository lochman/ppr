#include "MaskService.h"
#include <cstdint>
#include <vector>
#include "defs.h"

const uint8_t bit1 = 0x80; // 1000 0000
const uint8_t bit2 = 0x40; // 0100 0000
const uint8_t bit3 = 0x20; // 0010 0000
const uint8_t bit4 = 0x10; // 0001 0000
const uint8_t bit5 = 0x08; // 0000 1000
const uint8_t bit6 = 0x04; // 0000 0100
const uint8_t bit7 = 0x02; // 0000 0010
const uint8_t bit8 = 0x01; // 0000 0001

int is_bit_set(uint8_t mask, const uint8_t bit) {
	return (mask & bit) == bit;
}

HRESULT fill_lvl_with_values(IGlucoseLevels *level, std::vector<TGlucoseLevel> &gl_levels) {
	TGlucoseLevel *lvl;
	level->SetLevelsCount(gl_levels.size());
	level->GetLevels(&lvl);
	level->AddRef();
	memcpy(lvl, gl_levels.data(), gl_levels.size() * sizeof TGlucoseLevel);
	return S_OK;
}

MaskService::MaskService(IGlucoseLevels *levels) {
	TGlucoseLevel *gl_levels;
	levels->GetLevelsCount(&segment_size);
	levels->GetLevels(&gl_levels);
	this->levels.resize(MASK_COUNT);

	for (int i = MASK_COUNT; i > 0; i--) {
		get_masked_values(gl_levels, i);
	}
}

HRESULT MaskService::get_masked_values(const TGlucoseLevel *levels, uint8_t mask) {
	std::vector<TGlucoseLevel> gl_levels;
	int bits[8];
	size_t i = 1, bit;
	
	bits[0] = is_bit_set(mask, bit1);
	bits[1] = is_bit_set(mask, bit2);
	bits[2] = is_bit_set(mask, bit3);
	bits[3] = is_bit_set(mask, bit4);
	bits[4] = is_bit_set(mask, bit5);
	bits[5] = is_bit_set(mask, bit6);
	bits[6] = is_bit_set(mask, bit7);
	bits[7] = is_bit_set(mask, bit8);

	gl_levels.push_back(levels[0]); // always add first point

	while (i < segment_size - 1) {
		for (bit = 0; bit < 8; bit++) {
			if (bits[bit] && i < segment_size - 1) {
				gl_levels.push_back(levels[i]);
			}
			i++;
			if (i >= segment_size - 1) { break; }
		}
	}
	gl_levels.push_back(levels[segment_size - 1]); // always add last point
	this->levels[mask - 1] = new CGlucoseLevels();
	fill_lvl_with_values(this->levels[mask - 1], gl_levels);
	return S_OK;
}

HRESULT MaskService::get_mask(IGlucoseLevels **levels, uint8_t mask) {
	*levels = this->levels[mask - 1];
	return S_OK;
}

HRESULT MaskService::get_inverse_mask(IGlucoseLevels **levels, uint8_t mask) {
	if (mask == MASK_COUNT) { get_mask(levels, mask); }
	else { get_mask(levels, ~mask); }
	return S_OK;
}

MaskService::~MaskService() {
	for (size_t i = 0; i < levels.size(); i++) { if (levels[i] != NULL) levels[i]->Release(); }
}

HRESULT MaskService::get_segment_size(size_t *size) {
	*size = this->segment_size;
	return S_OK;
}

HRESULT MaskService::get_mask_count(size_t *size) {
	*size = levels.size();
	return S_OK;
}
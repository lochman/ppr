#include "MaskService.h"
#include <cstdint>
#include <vector>

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

void fill_lvl_with_values(IGlucoseLevels *level, std::vector<TGlucoseLevel> &gl_levels) {
	TGlucoseLevel *lvl;
	level->AddRef();
	level->SetLevelsCount(gl_levels.size());
	level->GetLevels(&lvl);
	memcpy(lvl, gl_levels.data(), gl_levels.size() * sizeof TGlucoseLevel);
}

MaskService::MaskService(TGlucoseLevel *levels, size_t const &size) : levels(levels), size(size) {
	std::vector<TGlucoseLevel> gl_levels;
	
	for (int i = 255; i > 0; i--) {
		get_masked_values(gl_levels, i);
		fill_lvl_with_values(&masks[i - 1], gl_levels);
		gl_levels.clear();
	}
}

void MaskService::get_masked_values(std::vector<TGlucoseLevel> &glucose_levels, uint8_t mask) const {
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

	glucose_levels.push_back(levels[0]); // always add first point

	while (i < size - 1) {
		for (bit = 0; bit < 8; bit++) {
			if (bits[bit] && i < size - 1) {
				glucose_levels.push_back(levels[i]);
			}
			i++;
			if (i >= size - 1) { break; }
		}
	}
	
	glucose_levels.push_back(levels[size - 1]); // always add last point
	//printf("Got mask with size %zd\n", glucose_levels.size());
}

void MaskService::get_mask(IGlucoseLevels **levels, uint8_t mask) {
	*levels = &masks[mask - 1];
	/*
	std::vector<TGlucoseLevel> gl_levels;
	IGlucoseLevels *level = new CGlucoseLevels();
	get_masked_values(gl_levels, mask);
	fill_lvl_with_values(level, gl_levels);
	*levels = level;
	*/
	/*
	size_t m;
	level->GetLevelsCount(&m);
	TGlucoseLevel *lvl;
	level->GetLevels(&lvl);

	printf("getting mask %d size is %zd\n", mask, m);
	printf("last is %f\n", lvl[m - 1].datetime);
	printf("222 %f, %zd\n", gl_levels[gl_levels.size() - 1].datetime, gl_levels.size());
	*/
}

void MaskService::get_inverse_mask(IGlucoseLevels **levels, uint8_t mask) {
	get_mask(levels, ~mask);
}

void MaskService::get_levels(TGlucoseLevel **levels) const {
	*levels = this->levels;
}

void MaskService::get_levels_size(size_t *size) const {
	*size = this->size;
}
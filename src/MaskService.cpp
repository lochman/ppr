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

MaskService::MaskService(IGlucoseLevels *mEnumeratedLevels) {
	mEnumeratedLevels->GetLevels(&levels);
	mEnumeratedLevels->GetLevelsCount(&size);
}

void MaskService::get_masked_values(std::vector<TGlucoseLevel *> *glucose_levels, uint8_t mask) {
	int bits[8];
	size_t i = 0, bit;
	
	bits[0] = is_bit_set(mask, bit1);
	bits[1] = is_bit_set(mask, bit2);
	bits[2] = is_bit_set(mask, bit3);
	bits[3] = is_bit_set(mask, bit4);
	bits[4] = is_bit_set(mask, bit5);
	bits[5] = is_bit_set(mask, bit6);
	bits[6] = is_bit_set(mask, bit7);
	bits[7] = is_bit_set(mask, bit8);

	while (i < size) {
		for (bit = 0; bit < 8; bit++) {
			if (bits[bit] && i < size) {
				glucose_levels->push_back(&levels[i]);
			}
			i++;
			if (i >= size) { return; } // all values have been masked
		}
	}
}

void MaskService::get_inverse_masked_values(std::vector<TGlucoseLevel *> *glucose_levels, uint8_t mask) {
	get_masked_values(glucose_levels, ~mask);
}
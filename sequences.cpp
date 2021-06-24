#include <cstring>
#include <cassert>
#include <algorithm>

#include "sequences.hpp"



using namespace std;


/* Bitshift to the left all the bits in the array with a maximum of 7 bits.
 * Overflow on the left will be set into the previous cell.
 */
void leftshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	for (uint64_t i=0 ; i<length-1 ; i++) {
		bitarray[i] = (bitarray[i] << bitshift) | (bitarray[i+1] >> (8-bitshift));
	}
	bitarray[length-1] <<= bitshift;
}

/* Similar to the previous function but on the right */
void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	for (uint64_t i=length-1 ; i>0 ; i--) {
		bitarray[i] = (bitarray[i-1] << (8-bitshift)) | (bitarray[i] >> bitshift);
	}
	bitarray[0] >>= bitshift;
}

/* Fusion to bytes into one.
 * The merge_index higher bits are from left_bits the others from right_bits
 */
uint8_t fusion8(uint8_t left_bits, uint8_t right_bits, size_t merge_index) {
	uint8_t mask = 0xFF << (8-merge_index);
	return (left_bits & mask) | (right_bits & ~mask);
}


uint KffSeqStream::next_sequence(uint8_t * & seq, uint8_t * & data) {
	if (this->reader.has_next()) {
		return this->reader.next_block(seq, data);
	}

	return 0;
}




// #include <iostream>
void subsequence(const uint8_t * sequence, const uint seq_size, uint8_t * extracted, const uint begin_nucl, const uint end_nucl) {
	// Extract the correct slice
	uint seq_left_offset = (4 - seq_size % 4) % 4;
	uint extract_start_byte = (seq_left_offset + begin_nucl) / 4;
	uint extract_stop_byte = (seq_left_offset + end_nucl) / 4;

	memcpy(extracted, sequence + extract_start_byte, extract_stop_byte - extract_start_byte + 1);

	// Align the bits
	uint extract_left_offset = (seq_left_offset + begin_nucl) % 4;
	uint extract_right_offset = (seq_size - end_nucl - 1) % 4;
	

	if (extract_right_offset < 4 - extract_left_offset) {
		rightshift8(extracted, extract_stop_byte - extract_start_byte + 1, extract_right_offset * 2);
	} else {
		leftshift8(extracted, extract_stop_byte - extract_start_byte + 1, (4 - extract_right_offset) * 2);
	}
}


vector<uint64_t> compute_mini_candidates(const uint8_t * seq, const uint size, const uint k, const uint m) {
	vector<uint64_t> candidates(size - m + 1);

	uint offset = (4 - (size % 4)) % 4;
	uint64_t current_value = 0;
	// Compute prefix of first candidate
	for (uint i=0 ; i<m-1 ; i++) {
		uint idx = offset + i;
		uint byte_idx = idx/4;
		uint nucl_shift = 3 - (idx % 4);

		uint nucl = (seq[byte_idx] >> (nucl_shift * 2)) & 0b11;
		current_value = (current_value << 2) + nucl;
	}

	// Compute minimizer candidates
	uint64_t m_mask = (1 << (m*2)) - 1;
	for (uint i=m-1 ; i<size ; i++) {
		uint idx = offset + i;
		uint byte_idx = idx/4;
		uint nucl_shift = 3 - (idx % 4);

		uint nucl = (seq[byte_idx] >> (nucl_shift * 2)) & 0b11;
		current_value = ((current_value << 2) + nucl) & m_mask;	
		candidates[i - m + 1] = current_value;
	}

	return candidates;
}

vector<pair<int, uint64_t> > compute_minizers(const uint8_t * seq, const uint size, const uint k, const uint m) {
	vector<pair<int, uint64_t> > minimizers;
	// Get all the candidates
	vector<uint64_t> candidates = compute_mini_candidates(seq, size, k, m);

	int prev_pos = k + 2;
	// Compute the minimizer of each sliding window of size k - m
	for (uint i=0 ; i<=size-k ; i++) {
		auto smallest = min_element(candidates.begin()+i, candidates.begin()+i+(k-m));
		int pos = smallest - candidates.begin();
		// New minimizer ?
		if (pos != prev_pos) {
			prev_pos = pos;
			minimizers.emplace_back(pos, *smallest);
		}
	}

	return minimizers;
}

std::vector<pair<uint, uint> > compute_skmers(const uint seq_size, const uint k, const uint m, std::vector<std::pair<int, uint64_t> > & minimizers) {
	// Superkmer list
	vector<pair<uint, uint> > skmers;

	uint current_begin = 0;
	for (uint i=0 ; i<minimizers.size()-1 ; i++) {
		pair<int, uint64_t> & current_mini = minimizers[i];
		pair<int, uint64_t> & next_mini = minimizers[i+1];

		uint end = 0;
		uint new_begin = 0;
		// First minimizer is dominant
		if (current_mini.second <= next_mini.second) {
			end = current_mini.first + k - 1;
			new_begin = current_mini.first + 1;
		}
		// Second minimizer is dominant
		else {
			new_begin = next_mini.first - k + m;
			end = current_mini.first + m - 1 + (new_begin - current_begin);
		}

		// Add the superkmer and save position
		skmers.emplace_back(current_begin, end);
		current_begin = new_begin;
	}

	// Save the last skmer
	skmers.emplace_back(current_begin, seq_size-1);

	return skmers;
}

std::vector<pair<uint, uint> > compute_skmers(const uint8_t * seq, const uint size, const uint k, const uint m) {
	std::vector<std::pair<int, uint64_t> > minimizers = compute_minizers(seq, size, k, m);

	return compute_skmers(size, k, m, minimizers);
}


void search_mini(uint8_t * seq, const uint size, const uint m, uint & minimizer, uint & minimizer_position) {
	// Datastructure prepare
	uint k_bytes = (size + 3) / 4;
	uint m_bytes = (m + 3) / 4;
	uint8_t * bin_copy = new uint8_t[k_bytes];
	// Mask to cover all the minimizer bytes except the highest one
	uint low_mask = (1 << (8 * (m_bytes - 1))) - 1;
	// Mask to cover usefull bits of the higher byte of the minimizer
	uint high_mask = (1 << (2 * (((m-1) % 4) + 1))) - 1;


	// Minimizer prepare (Do not use memcpy for endianess problems !!)
	minimizer = 0;
	for (uint i=0 ; i<m_bytes-1 ; i++) {
		minimizer += ((uint)seq[k_bytes - 1 - i]) << (8 * i);
	}
	minimizer += ((uint)seq[k_bytes - 1 - (m_bytes - 1)] & high_mask) << (8 * (m_bytes - 1));
	uint mini_candidate = minimizer;

	// Forward search
	memcpy(bin_copy, seq, k_bytes);
	for (int m_idx=size-m ; m_idx>=0 ; m_idx--) {
		// Update minimizer
		if (mini_candidate <= minimizer) {
			minimizer = mini_candidate;
			minimizer_position = m_idx;
		}

		// Shift everything
		rightshift8(bin_copy, k_bytes, 2);
		mini_candidate >>= 2;
		// remove first byte
		mini_candidate &= low_mask;
		// Update first byte
		mini_candidate += ((uint)bin_copy[k_bytes - 1 - (m_bytes - 1)] & high_mask) << (8 * (m_bytes - 1));
	}

	delete[] bin_copy;
}


uint64_t seq_to_uint(const uint8_t * seq, uint seq_size) {
	uint nucl_to_extract = seq_size;
	if (nucl_to_extract > 32)
		nucl_to_extract = 32;

	uint seq_offset = (4 - (seq_size % 4)) % 4;
	uint seq_bytes = (seq_size + 3) / 4;
	uint useless_seq_nucl = seq_size - nucl_to_extract;

	uint suff_offset = (4 - (nucl_to_extract % 4)) % 4;
	uint mask = (1 << (2 * (4 - suff_offset))) - 1;
	uint suff_first_byte = (seq_offset + useless_seq_nucl) / 4;

	uint64_t val = seq[suff_first_byte] & mask;
	for (uint idx=suff_first_byte+1 ; idx<seq_bytes ; idx++) {
		val <<= 8;
		val += seq[idx];
	}

	return val;
}


void uint_to_seq(uint seq, uint8_t * bin_seq, uint size) {
	uint seq_bytes = (size + 3) / 4;

	for (int idx=seq_bytes-1 ; idx>=0 ; idx--) {
		bin_seq[idx] = seq & 0b11111111;
		seq >>= 8;
	}
}

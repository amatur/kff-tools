#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <queue>

#include "encoding.hpp"
#include "sequences.hpp"
#include "compact.hpp"
#include "merge.hpp"


using namespace std;


Compact::Compact() {
	input_filename = "";
	output_filename = "";

	this->buffer_size = 1 << 10;
	this->next_free = 0;
	this->kmer_buffer = (uint8_t *)malloc(this->buffer_size);
	memset(this->kmer_buffer, 0, this->buffer_size);
}


Compact::~Compact() {
	free(this->kmer_buffer);
}


void Compact::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("compact", "Read a kff file and try to compact the kmers from minimizer sections. The available ram must be sufficent to load a complete minimizer section into memory.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", input_filename, "Input kff file to compact.");
	input_option->required();
	input_option->check(CLI::ExistingFile);
	CLI::Option * out_option = subapp->add_option("-o, --outfile", output_filename, "Kff to write (must be different from the input)");
	out_option->required();
}


void Compact::exec() {
	Kff_file infile(input_filename, "r");
	Kff_file outfile(output_filename, "w");

	outfile.write_encoding(infile.encoding);
	
	outfile.set_uniqueness(infile.uniqueness);
	outfile.set_canonicity(infile.canonicity);
	
	// Metadata transfer
  uint8_t * metadata = new uint8_t[infile.metadata_size];
  infile.read_metadata(metadata);
  outfile.write_metadata(infile.metadata_size, metadata);
  delete[] metadata;

	while (infile.tellp() != infile.end_position) {
		char section_type = infile.read_section_type();

		if (section_type == 'v') {
			Section_GV isgv(&infile);
			Section_GV osgv(&outfile);
			for (auto & p : isgv.vars)
				osgv.write_var(p.first, p.second);
			isgv.close();
			osgv.close();
		}
		else if (section_type == 'i') {
			Section_Index si(&infile);
			si.close();
		}
		else if (section_type == 'r') {
			Section_Raw sr(&infile);
			sr.close();
		}
		else if (section_type == 'm') {
			uint k = outfile.global_vars["k"];
			uint m = outfile.global_vars["m"];

			// Rewrite a value section if max is not sufficently large
			if (outfile.global_vars["max"] < (k-m) * 2) {
				unordered_map<string, uint64_t> values(outfile.global_vars);
				Section_GV sgv(&outfile);

				for (auto & p : values)
					if (p.first != "max")
						sgv.write_var(p.first, p.second);
				sgv.write_var("max", (k-m)*2);

				sgv.close();
			}

			// Compact and save the kmers
			Section_Minimizer sm(&infile);
			this->compact_section(sm, outfile);
			sm.close();
		}
	}

	infile.close();
	outfile.close();
}

void Compact::compact_section(Section_Minimizer & ism, Kff_file & outfile) {
	// General variables
	uint k = outfile.global_vars["k"];
	uint m = outfile.global_vars["m"];
	uint data_size = outfile.global_vars["data_size"];
	
	// 1 - Load the input section
	vector<vector<long> > kmers_per_index = this->prepare_kmer_matrix(ism);
	
	// 2 - Compact kmers
	vector<pair<uint8_t *, uint8_t *> > to_compact = this->greedy_assembly(kmers_per_index, k, m);
	vector<vector<uint8_t *> > paths = this->pairs_to_paths(to_compact);

	Section_Minimizer osm(&outfile);
	osm.write_minimizer(ism.minimizer);
	this->write_paths(paths, osm, k, m, data_size);
	osm.close();

	// Cleaning
}

// this->buffer_size = 1 << 10;
// 	this->next_free = 0;
// 	this->kmer_buffer

vector<vector<long> > Compact::prepare_kmer_matrix(Section_Minimizer & sm) {
	vector<vector<long> > kmer_matrix;
	kmer_matrix.resize(sm.k - sm.m + 1);
	
	uint64_t max_nucl = sm.k + sm.max - 1;
	uint64_t max_seq_bytes = (max_nucl + 3) / 4;
	uint64_t kmer_bytes = (sm.k - sm.m + 3) / 4;
	uint64_t mini_pos_size = (static_cast<uint>(ceil(log2(max_nucl))) + 7) / 8;

	uint8_t * seq_buffer = new uint8_t[max_seq_bytes];
	uint8_t * data_buffer = new uint8_t[sm.data_size * sm.max];

	// 1 - Load the input section
	for (uint n=0 ; n<sm.nb_blocks ; n++) {
		uint64_t mini_pos = 0xFFFFFFFFFFFFFFFF;
		// Read sequence
		uint nb_kmers = sm.read_compacted_sequence_without_mini(
			seq_buffer, data_buffer, mini_pos);

		// Add kmer by index
		for (uint kmer_idx=0 ; kmer_idx<nb_kmers ; kmer_idx++) {
			uint kmer_pos = sm.k - (uint)sm.m - mini_pos + kmer_idx;

			// Realloc if needed
			if (this->buffer_size - this->next_free < kmer_bytes + sm.data_size + mini_pos_size) {
				this->kmer_buffer = (uint8_t *) realloc((void *)this->kmer_buffer, this->buffer_size*2);
				memset(this->kmer_buffer + this->buffer_size, 0, this->buffer_size);
				this->buffer_size *= 2;
			}

			// Copy kmer sequence
			subsequence(seq_buffer, sm.k - sm.m + nb_kmers - 1, this->kmer_buffer + next_free, kmer_idx, kmer_idx + sm.k - sm.m - 1);
			// Copy data array
			memcpy(this->kmer_buffer + next_free + kmer_bytes, data_buffer + kmer_idx * sm.data_size, sm.data_size);
			// Write mini position
			uint kmer_mini_pos = mini_pos - kmer_idx;
			for (int b=mini_pos_size-1 ; b>=0 ; b--) {
				*(this->kmer_buffer + next_free + kmer_bytes + sm.data_size + b) = kmer_mini_pos & 0xFF;
				kmer_mini_pos >>= 8;
			}
			// Update
			kmer_matrix[kmer_pos].push_back(this->next_free);
			next_free += kmer_bytes + sm.data_size + mini_pos_size;
		}
	}

	delete[] seq_buffer;
	delete[] data_buffer;

	return kmer_matrix;
}

vector<pair<uint8_t *, uint8_t *> > Compact::greedy_assembly(vector<vector<long> > & positions, const uint k, const uint m) {
	uint nb_kmers = k - m + 1;
	vector<pair<uint8_t *, uint8_t *> > assembly;

	// Index kmers from the 0th set
	for (long kmer_pos : positions[0])
		assembly.emplace_back(nullptr, kmer_buffer + kmer_pos);

	for (uint i=0 ; i<nb_kmers-1 ; i++) {
		// Index kmers in ith set
		unordered_map<uint64_t, vector<uint8_t *> > index;
		
		for (long kmer_pos : positions[i]) {
			uint8_t * kmer = kmer_buffer + kmer_pos;
			// Get the suffix
			uint64_t val = subseq_to_uint(kmer, nb_kmers-1, 1, nb_kmers-2);
			// Add a new vector for this value
			if (index.find(val) == index.end())
				index[val] = vector<uint8_t *>();
			// Add the kmer to the value list
			index[val].push_back(kmer);
		}

		// link kmers from (i+1)th set to ith kmers.
		for (long kmer_pos : positions	[i+1]) {
			uint8_t * kmer = kmer_buffer + kmer_pos;
			uint64_t val = subseq_to_uint(kmer, nb_kmers-1, 0, nb_kmers-3);

			if (index.find(val) == index.end()) {
				// No kmer available for matching
				assembly.emplace_back(nullptr, kmer);
			} else {
				bool chaining_found = false;
				uint candidate_pos = 0;
				// verify complete matching for candidates kmers
				for (uint8_t * candidate : index[val]) {
					// If the kmers can be assembled
					if (sequence_compare(
								kmer, nb_kmers-1, 0, nb_kmers-3,
								candidate, nb_kmers-1, 1, nb_kmers-2
							) == 0) {
						// Update status
						chaining_found = true;
						assembly.emplace_back(candidate, kmer);

						// remove candidate from list
						index[val].erase(index[val].begin()+candidate_pos);
						// Quit candidate searching
						break;
					}

					candidate_pos += 1;
				}
				// If no assembly possible, create a new superkmer
				if (not chaining_found)
					assembly.emplace_back(nullptr, kmer);
			}
		}
	}

	// Index last kmers without compaction
	int assembly_idx = assembly.size()-1;
	for (auto it=positions[nb_kmers-1].end() ; it>positions[nb_kmers-1].begin() ; it--) {
		long pos = *(it-1);
		uint8_t * kmer = kmer_buffer + 	pos;

		if (assembly_idx < 0 or kmer != assembly[assembly_idx].second) {
			assembly.emplace_back(nullptr, kmer);
		} else {
			assembly_idx--;
		}
	}

	return assembly;
}

vector<vector<uint8_t *> > Compact::pairs_to_paths(const vector<pair<uint8_t *, uint8_t *> > & to_compact) {
	vector<vector<uint8_t *> > paths;
	unordered_map<uint8_t *, uint> path_registry;

	for (const pair<uint8_t *, uint8_t *> & p : to_compact) {
		// First element of a compaction path
		if (p.first == nullptr) {
			path_registry[p.second] = paths.size();
			paths.emplace_back(vector<uint8_t *>());

			uint vec_idx = path_registry[p.second];
			paths[vec_idx].push_back(p.second);
		}
		// Extending existing path
		else {
			uint vec_idx = path_registry[p.first];
			paths[vec_idx].push_back(p.second);
			path_registry.erase(p.first);
			path_registry[p.second] = vec_idx;
		}
	}

	return paths;
}

void Compact::write_paths(const vector<vector<uint8_t *> > & paths, Section_Minimizer & sm, const uint k, const uint m, const uint data_size) {
	uint kmer_bytes = (k - m + 3) / 4;
	uint kmer_offset = (4 - ((k - m) % 4)) % 4;
	uint mini_pos_size = (static_cast<uint>(ceil(log2(k - m + 1))) + 7) / 8;

	uint max_skmer_bytes = (2 * (k - m) + 3) / 4;
	uint8_t * skmer_buffer = new uint8_t[max_skmer_bytes + 1];
	uint data_bytes = (k - m + 1) * data_size;
	uint8_t * data_buffer = new uint8_t[data_bytes];

	// uint8_t encoding[] = {0, 1, 3, 2};
	// Stringifyer strif(encoding);

	// Write skmer per skmer
	for (const vector<uint8_t *> & path : paths) {
		// Cleaning previous skmers/data
		memset(skmer_buffer, 0, max_skmer_bytes + 1);
		memset(data_buffer, 0, data_bytes);

		// Get the skmer minimizer position
		uint mini_pos = 0;
		uint8_t * mini_pos_pointer = path[0] + kmer_bytes + data_size;
		for (uint b=0 ; b<mini_pos_size ; b++) {
			mini_pos <<= 8;
			mini_pos += *mini_pos_pointer;
			mini_pos_pointer += 1;
		}

		// Usefull variables
		uint skmer_size = k - m - 1 + path.size();
		// uint skmer_bytes = (skmer_size + 3) / 4;
		uint skmer_offset = (4 - (skmer_size % 4)) % 4;

		// Save the first kmer
		// cout << "+" << strif.translate(path[0], k-m);
		memcpy(skmer_buffer, path[0], kmer_bytes);
		leftshift8(skmer_buffer, kmer_bytes, 2 * kmer_offset);
		rightshift8(skmer_buffer, kmer_bytes + 1, 2 * skmer_offset);
		// Save the first data
		memcpy(data_buffer, path[0] + kmer_bytes, data_size);

		// Compact kmer+data one by one
		for (uint kmer_idx = 1 ; kmer_idx<path.size() ; kmer_idx++) {
			uint8_t * kmer = path[kmer_idx];

			// cout << " -> " << strif.translate(kmer, k-m);
			// Compute compaction position
			uint compact_nucl_pos = skmer_offset + k - m - 1 + kmer_idx;
			uint compact_byte = compact_nucl_pos / 4;
			uint compact_shift = 3 - (compact_nucl_pos % 4);
			// Compact the nucleotide
			uint8_t last_nucl = kmer[kmer_bytes - 1] & 0b11;
			skmer_buffer[compact_byte] |= last_nucl << (2 * compact_shift);
			// Copy data
			memcpy(data_buffer + kmer_idx * data_size, path[kmer_idx] + kmer_bytes, data_size);
		}
		// cout << endl;

		// Write everything in the file
		sm.write_compacted_sequence_without_mini(skmer_buffer, skmer_size, mini_pos, data_buffer);
		// cout << " " << strif.translate(skmer_buffer, skmer_size) << endl;
	}
	// cout << endl;

	delete[] skmer_buffer;
	delete[] data_buffer;
}


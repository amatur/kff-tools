#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "spss_to_kff.hpp"
#include "encoding.hpp"
#include "merge.hpp"
#include "sequences.hpp"


using namespace std;



SpssToKff::SpssToKff() {
	basename = "";
	output_filename = "";
	data_size = 0;
	k = 0;
	max_kmerseq = 255;
	delimiter = " ";
	data_delimiter = ",";
	minimizer_size=10;
}

void SpssToKff::cli_prepare(CLI::App * app) {
	this->subapp = app->add_subcommand("SpssToKff", "Convert a text kmer file or a text sequence file into a kff file. Kmers or sequences must be 1 per line. If data size is more than 0, then the delimiters are used to split each line.");
	CLI::Option * input_option = subapp->add_option("-i, --infile", basename, "basename .instr: A text file with one sequence per line (sequence omitted if its size < k). Empty data is added (size defined by -d option).");
	input_option->required();
	//input_option->check(CLI::ExistingFile);
	// CLI::Option * output_option = subapp->add_option("-o, --outfile", output_filename, "The kff output file name.");
	// output_option->required();
	CLI::Option * k_opt = subapp->add_option("-k, --kmer-size", k, "Mandatory kmer size");
	k_opt->required();
	CLI::Option * max_kmerseq_opt = subapp->add_option("-m, --max-kmer-seq", max_kmerseq, "The maximum number of kmer that can be inside of sequence in the output (default 255).");
	max_kmerseq_opt->required();

	CLI::Option * mini_size_opt = subapp->add_option("-s, --minimizer-size", minimizer_size, "Minimizer size");
	mini_size_opt->required();

	CLI::Option * data_size_opt = subapp->add_option("-d, --data-size", data_size, "Data size in Bytes (Default 0, max 8).");
	data_size_opt->required();
	//subapp->add_option("--delimiter", this->delimiter, "Character used as a delimiter between the sequence and the data (default ' ').");
	//subapp->add_option("--data-delimiter", this->data_delimiter, "Character used as a delimiter between two kmers data from the same sequence (default ',').");
}


/** Read txt sequence file (1 sequence per line)
 */
class TxtSeqStream : public SequenceStream {
private:
  std::fstream fs;
  Binarizer bz;

  uint buffer_size;
  uint8_t * seq_buffer;
  uint8_t * data_buffer;

  uint k;
  uint data_size;
  bool data_present;

  string delimiter;
  string data_delimiter;

public:
  TxtSeqStream(const std::string filename, const uint8_t encoding[4], uint k, uint data_size, string delim, string data_delim) 
      : fs(filename, std::fstream::in)
      , bz(encoding)
      , buffer_size(1024)
      , seq_buffer(new uint8_t[1024])
      , data_buffer(new uint8_t[(1024 * 4 - k + 1	) * data_size])
      , k(k)
      , data_size(data_size)// uint8_t * seq;
	// //uint8_t * sub_seq = new uint8_t[(max_kmerseq + 3) / 4];
	// uint8_t * data = new uint8_t[data_size * max_kmerseq];

	//uint8_t * byte_seq;
	// Write the minimizer
	//byte_seq = encode("AAATAACACA");
	//sm.write_minimizer(byte_seq);
	//	delete [] byte_seq;


	// std::fstream pos_file("~/projects/Blight/example/m_88095.pos", std::ios_base::in);
    // uint64_t minimizer_pos;
      , delimiter(delim)
      , data_delimiter(data_delim)
  {};
  ~TxtSeqStream() {
    this->fs.close();
    delete[] this->seq_buffer;
    delete[] this->data_buffer;
  }

	uint next_sequence(uint8_t * & seq, uint8_t * & data) {
		// Verify stream integrity
		if (not this->fs)
			return 0;

		// read next sequence
		string line;
		getline(this->fs, line);
		if (line.size() == 0)
			return 0;

		// Get the split limit
		size_t seq_size = 0;
		
		if (this->data_size == 0)
			seq_size = line.size();
		else {
			seq_size = line.find(this->delimiter);
			
			if (seq_size == string::npos) {
				cerr << "Delimiter not found in" << endl << "\t" << line << endl;
				exit(1);
			}
		}

		uint nb_kmers = seq_size - k + 1;

		// Update buffers
		if (seq_size > this->buffer_size * 4) {
			delete[] this->seq_buffer;
			this->buffer_size = (seq_size + 3) / 4;
			this->seq_buffer = new uint8_t[this->buffer_size];

			delete[] this->data_buffer;
			this->data_buffer = new uint8_t[(seq_size - k + 1) * this->data_size];
		}
		// convert/copy the sequence
		this->bz.translate(line, seq_size, this->seq_buffer);
		seq = this->seq_buffer;

		// If counts
		if (this->data_size > 0) {
			char * str_data = (char *)line.c_str() + seq_size + 1;

			for (uint i=0 ; i<nb_kmers ; i++) {
				unsigned long count = strtoul(str_data, &str_data, 10);
				str_data += 1;
				for (uint d=0 ; d<this->data_size ; d++) {
					data[data_size * i + (data_size-1-d)] = (uint8_t)count & 0xFF;
					count >>= 8;
				}
			}
		}

		return seq_size;
	}

	int next_sequence(uint8_t * & seq, uint max_seq_size, uint8_t * & data, uint max_data_size) {
		return this->next_sequence(seq, data);
  }
};

void SpssToKff::exec() {

	this->output_filename = basename+".ess.kff";
	Kff_file outfile(this->output_filename, "w");
	string minimizer_filename=basename+".minimizer";
	string position_filename=basename+".pos";
	string input_filename = basename+".instr";
	string encoding_filename = basename+".encoding";
	
	//uint minimizer_size = 10;

		// Write needed variables
	Section_GV sgv(&outfile);
	sgv.write_var("k", this->k);
	sgv.write_var("m", minimizer_size);
	sgv.write_var("data_size", this->data_size);
	//sgv.write_var("ordered", 0);
	sgv.write_var("max", this->max_kmerseq);
	sgv.close();


	std::fstream myfile(encoding_filename, std::ios_base::in);
    int a, b, c, d;
    myfile >> a >> b >> c >> d;

	const uint8_t encoding[4] = {uint8_t(a), uint8_t(b), uint8_t(c), uint8_t(d)};

	// Write the sequences inside of a minimizer section
	Section_Minimizer sm(&outfile);

	std::fstream pos_file(position_filename, std::ios_base::in);
    uint64_t minimizer_pos;
    
	uint8_t * seq_mini;
	uint8_t * data_mini = new uint8_t[0];
	uint seq_size = 0;

	TxtSeqStream stream_mini(minimizer_filename, encoding, minimizer_size, 0, " ", ",");
	while ((seq_size = stream_mini.next_sequence(seq_mini, data_mini)) > 0) {
		sm.write_minimizer(seq_mini);
		break;
	}
	
	uint8_t * seq;
	uint8_t * data = new uint8_t[this->data_size *  this->max_kmerseq];
	TxtSeqStream stream(input_filename, encoding, this->k ,this->data_size,  " ", ",");
	while ((seq_size = stream.next_sequence(seq, data)) > 0) {
		if (seq_size < k){
			cerr<<"Error: seq size is smaller than k"<<endl;
			exit(2);
		}
		uint nb_kmers = seq_size - this->k + 1;
		// Full sequence copy
		if (nb_kmers > this->max_kmerseq) {
			cerr<<"Error: nb kmers exceed max "<<nb_kmers<<">"<<this->max_kmerseq<<endl;
			exit(2);
		}

		if (!(pos_file >> minimizer_pos))
		{
			cerr<<"Error: Invalid position file"<<endl;
			exit(2);
		}
		
		uint nucl_size =  this->k + nb_kmers - 1;
		sm.write_compacted_sequence(seq, nucl_size, minimizer_pos, data);
	}
	sm.close();
	outfile.close();
}


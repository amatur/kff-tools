#include <string>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kfftools.hpp"


#ifndef SpssToKff_H
#define SpssToKff_H

class SpssToKff: public KffTool {
private:
	std::string basename;
	std::string output_filename;

	uint data_size;
	uint k;
	uint max_kmerseq;
	uint minimizer_size;

	std::string delimiter;
  std::string data_delimiter;

	void monofile();
	void multifile();

public:
	SpssToKff();
	void cli_prepare(CLI::App * subapp);
	void exec();
	void exec2();
};

#endif
//
// Created by zsc on 2023/11/6.
//

#ifndef RUNE_FASTA_IO_H
#define RUNE_FASTA_IO_H
#include <algorithm>
#include <cstring>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <zlib.h>

class fasta_io {
private:
  std::unordered_set<char> bases = {'A', 'T', 'G', 'C', 'N'};
  std::string file_name;
  bool check_seq(const std::string &seq) {
    return std::all_of(seq.begin(), seq.end(), [this](const char base) {
      return bases.find(base) != bases.end();
    });
  }
  std::unordered_map<std::string, std::string> read_text();
  std::unordered_map<std::string, std::string> read_gz();

public:
  explicit fasta_io(std::string file_name) {
    this->file_name = std::move(file_name);
  }
  // read fasta file and return an unordered_name:
  // sequence_id => sequence
  std::unordered_map<std::string, std::string> read();
};

#endif // RUNE_FASTA_IO_H

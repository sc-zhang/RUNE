//
// Created by zsc on 2023/11/11.
//

#ifndef RUNE_DUMPER_H
#define RUNE_DUMPER_H
#include <utility>

#include "bin_io.h"
#include "fasta_io.h"
#include "k_bin.h"
#include "msg.h"
#include "rune.h"

class dumper {
private:
  std::string fasta_file;
  std::string kmer_file;
  uint8_t k_size = 0;
  std::unordered_map<std::string, std::string> mp_seq;
  std::unordered_map<uint64_t, std::pair<uint32_t, uint64_t>> mp_kmer;
  std::unordered_map<std::string, uint32_t> id_sample;
  std::unordered_map<uint32_t, std::string> sample_id;

public:
  explicit dumper(std::string fasta_file, uint8_t k_size,
                  std::string kmer_file) {
    this->fasta_file = std::move(fasta_file);
    this->k_size = k_size;
    this->kmer_file = std::move(kmer_file);
  }
  void extract();
  void save();
  std::unordered_map<uint64_t, std::pair<uint32_t, uint64_t>> get_kmer_db();
  uint32_t get_sample_id(const std::string &sample);
  std::string get_sample_name(const uint32_t &id);
};

#endif // RUNE_DUMPER_H

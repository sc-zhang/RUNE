//
// Created by zsc on 2023/11/11.
//

#ifndef RUNE_LOADER_H
#define RUNE_LOADER_H
#include "bin_io.h"
#include "k_bin.h"
#include "msg.h"
#include "rune.h"
#include <fstream>

class loader {
private:
  std::string kmer_file;
  uint8_t k_size;
  std::unordered_map<uint64_t, std::pair<uint32_t, uint64_t>> mp_kmer;
  std::unordered_map<std::string, uint32_t> id_sample;
  std::unordered_map<uint32_t, std::string> sample_id;
  bin_io kmer_bin_io;

public:
  loader(std::string kmer_bin_file, uint8_t k_size)
      : kmer_bin_io(kmer_bin_file) {
    this->kmer_file = std::move(kmer_bin_file);
    this->k_size = k_size;
  };
  void load();
  void save(const std::string &output_file);
  std::unordered_map<uint64_t, uint64_t> get_kmer_db() const;
  uint32_t get_sample_id(const std::string &sample);
  std::string get_sample_name(const uint32_t &id);
};

#endif // RUNE_LOADER_H

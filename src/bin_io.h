//
// Created by zsc on 2023/11/7.
//

#ifndef RUNE_BIN_IO_H
#define RUNE_BIN_IO_H
#include <cstdint>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

class Header {
public:
  char magic[5] = {'R', 'U', 'N', 'E', '\0'};
  uint8_t k_size{};
  uint32_t sample_count{};
  uint64_t record_count{};
};

class Record {
public:
  uint64_t kbin;
  uint64_t kpos;
  uint32_t sample_idx;
};

class bin_io {
private:
  std::string file_name;
  std::fstream fs;

public:
  std::unordered_map<uint64_t, std::pair<uint32_t, uint64_t>> mp_kmer_records;
  std::unordered_map<uint32_t, std::string> id_name;
  std::unordered_map<std::string, uint32_t> sample_id;
  uint8_t k_size{};
  explicit bin_io(std::string file_name);
  bin_io(std::string file_name, uint8_t k_size);
  bool
  write(std::unordered_map<uint64_t, std::pair<uint32_t, uint64_t>> &mp_kmer,
        std::unordered_map<uint32_t, std::string> &mp_sample_id);
  bool read();
};

#endif // RUNE_BIN_IO_H

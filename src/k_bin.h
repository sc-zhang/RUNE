//
// Created by zsc on 2023/11/2.
//

#ifndef RUNE_K_BIN_H
#define RUNE_K_BIN_H
#include <algorithm>
#include <cstdint>
#include <string>
#include <unordered_map>

class k_bin {
private:
  uint8_t k_size;
  uint64_t pos;
  uint64_t kbin;
  uint64_t rbin;
  const uint64_t KMER_MASK;
  std::string seq;
  std::unordered_map<char, uint64_t> mp_base = {
      {'A', 0LL}, {'C', 1LL}, {'G', 2LL}, {'T', 3LL}};
  std::unordered_map<uint64_t, char> mp_bin = {
      {0LL, 'A'}, {1LL, 'C'}, {2LL, 'G'}, {3LL, 'T'}};

public:
  explicit k_bin(std::string seq, uint8_t k_size = 17)
      : KMER_MASK(k_size < 32 ? ((1LL << (k_size * 2)) - 1)
                              : 0xFFFFFFFFFFFFFFFF) {
    this->seq = std::move(seq);
    this->k_size = k_size;
    this->pos = 0;
    this->kbin = 0;
    this->rbin = 0;
  }
  uint64_t kmer2bin(const std::string &k_mer);
  std::string bin2kmer(uint64_t k_bin);
  uint64_t rev_bin(uint64_t k_bin) const;
  void get_kmer();
  uint64_t get_kbin() const;
  uint64_t get_rbin() const;
  uint64_t get_pos() const;
  void reset_pos();
};

#endif // RUNE_K_BIN_H

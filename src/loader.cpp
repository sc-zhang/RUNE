//
// Created by zsc on 2023/11/11.
//
#include "loader.h"

void loader::load() {
  msg message = msg(false);
  message.info("Loading k-mer");
  // k_bin kb = k_bin("", k_size);
  bool is_valid = kmer_bin_io.read();
  if (is_valid) {
    message.info(std::to_string(kmer_bin_io.mp_kmer_records.size()) +
                 " Unique kmer Loaded");
  } else {
    message.err("Invalid or incomplete binary file, exiting...");
    exit(-1);
  }
  this->k_size = kmer_bin_io.k_size;
}
void loader::save(const std::string &output_file) {
  msg message = msg(false);
  message.info("Kmer size: " + std::to_string(this->k_size));
  message.info("Exporting");

  std::ofstream ofs;
  std::ostream *fs = &std::cout;
  if (!output_file.empty()) {
    ofs.open(output_file);
    fs = &ofs;
  }
  k_bin kb = k_bin("", this->k_size);
  for (auto &it : this->get_kmer_db()) {
    if (((it.second & rune::MASK::ID_MASK) >> 32) == rune::FLAG::UNKNOWN) {
      continue;
    }
    (*fs) << kb.bin2kmer(it.first) << "\t"
          << this->get_sample_name((it.second & rune::MASK::ID_MASK) >> 32)
          << "\t" << (it.second & rune::MASK::POS_MASK) << "\n";
  }
  ofs.close();
  message.info("Exported");
}
std::unordered_map<uint64_t, uint64_t> loader::get_kmer_db() const {
  return this->kmer_bin_io.mp_kmer_records;
}
uint32_t loader::get_sample_id(const std::string &sample) {
  return this->kmer_bin_io.sample_id[sample];
}
std::string loader::get_sample_name(const uint32_t &id) {
  return this->kmer_bin_io.id_name[id];
}

//
// Created by zsc on 2023/11/11.
//

#include "dumper.h"

void dumper::extract() {
  msg message = msg(true);
  message.info("Loading fasta");
  fasta_io ff = fasta_io(this->fasta_file);
  this->mp_seq = ff.read();
  if (this->mp_seq.empty()) {
    message.err("Invalid fasta file");
    exit(-1);
  }
  // set sample index, 1-based
  uint32_t idx = 1;
  for (auto &it : this->mp_seq) {
    this->id_sample[it.first] = idx;
    this->sample_id[idx] = it.first;
    ++idx;
  }

  // get uniq kmer map
  // k_bin=>sample_index
  // value 0 means not unique
  uint32_t seq_cnt = 0, skip_seq_cnt = 0;
  message.info("Generating k-mers with " + std::to_string(this->k_size));
  for (auto &it : this->mp_seq) {
    // skip sequence shorter than the size of kmer
    if (it.second.size() < this->k_size) {
      ++skip_seq_cnt;
      continue;
    }
    ++seq_cnt;
    k_bin kb = k_bin(it.second, this->k_size);
    while (kb.get_pos() <= it.second.size() - this->k_size) {
      kb.get_kmer();
      uint64_t kbin = kb.get_kbin();
      uint64_t rbin = kb.get_rbin();
      if (!this->mp_kmer.count(kbin)) {
        this->mp_kmer[kbin] =
            std::make_pair(this->id_sample[it.first], kb.get_pos());
      } else {
        this->mp_kmer[kbin] = std::make_pair(rune::flag::unknown, 0);
      }
      if (!this->mp_kmer.count(rbin)) {
        this->mp_kmer[rbin] =
            std::make_pair(this->id_sample[it.first], kb.get_pos());
      } else {
        this->mp_kmer[rbin] = std::make_pair(rune::flag::unknown, 0);
      }
    }
  }
  message.info(std::to_string(seq_cnt) + " sequences loaded, " +
               std::to_string(skip_seq_cnt) + " skipped.");
}

void dumper::save() {
  msg message = msg(true);
  message.info("Writing k-mers");
  bin_io bio = bin_io(this->kmer_file, this->k_size);
  bio.write(this->mp_kmer, this->sample_id);
}

std::unordered_map<uint64_t, std::pair<uint32_t, uint64_t>>
dumper::get_kmer_db() {
  return this->mp_kmer;
}
uint32_t dumper::get_sample_id(const std::string &sample) {
  return this->id_sample[sample];
};
std::string dumper::get_sample_name(const uint32_t &id) {
  return this->sample_id[id];
}
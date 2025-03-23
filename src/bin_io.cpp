//
// Created by zsc on 2023/11/7.
//

#include "bin_io.h"
#include "rune.h"
#include <cstring>
#include <utility>

bin_io::bin_io(std::string file_name) {
  this->file_name = std::move(file_name);
}

bin_io::bin_io(std::string file_name, uint8_t k_size) {
  this->file_name = std::move(file_name);
  this->k_size = k_size;
}

bool bin_io::write(std::unordered_map<uint64_t, uint64_t> &mp_kmer,
                   std::unordered_map<uint32_t, std::string> &mp_sample_id) {
  fs.open(this->file_name, std::ios_base::out | std::ios_base::binary);
  if (fs) {
    Header header = Header();
    if (!fs.write((char *)&header, sizeof(header))) {
      fs.close();
      return false;
    }

    header.k_size = this->k_size;
    header.sample_count = mp_sample_id.size();
    header.record_count = 0;

    std::vector<std::string> samples(header.sample_count);

    for (auto &it : mp_sample_id) {
      samples[it.first - 1] = it.second;
    }
    for (auto &sample : samples) {
      uint16_t id_length = sample.size();
      if (!fs.write((char *)&id_length, sizeof(id_length))) {
        fs.close();
        return false;
      }
      if (!fs.write(sample.c_str(), id_length)) {
        fs.close();
        return false;
      }
    }
    uint64_t keep_record_count = 0;
    for (auto &it : mp_kmer) {
      Record record = Record();
      if (((it.second & rune::MASK::ID_MASK) >> 32) == rune::FLAG::UNKNOWN) {
        continue;
      }
      record.kbin = it.first;
      record.sample_idx_kpos = it.second;
      if (fs.write((char *)&record, sizeof(record))) {
        ++keep_record_count;
      } else {
        fs.close();
        return false;
      }
    }

    header.record_count = keep_record_count;
    fs.seekp(0, std::ios::beg);
    if (!fs.write((char *)&header, sizeof(header))) {
      fs.close();
      return false;
    }
    fs.close();
    return true;
  } else {
    return false;
  }
}

bool bin_io::read() {
  fs.open(this->file_name, std::ios_base::in | std::ios_base::binary);
  if (fs) {
    Header header;
    fs.read((char *)&header, sizeof(header));
    if (strcmp(header.magic, "RUNE") != 0) {
      fs.close();
      return false;
    }
    this->k_size = header.k_size;
    for (uint32_t i = 0; i < header.sample_count; ++i) {
      uint16_t id_length;
      fs.read((char *)&id_length, sizeof(id_length));
      std::string sample;
      char *buffer = new char[id_length + 1];
      fs.read(buffer, id_length);
      buffer[id_length] = '\0';
      sample = buffer;
      this->id_name[i + 1] = sample;
      this->sample_id[sample] = i + 1;
    }
    for (uint64_t i = 0; i < header.record_count; ++i) {
      Record record{};
      fs.read((char *)&record, sizeof(record));
      this->mp_kmer_records[record.kbin] = record.sample_idx_kpos;
    }
    fs.close();
    if (this->mp_kmer_records.size() != header.record_count) {
      return false;
    }
    return true;
  } else {
    return false;
  }
}

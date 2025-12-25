//
// Created by zsc on 24-11-14.
//

#include "args.h"

args::args(int argc, char **argv) {
  bool is_valid = true;
  // if no argument, show help
  if (argc == 1) {
    is_valid = false;
  }
  // load arguments
  for (int i = 1; i < argc; ++i) {
    if (i == 1) {
      if (strcmp(argv[i], "dump") != 0 && strcmp(argv[i], "load") != 0) {
        if (strcmp(argv[i], "-h") != 0 && strcmp(argv[i], "--help") != 0) {
          std::cout << "Invalid operation: " + std::string(argv[i])
                    << std::endl;
        }
        is_valid = false;
        break;
      } else {
        this->operate = argv[1];
      }
    } else {
      if (strcmp(argv[i], "-i") == 0) {
        if (i + 1 < argc) {
          this->input_file = argv[i + 1];
          ++i;
          continue;
        } else {
          std::cout << "Input file not found" << std::endl;
          is_valid = false;
          break;
        }
      } else if (strcmp(argv[i], "-o") == 0) {
        if (i + 1 < argc) {
          this->output_file = argv[i + 1];
          ++i;
          continue;
        } else {
          std::cout << "Output file not found" << std::endl;
          is_valid = false;
          break;
        }
      } else if (strcmp(argv[i], "-k") == 0) {
        if (i + 1 < argc) {
          std::string k_size_arg = argv[i + 1];
          this->k_size = 0;
          for (char &c : k_size_arg) {
            this->k_size *= 10;
            this->k_size += c - '0';
          }
          if (this->k_size < 3 || this->k_size > 32) {
            std::cout << "Invalid kmer size, should be in range of [3, 32]"
                      << std::endl;
            is_valid = false;
            break;
          }
        } else {
          std::cout << "Kmer size not set" << std::endl;
          is_valid = false;
          break;
        }
      }
    }
  }
  if (!is_valid || this->input_file.empty()) {
    std::string help_string =
        "usage: rune [-h] {dump,load} -i input_file -k "
        "kmer_size -o output_file\n"
        "options:\n"
        "  -h, --help        show this help message and exit\n"
        "  -i input_file     Input file\n"
        "  -k kmer_size      kmer size, only for dump, "
        "should be [3, 32], default=25\n"
        "  -o output_file    Output file\n";
    std::cout << help_string << std::endl;
    this->args_enough = false;
  }
}

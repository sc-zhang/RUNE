#include "args.h"
#include "dumper.h"
#include "loader.h"
#include "msg.h"
#include <bitset>

int main(int argc, char *argv[]) {
  args arguments = args(argc, argv);
  if (arguments.args_enough) {
    msg message = msg(false);
    message.info("Starting");
    if (arguments.operate == "load") {
      loader kmer_loader = loader(arguments.input_file, arguments.k_size);
      kmer_loader.load();
      kmer_loader.save(arguments.output_file);
    } else if (arguments.operate == "dump") {
      if (arguments.output_file.empty()) {
        message.err(
            "dump could not work if no output file provided, aborting...");
        exit(-1);
      }
      dumper kmer_dumper =
          dumper(arguments.input_file, arguments.k_size, arguments.output_file);
      kmer_dumper.extract();
      kmer_dumper.save();
    }
    message.info("Finished");
  }
  return 0;
}

#include <zlib.h>
#include <argparse/argparse.hpp>
#include "../lib/include/HyperLogLog.hpp"

extern "C" {
#include "kseq.h"
}

KSEQ_INIT(gzFile, gzread)

argparse::ArgumentParser get_parser();

int main(int argc, char* argv[])
{
    gzFile fp;
    kseq_t* seq;
    auto parser = get_parser();
    parser.parse_args(argc, argv);
    auto k = parser.get<std::size_t>("-k");
    auto b = parser.get<std::size_t>("-b");
    auto input_filename = parser.get<std::string>("--input");
    if (k > 32) throw std::runtime_error("0 < k <= 32");
    if (b >= 64) throw std::runtime_error("b < 64");
    HyperLogLog hll(k, b);

    fp = NULL;
    if (input_filename == "") {
        if ((fp = gzdopen(fileno(stdin), "r")) == NULL) {
            return 1;
        }
    } else {
        if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
            return 1;
        }
    }
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        hll.add(seq->seq.s, seq->seq.l);
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    std::cout << hll.count() << "\n";
    return 0;
}

argparse::ArgumentParser get_parser()
{
    argparse::ArgumentParser parser("khll");
    parser.add_description("Estimate number of distinct k-mers from a stream");
    parser.add_argument("-k")
        .help("k-mer size")
        .scan<'u', std::size_t>()
        .required();
    parser.add_argument("-b")
        .help("header size (number of msb bits used as index)")
        .scan<'u', std::size_t>()
        .default_value(std::size_t(12));
    parser.add_argument("-i", "--input")
        .help("input filename [stdin]")
        .default_value(std::string(""));
    return parser;
}

#include "../include/build.hpp"
#include "../../lib/include/HyperLogLog.hpp"
#include <filesystem>
#include <iostream>
#include <zlib.h>
extern "C" {
#include "../include/kseq.h"
}

KSEQ_INIT(gzFile, gzread)

int build_main(const argparse::ArgumentParser& parser)
{
    using namespace sketching;
    auto k = parser.get<std::size_t>("-k");
    auto b = parser.get<std::size_t>("-b");
    auto e = parser.get<double>("-e");
    auto passthrough = parser.get<bool>("--passthrough");
    auto input_filename = parser.get<std::string>("--input");
    auto sketch_filename = parser.get<std::string>("--sketch");

    gzFile fp = NULL;
    if (input_filename == "") {
        if ((fp = gzdopen(fileno(stdin), "r")) == NULL) {
            return 1;
        }
    } else {
        if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
            return 1;
        }
    }

    HyperLogLog hll;
    if (sketch_filename != "" and std::filesystem::exists(sketch_filename)) {
        hll = HyperLogLog::load(sketch_filename);
    } else if (e < 0) {
        hll = HyperLogLog(k, uint8_t(b));
    } else {
        hll = HyperLogLog(k, e);
    }

    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        hll.add(seq->seq.s, seq->seq.l);
        if (passthrough) {
            std::cout <<  ">" << std::string(seq->name.s, seq->name.l) << "\n";
            std::cout << std::string(seq->seq.s, seq->seq.l) << "\n";
            if (seq->qual.l != 0) {
                if (seq->qual.l != seq->seq.l) {
                    std::cerr << "sequence and its quality string do not match in length\n";
                }
                std::cout << "#\n";
                std::cout << std::string(seq->qual.s, seq->qual.l) << "\n";
            }
        }
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);

    if (sketch_filename != "") {
        hll.store(sketch_filename);
    }

    std::cerr << hll.count() << "\n";
    return 0;
}

argparse::ArgumentParser get_parser_build()
{
    argparse::ArgumentParser parser("build");
    parser.add_description("Build HyperLogLog from FastX files");
    parser.add_argument("-k")
        .help("k-mer size")
        .scan<'u', std::size_t>()
        .required();
    parser.add_argument("-b")
        .help("header size (number of msb bits used as index)")
        .scan<'u', std::size_t>()
        .default_value(std::size_t(12));
    parser.add_argument("-e")
        .help("error rate of the HLL sketch. Supersedes option -b if present")
        .scan<'f', double>()
        .default_value(-1.0);
    parser.add_argument("-p", "--passthrough")
        .help("forward records to stdout")
        .default_value(false)
        .implicit_value(true);
    parser.add_argument("-i", "--input")
        .help("input filename [stdin]")
        .default_value(std::string(""));
    parser.add_argument("-s", "--sketch")
        .help("hll sketch, create or update with stream depending on if it exists or not")
        .default_value("");
    return parser;
}

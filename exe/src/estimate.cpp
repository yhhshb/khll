#include "../include/estimate.hpp"
#include "../../lib/include/HyperLogLog.hpp"
#include <filesystem>
#include <iostream>

int estimate_main(const argparse::ArgumentParser& parser)
{
    using namespace sketching;
    auto sketch_filename = parser.get<std::string>("--sketch");
    HyperLogLog hll;
    if (std::filesystem::exists(sketch_filename)) hll = HyperLogLog::load(sketch_filename);
    else throw std::runtime_error("sketch does not exist");
    std::cout << hll.count() << "\n";
    return 0;
}

argparse::ArgumentParser get_parser_estimate()
{
    argparse::ArgumentParser parser("estimate");
    parser.add_description("Print sketch estimation");
    parser.add_argument("-s", "--sketch")
        .help("hll sketch to query")
        .required();
    return parser;
}
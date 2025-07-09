#include "../include/build.hpp"
#include "../../lib/include/HyperLogLog.hpp"
#include <filesystem>
#include <fstream>

int merge_main(const argparse::ArgumentParser& parser) 
{
    auto trim = [](std::string& s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
    };
    using namespace sketching;
    auto file_lists = parser.get<std::vector<std::string>>("--input-lists");
    auto sketches_filenames = parser.get<std::vector<std::string>>("sketches");
    auto output_filename = parser.get<std::string>("--output-sketch");

    for (auto const& list_filename : file_lists) {
        std::string buffer;
        std::ifstream flist(list_filename);
        while (std::getline(flist, buffer)) {
            trim(buffer);
            sketches_filenames.push_back(buffer);
        }
    }

    HyperLogLog hll;
    if (not sketches_filenames.empty()) {
        hll = HyperLogLog::load(sketches_filenames.back());
        sketches_filenames.pop_back();
        for (const auto& sketch_filename : sketches_filenames) {
            HyperLogLog other = HyperLogLog::load(sketch_filename);
            hll += other;
        }
    }

    if (output_filename != "") hll.store(output_filename);
    std::cerr << hll.count() << " / " << hll.size() << "\n";

    return 0;
}

argparse::ArgumentParser get_parser_merge()
{
    argparse::ArgumentParser parser("merge");
    parser.add_description("Merge HyperLogLogs");
    parser.add_argument("-i", "--input-lists")
        .help("file(s) listing sketches to be merged (1 sketch filename per row)")
        .nargs(argparse::nargs_pattern::any);
    parser.add_argument("sketches")
        .help("list of sketches to be merged")
        .nargs(argparse::nargs_pattern::any);
    parser.add_argument("-o", "--output-sketch")
        .help("output sketch (optional)")
        .default_value("");
    return parser;
}

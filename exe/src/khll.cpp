#include <iostream>
#include "../include/build.hpp"
#include "../include/estimate.hpp"
#include "../include/merge.hpp"

int main(int argc, char* argv[])
{
    auto build_parser = get_parser_build();
    auto estimate_parser = get_parser_estimate();
    auto merge_parser = get_parser_merge();
    argparse::ArgumentParser program(argv[0]);
    program.add_subparser(build_parser);
    program.add_subparser(estimate_parser);
    program.add_subparser(merge_parser);
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    if (program.is_subcommand_used(build_parser)) return build_main(build_parser);
    else if (program.is_subcommand_used(estimate_parser)) return estimate_main(estimate_parser);
    else if (program.is_subcommand_used(merge_parser)) return merge_main(merge_parser);
    else std::cerr << program << std::endl;
    return 0;
}



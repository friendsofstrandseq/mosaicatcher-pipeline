/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#include <iostream>

#include <htslib/sam.h>
#include "version.hpp"
#include "count.hpp"
#include "segmentation.hpp"
#include "simulate.hpp"
#include "calc_bins.hpp"
#include "sces.hpp"
#include "ploidy.hpp"

int main(int argc, char **argv)
{

    if (argc >= 2 && std::string(argv[1]) == "count") {
        return main_count(argc-1, argv+1);

    } else if (argc >= 2 && std::string(argv[1]) == "simulate") {
        return main_simulate(argc-1, argv+1);

    } else if (argc >= 2 && std::string(argv[1]) == "segment") {
        return main_segment(argc-1, argv+1);

    } else if (argc >= 2 && std::string(argv[1]) == "makebins") {
        return main_calc_bins(argc-1, argv+1);

    } else if (argc >= 2 && std::string(argv[1]) == "states") {
        return main_strand_states(argc-1, argv+1);

    } else if (argc >= 2 && std::string(argv[1]) == "hmm") {
        return main_hmm(argc-1, argv+1);

    } else if (argc >= 2 && std::string(argv[1]) == "--version") {
        std::cout << "Mosaicatcher " << STRINGIFYMACRO(MOSAIC_VERSION_MAJOR) ;
        std::cout << "." << STRINGIFYMACRO(MOSAIC_VERSION_MINOR) << " (";
		std::cout << STRINGIFYMACRO(MOSAIC_GIT_COMMIT_HASH) << "-";
        std::cout << STRINGIFYMACRO(MOSAIC_GIT_BRANCH) << ")" << std::endl;
        std::cout << "htslib " << STRINGIFYMACRO(HTSLIB_VERSION) << std::endl;
        std::cout << "boost version " << STRINGIFYMACRO(BOOST_VERSION) << std::endl;
        return 0;

    } else {
        std::cout << std::endl;
        std::cout << "Mosaicatcher " << STRINGIFYMACRO(MOSAIC_VERSION_MAJOR);
        std::cout << "." << STRINGIFYMACRO(MOSAIC_VERSION_MINOR) << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:   " << argv[0] << " [--version] [--help] <command> [options]" << std::endl << std::endl;
        std::cout << "Commands:" << std::endl;
        std::cout << "    count       Count binned reads in Strand-seq cells" << std::endl;
        //std::cout << "    hmm         Count binned reads in arbitrary ploidy" << std::endl;
        std::cout << "    segment     Find a segmentation across binned counts" << std::endl;
        std::cout << "    simulate    Simulate Strand-seq data" << std::endl;
        std::cout << "    makebins    Create variable-width bins based on real data" << std::endl;
        std::cout << "    states      Determine strand states including SCEs" << std::endl;
        std::cout << std::endl;

        if (argc >= 2 && std::string(argv[1]) == "--help") {
            std::cout << "Mosaicatcher calls structural variants (SVs) in Strand-seq data." << std::endl;
            std::cout << "It is currently under development." << std::endl;
            std::cout << "This software comes with absolutely no warranty." << std::endl << std::endl;
            return 0;
        }
    }
    return 1;
}


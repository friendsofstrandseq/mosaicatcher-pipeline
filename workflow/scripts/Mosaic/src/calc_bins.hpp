/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>


struct Conf_calc_bins {
    boost::filesystem::path f_in;
    boost::filesystem::path f_out;
    boost::filesystem::path f_excl;
    unsigned minMapQual;
    unsigned num_reads;
    unsigned window;
};

int main_calc_bins(int argc, char **argv)
{
    Conf_calc_bins conf;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
    ("help,?", "show help message")
    ("mapq,q", boost::program_options::value<unsigned>(&conf.minMapQual)->default_value(10), "min mapping quality")
    ("window,w", boost::program_options::value<unsigned>(&conf.window)->default_value(100000), "window size to approximate")
    ("numreads,n", boost::program_options::value<unsigned>(&conf.num_reads)->default_value(20), "sample 1/n of reads (reduce memory)")
    ("out,o", boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("bins.bed"), "output file for bins")
    ("exclude,x", boost::program_options::value<boost::filesystem::path>(&conf.f_excl), "Exclude chromosomes (no regions!)")
    ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&conf.f_in), "input bam file (WGS or merged libraries)")
    ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;

    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if (vm.count("help") || (!vm.count("input-file"))) {
        std::cout << "Usage: " << argv[0] << " [OPTIONS] <wgs.bam>" << std::endl;
        std::cout << visible_options << "\n";
        if(vm.count("help")) {
            std::cout << "Specify whole genome sequencing data (or a set of Strand-seq cells" << std::endl;
            std::cout << "merged into a single file) which were sequenced under equal conditions." << std::endl;
            std::cout << "This tool will create bins of variable width but which contian the same" << std::endl;
            std::cout << "number of reads. This way we hope to overcome mappability issues." << std::endl;
            std::cout <<  std::endl;
        }
        return vm.count("help") ? 0 : 1;
    }

    // Open bam file
    samFile* samfile = sam_open(conf.f_in.string().c_str(), "r");
    if (samfile == NULL) {
        std::cerr << "Fail to open file " << conf.f_in.string() << std::endl;
        return 1;
    }
    hts_idx_t* idx = sam_index_load(samfile, conf.f_in.string().c_str());
    if (idx == NULL) {
        std::cerr << "Fail to open index for " << conf.f_in.string() << std::endl;
        return 1;
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Exclude chromosomes:
    std::set<std::string> excl_chroms;
    if (vm.count("exclude")) {
        std::string line;
        std::ifstream excl_file(conf.f_excl.string().c_str());
        if (excl_file.is_open()) {
            while ( getline (excl_file, line) ) {
                boost::trim(line);
                if (line.size()>0)
                    excl_chroms.insert(line);
            }
            excl_file.close();
        } else {
            std::cerr << "Fail to open exclude file " << conf.f_excl.string() << std::endl;
            return 1;
        }
    }

    unsigned n_chroms = 0;
    for (int chrom = 0; chrom < hdr->n_targets; ++chrom)
        if (!excl_chroms.count(hdr->target_name[chrom]))
            ++n_chroms;
    boost::progress_display show_progress( n_chroms );


    // Save all positions in vector (sampled to every n-th position)
    std::vector<std::vector<uint32_t> > positions;
    std::vector<int> tids;
    for (int chrom = 0; chrom < hdr->n_targets; ++chrom) {

        // Skip unused chromosomes
        if (excl_chroms.count(hdr->target_name[chrom])) continue;
        positions.push_back(std::vector<uint32_t>());
        tids.push_back(chrom);

        ++show_progress;
        unsigned n_reads = 0;
        hts_itr_t* iter = sam_itr_queryi(idx, chrom, 0, hdr->target_len[chrom]);
        bam1_t* rec = bam_init1();
        while (sam_itr_next(samfile, iter, rec) >= 0) {

            if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP))
                continue;
            if ((rec->core.qual < conf.minMapQual) || (rec->core.tid<0))
                continue;
            if (++n_reads == conf.num_reads) {
                positions[chrom].push_back(rec->core.pos);
                n_reads = 0;
            }
        }
        bam_destroy1(rec);
        hts_itr_destroy(iter);
    }
    hts_idx_destroy(idx);
    sam_close(samfile);


    // 'positions' only contains non-exlcuded chromosomes
    // 'tids' contains matching tids.

    // Calculate average number of reads per window-size
    uint64_t total_num = 0;
    uint64_t total_len = 0;
    for (unsigned i = 0; i<positions.size(); ++i) {
        total_num    += positions[i].size() * conf.num_reads;
        total_len += hdr->target_len[ tids[i] ];
    }
    float avg_num_reads = total_num/(float)total_len * conf.window;
    unsigned jump = std::round(avg_num_reads/(float)conf.num_reads);
    if (jump < 1) jump = 1;


    std::cout << "Approximated window size wished: " << conf.window/1000 << "kb." << std::endl;
    std::cout << "Input data contains: " << avg_num_reads << " reads per " << conf.window/1000 << "kb." << std::endl;
    std::cout << "Genome will be split every " << jump * conf.num_reads << " reads (divergence due to sampling)." << std::endl;
    std::cout << "Writing to " << conf.f_out.string() << std::endl;

    // Prepare output file
    std::ofstream out(conf.f_out.string());
    if (!out.is_open()) {
        std::cerr << "Fail to open file " << conf.f_out.string() << std::endl;
        return 1;
    }

    // Split positions every ~avg_num_reads and output
    for (unsigned i=0; i<positions.size(); ++i) {
        unsigned prev_pos = 0;
        for (unsigned list_pos = jump; list_pos < positions[i].size(); list_pos+=jump)
        {
            out << hdr->target_name[ tids[i] ] << "\t" << prev_pos << "\t" << positions[i][list_pos] << std::endl;
            prev_pos = positions[i][list_pos];
        }
        // last bin:
        out << hdr->target_name[ tids[i] ] << "\t" << prev_pos << "\t" << hdr->target_len[ tids[i] ] << std::endl;
    }

    out.close();
    return 0;
}


/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <tuple>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <htslib/sam.h>

#include "version.hpp"
#include "intervals.hpp"
#include "counter.hpp"
#include "distribution.hpp"
#include "hmm.hpp"
#include "iocounts.hpp"


/**
 * @file
 * @defgroup count Bin, count and classify W/C reads.
 *
 * Summary of how Strand-seq data is binned, counted and classified.
 *
 * ## Strand-seq read counting
 *
 * @todo write documentation about counting.
*/


using interval::Interval;
using count::TGenomeCounts;
using count::Counter;


struct Conf {
    std::vector<boost::filesystem::path> f_in;
    boost::filesystem::path f_out;
    boost::filesystem::path f_bins;
    boost::filesystem::path f_excl;
    boost::filesystem::path f_info;
    boost::filesystem::path f_sample_info;
    boost::filesystem::path f_removed_bins;
    boost::filesystem::path f_segments;
    int minMapQual;
    unsigned int window;
    std::string mode;
};


/**
 *
 */
void run_standard_HMM(std::vector<TGenomeCounts> & counts,
                      std::vector<unsigned> const & good_cells,
                      std::vector<CellInfo>  & cells,
                      std::vector<unsigned> const & good_bins,
                      std::vector<int32_t> const & good_map,
                      std::unordered_map<std::string, SampleInfo> const & samples,
                      float p_trans)
{
    // Set up and run HMM:
    hmm::HMM<unsigned, hmm::MultiVariate<hmm::NegativeBinomial> > hmm({"CC", "WC", "WW"});
    hmm.set_initials({0.3333, 0.3333, 0.3333});
    hmm.set_transitions({1-2*p_trans, p_trans,     p_trans,    \
        p_trans,     1-2*p_trans, p_trans,    \
        p_trans,     p_trans,     1-2*p_trans});

    for (auto i = good_cells.begin(); i != good_cells.end(); ++i)
    {
        // set NB(n,p) parameters according to `p` of sample and mean of cell.
        float p = samples.at(cells[*i].sample_name).p;
        float n = (float)cells[*i].mean_bin_count * p / (1-p);
        float a = 0.1;
        cells[*i].nb_p = p;
        cells[*i].nb_r = n;
        cells[*i].nb_a = a;

        //std::cout << "NB parameters for cell <?>" << ": p=" << p << "\tn=" << n << "\tz=" << z << std::endl;

        hmm.set_emissions( {\
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p, (1-a)*n), hmm::NegativeBinomial(p, a*n)}), // CC
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p, n/2),     hmm::NegativeBinomial(p, n/2)}), // WC
            hmm::MultiVariate<hmm::NegativeBinomial>({hmm::NegativeBinomial(p, a*n),     hmm::NegativeBinomial(p, (1-a)*n)})  // WW
        });
        run_HMM(hmm, counts[*i], good_bins, good_map);
    }
}





int main_count(int argc, char **argv)
{

    // Command line options
    Conf conf;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
    ("help,?", "show help message")
    ("verbose,v", "Be more verbose in the output")
    ("mapq,q", boost::program_options::value<int>(&conf.minMapQual)->default_value(10), "min mapping quality")
    ("window,w", boost::program_options::value<unsigned int>(&conf.window)->default_value(500000), "window size of fixed windows")
    ("out,o", boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("out.txt.gz"), "output file for counts + strand state (gz)")
    ("bins,b", boost::program_options::value<boost::filesystem::path>(&conf.f_bins), "BED file with manual bins (disables -w). See also 'makebins'")
    ("exclude,x", boost::program_options::value<boost::filesystem::path>(&conf.f_excl), "Exclude chromosomes and regions")
    ("info,i", boost::program_options::value<boost::filesystem::path>(&conf.f_info), "Write info about samples")
    ("do-not-filter-by-WC", "When black-listing bins, only consider coverage and not WC/WW/CC states")
    ("do-not-blacklist-hmm", "Do not output a blacklist (None bins). Bins will be blacklisted for parameter estimation, but not during HMM")
    ;
    // Note: Currently the blacklisting is done after counting via R/norm.R. A better way would be to 
    // input the blacklist + normalization into MosaiCatcher during the counting, then the HMM could be run
    // only on non-blacklisted bins. Now it is run through many bad bin that potentially affect the quality of results.

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
    ("input-file", boost::program_options::value<std::vector<boost::filesystem::path> >(&conf.f_in), "input bam file(s)")
    ("sample_info,S", boost::program_options::value<boost::filesystem::path>(&conf.f_sample_info),   "write info per sample")
    ("removed_bins,R", boost::program_options::value<boost::filesystem::path>(&conf.f_removed_bins), "bins that were removed (bed file)")
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
    if (!vm["window"].defaulted() && vm.count("bins")) {
        std::cerr << "[Error] -w and -b cannot be specified together" << std::endl << std::endl;
        goto print_usage_and_exit;
    }
    if (vm.count("bins") && vm.count("exclude")) {
        std::cerr << "[Error] Exclude chromosomes (-x) have no effect when -b is specified. Stop" << std::endl << std::endl;
        goto print_usage_and_exit;
    }

    if (vm.count("help") || !vm.count("input-file"))
    {
    print_usage_and_exit:
        std::cout << std::endl;
        std::cout << "Mosaicatcher " << STRINGIFYMACRO(MOSAIC_VERSION_MAJOR);
        std::cout << "." << STRINGIFYMACRO(MOSAIC_VERSION_MINOR) << std::endl;
        std::cout << "> Count reads from Strand-seq BAM files." << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:   " << argv[0] << " [OPTIONS] <cell1.bam> <cell2.bam> ..." << std::endl << std::endl;
        std::cout << visible_options << std::endl;
        std::cout << "Notes:" << std::endl;
        std::cout << "  * writes a table of bin counts and state classifcation as a gzip file (default: out.txt.gz)" << std::endl;
        std::cout << "  * Reads are counted by start position" << std::endl;
        std::cout << "  * One cell per BAM file, including SM tag in header" << std::endl;
        std::cout << "  * For paired-end data, only read 1 is counted" << std::endl;
        return vm.count("help") ? 0 : 1;
    }


    /////////////////////////////////////////////////////////// global variables
    /* leave one BAM header open to get chrom names & lengths */
    bam_hdr_t* hdr = NULL;

    /* regarding each cell */
    std::vector<CellInfo>       cells(conf.f_in.size());
    std::vector<TGenomeCounts>  counts(conf.f_in.size());
    std::vector<unsigned>       good_cells;

    /* regarding each sample */
    std::unordered_map<std::string, SampleInfo> samples;

    /* regarding bins */
    std::vector<Interval>       bins;
    std::vector<int32_t>        chrom_map;
    std::vector<unsigned>       good_bins;
    std::vector<int32_t>        good_map;
    ////////////////////////////////////////////////////////////////////////////



    //
    // Chapter: Binning & counting
    // ===========================
    //

    // Read sample names from headers.
    // Keep one header throughout the program.
    if (vm.count("verbose")) std::cout << "[Info] Exploring SAM headers..." << std::endl;
    for(unsigned i = 0; i < conf.f_in.size(); ++i)
    {
        cells[i].id = (int32_t)i;
        cells[i].bam_file = conf.f_in[i].string();
        samFile* samfile = sam_open(conf.f_in[i].string().c_str(), "r");
        if (samfile == NULL) {
            std::cerr << "[Error] Fail to open file " << conf.f_in[i].string() << std::endl;
            return 1;
        }
        hdr = sam_hdr_read(samfile);
        if (!get_RG_tag("SM", hdr->text, cells[i].sample_name)) {
            std::cerr << "[Error] Each BAM file has to have exactly one RG tag. Group cells " << std::endl;
            std::cerr << "        belonging to the same sample by the SM tag." << std::endl;
            std::cerr << "        Problematic file: " << conf.f_in[i].string() << std::endl << std::endl;
            goto print_usage_and_exit;
        }
        if (!get_RG_tag("ID", hdr->text, cells[i].cell_name, /*allow_multiple_matches = */ true)) {
            std::cerr << "[Error] Each BAM file has to have exactly one RG tag." << std::endl;
            std::cerr << "        Problematic file: " << conf.f_in[i].string() << std::endl;
            goto print_usage_and_exit;
        }
        sam_close(samfile);
    }


    // Bin the genome
    unsigned median_binsize;
    chrom_map = std::vector<int32_t>(hdr->n_targets, -1);
    if (vm.count("bins"))
    {
        if (!read_dynamic_bins(bins,
                               chrom_map,
                               conf.f_bins.string().c_str(),
                               hdr))
            return 1;
        TMedianAccumulator<unsigned> med_acc;
        for (Interval const & b : bins)
            med_acc(b.end - b.start);
        median_binsize = boost::accumulators::median(med_acc);
        std::cout << "[Info] Reading " << bins.size() << " variable-width bins with median bin size of " << round(median_binsize/1000) << "kb" << std::endl;
    }
    else
    {
        std::vector<Interval> exclude;
        if (vm.count("exclude")) {
            read_exclude_file(conf.f_excl.string(), hdr, exclude, vm.count("verbose"));
            sort(exclude.begin(), exclude.end(), interval::invt_less);
        }
        std::cout << "[Info] Creating " << round(conf.window/1000) << "kb bins with " << exclude.size() << " excluded regions" << std::endl;
        create_fixed_bins(bins,
                          chrom_map,
                          conf.window,
                          exclude,
                          hdr->n_targets,
                          hdr->target_len);
        median_binsize = conf.window;
    }
    // add last element for easy calculation of number of bins
    chrom_map.push_back((int32_t)bins.size());


    // Count in bins. If A bam file cannot be read, the cell is ignored and
    //     the respective entry in `counts` and `cells` will be erased.
    std::cout << "[Info] Reading " << conf.f_in.size() <<  " BAM files...";
    boost::progress_display show_progress1(conf.f_in.size());
    for (unsigned i = 0, i_f = 0; i_f < conf.f_in.size(); ++i, ++i_f)
    {
        if (!count_sorted_reads(conf.f_in[i_f].string(),
                                bins,
                                chrom_map,
                                hdr,
                                conf.minMapQual,
                                counts[i],
                                cells[i]))
        {
            std::cerr << "[Warning] Ignoring cell " << conf.f_in[i_f].string() << std::endl;
            counts.erase(counts.begin()+i);
            cells.erase(cells.begin()+i);
            --i;
        }
        ++show_progress1;
    }





    //
    // Chapter: Filter cells and bins and estimate NB parameter p
    // ==========================================================
    //

    // median per cell
    count::set_median_per_cell(counts, cells);

    // filter cells with low counts
    good_cells = count::get_good_cells(counts, cells);

    // filter bins with abnormal counts
    if (good_cells.size() < 5) {
        std::cerr << "[Warning] Only few cells with sufficient coverage. I will not filter bad bins" << std::endl;
        good_bins.resize(bins.size());
        std::iota(good_bins.begin(), good_bins.end(), 0); // fill with 0,1,2,...
    } else {
        good_bins = count::get_good_bins(counts,
                                         cells,
                                         good_cells,
                                         vm.count("verbose"),
                                         !vm.count("do-not-filter-by-WC"));
        if (vm.count("verbose")) std::cout << "[Info] Filtered out " << bins.size() - good_bins.size() << " bad bins from " << bins.size() << std::endl;
    }

    // build chrom_map for good bins
    good_map = std::vector<int32_t>(chrom_map.size() - 1, -1);
    int32_t pos = 0;
    for (int32_t chr = 0; chr < static_cast<int32_t>(good_map.size()); ++chr) {
        while (pos < good_bins.size() && bins[good_bins[pos]].chr < chr)
            ++pos;
        // now goodit is either at first occurence of chr, or at the end.
        if (pos >= good_bins.size()) good_map[chr] = (int32_t)good_bins.size();
        else good_map[chr] = pos;
    }
    // add last element for easy calculation of number of bins
    good_map.push_back((int32_t)good_bins.size());



    // calculate cell means and cell variances, grouped by sample (not cell)
    calculate_new_cell_mean(samples, cells, counts, good_cells, good_bins);


    // Estimation of parameter p per sample (should work even with one cell only)
    for (auto it = samples.begin(); it != samples.end(); ++it) {
        SampleInfo & s = it->second;
        s.p = std::inner_product(s.means.begin(), s.means.end(), s.means.begin(), 0.0f) \
        / std::inner_product(s.means.begin(), s.means.end(), s.vars.begin(), 0.0f);
    }

    // Write sample information to file
    if (vm.count("sample_info")) {
        if (vm.count("verbose")) std::cout << "[Write] sample information: " << conf.f_sample_info.string() << std::endl;
        std::ofstream out(conf.f_sample_info.string());
        if (out.is_open()) {
            out << "sample\tcells\tp\tmeans\tvars" << std::endl;
            for (auto it = samples.begin(); it != samples.end(); ++it) {
                SampleInfo const & s = it->second;
                out << it->first << "\t" << s.means.size() << "\t" << s.p << "\t" << s.means[0];
                for (unsigned k=1; k<s.means.size(); ++k) out << "," << s.means[k];
                out << "\t" << s.vars[0];
                for (unsigned k=1; k<s.vars.size(); ++k) out << "," << s.vars[k];
                out << std::endl;
            }
        } else {
            std::cerr << "[Warning] Cannot write to " << conf.f_sample_info.string() << std::endl;
        }
    }




    //
    // Chapter: Run HMM
    // ================
    //
    if(vm.count("do-not-blacklist-hmm")) {
        if (vm.count("verbose")) std::cout << "[Info] Previous filters are not used during HMM phase" << std::endl;
        std::vector<unsigned> normal_bins(bins.size());
        std::iota(normal_bins.begin(), normal_bins.end(), 0);
        run_standard_HMM(counts,
                         good_cells,
                         cells,
                         normal_bins,
                         chrom_map,
                         samples,
                         10.0f / bins.size());
    } else {
        run_standard_HMM(counts,
                         good_cells,
                         cells,
                         good_bins,
                         good_map,
                         samples,
                         10.0f / good_bins.size());
    }



    // Print cell information:
    if (vm.count("info")) {
        if (vm.count("verbose")) std::cout << "[Write] Cell summary: " << conf.f_info.string() << std::endl;
        write_cell_info(conf.f_info.string(), cells);
    }


    // Write final counts + classification
    std::cout << "[Write] count table: " << conf.f_out.string() << std::endl;
    {
        // TODO: why do I pass vector<pair>? I could make it two separate vectors. Just check where else the io function is called.
        struct sample_cell_name_wrapper {
            std::vector<CellInfo> const & cells;
            sample_cell_name_wrapper(std::vector<CellInfo> const & cells) : cells(cells)
            {}
            std::pair<std::string,std::string> operator[](size_t i) const {
                return std::make_pair(cells[i].sample_name, cells[i].cell_name);
            }
        };

        if (!io::write_counts_gzip(conf.f_out.string(),
                                   counts,
                                   bins,
                                   hdr->target_name,
                                   sample_cell_name_wrapper(cells)) )
            return 1;
    }
    
    return 0;
}

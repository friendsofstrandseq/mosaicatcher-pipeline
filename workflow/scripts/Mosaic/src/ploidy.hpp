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

using interval::Interval;
using count::TGenomeCounts;
using count::Counter;

std::vector<std::string> get_states(unsigned ploidy)
{
    assert(ploidy>0 && ploidy < 6);
    std::vector<std::string> states;
    for (unsigned len=1; len < ploidy+1; ++len) {
        for (unsigned w = 0; w <= len; ++w) {

            unsigned c = len - w;
            std::string s;
            for (unsigned i = 0; i < w; ++i)
                s.push_back('W');
            for (unsigned i = 0; i < c; ++i)
                s.push_back('C');
            states.push_back(std::move(s));
        }
    }
    return states;
}


void setup_HMM_emissions(hmm::HMM<unsigned, hmm::MultiVariate<hmm::NegativeBinomial>> & hmm,
                         unsigned ploidy,
                         double nb_p,
                         double nb_r,
                         double nb_a,
                         double prior,
                         unsigned max_ploidy = 4)
{
    assert(max_ploidy>0 && max_ploidy<6);
    std::vector<hmm::MultiVariate<hmm::NegativeBinomial>> emissions;

    double prior_others = (1-prior)/static_cast<double>(max_ploidy);

    //std::cout << "MultiVariate<NegativeBinomial>:" << std::endl;
    for (unsigned len=1; len < max_ploidy+1; ++len) {
        for (unsigned w = 0; w <= len; ++w) {
            // e.g. WWWW (with N=4) should be  NB[p,a] x NB[p,(1-a)r]
            // or   WWWC                       NB[p,3/4r] x NB[p,1/4r]
            // but  WW ??        I decide for  NB[p,2/4r] x NB[p,a]
            unsigned c = len - w;
            double fac_w = (w==0 ? nb_a : w/static_cast<double>(max_ploidy) );
            double fac_c = (c==0 ? nb_a : c/static_cast<double>(max_ploidy) );
            hmm::MultiVariate<hmm::NegativeBinomial> distribution
                ({  hmm::NegativeBinomial(nb_p, fac_w * nb_r),
                    hmm::NegativeBinomial(nb_p, fac_c * nb_r)
                }, (len == ploidy ? prior : prior_others));
            emissions.push_back(std::move(distribution));

            // DEBUG
            //std::cout << std::setprecision(3) << "\t" << "p = " << nb_p << "\tcell_mean = " << nb_r/nb_p*(1-nb_p) << "\tState = ";
            //for (unsigned x = 0; x < w; ++x) std::cout << "W";
            //for (unsigned x = 0; x < c; ++x) std::cout << "C";
            //std::cout << "\tFactors = [" << fac_w << "\t" << fac_c << "]\t means = [" << fac_w * nb_r/nb_p*(1-nb_p) << "\t" <<  fac_c * nb_r/nb_p*(1-nb_p) << "]" << std::endl;

        }
    }
    //std::cout << "-" << std::endl;
    hmm.set_emissions(emissions);
}


void setup_HMM_emissions(hmm::HMM<unsigned, hmm::CombinedNegBinAndBinomial> & hmm,
                         unsigned ploidy,
                         double nb_p,
                         double nb_r,
                         double nb_a,
                         double prior,
                         unsigned max_ploidy)
{
    assert(max_ploidy>0 && max_ploidy<6);
    std::vector<hmm::CombinedNegBinAndBinomial> emissions;

    double prior_others = (1-prior)/static_cast<double>(max_ploidy);

    //std::cout << "CombinedNegBinAndBinomial:" << std::endl;
    for (unsigned len=1; len < max_ploidy+1; ++len) {
        for (unsigned w = 0; w <= len; ++w) {

            unsigned c = len - w;
            double ratio = (c == 0 || c == len ? (c == 0 ? 0.02 : 0.98) : static_cast<double>(c)/static_cast<double>(len) );
            double r     = (len == 0 ? nb_a : nb_r * static_cast<double>(len)/static_cast<double>(max_ploidy));
            emissions.push_back(hmm::CombinedNegBinAndBinomial(nb_p, r, ratio, (len == ploidy ? prior : prior_others)));

            // DEBUG
            //std::cout << std::setprecision(3) << "\t" << "p = " << nb_p << "\tcell_mean = " << nb_r/nb_p*(1-nb_p) << "\tState = ";
            //for (unsigned x = 0; x < w; ++x) std::cout << "W";
            //for (unsigned x = 0; x < c; ++x) std::cout << "C";
            //std::cout << "\t NB mean = " << r/nb_p*(1-nb_p) << "\tBinomial ratio = " << ratio << std::endl;

        }
    }
    //std::cout << "-" << std::endl;
    hmm.set_emissions(emissions);

    std::cout << "NB params for ploidy " << ploidy << ": p = " << nb_p << "\tr = " << nb_r << "\tmean = " << nb_r/nb_p*(1-nb_p) << std::endl;   
}






template <typename TDistribution>
typename hmm::HMM<unsigned, TDistribution> setup_HMM(unsigned max_ploidy, double p_trans)
{
    assert(max_ploidy>0 && max_ploidy<6);
    std::vector<std::string> states = get_states(max_ploidy);
    std::vector<double> initials(states.size(), 1/static_cast<double>(states.size()));
    std::vector<double> transitions(states.size() * states.size(), p_trans);
    for (unsigned i = 0; i < states.size(); ++i) {
        transitions[i * states.size() + i] = 1 - (states.size() - 1) * p_trans;
    }
    hmm::HMM<unsigned, TDistribution> hmm(states);
    hmm.set_initials(initials);
    hmm.set_transitions(transitions);
    return hmm;
}


template <typename TDistribution>
void run_generic_HMM(std::vector<TGenomeCounts> & counts,
                     std::vector<unsigned> const & good_cells,
                     std::vector<CellInfo>  & cells,
                     std::vector<unsigned> const & good_bins,
                     std::vector<int32_t> const & good_map,
                     std::unordered_map<std::string, SampleInfo> const & samples,
                     float p_trans,
                     double prior,
                     unsigned ploidy = 2,
                     unsigned max_ploidy = 4)
{
    // Set up and run HMM:
    hmm::HMM<unsigned, TDistribution> hmm = setup_HMM<TDistribution>(max_ploidy, p_trans);

    //std::cout << "HMM states: ";
    //for (auto x: hmm.state_labels) {std::cout << x << "\t";}
    //std::cout << std::endl;

    for (auto i = good_cells.begin(); i != good_cells.end(); ++i)
    {
        // set NB(n,p) parameters according to `p` of sample and mean of cell.
        float p = samples.at(cells[*i].sample_name).p;
        float r = (float)cells[*i].mean_bin_count * p / (1-p);
        float a = 0.05;

        setup_HMM_emissions(hmm, ploidy, p, r, a, prior, max_ploidy);
        run_HMM(hmm, counts[*i], good_bins, good_map);
    }

}








struct Conf_ploidy {
    std::vector<boost::filesystem::path> f_in;
    boost::filesystem::path f_out;
    boost::filesystem::path f_bins;
    boost::filesystem::path f_excl;
    boost::filesystem::path f_info;
    boost::filesystem::path f_segments;
    unsigned ploidy;
    unsigned max_ploidy;
    int minMapQual;
    unsigned int window;
    std::string model;
    double prior;
    Conf_ploidy() : max_ploidy(4)
    {}
};





int main_hmm(int argc, char **argv)
{

    // Command line options
    Conf_ploidy conf;
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
    ("ploidy,p", boost::program_options::value<unsigned>(&conf.ploidy)->default_value(2), "Assume cells have this ploidy level (max 4)")
    ("prior,P",  boost::program_options::value<double>(&conf.prior)->notifier(in_range(0,1,"prior")), "Prior probability for copy number <p> (penalize other copy number states. Default: All copy numbers are equally probable)")
    ("model,m", boost::program_options::value<std::string>(&conf.model)->default_value("multiNB"), "Models for HMM: multiNB, Binomial+NB")
    ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
    ("input-file", boost::program_options::value<std::vector<boost::filesystem::path> >(&conf.f_in), "input bam file(s)")
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
    if (conf.model != "multiNB" && conf.model != "Binomial+NB") {
        std::cerr << "[Error] Unknown --model for the HMM." << std::endl << std::endl;
        goto print_usage_and_exit;
    }
    if (vm.count("prior") && (conf.prior < 0 || conf.prior > 1)) {
        std::cerr << "[Error] Prior probability (for CN " << conf.ploidy << ") has to be between 0 and 1. By default all CNs will be equally likely" << std::endl << std::endl;
        goto print_usage_and_exit;
    }

    if (vm.count("help") || !vm.count("input-file"))
    {
    print_usage_and_exit:
        std::cout << std::endl;
        std::cout << "Mosaicatcher " << STRINGIFYMACRO(MOSAIC_VERSION_MAJOR);
        std::cout << "." << STRINGIFYMACRO(MOSAIC_VERSION_MINOR) << std::endl;
        std::cout << "> Count reads from Strand-seq BAM files..." << std::endl;
        std::cout << "  Now for different ploidy levels, too!" << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:   " << argv[0] << " --ploidy N [OPTIONS] <cell1.bam> <cell2.bam> ..." << std::endl << std::endl;
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
        if (vm.count("verbose")) std::cout << "[Info] Reading " << bins.size() << " variable-width bins with median bin size of " << round(median_binsize/1000) << "kb" << std::endl;
    }
    else
    {
        std::vector<Interval> exclude;
        if (vm.count("exclude")) {
            read_exclude_file(conf.f_excl.string(), hdr, exclude, vm.count("verbose"));
            sort(exclude.begin(), exclude.end(), interval::invt_less);
        }
        if (vm.count("verbose")) std::cout << "[Info] Creating " << round(conf.window/1000) << "kb bins with " << exclude.size() << " excluded regions" << std::endl;
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
    if (vm.count("verbose")) std::cout << "[Info] Reading " << conf.f_in.size() <<  " BAM files...";
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
        good_bins = count::get_good_bins(counts, cells, good_cells);
        if (vm.count("verbose")) std::cout << "[Info] Filtered out " << bins.size() - good_bins.size() << " bad bins" << std::endl;
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





    //
    // Chapter: Run HMM
    // ================
    //
    double prior = (vm.count("prior") ? conf.prior : 1.0/static_cast<double>(conf.max_ploidy));
    std::cout << "[Info] Running HMM with model " << conf.model << " and expected ploidy ";
    std::cout << std::setprecision(3) << conf.ploidy << " (prior = " << prior << ")" << std::endl;
    if (conf.model == "multiNB")
        run_generic_HMM<hmm::MultiVariate<hmm::NegativeBinomial>>(
                    counts,
                    good_cells,
                    cells,
                    good_bins,
                    good_map,
                    samples,
                    10.0f / bins.size(),
                    prior,
                    conf.ploidy,
                    conf.max_ploidy);
    if (conf.model == "Binomial+NB")
        run_generic_HMM<hmm::CombinedNegBinAndBinomial>(
                    counts,
                    good_cells,
                    cells,
                    good_bins,
                    good_map,
                    samples,
                    10.0f / bins.size(),
                    prior,
                    conf.ploidy,
                    conf.max_ploidy);



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

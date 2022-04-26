/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#ifndef simulate_hpp
#define simulate_hpp


#include <iostream>
#include <vector>
#include <random>
#include <utility>
#include <chrono>
#include <unordered_map>

//#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/multiprecision/random.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>

#include "version.hpp"
#include "program_options.hpp"
#include "intervals.hpp"
#include "counter.hpp"
#include "iocounts.hpp"

/**
 * @file
 * @defgroup simulation Simulation of Strand-seq data
 *
 * Summary of how Strand-seq data simulation works
 *
 * ## Strand-seq simulation
 *
 * The simulation of a single cell contains of three steps:
 *
 *   1. Sample basic counts for both haplotypes from a negative binomial distribution
 *   2. Insert SV according to the user specificaitons onto one of the haplotypes.
 *   3. Render cell as WW,WC or CC, potentially including SCEs
 *
 * ### 1: Basic counts for a cell
 * First, we randomly select a mean coverage (\f$c\f$) between `minCoverage` and
 * `maxCoverage` (user-specified) from a uniform distribution.
 *
 * Internally, a simulated cell is represented by two haplotypes (`h1`, `h2`), 
 * which both have reads on the main strand (I call it *plus*), and a few 
 * abnormal reads that are flipped in orientation (i.e. on the *minus* strand). 
 *
 *   * *plus* reads are sampled from a negative binomial with parameters \f$p\f$
 *     (user-specified) and \f$n = c/2 * p/(1-p)\f$
 *   * *minus* reads are sampled from a zero-inflated geometric distribution:
 *     with probability \f$1-\alpha\f$ the bin gets 0 reads, otherwise the read
 *     number is sampled from a geometric distribution
 *
 * ### 2: Insertion of SVs
 * SVs are specified in a config file and their start and end positions do **not
 * need to align with bins**.
 *
 * SVs are introduced by changing the *plus* and *minus* counts according to 
 * what you would expect from an SV. This can be either done for a single
 * haplotype (heterozygous) or both haplotypes (homozygous). Below the changes
 * to the strands are listed in detail
 *
 * SV type         | Haplotype *plus* counts | Haplotype *minus* counts
 * ---             | ---                     | ---
 * deletion  *HET* / *HOM*)     | set to 0                | no change
 * duplication (*HET* / *HOM*)  | multiply by 2           | no change
 * inversion (*HET* / *HOM*)    | switch with *minus*     | switch with *plus*
 * inverted duplication (*HET*) | no change               | set to *plus*
 * false_del (*HOM*)            | divide by 2             | no change
 *
 * > For now, heterozygous SV are **always introduced on haplotype 1**. In *HET*
 * > SVs the aforementioned changes are applied only to haplotype 1, in *HOM*
 * > SVs to both haplotypes.
 *
 * The changes are applied to all bins involved in the SV. However, if a bin is
 * only involved partially, the rule is applied only to a fraction of the
 * counts.
 *
 *
 * ### 3: Render cells
 *
 * At the last step, haplotypes are "inherited" as either Watson or Crick
 * strands. This is randomly chosen for each chromosome so that on average we
 * will obtain a 25:50:25 ratio of WW:WC:CC chromosomes.
 */



namespace simulator {


using interval::Interval;
using count::TGenomeCounts;
using count::Counter;


/**
 * Save read counts of a single bin for 2 haplotypes, separated by correct
 * strand (*plus*) and flipped strand (*minus*).
 *
 * @ingroup simulation
 */
struct HaploCount {
    unsigned h1_plus, h1_minus, h2_plus, h2_minus;
    HaploCount() : h1_plus(0), h1_minus(0), h2_plus(0), h2_minus(0)
    {}
};

/**
 * @ingroup simulation
 */
enum SV_type {
    het_inv,
    hom_inv,
    het_del,
    hom_del,
    het_dup,
    hom_dup,
    inv_dup,
    false_del
};
    std::string SV_type_to_string(SV_type x) {
        switch(x) {
            case het_inv: return ("het_inv");
            case hom_inv: return ("hom_inv");
            case het_del: return ("het_del");
            case hom_del: return ("hom_del");
            case het_dup: return ("het_dup");
            case hom_dup: return ("hom_dup");
            case inv_dup: return ("inv_dup");
            case false_del: return ("false_del");
        }
        return ("?");
    }

/**
 * @ingroup simulation
 */
struct SV {
    Interval where;
    SV_type type;
    float vaf;
    SV(Interval const & intvl, SV_type const & type) :
        where(intvl), type(type), vaf(1)
    {}
    SV(Interval const & intvl, SV_type const & type, float vaf) :
        where(intvl), type(type), vaf(vaf)
    {}
};


typedef std::vector<HaploCount> THapCount;
typedef std::vector<std::string> THapType;

struct phased_counts {
    uint8_t h1_w;
    uint8_t h1_c;
    uint8_t h2_w;
    uint8_t h2_c;
};

/**
 * Turn the list of HaploCount information into Strand-seq data.
 * @ingroup simulation
 *
 * Initially we decide for each haplotype on which strand (W or C) it is going
 * to be inherited. Then, while traversing along the chromosome, there is a 
 * small chance in every bin that these states change --> this is an SCE.
 *
 * **Update**: Now this function also simulates phased reads for both
 * haplotypes, which are drawn from a binomial distribution using `phased_frac`
 * as a probability.
 *
 * @param hapls Vector of haplotypes (THapl) for each cell, which shall be written as W/C counts.
 * @param chrom_map Chromosome boarders
 * @param sce_prob Probabiliy per bin to change strands
 * @param strand_states Vector of inherited strand states. Note that `Interval`s
 *        get mis-used by inputting bin numbers instead of chromosomal positions
 * @param phases Empty vector of `phased_counts` which will be filled to the
 *        same size as the returned `TGenomeCounts` with counts of haplotypes
 *        H1 and H2 on Watson and Crick strands. This simulates phase data.
 * @param phased_frac Fraction of reads that can be phased (used in binomial
 *        distribution).
 * @return Final Watson/Crick counts that can be plotted.
 */
template <typename TRandomDev>
TGenomeCounts render_cell(THapCount const & hapls,
                          std::vector<int32_t> const & chrom_map,
                          float sce_prob,
                          std::vector<std::pair<Interval, std::string>> & strand_states,
                          std::vector<phased_counts> & phases,
                          TRandomDev & rd_gen,
                          float phased_frac = 0.1)
{
    std::uniform_real_distribution<> rd_unif(0,1);

    // Final counts to be written
    TGenomeCounts counts(chrom_map.back());
    phases.resize(chrom_map.back());

    // Go through all chromosomes
    for (int32_t chrom = 0; chrom<chrom_map.size()-1; ++chrom)
    {
        if (chrom_map[chrom+1] - chrom_map[chrom] < 1) continue;

        // strand states for both haplotypes: true = Watson, false = Crick
        // Initially, choose states with equal prob.
        bool W_h1 = rd_unif(rd_gen) < 0.5;
        bool W_h2 = rd_unif(rd_gen) < 0.5;
        std::string state = std::string(W_h1?"W":"C") + std::string(W_h2?"W":"C");

        // Iterate over bins
        unsigned start_bin = chrom_map[chrom];
        for(unsigned bin = chrom_map[chrom]; bin < chrom_map[chrom+1]; ++bin)
        {
            // Small chance of an SCE:
            if(bin > chrom_map[chrom] && rd_unif(rd_gen) < sce_prob)
            {
                // Write down interval
                strand_states.push_back(std::make_pair(Interval(chrom, start_bin, bin-1), state));
                start_bin = bin;

                // change the state of one haplotype
                if (rd_unif(rd_gen) < 0.5) W_h1 = !W_h1;
                else                       W_h2 = !W_h2;
                state = std::string(W_h1?"W":"C") + std::string(W_h2?"W":"C");
            }

            // Fill counts
            counts[bin].watson_count = (W_h1  ? hapls[bin].h1_plus : hapls[bin].h1_minus) +
                                       (W_h2  ? hapls[bin].h2_plus : hapls[bin].h2_minus);
            counts[bin].crick_count  = (!W_h1 ? hapls[bin].h1_plus : hapls[bin].h1_minus) +
                                       (!W_h2 ? hapls[bin].h2_plus : hapls[bin].h2_minus);

            // Simulate phased reads (for each read there is a small chance to be phased)
            if (W_h1) {
                phases[bin].h1_w = std::binomial_distribution<>(hapls[bin].h1_plus, phased_frac)(rd_gen);
                phases[bin].h1_c = std::binomial_distribution<>(hapls[bin].h1_minus, phased_frac)(rd_gen);
            } else {
                phases[bin].h1_c = std::binomial_distribution<>(hapls[bin].h1_plus, phased_frac)(rd_gen);
                phases[bin].h1_w = std::binomial_distribution<>(hapls[bin].h1_minus, phased_frac)(rd_gen);
            }
            if (W_h2) {
                phases[bin].h2_w = std::binomial_distribution<>(hapls[bin].h2_plus, phased_frac)(rd_gen);
                phases[bin].h2_c = std::binomial_distribution<>(hapls[bin].h2_minus, phased_frac)(rd_gen);
            } else {
                phases[bin].h2_c = std::binomial_distribution<>(hapls[bin].h2_plus, phased_frac)(rd_gen);
                phases[bin].h2_w = std::binomial_distribution<>(hapls[bin].h2_minus, phased_frac)(rd_gen);
            }
        }

        // write down interval
        strand_states.push_back(std::make_pair(Interval(chrom, start_bin, chrom_map[chrom+1]-1), state));
    }
    return counts;
}




/** (Partially) flip haplotype counts according to a certail SV type
 *
 * For a het inversion for example, plus and minus on h1 are exchanged.
 * If a breakpoint of the SV does not perfectly align with bin boundaries,
 * you can specify an additional fraction f [0,1] to say that the flip should
 * only occur in f% of the bin.
 *
 * @param h Haplotype of a single bin (containing 4 counts).
 * @param sv_type Type of SV, only a few are possible.
 * @param f Apply the flip only to a portion of this bin.
 */
inline void flip_strand(HaploCount & h, SV_type sv_type, float f = 1)
{
    HaploCount x = h;
    switch(sv_type) {
        case het_inv:
            h.h1_plus  = (1-f) * h.h1_plus  + f * x.h1_minus;
            h.h1_minus = (1-f) * h.h1_minus + f * x.h1_plus;
            break;
        case hom_inv:
            h.h1_plus  = (1-f) * h.h1_plus  + f * x.h1_minus;
            h.h1_minus = (1-f) * h.h1_minus + f * x.h1_plus;
            h.h2_plus  = (1-f) * h.h2_plus  + f * x.h2_minus;
            h.h2_minus = (1-f) * h.h2_minus + f * x.h2_plus;
            break;
        case het_del:
            h.h1_plus  = (1-f) * h.h1_plus  + f * 0;
            h.h1_minus = (1-f) * h.h1_minus + f * 0;
            break;
        case hom_del:
            h.h1_plus  = (1-f) * h.h1_plus  + f * 0;
            h.h1_minus = (1-f) * h.h1_minus + f * 0;
            h.h2_plus  = (1-f) * h.h2_plus  + f * 0;
            h.h2_minus = (1-f) * h.h2_minus + f * 0;
            break;
        case het_dup:
            h.h1_plus  = (1-f) * h.h1_plus  + f * 2 * h.h1_plus;
            h.h1_minus = (1-f) * h.h1_minus + f * 2 * h.h1_minus; // should be 0
            break;
        case hom_dup:
            h.h1_plus  = (1-f) * h.h1_plus  + f * 2 * h.h1_plus;
            h.h1_minus = (1-f) * h.h1_minus + f * 2 * h.h1_minus; // should be 0
            h.h2_plus  = (1-f) * h.h2_plus  + f * 2 * h.h2_plus;
            h.h2_minus = (1-f) * h.h2_minus + f * 2 * h.h2_minus; // should be 0
            break;
        case inv_dup:
            h.h1_plus  = (1-f) * h.h1_plus  + f * (x.h1_minus + x.h1_plus);
            h.h1_minus = (1-f) * h.h1_minus + f * (x.h1_minus + x.h1_plus);
            break;
        case false_del:
            h.h1_plus  = (1-f) * h.h1_plus  + f * x.h1_plus/2;
            h.h1_minus = (1-f) * h.h1_minus + f * x.h1_minus/2;
            h.h2_plus  = (1-f) * h.h2_plus  + f * x.h2_plus/2;
            h.h2_minus = (1-f) * h.h2_minus + f * x.h2_minus/2;
    }
}









/** Read the SV config file (5 columns)
 *
 * this function also checks that the SV definition is valid and within the
 * chromosome boarders. Note that this is important because this check will not
 * be done during `locate_partial_bins`.
 *
 * @param filename file name.
 * @param chrom_names vector of names of the chromosomes.
 * @param chrom_size  vector of chromosome sizes. Must match `chrom_names` in size.
 * @param sv_list List of SVs to be written (push_back is used)
 */
bool read_SV_config_file(std::string const & filename,
                         std::vector<std::string> const & chrom_names,
                         std::vector<int32_t> const & chrom_size,
                         std::vector<SV> & sv_list)
{
    std::ifstream interval_file(filename.c_str(), std::ifstream::in);
    if (interval_file.is_open()) {
        while (interval_file.good()) {
            std::string line;
            getline(interval_file, line);
            typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
            boost::char_separator<char> sep(" \t,;");
            Tokenizer tokens(line, sep);
            Tokenizer::iterator tokIter = tokens.begin();
            if (tokIter!=tokens.end())
            {
                std::string chrName = *tokIter++;

                // get chromosome id
                //int32_t tid = bam_name2id(hdr, chrName.c_str());
                int32_t tid = -1;
                for (int32_t i = 0; i < chrom_names.size(); ++i) {
                    if (chrom_names[i] == chrName) {
                        tid = i;
                        break;
                    }
                }
                if (tid >= 0 && tid < chrom_names.size())
                {
                    Interval ivl;
                    ivl.chr = tid;

                    if (tokIter == tokens.end()) {
                        std::cerr << "Warning: Invalid line: " << line << std::endl;
                        continue;
                    }
                    ivl.start = boost::lexical_cast<int32_t>(*tokIter++);
                    if (tokIter == tokens.end()) {
                        std::cerr << "Warning: Invalid line: " << line << std::endl;
                        continue;
                    }
                    ivl.end   = boost::lexical_cast<int32_t>(*tokIter++);

                    // check whether interval makes sens
                    if (ivl.start < 0 || ivl.start > ivl.end || ivl.end > chrom_size[tid]) {
                        std::cerr << "[Warning] Interval out of bounds " << line << std::endl;
                        continue;
                    }

                    // Get SV type and VAF
                    if (tokIter == tokens.end()) {
                        std::cerr << "Warning: Invalid line: " << line << std::endl;
                        continue;
                    }
                    std::string type_in = *tokIter++;
                    SV_type type_out;
                    if      (type_in == "het_inv")
                        type_out = het_inv;
                    else if (type_in == "hom_inv")
                        type_out = hom_inv;
                    else if (type_in == "het_del")
                        type_out = het_del;
                    else if (type_in == "hom_del")
                        type_out = hom_del;
                    else if (type_in == "het_dup")
                        type_out = het_dup;
                    else if (type_in == "hom_dup")
                        type_out = hom_dup;
                    else if (type_in == "inv_dup")
                        type_out = inv_dup;
                    else if (type_in == "false_del")
                        type_out = false_del;
                    else {
                        std::cerr << "[Warning] Unknown SV type. Ignored: " << line << std::endl;
                    }


                    // Variant allele frequency
                    if (tokIter == tokens.end()) {
                        std::cerr << "Warning: Invalid line: " << line << std::endl;
                        continue;
                    }
                    float sv_vaf = boost::lexical_cast<float>(*tokIter++);
                    assert(tokIter == tokens.end());

                    sv_list.push_back(SV(ivl, type_out, sv_vaf));
                } else {
                    std::cerr << "Warning: Chromosome not found: " << chrName << " in \"" << line << "\"" << std::endl;
                }
            }
        }
        interval_file.close();
    } else {
        std::cerr << "[Error] SV config file cannot be read: " << filename << std::endl;
        return false;
    }
    return true;
}




/**
 * Insert SVs onto haplotype counts.
 * @ingroup simulator
 *
 * This function does a couple of steps at the same time **for each SV**:
 *    * smaple which cells are supposed to carry the SV with a probability
 *      of *vaf* for each cell (this is probabilistic, so the actual number can
 *      differ from the expected vaf).
 *    * Insert SV onto haplotypes (currently just on h1), including
 *      *fractional bins* at the ends.
 *    * Note down ids of cell carrying the SV into `inserted_SVs` (refers to 
 *      `cells`).
 *    * Calculate `optimal_breakpoints`, which are **always the right end of the
 *      bins**. The bin adapts to the left or right if the real SV breakpoint
 *      is < or > 50% of the bin.
 *
 * @param haplotypes Haplotype Counts that will be edited according to `flip_strands`.
 * @param inserted_SVs List of SVs and carriers that are inserted into haplotypes (initially empty).
 * @param cells List of `CellInfo`, which include cell and sample names etc.
 * @param optimal_breakpoints Set of optimal bin positions (e.g. optimal segmentation).
 * @param sv_list List of SVs to be inserted.
 * @param bins List of Intervals (bins) across the genome.
 * @param rd_gen Just pass the random generator here to flip a coin for each carrier.
 */
template <typename TRandomGenerator>
void simulate_SVs(std::vector<THapCount> & haplotypes,
                  std::vector<std::vector<unsigned>> & inserted_SVs,
                  std::set<unsigned> & optimal_breakpoints,
                  std::vector<SV> const & sv_list,
                  std::vector<Interval> const & bins,
                  std::vector<int32_t> const & chrom_map,
                  TRandomGenerator & rd_gen)
{
    // Need a uniform dist. to sample carriers of the SV.
    std::uniform_real_distribution<> rd_unif(0,1);

    // Prepare entries in inserted_SVs for each SV
    inserted_SVs.resize(sv_list.size());

    for (unsigned j = 0; j < sv_list.size(); ++j)
    {
        // SV position in bin coordinates
        SV const & sv = sv_list[j];
        auto   sv_bins = interval::locate_partial_bins(sv.where, bins, chrom_map);
        int32_t   binl = sv_bins.first.first;
        int32_t   binr = sv_bins.first.second;
        float       fl = sv_bins.second.first;
        float       fr = sv_bins.second.second;

        // Try for each cell
        for (unsigned i = 0; i < haplotypes.size(); ++i)
        {
            if (rd_unif(rd_gen) < sv.vaf)
            {
                // List which cells got which SV.
                inserted_SVs[j].push_back(i);

                // There is only one bin --> fraction is fl - (1-fr)
                if (binl == binr)
                {
                    flip_strand(haplotypes[i][binl], sv.type, fl - (1-fr));
                    optimal_breakpoints.insert(binl); // right
                    if (binl > chrom_map[sv.where.chr])
                        optimal_breakpoints.insert(binl-1); // left

                } else {

                    // Partially flip first bin
                    flip_strand(haplotypes[i][binl], sv.type, fl);
                    if (fl > 0.5 && binl > chrom_map[sv.where.chr]) // left
                        optimal_breakpoints.insert(binl-1);
                    else
                        optimal_breakpoints.insert(binl);

                    // Partially flip last bin
                    flip_strand(haplotypes[i][binr], sv.type, fr);
                    if (fr < 0.5 && binr > chrom_map[sv.where.chr]) // right
                        optimal_breakpoints.insert(binr-1);
                    else
                        optimal_breakpoints.insert(binr);

                    // Completely flip all bins in between
                    for (unsigned bin = binl+1; bin < binr; ++bin) {
                        flip_strand(haplotypes[i][bin], sv.type);
                    }
                }
            }
        }
    }
}





} /* namespace */

struct Conf_simul {
    bool verbose;
    unsigned n_cells;
    unsigned window;
    boost::filesystem::path f_sv;
    boost::filesystem::path f_out;
    boost::filesystem::path f_phases;
    boost::filesystem::path f_sce;
    boost::filesystem::path f_fai;
    boost::filesystem::path f_svs;
    boost::filesystem::path f_segment;
    boost::filesystem::path f_info;
    unsigned seed;
    double p, min_cov, max_cov, alpha, phased_frac;
    unsigned sce_num;
    std::string sample_name;
};


int main_simulate(int argc, char **argv)
{

    using interval::Interval;
    using count::TGenomeCounts;
    using count::Counter;
    using simulator::SV;
    using simulator::HaploCount;
    using simulator::SV_type;
    using simulator::THapCount;
    using simulator::THapType;


    // Command line options
    Conf_simul conf;
    boost::program_options::options_description po_generic("Generic options");
    po_generic.add_options()
    ("help,?", "show help message")
    ("verbose,v", "tell me more")
    ("seed", boost::program_options::value<unsigned>(&conf.seed), "Random generator seed")
    ("window,w", boost::program_options::value<unsigned>(&conf.window)->default_value(100000)->notifier(in_range(1000,10000000,"window")), "window size of fixed windows")
    ("numcells,n", boost::program_options::value<unsigned>(&conf.n_cells)->default_value(10)->notifier(in_range(0,500,"numcells")), "number of cells to simulate")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&conf.f_fai), "Chrom names & length file. Default: GRch38")
    ;

    boost::program_options::options_description po_out("Output options");
    po_out.add_options()
    ("out,o",         boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("out.txt.gz"), "output count file")
    ("phases,P",      boost::program_options::value<boost::filesystem::path>(&conf.f_phases), "output phased reads into a file")
    ("sceFile,S",     boost::program_options::value<boost::filesystem::path>(&conf.f_sce), "output the positions of SCEs")
    ("variantFile,V", boost::program_options::value<boost::filesystem::path>(&conf.f_svs), "output SVs and which cells they were simulated in")
    ("segmentFile,U", boost::program_options::value<boost::filesystem::path>(&conf.f_segment), "output optimal segmentation according to SVs and SCEs.")
    ("info,i",        boost::program_options::value<boost::filesystem::path>(&conf.f_info), "Write info about samples")
    ("sample-name",   boost::program_options::value<std::string>(&conf.sample_name)->default_value("simulated"), "Use this sample name in the output")
    ;

    boost::program_options::options_description po_rand("Radnomization parameters");
    po_rand.add_options()
    ("nbinom_p,p",    boost::program_options::value<double>(&conf.p)->default_value(0.8,"0.8")->notifier(in_range(0.01,0.99,"nbinom")), "p parameter of the NB distirbution")
    ("minCoverage,c", boost::program_options::value<double>(&conf.min_cov)->default_value(10)->notifier(in_range(1,500,"minCoverage")), "min. read coverage per bin")
    ("maxCoverage,C", boost::program_options::value<double>(&conf.max_cov)->default_value(60)->notifier(in_range(1,500,"maxCoverage")), "max. read coverage per bin")
    ("alpha,a",       boost::program_options::value<double>(&conf.alpha)->default_value(0.1,"0.1")->notifier(in_range(0,1,"alpha")), "noise added to all bins: mostly 0, but for a fraction alpha drawn from geometrix distribution")
    ("scesPerCell,s", boost::program_options::value<unsigned>(&conf.sce_num)->default_value(4)->notifier(in_range(0,200,"scesPerCell")), "Average number of SCEs per cell")
    ("phasedFraction,z", boost::program_options::value<double>(&conf.phased_frac)->default_value(0.1)->notifier(in_range(0,1,"scesPerCell")), "Average number of SCEs per cell")
    ;

    boost::program_options::options_description po_hidden("Hidden options");
    po_hidden.add_options()
    ("sv_config_file", boost::program_options::value<boost::filesystem::path>(&conf.f_sv), "Config file for SVs (see details)")
    ;

    boost::program_options::positional_options_description po_positional;
    po_positional.add("sv_config_file", -1);

    boost::program_options::options_description po_cmdline_options;
    po_cmdline_options.add(po_generic).add(po_out).add(po_rand).add(po_hidden);
    boost::program_options::options_description po_visible_options;
    po_visible_options.add(po_generic).add(po_out).add(po_rand);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(po_cmdline_options).positional(po_positional).run(), vm);
    boost::program_options::notify(vm);

    conf.verbose = vm.count("verbose");

    if (vm.count("help") || !vm.count("sv_config_file") || !file_exists(conf.f_sv.string()))
    {
        std::cout << std::endl;
        std::cout << "Mosaicatcher " << STRINGIFYMACRO(MOSAIC_VERSION_MAJOR);
        std::cout << "." << STRINGIFYMACRO(MOSAIC_VERSION_MINOR) << std::endl;
        std::cout << "> Simulate binned Strand-seq data." << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:   " << argv[0] << " [OPTIONS] SV-conf-file" << std::endl << std::endl;
        std::cout << po_visible_options << std::endl;

        if (vm.count("help")) {
            std::cout << "Simulate binned Strand-seq cells including structural variants (SVs)" << std::endl;
            std::cout << "and sister chromatid exchange events (SCEs). Type, size, position and" << std::endl;
            std::cout << "frequency of SVs are specified by a config file. To not include SVs" << std::endl;
            std::cout << "specify an empty file. The SV config file is a tab-separated file with" << std::endl;
            std::cout << "5 columns (chrom, start, end, SV type, avg. freuqency)." << std::endl;
            std::cout << "The allowed SV types are" << std::endl;
            std::cout << "  - het_del, hom_del" << std::endl;
            std::cout << "  - het_dup, hom_dup" << std::endl;
            std::cout << "  - het_inv, hom_inv" << std::endl;
            std::cout << "  - inv_dup" << std::endl;
            std::cout << "  - false_del (to simulate lower mappability region)" << std::endl;
            std::cout << "SV breakpoints do not need to align with bin boundaries." << std::endl;
            std::cout << std::endl;
            std::cout << "Note: The frequency (5-th colum) is interpreted as an expected" << std::endl;
            std::cout << "      fraction of cells carrying the SV - the exact number of cells" << std::endl;
            std::cout << "      can eventually differ form that expectation." << std::endl;
        }
        return vm.count("help") ? 0 : 1;
    }

    // Error when coverage & p are not applicable to boost negative_binomial
    double nb_n_param = conf.min_cov/2 * conf.p / (1-conf.p);
    if (nb_n_param < 1) {
        std::cerr << "[Error] The specified values for p and minCov will cause a problem" << std::endl;
        std::cerr << "        for the negative binomial distribution. This is currently" << std::endl;
        std::cerr << "        a limitation of boost::random::negative_binomial_distribution," << std::endl;
        std::cerr << "        which only allows integer values > 0 as a value for n" << std::endl;
        std::cerr << "        This behaviour is different in boost::math::negative_binomial" << std::endl;
        std::cerr << "        and again different in the implementation used in R." << std::endl;
        std::cerr << "        std::negative_binomial is not giving an error but produces" << std::endl;
        std::cerr << "        wrong numbers!" << std::endl;
        std::cerr << "        This limitation is also the reason why background reads cannot be" << std::endl;
        std::cerr << "        modelled by an NB for now." << std::endl;
        std::cerr << "Please choose a higher p or minCov!" << std::endl;
        return 0;
    }


    std::chrono::steady_clock::time_point t1, t2;

    // global vars
    std::vector<Interval>       bins;
    std::vector<int32_t>        chrom_map;
    std::vector<int32_t>        chrom_sizes;
    std::vector<std::string>    chrom_names;
    std::vector<THapCount>      haplotypes;
    std::vector<THapType>       chrom_states;
    std::vector<CellInfo>       cells;

    // Read genome or use GRch38 by default
    if (vm.count("genome")) {
        std::cerr << "[Error]: reading a genome file is not implemented. Leave out this flag to use GRch38" << std::endl;
        return 2;
    } else {
        chrom_sizes = { \
            248956422, 242193529, 198295559, 190214555,
            181538259, 170805979, 159345973, 145138636,
            138394717, 133797422, 135086622, 133275309,
            114364328, 107043718, 101991189,  90338345,
            83257441,  80373285,  58617616,  64444167,
            46709983,  50818468, 156040895,  57227415 };     // GRCh38
        chrom_names = { \
            "chr1",  "chr2",  "chr3",   "chr4",  "chr5",  "chr6",
            "chr7",  "chr8",  "chr9",  "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
            "chr19", "chr20", "chr21", "chr22",  "chrX",  "chrY" };
    }

    // Generage bins & chrom_map
    chrom_map.resize(chrom_sizes.size());
    create_fixed_bins(bins, chrom_map, conf.window, std::vector<Interval>(), (int32_t)chrom_sizes.size(), chrom_sizes);
    chrom_map.push_back((int32_t)bins.size());

    // Random generator
    std::random_device rd;
    std::mt19937 rd_gen(rd());
    boost::random::mt19937 rd_gen_boost;
    if (vm.count("seed")) {
        rd_gen_boost.seed(conf.seed);
        rd_gen      .seed(conf.seed);
        std::cout << "[Info] Using random seed " << conf.seed << std::endl;
    }


    // Generate basic haplotype counts of each cells, including random noise
    t1 = std::chrono::steady_clock::now();
    std::uniform_real_distribution<> rd_cov(conf.min_cov, conf.max_cov);
    std::uniform_real_distribution<> rd_unif(0,1);
    double p = conf.p;

    if (conf.verbose) std::cout << "Simulating  " << conf.n_cells << " cells" << std::endl;
    for (unsigned i = 0; i < conf.n_cells; ++i)
    {
        double cov_per_bin = rd_cov(rd_gen);
        std::geometric_distribution<>         rd_geom(5/(5+log2(cov_per_bin)));
        //std::negative_binomial_distribution<> rd_nb(cov_per_bin/2 * p/(1-p), p); // gives wrong results for small r !!
        boost::random::negative_binomial_distribution<> rd_nb_boost(cov_per_bin/2 * p/(1-p), p);
        boost::random::variate_generator<boost::mt19937&, boost::random::negative_binomial_distribution<> > rd_nb(rd_gen_boost, rd_nb_boost);


        CellInfo cell;
        cell.median_bin_count = static_cast<unsigned>(cov_per_bin);
        cell.sample_name = conf.sample_name;
        cell.cell_name   = std::string("cell_") + std::to_string(i);
        cell.bam_file    = "no_file";
        cell.nb_p        = p;
        cell.nb_a        = 0;
        cell.nb_r        = cov_per_bin * p / (1-p);

        THapCount count(bins.size());
        for (unsigned bin = 0; bin < bins.size(); ++bin)
        {
            count[bin].h1_plus = rd_nb();
            count[bin].h2_plus = rd_nb();
            count[bin].h1_minus = (rd_unif(rd_gen)<conf.alpha) ? rd_geom(rd_gen) : 0;
            count[bin].h2_minus = (rd_unif(rd_gen)<conf.alpha) ? rd_geom(rd_gen) : 0;
        }
        haplotypes.push_back(std::move(count));
        cells.push_back(std::move(cell));
    }


    // SV part
    std::vector<SV> sv_list;
    read_SV_config_file(conf.f_sv.string(), chrom_names, chrom_sizes, sv_list);

    // Insert SVs
    std::vector<std::vector<unsigned>> inserted_SVs;
    std::set<unsigned> optimal_breakpoints;
    simulate_SVs(haplotypes,
                 inserted_SVs,
                 optimal_breakpoints,
                 sv_list,
                 bins,
                 chrom_map,
                 rd_gen);


    t2 = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    if (conf.verbose) std::cout << "[Info] Simulation took " << time.count() << " sec." << std::endl;


    // Write SV information (especially which cells contain the SVs)
    if (vm.count("variantFile"))
    {
        std::cout << "[Write] Variant summary: " << conf.f_svs.string() << std::endl;
        std::ofstream out(conf.f_svs.string());
        if (out.is_open())
        {
            out << "chrom\tstart\tend\tSV_type\tsample\tcell" << std::endl;
            for (unsigned j = 0; j < sv_list.size(); ++j)
            {
                SV const & sv = sv_list[j];
                for (unsigned carrier_id : inserted_SVs[j])
                {
                    out << chrom_names[sv.where.chr] << "\t";
                    out << sv.where.start << "\t" << sv.where.end << "\t";
                    out << SV_type_to_string(sv.type) << "\t";
                    out << cells[carrier_id].sample_name << "\t";
                    out << cells[carrier_id].cell_name << std::endl;
                }
            }
        } else {
            std::cerr << "[Warning] Cannot write to " << conf.f_svs.string() << std::endl;
        }
    }

    // Update optimal_breakpoints by one breakpoint at the end of each chrom.
    for (int32_t chrom = 1; chrom < chrom_map.size(); ++chrom)
        optimal_breakpoints.insert(chrom_map[chrom]-1);

    // Write optimal segmentation file, which includes SV breakpoints and SCE breakpoints.
    if (vm.count("segmentFile"))
    {
        std::cout << "[Write] Segmentation file: " << conf.f_segment.string() << std::endl;
        std::ofstream out(conf.f_segment.string());
        if (out.is_open())
        {
            out << "k\tchrom\tbps" << std::endl;
            for (auto const & bin : optimal_breakpoints)
            {
                int32_t chrom = std::upper_bound(chrom_map.begin(), chrom_map.end(), bin) - chrom_map.begin() -1;
                out << 0 << "\t" << chrom_names[chrom] << "\t";
                out << bin - chrom_map[chrom] << std::endl;
            }
        } else {
            std::cerr << "[Warning] Cannot write to " << conf.f_segment.string() << std::endl;
        }
    }
    

    // Turn haplotypes into TGenomeCounts and simulate SCEs
    t1 = std::chrono::steady_clock::now();
    std::vector<TGenomeCounts> final_counts;
    std::vector<std::pair<Interval,std::string>> strand_states;
    std::vector<unsigned> str_states_cells;
    std::vector<std::vector<simulator::phased_counts>> phases;

    for (unsigned i = 0; i < conf.n_cells; ++i)
    {
        unsigned cell_pos = strand_states.size();
        std::vector<simulator::phased_counts> phase;
        final_counts.push_back(render_cell(haplotypes[i],   // simulated counts
                                           chrom_map,       // chromosome boundaries
                                           (float)conf.sce_num / bins.size(), // sce_probability
                                           strand_states, // Intervals with inherited strand states
                                           phase,
                                           rd_gen,
                                           conf.phased_frac));

        // note down which cells belong to these intervals
        for (; cell_pos < strand_states.size(); ++cell_pos)
            str_states_cells.push_back(i);

        // save phased reads as additional output
        phases.push_back(phase);
    }
    t2 = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    if (conf.verbose) std::cout << "[Info] Rendering cells took " << time.count() << " sec." << std::endl;



    // write down SCEs
    if (vm.count("sceFile"))
    {
        std::cout << "[Write] Inherited strand states (SCEs) to file " << conf.f_sce.string() << std::endl;
        std::ofstream out(conf.f_sce.string());
        if (out.is_open()) {
            out << "sample\tcell\tchrom\tstart\tend\tclass" << std::endl;
            for (unsigned i = 0; i < strand_states.size(); ++i) {
                out << conf.sample_name << "\t";
                out << "cell_" << std::to_string(str_states_cells[i]) << "\t";
                out << chrom_names[strand_states[i].first.chr] << "\t";
                out << bins[(strand_states[i].first).start].start << "\t";
                out << bins[(strand_states[i].first).end].end << "\t";
                out << strand_states[i].second << std::endl;
            }
        } else {
            std::cerr << "[Warning] Cannot write to " << conf.f_sce.string() << std::endl;
        }
    }


    // Get total number of reads per cell and print cell information:
    for (unsigned i = 0; i < final_counts.size(); ++i) {
        for (unsigned bin = 0; bin < bins.size(); ++bin) {
            cells[i].n_mapped += final_counts[i][bin].watson_count + final_counts[i][bin].crick_count;
        }
    }

    //
    // Chapter: Filter cells and bins and run HMM

    t1 = std::chrono::steady_clock::now();

    // median per cell
    count::set_median_per_cell(final_counts, cells);

    // filter cells with low counts and set pass_qc = false for bad cells;
    std::vector<unsigned> good_cells;
    good_cells = count::get_good_cells(final_counts, cells);
    for (auto c : cells) c.pass_qc = false;
    for (unsigned cid : good_cells) cells[cid].pass_qc = true;
    if (cells.size() > good_cells.size())
        std::cout << "[Info] " << cells.size() - good_cells.size() << "/" << cells.size()
                  << " cells were deemed QC fail by HMM (which is purely based on coverage for now)" << std::endl;

    // filter bins with abnormal counts (not happening here)
    std::vector<unsigned> good_bins(bins.size());
    std::iota(good_bins.begin(), good_bins.end(), 0); // fill with 0,1,2,...

    // calculate cell means and cell variances, grouped by sample (not cell)
    std::unordered_map<std::string, SampleInfo> samples;
    calculate_new_cell_mean(samples, cells, final_counts, good_cells, good_bins);

    // Estimation of parameter p per sample
    for (auto it = samples.begin(); it != samples.end(); ++it) {
        SampleInfo & s = it->second;
        s.p = std::inner_product(s.means.begin(), s.means.end(), s.means.begin(), 0.0f) \
        / std::inner_product(s.means.begin(), s.means.end(), s.vars.begin(), 0.0f);
    }

    // Chapter: Run HMM
    run_standard_HMM(final_counts,
                     good_cells,
                     cells,
                     good_bins,
                     chrom_map,
                     samples,
                     10.0f / bins.size());

    t2 = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    if (conf.verbose) std::cout << "[Info] Running HMM took " << time.count() << " sec." << std::endl;


    // Write info files
    if (vm.count("info")) {
        std::cout << "[Write] Cell info to " << conf.f_info.string() << ". Note that NB parameters are estimated before SVs are introduced." << std::endl;
        write_cell_info(conf.f_info.string(), cells);
    }




    // write down phases
    if (vm.count("phases"))
    {
        std::cout << "[Write] Phases to " << conf.f_phases.string() << std::endl;
        std::ofstream out(conf.f_phases.string());
        if (out.is_open())
        {
            out << "chrom\tstart\tend\tsample\tcell\th1_w\th1_c\th2_w\th2_c" << std::endl;
            for (unsigned i = 0; i < phases.size(); ++i)
            {
                for (unsigned j = 0; j < phases[i].size(); ++j)
                {
                    out << chrom_names[bins[j].chr] << "\t";
                    out << bins[j].start << "\t";
                    out << bins[j].end << "\t";
                    out << cells[i].sample_name << "\t";
                    out << cells[i].cell_name << "\t";
                    out << static_cast<int>(phases[i][j].h1_w) << "\t";
                    out << static_cast<int>(phases[i][j].h1_c) << "\t";
                    out << static_cast<int>(phases[i][j].h2_w) << "\t";
                    out << static_cast<int>(phases[i][j].h2_c) << std::endl;
                }
            }
        } else {
            std::cerr << "[Warning] Cannot write to " << conf.f_phases.string() << std::endl;
        }
    }

    


    // write down counts
    std::cout << "[Write] Count table " << conf.f_out.string() << std::endl;
    std::vector<std::pair<std::string, std::string>> sample_cell_names;
    for (unsigned i=0; i<conf.n_cells; ++i) {
        sample_cell_names.push_back(std::make_pair(conf.sample_name, "cell_" + std::to_string(i)));
    }
    io::write_counts_gzip(conf.f_out.string(), final_counts, bins, chrom_names, sample_cell_names);

    return 0;
}





#endif /* simulate_hpp */

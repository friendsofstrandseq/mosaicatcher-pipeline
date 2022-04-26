/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#ifndef counter_hpp
#define counter_hpp

#include <algorithm>
#include <unordered_map>

#include <htslib/sam.h>
#include "utils.hpp"
#include "intervals.hpp"

namespace count {



using interval::Interval;


/**
 * @ingroup count
 */
template <typename TInner>
struct Counter {
    //static const std::vector<std::string> label_names;
    //static const std::map<std::string, uint8_t> label_id;
    TInner watson_count, crick_count;
    std::string label;
    unsigned n_supplementary;

    Counter() : watson_count(static_cast<TInner>(0)), crick_count(static_cast<TInner>(0)), label("None"), n_supplementary(0)
    {}

    bool set_label(std::string const & s) {
        label = s;
        return true;
    }

    std::string const & get_label() const {
        return label;
    }
};

/**
 * @ingroup count
 */
typedef std::vector<Counter<unsigned>> TGenomeCounts;


/**
 * Calculate median bin count and return list of cells above.
 * @ingroup count
 *
 * **Important:** `cell.median_per_bin` must have been set beforehands!
 *
 * @param counts Matrix witch raw cell counts - labels are set to "None" if cell has too few reads.
 * @param min_median Cutoff for the minimum median read count (dafault: 3).
 */
    std::vector<unsigned> get_good_cells(std::vector<TGenomeCounts> & counts,
                                         std::vector<CellInfo> & cells,
                                         unsigned min_median = 3)
{
    assert(counts.size() == cells.size());

    std::vector<unsigned> good_cells;
    for (unsigned i=0; i<counts.size(); ++i) {
        if (cells[i].median_bin_count >= min_median)
            good_cells.push_back(i);
        else {
            cells[i].pass_qc = false;
            for (unsigned bin=0; bin < counts[i].size(); ++bin)
                counts[i][bin].set_label("None");
        }
    }
    return good_cells;
}


/**
 * Write median per cell into CellInfo.
 * @ingroup count
 *
 * @param counts matrix of counts.
 * @param cells CellInfo to be written.
 */
void set_median_per_cell(std::vector<TGenomeCounts> const & counts,
                         std::vector<CellInfo> & cells)
{
    assert(counts.size() == cells.size());

    for (unsigned i = 0; i < counts.size(); ++i) {
        TMedianAccumulator<unsigned int> med_acc;
        for (Counter<unsigned> const & count_bin : counts[i])
            med_acc(count_bin.watson_count + count_bin.crick_count);
        cells[i].median_bin_count = boost::accumulators::median(med_acc);
    }
}




/**
 * Filter bins across all cells.
 * @ingroup count
 *
 * @param counts Matrix witch raw cell counts.
 * @param cells Vector of CellInfos.
 * @param good_cells Ignore all other cells, e.g. the ones with too few reads.
 */
std::vector<unsigned> get_good_bins(std::vector<TGenomeCounts> const & counts,
                                    std::vector<CellInfo> const & cells,
                                    std::vector<unsigned> const & good_cells,
                                    bool verbose = 0,
                                    bool filter_by_WC = true)
{
    if (counts.size() < 1) return {};

    unsigned N = counts[0].size();

    // Median-Normalized counts
    std::vector<std::vector<std::tuple<float,float>>> norm_counts(counts.size());
    for (auto i = good_cells.begin(); i != good_cells.end(); ++i) {
        norm_counts[*i] = std::vector<std::tuple<float,float>>(N);
        for (unsigned bin = 0; bin < N; ++bin)
            norm_counts[*i][bin] = std::make_tuple(counts[*i][bin].watson_count/(float)cells[*i].median_bin_count,
                                                   counts[*i][bin].crick_count /(float)cells[*i].median_bin_count);
    }

    // mean + variance per bin
    std::vector<float> bin_means(N);
    std::vector<float> bin_variances(N);
    std::vector<float> bin_WC_fraction(N, 0.0);

    if (good_cells.size() < 12) {
        std::cerr << "[Warning] Too few cells, I will not filter bins by WC/CC/WW states." << std::endl;
        filter_by_WC = false;
    }

    for (unsigned bin = 0; bin < N; ++bin) {
        TMeanVarAccumulator<float> meanvar_acc;
        for (auto i = good_cells.begin(); i != good_cells.end(); ++i)
            meanvar_acc(std::get<0>(norm_counts[*i][bin]) + std::get<1>(norm_counts[*i][bin]));
        bin_means[bin]     = boost::accumulators::mean(meanvar_acc);
        bin_variances[bin] = boost::accumulators::variance(meanvar_acc);

        // check how many cells show a WC sginal in this bin
        // TODO: this is currently based on fixed cut-offs!
        if (filter_by_WC) {
            for (auto i = good_cells.begin(); i != good_cells.end(); ++i) {
                float bin_WC_frac = std::get<0>(norm_counts[*i][bin]) / (std::get<0>(norm_counts[*i][bin]) + std::get<1>(norm_counts[*i][bin]));
                if (bin_WC_frac > 0.15 & bin_WC_frac < 0.85) {
                    bin_WC_fraction[bin] += 1;
                }
            }
            bin_WC_fraction[bin] /= (float)good_cells.size();
        }
    }

    // finding good bins
    // to do: this could be extended by checking if bin is always in WC state!
    TMeanVarAccumulator<float> meanvar_acc;
    std::vector<char> reasons(N,'.');
    std::vector<unsigned> good_bins;

    for (unsigned bin = 0; bin < N; ++bin)
        meanvar_acc(bin_means[bin]);
    float mean_mean = boost::accumulators::mean(meanvar_acc);
    float mean_sd   = std::sqrt(boost::accumulators::variance(meanvar_acc));

    for (unsigned bin = 0; bin < N; ++bin)
    {
        if (bin_means[bin] < 0.1) {
            reasons[bin] = 'l'; // mean too low
        // } else if (bin_means[bin] > mean_mean + 4 * mean_sd) {
        //     reasons[bin] = 'h'; // mean too high
        } else if (filter_by_WC && bin_WC_fraction[bin] < 0.05) {
            reasons[bin] = 'w'; // never WC
        } else if (filter_by_WC && bin_WC_fraction[bin] > 0.95) {
            reasons[bin] = 'c'; // always WC
        } else {
            reasons[bin] = '_';
            good_bins.push_back(bin);
        }
    }

    if (verbose) {
        std::cout << "[Info] Reasons for filtered bins (each character is one bin):" << std::endl;
        std::cout << "[...]    l = low coverage across all cells" << std::endl;
        std::cout << "[...]    h = high coverage across all cells" << std::endl;
        std::cout << "[...]    c = always WC across cells" << std::endl;
        std::cout << "[...]    w = never WC across cells";
        unsigned x = 0;
        while(x < reasons.size()) {
            if (x % 80 == 0) std::cout << std::endl << "[...] ";
            std::cout << reasons[x++];
        }
        std::cout << std::endl;
    }
    return good_bins;
}


std::vector<unsigned> get_good_bins(std::vector<TGenomeCounts> const & counts,
                                    std::vector<CellInfo> const & cells,
                                    bool verbose = 0,
                                    bool filter_by_WC = true)
{
    // When no good_cells where supplied, use all cells
    std::vector<unsigned> good_cells(counts.size());
    for (unsigned i=0; i < good_cells.size(); ++i)
        good_cells[i] = i;

    return get_good_bins(counts, cells, good_cells);
}










/**
 * Count reads from BAM file.
 * @ingroup count
 *
 * Count start positions, which are expected to be sorted, into sorted bins.
 * Suitable for both fixed and variable-width bins.
 */
bool count_sorted_reads(std::string const & filename,
                        std::vector<Interval> const & bins,
                        std::vector<int32_t> const & chrom_map,
                        bam_hdr_t * hdr,
                        int min_map_qual,
                        TGenomeCounts & counts,
                        CellInfo & cell)
{

    // Open bam file
    samFile* samfile = sam_open(filename.c_str(), "r");
    if (samfile == NULL) {
        std::cerr << "Fail to open file " << filename << std::endl;
        return false;
    }
    hts_idx_t* idx = sam_index_load(samfile, filename.c_str());
    if (idx == NULL) {
        std::cerr << "Fail to open index for " << filename << std::endl;
        return false;
    }

    counts.resize(bins.size(), Counter<unsigned>());

    // access samfile chrom per chrom
    for (int32_t chrom = 0; chrom < hdr->n_targets; ++chrom) {

        // skip chromosomes with no bins
        if (chrom_map[chrom+1] - chrom_map[chrom] < 1)
            continue;

        int32_t bin = chrom_map[chrom];
        hts_itr_t* iter = sam_itr_queryi(idx, chrom, 0, hdr->target_len[chrom]);
        bam1_t* rec = bam_init1();
        while (sam_itr_next(samfile, iter, rec) >= 0) {

            if (rec->core.flag & BAM_FUNMAP) {
                continue;
            }
            ++cell.n_mapped;

            // expect pos to be sorted
            int32_t pos = rec->core.pos;

            // skip all bins left of this position.
            // Stop when all bins of the chromosome are done
            while (pos >= bins[bin].end)
                if (bin++ == chrom_map[chrom+1])
                    goto end_of_chromosome;

            // Ignore reads before until we reach the start of a bin
            if (pos < bins[bin].start)
                continue;

            assert(pos >= bins[bin].start && pos < bins[bin].end);

            // Don't read every RG tag because that might slow down BAM parsing.
            // auto x = bam_aux_get(rec, "RG");

            // Ignore certain alignments
            if (rec->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL)) {
                ++cell.n_supplementary;
                ++(counts[bin].n_supplementary);
                continue;
            } if (rec->core.flag & BAM_FDUP) {
                ++cell.n_pcr_dups;
                continue;
            } if ((rec->core.qual < min_map_qual) || (rec->core.tid<0)) {
                ++cell.n_low_mapq;
                continue;
            } if (rec->core.flag & BAM_FREAD2) {
                ++cell.n_read2s;
                continue;
            }
            ++cell.n_counted;
            if (rec->core.flag & BAM_FREVERSE)
                ++( counts[bin].watson_count );
            else
                ++( counts[bin].crick_count );
        }
    end_of_chromosome:
        bam_destroy1(rec);
        hts_itr_destroy(iter);
    }

    hts_idx_destroy(idx);
    sam_close(samfile);
    return true;
}


/**
 * Calculate mean and var. of all cells of a sample based on `good_bins` only.
 * @ingroup count
 *
 * Calculate mean and variance of reads per bin, but only use `good_bins` (this
 * make it different from `CellInfo.median_bin_count`). This value is written
 * into `CellInfo.mean_bin_count` and also appended to `SampleInfo.means`.
 * The values in `SampleInfo.means` and `SampleInfo.vars` will later be used 
 * to estimate the *p* parameter of the negative binomial.
 *
 * There is a version requiring `good_cells` and a version without (using all
 * cells)
 *
 * @param samples Map of cell name to sampleInfo --> variables `means` and 
 *        `vars` are updated.
 * @param cells List of `CellInfos`, which must have `sample_name` set. Also 
 *        member `mean_bin_count` will be set, but only for `good_cells`.
 * @param counts Matrix of raw counts.
 * @param good_cells Ignore all other cells.
 * @param good_bins Ignore all other bins. Also cell mean and var. is 
 *        re-calculated using only good_cells.
 */
bool calculate_new_cell_mean(std::unordered_map<std::string, SampleInfo> & samples,
                           std::vector<CellInfo> & cells,
                           std::vector<TGenomeCounts> const & counts,
                           std::vector<unsigned> const & good_cells,
                           std::vector<unsigned> const & good_bins)
{
    // calculate cell means and cell variances, grouped by sample (not cell)
    for (auto i = good_cells.begin(); i != good_cells.end(); ++i) {

        // Get mean and var for this cell, but only from good bins!
        TMeanVarAccumulator<float> acc;
        for (unsigned bini = 0; bini < good_bins.size(); ++bini) {
            acc(counts[*i][good_bins[bini]].crick_count + counts[*i][good_bins[bini]].watson_count);
        }
        // emplace finds key if existing and returns (it,false);
        // otherwise it inserts (key,value) and returns (it,true).
        auto it = samples.begin();
        std::tie(it, std::ignore) = samples.emplace(cells[*i].sample_name, SampleInfo());
        float cell_mean = boost::accumulators::mean(acc);
        float cell_var  = boost::accumulators::variance(acc);

        cells[*i].mean_bin_count = cell_mean;
        (it->second).means.push_back(cell_mean);
        (it->second).vars.push_back(cell_var);
    }
    return true;
}
bool calculate_new_cell_mean(std::unordered_map<std::string, SampleInfo> & samples,
                             std::vector<CellInfo> & cells,
                             std::vector<TGenomeCounts> const & counts,
                             std::vector<unsigned> const & good_bins)
{
    // calculate cell means and cell variances, grouped by sample (not cell)
    for (unsigned i = 0; i < counts.size(); ++i) {

        // Get mean and var for this cell, but only from good bins!
        TMeanVarAccumulator<float> acc;
        for (unsigned bini = 0; bini < good_bins.size(); ++bini) {
            acc(counts[i][good_bins[bini]].crick_count + counts[i][good_bins[bini]].watson_count);
        }
        // emplace finds key if existing and returns (it,false);
        // otherwise it inserts (key,value) and returns (it,true).
        auto it = samples.begin();
        std::tie(it, std::ignore) = samples.emplace(cells[i].sample_name, SampleInfo());
        float cell_mean = boost::accumulators::mean(acc);
        float cell_var  = boost::accumulators::variance(acc);

        cells[i].mean_bin_count = cell_mean;
        (it->second).means.push_back(cell_mean);
        (it->second).vars.push_back(cell_var);
    }
    return true;
}



} /* namespace */
#endif /* counter_hpp */

/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#ifndef segmentation_hpp
#define segmentation_hpp

#include <cassert>
#include <iostream>
#include <vector>
#include <limits>
#include <iomanip> // setw
#include <numeric> // partial_sum
#include <chrono>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>

#include "program_options.hpp"
#include "version.hpp"
#include "utils.hpp"
#include "iocounts.hpp"


using interval::Interval;
typedef std::vector<count::Counter<double>> TGenomeCountDouble;
using count::Counter;

/**
 * @file
 * @defgroup segmentation Segmentation algorithm
 * Summary of how segments across binned Strand-seq data of multiple cells are
 * found.
 *
 * ### Overview of the segmentation algorithm
 *
 * We apply an optimal multivariate segmentation algorithm to find segments
 * across all cells simultaneously.
 *
 * See `main_segment` for a description of which steps are done.
 *
 * ### Handling of None bins:
 *
 * There are three options to handle None bins
 *   - utilize-none   Treat None bins like any other bins (default)
 *   - penalize-none  Add a cost factor to the cost matrix during segmentation to
 *                    penalize segments that go through None bins. This should
 *                    force the algorithm to place change points at the boarders of
 *                    None stretches. It is currently a bit buggy.
 *   - remove-none    Remove None regions from the data prior to segmentation. The
 *                    results of the segmentation refer to the original bins again.
 *
 * @todo Write description
 */


template<class T>
using Matrix = std::vector<std::vector<T> >;

/**
 * @fn template <typename TMat> void print_mat(TMat const & G)
 * @ingroup segmentation
 * Print a Matrix (vector<vector>).
 *
 * @param G 2D matrix[row_index][column_index]
 */
template <typename TMat> void print_mat(TMat const & G) {
    for (auto row : G) {
        for (auto x : row)
            std::cout << std::setw(10) << x;
        std::cout << std::endl;
    }
}

template <typename Printable>
std::ostream& operator<< (std::ostream& stream, Matrix<Printable> const & m) {
    if (m.size() == 0 || m[0].size() == 0)
        return (stream << "[EMPTY MATRIX]" << std::endl);
    for (unsigned r = 0; r < m.size(); ++r) {
        //stream << std::setw(8) << std::setprecision(2) << m[r][0];
        for (unsigned c = 0; c < m[r].size(); ++c)
            stream << std::setw(10) << std::setprecision(4) << m[r][c];
        stream << std::endl;
    }
    return stream;
}



/**
 * @fn bool optimal_segment_dp(Matrix<double> const & cost, int max_cp, Matrix<int> & breakpoints, std::vector<double> & sse)
 * @ingroup segmentation
 * Find optimal segmentation based on a cost matrix.
 *  
 * This is a dynamic programming algorithm to find an optimal segmentation
 * in a cost matrix (see `calculate_cost_matrix`). It calculates a max_cp x N
 * matrix of optimal segmentation cost (internally) and breakpoints (returned
 * via `breakpoints`).
 *
 * ### Optimal cost matrix `dp` of size N x max_cp
 *
 * `dp[n, cp]` stores the cost of the optimal segmentation with
 * `cp` breakpoints of segment `[0, n]`. This is calculated via
 * dynamic programming from the optimal cost of a segmentation
 * of `[0, k-1]` with `cp-1` breakpoints (note that 0 <= k < n) plus
 * the actual cost of segment `[k, n]` from the "cost" matrix.
 *
 * The breakpoints listed in the final matrix are the right-most bin
 * of a segment in 1-based coordinates, or one right of the last one
 * in a 0-based coordinate system. In order to assign positions to the
 * change points, use:
 *    from = bins[ chrom_map[chrom] + (k == 0 ? 0 : breakpoints[cp][k-1]) ].start
 *    to   = bins[ chrom_map[chrom] + breakpoints[cp][k] - 1 ].end
 *
 * @param cost Cost matrix (see `calculate_cost_matrix`)
 * @param max_cp Maximum number of breakpoints (<= N)
 * @param breakpoints Breakpoints for k=1..max_cp will be written in here
 * @param sse sum of squared error value written here.
 */
bool optimal_segment_dp(Matrix<double> const & cost,
                  int max_cp,
                  Matrix<int> & breakpoints,
                  std::vector<double> & sse)
{
    // Determine max_k and N from cost matrix
    unsigned max_k  = (unsigned)cost.size();
    if (max_k < 2 || max_k > 5000) {
        std::cerr << "[Error] unreasonable value for max_k: " << max_k << std::endl;
        return false;
    }
    unsigned N      = (unsigned)cost[0].size();
    if (N < 2 || N > 25000) {
        std::cerr << "[Error] unreasonable value for N: " << N << std::endl;
        return false;
    }
    if (max_cp > N || max_k > N) {
        std::cerr << "[Error] max_cp or max_k cannot be larger than N: N=" << N << ", max_cp=" << max_cp << ", max_k=" << max_k << std::endl;
        return false;
    }

    // cost of optimal segmentation matrix
    Matrix<double> dp(N, std::vector<double>(max_cp, -1));

    // backtrack matrix mt
    Matrix<int>    mt(N, std::vector<int>(max_cp - 1, -1));

    // Initialization: dp[i][0] is just cost[i][0]
    for (unsigned r = 0; r < max_k; ++r)
        dp[r][0] = cost[r][0];
    for (unsigned r = max_k; r < N; ++r)
        dp[r][0] = -1;


    // Dynamic programming: column-wise (number of change points).
    for (unsigned cp=1; cp < max_cp; ++cp)
    {
        // row j = position j
        for (unsigned j=0; j<N; ++j)
        {
            double   z_min = -1;
            unsigned i_min =  j;
            int      j0    = (j<max_k) ? j : max_k;

            // Iterate over possible new mt between 0 and j
            for (unsigned k=0; k<j0; ++k)
            {
                // Best segmentation from 0 to j-k-1 with cp-1 change points
                double mI_prev_col = dp[j-k-1][cp-1];

                // Cost of segment from j-k to j
                double cost_seg    = cost[k][j-k];

                // Keep value with minimal cost and set mt matriv=x
                // Note that values < 0 represent positive infinity
                if (mI_prev_col >= 0 && cost_seg >= 0)
                {
                    if (z_min < 0 || mI_prev_col + cost_seg < z_min) {
                        z_min = mI_prev_col + cost_seg;
                        i_min = j-k;
                    }
                }
            } /* for k */
            dp[j][cp]   = z_min;
            mt[j][cp-1] = i_min;
        } /* for j */
    } /* for cp */


    // breakpoints: Row cp contains the breakpoints for a segmentation with cp breakpoints.
    // breakpoints[cp][cp] is always N
    sse     = std::vector<double>(max_cp);
    breakpoints = Matrix<int>(max_cp, std::vector<int>(max_cp, -1)); // important: initialize to -1

    for (unsigned cp = 0; cp < max_cp; ++cp)
    {
        // just write down SSE
        double z = dp[N-1][cp];
        //sse[cp] = (z > 0) ? -(double)N / 2.0 * (1 + log(2*M_PI) + log(z / N)) : -1e10;  // 1e10 as a really really low value
        sse[cp] = (z >= 0) ? z : 1e10;  // 1e10 as a really really high value

        // Backtrack to get breakpoints
        // i is always the bin right of the changepoint (i.e. last bin in interval is i-1)
        // For one segment (cp==1) --> this value is consequently N
        int i = (int)N;
        breakpoints[cp][cp] = N;
        for (int j = cp - 1; j >= 0; --j) {
            i = mt[i-1][j];
            breakpoints[cp][j] = i;
        }
    }
    return true;
}


/**
 * @fn Matrix<double> calculate_cost_matrix(std::vector<double> const & data, int max_k)
 * @ingroup segmentation
 * Calculate segmentation cost matrix on a single vector.
 *
 * Function re-implemented from Wolfgang Hubers tilingArray package
 * (https://github.com/Bioconductor-mirror/tilingArray/blob/master/R/costMatrix.R).
 * This is the version that runs on a single vector.
 *
 * @todo Profile this algorithm, to see where time is spent.
 * @todo Potentially replace calculate_cost_matrix_single by a faster version using
 *       boost UBLAS vectors and matrices.
 * 
 * @param data Single data vector.
 * @param max_k Maximum number of bins for each segment. max_k <= N.
 */
Matrix<double> calculate_cost_matrix(std::vector<double> const & data,
                                                int max_k)
{
    // See https://www.bioconductor.org/packages/devel/bioc/vignettes/tilingArray/inst/doc/costMatrix.pdf
    // for an explanation of the algebra

    // Rows 0 ... max_k => segment lengths (1 ... max_k+1)
    // Columns 0 ... N  => start position of segment

    double ZERO_THR = 1e-10;

    unsigned N = (unsigned)data.size();
    std::vector<double> cr(N);
    std::vector<double> cq(N);

    // cr = cumsum(data)
    std::partial_sum(data.begin(), data.end(), cr.begin());                             // O(N)

    // cq = cumsum(data^2)
    std::transform(data.begin(), data.end(), cq.begin(), [](double u){return u*u;});     // O(N)
    std::partial_sum(cq.begin(), cq.end(), cq.begin());                                 // O(N)

    // Cost matrix G
    Matrix<double> G(max_k, std::vector<double>(N, -1));                                  // O(N * maxk)

    // Initialization for segments of lengths 1...max_k at position 0.
    for (unsigned k = 0; k < max_k; ++k) {                                               // O(maxk)
        G[k][0] = cq[k] - cr[k]*cr[k]/(double)(k+1);
        if (G[k][0] < ZERO_THR) {
            if (G[k][0] < -ZERO_THR)
                std::cout << "[Internal] Large negative value in cost matrix: "
                          << "G[" << k << "][" << 0 << "] = " << G[k][0]
                          << " --> 0" << std::endl;
	    G[k][0] = 0;
	}
    }

    // Iterate per row (k)
    for (unsigned k = 1; k <= max_k && k <= N-1; ++k)                                   // O( maxk * ...
    {
        // values for a segment of length k
        std::vector<double> cqk(N-k, -1);                                                //           ... N
        std::vector<double> crk(N-k, -1);                                                //           ... N
        for (unsigned i = 0; i < N-k; ++i) {                                            //           ... N
            cqk[i] = cq[i+k] - cq[i];
            crk[i] = cr[i+k] - cr[i];
        }

        // Iterate columns j = 1...N-k
        for (unsigned j = 1; j < N-k+1; ++j) {                                          //           ... N
            unsigned i = j - 1;
            assert(cqk[i]>=0);
            assert(crk[i]>=0);
            G[k-1][j] = cqk[i] - crk[i]*crk[i]/(double)k;

            // force small values to zero
            if (G[k-1][j] < ZERO_THR) {
                if (G[k-1][j] < -ZERO_THR)
                    std::cout << "[Internal] Large negative value in cost matrix: "
                              << "G[" << k-1 << "][" << j << "] = " << G[k-1][j]
                              << " --> 0" << std::endl;
                G[k-1][j] = 0;
            }
        }
    }

    return G;                                                                           // TOTAL: O(N * maxk) <= O(N^2)
}



/**
 * @fn calculate_cost_matrix(Matrix<double> const & data, int max_k, Matrix<double> & G)
 * @ingroup segmentation
 * Calculate segmentation cost matrix.
 *
 * Function re-implemented from Wolfgang Hubers tilingArray
 * (https://github.com/Bioconductor-mirror/tilingArray/blob/master/R/costMatrix.R).
 * The difference to the tilingArray package is that, unlike there,
 * a separate piecewise-constant function is fitted to every sample
 * (in the original version all samples are averaged into one).
 * This is done by calculating the cost matrix seperately on each sample
 * and then averaging across them. Complexity is quite high now, with
 * O(N * maxk * J) where J is the number of strands or samples.
 *
 * @todo Test for correctness (e.g. write unit test)
 *
 * @param data Matrix of data values, each row representing a single sample 
 *             (or strand). Note that all must have the same length and should
 *             be normalized beforehands!
 * @param max_k Number of bins of longest possible segment. This reduces the
 *              complexity from O(N^2) to O(N*max_k).
 * @param G Cost matrix. Reference will be resized and filled.
 * 
 * @addtogroup segmentation
 */
bool calculate_cost_matrix(Matrix<double> const & data,                         // TOTAL O ( J * N * maxk)
                                    int max_k,
                                    Matrix<double> & G)
{

    unsigned J  = (unsigned)data.size();        // number of cells or strands
    if (J < 1 || J > 2000) {
        std::cerr << "[Warning] more than 2000 strands might be too much: " << J << std::endl;
        return false;
    }
    unsigned N      = (unsigned)data[0].size(); // length of signal
    if (N < 2 || N > 25000) {
        std::cerr << "[Error] unreasonable value for N: " << N << std::endl;
        return false;
    }
    if (max_k > N) {
        std::cerr << "[Error] max_k cannot be larger than N: N=" << N << ", max_k=" << max_k << std::endl;
        return false;
    }


    G = Matrix<double>(max_k, std::vector<double>(N, 0));
    for (unsigned s = 0; s < J; ++s) {

        // Get cost of this one sample
        Matrix<double> cost_j = calculate_cost_matrix(data[s], max_k);

        // Add to total cost
        for (unsigned i=0; i<max_k; ++i)
            for (unsigned j=0; j<N; ++j)
                G[i][j] += cost_j[i][j];

    }

    // Finally, divide by J
    for (unsigned i=0; i<max_k; ++i)
        for (unsigned j=0; j<N; ++j)
            G[i][j] /= (double)J;

    return true;
}


struct Conf_segment {
    boost::filesystem::path f_in;
    boost::filesystem::path f_out;
    boost::filesystem::path f_cost_mat;
    float max_bp_per_Mb;
    unsigned max_bp_intercept;
    float none_penalty;
    float merge_threshold;
    unsigned max_segment_length;
    bool remove_bad_cells;
    unsigned min_num_bins_per_segment;
    std::string cm_chrom;
};


/**
 * @fn main_segment(int argc, char** argv)
 * @ingroup segmentation
 * Main program for segmentation
 *
 * 1. Read Strand-seq counts from table, incl. 'class' column
 * 2. Remove cells which are all 'none' (opt-out).
 * 3. Detect stretches of removed bins (from 'None' in the counts table) and ...
 *    1. do nothing special,
 *    2. penalize these regions in the cost matrix, or
 *    3. remove these regions from the data before segmentation.
 * 4. Finally, run dynamic programming segmentation and report optimal segments
 *    for a different numbers of change points.
 */
int main_segment(int argc, char** argv) {

    Conf_segment conf;
    boost::program_options::options_description po_generic("Generic options");
    po_generic.add_options()
    ("help,?", "show help message")
    ("out,o", boost::program_options::value<boost::filesystem::path>(&conf.f_out)->default_value("segments.txt"), "output file for counts")
    ;

    boost::program_options::options_description po_segmentation("Segmentation options");
    po_segmentation.add_options()
    ("max_bp_per_Mb,m", boost::program_options::value<float>(&conf.max_bp_per_Mb)->default_value(0.25)->notifier(in_range(0,1,"max_bp_per_Mb")), "max. number of change points per Mb")
	("max_bp_intercept,i", boost::program_options::value<unsigned>(&conf.max_bp_intercept)->default_value(15)->notifier(in_range(0,50,"max_bp_intercept")), "max. number of cp, add this constant")
    ("max_segment,M", boost::program_options::value<unsigned>(&conf.max_segment_length)->default_value(100000000), "maximum segment length")
    ("penalize-none", boost::program_options::value<float>(&conf.none_penalty)->implicit_value(100), "Penalize segments through removed bins (which are marked by 'None' in the counts table).")
    ("forbid-small-segments", boost::program_options::value<unsigned>(&conf.min_num_bins_per_segment)->default_value(1)->notifier(in_range(1,20,"forbid-small-segments")), "Penalize segments shorter that this number of bins")
    ("remove-none", "Remove segments through removed bins before segmentation. Mutually exclusive with --penalize-none.")
    ("do-not-normalize-cells", "Instead of using cell-normalized counts for each strand, use the raw numbers.")
    ("do-not-remove-bad-cells", "Keep all cells (by default, cells which are marked 'None' in all bins get removed")
    ;

    boost::program_options::options_description po_hidden("Hidden options");
    po_hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&conf.f_in), "mosaicatcher count file")
    ("cost-matrix,c", boost::program_options::value<boost::filesystem::path>(&conf.f_cost_mat), "write cost matrix to file")
    ("cost-matrix-chrom,C", boost::program_options::value<std::string>(&conf.cm_chrom), "Chromosome for which cost matrix should be written (by default first chrom)")
    ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", 1);


    boost::program_options::options_description cmdline_options;
    cmdline_options.add(po_generic).add(po_segmentation).add(po_hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(po_generic).add(po_segmentation);
    boost::program_options::variables_map vm;

    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);


    // Check arguments
    if (vm.count("help") || !vm.count("input-file")) {

        std::cout << std::endl;
        std::cout << "Mosaicatcher " << STRINGIFYMACRO(MOSAIC_VERSION_MAJOR);
        std::cout << "." << STRINGIFYMACRO(MOSAIC_VERSION_MINOR) << std::endl;
        std::cout << "> Find a segmentation across multiple cells." << std::endl;
        std::cout << std::endl;
        std::cout << "Usage:   " << argv[0] << " [OPTIONS] counts.txt.gz" << std::endl << std::endl;
        std::cout << visible_options << std::endl;
        return vm.count("help") ? 0 : 1;
    }
    conf.remove_bad_cells = true;
    if (vm.count("do-not-remove-bad-cells"))
        conf.remove_bad_cells = false;
    if (vm.count("penalize-none") && vm.count("remove-none")) {
        std::cerr << "[Error] --penalize-none and --remove-none are mutually exclusive." << std::endl;
        return 1;
    }



    /////////////////////////////////////////////////////////// global variables
    /* counts/cells */
    std::vector<std::pair<std::string,std::string>> sample_cell_names;
    std::vector<TGenomeCountDouble> counts;

    /* bins */
    std::vector<std::string> chromosomes;
    std::vector<Interval> bins;
    std::vector<int32_t>  chrom_map;
    /////////////////////////////////////////////////////////// global variables




    // 1. Reading count file
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::cout << "[Info] Reading file " << conf.f_in.string() << " (this is a bit slow)..." << std::endl;
    if (!io::read_counts_gzip(conf.f_in.string(),
                              counts,
                              chromosomes,
                              sample_cell_names,
                              bins))
        return 1;
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "[Info] " << counts.size() << " cells found." << std::endl;
    std::cout << "[Time] Reading took " << time.count() << " seconds." << std::endl;


    // Catch case of empty table
    if (chromosomes.empty())
        return 0;
    if (conf.cm_chrom.empty()) conf.cm_chrom = chromosomes[0];


    // If more than one sample - give a warning
    std::set<std::string> samples;
    for (auto s_c : sample_cell_names) samples.insert(s_c.first);
    if (samples.size()>1)
        std::cerr << "[Warning] Found > 1 sample, but Samples will be ignored for now. All cells will be treated as being from the first sample" << std::endl;


    // Create chromosome map
    chrom_map = std::vector<int32_t>(chromosomes.size());
    if (!make_chrom_map(bins, chrom_map)) {
        std::cerr << "[Error] Something is wrong with the intervals" << std::endl;
        return 2;
    }
    chrom_map.push_back((int32_t)bins.size());
    std::cout << "[Info] " << bins.size() << " intervals across " << chromosomes.size() << " chromosomes." << std::endl;



    // 2. remove bad cells (i.e. with all None).
    std::vector<unsigned> good_cells = io::get_good_cells(counts);
    if (conf.remove_bad_cells) {
        unsigned prev = counts.size();

        for (unsigned i = 0; i < good_cells.size(); ++i) {
            if (i != good_cells[i]) {
                counts[i] = counts[good_cells[i]];
                sample_cell_names[i] = sample_cell_names[good_cells[i]];
            }
        }
        counts.resize(good_cells.size()); // should
        std::cout << "[Info] Removed " << prev - counts.size() << " bad cells, leaving " << counts.size() << " good ones" << std::endl;
    }



    // Determine window size & Mean count per sample (ignoring None bins)
    unsigned window_size;
    std::vector<float> mean_per_cell(counts.size());
    {
        TMeanVarAccumulator<float> mean_acc;
        for (auto bin : bins)
            mean_acc((float)(bin.end - bin.start));
        window_size = (unsigned) boost::accumulators::mean(mean_acc);

        for (unsigned i=0; i<counts.size(); ++i) {

            // Calculate mean without using 'None' bins
            TMeanVarAccumulator<float> mean_acc;
            for (auto j = 0; j < bins.size(); ++j)
                // None bins are ignored during calculation of mean counts!
                if (counts[i][j].label != "None")
                    mean_acc((counts[i][j]).watson_count + (counts[i][j]).crick_count);

            // if there was at least one value write down the mean, otherwise set it to 1.0
            if (boost::accumulators::count(mean_acc)==0) {
                std::cerr << "[Warning] Calculating mean of cell with only black-listed bins" << std::endl;
                mean_per_cell[i] = 1.0;
            } else {
                mean_per_cell[i] = boost::accumulators::mean(mean_acc);
            }
        }
    }




    // 3. Determine `good_bins` from stretches of 'None' in the data.
    // This is more flexible than using low-counts, because it allows
    // the user to alter the labels of certain regions.
    //
    std::vector<std::vector<std::pair<unsigned,unsigned>>> none_reg_local(chromosomes.size());
    std::vector<unsigned> good_bins;
    //std::vector<int32_t> good_map;
    {

        // Finding 'none' stretches
        unsigned num_none_bins = 0;
        for (int32_t chrom = 0; chrom < chromosomes.size(); ++chrom) {

            unsigned pos = chrom_map[chrom];
            while (pos < chrom_map[chrom+1]) {
                if (counts[0][pos].label != "None") {
                    good_bins.push_back(pos);
                    ++pos;
                    continue;
                }
                // iterate through consecutive stretch of None
                unsigned start = pos;
                while (pos < chrom_map[chrom+1] && counts[0][pos].label == "None")
                    ++pos;
                none_reg_local[chrom].push_back(std::make_pair(start - chrom_map[chrom], pos - 1 - chrom_map[chrom]));
                num_none_bins += pos - start;

                // adjust boundaries before removing None stretches...
                if (vm.count("remove-none")) {
                    // first make sure that the current None stretch does not span the entire chromosome
                    if (start == chrom_map[chrom] && pos == chrom_map[chrom+1]) {
                        std::cerr << "None stretch seems to span the whole chromosome" << std::endl;
                    }
                    // if this None stretch is not at the start of the chromosome ...
                    if (start != chrom_map[chrom]) {
                        // add this None stretch to the interval left of it
                        bins[start-1].end = bins[pos-1].end;
                        std::cout << "[Debug] merge " << pos - start << " bins to " << chromosomes[bins[start].chr] << "[" << bins[start-1].start << "-" << bins[start-1].end << "]" << std::endl;
                    }
                    // If it is at the start, though, count it th the interval on the right
                    else {
                        bins[pos].start = bins[start].start;
                        std::cout << "[Debug] merge " << pos - start << " bins to " << chromosomes[bins[pos].chr] << "[" << bins[pos].start << "-" << bins[pos].end << "]" << std::endl;
                    }
                    // todo: when outputting the adjusted bins, the outputted indices refer to this subset of bins.
                    // They should, however, refer to the original bins (prior to `None` removal.
                }
            }
        std::cout << "[Info] Found " << none_reg_local[chrom].size() << " 'None' stretches on " << chromosomes[chrom] << ", in total " << num_none_bins << " bins." << std::endl;
        }


        // Alternative 3.3
        //     "--remove-none"
        // Remove 'none' bins from data before running segmentation.
        if (vm.count("remove-none")) {

            // remove bins
            for (unsigned j = 0; j < good_bins.size(); ++j)
                bins[j] = bins[good_bins[j]];
            bins.resize(good_bins.size());

            // remove data values
            for (unsigned i = 0; i < counts.size(); ++i) {
                for (unsigned j = 0; j < good_bins.size(); ++j)
                    counts[i][j] = counts[i][good_bins[j]];
                counts[i].resize(good_bins.size());
            }

            // remake chrom_map
            chrom_map.resize(good_bins.size());
            make_chrom_map(bins, chrom_map);
            chrom_map.push_back(good_bins.size());
        }

    }





    // prepare OUTPUT file
    std::ofstream out(conf.f_out.string());
    out << "# Breakpoint file generated via this command: " << std::endl;
    out << "# $ ";
    for (unsigned x = 0; x < argc; ++x) out << argv[x] << " ";
    out << std::endl;
    out << "# sample       Sample name" << std::endl;
    out << "# cells        Number of cells used in segmentation" << std::endl;
    out << "# . chrom      Segmented chromosome" << std::endl;
    out << "# . bins       Number of bins during segmentation. This can be lower " << std::endl;
    out << "#              than total number of bins on this chromosome (--remove-none)" << std::endl;
    out << "# . maxcp      Max. number of change points used (max_bp_intercept + chrom_size/1e6 * max_bp_per_Mb)" << std::endl;
    out << "# . maxseg     Max. segment length in number of bins (--max_segment)" << std::endl;
    out << "# . none_bins  Number of `None` bins in this chrom" << std::endl;
    out << "# . none_regs  Number of consecutive `None` regions in this chrom" << std::endl;
    out << "# . action     How to treat `None` bins (--remove-none, --penalize-none)" << std::endl;
    out << "# ... k        Current max. number of change points" << std::endl;
    out << "# ... sse      Standard squared error associated with the current segmentation" << std::endl;
    out << "# ..... bps    Index of last bin in the segment (0-based)" << std::endl;
    out << "# ..... start  Start position [bp] of the segment defined by bps" << std::endl;
    out << "# ..... end    End position [bp] of the segment defined by bps" << std::endl;

    out << "sample\t";
    out << "cells\t";
    out << "chrom\t";
    out << "bins\t";
    out << "maxcp\t";
    out << "maxseg\t";
    out << "none_bins\t";
    out << "none_regs\t";
    out << "action\t";
    out << "k\t";
    out << "sse\t";
    out << "bps\t";
    out << "start\t";
    out << "end" << std::endl;


    // Segmentation chromosome per chromosome
    // This loop can later be parallelized
    for (int32_t chrom=0; chrom < chromosomes.size(); ++chrom)
    {
        std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();

        // Put both strands of each cell into data.
        // Normalize each cell by its mean coverage first, by default,
        // but this can be switched off.
        unsigned N = chrom_map[chrom+1] - chrom_map[chrom];
        if (N <= 1) continue;

        Matrix<double> data;
        for (unsigned i = 0; i < counts.size(); ++i)
        {
            if (!vm.count("do-not-normalize-cells"))
            {
                if (i==0 && chrom == 0) std::cout << "[Info] Normalize cells by mean counts (excl. None regions)" << std::endl;
                double cell_mean = mean_per_cell[i];
                std::vector<double> tmp(N);
                std::transform(counts[i].begin() + chrom_map[chrom],
                               counts[i].begin() + chrom_map[chrom] + N,
                               tmp.begin(),
                               [cell_mean](Counter<double> const & c){return (double) c.watson_count / cell_mean;});
                data.push_back(tmp);
                std::transform(counts[i].begin() + chrom_map[chrom],
                               counts[i].begin() + chrom_map[chrom] + N,
                               tmp.begin(),
                               [cell_mean](Counter<double> const & c){return (double) c.crick_count / cell_mean;});
                data.push_back(tmp);
            } else {
                std::vector<double> tmp(N);
                std::transform(counts[i].begin() + chrom_map[chrom],
                               counts[i].begin() + chrom_map[chrom] + N,
                               tmp.begin(),
                               [](Counter<double> const & c){return (double) c.watson_count;});
                data.push_back(tmp);
                std::transform(counts[i].begin() + chrom_map[chrom],
                               counts[i].begin() + chrom_map[chrom] + N,
                               tmp.begin(),
                               [](Counter<double> const & c){return (double) c.crick_count;});
                data.push_back(tmp);
            }
        }


        // parameters:
        // max_k = longest allowed segment
        // max_cp = max. number of change points, at least 10
        unsigned chrom_size = bins[chrom_map[chrom+1]-1].end - bins[chrom_map[chrom]].start;
        unsigned max_k  = std::min(std::max(static_cast<unsigned>(10),
                                            static_cast<unsigned>(conf.max_segment_length/window_size)),
                                   static_cast<unsigned>(N));
        unsigned max_cp = std::min(conf.max_bp_intercept + static_cast<unsigned>(ceil( (float)chrom_size/1e6 * conf.max_bp_per_Mb )),  N-1);


        // New Cost matrix
        Matrix<double> new_cost;
        if (!calculate_cost_matrix(data, max_k, new_cost)) {
            std::cout << "[Error] Cost matrix calculation failed on " << chromosomes[chrom] << std::endl;
            continue;
        }

        // also, check that values make sense
        for (unsigned k=0; k<max_k; ++k)
            for (unsigned j=0; j<N-k; ++j)
                assert(new_cost[k][j] >= 0);
        for (unsigned k=0; k<max_k; ++k)
            for (unsigned j=N-k; j<N; ++j)
                assert(new_cost[k][j] <= 0);


        // Penalize short segments
        if (conf.min_num_bins_per_segment > 1) {
            for (unsigned k = 0; k + 1 < conf.min_num_bins_per_segment && k < max_k; ++k) {
                for (unsigned j = 0; j < N-k; ++j) {
                    new_cost[k][j] += static_cast<float>(100);
                }
            }
        }



        // Alternative 3.2
        // Force segmentation to use "None" intervals.
        // This is done by setting the cost for such intervals to 0 and adding
        // a penalty to intervals violating the boarders.
        if (vm.count("penalize-none"))
        {
            if (chrom == 0)
                std::cout << "[Warning] --penalize-none is currently buggy - some None stretches are being ignored. This might be due to repeated overwriting of values in the breakpoint table, but I am not sure yet" << std::endl;

            for (auto stretch : none_reg_local[chrom])
            {
                // Penalize all segments that violate these boarders.
                // Note that this penalty becomes weaker the longer the segment!
                for (unsigned k = 0; k < max_k; ++k)
                    for (unsigned j = (k>stretch.first ? 0 : stretch.first - k); j < stretch.second && j < N-k; ++j)
                        new_cost[k][j] += conf.none_penalty;

                // Finally, favor the usage of exactly the consecutive None segment
                new_cost[stretch.second - stretch.first][stretch.first] = 0;
            }
        } // vm.count("penalize-none")

        if (N<1) {
            std::cout << "[Warning]: Chromosome " << chromosomes[chrom] << " was completely removed via None bins" << std::endl;
            continue;
        }


        // print cost matrix for a chromosome
        if (vm.count("cost-matrix") && chromosomes[chrom] == conf.cm_chrom) {
            std::ofstream out(conf.f_cost_mat.string());
            if (out.is_open()) {
                std::cout << "[Write] cost matrix for chromosome '" << conf.cm_chrom << "' to " << conf.f_cost_mat.string() << std::endl;
                out << new_cost;
            } else {
                std::cerr << "[Warning] Cannot write to " << conf.f_cost_mat.string() << std::endl;
            }
        }


        // Find optimal segmentation
        Matrix<int> breakpoints;
        std::vector<double> sse;
        if (!optimal_segment_dp(new_cost, max_cp, breakpoints, sse)) {
            std::cerr << "[Warning] Segmentation failed on " << chromosomes[chrom] << std::endl;
            continue;
        }


        // number of none bins (for output)
        unsigned num_none = 0;
        for (auto stretch : none_reg_local[chrom])
            num_none += stretch.second - stretch.first + 1;
        unsigned num_none_regions = none_reg_local[chrom].size();

        // Output of breakpoints;
        if (out.is_open()) {
            std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
            auto time2 = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3).count();
            std::cout << "[Write] chromosome '" << chromosomes[chrom] << "' (took " << time2 << " sec) to file " << conf.f_out.string() << std::endl;
            for (unsigned cp = 0; cp < max_cp; ++cp) {
                for (unsigned k = 0; k <= cp; ++k) {
                    out << *samples.begin() <<  "\t";
                    out << counts.size() << "\t";
                    out << chromosomes[chrom] << "\t";
                    out << chrom_map[chrom+1] - chrom_map[chrom] << "\t";
                    out << max_cp << "\t";
                    out << max_k << "\t";
                    out << num_none << "\t";
                    out << num_none_regions << "\t";
                    if (vm.count("penalize-none"))
                        out << "penalize=" << conf.none_penalty << "\t";
                    else if (vm.count("remove-none"))
                        out << "remove-none" << "\t";
                    else
                        out << "utilize-none" << "\t";
                    out << cp+1 << "\t";
                    out << sse[cp] << "\t";
                    unsigned from = k==0 ? 0 : breakpoints[cp][k-1];
                    unsigned to   = breakpoints[cp][k]-1;
                    out << to << "\t";
                    out << bins[chrom_map[chrom] + from].start << "\t";
                    out << bins[chrom_map[chrom] + to].end << std::endl;
                }
            }
        } else {
            std::cerr << "[Warning] Cannot write to " << conf.f_out.string() << std::endl;
        }

    } // for chrom

    return 0;
} // main



#endif /* segmentation_hpp */

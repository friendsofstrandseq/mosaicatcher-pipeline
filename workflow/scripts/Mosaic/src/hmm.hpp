/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#ifndef hmm_hpp
#define hmm_hpp

#include <iostream>
#include <vector>
#include <cassert>
#include "utils.hpp"

// todo: remove if std::setprecision no longer needed
#include <iomanip>

namespace hmm {

    using count::TGenomeCounts;

    static const double INPUT_PRECISION = 0.0001;
    static const double MIN_HMM_LOG = -1e8;

    // Note: Do not copy this class (no copy constructor defined)
    template <typename TEmission, typename TDistribution>
    class HMM
    {
        /** 
         *  N = number of states
         */
    public:

        explicit HMM(std::vector<std::string> const & names) :
            N(names.size()), state_labels(names)
        {
            assert(N>0);
            assert(names.size() == N);

            pi = new double[N];
            log_trans = new double*[N];
            for (uint8_t i=0; i<N; ++i)
                log_trans[i] = new double[N];
        }

        ~HMM()
        {
            delete[] pi;
            for (uint8_t i=0; i<N; ++i)
                delete[] log_trans[i];
            delete[] log_trans;
        }

        void set_state_names(std::vector<std::string> const & names) {
            assert(names.size() == N);
            state_labels = names;
        }

        void set_transitions(std::vector<double> const & t)
        {
            assert(t.size() == N*N);
            for (uint8_t i=0; i<N; ++i) {
                double sum = 0;
                for (uint8_t j=0; j<N; ++j) {
                    sum += t[N*i+j];
                    log_trans[i][j] = log(t[N*i+j]);
                }
                assert(std::abs(1 - sum) < INPUT_PRECISION);
            }
        }

        void set_initials(std::vector<double> const & ini) {
            assert(ini.size() == N);
            assert(fabs(1 - sum(ini)) < INPUT_PRECISION);
            for (uint8_t i=0; i<N; ++i)
                pi[i] = ini[i];
        }

        void set_emissions(typename std::vector<TDistribution> const & e)
        {
            assert(e.size() == N);
            // assert that all distributions have the same dimension
            int dim = e[0].dim;
            for (auto x : e) {
                assert(x.dim == dim);
            }
            distributions = e; // copy!
        }

        std::vector<uint8_t> const & get_path() const
        {
            return path;
        }

        std::vector<std::string> get_path_labels() const {
            std::vector<std::string> labels;
            for (unsigned i=0; i<path.size(); i++) {
                assert(path[i] < N);
                labels.push_back(state_labels[path[i]]);
            }
            return labels;
        }


        // Algorithms
        double viterbi(std::vector<TEmission> const &); // writes "path"


//private:

        std::vector<std::string> state_labels;
        uint8_t N;           // Number of states
        double* pi;          // initial transition probabilites
        double** log_trans;  // transition probabilites (log scale)

        std::vector<double> calc_log_emissions(typename std::vector<TEmission>::const_iterator iter);
        typename std::vector<TDistribution> distributions;   // emission parameters - one "FakeMultivariateGaussian" per state
        std::vector<uint8_t> path;      // last path calculated;
    };


    /**
     *  Given an observation (which consists of multiple doubles for a 
     *  multivariate model), calculate the PDFs for all states.
     *  Result should be move-constructed efficiently
     *  TODO: how can I whether it is really moved and not copied?
     */
    template <typename TEmission, typename TDistribution>
    std::vector<double> HMM<TEmission, TDistribution>::calc_log_emissions(typename std::vector<TEmission>::const_iterator iter)
    {
        typename std::vector<TEmission>::const_iterator start = iter;
        std::vector<double> log_probs(N);

        for (uint8_t state = 0; state < N; ++state)
        {
            // multivariate emissions are encoded succeedingly in one long vector,
            // so "dist" consumed dist.dim many emissions.
            TDistribution const& dist = distributions[state];
            log_probs[state] = dist.calc_log_emission(iter);
            iter = start;
        }
        return log_probs;
    }




    /**
     *  Viterbi algorithm
     *  -----------------
     *  Sequence must be flattened for multivariate emissions, i.e. a 2-dim
     *  obervation across T=4 time points, e.g. [(0,1), (2,4), (1,3), (4,5)],
     *  would be written [0,1,2,4,1,3,4,5].
     */
    template <typename TEmission, typename TDistribution>
    double HMM<TEmission, TDistribution>::viterbi(std::vector<TEmission> const & seq)
    {
        // length of the sequence
        unsigned dim = distributions[0].dim;
        unsigned T = (unsigned)seq.size()/dim;

        // Nothing to do for empty sequences
        if (seq.size() < dim) {
            path = std::vector<uint8_t>();
            return 0;
        }


        // allocate a matrix for probabilities:
        // matrix[T=columns][N=rows]
        std::vector<std::vector<double> > matrix(T, std::vector<double>(N, 1));

        // allocate a matrix for the path through the matrix
        std::vector<std::vector<uint8_t> > trace(T, std::vector<uint8_t>(N, N));

        // Step 1: Initiation
        // Prepare emissions for first observation:
        typename std::vector<TEmission>::const_iterator it = seq.begin();
        std::vector<double> log_emission = calc_log_emissions(it);
        it += dim; // First observation was already processed

        for (uint8_t state = 0; state < N; ++state)
            matrix[0][state] = log(pi[state]) + log_emission[state];


        // Step 2: Iteration
        for (unsigned i = 1; i < T; ++i, it+=dim)
        {
            // get emission for observation(s) i
            log_emission = calc_log_emissions(it);


            for (uint8_t state = 0; state < N; ++state)
            {
                double max_log = MIN_HMM_LOG; // todo: must actually be smaller than all log(p)
                uint8_t max_path = N;

                for (uint8_t prev_state = 0; prev_state < N; ++prev_state) {
                    double log_prob = log_trans[prev_state][state] + matrix[i-1][prev_state];
                    if (log_prob > max_log) {
                        max_log = log_prob;
                        max_path = prev_state;
                    }
                }
                matrix[i][state] = max_log + log_emission[state];
                trace[i][state] = max_path;
            }
        }


        // Get maximum in last column
        double max_log = MIN_HMM_LOG;
        uint8_t max_path = N;
        for (uint8_t state = 0; state < N; ++state) {
            if (matrix[T-1][state] >= max_log) {
                max_log = matrix[T-1][state];
                max_path = state;
            }
        }


        // Step 3: Traceback
        // overwrite the member "path"
        path = std::vector<uint8_t>(T, N);
        for (unsigned i=T; i>0; ) {
            --i;
            path[i] = max_path;
            max_path = trace[i][max_path];
        }


        // print
        if (false) {
                        it = seq.begin();
                        log_emission = calc_log_emissions(it);
                        it += dim; // First observation was already processed

                        // header
                        std::cout << std::endl;
                        std::cout << std::setw(6) << "bin";
                        std::cout << std::setw(10*(int)N) << "matrix";
                        std::cout << std::setw(6*(int)N)  << "path";
                        std::cout << std::setw(10*(int)N) << "emission";
                        std::cout << std::setw(6) << "label";
                        std::cout << std::endl;

                        std::cout << std::setw(6) << "-----";
                        std::cout << " "; for (unsigned u=1; u<10*(unsigned)N; ++u) std::cout << "-";
                        std::cout << " "; for (unsigned u=1; u<6*(unsigned)N; ++u) std::cout << "-";
                        std::cout << " "; for (unsigned u=1; u<10*(unsigned)N; ++u) std::cout << "-";
                        std::cout << std::setw(6) << "-----";
                        std::cout << std::endl;

                        // matrix and path
                        for (unsigned i = 0; i < T; ++i)
                        {
                            log_emission = calc_log_emissions(it);
                            it += dim;

                            std::cout << std::setw(6) << (int)i;

                            // first set of cols: matrix
                            for (uint8_t state = 0; state < N; ++state)
                                std::cout << std::setw(10) << std::setprecision(1) << std::fixed << matrix[i][state];
                            // second set of cols: path
                            for (uint8_t state = 0; state < N; ++state)
                                std::cout << std::setw(6) << (int)trace[i][state];
                            // third set of cols: emisison probs
                            for (uint8_t state = 0; state < N; ++state)
                                std::cout << std::setw(10) << std::setprecision(1) << std::fixed << log_emission[state];
                            // last col: label
                            if (path[i] < N)
                                std::cout << std::setw(6) << state_labels[path[i]];
                            std::cout << std::endl;

                        }
        } // end of print

        return max_log;
    }

    //template HMM<double, MultiVariate<Gaussian> >;









    template <typename THMM>
    void run_HMM(THMM & hmm, TGenomeCounts & counts, std::vector<unsigned> const & good_bins, std::vector<int32_t> const & chrom_map)
    {
        for (int32_t chrom = 0; chrom < chrom_map.size()-1; ++chrom) {

            // skip empty chromosomes
            if (chrom_map[chrom+1] - chrom_map[chrom] < 1)
                continue;

            // Order: crick, watson, crick, watson, ...
            // but only use good_bins !!
            std::vector<unsigned> seq;
            for (unsigned bin = chrom_map[chrom]; bin < chrom_map[chrom+1]; ++bin) {
                seq.push_back(counts[good_bins[bin]].crick_count);
                seq.push_back(counts[good_bins[bin]].watson_count);
            }

            // Run viterbi
            hmm.viterbi(seq);
            std::vector<std::string> path = hmm.get_path_labels();
            assert(path.size() == (chrom_map[chrom+1] - chrom_map[chrom]));

            // write classification into Counter
            unsigned bin_in_path = 0;
            for (unsigned bin = chrom_map[chrom]; bin < chrom_map[chrom+1]; ++bin)
                counts[good_bins[bin]].set_label(path[bin_in_path++]);
        }
    }






}
#endif /* hmm_hpp */

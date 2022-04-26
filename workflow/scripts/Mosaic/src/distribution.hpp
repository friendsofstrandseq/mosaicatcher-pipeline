/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#ifndef distribution_h
#define distribution_h

#include <iostream>
#include <vector>
#include <cassert>
#include "utils.hpp"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/binomial.hpp>


namespace hmm {


    static const double MIN_EMISSION_LOG = -1e4;
    // TODO: check for p=0 or 1 in all pdfs --> take special care for log then!


    /**
     *  NegativeBinomial
     *  ----------------
     *  p(k) = {\Gamma(n + k) \over \Gamma(k+1) \Gamma(n) } p^n (1-p)^k
     *  Number of failures k before n-th success.
     *  Mean :=  n * (1-p)/p
     *  Variance := Mean/p
     */
    class NegativeBinomial
    {
    public:
        typedef unsigned TEmission;
        static const uint8_t dim = 1;
        boost::math::negative_binomial nb;


        NegativeBinomial(double p, double n) : nb(n,p)
        {}

        inline double calc_log_emission(std::vector<unsigned>::const_iterator iter) const
        {
            double pr = log(calc_emission(iter));
            if (pr < MIN_EMISSION_LOG)
                pr = MIN_EMISSION_LOG;
            return pr;
        }

        inline double calc_emission(std::vector<unsigned>::const_iterator iter) const
        {
            return boost::math::pdf(nb, *iter);
        }
    };

    std::ostream& operator<<(std::ostream& os, const NegativeBinomial& obj)
    {
        os << "Negative Binomial Distribution:" << std::endl;
        os << "       p = " << obj.nb.success_fraction() << std::endl;
        os << "       r = " << obj.nb.successes() << std::endl;
        os << "    mean = " << mean(obj.nb) << std::endl;
        os << "     var = " << variance(obj.nb) << std::endl;
        return os;
    }



    /**
     *  MultiVariate
     *  ------------
     *  A quick wrapper to sum up multiple 1-dim distirbutions of the same type into
     *  one "multivariate" one, in which components are independent.
     */
    template <typename TDistribution>
    class MultiVariate
    {
    public:
        uint8_t dim;
        std::vector<TDistribution> inner;
        double log_prior;

        MultiVariate(std::vector<TDistribution> const & distributions, double _prior = 1) :
        inner(distributions), dim(distributions.size()), log_prior(log(_prior))
        {
            assert(distributions.size() > 0 && distributions.size() < 256);
            for (auto dist : distributions)
                assert(dist.dim == (uint8_t)1);
        }

        inline double calc_log_emission(typename std::vector<typename TDistribution::TEmission>::const_iterator iter) const
        {
            double logp = 0;
            for (unsigned i=0; i<dim; ++i, ++iter)
                logp += inner[i].calc_log_emission(iter);
            logp += log_prior;
            if (logp < MIN_EMISSION_LOG)
                logp = MIN_EMISSION_LOG;
            return logp;
        }

        inline double calc_emission(std::vector<unsigned>::const_iterator iter) const
        {
            return exp(calc_log_emission(iter));
        }
    };

    template <typename TDistribution>
    std::ostream& operator<<(std::ostream& os, const MultiVariate<TDistribution>& obj)
    {
        os << "Multivariate distribution with " << (int)obj.dim << " dimensions:" << std::endl;
        for (unsigned i=0; i<obj.inner.size(); ++i)
            os << ">> " <<  i << ") " <<  obj.inner[i];
        os << std::endl;
        return os;
    }




    /**
     *  CombinedNegBinAndBinomial
     *  ----------------
     *  Model total coverage (w+c) by Negative Binomial and fraction c/(w+c) by
     *  a binomial distribution.
     */
    class CombinedNegBinAndBinomial
    {
    public:
        typedef unsigned TEmission;
        static const uint8_t dim = 2;
        boost::math::negative_binomial nb;
        double ratio;
        double prior;

        CombinedNegBinAndBinomial(double nb_p, double nb_r, double strand_ratio, double prior = 1) :
        nb(nb_r,nb_p), ratio(strand_ratio), prior(prior)
        {}

        inline double calc_log_emission(std::vector<unsigned>::const_iterator iter) const
        {
            double pr = log(calc_emission(iter));
            if (pr < MIN_EMISSION_LOG)
                pr = MIN_EMISSION_LOG;
            return pr;
        }

        inline double calc_emission(std::vector<unsigned>::const_iterator iter) const
        {
            unsigned c = *iter++;
            unsigned w = *iter++;
            boost::math::binomial_distribution<> binomial(c+w, ratio);
            return prior * boost::math::pdf(nb, w+c) * boost::math::pdf(binomial, w);
        }
    };


    
}
#endif /* distribution_h */

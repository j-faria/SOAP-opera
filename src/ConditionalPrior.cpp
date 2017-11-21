#include "ConditionalPrior.h"
#include "DNest4.h"
#include "Utils.h"
#include <cmath>
#include <typeinfo>

using namespace std;
using namespace DNest4;


extern ContinuousDistribution *log_muP_prior;
extern ContinuousDistribution *wP_prior;
extern ContinuousDistribution *log_muK_prior;

extern ContinuousDistribution *Pprior;
extern ContinuousDistribution *Kprior;
extern ContinuousDistribution *eprior;
extern ContinuousDistribution *phiprior;
extern ContinuousDistribution *wprior;


PlanetConditionalPrior::PlanetConditionalPrior()
{}


void PlanetConditionalPrior::from_prior(RNG& rng)
{
    if(hyperpriors)
    {
        center = log_muP_prior->rvs(rng);
        width = wP_prior->rvs(rng);
        muK = exp(log_muK_prior->rvs(rng));
    }
}

double PlanetConditionalPrior::perturb_hyperparameters(RNG& rng)
{
    double logH = 0.;

    if(hyperpriors)
    {
        int which = rng.rand_int(3);

        if(which == 0)
        {
            logH -= log_muP_prior->log_pdf(center);
            center = log_muP_prior->rvs(rng);
            logH += log_muP_prior->log_pdf(center);
        }
        else if(which == 1)
            width = wP_prior->rvs(rng);
        else
        {
            muK = log(muK);

            logH -= log_muK_prior->log_pdf(muK);
            muK = log_muK_prior->rvs(rng);
            logH += log_muK_prior->log_pdf(muK);

            muK = exp(muK);
        }
    }

    return logH;
}

// vec[0] = period
// vec[1] = amplitude
// vec[2] = phase
// vec[3] = ecc
// vec[4] = viewing angle

double PlanetConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
    //cout << "type of Pprior:" << typeid(Pprior).name() << endl;
    if(hyperpriors)
    {
        if(vec[2] < 0. || vec[2] > 2.*M_PI ||
           vec[3] < 0. || vec[3] >= 1.0 ||
           vec[4] < 0. || vec[4] > 2.*M_PI)
             return -1E300;

        Pprior = new Laplace(center, width);
        Kprior = new Exponential(muK);
    }
    else
    {
        if(vec[0] < 1. || vec[0] > 1E4 ||
           vec[1] < 0. ||
           vec[2] < 0. || vec[2] > 2.*M_PI ||
           vec[3] < 0. || vec[3] >= 1.0 ||
           vec[4] < 0. || vec[4] > 2.*M_PI)
             return -1E300;
    }

    return Pprior->log_pdf(vec[0]) + 
           Kprior->log_pdf(vec[1]) + 
           phiprior->log_pdf(vec[2]) + 
           eprior->log_pdf(vec[3]) + 
           wprior->log_pdf(vec[4]);
}

void PlanetConditionalPrior::from_uniform(std::vector<double>& vec, int id) const
{
    if(hyperpriors)
    {
        Pprior = new Laplace(center, width);
        Kprior = new Exponential(muK);
    }
    vec[0] = Pprior->cdf_inverse(vec[0]);
    vec[1] = Kprior->cdf_inverse(vec[1]);
    vec[2] = phiprior->cdf_inverse(vec[2]); //2.*M_PI*vec[2];
    vec[3] = eprior->cdf_inverse(vec[3]);
    vec[4] = wprior->cdf_inverse(vec[4]); //2.*M_PI*vec[4];
}

void PlanetConditionalPrior::to_uniform(std::vector<double>& vec, int id) const
{
    if(hyperpriors)
    {
        Pprior = new Laplace(center, width);
        Kprior = new Exponential(muK);
    }
    vec[0] = Pprior->cdf(vec[0]);
    vec[1] = Kprior->cdf(vec[1]);
    vec[2] = phiprior->cdf(vec[2]); //vec[2]/(2.*M_PI);
    vec[3] = eprior->cdf(vec[3]);
    vec[4] = wprior->cdf(vec[4]); //vec[4]/(2.*M_PI);
}

void PlanetConditionalPrior::print(std::ostream& out) const
{
    if(hyperpriors)
        out<<center<<' '<<width<<' '<<muK<<' ';
}


void PlanetConditionalPrior::print0(std::ostream& out) const
{
    out<<0.<<' '<<0.<<' '<<0.<<' ';
}





// conditional prior for the active regions

extern ContinuousDistribution *size_prior;
extern ContinuousDistribution *longitude_prior;
extern ContinuousDistribution *latitude_prior;

SOAPConditionalPrior::SOAPConditionalPrior()
{}


void SOAPConditionalPrior::from_prior(RNG& rng)
{}

double SOAPConditionalPrior::perturb_hyperparameters(RNG& rng)
{
    double logH = 0.;
    return logH;
}

// vec[0] = size
// vec[1] = longitude
// vec[2] = latitude

double SOAPConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
    if(vec[0] < 0. || vec[0] > 1. ||
       vec[1] < 0. || vec[1] > 2.*M_PI ||
       vec[2] < 0. || vec[2] > 2.*M_PI)
         return -1E300;

    return 0.;
    // return size_prior->log_pdf(vec[0]) + 
    //        longitude_prior->log_pdf(vec[1]) + 
    //        latitude_prior->log_pdf(vec[2]);
}

void SOAPConditionalPrior::from_uniform(std::vector<double>& vec, int id) const
{
    vec[0] = size_prior->cdf_inverse(vec[0]);
    vec[1] = longitude_prior->cdf_inverse(vec[1]);
    vec[2] = latitude_prior->cdf_inverse(vec[2]); //2.*M_PI*vec[2];
}

void SOAPConditionalPrior::to_uniform(std::vector<double>& vec, int id) const
{
    vec[0] = size_prior->cdf(vec[0]);
    vec[1] = longitude_prior->cdf(vec[1]);
    vec[2] = phiprior->cdf(vec[2]); //vec[2]/(2.*M_PI);
}

void SOAPConditionalPrior::print(std::ostream& out) const
{}


void SOAPConditionalPrior::print0(std::ostream& out) const
{}





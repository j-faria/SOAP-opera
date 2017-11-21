#ifndef DNest4_Model
#define DNest4_Model

#include <vector>
#include "ConditionalPrior.h"
#include "RJObject/RJObject.h"
#include "RNG.h"
#include "Data.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

// whether the model includes a GP component
extern const bool GP;

// whether there are observations after the change in HARPS fibers
extern const bool obs_after_HARPS_fibers;

// whether the model includes a linear trend
extern const bool trend;

#define RADSUN 696000;
#define lightC 299792458.  // speed of light in m/s

class Star
{
    public:
        double prot, vrot, incl, limba1, limba2, psi, \
               rad, Temp, Temp_diff_spot;
};

class CCF
{
    public:
        std::vector<double> rv, intensity;
        double width, step, v_interval;
        int n, n_v;
};

// class ActiveRegion
// {
//     public:
//         double longi, lati, size;
//         int active_region_type; // 0 spot, 1 plage
//         int check; // 0 false, 1 true
// };


class Model
{
    private:
        DNest4::RJObject<PlanetConditionalPrior> planets;
        DNest4::RJObject<SOAPConditionalPrior> active_regions;

        double background;
        //std::vector<double> offsets;
        double slope, quad;
        double fiber_offset;

        double extra_sigma;

        // Parameters for the quasi-periodic extra noise
        double eta1, eta2, eta3, eta4, eta5;
        double a,b,c,P;

        /*celerite::solver::CholeskySolver<double> solver;
        Eigen::VectorXd alpha_real,
                 beta_real,
                 alpha_complex_real,
                 alpha_complex_imag,
                 beta_complex_real,
                 beta_complex_imag;*/

        // The signal
        std::vector<long double> mu;
        void calculate_mu();




        // eccentric and true anomalies
        double ecc_anomaly(double time, double prd, double ecc, double peri_pass);
        double eps3(double e, double M, double x);
        double keplerstart3(double e, double M);
        double true_anomaly(double time, double prd, double ecc, double peri_pass);

        // The covariance matrix for the data
        Eigen::MatrixXd C;
        void calculate_C();

        //QPkernel *kernel;
        //HODLR_Tree<QPkernel> *A;

        unsigned int staleness;

    public:
        Model();

        void setupHODLR();

        // The star, ccf and active regions
        Star star;
        CCF ccf, ccf_active_region;

        // struct star
        // {
        //   double prot, vrot, incl, limba1, limba2, psi, \
        //          rad_sun, rad, Temp, Temp_diff_spot;
        // };

        // Generate the point from the prior
        void from_prior(DNest4::RNG& rng);

        // Metropolis-Hastings proposals
        double perturb(DNest4::RNG& rng);

        // Likelihood function
        double log_likelihood() const;

        // Print to stream
        void print(std::ostream& out) const;

        // Return string with column information
        std::string description() const;
};

#endif


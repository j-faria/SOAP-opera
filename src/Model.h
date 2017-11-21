#ifndef DNest4_Model
#define DNest4_Model

#include <math.h>
#include <vector>
#include <algorithm> // for std::copy
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

#define RADSUN 696000
#define lightC 299792458.  // speed of light in m/s


class Star
{
    public:
        //Star(double prot=25.05);
        double prot=25.05;
        //double vrot=(2.0*M_PI * 1.0 * RADSUN)/(prot*86400.);
        double incl=90.;
        double limba1=0.29, limba2=0.34;
        double psi=0.0;
        double rad=1.0;
        double Temp=5778, Temp_diff_spot=663;

        double vrot() const {return (2.0*M_PI * rad * RADSUN)/(prot * 86400.); }
};

class CCF
{
    public:
        double rv[401] = {-20.,-19.9,-19.8,-19.7,-19.6,-19.5,-19.4,-19.3,-19.2,-19.1,-19.,-18.9,-18.8,-18.7,-18.6,-18.5,-18.4,-18.3,-18.2,-18.1,-18.,-17.9,-17.8,-17.7,-17.6,-17.5,-17.4,-17.3,-17.2,-17.1,-17.,-16.9,-16.8,-16.7,-16.6,-16.5,-16.4,-16.3,-16.2,-16.1,-16.,-15.9,-15.8,-15.7,-15.6,-15.5,-15.4,-15.3,-15.2,-15.1,-15.,-14.9,-14.8,-14.7,-14.6,-14.5,-14.4,-14.3,-14.2,-14.1,-14.,-13.9,-13.8,-13.7,-13.6,-13.5,-13.4,-13.3,-13.2,-13.1,-13.,-12.9,-12.8,-12.7,-12.6,-12.5,-12.4,-12.3,-12.2,-12.1,-12.,-11.9,-11.8,-11.7,-11.6,-11.5,-11.4,-11.3,-11.2,-11.1,-11.,-10.9,-10.8,-10.7,-10.6,-10.5,-10.4,-10.3,-10.2,-10.1,-10.,-9.9,-9.8,-9.7,-9.6,-9.5,-9.4,-9.3,-9.2,-9.1,-9.,-8.9,-8.8,-8.7,-8.6,-8.5,-8.4,-8.3,-8.2,-8.1,-8.,-7.9,-7.8,-7.7,-7.6,-7.5,-7.4,-7.3,-7.2,-7.1,-7.,-6.9,-6.8,-6.7,-6.6,-6.5,-6.4,-6.3,-6.2,-6.1,-6.,-5.9,-5.8,-5.7,-5.6,-5.5,-5.4,-5.3,-5.2,-5.1,-5.,-4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9.,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12.,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13.,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14.,14.1,14.2,14.3,14.4,14.5,14.6,14.7,14.8,14.9,15.,15.1,15.2,15.3,15.4,15.5,15.6,15.7,15.8,15.9,16.,16.1,16.2,16.3,16.4,16.5,16.6,16.7,16.8,16.9,17.,17.1,17.2,17.3,17.4,17.5,17.6,17.7,17.8,17.9,18.,18.1,18.2,18.3,18.4,18.5,18.6,18.7,18.8,18.9,19.,19.1,19.2,19.3,19.4,19.5,19.6,19.7,19.8,19.9,20.};
        double intensity[401] = {1.,1.,1.,1.,1.,1.,0.999997,0.999987,0.999969,0.999946,0.999914,0.999873,0.999824,0.999770,0.999712,0.999646,0.999565,0.999467,0.999363,0.999256,0.999144,0.999007,0.998853,0.998707,0.998587,0.998473,0.998363,0.998212,0.998016,0.997813,0.997606,0.997400,0.997183,0.996980,0.996792,0.996615,0.996475,0.996329,0.996146,0.995939,0.995630,0.995291,0.994955,0.994696,0.994499,0.994347,0.994125,0.993960,0.993763,0.993576,0.993308,0.992977,0.992587,0.992159,0.991826,0.991566,0.991394,0.991203,0.991100,0.990897,0.990719,0.990564,0.990412,0.990156,0.989887,0.989589,0.989304,0.989158,0.989121,0.989000,0.988843,0.988623,0.988432,0.988528,0.988838,0.988901,0.988589,0.988242,0.988125,0.988319,0.988456,0.988268,0.987867,0.987870,0.988523,0.989242,0.989498,0.989436,0.989436,0.989718,0.990232,0.990369,0.990243,0.990289,0.990740,0.991413,0.992010,0.992503,0.992836,0.993384,0.993906,0.994086,0.993972,0.994138,0.994707,0.995474,0.996239,0.996464,0.996451,0.996664,0.997101,0.997342,0.997075,0.996495,0.995906,0.995707,0.995824,0.995999,0.995926,0.995851,0.995530,0.994877,0.994075,0.993292,0.992458,0.991459,0.990311,0.988827,0.987561,0.986408,0.985270,0.983787,0.982087,0.980110,0.977944,0.975737,0.973192,0.970506,0.967436,0.964180,0.960548,0.956900,0.953096,0.949075,0.944613,0.939960,0.935270,0.930378,0.925143,0.918987,0.912432,0.905501,0.898532,0.890977,0.882767,0.874185,0.865590,0.856954,0.848009,0.838299,0.828086,0.817697,0.807004,0.795838,0.783927,0.771817,0.759510,0.747148,0.734578,0.721594,0.708664,0.695641,0.682500,0.669051,0.655480,0.642190,0.629008,0.616012,0.602902,0.589843,0.577070,0.564816,0.553061,0.541524,0.530339,0.519612,0.509466,0.499994,0.491271,0.482997,0.475209,0.467924,0.461424,0.455704,0.450760,0.446442,0.442709,0.439744,0.437645,0.436381,0.435681,0.435726,0.436643,0.438486,0.441322,0.444862,0.449187,0.454290,0.460341,0.467379,0.475094,0.483560,0.492589,0.502501,0.513127,0.524429,0.536296,0.548498,0.561088,0.574075,0.587437,0.600949,0.614518,0.628263,0.642166,0.656286,0.670336,0.684196,0.697933,0.711569,0.725019,0.737943,0.750404,0.762499,0.774549,0.786308,0.797660,0.808542,0.819137,0.829603,0.839667,0.849133,0.858093,0.866652,0.875001,0.882857,0.890268,0.897163,0.904086,0.910890,0.917176,0.923064,0.928494,0.933710,0.938797,0.943597,0.948045,0.952163,0.956077,0.959752,0.963169,0.966586,0.969730,0.972528,0.974963,0.977306,0.979622,0.981972,0.984063,0.985695,0.987163,0.988779,0.990388,0.991791,0.992825,0.993676,0.994643,0.995759,0.996635,0.997173,0.997624,0.998121,0.998558,0.998756,0.998829,0.998932,0.999114,0.999471,0.999726,0.999886,0.999941,1.,0.999905,0.999711,0.999497,0.999235,0.998965,0.998792,0.998652,0.998623,0.998402,0.998231,0.997979,0.997788,0.997673,0.997189,0.996631,0.996106,0.995952,0.995829,0.995604,0.995154,0.994628,0.994199,0.993973,0.993790,0.993558,0.993261,0.992942,0.992610,0.992343,0.992156,0.992107,0.991981,0.991715,0.991260,0.990982,0.990884,0.990872,0.990873,0.990625,0.990240,0.990035,0.990181,0.990441,0.990637,0.990568,0.990440,0.990410,0.990633,0.990905,0.991024,0.990921,0.990893,0.991059,0.991364,0.991663,0.991866,0.992063,0.992306,0.992566,0.992819,0.993001,0.993267,0.993580,0.993909,0.994145,0.994326,0.994517,0.994791,0.995117,0.995375,0.995578,0.995797,0.996063,0.996393,0.996687,0.996930,0.997103,0.997267,0.997463,0.997660,0.997855,0.998037,0.998206,0.998384,0.998576,0.998750,0.998910,0.999045,0.999160,0.999261,0.999357,0.999457,0.999547,0.999625,0.999694,0.999759,0.999818,0.999871,0.999912,0.999944,0.999968,0.999986,0.999997,1.,1.,1.,1.,1.,1.};
        
        double width = *std::max_element(std::begin(rv), std::end(rv));
        double step = rv[1] - rv[0];
        int n = sizeof rv / sizeof *rv;;
        int n_v = 445;  // need to set this right, but it needs star.vrot!
        double v_interval = step * (n_v - 1) / 2.;
};

class CCFar
{
    public:
        double rv[401] = {-20.,-19.9,-19.8,-19.7,-19.6,-19.5,-19.4,-19.3,-19.2,-19.1,-19.,-18.9,-18.8,-18.7,-18.6,-18.5,-18.4,-18.3,-18.2,-18.1,-18.,-17.9,-17.8,-17.7,-17.6,-17.5,-17.4,-17.3,-17.2,-17.1,-17.,-16.9,-16.8,-16.7,-16.6,-16.5,-16.4,-16.3,-16.2,-16.1,-16.,-15.9,-15.8,-15.7,-15.6,-15.5,-15.4,-15.3,-15.2,-15.1,-15.,-14.9,-14.8,-14.7,-14.6,-14.5,-14.4,-14.3,-14.2,-14.1,-14.,-13.9,-13.8,-13.7,-13.6,-13.5,-13.4,-13.3,-13.2,-13.1,-13.,-12.9,-12.8,-12.7,-12.6,-12.5,-12.4,-12.3,-12.2,-12.1,-12.,-11.9,-11.8,-11.7,-11.6,-11.5,-11.4,-11.3,-11.2,-11.1,-11.,-10.9,-10.8,-10.7,-10.6,-10.5,-10.4,-10.3,-10.2,-10.1,-10.,-9.9,-9.8,-9.7,-9.6,-9.5,-9.4,-9.3,-9.2,-9.1,-9.,-8.9,-8.8,-8.7,-8.6,-8.5,-8.4,-8.3,-8.2,-8.1,-8.,-7.9,-7.8,-7.7,-7.6,-7.5,-7.4,-7.3,-7.2,-7.1,-7.,-6.9,-6.8,-6.7,-6.6,-6.5,-6.4,-6.3,-6.2,-6.1,-6.,-5.9,-5.8,-5.7,-5.6,-5.5,-5.4,-5.3,-5.2,-5.1,-5.,-4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9.,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12.,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13.,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14.,14.1,14.2,14.3,14.4,14.5,14.6,14.7,14.8,14.9,15.,15.1,15.2,15.3,15.4,15.5,15.6,15.7,15.8,15.9,16.,16.1,16.2,16.3,16.4,16.5,16.6,16.7,16.8,16.9,17.,17.1,17.2,17.3,17.4,17.5,17.6,17.7,17.8,17.9,18.,18.1,18.2,18.3,18.4,18.5,18.6,18.7,18.8,18.9,19.,19.1,19.2,19.3,19.4,19.5,19.6,19.7,19.8,19.9,20.};
        double intensity[401] = {1.,1.,1.,1.,1.,1.,0.999999,0.999996,0.999990,0.999982,0.999971,0.999954,0.999934,0.999913,0.999894,0.999874,0.999844,0.999801,0.999759,0.999721,0.999684,0.999624,0.999550,0.999488,0.999455,0.999437,0.999425,0.999373,0.999278,0.999178,0.999079,0.998984,0.998877,0.998777,0.998691,0.998619,0.998591,0.998550,0.998477,0.998377,0.998174,0.997945,0.997706,0.997545,0.997446,0.997400,0.997263,0.997172,0.997052,0.996935,0.996744,0.996474,0.996138,0.995757,0.995473,0.995287,0.995149,0.995010,0.994956,0.994792,0.994642,0.994498,0.994368,0.994110,0.993862,0.993560,0.993252,0.993087,0.993002,0.992821,0.992596,0.992283,0.991956,0.991898,0.992039,0.991955,0.991492,0.990974,0.990662,0.990636,0.990575,0.990188,0.989540,0.989252,0.989614,0.990037,0.989985,0.989584,0.989206,0.989124,0.989301,0.989075,0.988546,0.988136,0.988155,0.988406,0.988595,0.988628,0.988478,0.988576,0.988669,0.988459,0.987947,0.987730,0.987950,0.988384,0.988813,0.988703,0.988388,0.988325,0.988521,0.988479,0.987910,0.987071,0.986290,0.985938,0.985856,0.985801,0.985465,0.985129,0.984602,0.983697,0.982585,0.981428,0.980206,0.978852,0.977351,0.975515,0.973806,0.972136,0.970473,0.968479,0.966286,0.963787,0.961060,0.958353,0.955347,0.952252,0.948777,0.945100,0.941043,0.936948,0.932717,0.928239,0.923383,0.918382,0.913403,0.908314,0.902904,0.896730,0.890229,0.883414,0.876574,0.869226,0.861407,0.853332,0.845380,0.837479,0.829376,0.820685,0.811687,0.802648,0.793368,0.783703,0.773422,0.763068,0.752666,0.742308,0.731849,0.721075,0.710446,0.699846,0.689229,0.678352,0.667407,0.656800,0.646344,0.636120,0.625814,0.615549,0.605572,0.596071,0.586993,0.578077,0.569391,0.561063,0.553288,0.546101,0.539460,0.533114,0.527089,0.521421,0.516448,0.512077,0.508185,0.504585,0.501405,0.498748,0.496753,0.495352,0.494199,0.493393,0.493055,0.493293,0.494049,0.495182,0.496656,0.498422,0.500611,0.503352,0.506407,0.509871,0.513600,0.517804,0.522440,0.527567,0.533068,0.538805,0.544879,0.551299,0.558172,0.565314,0.572710,0.580400,0.588460,0.597014,0.605839,0.614911,0.624188,0.633690,0.643397,0.653093,0.662787,0.672550,0.682636,0.692805,0.703063,0.713286,0.723627,0.734133,0.744620,0.754884,0.764959,0.774971,0.784955,0.794710,0.804229,0.813423,0.822775,0.832134,0.841126,0.849759,0.857942,0.865919,0.873788,0.881416,0.888662,0.895481,0.902032,0.908339,0.914350,0.920324,0.925910,0.931049,0.935724,0.940241,0.944678,0.949017,0.953028,0.956492,0.959703,0.962990,0.966219,0.969194,0.971701,0.973957,0.976296,0.978722,0.980862,0.982607,0.984208,0.985824,0.987350,0.988602,0.989644,0.990651,0.991709,0.992915,0.993963,0.994856,0.995595,0.996283,0.996800,0.997223,0.997612,0.997910,0.998158,0.998444,0.998759,0.999197,0.999433,0.999640,0.999678,0.999781,1.,0.999886,0.999662,0.999381,0.999441,0.999556,0.999622,0.999416,0.999058,0.998759,0.998674,0.998699,0.998660,0.998502,0.998266,0.998018,0.997848,0.997762,0.997773,0.997658,0.997382,0.996935,0.996679,0.996606,0.996594,0.996566,0.996255,0.995797,0.995509,0.995567,0.995706,0.995756,0.995551,0.995275,0.995119,0.995227,0.995381,0.995362,0.995121,0.994965,0.995012,0.995206,0.995396,0.995472,0.995556,0.995696,0.995870,0.996030,0.996093,0.996233,0.996412,0.996609,0.996714,0.996763,0.996819,0.996971,0.997185,0.997332,0.997423,0.997539,0.997707,0.997943,0.998137,0.998276,0.998343,0.998412,0.998524,0.998640,0.998753,0.998850,0.998947,0.999061,0.999193,0.999309,0.999409,0.999487,0.999550,0.999604,0.999657,0.999715,0.999765,0.999808,0.999845,0.999880,0.999912,0.999941,0.999962,0.999976,0.999986,0.999994,0.999999,1.,1.,1.,1.,1.,1.};
        
        double width = *std::max_element(std::begin(rv), std::end(rv));
        double step = rv[1] - rv[0];
        int n = sizeof rv / sizeof *rv;
        int n_v = 445;  // need to set this right, but it needs star.vrot!
        double v_interval = step * (n_v - 1) / 2.;
};


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

        // The signal
        Eigen::MatrixXd ccfmodel;
        std::vector<long double> mu;
        void calculate_model();

        // eccentric and true anomalies
        double ecc_anomaly(double time, double prd, double ecc, double peri_pass);
        double eps3(double e, double M, double x);
        double keplerstart3(double e, double M);
        double true_anomaly(double time, double prd, double ecc, double peri_pass);

        // The covariance matrix for the data
        Eigen::MatrixXd C;
        void calculate_C();

        unsigned int staleness;

    public:
        Model();
        //~Model();

        int ngrid=300, nrho=20;

        // The star, ccf and active regions
        Star star;
        CCF ccf;
        CCFar ccf_active_region;
        double *CCFstar_quiet, *FLUXstar_quiet;
        double **xyz;
        double **f_spot_flux, **f_spot_bconv;
        double **f_spot_tot, *sum_spot;
        void initialize_star_quiet();


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


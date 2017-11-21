#include "Model.h"
#include "ConditionalPrior.h"
#include "DNest4.h"
#include "RNG.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>
#include <fstream>
#include <chrono>

extern "C" {
#include "starspot.h"
}

using namespace std;
using namespace Eigen;
using namespace DNest4;

#define TIMING false

extern ContinuousDistribution *Cprior; // systematic velocity, m/s
extern ContinuousDistribution *Jprior; // additional white noise, m/s

extern ContinuousDistribution *slope_prior; // m/s/day

extern ContinuousDistribution *log_eta1_prior;
extern ContinuousDistribution *log_eta2_prior;
extern ContinuousDistribution *eta3_prior;
extern ContinuousDistribution *log_eta4_prior;


void Model::from_prior(RNG& rng)
{
    planets.from_prior(rng);
    planets.consolidate_diff();
    
    active_regions.from_prior(rng);
    active_regions.consolidate_diff();

    background = Cprior->rvs(rng);
    extra_sigma = Jprior->rvs(rng);

    if(obs_after_HARPS_fibers)
        // between 0 m/s and 50 m/s
        fiber_offset = 50*rng.rand();

    if(trend)
        slope = slope_prior->rvs(rng);

    if(GP)
    {
        eta1 = exp(log_eta1_prior->rvs(rng)); // m/s

        eta2 = exp(log_eta2_prior->rvs(rng)); // days

        eta3 = eta3_prior->rvs(rng); // days

        eta4 = exp(log_eta4_prior->rvs(rng));
    }

    calculate_model();

    if(GP) calculate_C();

}

void Model::calculate_C()
{

    // Get the data
    const vector<double>& t = Data::get_instance().get_t();
    const vector<double>& sig = Data::get_instance().get_sig();


    int N = Data::get_instance().get_t().size();
    // auto begin = std::chrono::high_resolution_clock::now();  // start timing

    for(size_t i=0; i<N; i++)
    {
        for(size_t j=i; j<N; j++)
        {
            //C(i, j) = eta1*eta1*exp(-0.5*pow((t[i] - t[j])/eta2, 2) );
            C(i, j) = eta1*eta1*exp(-0.5*pow((t[i] - t[j])/eta2, 2) 
                       -2.0*pow(sin(M_PI*(t[i] - t[j])/eta3)/eta4, 2) );

            if(i==j)
                C(i, j) += sig[i]*sig[i] + extra_sigma*extra_sigma; //+ eta5*t[i]*t[i];
            else
                C(j, i) = C(i, j);
        }
    }

    // auto end = std::chrono::high_resolution_clock::now();
    // cout << "old GP: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << " ns" << "\t"; // << std::endl;

}


void Model::initialize_star_quiet()
{
    CCFstar_quiet = new double[ccf.n];
    FLUXstar_quiet = new double[ccf.n];

    xyz = (double **)malloc(sizeof(double *)*nrho);
    for (int j=0; j<nrho; j++) 
        xyz[j] = (double *)malloc(sizeof(double)*3);

    int npsi = 10;
    f_spot_flux = (double **)malloc(sizeof(double *)*npsi);
    f_spot_bconv = (double **)malloc(sizeof(double *)*npsi);
    f_spot_tot = (double **)malloc(sizeof(double *)*npsi);
    for (int j=0; j<npsi; j++) {
        f_spot_flux[j] = (double *)malloc(sizeof(double)*ccf.n_v);
        f_spot_bconv[j] = (double *)malloc(sizeof(double)*ccf.n_v);
        f_spot_tot[j] = (double *)malloc(sizeof(double)*ccf.n_v);
    }

    sum_spot = new double[npsi];

    // Calculates the flux and CCF in each cell of the grid and integrate
    // over the entire stellar disc to have the 
    // integrated flux (FLUXstar) and CCF (CCFstar) for the quiet star
    itot(star.vrot(), star.incl, star.limba1, star.limba2, 0., 0., 0., ngrid,
         ccf.rv, ccf.intensity, ccf.v_interval, ccf.n_v, ccf.n,
         CCFstar_quiet, FLUXstar_quiet);

}

void Model::calculate_model()
{
    // Get the times from the data
    const vector<double>& t = Data::get_instance().get_t();

    // Update or from scratch?
    bool update = (planets.get_added().size() < planets.get_components().size()) &&
            (staleness <= 10);

    // Get the planet components
    const vector< vector<double> >& pcomponents = (update)?(planets.get_added()):
                (planets.get_components());
    // at this point, pcomponents has:
    //  if updating: only the added planets' parameters
    //  if from scratch: all the planets' parameters

    // Zero the signal
    if(!update) // not updating, means recalculate everything
    {
        mu.assign(mu.size(), background);
        staleness = 0;
        if(trend)
        {
            for(size_t i=0; i<t.size(); i++)
                mu[i] += slope*(t[i] - t[0]); //+ quad*(t[i] - t[0])*(t[i] - t[0]);
        }

        if(obs_after_HARPS_fibers)
        {
            for(size_t i=Data::get_instance().index_fibers; i<t.size(); i++)
                //if (i>=Data::get_instance().index_fibers) mu[i] += fiber_offset;
                mu[i] += fiber_offset;
        }


    }
    else // just updating (adding) planets
        staleness++;

    #if TIMING
    auto begin = std::chrono::high_resolution_clock::now();  // start timing
    #endif

    double P, K, phi, ecc, viewing_angle, f, v, ti;
    for(size_t j=0; j<pcomponents.size(); j++)
    {
        if(hyperpriors)
            P = exp(pcomponents[j][0]);
        else
            P = pcomponents[j][0];
        
        K = pcomponents[j][1];
        phi = pcomponents[j][2];
        ecc = pcomponents[j][3];
        viewing_angle = pcomponents[j][4];

        for(size_t i=0; i<t.size(); i++)
        {
            ti = t[i];
            f = true_anomaly(ti, P, ecc, t[0]-(P*phi)/(2.*M_PI));
            v = K*(cos(f+viewing_angle) + ecc*cos(viewing_angle));
            mu[i] += v;
        }
    }

    /*************************************************************************/
    // Get the active region components
    const vector< vector<double> >& arcomponents = active_regions.get_components();
    
    int npsi = 10;
    double psi[npsi];
    // psi[0] = 0.; psi[1] = 0.1;
    for (int i=0; i<npsi; i++)
        psi[i] = 0.1*i;


    double s, lat, lon;

    for(size_t j=0; j<arcomponents.size(); j++)
    {
        s = arcomponents[j][0];
        lon = arcomponents[j][1];
        lat = arcomponents[j][2];
        //cout << s << " " << lon << "  " << lat << endl;

        // Calculates the position of this spot initialized at the disc center
        spot_init(s, lon, lat, star.incl, nrho, xyz);

        // Scans the yz-area where this spot is for different phases (psi) and
        // returns the spot's "non-contribution" to the total flux and its
        // "non-contribution" to the ccf, for each phase.
        spot_scan_npsi(xyz, nrho, psi, npsi, 
                       star.vrot(), star.incl, star.limba1, star.limba2, 
                       0., 0., 0., ngrid,
                       ccf.rv, ccf.intensity, ccf_active_region.intensity,
                       ccf.v_interval, ccf.n_v, ccf.n,
                       s, lon, lat,
                       f_spot_flux, f_spot_bconv, 
                       f_spot_tot, sum_spot,
                       0, star.Temp, star.Temp_diff_spot);
    

        //FLUXstar      = FLUXstar - FLUXactive_region           # total flux of the star affected by active regions
        //CCFstar_flux  = CCFstar_flux - CCFactive_region_flux   # CCF of the star affected by the flux effect of active regions
        //CCFstar_bconv = CCFstar_bconv - CCFactive_region_bconv # CCF of the star affected by the convective blueshift effect of active regions
        //CCFstar_tot   = CCFstar_tot - CCFactive_region_tot     # CCF of the star affected by the total effect of active regions
    }


    #if TIMING
    auto end = std::chrono::high_resolution_clock::now();
    cout << "Model eval took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()*1E-6 << " ms" << std::endl;
    #endif    

}

double Model::perturb(RNG& rng)
{
    const vector<double>& t = Data::get_instance().get_t();
    double logH = 0.;
    double fraction = 0.5;

    logH += active_regions.perturb(rng);
    active_regions.consolidate_diff();
    calculate_model();


    if(rng.rand() <= 0.5)
    {
        logH += planets.perturb(rng);
        planets.consolidate_diff();
        calculate_model();
    }

    if(GP) 
    {
        if(rng.rand() <= 0.25)
        {
            if(rng.rand() <= 0.25)
            {
                eta1 = exp(log_eta1_prior->rvs(rng)); // m/s
                //eta1 = log(eta1);
                //eta1 += log(1E4)*rng.randh(); // range of prior support
                //wrap(eta1, log(1E-5), log(1E-1)); // wrap around inside prior
                //eta1 = exp(eta1);
            }
            else if(rng.rand() <= 0.33330)
            {
                eta2 = exp(log_eta2_prior->rvs(rng)); // days
                //eta2 = log(eta2);
                //eta2 += log(1E12)*rng.randh(); // range of prior support
                //wrap(eta2, log(1E-6), log(1E6)); // wrap around inside prior
                //eta2 = exp(eta2);
            }
            else if(rng.rand() <= 0.5)
            {
                eta3 = eta3_prior->rvs(rng);
                //eta3 += 35.*rng.randh(); // range of prior support
                //wrap(eta3, 15., 50.); // wrap around inside prior
            }
            else
            {
                // eta4 = 1.0;
                eta4 = exp(log_eta4_prior->rvs(rng));
                //eta4 = log(eta4);
                //eta4 += log(1E10)*rng.randh(); // range of prior support
                //wrap(eta4, log(1E-5), log(1E5)); // wrap around inside prior
                //eta4 = exp(eta4);
            }

            calculate_C();
        }
    } // GP


    if(GP) fraction = 0.125;

    if(rng.rand() <= fraction)
    {
        // need to change logH
        logH -= Jprior->log_pdf(extra_sigma);
        extra_sigma = Jprior->rvs(rng);
        logH += Jprior->log_pdf(extra_sigma);

        #if GP
            calculate_C();
        #endif
    }

    if(rng.rand() <= fraction)
    {

        for(size_t i=0; i<mu.size(); i++)
        {
            mu[i] -= background;
            if(trend)
                mu[i] -= slope*(t[i]-t[0]);

            if (obs_after_HARPS_fibers)
                if (i >= Data::get_instance().index_fibers) mu[i] -= fiber_offset;
        }

        background = Cprior->rvs(rng);

        // propose new fiber offset
        if (obs_after_HARPS_fibers) 
        {
            fiber_offset += 50*rng.randh();
            wrap(fiber_offset, 0., 50);
        }

        // propose new slope
        if(trend)
            slope = slope_prior->rvs(rng);

        for(size_t i=0; i<mu.size(); i++)
        {
            mu[i] += background;
            if(slope)
                mu[i] += slope*(t[i]-t[0]);

            if (obs_after_HARPS_fibers)
                if (i >= Data::get_instance().index_fibers) mu[i] += fiber_offset;
        }
    }

    return logH;
}


double Model::log_likelihood() const
{
    int N = Data::get_instance().get_y().size();
    double logL = 0.;
    /** The following code calculates the log likelihood in the case of a GP model */

    // Get the data
    const vector<double>& y = Data::get_instance().get_y();

    #if TIMING
    auto begin = std::chrono::high_resolution_clock::now();  // start timing
    #endif

    if(GP)
    {
        // residual vector (observed y minus model y)
        VectorXd residual(y.size());
        for(size_t i=0; i<y.size(); i++)
            residual(i) = y[i] - mu[i];


        // perform the cholesky decomposition of C
        Eigen::LLT<Eigen::MatrixXd> cholesky = C.llt();
        // get the lower triangular matrix L
        MatrixXd L = cholesky.matrixL();

        double logDeterminant = 0.;
        for(size_t i=0; i<y.size(); i++)
            logDeterminant += 2.*log(L(i,i));

        VectorXd solution = cholesky.solve(residual);

        // y*solution
        double exponent = 0.;
        for(size_t i=0; i<y.size(); i++)
            exponent += residual(i)*solution(i);

        logL = -0.5*y.size()*log(2*M_PI)
               - 0.5*logDeterminant - 0.5*exponent;

        // calculate C^-1*(y-mu)
        // auto begin = std::chrono::high_resolution_clock::now();  // start timing
        // auto end = std::chrono::high_resolution_clock::now();
        // cout << "solve took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << " ns \t";// <<  std::endl;

        // auto begin1 = std::chrono::high_resolution_clock::now();  // start timing
        // auto end1 = std::chrono::high_resolution_clock::now();
        // cout << "solve took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end1-begin1).count() << " ns" << std::endl;

    } 
    else
    {

        /** The following code calculates the log likelihood in the case of a t-Student model without correlated noise*/
        //  for(size_t i=0; i<y.size(); i++)
        //  {
        //      var = sig[i]*sig[i] + extra_sigma*extra_sigma;
        //      logL += gsl_sf_lngamma(0.5*(nu + 1.)) - gsl_sf_lngamma(0.5*nu)
        //          - 0.5*log(M_PI*nu) - 0.5*log(var)
        //          - 0.5*(nu + 1.)*log(1. + pow(y[i] - mu[i], 2)/var/nu);
        //  }

        /** The following code calculates the log likelihood in the case of a Gaussian likelihood*/
        const vector<double>& sig = Data::get_instance().get_sig();

        double halflog2pi = 0.5*log(2.*M_PI);
        double var;
        for(size_t i=0; i<y.size(); i++)
        {
            var = sig[i]*sig[i] + extra_sigma*extra_sigma;
            logL += - halflog2pi - 0.5*log(var)
                    - 0.5*(pow(y[i] - mu[i], 2)/var);
        }

    } //GP

    #if TIMING
    auto end = std::chrono::high_resolution_clock::now();
    cout << "Likelihood took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()*1E-6 << " ms" << std::endl;
    #endif

    if(std::isnan(logL) || std::isinf(logL))
        logL = -1E300;
    return logL;
}

void Model::print(std::ostream& out) const
{
    // output precision
    out.setf(ios::fixed,ios::floatfield);
    out.precision(8);

    out<<extra_sigma<<'\t';
    
    if (obs_after_HARPS_fibers)
        out<<fiber_offset<<'\t';

    if(trend)
        out<<slope<<'\t'; //<<quad*1E8<<'\t';

    if(GP)
        out<<eta1<<'\t'<<eta2<<'\t'<<eta3<<'\t'<<eta4<<'\t';
  
    if (planets.get_fixed() and planets.get_components().size()==0)
        planets.print0(out);
    else
        planets.print(out);

    if (active_regions.get_fixed() and active_regions.get_components().size()==0)
        active_regions.print0(out);
    else
        active_regions.print(out);

    out<<' '<<staleness<<' ';
    out<<background<<' ';
}

string Model::description() const
{
    if(GP)
        return string("extra_sigma   eta1   eta2   eta3   eta4  planets.print   ar.print   staleness   background");
    else
        return string("extra_sigma   planets.print   ar.print   staleness   background");
}



/**
    Calculates the eccentric anomaly at time t by solving Kepler's equation.
    See "A Practical Method for Solving the Kepler Equation", Marc A. Murison, 2006

    @param t the time at which to calculate the eccentric anomaly.
    @param period the orbital period of the planet
    @param ecc the eccentricity of the orbit
    @param t_peri time of periastron passage
    @return eccentric anomaly.
*/
double Model::ecc_anomaly(double t, double period, double ecc, double time_peri)
{
    double tol;
    if (ecc < 0.8) tol = 1e-14;
    else tol = 1e-13;

    double n = 2.*M_PI/period;  // mean motion
    double M = n*(t - time_peri);  // mean anomaly
    double Mnorm = fmod(M, 2.*M_PI);
    double E0 = keplerstart3(ecc, Mnorm);
    double dE = tol + 1;
    double E;
    int count = 0;
    while (dE > tol)
    {
        E = E0 - eps3(ecc, Mnorm, E0);
        dE = abs(E-E0);
        E0 = E;
        count++;
        // failed to converge, this only happens for nearly parabolic orbits
        if (count == 100) break;
    }
    return E;
}


/**
    Provides a starting value to solve Kepler's equation.
    See "A Practical Method for Solving the Kepler Equation", Marc A. Murison, 2006

    @param e the eccentricity of the orbit
    @param M mean anomaly (in radians)
    @return starting value for the eccentric anomaly.
*/
double Model::keplerstart3(double e, double M)
{
    double t34 = e*e;
    double t35 = e*t34;
    double t33 = cos(M);
    return M + (-0.5*t35 + e + (t34 + 1.5*t33*t35)*t33)*sin(M);
}


/**
    An iteration (correction) method to solve Kepler's equation.
    See "A Practical Method for Solving the Kepler Equation", Marc A. Murison, 2006

    @param e the eccentricity of the orbit
    @param M mean anomaly (in radians)
    @param x starting value for the eccentric anomaly
    @return corrected value for the eccentric anomaly
*/
double Model::eps3(double e, double M, double x)
{
    double t1 = cos(x);
    double t2 = -1 + e*t1;
    double t3 = sin(x);
    double t4 = e*t3;
    double t5 = -x + t4 + M;
    double t6 = t5/(0.5*t5*t4/t2+t2);

    return t5/((0.5*t3 - 1/6*t1*t6)*e*t6+t2);
}



/**
    Calculates the true anomaly at time t.
    See Eq. 2.6 of The Exoplanet Handbook, Perryman 2010

    @param t the time at which to calculate the true anomaly.
    @param period the orbital period of the planet
    @param ecc the eccentricity of the orbit
    @param t_peri time of periastron passage
    @return true anomaly.
*/
double Model::true_anomaly(double t, double period, double ecc, double t_peri)
{
    double E = ecc_anomaly(t, period, ecc, t_peri);
    double f = acos( (cos(E)-ecc)/( 1-ecc*cos(E) ) );
    //acos gives the principal values ie [0:PI]
    //when E goes above PI we need another condition
    if(E>M_PI)
      f=2*M_PI-f;

    return f;
}

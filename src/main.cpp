#include <iostream>
#include <typeinfo>
#include "DNest4.h"
#include "Data.h"
#include "Model.h"
#include "ConditionalPrior.h"
#include <Eigen/Core>


using namespace std;
using namespace DNest4;

#include "default_priors.h"

const bool obs_after_HARPS_fibers = false;
const bool GP = false;
const bool hyperpriors = false;
const bool trend = false;


Model::Model()
:planets(5, 0, true, PlanetConditionalPrior())
,active_regions(3, 1, true, SOAPConditionalPrior())
,mu(Data::get_instance().get_t().size())
,C(Data::get_instance().get_t().size(), Data::get_instance().get_t().size())
,star() //erijiriejr
,ccf(),ccf_active_region()
{
    double ymin = Data::get_instance().get_y_min();
    double ymax = Data::get_instance().get_y_max();
    double topslope = Data::get_instance().topslope();

    Cprior = new Uniform(ymin, ymax);
    if(trend)
    	slope_prior = new Uniform(-topslope, topslope);

    star.prot = 25.05;

    Eigen::MatrixXd c = Data::get_instance().get_ccfdata();
    //cout << c.rows() << "x" << c.cols() << endl;

    initialize_star_quiet();

}

int main(int argc, char** argv)
{
	Data::get_instance().load("BL2009_dataset1.kms.rv", "kms");

	Sampler<Model> sampler = setup<Model>(argc, argv);
	//sampler.run();


	return 0;
}

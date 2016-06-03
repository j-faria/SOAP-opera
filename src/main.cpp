#include <iostream>
#include "Lookup.h"
#include "Data.h"
#include "MultiSite2.h"
#include "Start.h"
#include "MyModel.h"

using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{

	DataSet& full = DataSet::getRef("full");
	full.load("data/k10_harpsn.dat", "kms", 0);

	//DataSet& harps = DataSet::getRef("HARPS");
	//DataSet& sophie = DataSet::getRef("SOPHIE");
	//harps.load("test_offsets1.rdb", "ms", 1);
	//sophie.load("test_offsets2.rdb", "ms", 2);

	//cout<<DataSet().nsites()<<endl;
	//cout<<DataSet::getRef("full").N;
	//cout<<endl;



    //return 0;


	//Lookup::get_instance().load();
	start<MyModel>(argc, argv);
	return 0;
}

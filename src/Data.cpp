#include "Data.h"
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>

using namespace std;

Data Data::instance;

Data::Data()
{

}

void Data::load(const char* filename)
{
	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# Error. Couldn't open file "<<filename<<endl;
		return;
	}

	// Empty the vectors
	t.clear();
	y.clear();
	sig.clear();

	int it = 0;
	double temp1, temp2, temp3;
	while(fin>>temp1 && fin>>temp2 && fin>>temp3)
	{
		/*if (it==0 || it==1) {
			it++;
			continue;
		}*/
		t.push_back(temp1);
		y.push_back(temp2);
		sig.push_back(temp3);
		it++;
	}
	cout<<"# Loaded "<<t.size()<<" data points from file "
			<<filename<<endl;
	fin.close();
	//cout<<it<<endl;

	double mean = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
	//cout<<mean<<endl;

	// this is probably a stupid way to substract the mean and convert to m/s
	std::transform( y.begin(), y.end(), y.begin(), std::bind2nd( minus<double>(), mean ) );
	std::transform( y.begin(), y.end(), y.begin(), std::bind2nd( multiplies<double>(), 1000. ) );
	std::transform( y.begin(), y.end(), y.begin(), std::bind2nd( plus<double>(), mean ) );

	// the errorbars just need to be converted to m/s
	std::transform( sig.begin(), sig.end(), sig.begin(), std::bind2nd( multiplies<double>(), 1000. ) );
	
	//for (std::vector<double>::const_iterator i = sig.begin(); i != sig.end(); ++i)
    //std::cout << *i << '\n';
	//std::cout << '\n';


}


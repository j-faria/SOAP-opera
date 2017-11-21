#ifndef DNest4_Data
#define DNest4_Data

#include <vector>
#include <algorithm>
#include <Eigen/Core>

class Data
{
	private:
		std::vector<double> t, y, sig;
		Eigen::MatrixXd ccfdata;

	public:
		Data();
		//void load(const char* filename);
		void load(const char* filename, const char* units, int skip=0);
		int index_fibers;

		// Getters
		int N() const {return t.size();}
		const std::vector<double>& get_t() const { return t; }
		const Eigen::MatrixXd& get_ccfdata() const { return ccfdata; }
		double get_t_min() const { return *min_element(t.begin(), t.end()); }
		double get_t_max() const { return *max_element(t.begin(), t.end()); }

		const std::vector<double>& get_y() const { return y; }
		double get_y_min() const { return *min_element(y.begin(), y.end()); }
		double get_y_max() const { return *max_element(y.begin(), y.end()); }
		
		const std::vector<double>& get_sig() const { return sig; }
		// double get_y_span() const {return abs(get_y_max() - get_y_min());}
		double topslope() const {return abs(get_y_max() - get_y_min()) / (t.back() - t.front()); }

	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

#endif


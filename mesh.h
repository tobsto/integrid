// ****************************************************************************
// mesh 
//
// Description:	Like 'multigrid, but without inverse mapping'
//		M:			length of the grid
//		omega_min:		minimum	
//		omega_max:		maximum	
//		omega:			container for the grid	
//		domega:			container for the weigths	
//
// Tobias Stollenwerk, Last Modification: 29.2.2012
// ****************************************************************************
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<complex>

class mesh
{
	public:
	int M;
	double omega_min;
	double omega_max;
	vector <double> omega;
	vector <double> domega;

	vector<double> spoints;
	bool lendpoint;
	bool rendpoint;
	double omegal;
	double omegar;
	
	// constructor
	mesh();
	// destructor
	~mesh(){};
	// member function
	void add_spoint(double omega0);
	void add_lendpoint(double omegal);
	void add_rendpoint(double omegar);
	void create();
	// for testing purposes
	void testWeights();
	void testMonotony();

	// wrapper for multigrid member functions
	void add_gr_equi(int N, double omega_min, double omega_max, double omega_c);
	void add_gr_equi(int N, double omega_min, double omega_max, double omega_c, string id);
	void add_gr_tan(int N, double omega_min, double omega_max, double omega_c, double c);
	void add_gr_tan(int N, double omega_min, double omega_max, double omega_c, double c, string id);
	void add_gr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0);
	void add_gr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0, string id);

	void add_sgr_equi(int N, double omega_min, double omega_max, double omega_c);
	void add_sgr_equi(int N, double omega_min, double omega_max, double omega_c, string id);
	void add_sgr_tan(int N, double omega_min, double omega_max, double omega_c, double c);
	void add_sgr_tan(int N, double omega_min, double omega_max, double omega_c, double c, string id);
	void add_sgr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0);
	void add_sgr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0, string id);

	void replace_sgr_equi(int N, double omega_min, double omega_max, double omega_c, string id);
	void replace_sgr_tan(int N, double omega_min, double omega_max, double omega_c, double c, string id);
	void replace_sgr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0, string id);
	void replace_gr_equi(int N, double omega_min, double omega_max, double omega_c, string id);
	void replace_gr_tan(int N, double omega_min, double omega_max, double omega_c, double c, string id);
	void replace_gr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0, string id);

	// wrapper for the equi-Gridregion functions
	void add_gr_equi(double omega_min, double omega_max, double domega);
	void add_gr_equi(double omega_min, double omega_max, double domega, string id);
	void add_sgr_equi(double omega_min, double omega_max, double domega);
	void add_sgr_equi(double omega_min, double omega_max, double domega, string id);
	void replace_gr_equi(double omega_min, double omega_max, double domega, string id);
	void replace_sgr_equi(double omega_min, double omega_max, double domega, string id);

	// wrapper for the log-Gridregion functions
	void add_gr_log(double omega_c, double omega1, double domega_min, double domega_max);
	void add_gr_log(double omega_c, double omega1, double domega_min, double domega_max, string id);
	void add_sgr_log(double omega_c, double omega1, double domega_min, double domega_max);
	void add_sgr_log(double omega_c, double omega1, double domega_min, double domega_max, string id);
	void replace_gr_log(double omega_c, double omega1, double domega_min, double domega_max, string id);
	void replace_sgr_log(double omega_c, double omega1, double domega_min, double domega_max, string id);

	void add_sgr(string id, gridRegion & gr);
	void add_gr(string id, gridRegion & gr);
	void replace_sgr(string id, gridRegion & gr);
	void replace_gr(string id, gridRegion & gr);

	void rem_sgr(string id);
	void rem_sgr(unsigned int n);
	void rem_gr(string id);
	void rem_gr(unsigned int n);
	bool sgr_exists(string id);
	bool gr_exists(string id);

	private:
	multigrid mgr;
};

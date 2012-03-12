// Tobias Stollenwerk, Last Modification: 29.2.2012
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<complex>
#include <boost/archive/text_oarchive.hpp> 
#include <boost/archive/text_iarchive.hpp> 
#include <boost/serialization/string.hpp> 
#include <boost/serialization/vector.hpp>


using namespace std;

// contains information about the desired sub grid before it is included into the multigrid
class gridRegion
{
	public:
	double omega_l;
	double omega_c;
	double omega_r;
	double omega_m;
	double omega_p;
	double domega_min_l;
	double domega_min_r;
	int N_d;
	int N_i;
	vector<double> para_d ;
	vector<int> para_i;
	string type;
	string id;
	gridRegion (){};
	// constructor
	gridRegion ( double OMEGA_L, double OMEGA_C, double OMEGA_R, double OMEGA_M, double OMEGA_P, double DOMEGA_MIN_L, double DOMEGA_MIN_R, int N_D, int N_I, vector<double> PARA_D, vector <int> PARA_I, string TYPE, string ID);
	// copy constructor
	gridRegion (const gridRegion & sgr);
	// '=' operator
	gridRegion& operator=(const gridRegion& sgr);
};

// contains information the acutal subgrid after it has been included into the multigrid 
class subgrid
{
	public:
	double omega_l;
	double omega_r;
	int i_l;
	int i_r;
	int N_replaced;
	int N_d;
	int N_i;
	vector<double> para_d ;
	vector<int> para_i;
	string type;
	subgrid (){};
	subgrid ( double OMEGA_L, double OMEGA_R, int I_L, int I_R, int N_REPLACED, int N_D, int N_I, vector<double> PARA_D, vector <int> PARA_I, string TYPE);
};

class multigrid : public grid
{
	public:
	// default constructor
	multigrid ();
	// copy constructor
	multigrid (const multigrid & mgrid);
	// destructor
	~multigrid();
	// '=' operator
	multigrid& operator=(const multigrid& mgrid);

	// member variables
	// container for fundental subgrids
	vector<subgrid> fsubgrids;			
	// container for special subgrids
	vector<subgrid> ssubgrids;			
	// container for fundamental grid regions
	vector<gridRegion> fgridRegions;			
	// container for special grid regions
	vector<gridRegion> sgridRegions;			

	// member functions
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

	void create();
	void createVerbose();

	gridRegion get_sgr(string ID);
	gridRegion get_gr(string ID);
	gridRegion get_sgr(unsigned int n);
	gridRegion get_gr(unsigned int n);

	// for testing purposes
	void list();
	void plot();
	void testMonotony();
	void testInverse();
	void testWeights();

	// inverse
	unsigned int inverse (double omega);
	unsigned int finverse (double omega);

	// exception
	class xBadValues
	{
	};	

	private:
	void addFSubgrid_equi(int N, double omega_min, double omega_max);
	void addFSubgrid_tan(int N, double omega_min, double omega_max, double omega_c, double c);
	void addFSubgrid_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0);
	void addSSubgrid_equi(int N, double omega_min, double omega_max);
	void addSSubgrid_tan(int N, double omega_min, double omega_max, double omega_c, double c);
	void addSSubgrid_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0);
};


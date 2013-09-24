// ****************************************************************************
// This file is part of Integrid.
//
// Integrid is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Integrid is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Integrid.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright 2012 Tobias Stollenwerk
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

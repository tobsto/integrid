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

using namespace std;

class grid
{
	public:
	int M;
	double omega_min;
	double omega_max;
	vector <double> omega;
	vector <double> domega;
	virtual unsigned int inverse (double omega)=0;
};

void plot(grid & g);
void testInverse(grid & g);

class tangrid : public grid
{
	public:
	// grid parameters
	double omega_c;
	double c;

	// implicit grid parameter
	double u_0;
	double u_1;
	double du;

	// constructor
	tangrid (int N_T, double OMEGA_min, double OMEGA_max, double OMEGA_c, double C);
	// inverse
	unsigned int inverse (double omega);

	// exception for bad initial values
	class xBadValues
	{
	};	
};
class equigrid : public grid
{
	public:
	double delta_omega;
	// constructor
	equigrid (int N_e, double OMEGA_min, double OMEGA_max);
	// inverse
	unsigned int inverse (double omega);

	// exception for bad initial values
	class xBadValues
	{
	};	
};
class loggrid : public grid
{
	public:
	// grid parameters
	int N_l;
	int N_r;
	int Nk;
	double omegak;
	double omegak_0;

	// constructed grid parameters
	double domegak_min;
	double domegak;
	// implicit grid parameters;
	double omegak_0m;
	double omegak_0p;
	double c_1;
	double i_1;
	double c_2;
	double i_2;

	// default constructor
	loggrid (){};
	// constructor
	loggrid (int N_L, int N_R, double omega_minp, double omega_maxp, double omegakp, double omegak_0p);
	// destructor
	~loggrid();

	//"=" Operator
	loggrid& operator=(const loggrid& ) ;
	// inverse
	unsigned int inverse (double omega);
	// Testfunktionen für die inversen Abbildungen
	void test_inverse();

	// exception for bad initial values
	class xBadValues
	{
	};	
};

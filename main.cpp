// ****************************************************************************
// Tobias Stollenwerk, stollenwerk@th.physik.uni-bonn.de
// ****************************************************************************

#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<fstream>
#include<vector>
#include<string>
#include<complex>
#include<limits>
#include"grid.h"
#include"multigrid.h"

using namespace std;

void save(vector <double> func, multigrid & mgrid, string filename)
{
	ofstream out;
	out.open(filename.c_str());
	for (int j=0; j!=mgrid.M+1; j++)
	{
		out << scientific << setprecision(15) << mgrid.omega[j] << "\t" <<  func[j] << endl;
	}
	out.close();
}
void saveGrid(grid & mgrid, string filename)
{
	ofstream out;
	out.open(filename.c_str());
	for (int j=0; j!=mgrid.M+1; j++)
	{
		out << scientific << setprecision(15) << j << "\t" << mgrid.omega[j] << endl;
	}
	out.close();
}
double lorentz(double omega, double width, double center)
{
	return width/(M_PI*(pow(omega-center,2.0)+pow(width,2.0)));
}
int main(int argc, char * argv[])
{
	system ("mkdir -p output");

	/*
	// ***********************************************
	// ***********************************************
	// ********* First Example ***********************
	// ***********************************************
	// ***********************************************
	// initialize multigrid named mgrid
	multigrid mgrid;
	// create equidistant grid region from -4 to 4 with resolution 0.01
	mgrid.add_gr_equi(-4.0, 4.0, 0.01);
	// create logarithmically dense grid from with peak point 0.3, half width 0.5, maximal resolution 0.001 and minimal resolution 0.01 
	mgrid.add_gr_log(0.3, 0.5, 0.001, 0.01);
	// create logarithmically dense grid from with peak point 0.6, half width 0.5, maximal resolution 0.001 and minimal resolution 0.01 
	mgrid.add_gr_log(0.6, 0.5, 0.001, 0.01);
	// create the multigrid out of all previously defined grid regions
	mgrid.create();

	// save grid
	saveGrid(mgrid, "output/grid.dat");

	// ***********************************************
	// ********* Calculate Integral ******************
	// ***********************************************

	// calculate a Lorentz function on this grid
	vector<double> func(mgrid.M+1);
	double width=0.1;
	double center=0.3;
	for (int i=0; i<=mgrid.M; i++)
	{
		func[i]=lorentz(mgrid.omega[i], width, center);
	}
	save(func, mgrid, "output/lorentz_function.dat");

	// calculate the integral over this function
	double I=0;
	for (int i=0; i<=mgrid.M; i++)
	{
		I+=func[i]*mgrid.domega[i];
	}
	cout << "Integral over lorentz curve is " << I << endl;
	*/

	/*
	// ***********************************************
	// ***********************************************
	// ********* Grid types **************************
	// ***********************************************
	// ***********************************************
	equigrid egrid(100, -4, 4);
	saveGrid(egrid, "output/equigrid.dat");
	tangrid tgrid(100, -4, 4, 1, 0.5);
	saveGrid(tgrid, "output/tangrid.dat");
	loggrid lgrid(30, 30, -4, 4, 1, 0.4);
	saveGrid(lgrid, "output/loggrid.dat");
	*/

	// ***********************************************
	// ***********************************************
	// ********* Example for adding GR ***************
	// ***********************************************
	// ***********************************************
	multigrid mgrid;
	mgrid.add_gr_equi(100, -1, 1, 0);
	mgrid.add_gr_tan(100, 0.2, 0.5, 0.3, 0.01);
	mgrid.add_gr_log(100, 100, 0.4, 0.7, 0.6, 1E-6, "some_loggrid");
	mgrid.create();
	saveGrid(mgrid, "output/example_adding_gr.dat");

	
	/*
	// mesh: multigrid without inverse mapping
	mesh amesh;
	amesh.add_gr_equi(10, 0.0, 1, 0.5);
	amesh.add_gr_log(0.75, 0.1, 1E-1, 1E-2);
	amesh.add_spoint(0.76);
	for (int j=0; j!=100; j++)
	{
		amesh.add_spoint(0.6+j/100.0);
	}
	amesh.add_spoint(0.101);
	amesh.add_spoint(0.199);
	amesh.add_lendpoint(0.15);
	amesh.add_rendpoint(0.85);
	amesh.create();
	
	ofstream out;
	out.open("output/mesh.dat");
	for (int i=0; i!=amesh.M+1; i++)
	{
		out << i << "\t" << amesh.omega[i] << endl;
	}
	out.close();
	*/
}


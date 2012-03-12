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
#include"mesh.h"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>

using namespace std;
namespace mpi=boost::mpi;

int main(int argc, char * argv[])
{
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
	out.open("mesh.dat");
	for (int i=0; i!=amesh.M+1; i++)
	{
		out << i << "\t" << amesh.omega[i] << endl;
	}
	out.close();
}

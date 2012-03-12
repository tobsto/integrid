// ****************************************************************************
// grid 
//
// Description:	
//		Grid is the basis class for other grids. It has only some 
//		basic properties like:
//		M:			length of the grid
//		omega_min:		minimum	
//		omega_max:		maximum	
//		omega:			container for the grid	
//		domega:			container for the weigths	
//		inverse mapping:	inverse mapping for interpolation	
//
// Programm structure: 
//		'grid' is a abstract class (no instance of this class 
//		possible.
// 
// Tobias Stollenwerk, Last Modification: 29.2.2012
// ****************************************************************************
#include<iostream>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<limits>
#include<vector>
#include<string>
#include<complex>
#include<algorithm>
#include"grid.h"

using namespace std;

void plot(grid & g) 
{
	ofstream out("plotGrid.dat");
	for (int j=0; j!=g.M+1; j++)
	{
		out << scientific << setprecision(15) << j << "\t" << g.omega[j] << endl;
	}
	out.close();
}
void testInverse(grid & g)
{
	for (int i=0; i!=g.M+1; i++)
	{
		if (i!= g.inverse(g.omega[i]))
		{
			cerr << endl;
			cerr << "ERROR: inverse mapping failed at: i=" << i << "\t(";
			cerr << g.inverse(g.omega[i]);
			cerr  << " !=i)" << endl;
			cerr << "Break." << endl;
			//exit(1);
		}
	}
}

tangrid::tangrid (int N_T, double OMEGA_min, double OMEGA_max, double OMEGA_c, double C) : grid()
{
	// set the member variables
	this->M=N_T;
	this->omega_min=OMEGA_min;
	this->omega_max=OMEGA_max;
	this->omega_c=OMEGA_c;
	this->c=C;

	// check
	if (M<3  || omega_c<=omega_min || omega_c>=omega_max || c==0 || omega_min>=omega_max)
	{
		cerr << endl;
		cerr << "Error: tangrid constructor: Bad values." << endl;
		cerr << endl;
		cerr << "Values\t\t\tBreak conditions" << endl;
		cerr << "omega_min:\t"  << omega_min<< "\t(omega_min>=omega_max)" << endl;
		cerr << "omega_max:\t"  << omega_max<< "\t(omega_min>=omega_max)" << endl;
		cerr << "M:\t\t"        << M        << "\t(M<3)" << endl;
		cerr << "omega_c:\t"    << omega_c  << "\t(omega_c<=omega_min) or (omega_c>=omega_max)" << endl;
		cerr << "c:\t\t"        << c        << "\t(c==0)" << endl;
		cerr << endl;
		cerr << "Break." << endl;
		throw xBadValues();
	}
	// initialize tangential grid
	this->omega=vector <double> (M+1);	
	this->domega=vector <double> (M+1);	
	this-> u_0=atan((omega_min-omega_c)/c);
	this-> u_1=atan((omega_max-omega_c)/c);
	this-> du=(u_1-u_0)/double(M);
	for (int j=1; j!=M; j++)
	{
		omega[j]=c*tan(u_0 + du * j) + omega_c;	
	}
	omega[0]=omega_min;	
	omega[M]=omega_max;	

	// weights from trapez rule
	for (int i=1; i<=M-1; i++)
	{
		domega[i]=0.5*(omega[i+1]-omega[i-1]);	
	}
	domega[0]=0.5*(omega[1]-omega[0]);	
	domega[M]=0.5*(omega[M]-omega[M-1]);	
}
// inverse
unsigned int tangrid::inverse (double omega)
{
	// rounding by adding 0.5
	if (omega>=omega_min && omega<=omega_max )
	{
		return int((atan((omega-omega_c)/c)-u_0)/du+0.5);
	}
	else
	{
		cerr << endl;
		cerr << "ERROR: tangrid inverse: argument out of range. Break." << endl;
		cerr << scientific << setprecision(30) << omega << " not inside [" << omega_min << ", " << omega_max << "]" << endl;
		exit(1);
	}
}
equigrid::equigrid (int N_e, double OMEGA_min, double OMEGA_max) : grid()
{
	// set the member variables
	this->M=N_e;
	this->omega_min=OMEGA_min;
	this->omega_max=OMEGA_max;

	// check
	if (M<3  || omega_min>=omega_max)
	{
		cerr << endl;
		cerr << "Error: equigrid constructor: Bad values." << endl;
		cerr << endl;
		cerr << "Values\t\t\tBreak conditions" << endl;
		cerr << "omega_min:\t"  << omega_min<< "\t(omega_min>=omega_max)" << endl;
		cerr << "omega_max:\t"  << omega_max<< "\t(omega_min>=omega_max)" << endl;
		cerr << "M:\t\t"        << M        << "\t(M<3)" << endl;
		cerr << endl;
		cerr << "Break." << endl;
		throw xBadValues();
	}
	this->delta_omega=(omega_max-omega_min)/double(M);
	// initialize tangential grid
	this->omega=vector <double> (M+1);	
	this->domega=vector <double> (M+1);	
	omega[0]=omega_min;
	for (int j=1; j!=M; j++)
	{
		omega[j]=omega_min + j*delta_omega;
	}
	omega[M]=omega_max;

	// weights from trapez rule
	for (int i=1; i<=M-1; i++)
	{
		domega[i]=0.5*(omega[i+1]-omega[i-1]);	
	}
	domega[0]=0.5*(omega[1]-omega[0]);	
	domega[M]=0.5*(omega[M]-omega[M-1]);	
}
// inverse
unsigned int equigrid::inverse (double omega)
{
	// rounding by adding 0.5
	if (omega>=omega_min && omega<=omega_max )
	{
		return int((omega-omega_min)/delta_omega+0.5);
	}
	else
	{
		cerr << endl;
		cerr << "ERROR: equigrid inverse: argument out of range. Break." << endl;
		cerr << scientific << setprecision(30) << omega << " not inside [" << omega_min << ", " << omega_max << "]" << endl;
		exit(1);
	}
}
// constructor
loggrid::loggrid (int N_L, int N_R, double OMEGA_min, double OMEGA_max, double OMEGAK, double OMEGAK_0) : grid()
{
	// set the member variables
	this->N_l=N_L;
	this->N_r=N_R;
	this->omega_min=OMEGA_min;
	this->omega_max=OMEGA_max;
	this->omegak=OMEGAK;
	this->omegak_0=OMEGAK_0;

	// boundaries for the grid
	this->omegak_0m=omegak-omegak_0;
	this->omegak_0p=omegak+omegak_0;
	
	// check
	if (N_l<3 || N_r<3 ||  omegak<=omega_min || omegak>=omega_max || omegak_0<1E-15 || omega_min>=omega_max)
	{
		cerr << endl;
		cerr << "Error: loggrid constructor: Bad values." << endl;
		cerr << endl;
		cerr << "Values\t\t\tBreak conditions" << endl;
		cerr << "omega_min:\t"  << omega_min<< "\t(omega_min>=omega_max)" << endl;
		cerr << "omega_max:\t"  << omega_max<< "\t(omega_min>=omega_max)" << endl;
		cerr << "N_l:\t\t"       << N_l       << "\t(N_l<3)" << endl;
		cerr << "N_r:\t\t"       << N_r       << "\t(N_r<3)" << endl;
		cerr << "omegak:\t\t"   << omegak   << "\t(omegak<=omega_min) or (omegak>=omega_max)" << endl;
		cerr << "omegak_0:\t"   << omegak_0 << "\t(omegak_0<1E-15))" << endl;
		cerr << endl;
		cerr << "Break." << endl;
		throw xBadValues();
	}
	
	// **************************************************************************
	// **************************************************************************
	// *********    initialize singly peaked grid and weights   *****************
	// **************************************************************************
	// **************************************************************************

	// get grid parameters (from boundary conditions)
	
	this->c_1 = log((omegak-omega_min)/(omegak_0)) / (N_l-1);
	this->i_1 =  log(omegak-omega_min)/c_1;

	this->c_2 = log((omega_max-omegak)/(omegak_0)) / (N_r);
	this->i_2 = -log(omega_max-omegak)/c_2 + N_r;

	// calculate maximal resolution in logarithmic region
	this->domegak_min=min(exp(-c_1*(N_l-1-i_1))*(exp(c_1)-1), exp(-c_2*i_2)*(exp(c_2)-1));
	// check if grid difference is beyond machine precision
	if (min(fabs(domegak_min/(omegak+omegak_0)), fabs(domegak_min/(omegak-omegak_0)))<=numeric_limits<double>::epsilon())
	{

		cerr << endl;
		cerr << "Error: loggrid constructor error: Smallest grid difference is beyond machine precision. Break." << endl;
		cerr << endl;
		cerr << scientific << setprecision(numeric_limits<double>::digits10);
		cerr << "omegak:\t\t\t\t\t" << omegak  << endl;
		cerr << "omegak_0:\t\t\t\t" << omegak_0 << endl;
		cerr << "domegak_min:\t\t\t\t" << domegak_min << endl;
		cerr << "min(fabs(domegak_min/(omegak+omegak_0)),\nfabs(domegak_min/(omegak-omegak_0))):\t\t" << min(fabs(domegak_min/(omegak+omegak_0)), fabs(domegak_min/(omegak-omegak_0))) << endl;
		cerr << endl;
		cerr << "machine precision:\t\t\t" << numeric_limits<double>::epsilon()  << endl;
		throw xBadValues();
	}
	// get number of grid points in peak regions 
	this->Nk=max(int((omegak_0p-omegak_0m)/domegak_min)-1, 0);
	this->domegak=(omegak_0p-omegak_0m)/double(Nk+1);
	//cout << "peak grid at: " << N1+N2 << " to " << max(N1+N2+Nk-1,N1+N2) << endl;

	// allocate memory for
	// grid and weights
	this->M=N_l+N_r+Nk;
	this->omega=vector <double> (M+1);	
	this->domega=vector <double> (M+1);	

	//calculate grid and weights
	// Region I
	omega[0]=omega_min;
	for (int i=1; i!=N_l-1; i++)
	{
		omega[i]=-exp(-c_1*(i-i_1))+omegak;
	}
	omega[N_l-1]=omegak_0m;
	// Peak Region 
	if (Nk>0)
	{
		for (int i=1; i<=Nk; i++)
		{
			omega[N_l-1+i]=omegak_0m+domegak * i;
		}
	}
	// Region II
	omega[N_l+Nk]=omegak_0p;
	for (int i=N_l+Nk+1; i!=Nk+N_l+N_r; i++)
	{
		omega[i]=exp(c_2*(i-i_2-N_l-Nk))+omegak;
	}
	omega[N_l+N_r+Nk]=omega_max;

	// trapez rule
	for (int j=1; j<M; j++)
	{
		domega[j]=0.5*(omega[j+1]-omega[j-1]);	
	}
	domega[0]=0.5*(omega[1]-omega[0]);	
	domega[M]=0.5*(omega[M]-omega[M-1]);	
}
// destructor
loggrid::~loggrid() 
{
}
// inverse of singly peaked omega grid
unsigned int loggrid::inverse (double omega)
{
	// rounding by adding 0.5
	if (omega>=omega_min && omega<=omegak_0m )
	{
		return int(i_1-log(omegak-omega)/c_1+0.5);
	}
	else if (omega>omegak_0m && omega<omegak_0p )
	{
		if (Nk>0)
		{
			return int((omega-omegak_0m)/domegak +0.5)+N_l-1;
		}
		else
		{
			return N_l;
		}
	}
	else if (omega>=omegak_0p && omega<=omega_max )
	{
		return int(i_2+log(omega-omegak)/c_2+0.5)+N_l+Nk;
	}
	else
	{
		cerr << endl;
		cerr << "ERROR loggrid: single-peaked inverse: argument out of range. Break." << endl;
		cerr << scientific << setprecision(30) << omega << " not inside [" << omega_min << ", " << omega_max << "]" << endl;
		exit(1);
	}
}
void loggrid::test_inverse()
{
	for (int i=0; i!=M+1; i++)
	{
		if (i!= inverse(omega[i]))
		{
			cerr << endl;
			cerr << "ERROR loggrid: inverse mapping failed at: i=" << i << "\t(";
			cerr << inverse(omega[i]);
			cerr  << " !=i)" << endl;
			cerr << "Break." << endl;
			//exit(1);
		}
	}
}
// "=" Operator
loggrid& loggrid::operator=(const loggrid& A) 
{
        this -> N_l=A.N_l;
        this -> N_r=A.N_r;
        this -> Nk=A.Nk;
        this -> omega_min=A.omega_min;
        this -> omega_max=A.omega_max;
        this -> omegak=A.omegak;
        this -> omegak_0=A.omegak_0;
        this -> domegak_min=A.domegak_min;
        this -> domegak=A.domegak;
        this -> omegak_0m=A.omegak_0m;
        this -> omegak_0p=A.omegak_0p;
        this -> c_1=A.c_1;
        this -> i_1=A.i_1;
        this -> c_2=A.c_2;
        this -> i_2=A.i_2;
	this->omega.resize(A.M+1);
	this->domega.resize(A.M+1);
        this -> M=A.M;

	for (int j=0; j!=A.M+1; j++)
	{
		this->omega[j]=A.omega[j];
		this->domega[j]=A.domega[j];
	}
        return *this;
}

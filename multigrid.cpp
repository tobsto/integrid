// ****************************************************************************
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
#include"multigrid.h"

using namespace std;

gridRegion::gridRegion ( double OMEGA_L, double OMEGA_C, double OMEGA_R, double OMEGA_M, double OMEGA_P, double DOMEGA_MIN_L, double DOMEGA_MIN_R, int N_D, int N_I, vector<double> PARA_D, vector <int> PARA_I, string TYPE, string ID)
{
	this->type=TYPE;
	this->id=ID;
	this->omega_l=OMEGA_L;
	this->omega_c=OMEGA_C;
	this->omega_r=OMEGA_R;
	this->omega_m=OMEGA_M;
	this->omega_p=OMEGA_P;
	this->domega_min_l=DOMEGA_MIN_L;
	this->domega_min_r=DOMEGA_MIN_R;
	this->N_d=N_D;
	this->N_i=N_I;
	for (int n=0; n!=N_d; n++)
	{
		this->para_d.push_back(PARA_D[n]);	
	}
	for (int n=0; n!=N_i; n++)
	{
		this->para_i.push_back(PARA_I[n]);	
	}
}
// < operator for special grid regions (compares the centers of the grid regions)
bool lesserLeft( const gridRegion& a, const gridRegion& b)
{
	return a.omega_l < b.omega_l;
}
bool lesserCenter( const gridRegion& a, const gridRegion& b)
{
	return a.omega_c < b.omega_c;
}
subgrid::subgrid ( double OMEGA_L, double OMEGA_R, int I_L, int I_R, int N_REPLACED, int N_D, int N_I, vector<double> PARA_D, vector <int> PARA_I, string TYPE)
{
	this->type=TYPE;
	this->omega_l=OMEGA_L;
	this->omega_r=OMEGA_R;
	this->i_l=I_L;
	this->i_r=I_R;
	this->N_d=N_D;
	this->N_i=N_I;
	this->N_replaced=N_REPLACED;
	for (int n=0; n!=N_d; n++)
	{
		this->para_d.push_back(PARA_D[n]);	
	}
	for (int n=0; n!=N_i; n++)
	{
		this->para_i.push_back(PARA_I[n]);	
	}
}
// constructor
multigrid::multigrid()
{
}
// destructor
multigrid::~multigrid() 
{
}

void multigrid::list() 
{
	cout << "# FUNDAMENTAL GRID REGIONS:" << endl;
	for (int j=0; j!=this->fgridRegions.size(); j++)
	{
		cout << "# type: " << this->fgridRegions[j].type << "\tomega_l =" << this->fgridRegions[j].omega_l << "\tomega_r =" << this->fgridRegions[j].omega_r <<   "\tomega_c =" << this->fgridRegions[j].omega_c<<   "\tid: " << this->fgridRegions[j].id;
		if (this->fgridRegions[j].type=="equi")
		{
			cout << "\tN=" << this->fgridRegions[j].para_i[0] <<   "\tdomega=" << this->fgridRegions[j].para_d[0];
		}
		else if (this->fgridRegions[j].type=="tan")
		{
			cout << "\tN=" << this->fgridRegions[j].para_i[0] <<   "\tc=" << this->fgridRegions[j].para_d[1];
		}
		else if (this->fgridRegions[j].type=="log")
		{
			cout << "\tNl=" << this->fgridRegions[j].para_i[0] <<   "\tNr=" << this->fgridRegions[j].para_i[1] <<   "\tomega0=" << this->fgridRegions[j].para_d[1];
		}
		cout << endl;
	}
	cout << "# SPECIAL GRID REGIONS:" << endl;
	for (int j=0; j!=this->sgridRegions.size(); j++)
	{
		cout << "# type: " << this->sgridRegions[j].type << "\tomega_l =" << this->sgridRegions[j].omega_l << "\tomega_r =" << this->sgridRegions[j].omega_r <<   "\tomega_c =" << this->sgridRegions[j].omega_c<<   "\tid: " << this->sgridRegions[j].id;
		if (this->sgridRegions[j].type=="equi")
		{
			cout << "\tN=" << this->sgridRegions[j].para_i[0] <<   "\tdomega=" << this->sgridRegions[j].para_d[0];
		}
		else if (this->sgridRegions[j].type=="tan")
		{
			cout << "\tN=" << this->sgridRegions[j].para_i[0] <<   "\tc=" << this->sgridRegions[j].para_d[1];
		}
		else if (this->sgridRegions[j].type=="log")
		{
			cout << "\tNl=" << this->sgridRegions[j].para_i[0] <<   "\tNr=" << this->sgridRegions[j].para_i[1] <<   "\tomega0=" << this->sgridRegions[j].para_d[1];
		}
		cout << endl;
	}
	cout << "# FUNDAMENTAL SUBGRIDS:" << endl;
	for (int j=0; j!=this->fsubgrids.size(); j++)
	{
		cout << "# type: " << this->fsubgrids[j].type << "\tomega_l =" << this->fsubgrids[j].omega_l << "\tomega_r =" << this->fsubgrids[j].omega_r;
		if (this->fsubgrids[j].type=="equi")
		{
			cout << "\tN=" << this->fsubgrids[j].para_i[0] <<   "\tdomega=" << this->fsubgrids[j].para_d[0];
		}
		else if (this->fsubgrids[j].type=="tan")
		{
			cout << "\tN=" << this->fsubgrids[j].para_i[0] <<   "\tc=" << this->fsubgrids[j].para_d[1];
		}
		else if (this->fsubgrids[j].type=="log")
		{
			cout << "\tNl=" << this->fsubgrids[j].para_i[1] <<   "\tNr=" << this->fsubgrids[j].para_i[3] <<   "\tomega0=" << this->fsubgrids[j].para_d[1];
		}
		cout << endl;
	}
	cout << "# SPECIAL SUBGRIDS:" << endl;
	for (int j=0; j!=this->ssubgrids.size(); j++)
	{
		cout << "# type: " << this->ssubgrids[j].type << "\tomega_l =" << this->ssubgrids[j].omega_l << "\tomega_r =" << this->ssubgrids[j].omega_r;
		if (this->ssubgrids[j].type=="equi")
		{
			cout << "\tN=" << this->ssubgrids[j].para_i[0] <<   "\tdomega=" << this->ssubgrids[j].para_d[0];
		}
		else if (this->ssubgrids[j].type=="tan")
		{
			cout << "\tN=" << this->ssubgrids[j].para_i[0] <<   "\tc=" << this->ssubgrids[j].para_d[1];
		}
		else if (this->ssubgrids[j].type=="log")
		{
			cout << "\tNl=" << this->ssubgrids[j].para_i[1] <<   "\tNr=" << this->ssubgrids[j].para_i[3] <<   "\tomega0=" << this->ssubgrids[j].para_d[1];
		}
		cout << endl;
	}
}
void multigrid::plot() 
{
	ofstream out("plotGrid.dat");
	for (int j=0; j!=this->M+1; j++)
	{
		out << scientific << setprecision(15) << j << "\t" << this->omega[j] << endl;
	}
	out.close();
}
void multigrid::testMonotony()
{
	for (int j=0; j<=M-1; j++)
	{
		if (omega[j+1]<omega[j])
		{
			cerr << endl;
			cerr << "Error. Keine Monotonie im Integrationsgitter." << endl;
			cerr << "j: " << j <<  endl;	
			throw 1;
		}
		if (omega[j+1]==omega[j])
		{
			cerr << endl;
			cerr << "Error. Keine Strenge Monotonie im Integrationsgitter." << endl;
			cerr << "j: " << j <<  endl;	
			throw 1;
		}
	}
}
void multigrid::testInverse()
{
	for (int j=0; j!=M+1; j++)
	{
		if (j!= inverse(omega[j]))
		{
			cerr << endl;
			cerr << "ERROR: inverse mapping failed at: j=" << j << "\t(";
			cerr << inverse(omega[j]);
			cerr  << " !=j)" << endl;
			cerr << " at omega= " << scientific << setprecision(15)  << omega[j] << endl;
			cerr << "Break." << endl;
			throw 1;
		}
	}
}
void multigrid::testWeights()
{
	for (int j=0; j!=M+1; j++)
	{
		if (domega[j]<=numeric_limits<double>::epsilon())
		{
			cerr << endl;
			cerr << "ERROR: test multigrid weights: Smallest grid difference is beyond machine precision at: j=" << j ;
			cerr << " , domega= " << scientific << setprecision(15)  << omega[j] << endl;
			cerr << "Break." << endl;
			throw 1;
		}
	}
}

// inverse of the fundamental grid
unsigned int multigrid::finverse (double omega)
{
	// calculate number of grid points from other grid regions which are below omega
	int nShift=0;

	for (int n=int(fsubgrids.size()-1); n>=1; n--)
	{
		if (omega>=fsubgrids[n].omega_l)
		{
			if (omega<=fsubgrids[n].omega_r)
			{
				if (fsubgrids[n].type=="equi")
				{
					// subgrid.para_d[0] is delta_omega
					//cout << omega << "\tin region " << n << "\t" << int((omega-fsubgrids[n].omega_l)/fsubgrids[n].para_d[0]+0.5) << "\t" << fsubgrids[n].i_l << "\t" << int((omega-fsubgrids[n].omega_l)/fsubgrids[n].para_d[0]+0.5)+fsubgrids[n].i_l << endl;
					return int((omega-fsubgrids[n].omega_l)/fsubgrids[n].para_d[0]+0.5)+fsubgrids[n].i_l;
				}
				else if (fsubgrids[n].type=="tan")
				{
					// subgrid.para_d[0] is omega_c
					// subgrid.para_d[1] is c
					return int((atan((omega-fsubgrids[n].para_d[0])/fsubgrids[n].para_d[1])-fsubgrids[n].para_d[2])/fsubgrids[n].para_d[4]+0.5)+fsubgrids[n].i_l;
				}
				else if (fsubgrids[n].type=="log")
				{
					// subgrid.para_i[0] is M
					// subgrid.para_i[1] is N_l
					// subgrid.para_i[2] is Nk
					// subgrid.para_i[3] is N_r
					// subgrid.para_d[0] is omegak
					// subgrid.para_d[1] is omegak_0
					// subgrid.para_d[2] is domegak_min
					// subgrid.para_d[3] is domegak
					// subgrid.para_d[4] is omegak_0m
					// subgrid.para_d[5] is omegak_0p
					// subgrid.para_d[6] is c_1
					// subgrid.para_d[7] is i_1
					// subgrid.para_d[8] is c_2
					// subgrid.para_d[9] is i_2
					if (omega>=fsubgrids[n].omega_l && omega<=fsubgrids[n].para_d[4] )
					{
						return int(fsubgrids[n].para_d[7]-log(fsubgrids[n].para_d[0]-omega)/fsubgrids[n].para_d[6]+0.5)+fsubgrids[n].i_l;
					}
					else if (omega>fsubgrids[n].para_d[4] && omega<fsubgrids[n].para_d[5] )
					{
						if (fsubgrids[n].para_i[2]>0)
						{
							return int((omega-fsubgrids[n].para_d[4])/fsubgrids[n].para_d[3] +0.5)+fsubgrids[n].para_i[1]-1+fsubgrids[n].i_l;
						}
						else
						{
							return fsubgrids[n].para_i[1]+fsubgrids[n].i_l;
						}
					}
					else if (omega>=fsubgrids[n].para_d[5] && omega<=fsubgrids[n].omega_r )
					{
						return int(fsubgrids[n].para_d[9]+log(omega-fsubgrids[n].para_d[0])/fsubgrids[n].para_d[8]+0.5)+fsubgrids[n].para_i[1]+fsubgrids[n].para_i[2]+fsubgrids[n].i_l;
					}
				}
				else
				{
					cerr << "multigrid: finverse: An serious error has occured. Subgrid seems not to be of known type. Break." << endl;
					exit(1);
				}
			}
			else
			{
				nShift+=fsubgrids[n].N_replaced;
			}	
		}
	}
	// if omega is not inside one of the additional grid regions, it has to be in the first one
	if (omega>=fsubgrids[0].omega_l && omega<=fsubgrids[0].omega_r)
	{
		if (fsubgrids[0].type=="equi")
		{
			return int((omega-fsubgrids[0].omega_l)/fsubgrids[0].para_d[0]+0.5)+nShift;
		}
		else if (fsubgrids[0].type=="tan")
		{
			return int((atan((omega-fsubgrids[0].para_d[0])/fsubgrids[0].para_d[1])-fsubgrids[0].para_d[2])/fsubgrids[0].para_d[4]+0.5)+nShift;
		}
		else if (fsubgrids[0].type=="log")
		{
			if (omega>=fsubgrids[0].omega_l && omega<=fsubgrids[0].para_d[4] )
			{
				return int(fsubgrids[0].para_d[7]-log(fsubgrids[0].para_d[0]-omega)/fsubgrids[0].para_d[6]+0.5)+nShift;
			}
			else if (omega>fsubgrids[0].para_d[4] && omega<fsubgrids[0].para_d[5] )
			{
				if (fsubgrids[0].para_i[2]>0)
				{
					return int((omega-fsubgrids[0].para_d[4])/fsubgrids[0].para_d[3] +0.5)+fsubgrids[0].para_i[1]-1+nShift;
				}
				else
				{
					return fsubgrids[0].para_i[1]+nShift;
				}
			}
			else if (omega>=fsubgrids[0].para_d[5] && omega<=fsubgrids[0].omega_r )
			{
				return int(fsubgrids[0].para_d[9]+log(omega-fsubgrids[0].para_d[0])/fsubgrids[0].para_d[8]+0.5)+fsubgrids[0].para_i[1]+fsubgrids[0].para_i[2]+nShift;
			}
		}
		else
		{
			cerr << "multigrid: finverse: An serious error has occured. Subgrid seems not to be of known type. Break." << endl;
			exit(1);
		}
	}
	// if omega is not inside any grid region
	else	
	{
		cerr << endl;
		cerr << "ERROR: multigrid finverse: argument out of range. Break." << endl;
		cerr << scientific << setprecision(30) << omega << " not inside [" << this->omega_min << ", " << this->omega_max << "]" << endl;
		exit(1);
	}
}
unsigned int multigrid::inverse (double omega)
{
	// calculate number of grid points from other grid regions which are below omega
	int nShift=0;

	for (int n=int(ssubgrids.size()-1); n>=0; n--)
	{
		if (omega>=ssubgrids[n].omega_l)
		{
			if (omega<=ssubgrids[n].omega_r)
			{
				if (ssubgrids[n].type=="equi")
				{
					// subgrid.para_d[0] is delta_omega
					//cout << omega << "\tin region " << n << "\t" << int((omega-ssubgrids[n].omega_l)/ssubgrids[n].para_d[0]+0.5) << "\t" << ssubgrids[n].i_l << "\t" << int((omega-ssubgrids[n].omega_l)/ssubgrids[n].para_d[0]+0.5)+ssubgrids[n].i_l << endl;
					return int((omega-ssubgrids[n].omega_l)/ssubgrids[n].para_d[0]+0.5)+ssubgrids[n].i_l;
				}
				else if (ssubgrids[n].type=="tan")
				{
					// subgrid.para_d[0] is omega_c
					// subgrid.para_d[1] is c
					return int((atan((omega-ssubgrids[n].para_d[0])/ssubgrids[n].para_d[1])-ssubgrids[n].para_d[2])/ssubgrids[n].para_d[4]+0.5)+ssubgrids[n].i_l;
				}
				else if (ssubgrids[n].type=="log")
				{
					// subgrid.para_i[0] is M
					// subgrid.para_i[1] is N_l
					// subgrid.para_i[2] is Nk
					// subgrid.para_i[3] is N_r
					// subgrid.para_d[0] is omegak
					// subgrid.para_d[1] is omegak_0
					// subgrid.para_d[2] is domegak_min
					// subgrid.para_d[3] is domegak
					// subgrid.para_d[4] is omegak_0m
					// subgrid.para_d[5] is omegak_0p
					// subgrid.para_d[6] is c_1
					// subgrid.para_d[7] is i_1
					// subgrid.para_d[8] is c_2
					// subgrid.para_d[9] is i_2
					if (omega>=ssubgrids[n].omega_l && omega<=ssubgrids[n].para_d[4] )
					{
						return int(ssubgrids[n].para_d[7]-log(ssubgrids[n].para_d[0]-omega)/ssubgrids[n].para_d[6]+0.5)+ssubgrids[n].i_l;
					}
					else if (omega>ssubgrids[n].para_d[4] && omega<ssubgrids[n].para_d[5] )
					{
						if (ssubgrids[n].para_i[2]>0)
						{
							return int((omega-ssubgrids[n].para_d[4])/ssubgrids[n].para_d[3] +0.5)+ssubgrids[n].para_i[1]-1+ssubgrids[n].i_l;
						}
						else
						{
							return ssubgrids[n].para_i[1]+ssubgrids[n].i_l;
						}
					}
					else if (omega>=ssubgrids[n].para_d[5] && omega<=ssubgrids[n].omega_r )
					{
						return int(ssubgrids[n].para_d[9]+log(omega-ssubgrids[n].para_d[0])/ssubgrids[n].para_d[8]+0.5)+ssubgrids[n].para_i[1]+ssubgrids[n].para_i[2]+ssubgrids[n].i_l;
					}
				}
				else
				{
					cerr << "multigrid: inverse: An serious error has occured. Subgrid seems not to be of known type. Break." << endl;
					exit(1);
				}
			}
			else
			{
				nShift+=ssubgrids[n].N_replaced;
			}	
		}
	}
	// if omega is not inside one of the special grid regions, it has to be in the fundamental grid 
	if (omega>=this->omega_min && omega<=this->omega_max)
	{
		return finverse(omega)+nShift;
	}
	// if omega is not inside any grid region
	else	
	{
		cerr << endl;
		cerr << "ERROR: multigrid inverse: argument out of range. Break." << endl;
		cerr << scientific << setprecision(30) << omega << " not inside [" << this->omega_min << ", " << this->omega_max << "]" << endl;
		exit(1);
	}
}


// add fundamental equidistant grid region with id
void multigrid::add_gr_equi(int N_e, double omega_l, double omega_r, double omega_c, string id)
{
	double ntol=numeric_limits<double>::epsilon();
	if (N_e<3  || omega_l+ntol*N_e>=omega_r)
	{
		cerr << endl;
		cerr << "Error: multigrid: add equigrid: Bad values." << endl;
		cerr << endl;
		cerr << "Values\t\t\tBreak conditions" << endl;
		cerr << "omega_l:\t"  << omega_l<< "\t(omega_l+ntol*N_e>=omega_r)" << endl;
		cerr << "omega_r:\t"  << omega_r<< "\t(omega_l+ntol*N_e>=omega_r)" << endl;
		cerr << "ntol:\t"  << ntol<< "\t(omega_l+ntol*N_e>=omega_r)" << endl;
		cerr << "N_e:\t\t"        << N_e        << "\t(N_e<3)" << endl;
		cerr << endl;
		cerr << "Break." << endl;
		throw xBadValues();
	}
	// container for grid region parameter
	int N_i=1;
	int N_d=1;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// save grid parameter
	para_i[0]=N_e;
	double delta_omega=(omega_r-omega_l)/double(N_e);
	para_d[0]=delta_omega;
	double domega_min_l=delta_omega;
	double domega_min_r=delta_omega;
	double omega_m=omega_c-3*delta_omega;
	double omega_p=omega_c+3*delta_omega;
	
	this->fgridRegions.push_back( gridRegion(omega_l, omega_c, omega_r, omega_m, omega_p, domega_min_l, domega_min_r, N_d, N_i, para_d, para_i, "equi", id));	
}
// add special equidistant grid region
void multigrid::add_sgr_equi(int N_e, double omega_l, double omega_r, double omega_c, string id)
{
	double ntol=numeric_limits<double>::epsilon();
	if (N_e<3  || omega_l+ntol*N_e>=omega_r)
	{
		cerr << endl;
		cerr << "Error: multigrid: add equigrid: Bad values." << endl;
		cerr << endl;
		cerr << "Values\t\t\tBreak conditions" << endl;
		cerr << "omega_l:\t"  << omega_l<< "\t(omega_l+ntol*N_e>=omega_r)" << endl;
		cerr << "omega_r:\t"  << omega_r<< "\t(omega_l+ntol*N_e>=omega_r)" << endl;
		cerr << "ntol:\t"  << ntol<< "\t(omega_l+ntol*N_e>=omega_r)" << endl;
		cerr << "N_e:\t\t"        << N_e        << "\t(N_e<3)" << endl;
		cerr << endl;
		cerr << "Break." << endl;
		throw xBadValues();
	}
	// container for grid region parameter
	int N_i=1;
	int N_d=1;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// save grid parameter
	para_i[0]=N_e;
	double delta_omega=(omega_r-omega_l)/double(N_e);
	para_d[0]=delta_omega;
	double domega_min_l=delta_omega;
	double domega_min_r=delta_omega;
	double omega_m=omega_c-3*delta_omega;
	double omega_p=omega_c+3*delta_omega;
	
	this->sgridRegions.push_back( gridRegion(omega_l, omega_c, omega_r, omega_m, omega_p, domega_min_l, domega_min_r, N_d, N_i, para_d, para_i, "equi", id));	
}


// add fundamental tangential grid region
void multigrid::add_gr_tan(int N_t, double omega_l, double omega_r, double omega_c, double c, string id) 
{
	double ntol=numeric_limits<double>::epsilon();
	if (N_t<3  || omega_c<=omega_l || omega_c>=omega_r || c==0 || omega_l+ntol*N_t>=omega_r)
	{
		cerr << endl;
		cerr << "Error: multigrid: add tangrid: Bad values." << endl;
		cerr << endl;
		cerr << "Values\t\t\tBreak conditions" << endl;
		cerr << "omega_l:\t"  << omega_l<< "\t(omega_l+ntol*N_t>=omega_r)" << endl;
		cerr << "omega_r:\t"  << omega_r<< "\t(omega_l+ntol*N_t>=omega_r)" << endl;
		cerr << "ntol:\t"  << ntol<< "\t(omega_l+ntol*N_t>=omega_r)" << endl;
		cerr << "N_t:\t\t"        << N_t        << "\t(N_t<3)" << endl;
		cerr << "omega_c:\t"    << omega_c  << "\t(omega_c<=omega_l) or (omega_c>=omega_r)" << endl;
		cerr << "c:\t\t"        << c        << "\t(c==0)" << endl;
		cerr << endl;
		cerr << "Break." << endl;
		throw xBadValues();
	}
	// container for grid region parameter
	int N_i=1;
	int N_d=3;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// save grid parameter
	para_i[0]=N_t;
	para_d[0]=omega_c;
	para_d[1]=c;
	double u_0=atan((omega_l-omega_c)/c);
	double u_1=atan((omega_r-omega_c)/c);
	double du=(u_1-u_0)/double(N_t);
	para_d[2]=du;
	double domega_min_l=c*du;
	double domega_min_r=c*du;
	double omega_m=omega_c-3*c*du;
	double omega_p=omega_c+3*c*du;
	
	this->fgridRegions.push_back( gridRegion(omega_l, omega_c, omega_r, omega_m, omega_p, domega_min_l, domega_min_r, N_d, N_i, para_d, para_i, "tan", id));
}
// add special tangential grid region
void multigrid::add_sgr_tan(int N_t, double omega_l, double omega_r, double omega_c, double c, string id) 
{
	double ntol=numeric_limits<double>::epsilon();
	if (N_t<3  || omega_c<=omega_l || omega_c>=omega_r || c==0 || omega_l+ntol*N_t>=omega_r)
	{
		cerr << endl;
		cerr << "Error: multigrid: add tangrid: Bad values." << endl;
		cerr << endl;
		cerr << "Values\t\t\tBreak conditions" << endl;
		cerr << "omega_l:\t"  << omega_l<< "\t(omega_l+ntol*N_t>=omega_r)" << endl;
		cerr << "omega_r:\t"  << omega_r<< "\t(omega_l+ntol*N_t>=omega_r)" << endl;
		cerr << "ntol:\t"  << ntol<< "\t(omega_l+ntol*N_t>=omega_r)" << endl;
		cerr << "N_t:\t\t"        << N_t        << "\t(N_t<3)" << endl;
		cerr << "omega_c:\t"    << omega_c  << "\t(omega_c<=omega_l) or (omega_c>=omega_r)" << endl;
		cerr << "c:\t\t"        << c        << "\t(c==0)" << endl;
		cerr << endl;
		cerr << "Break." << endl;
		throw xBadValues();
	}
	// container for grid region parameter
	int N_i=1;
	int N_d=3;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// save grid parameter
	para_i[0]=N_t;
	para_d[0]=omega_c;
	para_d[1]=c;
	double u_0=atan((omega_l-omega_c)/c);
	double u_1=atan((omega_r-omega_c)/c);
	double du=(u_1-u_0)/double(N_t);
	para_d[2]=du;
	double domega_min_l=c*du;
	double domega_min_r=c*du;
	double omega_m=omega_c-3*c*du;
	double omega_p=omega_c+3*c*du;
	
	this->sgridRegions.push_back( gridRegion(omega_l, omega_c, omega_r, omega_m, omega_p, domega_min_l, domega_min_r, N_d, N_i, para_d, para_i, "tan", id));
}

// add fundamental logarithmic grid region
void multigrid::add_gr_log(int N_l, int N_r, double omega_l, double omega_r, double omegak, double omegak_0, string id) 
{
	double ntol=numeric_limits<double>::epsilon();
	if (N_l<3 || N_r<3 ||  omegak<=omega_l || omegak>=omega_r || omegak_0<1E-15 || omegak-omegak_0<=omega_l+ntol*N_l || omegak+omegak_0>=omega_r-ntol*N_r || omega_l+ntol*(N_l+N_r)>=omega_r)
	{
		cerr << endl;
		cerr << "Error: multigrid: add loggrid: Bad values." << endl;
		cerr << endl;
		cerr << "Values\t\t\tBreak conditions" << endl;
		cerr << "omega_l:\t"  << omega_l<< "\t(omega_l+ntol*(N_l+N_r)>=omega_r)" << endl;
		cerr << "omega_r:\t"  << omega_r<< "\t(omega_l+ntol*(N_l+N_r)>=omega_r)" << endl;
		cerr << "ntol:\t"  << ntol<< "\t(omega_l+ntol*(N_l+N_r)>=omega_r)" << endl;
		cerr << "N_l:\t\t"       << N_l       << "\t(N_l<3)" << endl;
		cerr << "N_r:\t\t"       << N_r       << "\t(N_r<3)" << endl;
		cerr << "omegak:\t\t"   << omegak   << "\t(omegak<=omega_l) or (omegak>=omega_r)" << endl;
		cerr << "omegak_0:\t"   << omegak_0 << "\t(omegak_0<1E-15 || omegak-omegak_0<=omega_l+ntol*N_l || omegak+omegak_0>=omega_r-ntol*N_r)" << endl;
		cerr << endl;
		cerr << "Break." << endl;
		throw xBadValues();
	}
	// container for grid region parameter
	int N_i=2;
	int N_d=8;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// save grid parameter
	para_i[0]=N_l;
	para_i[1]=N_r;
	para_d[0]=omegak;
	para_d[1]=omegak_0;
	double c_1 = log((omegak-omega_l)/(omegak_0)) / (N_l-1);
	double i_1 =  log(omegak-omega_l)/c_1;
	double c_2 = log((omega_r-omegak)/(omegak_0)) / (N_r);
	double i_2 = -log(omega_r-omegak)/c_2 + N_r;
	para_d[2]=c_1;
	para_d[3]=c_2;
	// grid at i=N_l-3
	double omega_m=-exp(-c_1*(N_l-3-i_1))+omegak;
	// grid at i=N_l+Nk+3
	double omega_p=exp(c_2*(3-i_2))+omegak;
	// minima domega on the left
	double domega_min_l=exp(-c_1*(N_l-1-i_1))*(exp(c_1)-1);
	// minima domega on the right
	double domega_min_r=exp(-c_2*i_2)*(exp(c_2)-1);
	
	this->fgridRegions.push_back( gridRegion(omega_l, omegak, omega_r, omega_m, omega_p, domega_min_l, domega_min_r, N_d, N_i, para_d, para_i, "log", id));
}
// add special logarithmic grid region
void multigrid::add_sgr_log(int N_l, int N_r, double omega_l, double omega_r, double omegak, double omegak_0, string id) 
{
	double ntol=numeric_limits<double>::epsilon();
	if (N_l<3 || N_r<3 ||  omegak<=omega_l || omegak>=omega_r || omegak_0<1E-15 || omegak-omegak_0<=omega_l+ntol*N_l || omegak+omegak_0>=omega_r-ntol*N_r || omega_l+ntol*(N_l+N_r)>=omega_r)
	{
		cerr << endl;
		cerr << "Error: multigrid: add loggrid: Bad values." << endl;
		cerr << endl;
		cerr << "omega_l:\t"  << omega_l<< "\t(omega_l+ntol*(N_l+N_r)>=omega_r)" << endl;
		cerr << "omega_r:\t"  << omega_r<< "\t(omega_l+ntol*(N_l+N_r)>=omega_r)" << endl;
		cerr << "ntol:\t"  << ntol<< "\t(omega_l+ntol*(N_l+N_r)>=omega_r)" << endl;
		cerr << "N_l:\t\t"       << N_l       << "\t(N_l<3)" << endl;
		cerr << "N_r:\t\t"       << N_r       << "\t(N_r<3)" << endl;
		cerr << "omegak:\t\t"   << omegak   << "\t(omegak<=omega_l) or (omegak>=omega_r)" << endl;
		cerr << "omegak_0:\t"   << omegak_0 << "\t(omegak_0<1E-15 || omegak_0>=omegak-(omega_l+ntol*N_l) || omegak_0>=omega_r-(omegak+ntol*N_r))" << endl;
		cerr << endl;
		cerr << "Break." << endl;
		throw xBadValues();
	}
	// container for grid region parameter
	int N_i=2;
	int N_d=8;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// save grid parameter
	para_i[0]=N_l;
	para_i[1]=N_r;
	para_d[0]=omegak;
	para_d[1]=omegak_0;
	double c_1 = log((omegak-omega_l)/(omegak_0)) / (N_l-1);
	double i_1 =  log(omegak-omega_l)/c_1;
	double c_2 = log((omega_r-omegak)/(omegak_0)) / (N_r);
	double i_2 = -log(omega_r-omegak)/c_2 + N_r;
	para_d[2]=c_1;
	para_d[3]=c_2;
	// grid at i=N_l-3
	double omega_m=-exp(-c_1*(N_l-3-i_1))+omegak;
	// grid at i=N_l+Nk+3
	double omega_p=exp(c_2*(3-i_2))+omegak;
	// minima domega on the left
	double domega_min_l=exp(-c_1*(N_l-1-i_1))*(exp(c_1)-1);
	// minima domega on the right
	double domega_min_r=exp(-c_2*i_2)*(exp(c_2)-1);
	
	this->sgridRegions.push_back( gridRegion(omega_l, omegak, omega_r, omega_m, omega_p, domega_min_l, domega_min_r, N_d, N_i, para_d, para_i, "log", id));
}
// add fundamental equidistant grid region (without id)
void multigrid::add_gr_equi(int N_e, double omega_l, double omega_r, double omega_c)
{
	this->add_gr_equi(N_e, omega_l, omega_r, omega_c, "");
}
// add special equidistant grid region (without id)
void multigrid::add_sgr_equi(int N_e, double omega_l, double omega_r, double omega_c)
{
	this->add_sgr_equi(N_e, omega_l, omega_r, omega_c, "");
}
// add fundamental tangential grid region (without id)
void multigrid::add_gr_tan(int N_t, double omega_l, double omega_r, double omega_c, double c) 
{
	this->add_gr_tan(N_t, omega_l, omega_r, omega_c, c, "");
}
// add special tangential grid region (without id)
void multigrid::add_sgr_tan(int N_t, double omega_l, double omega_r, double omega_c, double c) 
{
	this->add_sgr_tan(N_t, omega_l, omega_r, omega_c, c, "");
}
// add fundamental logarithmic grid region (without id) 
void multigrid::add_gr_log(int N_l, int N_r, double omega_l, double omega_r, double omegak, double omegak_0) 
{
	this->add_gr_log(N_l, N_r, omega_l, omega_r, omegak, omegak_0, "");
}
// add special logarithmic grid region (without id)
void multigrid::add_sgr_log(int N_l, int N_r, double omega_l, double omega_r, double omegak, double omegak_0) 
{
	this->add_sgr_log(N_l, N_r, omega_l, omega_r, omegak, omegak_0, "");
}
// replace fundamental equidistant grid region 
void multigrid::replace_gr_equi(int N_e, double omega_l, double omega_r, double omega_c, string ID)
{
	if (this->gr_exists(ID))
	{
		this->rem_gr(ID);
	}
	this->add_gr_equi(N_e, omega_l, omega_r, omega_c, ID);
}
// replace special equidistant grid region 
void multigrid::replace_sgr_equi(int N_e, double omega_l, double omega_r, double omega_c, string ID)
{
	if (this->sgr_exists(ID))
	{
		this->rem_sgr(ID);
	}
	this->add_sgr_equi(N_e, omega_l, omega_r, omega_c, ID);
}
// replace fundamental tangential grid region 
void multigrid::replace_gr_tan(int N_t, double omega_l, double omega_r, double omega_c, double c, string ID) 
{
	if (this->gr_exists(ID))
	{
		this->rem_gr(ID);
	}
	this->add_gr_tan(N_t, omega_l, omega_r, omega_c, c, ID);
}
// replace special tangential grid region 
void multigrid::replace_sgr_tan(int N_t, double omega_l, double omega_r, double omega_c, double c, string ID) 
{
	if (this->sgr_exists(ID))
	{
		this->rem_sgr(ID);
	}
	this->add_sgr_tan(N_t, omega_l, omega_r, omega_c, c, ID);
}
// replace fundamental logarithmic grid region  
void multigrid::replace_gr_log(int N_l, int N_r, double omega_l, double omega_r, double omegak, double omegak_0, string ID) 
{
	if (this->gr_exists(ID))
	{
		this->rem_gr(ID);
	}
	this->add_gr_log(N_l, N_r, omega_l, omega_r, omegak, omegak_0, ID);
}
// replace special logarithmic grid region 
void multigrid::replace_sgr_log(int N_l, int N_r, double omega_l, double omega_r, double omegak, double omegak_0, string ID) 
{
	if (this->sgr_exists(ID))
	{
		this->rem_sgr(ID);
	}
	this->add_sgr_log(N_l, N_r, omega_l, omega_r, omegak, omegak_0, ID);
}

// wrapper for the equi-Gridregion functions
void multigrid::add_gr_equi(double omega_l, double omega_r, double domegap)
{
	int N=max(int((omega_r-omega_l)/domegap), 3);
	double omega_c=0.5*(omega_r+omega_l);
	this->add_gr_equi(N, omega_l, omega_r, omega_c);
}
void multigrid::add_gr_equi(double omega_l, double omega_r, double domegap, string ID)
{
	int N=max(int((omega_r-omega_l)/domegap), 3);
	double omega_c=0.5*(omega_r+omega_l);
	this->add_gr_equi(N, omega_l, omega_r, omega_c, ID);
}
void multigrid::add_sgr_equi(double omega_l, double omega_r, double domegap)
{
	int N=max(int((omega_r-omega_l)/domegap), 3);
	double omega_c=0.5*(omega_r+omega_l);
	this->add_sgr_equi(N, omega_l, omega_r, omega_c);
}
void multigrid::add_sgr_equi(double omega_l, double omega_r, double domegap, string ID)
{
	int N=max(int((omega_r-omega_l)/domegap), 3);
	double omega_c=0.5*(omega_r+omega_l);
	this->add_sgr_equi(N, omega_l, omega_r, omega_c, ID);
}
void multigrid::replace_gr_equi(double omega_l, double omega_r, double domegap, string ID)
{
	if (this->gr_exists(ID))
	{
		this->rem_gr(ID);
	}
	this->add_gr_equi(omega_l, omega_r, domegap, ID);
}
void multigrid::replace_sgr_equi(double omega_l, double omega_r, double domegap, string ID)
{
	if (this->sgr_exists(ID))
	{
		this->rem_sgr(ID);
	}
	this->add_sgr_equi(omega_l, omega_r, domegap, ID);
}

// wrapper for logarithmic grid regions
void multigrid::add_gr_log(double omega_c, double omega1, double domega_min, double domega_max)
{
	double omega0 = min(  max(omega1 *domega_min/domega_max, 1E-10),  omega1/1.1);
	int Nlog=max(int(log(omega1/omega0)*omega1/domega_max+0.5), 3);
	this->add_gr_log(Nlog, Nlog+1, omega_c-omega1, omega_c+omega1, omega_c, omega0,"");
}
void multigrid::add_gr_log(double omega_c, double omega1, double domega_min, double domega_max, string id)
{
	double omega0 = min(  max(omega1 *domega_min/domega_max, 1E-10),  omega1/1.1);
	int Nlog=max(int(log(omega1/omega0)*omega1/domega_max+0.5), 3);
	this->add_gr_log(Nlog, Nlog+1, omega_c-omega1, omega_c+omega1, omega_c, omega0, id);
}
void multigrid::add_sgr_log(double omega_c, double omega1, double domega_min, double domega_max)
{
	double omega0 = min(  max(omega1 *domega_min/domega_max, 1E-10),  omega1/1.1);
	int Nlog=max(int(log(omega1/omega0)*omega1/domega_max+0.5), 3);
	this->add_sgr_log(Nlog, Nlog+1, omega_c-omega1, omega_c+omega1, omega_c, omega0,"");
}
void multigrid::add_sgr_log(double omega_c, double omega1, double domega_min, double domega_max, string id)
{
	double omega0 = min(  max(omega1 *domega_min/domega_max, 1E-10),  omega1/1.1);
	int Nlog=max(int(log(omega1/omega0)*omega1/domega_max+0.5), 3);
	this->add_sgr_log(Nlog, Nlog+1, omega_c-omega1, omega_c+omega1, omega_c, omega0, id);
}
void multigrid::replace_gr_log(double omega_c, double omega1, double domega_min, double domega_max, string ID)
{
	if (this->gr_exists(ID))
	{
		this->rem_gr(ID);
	}
	this->add_gr_log(omega_c, omega1, domega_min, domega_max, ID);
}
void multigrid::replace_sgr_log(double omega_c, double omega1, double domega_min, double domega_max, string ID)
{
	if (this->sgr_exists(ID))
	{
		this->rem_sgr(ID);
	}
	this->add_sgr_log(omega_c, omega1, domega_min, domega_max, ID);
}
// add special grid region from existing grid region 
void multigrid::add_sgr(string ID, gridRegion & gr)
{
	if (gr.type=="equi")
	{
		this->add_sgr_equi(gr.para_i[0], gr.omega_l, gr.omega_r, gr.omega_c, ID);
	}
	else if (gr.type=="tan")
	{
		this->add_sgr_tan(gr.para_i[0], gr.omega_l, gr.omega_r, gr.omega_c, gr.para_d[1], ID);
	}
	else if (gr.type=="log")
	{
		this->add_sgr_log(gr.para_i[0], gr.para_i[1], gr.omega_l, gr.omega_r, gr.para_d[0],  gr.para_d[1], ID);
	}
}
// add fundamental grid region from existing grid region 
void multigrid::add_gr(string ID, gridRegion & gr)
{
	if (gr.type=="equi")
	{
		this->add_gr_equi(gr.para_i[0], gr.omega_l, gr.omega_r, gr.omega_c, ID);
	}
	else if (gr.type=="tan")
	{
		this->add_gr_tan(gr.para_i[0], gr.omega_l, gr.omega_r, gr.omega_c, gr.para_d[1], ID);
	}
	else if (gr.type=="log")
	{
		this->add_gr_log(gr.para_i[0], gr.para_i[1], gr.omega_l, gr.omega_r, gr.para_d[0],  gr.para_d[1], ID);
	}
}
// replace special grid region from existing grid region 
void multigrid::replace_sgr(string ID, gridRegion & gr)
{
	if (gr.type=="equi")
	{
		this->replace_sgr_equi(gr.para_i[0], gr.omega_l, gr.omega_r, gr.omega_c, ID);
	}
	else if (gr.type=="tan")
	{
		this->replace_sgr_tan(gr.para_i[0], gr.omega_l, gr.omega_r, gr.omega_c, gr.para_d[1], ID);
	}
	else if (gr.type=="log")
	{
		this->replace_sgr_log(gr.para_i[0], gr.para_i[1], gr.omega_l, gr.omega_r, gr.para_d[0],  gr.para_d[1], ID);
	}
}
// replace fundamental grid region from existing grid region 
void multigrid::replace_gr(string ID, gridRegion & gr)
{
	if (gr.type=="equi")
	{
		this->replace_gr_equi(gr.para_i[0], gr.omega_l, gr.omega_r, gr.omega_c, ID);
	}
	else if (gr.type=="tan")
	{
		this->replace_gr_tan(gr.para_i[0], gr.omega_l, gr.omega_r, gr.omega_c, gr.para_d[1], ID);
	}
	else if (gr.type=="log")
	{
		this->replace_gr_log(gr.para_i[0], gr.para_i[1], gr.omega_l, gr.omega_r, gr.para_d[0],  gr.para_d[1], ID);
	}
}

// add fundamental equidistant sub grid
void multigrid::addFSubgrid_equi(int N_e, double omega_l, double omega_r) 
{
	// container for grid region parameter
	int N_i=1;
	int N_d=1;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// create equidistant grid
	equigrid egrid(N_e, omega_l, omega_r);
	// save grid parameter
	para_i[0]=egrid.M;
	para_d[0]=egrid.delta_omega;
	int i_l;
	int i_r;
	int N_replaced;
	
	// machine precision
	double ntol=10*numeric_limits<double>::epsilon();

	// if this is the first grid region
	if (fsubgrids.size()==0)
	{
		this->omega_min=egrid.omega_min;
		this->omega_max=egrid.omega_max;
		this->M=egrid.M;
		this->omega=vector<double> (M+1);
		this->domega=vector<double> (M+1);
		for (int j=0; j!=M+1; j++)
		{
			omega[j]=egrid.omega[j];
			domega[j]=egrid.domega[j];
		}
		i_l=0;
		i_r=M+1;
		N_replaced=M;
	}
	else
	{
		if (omega_l<this->omega_min || omega_l>=this->omega_max || omega_r<=this->omega_min || omega_r>this->omega_max )
		{
			cerr << "Error: Multigrid: addFSubgrid_equi: Boundaries of subgrid exceeded boundaries of base grid. Break." << endl;
			cerr << "omega_l: " << omega_l << "  omega_r: " << omega_r << endl;
			cerr << "omega_m: " << omega_min << "  omega_m: " << omega_max << endl;
			throw 1;
		}
		// check if there is no intersection between the subgrids
		for (int n=1; n!=fsubgrids.size(); n++)
		{
			if ((omega_l<fsubgrids[n].omega_r && fsubgrids[n].omega_r<=omega_r) || (omega_l<=fsubgrids[n].omega_l && fsubgrids[n].omega_l<omega_r) || (fsubgrids[n].omega_l<omega_r && omega_r<=fsubgrids[n].omega_r) || (fsubgrids[n].omega_l<=omega_l && omega_l<fsubgrids[n].omega_r))
			{
				cerr << "Error: Multigrid: Intersection of subgrids is not allowed. Break." << endl;
				cerr << "While trying to create subgrid between " << scientific << setprecision(15) << omega_l << " and " << omega_r << endl;
				cerr << "Subgrid " << n << " is already between " << fsubgrids[n].omega_l << " and " << fsubgrids[n].omega_r << endl;
				throw 1;
			}
		}
		// find index of left and right boundary of new grid region
		i_l=0;
		while (i_l<=this->M && this->omega[i_l]<=omega_l-ntol*max(fabs(omega_l),1.0))
		{
			i_l++;
		}
		i_r=this->M;
		while (i_r>=0 && this->omega[i_r]>=omega_r+ntol*max(fabs(omega_r),1.0))
		{
			i_r--;
		}

		// erase region between i_l and i_r in old grid
		if (i_l<=i_r)
		{
			this->omega.erase(this->omega.begin() + i_l, this->omega.begin() + i_r+1);
		}
		// put in new grid region 
		this->omega.insert(this->omega.begin() + i_l, egrid.omega.begin(), egrid.omega.end());
		this->M=omega.size()-1;

		// number of inserted - replaced grid points
		N_replaced=egrid.M+1-max((i_r-i_l)+1,0);
		// adjust boundaries indices for other grid regions
		for (int n=0; n!=fsubgrids.size(); n++)
		{
			// if the grid region is above the current one
			if (omega_r<=fsubgrids[n].omega_l)
			{
				fsubgrids[n].i_l+=N_replaced;
				fsubgrids[n].i_r+=N_replaced;
			}
		}
	}
	
	fsubgrids.push_back( subgrid(omega_l, omega_r, i_l, i_r, N_replaced, N_d, N_i, para_d, para_i, "equi"));	
}
// add fundamental tangential sub grid
void multigrid::addFSubgrid_tan(int N_t, double omega_l, double omega_r, double omega_c, double c) 
{
	// container for grid region parameter
	int N_i=1;
	int N_d=5;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// create tangential grid
	tangrid tgrid(N_t, omega_l, omega_r, omega_c, c);

	// save grid parameter
	para_i[0]=tgrid.M;
	para_d[0]=tgrid.omega_c;
	para_d[1]=tgrid.c;
	para_d[2]=tgrid.u_0;
	para_d[3]=tgrid.u_1;
	para_d[4]=tgrid.du;

	int i_l;
	int i_r;
	int N_replaced;
	
	// machine precision
	double ntol=10*numeric_limits<double>::epsilon();

	// if this is the first grid region
	if (fsubgrids.size()==0)
	{
		this->omega_min=tgrid.omega_min;
		this->omega_max=tgrid.omega_max;
		this->M=tgrid.M;
		this->omega=vector<double> (M+1);
		this->domega=vector<double> (M+1);
		for (int j=0; j!=M+1; j++)
		{
			omega[j]=tgrid.omega[j];
			domega[j]=tgrid.domega[j];
		}
		i_l=0;
		i_r=M+1;
		N_replaced=M;
	}
	else
	{
		if (omega_l<this->omega_min || omega_l>=this->omega_max || omega_r<=this->omega_min || omega_r>this->omega_max )
		{
			cerr << "Error: Multigrid: addFSubgrid_tan: Boundaries of subgrid exceeded boundaries of base grid. Break." << endl;
			throw 1;
		}
		// check if there is no intersection between the fsubgrids
		for (int n=1; n!=fsubgrids.size(); n++)
		{
			if ((omega_l<fsubgrids[n].omega_r && fsubgrids[n].omega_r<=omega_r) || (omega_l<=fsubgrids[n].omega_l && fsubgrids[n].omega_l<omega_r) || (fsubgrids[n].omega_l<omega_r && omega_r<=fsubgrids[n].omega_r) || (fsubgrids[n].omega_l<=omega_l && omega_l<fsubgrids[n].omega_r))
			{
				cerr << "Error: Multigrid: Intersection of subgrids is not allowed. Break." << endl;
				cerr << "While trying to create subgrid between " << scientific << setprecision(15) << omega_l << " and " << omega_r << endl;
				cerr << "Subgrid " << n << " is already between " << fsubgrids[n].omega_l << " and " << fsubgrids[n].omega_r << endl;
				throw 1;
			}
		}
		// find index of left and right boundary of new grid region
		i_l=0;
		while (i_l<=this->M && this->omega[i_l]<=omega_l-ntol*max(fabs(omega_l),1.0))
		{
			i_l++;
		}
		i_r=this->M;
		while (i_r>=0 && this->omega[i_r]>=omega_r+ntol*max(fabs(omega_r),1.0))
		{
			i_r--;
		}

		// erase region between i_l and i_r in old grid
		if (i_l<=i_r)
		{
			this->omega.erase(this->omega.begin() + i_l, this->omega.begin() + i_r+1);
		}
		// put in new grid region 
		this->omega.insert(this->omega.begin() + i_l, tgrid.omega.begin(), tgrid.omega.end());
		this->M=omega.size()-1;

		// number of inserted - replaced grid points
		N_replaced=tgrid.M+1-max((i_r-i_l)+1,0);
		// adjust boundaries indices for other grid regions
		for (int n=0; n!=fsubgrids.size(); n++)
		{
			// if the grid region is above the current one
			if (omega_r<fsubgrids[n].omega_l)
			{
				fsubgrids[n].i_l+=N_replaced;
				fsubgrids[n].i_r+=N_replaced;
			}
		}
	}
	fsubgrids.push_back( subgrid(omega_l, omega_r, i_l, i_r, N_replaced, N_d, N_i, para_d, para_i, "tan"));	
}

// add fundamental logarithmic sub grid
void multigrid::addFSubgrid_log(int N_l, int N_r, double omega_l, double omega_r, double omegak, double omegak_0) 
{
	// container for grid region parameter
	int N_i=4;
	int N_d=10;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// create tangential grid
	loggrid lgrid(N_l, N_r, omega_l, omega_r, omegak, omegak_0);
	// save grid parameter
	para_i[0]=lgrid.M;
	para_i[1]=lgrid.N_l;
	para_i[2]=lgrid.Nk;
	para_i[3]=lgrid.N_r;
	para_d[0]=lgrid.omegak;
	para_d[1]=lgrid.omegak_0;
	para_d[2]=lgrid.domegak_min;
	para_d[3]=lgrid.domegak;
	para_d[4]=lgrid.omegak_0m;
	para_d[5]=lgrid.omegak_0p;
	para_d[6]=lgrid.c_1;
	para_d[7]=lgrid.i_1;
	para_d[8]=lgrid.c_2;
	para_d[9]=lgrid.i_2;

	int i_l;
	int i_r;
	int N_replaced;
	
	// machine precision
	double ntol=10*numeric_limits<double>::epsilon();

	// if this is the first grid region
	if (fsubgrids.size()==0)
	{
		this->omega_min=lgrid.omega_min;
		this->omega_max=lgrid.omega_max;
		this->M=lgrid.M;
		this->omega=vector<double> (M+1);
		this->domega=vector<double> (M+1);
		for (int j=0; j!=M+1; j++)
		{
			omega[j]=lgrid.omega[j];
			domega[j]=lgrid.domega[j];
		}
		i_l=0;
		i_r=M+1;
		N_replaced=M;
	}
	else
	{
		if (omega_l<this->omega_min || omega_l>=this->omega_max || omega_r<=this->omega_min || omega_r>this->omega_max )
		{
			cerr << "Error: Multigrid: addFSubgrid_log: Boundaries of subgrid exceeded boundaries of base grid. Break." << endl;
			throw 1;
		}
		// check if there is no intersection between the subgrids
		for (int n=1; n!=fsubgrids.size(); n++)
		{
			if ((omega_l<fsubgrids[n].omega_r && fsubgrids[n].omega_r<=omega_r) || (omega_l<=fsubgrids[n].omega_l && fsubgrids[n].omega_l<omega_r) || (fsubgrids[n].omega_l<omega_r && omega_r<=fsubgrids[n].omega_r) || (fsubgrids[n].omega_l<=omega_l && omega_l<fsubgrids[n].omega_r))
			{
				cerr << "Error: Multigrid: Intersection of subgrids is not allowed. Break." << endl;
				cerr << "While trying to create subgrid between " << scientific << setprecision(15) << omega_l << " and " << omega_r << endl;
				cerr << "Subgrid " << n << " is already between " << fsubgrids[n].omega_l << " and " << fsubgrids[n].omega_r << endl;
				throw 1;
			}
		}
		// find index of left and right boundary of new grid region
		i_l=0;
		while (i_l<=this->M && this->omega[i_l]<=omega_l-ntol*max(fabs(omega_l),1.0))
		{
			i_l++;
		}
		i_r=this->M;
		while (i_r>=0 && this->omega[i_r]>=omega_r+ntol*max(fabs(omega_r),1.0))
		{
			i_r--;
		}

		// erase region between i_l and i_r in old grid
		if (i_l<=i_r)
		{
			this->omega.erase(this->omega.begin() + i_l, this->omega.begin() + i_r+1);
		}
		// put in new grid region 
		this->omega.insert(this->omega.begin() + i_l, lgrid.omega.begin(), lgrid.omega.end());
		this->M=omega.size()-1;

		// number of inserted - replaced grid points
		N_replaced=lgrid.M+1-max((i_r-i_l)+1,0);
		// adjust boundaries indices for other grid regions
		for (int n=0; n!=fsubgrids.size(); n++)
		{
			// if the grid region is above the current one
			if (omega_r<fsubgrids[n].omega_l)
			{
				fsubgrids[n].i_l+=N_replaced;
				fsubgrids[n].i_r+=N_replaced;
			}
		}
	}
	fsubgrids.push_back( subgrid(omega_l, omega_r, i_l, i_r, N_replaced, N_d, N_i, para_d, para_i, "log"));
}
// add special equidistant sub grid
void multigrid::addSSubgrid_equi(int N_e, double omega_l, double omega_r) 
{
	// container for grid region parameter
	int N_i=1;
	int N_d=1;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// create equidistant grid
	equigrid egrid(N_e, omega_l, omega_r);
	// save grid parameter
	para_i[0]=egrid.M;
	para_d[0]=egrid.delta_omega;
	int i_l;
	int i_r;
	int N_replaced;
	
	// machine precision
	double ntol=10*numeric_limits<double>::epsilon();

	if (fsubgrids.size()==0)
	{
		cerr << "Error: Multigrid: It is not allowed to create a special grid if there is no fundamental grid. Break." << endl;
		throw 1;
	}
	if (omega_l<this->omega_min || omega_l>=this->omega_max || omega_r<=this->omega_min || omega_r>this->omega_max )
	{
		cerr << "Error: Multigrid: addSSubgrid_equi: Boundaries of subgrid exceeded boundaries of base grid. Break." << endl;
		throw 1;
	}
	// check if there is no intersection between the subgrids
	for (int n=0; n<ssubgrids.size(); n++)
	{
		if ((omega_l<ssubgrids[n].omega_r && ssubgrids[n].omega_r<=omega_r) || (omega_l<=ssubgrids[n].omega_l && ssubgrids[n].omega_l<omega_r) || (ssubgrids[n].omega_l<omega_r && omega_r<=ssubgrids[n].omega_r) || (ssubgrids[n].omega_l<=omega_l && omega_l<ssubgrids[n].omega_r))
		{
			cerr << "Error: Multigrid: Intersection of subgrids is not allowed. Break." << endl;
			cerr << "While trying to create subgrid between " << scientific << setprecision(15) << omega_l << " and " << omega_r << endl;
			cerr << "Subgrid " << n << " is already between " << ssubgrids[n].omega_l << " and " << ssubgrids[n].omega_r << endl;
			throw 1;
		}
	}
	// find index of left and right boundary of new grid region
	i_l=0;
	while (i_l<=this->M && this->omega[i_l]<=omega_l-ntol*max(fabs(omega_l),1.0))
	{
		i_l++;
	}
	i_r=this->M;
	while (i_r>=0 && this->omega[i_r]>=omega_r+ntol*max(fabs(omega_r),1.0))
	{
		i_r--;
	}

	// erase region between i_l and i_r in old grid
	if (i_l<=i_r)
	{
		this->omega.erase(this->omega.begin() + i_l, this->omega.begin() + i_r+1);
	}
	// put in new grid region 
	this->omega.insert(this->omega.begin() + i_l, egrid.omega.begin(), egrid.omega.end());
	this->M=omega.size()-1;

	// number of inserted - replaced grid points
	N_replaced=egrid.M+1-max((i_r-i_l)+1,0);
	// adjust boundaries indices for other grid regions
	for (int n=0; n!=ssubgrids.size(); n++)
	{
		// if the grid region is above the current one
		if (omega_r<=ssubgrids[n].omega_l)
		{
			ssubgrids[n].i_l+=N_replaced;
			ssubgrids[n].i_r+=N_replaced;
		}
	}
	ssubgrids.push_back( subgrid(omega_l, omega_r, i_l, i_r, N_replaced, N_d, N_i, para_d, para_i, "equi"));	
}
// add special tangential sub grid
void multigrid::addSSubgrid_tan(int N_t, double omega_l, double omega_r, double omega_c, double c) 
{
	// container for grid region parameter
	int N_i=1;
	int N_d=5;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// create tangential grid
	tangrid tgrid(N_t, omega_l, omega_r, omega_c, c);

	// save grid parameter
	para_i[0]=tgrid.M;
	para_d[0]=tgrid.omega_c;
	para_d[1]=tgrid.c;
	para_d[2]=tgrid.u_0;
	para_d[3]=tgrid.u_1;
	para_d[4]=tgrid.du;

	int i_l;
	int i_r;
	int N_replaced;
	
	// machine precision
	double ntol=10*numeric_limits<double>::epsilon();

	if (fsubgrids.size()==0)
	{
		cerr << "Error: Multigrid: It is not allowed to create a special grid if there is no fundamental grid. Break." << endl;
		throw 1;
	}
	if (omega_l<this->omega_min || omega_l>=this->omega_max || omega_r<=this->omega_min || omega_r>this->omega_max )
	{
		cerr << "Error: Multigrid: addSSubgrid_tan: Boundaries of subgrid exceeded boundaries of base grid. Break." << endl;
		throw 1;
	}
	// check if there is no intersection between the ssubgrids
	for (int n=0; n<ssubgrids.size(); n++)
	{
		if ((omega_l<ssubgrids[n].omega_r && ssubgrids[n].omega_r<=omega_r) || (omega_l<=ssubgrids[n].omega_l && ssubgrids[n].omega_l<omega_r) || (ssubgrids[n].omega_l<omega_r && omega_r<=ssubgrids[n].omega_r) || (ssubgrids[n].omega_l<=omega_l && omega_l<ssubgrids[n].omega_r))
		{
			cerr << "Error: Multigrid: Intersection of subgrids is not allowed. Break." << endl;
			cerr << "While trying to create subgrid between " << scientific << setprecision(15) << omega_l << " and " << omega_r << endl;
			cerr << "Subgrid " << n << " is already between " << ssubgrids[n].omega_l << " and " << ssubgrids[n].omega_r << endl;
			throw 1;
		}
	}
	// find index of left and right boundary of new grid region
	i_l=0;
	while (i_l<=this->M && this->omega[i_l]<=omega_l-ntol*max(fabs(omega_l),1.0))
	{
		i_l++;
	}
	i_r=this->M;
	while (i_r>=0 && this->omega[i_r]>=omega_r+ntol*max(fabs(omega_r),1.0))
	{
		i_r--;
	}

	// erase region between i_l and i_r in old grid
	if (i_l<=i_r)
	{
		this->omega.erase(this->omega.begin() + i_l, this->omega.begin() + i_r+1);
	}
	// put in new grid region 
	this->omega.insert(this->omega.begin() + i_l, tgrid.omega.begin(), tgrid.omega.end());
	this->M=omega.size()-1;

	// number of inserted - replaced grid points
	N_replaced=tgrid.M+1-max((i_r-i_l)+1,0);
	// adjust boundaries indices for other grid regions
	for (int n=0; n!=ssubgrids.size(); n++)
	{
		// if the grid region is above the current one
		if (omega_r<ssubgrids[n].omega_l)
		{
			ssubgrids[n].i_l+=N_replaced;
			ssubgrids[n].i_r+=N_replaced;
		}
	}
	ssubgrids.push_back( subgrid(omega_l, omega_r, i_l, i_r, N_replaced, N_d, N_i, para_d, para_i, "tan"));	
}

// add special logarithmic sub grid
void multigrid::addSSubgrid_log(int N_l, int N_r, double omega_l, double omega_r, double omegak, double omegak_0) 
{
	// container for grid region parameter
	int N_i=4;
	int N_d=10;
	vector<int> para_i(N_i);
	vector<double> para_d(N_d);
	// create tangential grid
	loggrid lgrid(N_l, N_r, omega_l, omega_r, omegak, omegak_0);
	// save grid parameter
	para_i[0]=lgrid.M;
	para_i[1]=lgrid.N_l;
	para_i[2]=lgrid.Nk;
	para_i[3]=lgrid.N_r;
	para_d[0]=lgrid.omegak;
	para_d[1]=lgrid.omegak_0;
	para_d[2]=lgrid.domegak_min;
	para_d[3]=lgrid.domegak;
	para_d[4]=lgrid.omegak_0m;
	para_d[5]=lgrid.omegak_0p;
	para_d[6]=lgrid.c_1;
	para_d[7]=lgrid.i_1;
	para_d[8]=lgrid.c_2;
	para_d[9]=lgrid.i_2;

	int i_l;
	int i_r;
	int N_replaced;
	
	// machine precision
	double ntol=10*numeric_limits<double>::epsilon();
	if (fsubgrids.size()==0)
	{
		cerr << "Error: Multigrid: It is not allowed to create a special grid if there is no fundamental grid. Break." << endl;
		throw 1;
	}
	if (omega_l<this->omega_min || omega_l>=this->omega_max || omega_r<=this->omega_min || omega_r>this->omega_max )
	{
		cerr << "Error: Multigrid: addSSubgrid_log: Boundaries of subgrid exceeded boundaries of base grid. Break." << endl;
		throw 1;
	}
	// check if there is no intersection between the subgrids
	for (int n=0; n<ssubgrids.size(); n++)
	{
		if ((omega_l<ssubgrids[n].omega_r && ssubgrids[n].omega_r<=omega_r) || (omega_l<=ssubgrids[n].omega_l && ssubgrids[n].omega_l<omega_r) || (ssubgrids[n].omega_l<omega_r && omega_r<=ssubgrids[n].omega_r) || (ssubgrids[n].omega_l<=omega_l && omega_l<ssubgrids[n].omega_r))
		{
			cerr << "Error: Multigrid: Intersection of subgrids is not allowed. Break." << endl;
			cerr << "While trying to create subgrid between " << scientific << setprecision(15) << omega_l << " and " << omega_r << endl;
			cerr << "Subgrid " << n << " is already between " << ssubgrids[n].omega_l << " and " << ssubgrids[n].omega_r << endl;
			throw 1;
		}
	}
	// find index of left and right boundary of new grid region
	i_l=0;
	while (i_l<=this->M && this->omega[i_l]<=omega_l-ntol*max(fabs(omega_l),1.0))
	{
		i_l++;
	}
	i_r=this->M;
	while (i_r>=0 && this->omega[i_r]>=omega_r+ntol*max(fabs(omega_r),1.0))
	{
		i_r--;
	}

	// erase region between i_l and i_r in old grid
	if (i_l<=i_r)
	{
		this->omega.erase(this->omega.begin() + i_l, this->omega.begin() + i_r+1);
	}
	// put in new grid region 
	this->omega.insert(this->omega.begin() + i_l, lgrid.omega.begin(), lgrid.omega.end());
	this->M=omega.size()-1;

	// number of inserted - replaced grid points
	N_replaced=lgrid.M+1-max((i_r-i_l)+1,0);
	// adjust boundaries indices for other grid regions
	for (int n=0; n!=ssubgrids.size(); n++)
	{
		// if the grid region is above the current one
		if (omega_r<ssubgrids[n].omega_l)
		{
			ssubgrids[n].i_l+=N_replaced;
			ssubgrids[n].i_r+=N_replaced;
		}
	}
	ssubgrids.push_back( subgrid(omega_l, omega_r, i_l, i_r, N_replaced, N_d, N_i, para_d, para_i, "log"));
}


// get intersection points of the subgrids depending on the type of subgrid
bool getIntersectionPoint(double & omegas, gridRegion & grl, gridRegion & grr)
{
	// loggrid:
	// para_i[0] is N_l
	// para_i[1] is N_r
	// para_d[0] is omegak
	// para_d[1] is omegak_0
	// para_d[2] is c_1
	// para_d[3] is c_2

	// tangrid:
	// para_i[0] is N_t
	// para_d[0] is omega_c
	// para_d[1] is c
	// para_d[2] is du 

	// equigrid:
	// para_i[0] is N_e
	// para_d[0] is delta_omega


	if (grl.type=="log" && grr.type=="log")
	{
		omegas=(grl.para_d[3]*grl.para_d[0]+grr.para_d[2]*grr.para_d[0])/(grl.para_d[3]+grr.para_d[2]);
		return true;
	}
	else if (grl.type=="log" && grr.type=="tan")
	{
		double phalf=grr.para_d[0] + (grl.para_d[3]*grr.para_d[1])/(2*grr.para_d[2]);
		double q=pow(grr.para_d[0], 2)+(grl.para_d[3]*grr.para_d[1]*grl.para_d[0])/grr.para_d[2] + pow(grr.para_d[1],2);
		double discriminant=pow(phalf,2)-q;
		// in no intersection exists return error code
		if (discriminant<0.0)
		{
			return false;
		}
		else
		{
			omegas=phalf - sqrt(discriminant);
			return true;
		}
	}
	else if (grl.type=="tan" && grr.type=="log")
	{
		double phalf=grl.para_d[0] - (grr.para_d[2]*grl.para_d[1])/(2*grl.para_d[2]);
		double q=pow(grl.para_d[0],2)-(grr.para_d[2]*grl.para_d[1]*grr.para_d[0])/grl.para_d[2] + pow(grl.para_d[1],2);
		double discriminant=pow(phalf,2)-q;
		// in no intersection exists return error code
		if (discriminant<0.0)
		{
			return false;
		}
		else 
		{
			omegas=phalf + sqrt(discriminant);
			return true;
		}
	}
	else if (grl.type=="tan" && grr.type=="tan")
	{
		// get intersection point
		double phalf=(grl.para_d[2]*grl.para_d[0]*grr.para_d[1]-grr.para_d[2]*grr.para_d[0]*grl.para_d[1])/(grl.para_d[2]*grr.para_d[1]-grr.para_d[2]*grl.para_d[1]);
		double q=(grr.para_d[1]*grl.para_d[2]*pow(grl.para_d[0],2)-grl.para_d[1]*grr.para_d[2]*pow(grr.para_d[0],2)+grl.para_d[2]*grr.para_d[1]*pow(grl.para_d[1],2)-grr.para_d[2]*grl.para_d[1]*pow(grr.para_d[1],2))/(grl.para_d[2]*grr.para_d[1]-grr.para_d[2]*grl.para_d[1]);
		double discriminant=pow(phalf,2)-q;
		// in no intersection exists retrun error code
		if (discriminant<0.0)
		{
			//cout << "\t\t# " << discriminant << endl;
			return false;
		}
		// intersection point has to be inside the grids
		else
		{
			double omegas_p=phalf + sqrt(discriminant);
			double omegas_m=phalf - sqrt(discriminant);
			//cout << "\t\t# " << omegas_p << "\t" << omegas_p << "\t" << grl.para_d[0] << "\t" << grr.para_d[0] << endl;
			if (omegas_p<=grr.para_d[0] && omegas_p>=grl.para_d[0])
			{
				omegas=omegas_p;
				return true;
			}
			else if (omegas_m<=grr.para_d[0] && omegas_m>=grl.para_d[0])
			{
				omegas=omegas_m;
				return true;
			}
			else
			{
				return false;
			}	
		}
	}
	else if (grl.type=="log" && grr.type=="equi")
	{
		omegas=grl.para_d[0]+grr.para_d[0]/grl.para_d[3];
		return true;
	}
	else if (grl.type=="equi" && grr.type=="log")
	{
		omegas=grr.para_d[0]-grl.para_d[0]/grr.para_d[2];
		return true;
	}
	else if (grl.type=="tan" && grr.type=="equi")
	{
		double discriminant=grr.para_d[0]/(grl.para_d[1]*grl.para_d[2]) -1.0;
		if (discriminant<0.0)
		{
			return false;
		}
		else
		{
			omegas=grl.para_d[0]+grl.para_d[1]*sqrt(discriminant);
			return true;
		}
	}
	else if (grl.type=="equi" && grr.type=="tan")
	{
		double discriminant=grl.para_d[0]/(grr.para_d[1]*grr.para_d[2]) -1.0;
		if (discriminant<0.0)
		{
			return false;
		}
		else
		{
			omegas=grr.para_d[0]-grr.para_d[1]*sqrt(discriminant);
			return true;
		}
	}
	else if (grl.type=="equi" && grr.type=="equi")
	{
		if (grr.para_d[0]<grl.para_d[0])
		{
			omegas=grr.omega_l;
		}
		else if (grr.para_d[0]>grl.para_d[0])
		{
			omegas=grl.omega_r;
		}
		else
		{
			omegas=0.5*(grl.omega_r+grr.omega_l);
		}
		return true;
	}
	else
	{
		cerr << "Error: multigrid: getIntersectionPoint routine. Break" << endl;
		throw 1;
	}	
}


// get intersection points of the subgrids depending on the type of subgrid
void adjustGridParameters(double & omegas, gridRegion & grl, gridRegion & grr)
{
	// loggrid:
	// para_i[0] is N_l
	// para_i[1] is N_r
	// para_d[0] is omegak
	// para_d[1] is omegak_0
	// para_d[2] is c_1
	// para_d[3] is c_2

	// tangrid:
	// para_i[0] is N_t
	// para_d[0] is omega_c
	// para_d[1] is c
	// para_d[2] is du 

	if (grl.type=="log")
	{
		grl.para_i[1]=max(int(log((omegas-grl.para_d[0])/(grl.para_d[1]))/grl.para_d[3] + 0.5), 3);
		grl.omega_r=omegas;
	}
	else if (grl.type=="tan")	
	{
		grl.para_i[0]=max(int((atan2(omegas-grl.para_d[0], grl.para_d[1])-atan2(grl.omega_l-grl.para_d[0], grl.para_d[1]))/grl.para_d[2]+0.5), 3);
		grl.omega_r=omegas;
	}
	else if (grl.type=="equi")
	{
		grl.para_i[0]=max(int((omegas-grl.omega_l)/(grl.omega_r-grl.omega_l)*grl.para_i[0]+0.5), 3);
		grl.omega_r=omegas;
	}
	else
	{
		cerr << "Error: multigrid: adjustGridParameters routine. Break" << endl;
		throw 1;
	}	

	if (grr.type=="log")
	{
		grr.para_i[0]=max(int(log((grr.para_d[0]-omegas)/grr.para_d[1])/grr.para_d[2] + 0.5)+1, 3);
		grr.omega_l=omegas;
	}
	else if (grr.type=="tan")
	{
		grr.para_i[0]=max(int((atan2(grr.omega_r-grr.para_d[0], grr.para_d[1])-atan2(omegas-grr.para_d[0], grr.para_d[1]))/grr.para_d[2]+0.5), 3);
		grr.omega_l=omegas;
	}
	else if (grr.type=="equi")
	{
		grr.para_i[0]=max(int((grr.omega_r-omegas)/(grr.omega_r-grr.omega_l)*grr.para_i[0]+0.5), 3);
		grr.omega_l=omegas;
	
	}
	else
	{
		cerr << "Error: multigrid: adjustGridParameters routine. Break" << endl;
		throw 1;
	}
}

// get intersection points of the subgrids depending on the type of subgrid and decide if a gridregion is skipped
int getIntersection(gridRegion & grl, gridRegion & grr)
{
	double omegas;
	bool spExists;
	spExists=getIntersectionPoint(omegas, grl, grr);
	if (!spExists)
	{
		// if no intersection point exists, keep the grid with higher resolution
		if (grl.domega_min_r>grr.domega_min_l)
		{
			// but outside of the right grid's boundaries there may be space for an implementation of the left grid
			if (grl.omega_p<=grr.omega_l)
			{
				//cout << "\t + case E-1" << endl;
				omegas=grr.omega_l;
			}
			// if there is not enough space, skip the left grid
			else
			{
				//cout << "\t + case E-2" << endl;
				return 1;
			}
		}
		else
		{
			// but outside of the left grid's boundaries there may be space for an implementation of the right grid
			if (grr.omega_m>=grl.omega_r)
			{
				//cout << "\t + case F-1" << endl;
				omegas=grl.omega_r;
			}
			// if there is not enough space, skip the right grid
			else
			{
				//cout << "\t + case F-2" << endl;
				return 2;
			}
		}
	}
	// case B, maximal resolution of the left grid is to small -> skip it
	if (omegas<grl.omega_p && omegas >= grr.omega_l)
	{
		// case D: intersection point is inside BOTH puffer-regions, skip the grid with the lower maximal resolution
		if (omegas>grr.omega_m && omegas <= grl.omega_r)
		{
			if(grr.domega_min_l<grl.domega_min_r)
			{
				//cout << "\t + case D-1" << endl;
				return 1;
			}	
			else
			{
				//cout << "\t + case D-2" << endl;
				return 2;
			}	
		}
		else
		{
			//cout << "\t + case B" << endl;
			return 1;
		}
	}
	// case E, there is no intersection point because left grid resolution is to small,
	else if (omegas<grr.omega_l)
	{
		// but outside of the right grid's boundaries there may be space for an implementation of the left grid
		if (grl.omega_p<=grr.omega_l)
		{
			//cout << "\t + case E-1" << endl;
			omegas=grr.omega_l;
		}
		// if there is not enough space, skip the left grid
		else
		{
			//cout << "\t + case E-2" << endl;
			return 1;
		}
	}
	// case C, maximal resolution of the right grid is to small -> skip it
	else if (omegas>grr.omega_m && omegas <= grl.omega_r)
	{
		//cout << "\t + case C" << endl;
		return 2;
	}
	// case F, there is no intersection point because right grid resolution is to small,
	else if (omegas>grl.omega_r)
	{
		// but outside of the left grid's boundaries there may be space for an implementation of the right grid
		if (grr.omega_m>=grl.omega_r)
		{
			//cout << "\t + case F-1" << endl;
			omegas=grl.omega_r;
		}
		// if there is not enough space, skip the right grid
		else
		{
			//cout << "\t + case F-2" << endl;
			return 2;
		}
	}
	////cout << "\t + Intersection point: " << scientific << setprecision(15) << omegas << endl;
	adjustGridParameters(omegas, grl, grr);
	return 0;
}

// cut grid regions from left or right (e.g. if they are out ouf bounds)
void cutGridRegion_right(double & omegas, gridRegion & grl)
{
	// loggrid:
	// para_i[0] is N_l
	// para_i[1] is N_r
	// para_d[0] is omegak
	// para_d[1] is omegak_0
	// para_d[2] is c_1
	// para_d[3] is c_2

	// tangrid:
	// para_i[0] is N_t
	// para_d[0] is omega_c
	// para_d[1] is c
	// para_d[2] is du 

	if (grl.type=="log")
	{
		grl.para_i[1]=max(int(log((omegas-grl.para_d[0])/(grl.para_d[1]))/grl.para_d[3] + 0.5), 3);
		grl.omega_r=omegas;
	}
	else if (grl.type=="tan")	
	{
		grl.para_i[0]=max(int((atan2(omegas-grl.para_d[0], grl.para_d[1])-atan2(grl.omega_l-grl.para_d[0], grl.para_d[1]))/grl.para_d[2]+0.5), 3);
		grl.omega_r=omegas;
	}
	else if (grl.type=="equi")
	{
		grl.para_i[0]=max(int((omegas-grl.omega_l)/(grl.omega_r-grl.omega_l)*grl.para_i[0]+0.5), 3);
		grl.omega_r=omegas;
	}
	else
	{
		cerr << "Error: multigrid: cutGridRegion_right routine. Break" << endl;
		throw 1;
	}	
}
void cutGridRegion_left(double & omegas, gridRegion & grr)
{
	// loggrid:
	// para_i[0] is N_l
	// para_i[1] is N_r
	// para_d[0] is omegak
	// para_d[1] is omegak_0
	// para_d[2] is c_1
	// para_d[3] is c_2

	// tangrid:
	// para_i[0] is N_t
	// para_d[0] is omega_c
	// para_d[1] is c
	// para_d[2] is du 

	if (grr.type=="log")
	{
		grr.para_i[0]=max(int(log((grr.para_d[0]-omegas)/grr.para_d[1])/grr.para_d[2] + 0.5)+1, 3);
		grr.omega_l=omegas;
	}
	else if (grr.type=="tan")
	{
		grr.para_i[0]=max(int((atan2(grr.omega_r-grr.para_d[0], grr.para_d[1])-atan2(omegas-grr.para_d[0], grr.para_d[1]))/grr.para_d[2]+0.5), 3);
		grr.omega_l=omegas;
	}
	else if (grr.type=="equi")
	{
		grr.para_i[0]=max(int((grr.omega_r-omegas)/(grr.omega_r-grr.omega_l)*grr.para_i[0]+0.5), 3);
		grr.omega_l=omegas;
	
	}
	else
	{
		cerr << "Error: multigrid: cutGridRegion_left routine. Break" << endl;
		throw 1;
	}	
}
// check if grid regions shall be skipped if they are out of bounds
bool skipGridRegion(gridRegion & gr, gridRegion & fgr, double epsilon)
{
	if ( gr.omega_r<=fgr.omega_l+epsilon || gr.omega_l>=fgr.omega_r-epsilon)
	{
		return true;
	}
	else if ( gr.omega_c<=fgr.omega_l+epsilon || gr.omega_c>=fgr.omega_r-epsilon)
	{
		return true;
	}
	else if (gr.type=="log")
	{
		if ( gr.omega_c+gr.para_d[1]>=fgr.omega_r-epsilon || gr.omega_c-gr.para_d[1]<=fgr.omega_l+epsilon )
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}
// createMultigrid
void multigrid::create() 
{
	// ***************************************************************************
	// check if grid regions intersect with boundaries of the first fundamental grid
	// -> skip or cut them
	// ***************************************************************************
	// get length of first fundamental grid region
	double ntol=numeric_limits<double>::epsilon();
	int Mfirst=0;
	for (int k=0; k<fgridRegions[0].para_i.size(); k++)
	{
		Mfirst+=fgridRegions[0].para_i[k];
	}
	vector<int> eraseCandidates;
	// find out which regions shall be skipped (or cutted)
	for (int n=1; n<fgridRegions.size(); n++)
	{
		//cout << "list region\t" << fgridRegions[n].id << "\tn=" << n << endl;
		// skip
		if (skipGridRegion(fgridRegions[n], fgridRegions[0], ntol*Mfirst))
		{
			//cout << "skip region " << fgridRegions[n].id << "\tn=" << n << endl;
			eraseCandidates.push_back(n);
		}
		// cut
		else if (fgridRegions[n].omega_r>fgridRegions[0].omega_r || fgridRegions[n].omega_l<fgridRegions[0].omega_l)
		{
			//cout << "cut region " << fgridRegions[n].id << "\tn=" << n << endl;
			if (fgridRegions[n].omega_r>fgridRegions[0].omega_r)
			{
				cutGridRegion_right(fgridRegions[0].omega_r, fgridRegions[n]);
			}
			if (fgridRegions[n].omega_l<fgridRegions[0].omega_l)
			{
				cutGridRegion_left(fgridRegions[0].omega_l, fgridRegions[n]);
			}
		}
	}
	// delete regions 
	int counter=0;
	for (int j=0; j<eraseCandidates.size(); j++)
	{
		//cout << "erase fgridregion " << eraseCandidates[j] << endl;
		fgridRegions.erase(fgridRegions.begin()+eraseCandidates[j]-counter);	
		counter++;
	}
	eraseCandidates.erase(eraseCandidates.begin(), eraseCandidates.end());
	for (int n=0; n<sgridRegions.size(); n++)
	{
		//cout << "list region\t" << sgridRegions[n].id << "\tn=" << n << endl;
		if (skipGridRegion(sgridRegions[n], fgridRegions[0], ntol*Mfirst))
		{
			//cout << "skip region " << sgridRegions[n].id << "\tn=" << n << endl;
			eraseCandidates.push_back(n);
		}
		else if (sgridRegions[n].omega_r>fgridRegions[0].omega_r || sgridRegions[n].omega_l<fgridRegions[0].omega_l)
		{
			//cout << "cut region " << sgridRegions[n].id << "\tn=" << n << endl;
			if (sgridRegions[n].omega_r>fgridRegions[0].omega_r)
			{
				cutGridRegion_right(fgridRegions[0].omega_r, sgridRegions[n]);
			}
			if (sgridRegions[n].omega_l<fgridRegions[0].omega_l)
			{
				cutGridRegion_left(fgridRegions[0].omega_l, sgridRegions[n]);
			}
		}
	}
	counter=0;
	for (int j=0; j<eraseCandidates.size(); j++)
	{
		//cout << "erase sgridregion " << eraseCandidates[j] << endl;
		sgridRegions.erase(sgridRegions.begin()+eraseCandidates[j]-counter);	
		counter++;
	}

	// clean up if there was a multigrid created before
	omega.erase(omega.begin(), omega.end());
	fsubgrids.erase(fsubgrids.begin(), fsubgrids.end());
	ssubgrids.erase(ssubgrids.begin(), ssubgrids.end());

	// *************************************
	// ***  create fundamental grid  *******
	// *************************************
	// sort grid regions
	sort(fgridRegions.begin(), fgridRegions.end(), lesserLeft);
	// check if there are grid regions with the same id
	for (int n=0; n<fgridRegions.size(); n++)
	{
		for (int m=n+1; m<fgridRegions.size(); m++)
		{
			if (fgridRegions[n].id==fgridRegions[m].id && fgridRegions[m].id!="")
			{
				cerr << "Error: multigrid: create: Two fundamental grid regions are not allowed to have the same name: " << fgridRegions[n].id << endl;
				cerr << "Break." << endl;
				throw xBadValues();
			}
		}
	}
	// create first fundamental subgrid
	if (fgridRegions[0].type=="equi")
	{
		addFSubgrid_equi(fgridRegions[0].para_i[0], fgridRegions[0].omega_l, fgridRegions[0].omega_r);
	}
	else if (fgridRegions[0].type=="tan")
	{
		addFSubgrid_tan(fgridRegions[0].para_i[0], fgridRegions[0].omega_l, fgridRegions[0].omega_r, fgridRegions[0].para_d[0], fgridRegions[0].para_d[1]);
	}
	else if (fgridRegions[0].type=="log")
	{
		addFSubgrid_log(fgridRegions[0].para_i[0], fgridRegions[0].para_i[1], fgridRegions[0].omega_l, fgridRegions[0].omega_r, fgridRegions[0].para_d[0], fgridRegions[0].para_d[1]);
	}

	// copy the the other fundamental grid regions to tfgridRegions
	vector<gridRegion> tfgridRegions;
	for (int n=1; n<fgridRegions.size(); n++)
	{
		tfgridRegions.push_back(fgridRegions[n]);
	}

	// flag decides if a grid region is skipped due to too much overlap with another grid region
	int flag;
	bool skipflag=true;
	int skipindex;
	vector<gridRegion> tfgridRegions_copy;
	while (skipflag && tfgridRegions.size()>1)
	{
		//cout << endl;
		//cout << "List of grid regions: " << endl;
		//for (int n=0; n<tfgridRegions.size(); n++)
		//{
			//cout << scientific << setprecision(5) <<  n << "\t" << tfgridRegions[n].type <<  "\tfrom\t" << tfgridRegions[n].omega_l <<"\tover\t" <<  tfgridRegions[n].omega_c <<"\tto\t" <<  tfgridRegions[n].omega_r  <<"\twith\t" <<  tfgridRegions[n].domega_min_l<<"\tand\t" <<  tfgridRegions[n].domega_min_r << endl;
		//}

		skipflag=false;
		// make a copy of tfgridRegions for the case that a region is decided to be skipped and therefore all regions
		// to the left have to be reevaluated
		tfgridRegions_copy.erase(tfgridRegions_copy.begin(), tfgridRegions_copy.end());
		for (int n=0; n<int(tfgridRegions.size()-1); n++)
		{
			tfgridRegions_copy.push_back(tfgridRegions[n]);
		}
		// get intersection points and adjust grid parameters accordingly
		//cout << endl;
		for (int n=0; n<int(tfgridRegions.size()-1); n++)
		{
			//cout << "Region " << n << " and " << n+1 << endl;
			if (tfgridRegions[n].omega_r>tfgridRegions[n+1].omega_l)
			{
				// get intersection point
				flag=getIntersection(tfgridRegions[n], tfgridRegions[n+1]);
				// if flag==0: Both regions will be realized. There may be changes to the grid parameters
				// if flag==1: the left grid region will be skipped, only the right one will be realized
				if (flag==1)
				{
					skipflag=true;
					skipindex=n;
					//cout << "skip region " << tfgridRegions[n].id << "\tn=" << n << endl;
					break;
				}
				// if flag==2: the right grid region will be skipped, only the left one will be realized
				else if (flag==2)
				{
					skipflag=true;
					skipindex=n+1;
					//cout << "skip region " << tfgridRegions[n+1].id << "\tn=" << n+1 << endl;
					break;
				}
			}
		}
		if (skipflag)
		{
			// if a region was decided to be skipped
			// overwrite all grid regions to the left of the skipped grid regions by their original implementation
			// so that they can be reevaluated in the next iteration
			tfgridRegions.erase(tfgridRegions.begin(), tfgridRegions.begin()+skipindex);
			tfgridRegions.insert(tfgridRegions.begin(), tfgridRegions_copy.begin(), tfgridRegions_copy.begin()+skipindex);

			// and erase the region to be skipped
			tfgridRegions.erase(tfgridRegions.begin()+skipindex);
		}
	}

	// create fundamental subgrids
	for (int n=0; n!=tfgridRegions.size(); n++)
	{
		if (tfgridRegions[n].type=="equi")
		{
			addFSubgrid_equi(tfgridRegions[n].para_i[0], tfgridRegions[n].omega_l, tfgridRegions[n].omega_r);
		}
		else if (tfgridRegions[n].type=="tan")
		{
			addFSubgrid_tan(tfgridRegions[n].para_i[0], tfgridRegions[n].omega_l, tfgridRegions[n].omega_r, tfgridRegions[n].para_d[0], tfgridRegions[n].para_d[1]);
		}
		else if (tfgridRegions[n].type=="log")
		{
			addFSubgrid_log(tfgridRegions[n].para_i[0], tfgridRegions[n].para_i[1], tfgridRegions[n].omega_l, tfgridRegions[n].omega_r, tfgridRegions[n].para_d[0], tfgridRegions[n].para_d[1]);
		}
	}

	// *************************************
	// *************************************
	// ***  insert special grid regions  ***
	// *************************************
	// *************************************
	// sort grid regions
	sort(sgridRegions.begin(), sgridRegions.end(), lesserCenter);
	// check if there are grid regions with the same id
	for (int n=0; n<sgridRegions.size(); n++)
	{
		for (int m=n+1; m<sgridRegions.size(); m++)
		{
			if (sgridRegions[n].id==sgridRegions[m].id && sgridRegions[m].id!="")
			{
				cerr << "Error: multigrid: create: Two special grid regions are not allowed to have the same name: " << sgridRegions[n].id << endl;
				cerr << "Break." << endl;
				throw xBadValues();
			}
		}
	}
	// copy the special grid regions to tsgridRegions
	vector <gridRegion> tsgridRegions;
	for (int n=0; n<sgridRegions.size(); n++)
	{
		tsgridRegions.push_back(sgridRegions[n]);
	}

	// flag decides if a grid region is skipped due to too much overlap with another grid region
	skipflag=true;
	skipindex;
	vector<gridRegion> tsgridRegions_copy;
	while (skipflag && tsgridRegions.size()>1)
	{
		//cout << endl;
		//cout << "List of grid regions: " << endl;
		//for (int n=0; n<tsgridRegions.size(); n++)
		//{
			//cout << scientific << setprecision(5) <<  n << "\t" << tsgridRegions[n].type <<  "\tfrom\t" << tsgridRegions[n].omega_l <<"\tover\t" <<  tsgridRegions[n].omega_c <<"\tto\t" <<  tsgridRegions[n].omega_r  <<"\twith\t" <<  tsgridRegions[n].domega_min_l<<"\tand\t" <<  tsgridRegions[n].domega_min_r << endl;
		//}

		skipflag=false;
		// make a copy of tsgridRegions for the case that a region is decided to be skipped and therefore all regions
		// to the left have to be reevaluated
		tsgridRegions_copy.erase(tsgridRegions_copy.begin(), tsgridRegions_copy.end());
		for (int n=0; n<int(tsgridRegions.size()-1); n++)
		{
			tsgridRegions_copy.push_back(tsgridRegions[n]);
		}
		// get intersection points and adjust grid parameters accordingly
		//cout << endl;
		for (int n=0; n<int(tsgridRegions.size()-1); n++)
		{
			//cout << "Region " << n << " and " << n+1 << endl;
			if (tsgridRegions[n].omega_r>tsgridRegions[n+1].omega_l)
			{
				// get intersection point
				flag=getIntersection(tsgridRegions[n], tsgridRegions[n+1]);
				// if flag==0: Both regions will be realized. There may be changes to the grid parameters
				// if flag==1: the left grid region will be skipped, only the right one will be realized
				if (flag==1)
				{
					skipflag=true;
					skipindex=n;
					//cout << "skip region " << tsgridRegions[n].id << "\tn=" << n << endl;
					break;
				}
				// if flag==2: the right grid region will be skipped, only the left one will be realized
				else if (flag==2)
				{
					skipflag=true;
					skipindex=n+1;
					//cout << "skip region " << tsgridRegions[n+1].id << "\tn=" << n+1 << endl;
					break;
				}
			}
		}
		if (skipflag)
		{
			// if a region was decided to be skipped
			// overwrite all grid regions to the left of the skipped grid regions by their original implementation
			// so that they can be reevaluated in the next iteration
			tsgridRegions.erase(tsgridRegions.begin(), tsgridRegions.begin()+skipindex);
			tsgridRegions.insert(tsgridRegions.begin(), tsgridRegions_copy.begin(), tsgridRegions_copy.begin()+skipindex);

			// and erase the region to be skipped
			tsgridRegions.erase(tsgridRegions.begin()+skipindex);
		}
	}

	// create subgrids
	for (int n=0; n!=tsgridRegions.size(); n++)
	{
		if (tsgridRegions[n].type=="equi")
		{
			addSSubgrid_equi(tsgridRegions[n].para_i[0], tsgridRegions[n].omega_l, tsgridRegions[n].omega_r);
		}
		else if (tsgridRegions[n].type=="tan")
		{
			addSSubgrid_tan(tsgridRegions[n].para_i[0], tsgridRegions[n].omega_l, tsgridRegions[n].omega_r, tsgridRegions[n].para_d[0], tsgridRegions[n].para_d[1]);
		}
		else if (tsgridRegions[n].type=="log")
		{
			addSSubgrid_log(tsgridRegions[n].para_i[0], tsgridRegions[n].para_i[1], tsgridRegions[n].omega_l, tsgridRegions[n].omega_r, tsgridRegions[n].para_d[0], tsgridRegions[n].para_d[1]);
		}
	}

	this->M=this->omega.size()-1;
	this->domega.resize(M+1);
	// calculate weights by trapez rule
	for (int j=1; j<M; j++)
	{
		domega[j]=0.5*(omega[j+1]-omega[j-1]);	
	}
	domega[0]=0.5*(omega[1]-omega[0]);	
	domega[M]=0.5*(omega[M]-omega[M-1]);	
}

// copy constructor for multigrid
multigrid::multigrid (const multigrid & mgrid)
{
	this->M=mgrid.M;
	this->omega_min=mgrid.omega_min;
	this->omega_max=mgrid.omega_max;
	for (int j=0; j!=M+1; j++)
	{
		this->omega.push_back(mgrid.omega[j]);
		this->domega.push_back(mgrid.domega[j]);
	}
	for (int n=0; n<int(mgrid.fgridRegions.size()); n++)
	{
		this->fgridRegions.push_back(mgrid.fgridRegions[n]);
	}
	for (int n=0; n<int(mgrid.sgridRegions.size()); n++)
	{
		this->sgridRegions.push_back(mgrid.sgridRegions[n]);
	}
	for (int n=0; n<int(mgrid.fsubgrids.size()); n++)
	{
		this->fsubgrids.push_back(mgrid.fsubgrids[n]);
	}
	for (int n=0; n<int(mgrid.ssubgrids.size()); n++)
	{
		this->ssubgrids.push_back(mgrid.ssubgrids[n]);
	}
}


// '=' operator for multigrid
multigrid& multigrid::operator=(const multigrid& mgrid) 
{
	this->M=mgrid.M;
	this->omega_min=mgrid.omega_min;
	this->omega_max=mgrid.omega_max;
	this->omega.erase(omega.begin(), omega.end());
	this->domega.erase(domega.begin(), domega.end());
	for (int j=0; j!=M+1; j++)
	{
		this->omega.push_back(mgrid.omega[j]);
		this->domega.push_back(mgrid.domega[j]);
	}
	this->fgridRegions.erase(fgridRegions.begin(), fgridRegions.end());
	for (int n=0; n<int(mgrid.fgridRegions.size()); n++)
	{
		this->fgridRegions.push_back(mgrid.fgridRegions[n]);
	}
	this->sgridRegions.erase(sgridRegions.begin(), sgridRegions.end());
	for (int n=0; n<int(mgrid.sgridRegions.size()); n++)
	{
		this->sgridRegions.push_back(mgrid.sgridRegions[n]);
	}
	this->fsubgrids.erase(fsubgrids.begin(), fsubgrids.end());
	for (int n=0; n<int(mgrid.fsubgrids.size()); n++)
	{
		this->fsubgrids.push_back(mgrid.fsubgrids[n]);
	}
	this->ssubgrids.erase(ssubgrids.begin(), ssubgrids.end());
	for (int n=0; n<int(mgrid.ssubgrids.size()); n++)
	{
		this->ssubgrids.push_back(mgrid.ssubgrids[n]);
	}
	return *this;
}

// copy constructor for gridRegion
gridRegion::gridRegion (const gridRegion & sgr)
{
	this->omega_l=sgr.omega_l;
	this->omega_r=sgr.omega_r;
	this->omega_c=sgr.omega_c;
	this->omega_m=sgr.omega_m;
	this->omega_p=sgr.omega_p;
	this->domega_min_l=sgr.domega_min_l;
	this->domega_min_r=sgr.domega_min_r;
	this->N_d=sgr.N_d;
	this->N_i=sgr.N_i;
	this->type=sgr.type;
	this->id=sgr.id;
	for (int n=0; n!=this->N_d; n++)
	{
		this->para_d.push_back(sgr.para_d[n]);
	}
	for (int n=0; n!=this->N_i; n++)
	{
		this->para_i.push_back(sgr.para_i[n]);
	}
}

// '=' operator for gridRegion
gridRegion& gridRegion::operator=(const gridRegion& sgr) 
{
	this->omega_l=sgr.omega_l;
	this->omega_r=sgr.omega_r;
	this->omega_c=sgr.omega_c;
	this->omega_m=sgr.omega_m;
	this->omega_p=sgr.omega_p;
	this->domega_min_l=sgr.domega_min_l;
	this->domega_min_r=sgr.domega_min_r;
	this->N_d=sgr.N_d;
	this->N_i=sgr.N_i;
	this->type=sgr.type;
	this->id=sgr.id;
	this->para_d.erase(para_d.begin(), para_d.end());
	this->para_i.erase(para_i.begin(), para_i.end());
	for (int n=0; n!=this->N_d; n++)
	{
		this->para_d.push_back(sgr.para_d[n]);
	}
	for (int n=0; n!=this->N_i; n++)
	{
		this->para_i.push_back(sgr.para_i[n]);
	}
	return *this;
}
// remove a single grid region either by id or by index
void multigrid::rem_sgr(string ID)
{
	for (int n=0; n!=this->sgridRegions.size(); n++)
	{
		if (this->sgridRegions[n].id==ID)
		{
			this->sgridRegions.erase(this->sgridRegions.begin()+n);
			return;
		}
	}
	cerr << "multigrid: rem_sgr: There is no special grid region with the id: " << ID << endl;
	cerr << "Break" << endl;
	throw 1;
}
void multigrid::rem_gr(string ID)
{
	for (int n=0; n!=this->fgridRegions.size(); n++)
	{
		if (this->fgridRegions[n].id==ID)
		{
			this->fgridRegions.erase(this->fgridRegions.begin()+n);
			return;
		}
	}
	cerr << "multigrid: rem_gr: There is no fundamental grid region with the id: " << ID << endl;
	cerr << "Break" << endl;
	throw 1;
}
void multigrid::rem_sgr(unsigned int n)
{
	if (n<=this->sgridRegions.size())
	{
		this->sgridRegions.erase(this->sgridRegions.begin()+n);
		return;
	}
	else
	{
		cerr << "multigrid: rem_sgr: There is no special grid region with the index: " << n << endl;
		cerr << "There are only " << this->sgridRegions.size() << " special grid regions" << endl;
		cerr << "Break" << endl;
		throw 1;
	}
}
void multigrid::rem_gr(unsigned int n)
{
	if (n<=this->fgridRegions.size())
	{
		this->fgridRegions.erase(this->fgridRegions.begin()+n);
		return;
	}
	else
	{
		cerr << "multigrid: rem_gr: There is no fundamental grid region with the index: " << n << endl;
		cerr << "There are only " << this->sgridRegions.size() << " special grid regions" << endl;
		cerr << "Break" << endl;
		throw 1;
	}
}
// check if grid region exists
bool multigrid::sgr_exists(string ID)
{
	for (int n=0; n!=this->sgridRegions.size(); n++)
	{
		if (this->sgridRegions[n].id==ID)
		{
			return true;
		}
	}
	return false;
}
bool multigrid::gr_exists(string ID)
{
	for (int n=0; n!=this->fgridRegions.size(); n++)
	{
		if (this->fgridRegions[n].id==ID)
		{
			return true;
		}
	}
	return false;
}
// obtain a single grid region either by id or by index
gridRegion multigrid::get_sgr(string ID)
{
	for (int n=0; n!=this->sgridRegions.size(); n++)
	{
		if (sgridRegions[n].id==ID)
		{
			return sgridRegions[n];
		}
	}
	cerr << "multigrid: get_sgr: There is no special grid region with the id: " << ID << endl;
	cerr << "Break" << endl;
	throw 1;
}
gridRegion multigrid::get_gr(string ID)
{
	for (int n=0; n!=this->fgridRegions.size(); n++)
	{
		if (fgridRegions[n].id==ID)
		{
			return fgridRegions[n];
		}
	}
	cerr << "multigrid: get_gr: There is no fundamental grid region with the id: " << ID << endl;
	cerr << "Break" << endl;
	throw 1;
}
gridRegion multigrid::get_sgr(unsigned int n)
{
	if (n<=this->sgridRegions.size())
	{
		return sgridRegions[n];
	}
	else
	{
		cerr << "multigrid: get_sgr: There is no special grid region with the index: " << n << endl;
		cerr << "There are only " << this->sgridRegions.size() << " special grid regions" << endl;
		cerr << "Break" << endl;
		throw 1;
	}
}
gridRegion multigrid::get_gr(unsigned int n)
{
	if (n<this->fgridRegions.size())
	{
		return fgridRegions[n];
	}
	else
	{
		cerr << "multigrid: get_gr: There is no fundamental grid region with the index: " << n << endl;
		cerr << "There are only " << this->fgridRegions.size() << " fundamental grid regions" << endl;
		cerr << "Break" << endl;
		throw 1;
	}
}

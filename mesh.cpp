// ****************************************************************************
// mesh 
//
// Description:	Like multigrid, but without inverse mapping
//		M:			length of the grid
//		omega_min:		minimum	
//		omega_max:		maximum	
//		omega:			container for the grid	
//		domega:			container for the weigths	
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
#include"multigrid.h"
#include"mesh.h"

using namespace std;

void mesh::testMonotony()
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
void mesh::testWeights()
{
	for (int j=0; j!=this->M+1; j++)
	{
		if (this->domega[j]<=numeric_limits<double>::epsilon())
		{
			cerr << endl;
			cerr << "ERROR: test mesh weights: Smallest grid difference is beyond machine precision ("<< numeric_limits<double>::epsilon()<<") at: j=" << j ;
			cerr << " , domega= " << scientific << setprecision(15)  << this->domega[j] << endl;
			cerr << "Break." << endl;
			throw 1;
		}
	}
}

mesh::mesh()
{
	this->lendpoint=false;
	this->rendpoint=false;
}

void mesh::add_spoint(double omega0)
{
	this->spoints.push_back(omega0);	
}
void mesh::add_lendpoint(double OMEGAL)
{
	lendpoint=true;
	this->omegal=OMEGAL;
}
void mesh::add_rendpoint(double OMEGAR)
{
	rendpoint=true;
	this->omegar=OMEGAR;
}
void mesh::create()
{
	// create multigrid
	this->mgr.create();
	// copy multigrid omega
	this->omega=vector<double> (this->mgr.M+1);
	for (int j=0; j!=this->mgr.M+1; j++)
	{
		this->omega[j]=mgr.omega[j];
	}
	this->omega_min=this->mgr.omega_min;
	this->omega_max=this->mgr.omega_max;
	double ntol=numeric_limits<double>::epsilon();
	// sort out points endpoints which are below/above the endpoints of the multigrid
	if (lendpoint && (this->omegal>=this->omega_max || this->omegal<=this->omega_min))
	{
		lendpoint=false;
	}
	if (rendpoint && (this->omegar>=this->omega_max || this->omegar<=this->omega_min))
	{
		rendpoint=false;
	}
	if (lendpoint && rendpoint && (omegar-omegal)<=ntol)
	{
		cerr << "Error: mesh: Left endpoint (" << omegal << ") must be greater (plus machineprecision) than right endpoint (" << omegar << "). Break." << endl;
		exit(1);
	}
	// sort out points which are too narrow to each other or above/below the endpoints or the endpoints of the multigrid
	vector<int> eraseCandidates;
	for (int i=0; i<this->spoints.size(); i++)
	{
		if (lendpoint && this->spoints[i]<=omegal)
		{
			//cout << "skip point " << i << " : " << this->spoints[i] << endl;
			eraseCandidates.push_back(i);
		}
		else if (rendpoint && this->spoints[i]>=omegar)
		{
			//cout << "skip point " << i << " : " << this->spoints[i] << endl;
			eraseCandidates.push_back(i);
		}
		else if (this->spoints[i]<=this->omega_min)
		{
			//cout << "skip point " << i << " : " << this->spoints[i] << endl;
			eraseCandidates.push_back(i);
		}
		else if (this->spoints[i]>=this->omega_max)
		{
			//cout << "skip point " << i << " : " << this->spoints[i] << endl;
			eraseCandidates.push_back(i);
		}
		else
		{
			for (int j=i+1; j<this->spoints.size(); j++)
			{
				if (fabs(this->spoints[i]-this->spoints[j])<=numeric_limits<double>::epsilon())
				{
					//cout << "skip point " << i << " : " << this->spoints[i] << endl;
					eraseCandidates.push_back(i);
					break;
				}
			}
		}
	}
	int counter=0;
	for (int j=0; j<eraseCandidates.size(); j++)
	{
		this->spoints.erase(this->spoints.begin()+eraseCandidates[j]-counter);	
		counter++;
	}

	// sort out points which are too narrow to the end points (if they exist) 
	eraseCandidates.clear();
	if (lendpoint || rendpoint)	
	{
		for (int j=0; j<this->spoints.size(); j++)
		{
			if (lendpoint && fabs(this->omegal-this->spoints[j])<=numeric_limits<double>::epsilon())
			{
				//cout << "skip point " << j << " : " << this->spoints[j] << endl;
				eraseCandidates.push_back(j);
			}
			else if (rendpoint && fabs(this->omegar-this->spoints[j])<=numeric_limits<double>::epsilon())
			{
				//cout << "skip point " << j << " : " << this->spoints[j] << endl;
				eraseCandidates.push_back(j);
			}
		}
		counter=0;
		for (int j=0; j<eraseCandidates.size(); j++)
		{
			this->spoints.erase(this->spoints.begin()+eraseCandidates[j]-counter);	
			counter++;
		}
	}

	// sort spoints and add endpoints to spoints container
	sort(spoints.begin(), spoints.end());
	if (lendpoint)
	{
		spoints.insert(spoints.begin(), omegal);
	}
	if (rendpoint)
	{
		spoints.push_back(omegar);
	}

	// get inverse of special points
	int I;
	vector<int> indices;
	for (int i=0; i!=this->spoints.size(); i++)
	{
		I=this->mgr.inverse(this->spoints[i]);
		if (this->spoints[i]>=this->omega[I])
		{
			//cout << "point " << spoints[i] << " at " << I << endl;
			indices.push_back(I);
		}
		else
		{
			//cout << "point " << spoints[i] << " at " << I-1 << endl;
			indices.push_back(I-1);
		}
	}

	// sort out original grid points which are too narrow to the special points
	for (int i=0; i!=this->spoints.size(); i++)
	{
		I=indices[i];
		if (fabs(this->spoints[i]-this->omega[I])<=ntol)
		{
			eraseCandidates.push_back(I);
			//cout << "skip orig. point " << I << " : " << this->omega[I] << endl;
		}
		if (fabs(this->spoints[i]-this->omega[I+1])<=ntol)
		{
			eraseCandidates.push_back(I+1);
			//cout << "skip orig. point " << I+1 << " : " << this->omega[I+1] << endl;
		}
	}

	// remove duplicates
	sort(eraseCandidates.begin(), eraseCandidates.end());
        eraseCandidates.erase(unique(eraseCandidates.begin(), eraseCandidates.end()), eraseCandidates.end());

	// calculate number of erase candidate points left of insert points
	vector<int> nleft(indices.size(), 0);
	for (int i=0; i!=indices.size(); i++)
	{
		for (int j=0; j<eraseCandidates.size(); j++)
		{
			if (eraseCandidates[j]<=indices[i])
			{
				nleft[i]++;
			}
		}		
		//cout << "nleft[" << i <<"] = " << nleft[i] << endl;
	}

	// erase eraseCandidates
	counter=0;
	for (int j=0; j<eraseCandidates.size(); j++)
	{
		this->omega.erase(this->omega.begin()+eraseCandidates[j]-counter);
		counter++;
	}

	// insert special points
	counter=0;
	int lindex, rindex;
	for (int i=0; i!=this->spoints.size(); i++)
	{
		// save indices of left and right endpoints
		if (lendpoint && i==0)
		{
			lindex=indices[i]+1-nleft[i]+counter;
		}
		if (rendpoint && i==this->spoints.size()-1)
		{
			rindex=indices[i]+1-nleft[i]+counter;
		}
		//cout << "insert grid point " << indices[i]+1-nleft[i]+counter << " : " << this->spoints[i] << endl;
		this->omega.insert(this->omega.begin()+indices[i]+1-nleft[i]+counter, this->spoints[i]);
		counter++;	
	}

	// cut all original grid points left from omegal and right from omegar
	if (rendpoint)
	{
		this->omega.erase(this->omega.begin()+rindex+1,         this->omega.end()          );
	}
	if (lendpoint)
	{
		this->omega.erase(this->omega.begin()         ,         this->omega.begin()+lindex);
	}

	// trapez rule
	this->M=this->omega.size()-1;
	this->domega.resize(this->omega.size());
	for (int j=1; j<this->M; j++)
	{
		this->domega[j]=0.5*(this->omega[j+1]-this->omega[j-1]);	
	}
	this->domega[      0]=0.5*(this->omega[      1]-this->omega[        0]);	
	this->domega[this->M]=0.5*(this->omega[this->M]-this->omega[this->M-1]);	
}


// wrapper for multigrid member functions

void mesh::add_gr_equi(int N, double omega_min, double omega_max, double omega_c)
{
	this->mgr.add_gr_equi(N, omega_min, omega_max, omega_c);
}
void mesh::add_gr_equi(int N, double omega_min, double omega_max, double omega_c, string id)
{
	this->mgr.add_gr_equi(N, omega_min, omega_max, omega_c, id);
}
void mesh::add_gr_tan(int N, double omega_min, double omega_max, double omega_c, double c)
{
	this->mgr.add_gr_tan(N, omega_min, omega_max, omega_c, c);
}
void mesh::add_gr_tan(int N, double omega_min, double omega_max, double omega_c, double c, string id)
{
	this->mgr.add_gr_tan(N, omega_min, omega_max, omega_c, c, id);
}
void mesh::add_gr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0)
{
	this->mgr.add_gr_log(N, M, omega_min, omega_max, omegak, omegak_0);
}
void mesh::add_gr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0, string id)
{
	this->mgr.add_gr_log(N, M, omega_min, omega_max, omegak, omegak_0, id);
}

void mesh::add_sgr_equi(int N, double omega_min, double omega_max, double omega_c)
{
	this->mgr.add_sgr_equi(N, omega_min, omega_max, omega_c);
}
void mesh::add_sgr_equi(int N, double omega_min, double omega_max, double omega_c, string id)
{
	this->mgr.add_sgr_equi(N, omega_min, omega_max, omega_c, id);
}
void mesh::add_sgr_tan(int N, double omega_min, double omega_max, double omega_c, double c)
{
	this->mgr.add_sgr_tan(N, omega_min, omega_max, omega_c, c);
}
void mesh::add_sgr_tan(int N, double omega_min, double omega_max, double omega_c, double c, string id)
{
	this->mgr.add_sgr_tan(N, omega_min, omega_max, omega_c, c, id);
}
void mesh::add_sgr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0)
{
	this->mgr.add_sgr_log(N, M, omega_min, omega_max, omegak, omegak_0);
}
void mesh::add_sgr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0, string id)
{
	this->mgr.add_sgr_log(N, M, omega_min, omega_max, omegak, omegak_0, id);
}

void mesh::replace_sgr_equi(int N, double omega_min, double omega_max, double omega_c, string id)
{
	this->mgr.replace_sgr_equi(N, omega_min, omega_max, omega_c, id);
}
void mesh::replace_sgr_tan(int N, double omega_min, double omega_max, double omega_c, double c, string id)
{
	this->mgr.replace_sgr_tan(N, omega_min, omega_max, omega_c, c, id);
}
void mesh::replace_sgr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0, string id)
{
	this->mgr.replace_sgr_log(N, M, omega_min, omega_max, omegak, omegak_0, id);
}
void mesh::replace_gr_equi(int N, double omega_min, double omega_max, double omega_c, string id)
{
	this->mgr.replace_gr_equi(N, omega_min, omega_max, omega_c, id);
}
void mesh::replace_gr_tan(int N, double omega_min, double omega_max, double omega_c, double c, string id)
{
	this->mgr.replace_gr_tan(N, omega_min, omega_max, omega_c, c, id);
}
void mesh::replace_gr_log(int N, int M, double omega_min, double omega_max, double omegak, double omegak_0, string id)
{
	this->mgr.replace_gr_log(N, M, omega_min, omega_max, omegak, omegak_0, id);
}

// wrapper for the equi-gridregion functions
void mesh::add_gr_equi(double omega_l, double omega_r, double domegap)
{
	this->mgr.add_gr_equi(omega_l, omega_r, domegap);
}
void mesh::add_gr_equi(double omega_l, double omega_r, double domegap, string ID)
{
	this->mgr.add_gr_equi(omega_l, omega_r, domegap, ID);
}
void mesh::add_sgr_equi(double omega_l, double omega_r, double domegap)
{
	this->mgr.add_sgr_equi(omega_l, omega_r, domegap);
}
void mesh::add_sgr_equi(double omega_l, double omega_r, double domegap, string ID)
{
	this->mgr.add_sgr_equi(omega_l, omega_r, domegap, ID);
}
void mesh::replace_gr_equi(double omega_l, double omega_r, double domegap, string ID)
{
	this->mgr.replace_gr_equi(omega_l, omega_r, domegap, ID);
}
void mesh::replace_sgr_equi(double omega_l, double omega_r, double domegap, string ID)
{
	this->mgr.replace_sgr_equi(omega_l, omega_r, domegap, ID);
}


// wrapper for the log-Gridregion functions
void mesh::add_gr_log(double omega_c, double omega1, double domega_min, double domega_max)
{
	this->mgr.add_gr_log(omega_c, omega1, domega_min, domega_max);
}
void mesh::add_gr_log(double omega_c, double omega1, double domega_min, double domega_max, string id)
{
	this->mgr.add_gr_log(omega_c, omega1, domega_min, domega_max, id);
}
void mesh::add_sgr_log(double omega_c, double omega1, double domega_min, double domega_max)
{
	this->mgr.add_sgr_log(omega_c, omega1, domega_min, domega_max);
}
void mesh::add_sgr_log(double omega_c, double omega1, double domega_min, double domega_max, string id)
{
	this->mgr.add_sgr_log(omega_c, omega1, domega_min, domega_max, id);
}
void mesh::replace_gr_log(double omega_c, double omega1, double domega_min, double domega_max, string id)
{
	this->mgr.replace_gr_log(omega_c, omega1, domega_min, domega_max, id);
}
void mesh::replace_sgr_log(double omega_c, double omega1, double domega_min, double domega_max, string id)
{
	this->mgr.replace_sgr_log(omega_c, omega1, domega_min, domega_max, id);
}

void mesh::add_sgr(string id, gridRegion & gr)
{
	this->mgr.add_sgr(id, gr);
}
void mesh::add_gr(string id, gridRegion & gr)
{
	this->mgr.add_gr(id, gr);
}
void mesh::replace_sgr(string id, gridRegion & gr)
{
	this->mgr.replace_sgr(id, gr);
}
void mesh::replace_gr(string id, gridRegion & gr)
{
	this->mgr.replace_gr(id, gr);
}

void mesh::rem_sgr(string id)
{
	this->mgr.rem_sgr(id);
}
void mesh::rem_sgr(unsigned int n)
{
	this->mgr.rem_sgr(n);
}
void mesh::rem_gr(string id)
{
	this->mgr.rem_gr(id);
}
void mesh::rem_gr(unsigned int n)
{
	this->mgr.rem_gr(n);
}
bool mesh::sgr_exists(string id)
{
	return this->mgr.sgr_exists(id);
}
bool mesh::gr_exists(string id)
{
	return this->mgr.gr_exists(id);
}

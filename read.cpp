#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "read.hpp"

using namespace std;

ReadC_O::ReadC_O(){
    ifstream f("oxygen_concentration.txt");
    string line;
    vector<double> T;
    vector<double> oversaturation;
    int N = 0;

    while (getline(f, line)) {
        istringstream ss(line);
        string col1, col2;
        ss >> col1 >> col2;
        T.push_back(double(atof(col1.c_str())));
        oversaturation.push_back(double(atof(col2.c_str())));
        N++;
    }

    double *T_for_interp;
    double *oversaturation_for_interp;
    T_for_interp = new double [N];
    oversaturation_for_interp = new double [N];
    for (int i=0; i<N; i++){
        T_for_interp[i] = T[N-1-i];
        oversaturation_for_interp[i] = oversaturation[N-1-i];
    }

    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_cspline, N);

    gsl_spline_init (spline, T_for_interp, oversaturation_for_interp, N);
    
    delete [] T_for_interp;
    delete [] oversaturation_for_interp;
}



ReadC_O::~ReadC_O(){
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
}


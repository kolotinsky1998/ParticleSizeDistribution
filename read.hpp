#pragma once

class ReadC_O{
    public:
        ReadC_O();
	~ReadC_O();
	gsl_spline *spline;
	gsl_interp_accel *acc;
};

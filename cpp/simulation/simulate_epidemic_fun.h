#ifndef SIMULATE_EPIDEMIC_FUN_PERF_PKG__H
#define SIMULATE_EPIDEMIC_FUN_PERF_PKG__H

#include <Rcpp.h>
using namespace Rcpp;

extern double mIncub;
extern double sdIncub;
extern double maxPCRDetectability;
extern double m_incub_g;
extern double sd_incub_g;
extern double shape_incub_g;
extern double scale_incub_g;

extern double mGamma;
extern double vGamma;
extern double shift;
extern double shape_infectivity;
extern double scale_infectivity;

extern double shape_gt;
extern double scale_gt;

extern double mu_f;
extern double sigma_f; 
extern double alpha_f;
extern double tau;

extern double cf;


RcppExport double dtost(double tost, double inc);
RcppExport double ptost(double tost, double inc);
RcppExport void test();
RcppExport double incubPeriod_g();
RcppExport double incubPeriod();

RcppExport NumericVector foi(
	double t, 
	double dt, 
	double lastDate,  
	int ageCurr,
	NumericVector incub_periods,
	NumericVector symptomOnset, 
	NumericVector infectionDate,
	IntegerVector age, 
	IntegerVector infStatus, 
	double alpha,
	double beta,
	double rInfC,
	double hhsize,
  double mainHHSize,
	std::string contact_pattern = "homogeneous",
	int display=0
  );

#endif

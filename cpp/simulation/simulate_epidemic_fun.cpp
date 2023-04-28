#include <Rcpp.h>
#include "simulate_epidemic_fun.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

double mIncub = 1.63;
double sdIncub = 0.5;
double maxPCRDetectability = 10.0;

double m_incub_g = 4.07;
double sd_incub_g = 2.12;
double shape_incub_g = pow(m_incub_g,2) / pow(sd_incub_g, 2);
double scale_incub_g = pow(sd_incub_g,2) / m_incub_g;

double mGamma = 26.1;
double vGamma = 7.0;
double shift = 25.6;
double shape_infectivity = pow(mGamma,2) / vGamma;
double scale_infectivity = vGamma / mGamma;

double shape_gt = 4.0;
double scale_gt = 5.0/4.0;

double mu_f = -4.0;
double sigma_f = 1.85; 
double alpha_f = 5.85;
double tau = exp(mIncub+0.5*pow(0.5, 2)); 

double cf = alpha_f / (sigma_f*(1- pow(1+exp((mu_f+tau)/sigma_f), -alpha_f)));
  

//////////////////////////////////////////////
// [[Rcpp::export]]
double dtost(double tost, double inc) {
	double out(0.0);

	if (tost<0) {
		out = exp(-(tost*tau/inc-mu_f)/sigma_f) / pow(1+exp(-(tost*tau/inc-mu_f)/sigma_f), alpha_f+1);
	} else {
		out = exp(-(tost-mu_f)/sigma_f) / pow(1+exp(-(tost-mu_f)/sigma_f),alpha_f+1);
	}

	out *= cf;
	return out;
}

//////////////////////////////////////////////
// [[Rcpp::export]]
double ptost(double tost, double inc) {
	double out(0.0);

	if (tost<0) {
		out = inc / tau * ( pow(1+exp(-(tost*tau/inc-mu_f)/sigma_f), -alpha_f) - pow(1+exp((mu_f+tau)/sigma_f), -alpha_f) );
	
	} else {
		// Integrate from -inc to 0
		out = inc / tau * ( pow(1+exp(mu_f/sigma_f), -alpha_f) - pow(1+exp((mu_f+tau)/sigma_f), -alpha_f) );

		// Integrate from 0 to tost
		out +=  pow(1+exp(-(tost-mu_f)/sigma_f), -alpha_f) - pow(1+exp(mu_f/sigma_f), -alpha_f) ;
	}

	out *= cf * sigma_f / alpha_f;
	return out;
}

//////////////////////////////////////////////
// [[Rcpp::export]]
void test() {
	Rcout << "Mean incubation period: " << tau << endl;
	Rcout << "Cumulative tost: " << cf * sigma_f / alpha_f * pow(1+exp(-(40-mu_f)/sigma_f), -alpha_f) << endl;
}


//////////////////////////////////////////////
// [[Rcpp::export]]
double incubPeriod_g() {
	double d = rgamma(1, shape_incub_g, scale_incub_g)[0];
	//while(d<1) d = rgamma(1, shape_incub_g, scale_incub_g)[0];
	return d;
}


//////////////////////////////////////////////
// [[Rcpp::export]]
double incubPeriod() 
  {
  return rlnorm(1, mIncub, sdIncub)[0];
}

//////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector foi(
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
	std::string contact_pattern,
	int display
  ) {
	// Force of infection from 
	// 	0: community
	// 	1 to #infected: infected household contacts
	NumericVector fois(symptomOnset.size() + 1);

	//Infection by the community
	fois[0] = alpha*dt;

  	// Infection by infected individuals within the same household
	for (int index = 0; index < symptomOnset.size(); ++index) {

		double k = 0.0;

		if ( t >= infectionDate[index] && infectionDate[index] < lastDate ) {
			k+= ptost(t+dt-symptomOnset[index], incub_periods[index]) - ptost(t-symptomOnset[index], incub_periods[index]);
		} 
		
		// Foi depends on household size - non parametric formulation (reference = child-child in households of size 4)
		fois[index + 1] = beta * k / (hhsize / mainHHSize);
		
		// Relative infectivity of adults
		if ( age[index] < 2 ) fois[index + 1] /= rInfC;

		// Contact rate between infector - infectee pairs
		if (contact_pattern == "heterogeneous") {
			// Child-Mother
			if ( (age[index] == 0 && ageCurr == 2) || (age[index] == 2 && ageCurr == 0) ) fois[index+1] *= 1.17;

			// Child-Father
			if ( (age[index] == 2 && ageCurr == 1) || (age[index] == 1 && ageCurr == 2) ) fois[index+1] *= 0.55;
            
      // Mother-Father
			if ( (age[index] == 0 && ageCurr == 1) || (age[index] == 1 && ageCurr == 0) ) fois[index+1] *= 1.31;
			
		}

	} 
  
  return fois;
}

#include <cmath>
#include <random>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/gamma.hpp>
#include "utils.h"

using namespace std;


double mIncub = 1.63;
double sdIncub = 0.5;

double mu_f = -4.0;
double sigma_f = 1.85; 
double alpha_f = 5.85;
double tau = exp(mIncub+0.5*pow(0.5, 2));

double cf = alpha_f / (sigma_f*(1- pow(1+exp((mu_f+tau)/sigma_f), -alpha_f)));


//----------------------------------------------------------------------
// Cumulative probability density functions
//----------------------------------------------------------------------
double plnorm(double x, double mean, double sd) {
   boost::math::lognormal_distribution<> d(mean, sd);
   return cdf(d, x);
}

double pgamma(double x, double shape, double scale) {
   boost::math::gamma_distribution<> d(shape, scale);
   return cdf(d, x);
}

double ptost(double tost, double inc) {
    double out(0.0);

    if (tost<0) {
        out = inc / tau * ( pow(1+exp(-(tost*tau/inc-mu_f)/sigma_f), -alpha_f) - pow(1+exp((mu_f+tau)/sigma_f), -alpha_f) );
    
    } else {
        // Integrate from -inc to 0
        out = inc / tau * ( pow(1+exp(mu_f/sigma_f), -alpha_f) - pow(1+exp((mu_f+tau)/sigma_f), -alpha_f) );

        // Integrate from 0 to tost
        out += pow(1+exp(-(tost-mu_f)/sigma_f), -alpha_f) - pow(1+exp(mu_f/sigma_f), -alpha_f) ;
    }

    out *= cf * sigma_f / alpha_f;
    return out;
}

//----------------------------------------------------------------------
// Probability density functions
//----------------------------------------------------------------------
double dlnorm(double x, double mean, double sd) {
    return exp( - pow(log(x) - mean, 2) / ( 2 * pow(sd,2) ) ) / (x * sd * 2.50662827);
}

double dgamma(double x, double shape, double scale) {
   boost::math::gamma_distribution<> d(shape, scale);
   return pdf(d, x);
}

// Log probability density distributions
double logdexp(double x, double rate) {
    return log(rate) - rate * x;
}

double logdlnorm(double x, double mean, double sd) {
    return -log(x) - log(sd) - 0.918938533 - pow(log(x)-mean, 2) / (2.0 * pow(sd, 2));
}

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

//----------------------------------------------------------------------
// Random generator
//----------------------------------------------------------------------
double rnorm(std::mt19937_64& gen, double mean, double sd) {
    std::normal_distribution<double> dist(mean, sd);
    return dist(gen);
}

double runif(std::mt19937_64& gen, double lower_value, double upper_value) {
    std::uniform_real_distribution<double> dist(lower_value, upper_value);
    return dist(gen);
}

double rlnorm(std::mt19937_64& gen, double mean, double sd) {
    std::lognormal_distribution<double> dist(mean, sd);
    return dist(gen);
}

double rgamma(std::mt19937_64& gen, double shape, double scale) {
    std::gamma_distribution<double> dist(shape, scale);
    return dist(gen);
}

double rexp(std::mt19937_64& gen, double rate) {
    std::exponential_distribution<double> dist(rate);
    return dist(gen);
}


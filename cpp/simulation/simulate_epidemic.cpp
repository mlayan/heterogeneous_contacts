#include <Rcpp.h>
#include "simulate_epidemic_fun.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// Mother: adult = 0
// Father: adult = 1
// Child: adult = 2

////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
DataFrame hhEpidemic(
	DataFrame H,
  double alpha,
 	double beta,
	double rSC,
	double rInfC,
  int followUp,
  double mainHHSize,
  std::string contact_pattern,
	double dt
  )
  {

  // Simulate epidemic within a household

  //Rcout << "Alpha: " << alpha << endl; 
  //Rcout << "Beta 4: " << beta << endl; 
  //Rcout << "rSC: " << rSC << endl;
  //Rcout << "rInfSC: " << rInfC << endl;  
  //Rcout << "Main hh size: " << mainHHSize << endl;
  //Rcout << "Contact pattern: " << contact_pattern << endl;  

  // Data info
  int hhsize = H.nrows();
  IntegerVector adult = H["adult"];

  // Initialize transmission vectors
  IntegerVector infectionStatus = H["infectionStatus"];
  CharacterVector origins(hhsize);
  NumericVector dds_temp = H["dds"];
  NumericVector dds = clone(dds_temp);
  NumericVector ddi = rep(1000.0, hhsize);
  NumericVector incub_period = rep(1000.0, hhsize);

  // Susceptible and index cases at the start of the epidemic
  IntegerVector infectors(0);
  IntegerVector sus(0);
  double lastDate = 0;
  double firstInfection = 1000.0;
  double startFollowUp = 1000.0;

  for (int ind=0; ind<hhsize; ind++) {

  	if (dds[ind] != 1000.0) { // Index cases

  	infectors.push_back(ind);

    	if ( infectionStatus[ind] == 1 ) { // Symptomatic cases 
    	//ddi[ind] = rIncub(dds[ind]);
        //ddi[ind] = dds[ind] - incubPeriod();
        incub_period[ind] = incubPeriod();

    	} else { // Asymptomatic index cases are detected 5 days after infection
    	//ddi[ind] = dds[ind] - detectionPeriod2(0.0, 1.0);
        //ddi[ind] = dds[ind] - runif(1, 2.0, 7.0)[0];
        incub_period[ind] = runif(1, 2.0, 7.0)[0];        

    	}
      
      ddi[ind] = dds[ind] - incub_period[ind];
      startFollowUp = min(startFollowUp, dds[ind]);
    	firstInfection = min(firstInfection, ddi[ind]);

    } else { // Susceptible household members
    	sus.push_back(ind);
    }
  }

  lastDate = startFollowUp + followUp;


  // Simulate Epidemic within household
  for (int n = 0; n <= (lastDate-firstInfection)/dt; ++n) { 

   // Current time
    double curr_time = firstInfection + dt * n;  
    
    // Draw new infections
    IntegerVector newInfected;
    NumericVector infected = runif(sus.size());
    
    // Loop on susceptibles 
    for (int s = 0; s < sus.size(); ++s) {
      // Individual
      int ind = sus[s];

      // Probability of getting infected
      int display=0;
      //if (n==0 && s==0) display = 1;

      NumericVector FOIS = foi(
        curr_time, 
        dt, 
        lastDate,
        adult[ind],
        incub_period[infectors],
        dds[infectors], 
        ddi[infectors], 
        adult[infectors], 
        infectionStatus[infectors], 
        alpha, 
        beta, 
        rInfC,
        hhsize,
        mainHHSize,
        contact_pattern, 
        display
        );
      
      // Update foi according to age of contact 
      NumericVector fois = clone(FOIS);

      if (fois.size() > 1) { // fois[0] corresponds to the infection within the community
        if ( adult[ind] < 2 ) { // Relative susceptibility of adults
          for (int i = 1; i<fois.size(); ++i) fois[i] /= rSC;

        }
      }
      
      // Is ind infected at curr_time?
      double pInfection = sum(fois); 

      if ( infected[s] < pInfection) {
        // Update ddi and dds upon infection
        ddi[ind] = curr_time + dt*runif(1)[0];
        infectionStatus[ind] = 1;
        incub_period[ind] = incubPeriod();
        dds[ind] = ddi[ind] + incub_period[ind];
             
        // Update new infected individuals
        newInfected.push_back(ind);
      }
    }
    
    // Update sus and infectors
    for (int index = 0; index < newInfected.size(); ++index) {
      
      infectors.push_back(newInfected[index]);

      for (int i = 0; i < sus.size(); ++i) {
        if (newInfected[index] == sus[i]) {
          sus.erase(i);
        }
      }
    }
  }
  
  // Output
  return DataFrame::create(
  	_("indid") = H["indid"],
  	_("hhid") = H["hhid"],
  	_("hhsize") = H["hhsize"],
  	_("ddi") = ddi,
  	_("dds") = dds,
  	_("infectionStatus") = infectionStatus,
    _("startFollowUp") = rep(startFollowUp, hhsize),
  	_("studyPeriod") = rep(lastDate, hhsize),
  	_("adult") = adult,
  	_("index") = H["index"],
    _("origins") = origins    
  	);
}





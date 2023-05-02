#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include "Household.h"
#include "utils.h"

using namespace std;

// Constant parameters
double mGamma = 26.1;
double vGamma = 7.0;
double shift = 25.6;
double shapeInf = pow(mGamma, 2) / vGamma;
double scaleInf = vGamma / mGamma;

double mIncub_g=4.07;
double sdIncub_g=2.12;
double shape_incub=pow(mIncub_g,2)/pow(sdIncub_g,2);
double scale_incub=pow(sdIncub_g,2)/mIncub_g;

double shape_gt = 4.0;
double scale_gt = 5.0/4.0;


//----------------------------------------------------------------------
// Infectivity profile
//----------------------------------------------------------------------
double infectivityProfile(
                          double origin,
                          double infectionDate,
                          int infectionStatus,
                          double tinf,
                          int studyPeriod,
                          int display
                          ) {

  // Density
  double out = 0.0;

  if (infectionStatus == 1) { // Symptomatic infector
    if (infectionDate < studyPeriod && infectionDate < tinf) {
      out = dtost(tinf-origin, origin-infectionDate); 
    }
  } 

  return out;
}


double cumulativeInfectivity(
                             double origin,
                             double infectionDate,
                             int infectionStatus,
                             double tinf,
                             int studyPeriod,
                             int display
                             )
{

  // Cumulative infectivity profile
  double out = 0.0;

  // Cumulative density
  if (infectionStatus == 1) { // Symptomatic infector

    if (infectionDate < studyPeriod && infectionDate < tinf) {
      out = ptost(tinf-origin, origin-infectionDate); 
    }
  }

  return out;
}

//----------------------------------------------------------------------
// Household class - Methods
//----------------------------------------------------------------------
Household::Household() : m_size(0), m_notInfected(0), m_startFollowUp(-1) {

  m_indid.resize(0);
  m_onsetTime.resize(0);
  m_augmentedOnsetTime.resize(0);
  m_infected.resize(0);
  m_confCase.resize(0);
  m_studyPeriod.resize(0);
  m_age.resize(0);
  m_infTime.resize(0);
  m_cumLambda.resize(0);
  m_instLambda.resize(0); 

  m_contactPattern= "homogeneous"; 
}



Household::Household(std::string contact_pattern) {

  m_size=0;
  m_notInfected=0;
  m_startFollowUp=-1;

  m_indid.resize(0);
  m_onsetTime.resize(0);
  m_augmentedOnsetTime.resize(0);
  m_infected.resize(0);
  m_confCase.resize(0);
  m_studyPeriod.resize(0);
  m_age.resize(0);
  m_infTime.resize(0);
  m_cumLambda.resize(0);
  m_instLambda.resize(0); 

  m_contactPattern=contact_pattern; 
}




std::vector<int> Household::getSpInfected(std::vector<int> index) const{
    unsigned i, k;
    std::vector<int> output;
    output.resize(0);

    for(i=0; i<index.size();i++) {
        k = index[i];
        output.push_back(m_infected[k]);
    }

  return output;
}

size_t Household::nSymptomatic() const{
  size_t nSymptomatic = 0;
  for (size_t i=0; i<m_size; i++) {
    if (m_onsetTime[i] < 1000.0) nSymptomatic += 1;
  }

  return nSymptomatic;
}


std::string Household::getindid() const {
  std::string indid;
  for (size_t i=0; i<m_size;i++) {
    indid += std::to_string(m_indid[i]) + " "; 
  }

  return indid;
}


// Add individual to household
void Household::addIndividual(
  int indid, 
  double onsetTime, 
  int isCase, 
  int age, 
  int startFollowUp,
  int studyPeriod, 
  double ddi
  ) {

  m_size += 1; // Update size

  // Update vector attributes
  m_indid.push_back(indid);
  m_onsetTime.push_back(onsetTime);
  m_augmentedOnsetTime.push_back(onsetTime);
  m_infected.push_back(isCase);
  m_studyPeriod.push_back(studyPeriod);
  m_age.push_back(age);
  
  m_startFollowUp = startFollowUp;

  // Add confirmed cases
    if (isCase > 0) { // 1: symptomatic cases, 2: asymptomatic cases, 3: symptomatic cases with unknown symptom onset
        m_confCase.push_back(m_size-1);
    } else {
        m_notInfected += 1;
    }

  // Initialize infection time
    m_infTime.push_back(ddi);
}


// Change infection time or symptom onset
void Household::setOnsetTime(int index, double onsetTime) {
  m_augmentedOnsetTime[index] = onsetTime;
}

void Household::setInfTime(int index, double infTime) {
  m_infTime[index] = infTime;
}


// Reinitialize object
void Household::newHousehold() {
  m_size = 0;
  m_notInfected = 0;  
  m_startFollowUp = -1;
  m_indid.clear();
  m_onsetTime.clear();
  m_augmentedOnsetTime.clear();
  m_infected.clear();
  m_confCase.clear();
  m_studyPeriod.clear();
  m_age.clear();
  m_infTime.clear();
}


// Initialize augmented data
void Household::initialOnsetTime(int index, std::mt19937_64& gen, int display) {

    if ( m_infected[index] > 0 ) {   // Infected household members 
      m_augmentedOnsetTime[index] = m_onsetTime[index] + runif(gen, 0.0, 1.0);
    } else {
      m_augmentedOnsetTime[index] = 1000.0;
    }

    if (display) cout << "Observed onset time: " << m_onsetTime[index] << " - New onset time" << m_augmentedOnsetTime[index] << endl;
}


void Household::initialInfTime(int index, double maxPCRDetectability, std::mt19937_64& gen, int display) {

    double incubPeriod(0.0);

    if ( m_infected[index] == 1 ) {   // Symptomatic case 
      incubPeriod = rlnorm(gen, mIncub, sdIncub);
      m_infTime[index] = m_augmentedOnsetTime[index] - incubPeriod;

    } else {    // Covid-free household members or household members with unknown final outcome
      m_infTime[index] = 1000.0; 
    }

    if (display) cout << "Inf time " << m_indid[index] << ": " << m_infTime[index] << endl;
}




// Compute the risk of infection
void Household::compute_infectivity_profiles(int display) {
  
  size_t infector, infectee;
  double t1;

  // Set m_cumLambda and m_instLambda size
  m_cumLambda.resize(m_size);
  m_instLambda.resize(m_size);
  for (infector=0; infector<m_size; infector++) {
    m_cumLambda[infector].resize(m_size);
    m_instLambda[infector].resize(m_size);
  }

  // Compute cumulative and instantaneous person-to-person transmission rates
  for (infector=0; infector<m_size; infector++) {
    for (infectee=0; infectee<m_size; infectee++) {

      if ( m_infected[infectee] == 0 ) {
        t1 = m_studyPeriod[infectee];
      }else{
        t1 = m_infTime[infectee];
      }

      if (display) cout << "Infector: " << infector << " " << m_infTime[infector]  << " " << m_augmentedOnsetTime[infector] << " Infectee: " << infectee << " " << t1 << " " << m_infected[infectee] << endl;

      // Cumulative transmission rate for household contacts
      if (m_infected[infector] > 0 && m_infTime[infector] < t1 && m_infected[infectee] >=0) {

        m_cumLambda[infector][infectee] += cumulativeInfectivity(
          m_augmentedOnsetTime[infector], 
          m_infTime[infector], 
          m_infected[infector], 
          t1, 
          m_studyPeriod[infectee],
          display
        );

      }

      // Instantaneous transmission rate for secondary cases
      if (m_infected[infector] > 0 && m_infTime[infector] < m_infTime[infectee] && m_infected[infectee] > 0) {

        m_instLambda[infector][infectee] += infectivityProfile(
          m_augmentedOnsetTime[infector], 
          m_infTime[infector], 
          m_infected[infector],
          m_infTime[infectee], 
          m_studyPeriod[infectee],
          display
        );
      }

    }
  }

   if (display) {
    // Cumulative transmission rate
    for (size_t i=0; i<m_size;i++) {
      for (size_t j=0; j<m_size; j++) {
        cout << m_cumLambda[i][j] << " "; 
      }
      cout << endl;
    }
    cout << endl;

    // Instantaneous transmission rate
    for (size_t i=0; i<m_size;i++) {
      for (size_t j=0; j<m_size; j++) {
        cout << m_instLambda[i][j] << " "; 
      }
      cout << endl;
    }
    cout << endl;

    // Infection time 
    for (size_t i=0; i<m_size; i++) {
      cout << m_infTime[i] << " ";
    }
    cout << endl;
  }

}

// Update m_cumLambda and m_instLambda matrices during data augmentation
void Household::update_infectivity_profiles(int ind, int display) {
  
  size_t infector, infectee;
  double t1;

  // Ind is the infector
  for (infectee=0; infectee<m_size; infectee++) {
    // Reinitialize
    m_cumLambda[ind][infectee] = 0;
    m_instLambda[ind][infectee] = 0;

    // Cumulative transmission rate for all household contacts
    if ( m_infected[infectee] == 0 ) {
      t1 = m_studyPeriod[infectee];
    }else{
      t1 = m_infTime[infectee];
    }

    if (m_infected[ind] > 0 && m_infTime[ind] < t1 && m_infected[infectee] >=0) {

      m_cumLambda[ind][infectee] = cumulativeInfectivity(
      m_augmentedOnsetTime[ind], 
      m_infTime[ind], 
      m_infected[ind],
      t1, 
      m_studyPeriod[infectee]
      );
    }

    // Instantaneous transmission rate for secondary cases
    if (m_infected[ind] > 0 && m_infTime[ind] < m_infTime[infectee] && m_infected[infectee] > 0) {
      m_instLambda[ind][infectee] = infectivityProfile(
      m_augmentedOnsetTime[ind], 
      m_infTime[ind], 
      m_infected[ind], 
      m_infTime[infectee], 
      m_studyPeriod[infectee]
      );
    }
  }


  // Ind is an infectee
  if (m_infected[ind] == 0) {
    t1 = m_studyPeriod[ind];
  } else {
    t1 = m_infTime[ind];    
  }

  for (infector=0; infector<m_size;infector++) {    
    // Reinitialize
    m_cumLambda[infector][ind] = 0;
    m_instLambda[infector][ind] = 0;

    // Cumulative transmission rate for all household contacts
    if (m_infected[infector] > 0 && m_infTime[infector] < t1 && m_infected[ind] >=0) {

      m_cumLambda[infector][ind] = cumulativeInfectivity(
      m_augmentedOnsetTime[infector], 
      m_infTime[infector], 
      m_infected[infector], 
      t1, 
      m_studyPeriod[ind]
      );
    }

    // Instantaneous transmission rate for secondary cases
    if (m_infected[infector] > 0 && m_infTime[infector] < m_infTime[ind] && m_infected[ind] > 0) {
      m_instLambda[infector][ind] = infectivityProfile(
      m_augmentedOnsetTime[infector], 
      m_infTime[infector], 
      m_infected[infector], 
      m_infTime[ind], 
      m_studyPeriod[ind]
      );
    }
  }

  if (display) {
    for (size_t i=0; i<m_size;i++) {
      for (size_t j=0; j<m_size; j++) {
        cout << m_cumLambda[i][j] << " "; 
      }
      cout << endl;
    }
    
    for (size_t i=0; i<m_size; i++) {
      cout << m_infTime[i] << " ";
    }
    cout << endl;
  }

}



// Propose new value for infection time and symptom onset time
double Household::newOnsetTime(int index, std::mt19937_64& gen) {

    double onsetTime(0.0);

    if ( m_infected[index] > 0 ) {           // Infected household member 
      onsetTime = m_onsetTime[index] + runif(gen, 0.0, 1.0);

    } else {    // Covid-free household members or household members with unknown final outcome
      onsetTime = 1000.0;
    }

    return onsetTime;
}



double Household::newInfTime(int index, double maxPCRDetectability, std::mt19937_64& gen) {

    double infDate(0.0), incubPeriod(0.0);

    if ( m_infected[index] == 1 ) {           // Symptomatic case 

      incubPeriod = rlnorm(gen, mIncub, sdIncub);
      infDate += m_augmentedOnsetTime[index] - incubPeriod;

    } else {    // Covid-free household members or household members with unknown final outcome
      infDate = 1000.0;
    }

    return infDate;
}





double Household::log_dIncub(int index, double maxPCRDetectability, int display) {

    double out = 0.0;

    if ( m_infected[index] == 1 ) {       // Symptomatic case with known symptom onset 
      double incubPeriod = m_augmentedOnsetTime[index] - m_infTime[index];
      out += logdlnorm(incubPeriod, mIncub, sdIncub);
    }

    if (display) cout << "log_dIncub: " << out << endl;

    return out;
}





double Household::log_dIncub_g(int index, double maxPCRDetectability, int display) {

    double out = 0.0;

    if ( m_infected[index] == 1 ) {       // Symptomatic case with known symptom onset 
      double incubPeriod = m_augmentedOnsetTime[index] - m_infTime[index];
      out += log(dgamma(incubPeriod, shape_incub, scale_incub));
    }

    if (display) cout << "log_dIncub: " << out << endl;

    return out;
}







double Household::compute_log_lik(
  std::vector<double> parameter, 
  std::vector<int> selectedParam, 
  double maxPCRDetectability,
  int display
  ) {

  // Identify index case
  size_t firstCase(0);
  auto it = std::min_element(m_infTime.begin(), m_infTime.end());
  double t0 = *it;
  for (size_t i =0; i < m_size; i++) {
    if ( m_infTime[i] == t0 ) firstCase = i;
  }

  // Log Likelihood
  double LL = 0.0;

  for (size_t i = 0; i < m_size; i++) { // All individuals

    if (display) {
      cout << i << " " << firstCase << " " << m_augmentedOnsetTime[i] << " " << m_infTime[i] << " " << m_infected[i] << endl; 
    }

    if ( i == firstCase ) {                                           		// Incubation period of 1st infected
      LL += log_dIncub(i, maxPCRDetectability, display); 

    } else if ( m_infected[i] > 0 && i != firstCase) {                 // Contribution of the secondary cases
      LL += log_S(i, t0, parameter, selectedParam, display);
    	LL += log_pInf(i, parameter, selectedParam, display);
    	LL += log_dIncub(i, maxPCRDetectability, display);

    } else {																// Non infected individuals tha contribute to the likelihood
      if (m_infected[i] ==0) LL += log_S(i, t0, parameter, selectedParam, display); 

    }
    if (display) cout << "end individual" << endl;

  }

  if (display) cout << "LL: " << LL << "\n\n";

  return LL;
}



// Probability of infection at t1
double Household::log_pInf(
  size_t curr, 
  std::vector<double> parameter, 
  std::vector<int> selectedParam, 
  int display
  ){

    // Instaneous risk of infection from community
    double alpha = parameter[0];

    // Instantaneous risk of infection from other household members
    double beta_i(0.0), beta(0.0);
    for (size_t ind=0; ind<m_size; ++ind) {
      if (curr != ind) {
        beta_i = m_instLambda[ind][curr];

        if (display) cout << "pInf: " << curr << " " << ind << " " << beta_i << endl;

        // Relative infectivity of symptomatic adult cases versus symptomatic child cases
        if ( m_infected[ind] == 1 && m_age[ind] != 2 ) beta_i *= parameter[3];

        // Contat pattern
        if (m_contactPattern == "heterogeneous") { // Reference = child-child pair in household of size 4 
          // Child (2) - Mother (0)
          if ( (m_age[curr] == 0 && m_age[ind] == 2) || (m_age[curr] == 2 && m_age[ind] == 0) ) beta_i *= 1.17;

          // Child (2) - Father (1)
          if ( (m_age[curr] == 1 && m_age[ind] == 2) || (m_age[curr] == 2 && m_age[ind] == 1) ) beta_i *= 0.55;

          // Mother (0) - Father (1)
          if ( (m_age[curr] == 0 && m_age[ind] == 1) || (m_age[curr] == 1 && m_age[ind] == 0) ) beta_i *= 1.31;
        }

        beta += beta_i;
      }
    }

    //if (display) cout << beta << endl;

    // Relative susceptibility of adults
    if ( m_age[curr] != 2 ) beta *= parameter[2];

    // Instantaneous per capita transmission rate 
    beta *= parameter[1] / (m_size / 4.0);

    if (display) cout << "log_pInf: " << log( (alpha+beta) ) << endl;

    return log( (alpha + beta) );
}




// Survival from t0 to t1
double Household::log_S(
  size_t curr, 
  double t0, 
  std::vector<double> parameter, 
  std::vector<int> selectedParam, 
  int display
  ) {

  // If the individual is not infected, it's infection
  double t1;
  if ( m_infected[curr] == 0 ) {
  	t1 = m_studyPeriod[curr];
  }else{
  	t1 = m_infTime[curr];
  }
  
  // Cumulative risk of infection from community
  double alpha = parameter[0] * (t1 - t0);

  // Cumulative risk of infection from other household members
  double beta(0.0), beta_i(0.0);
	for (size_t ind=0; ind<m_size; ++ind) {
    if (curr != ind) {
      
      beta_i = m_cumLambda[ind][curr];

      if (display) cout << "S: " << curr << " " << ind << " " << beta_i << endl;

      // Relative infectivity of symptomatic adult cases versus symptomatic child cases
      if ( m_infected[ind] == 1 && m_age[ind] != 2 ) beta_i *= parameter[3];

      // Contat pattern
      if (m_contactPattern == "heterogeneous") {
        // Child (2) - Mother (0)
        if ( (m_age[curr] == 0 && m_age[ind] == 2) || (m_age[curr] == 2 && m_age[ind] == 0) ) beta_i *= 1.17;

        // Child (2) - Father (1)
        if ( (m_age[curr] == 1 && m_age[ind] == 2) || (m_age[curr] == 2 && m_age[ind] == 1) ) beta_i *= 0.55;

        // Mother (0) - Father (1)
        if ( (m_age[ind] == 0 && m_age[curr] == 1) || (m_age[ind] == 1 && m_age[curr] == 0) ) beta_i *= 1.31;
      }

      beta += beta_i;
    }
	}

  // Relative susceptibility of adults
  if ( m_age[curr] != 2 ) beta *= parameter[2];

  // Instantaneous per capita transmission rate
  beta *= parameter[1] / (m_size / 4.0);

  if (display) cout << beta << " " << alpha << endl;
  
  if (display) cout << "log_S: " << -(alpha+beta) << endl;

	return -(alpha + beta );
}


#ifndef MCMCOBJECT__H
#define MCMCOBJECT__H

#include <random>
#include <vector>
#include "Household.h"


// Class
class McmcObject
{
public:
	McmcObject();
	McmcObject(
            //unsigned seed,
            int nIterations,
            int nAdaptiveSteps,
            std::vector<Household> data,
            std::vector<double> parameter,
            std::vector<int> selectedParameter,
            std::vector<double> rateForRandomWalk,
            int nIterTimeInfection = 1, 
            double sdLNormInfPrior = 1.0,
            int lnormInfPrior = 1,
            double sdLNormSPrior = 1.0,
            int lnormSPrior = 1,
            double maxPCRDetectability = 10.0
            );
	~McmcObject() {};

	size_t getNumbHH() const { return m_nHH; }
	int iteration() const { return m_iterations; };
	int getNIterTimeInf() const {return m_iterTimeInfection; };
	int proposedMove(int index) const {return m_numberOfMoveProposed[index]; };
	int acceptedMove(int index) const {return m_numberOfMoveAccepted[index]; };
	int proposedMoveData() const {return m_proposedMoveData; };
	int acceptedMoveData() const {return m_acceptedMoveData; };
	double globalLogLik() const { return m_globalLogLik; };
	double hhLogLik(int index) const { return m_hhLogLik[index]; };
	std::vector<double> hhLogLik() const { return m_hhLogLik; };
	double rateRandomWalk(int index) const { return m_rateForRandomWalk[index]; };
	double parameter(int index) const { return m_parameter[index]; };
	size_t nParameters() const { return m_parameter.size(); };
	
	std::string getAugmentedDataHH(int house);
	std::string getIndId(int house);

	void resetMoves();
	void initial_param_values();
	void initialize_augmented_data();
	void initial_log_lik();
	void update_parameter(int parID, double step);
	void update_augmented_infection_times();
	void update_augmented_symptom_onset();

private:
	int m_iterations;
	int m_iterTimeInfection;
	double m_adaptive_steps;
	std::mt19937_64 m_gen;
	double m_globalLogLik;
	std::vector<double> m_hhLogLik;
	int m_acceptedMoveData;
	int m_proposedMoveData;
	std::vector<int> m_numberOfMoveAccepted;
	std::vector<int> m_numberOfMoveProposed;
	std::vector<double> m_accept;
	std::vector<double> m_parameter;
	std::vector<int>  m_selectedParam;
	std::vector<double> m_rateForRandomWalk;
	double m_sdLNormInfPrior;
	int m_lnormSPrior;
	double m_sdLNormSPrior;
	int m_lnormInfPrior;
	double m_maxPCRDetectability;
	double const m_opt_acceptance = 0.25;
	std::vector<Household> m_data;
	size_t m_nHH;
};

#endif

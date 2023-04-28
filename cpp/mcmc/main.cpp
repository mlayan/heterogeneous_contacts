#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <random>
#include <chrono>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <experimental/filesystem>
#include "Household.h"
#include "McmcObject.h"

using namespace std;
namespace fs = std::experimental::filesystem;

//----------------------------------------------------------------------
// Load data
//----------------------------------------------------------------------
std::vector<Household> buildData(std::string dataFile, std::string contact_pattern)
{

    // Variable with all the data
    std::vector<Household> output;
    output.resize(0);

    // Read file
    /*
    The data file should be a space-separated table with household ids
    sorted in ascending order and numbered from 0 to n-1 households
    */
    std::ifstream infile(dataFile.c_str());
    if(infile)
    {
        cout << "=============DATA===============" << endl;
        cout << "Read data file: " << dataFile << endl;

        int indid(0), lasthhid(-1), currhhid(0), hhsize(0), isCase(0), startFollowUp(0), studyPeriod(0), age(0);
        double onsetTime = 0.0;
        double ddi = 0.0;
        int numberOfCase(0), numberOfHousehold(0), numberOfSubject(0), numberOfDay(0);
        // 0: indid, 1: hhid, 2: hhsize, 3: dds, 4: case, 5: startFollowUp, 6: studyPeriod, 7: adult
        //int nSymptomatic(0);

        Household currHH(contact_pattern);

        // while(std::getline(infile, line))
        for ( std::string line; std::getline(infile, line); )
        {
            // Create a stringstream of the current line
            std::istringstream in(line);
            
            // Store information in variables 
            in >> indid >> currhhid >> hhsize >> onsetTime >> isCase >> startFollowUp >> studyPeriod >> age >> ddi;

            // Append previous household to the list of households
            if (currhhid != lasthhid)
                {
                    if (numberOfHousehold > 0) {
                        numberOfHousehold++;
                        output.push_back(currHH);
                        currHH.newHousehold();
                        lasthhid = currhhid;
                    } else {
                        numberOfHousehold++;
                        lasthhid = currhhid;
                    }
                }

            // Update household object and information parameters
            currHH.addIndividual(indid,
                                 onsetTime,
                                 isCase,
                                 age, 
                                 startFollowUp,
                                 studyPeriod,
                                 ddi
                                 );
            numberOfSubject++;
            if (isCase != 0) numberOfCase++;
            numberOfDay = max(numberOfDay, studyPeriod);
        }

        // Add last household
        output.push_back(currHH);

        // General information
        cout << "Number of households: " << numberOfHousehold << endl;
        cout << "Number of individuals: " << numberOfSubject << endl;
        cout << "Number of cases: " << numberOfCase<< endl;
        cout << "Number of days: " << numberOfDay << "\n\n";

        return(output);

    }else{
        cout << "ERROR: Cannot open the file." << endl;
        return(output);
    }
}


//----------------------------------------------------------------------
// Run mcmc
//----------------------------------------------------------------------
void runMCMC(McmcObject mcmc,
             std::string outputFile,
             std::string augmentedDataFile,
             int pas,
             /*int pasPrinting,*/
             std::vector<int> idOfSelectedParameter,
             std::vector<std::string> paramNames
)
{
	ofstream output(outputFile.c_str());
    //ofstream augmentedData(augmentedDataFile.c_str());

	int iteration, iter, parID;
	int numberOfIteration = int(mcmc.iteration() / pas);
	int nIterTimeInfection = mcmc.getNIterTimeInf();

	// Column names
	std::string colNames="iteration logLik ";
	for (size_t p=0; p < paramNames.size(); p++) colNames += paramNames[p] + " " + paramNames[p] + "_p " + paramNames[p] + "_a ";
	colNames += "data_p data_a";
	output << colNames << endl;

	// MCMC chain
	cout << "=============MCMC===============" << endl;

    // Initial state
    mcmc.initial_param_values();
	mcmc.initialize_augmented_data(); // Initialize infection time and symptom onset for all infected individuals
	mcmc.initial_log_lik(); 
	cout << "Initial log likelihood: " << mcmc.globalLogLik() << endl;

    output << "0 " << mcmc.globalLogLik() << " "; // Log likelihood
    for (size_t i = 0; i < mcmc.nParameters(); i++)
        output << mcmc.parameter(i) << " 0 0 ";
    output << "0 0" << endl;

    // for (size_t i = 0; i < mcmc.getNumbHH(); i++)
    //     augmentedData << mcmc.getIndId(i);
    // augmentedData << endl;

    // for (size_t i = 0; i < mcmc.getNumbHH(); i++)
    //     augmentedData << mcmc.getAugmentedDataHH(i);
    // augmentedData << endl;


    // Chain
	for (iteration = 0; iteration < numberOfIteration; iteration++)
	{
	    mcmc.resetMoves();

		for (iter = 0; iter < pas; iter++)
		{
			for (size_t selectedParameter = 0; selectedParameter < idOfSelectedParameter.size(); selectedParameter++)
			{
				parID = idOfSelectedParameter[selectedParameter];
				mcmc.update_parameter(parID, iteration*pas+iter+1);
			}

			// Augmented data
			for (int i=0; i < nIterTimeInfection; i++) {
                mcmc.update_augmented_infection_times(); // Update loglik at the household level
                mcmc.update_augmented_symptom_onset();
			}
            //cout << iteration * pas + iter << " " << mcmc.globalLogLik() << endl;
		}

		// Write log likelihood, parameter values, number of proposed/accepted move per parameter in the output file
		output << (iteration + 1) * pas << " " << mcmc.globalLogLik() << " "; // Log likelihood
        //cout << (iteration + 1) * pas << " " << mcmc.globalLogLik() << " "; // Log likelihood

		for (size_t i = 0; i < mcmc.nParameters(); i++)
            output << mcmc.parameter(i) << " " << mcmc.proposedMove(i) << " " << mcmc.acceptedMove(i) << " ";

        output << mcmc.proposedMoveData() << " " << mcmc.acceptedMoveData() << endl;

        // for (size_t i = 0; i <mcmc.getNumbHH(); i++)
        //     augmentedData << mcmc.getAugmentedDataHH(i);
        // augmentedData << endl;
	}

	output.close();
    // augmentedData.close();
}



//----------------------------------------------------------------------
// Main function
//----------------------------------------------------------------------
int main(int argc, char **argv)
{
    
    // Arguments passed to main
    std::string nSim = argv[1];
    std::string data_contact = argv[2];
    std::string inference_contact = argv[3];
    std::string param_combination = argv[4];

    double maxPCRDetectability = 10.0;     
    double sdrInf = 1.0; 
    double sdrS = 1.0; 
    int lnormInfPrior = 1; //0: uniform prior, 1: lognormal prior
    int lnormSPrior = 1; //0: uniform prior, 1: lognormal prior
    std::string chainID = "1";

    //==========Model parameters==========
    // Initial values
    std::vector<std::string> parameterNames = {"alpha", "beta", "delta", "rSC", "rInfAA", "rInfAC", "rInfSC"};
    //                                            0        1        2      3        4        5           6
    
    int numberOfParameters = parameterNames.size();
    std::vector<double> parameter(numberOfParameters);

    // Parameters to infer according to model
    std::vector<int> selectedParameter(numberOfParameters, 1);
    
    //parameter[0] = 0.001;
    //selectedParameter[0] = 0;
    
    //parameter[1] = 0.6;
    //selectedParameter[1] = 0;

    parameter[2] = 1.0;
    selectedParameter[2] = 0;

    parameter[4] = 1;
    selectedParameter[4] = 0;
    
    parameter[5] = 1;
    selectedParameter[5] = 0;

    //parameter[6] = 1;
    //selectedParameter[5] = 0;

	int parameterNumber;
	std::vector<int> idOfSelectedParameter(0);
    for(parameterNumber=0; parameterNumber<numberOfParameters; parameterNumber++)
    {
        if(selectedParameter[parameterNumber]==1) idOfSelectedParameter.push_back(parameterNumber);
    }

    cout << "ID of selected parameters: "; 
    for (auto i = idOfSelectedParameter.begin(); i != idOfSelectedParameter.end(); ++i)
    	std::cout << parameterNames[*i] << ' ';
    cout << endl;


    //==========MCMC parameters==========
    /*    
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed = d.count();
    unsigned seed;
    if ( chainID == "1" ) seed=20201208;
    if ( chainID == "2" ) seed=20210415;
    if ( chainID == "3" ) seed=20210601;
    cout << "Seed: " << seed << endl;
    */

    int pas = 40; 
    int numberOfIteration = 50000;
    //int nAdaptiveSteps = pas * numberOfIteration / 4;   // Doesn't work much...
    int nAdaptiveSteps = 0;
    int numberOfIterationTimeInfection = 1;

    // Variance of random walk
    std::vector<double> rateForRandomWalk(numberOfParameters);

    if (param_combination == "flu") {
        rateForRandomWalk = { 1.0,  0.2,  0.3,   0.3,    1.5,    0.5,    0.25};
        //                  alpha   beta delta   rSC    rInfAA  rInfAC  rInfSC   
    }

    if (param_combination == "covid") {
        rateForRandomWalk = { 1.0,  0.2,  0.3,   0.22,    1.5,    0.5,    0.22};
        //                  alpha   beta delta   rSC    rInfAA  rInfAC  rInfSC
    }
    
    

    //==========Output files==========
    //Paths
    std::string condition;
    std::string dataFile, outputFile, augmentedDataFile; 


    if (nSim == "recover") {
        //Data file
        // File structure :         0: indid, 1: hhid, 2: hhsize, 3: dds, 4: case, 5: studyPeriod, 6: adult
        dataFile="/pasteur/appa/homes/maylayan/MMMICovid/Contacts/data/recover.txt";        //Input file 
        outputFile="/pasteur/appa/homes/maylayan/MMMICovid/Contacts/results/mcmc_" + nSim + "_" + inference_contact + ".txt"; //Output file
        augmentedDataFile="/pasteur/appa/homes/maylayan/MMMICovid/Contacts/results/augmented_data_" + nSim + "_" + inference_contact + ".txt"; //Output file


    } else {

        if (data_contact == "homogeneous" && inference_contact == "homogeneous") condition = "homo_homo";
        if (data_contact == "homogeneous" && inference_contact == "heterogeneous") condition = "homo_hetero";
        if (data_contact == "heterogeneous" && inference_contact == "homogeneous") condition = "hetero_homo";
        if (data_contact == "heterogeneous" && inference_contact == "heterogeneous") condition = "hetero_hetero";

        std::string pathData="/pasteur/appa/homes/maylayan/MMMICovid/Contacts/data/" + param_combination + "/" + data_contact + "/";
        std::string pathOutput="/pasteur/appa/homes/maylayan/MMMICovid/Contacts/results/" + param_combination + "/" + condition + "/";
        
        fs::create_directories(pathOutput);
        
        // Data file
        // File structure :         0: indid, 1: hhid, 2: hhsize, 3: dds, 4: case, 5: studyPeriod, 6: adult
        dataFile=pathData + "data_" + nSim + ".txt";                         //Input file
        outputFile=pathOutput + "mcmc_" + nSim + ".txt";                    //Output file
        augmentedDataFile=pathOutput + "augmented_data_" + nSim + ".txt";                    //Output file
    }

    
    cout << "Input file: " << dataFile << endl;
    cout << "Output file: " << outputFile << "\n\n";
    
    
    //==========Build data==========
    // Load data
    std::vector<Household> hhData = buildData(dataFile, inference_contact);

    // Initialize MCMC object
    McmcObject mcmc(
        //seed, 
        numberOfIteration, 
        nAdaptiveSteps,
        hhData, 
        parameter, 
        selectedParameter, 
        rateForRandomWalk, 
        numberOfIterationTimeInfection, 
        sdrInf, 
        lnormInfPrior,
        sdrS,
        lnormSPrior,
        maxPCRDetectability          
        );

    //==========MCMC==========
    runMCMC(
        mcmc,
        outputFile,
        augmentedDataFile,
        pas,
        idOfSelectedParameter,
        parameterNames
        );

    return 0;
}

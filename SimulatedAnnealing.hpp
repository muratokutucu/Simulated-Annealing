#ifndef SIMULATEDANNEALING_H
#define SIMULATEDANNEALING_H

#include <utility> //pair, make_pair
#include <vector>
#include <random>
#include <ctime>

class SimulatedAnnealing{
private:
	int m_numberOfParameters;
	std::vector<std::pair<double, double>> m_constraints;

public:

	/**
	 * Constructor. Takes a vector containing the range of each parameter.
	 * @param constraints - Range of each parameter of the problem. Example : {[xmin, xmax], [ymin, ymax], ...}
	 */
	SimulatedAnnealing(std::vector<std::pair<double, double>> constraints){
		m_constraints = constraints;
		m_numberOfParameters = constraints.size();
	}

	double sigmoideFunction(double T, double T_MAX){
		double x = -10+(T/T_MAX)*20;
		return 1.0/(1.0 + exp(-x));
	}

	/**
	 * Returns a neighbor solution.
	 * @param steps - half length of search zone's. 
	 * @param from - value of each parameter of current solution.
	 */
	std::vector<double> getRandomNeighbor(std::vector<double> steps, std::vector<double> from){
		std::vector<double> res(m_numberOfParameters);
		for(int i = 0; i < m_numberOfParameters; i++){
			double step = steps[i];
			double expectedMin = from[i]-step;
			double expectedMax = from[i]+step;

			std::pair<double, double> range = m_constraints[i];
			double low = range.first, high = range.second;

			double min = (expectedMin < low) ? low : expectedMin;
			double max = (expectedMax > high) ? high : expectedMax;
		    //std::cout << "from " << from[i] << " min " << min << " max " << std::endl;
	        std::uniform_real_distribution<double> dist(min, max); //Type of random number distribution
	        std::mt19937 rng; //Mersenne Twister: Good quality random number generator
	        rng.seed(std::random_device{}()); //Initialize with non-deterministic seeds

			res[i] = dist(rng); //neighboor
		}
		return res;
	}

	/**
	 * Apply the simulated annealing algorithm in order to find the global optimum of the function 'evaluateScore' (passed as argument).
	 * @param initialTemperature - initial temperature. Determined empirically. 
	 * @param alpha - Temperature's decrease speed. Recommended value : 0.99
	 */
	template<typename Function>
	std::pair<std::vector<double>, std::vector<double>> run(Function evaluateScore, uint64_t iteration, double alpha, double initialTemperature){
		double T = initialTemperature;
		double T_MAX = initialTemperature;
		// definition of step for each parameter
		std::vector<double> steps(m_numberOfParameters);
		for(int i = 0; i < m_numberOfParameters; i++){
			std::pair<double, double> range = m_constraints[i];
			double low = range.first, high = range.second;
			steps[i] = (abs(high) - abs(low));
		}

		std::vector<double> best_solution(m_numberOfParameters), new_solution(m_numberOfParameters), solution(m_numberOfParameters);
		double best_score, new_score, score;
		
		//generate a first random solution
		for(int i = 0; i < m_numberOfParameters; i++){
			std::pair<double, double> range = m_constraints[i];
			double min = range.first, max = range.second;
			std::uniform_real_distribution<double> dist(min, max);  
       		std::mt19937 rng; 
			solution[i] = dist(rng);
		}

		score = evaluateScore(solution);
		best_solution = solution;
		best_score = score;

		//generator of 'double' between 0 and 1
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        std::mt19937 rng; 
        rng.seed(std::random_device{}()); 

		while (iteration--){
			new_solution = getRandomNeighbor(steps, solution); 
			new_score = evaluateScore(new_solution);

			if(new_score > score){
				solution = new_solution;
				score = new_score;
				if(new_score > best_score){ best_solution = new_solution; best_score = new_score;}
			} else if(dist(rng) < sigmoideFunction(T, T_MAX)){ //exp(-abs(score-new_score)/T);
				solution = new_solution;
				score = new_score;
			}

			T *= alpha; //reduce temperature
		}		
		return {best_solution, solution};
	}
};

#endif //SIMULATEDANNEALING_H

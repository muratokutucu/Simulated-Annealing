#ifndef SimulatedAnnealing_SimulatedAnnealing_H
#define SimulatedAnnealing_SimulatedAnnealing_H

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

	/**
	 * Apply the simulated annealing algorithm in order to find the global optimum of the function 'evaluateScore' (passed as argument).
	 * @param initialTemperature - initial temperature. Determined empirically. 
	 * @param alpha - Temperature's decrease speed. Recommended value : 0.99
	 */
	template<typename Function>
	std::vector<double> run(double initialTemperature, double alpha, Function evaluateScore){
		double T = initialTemperature;
		// definition of step for each parameter
		std::vector<double> steps(m_numberOfParameters);
		for(int i = 0; i < m_numberOfParameters; i++){
			std::pair<double, double> range = m_constraints[i];
			double low = range.first, high = range.second;
			steps[i] = (high - low);
		}


		std::vector<double> best_solution(m_numberOfParameters), new_solution(m_numberOfParameters), solution(m_numberOfParameters);
		double best_score, new_score, score;
		
		//generate a first random solution
		for(int i = 0; i < m_numberOfParameters; i++){
			std::pair<double, double> range = m_constraints[i];
			double low = range.first, high = range.second;
			solution[i] = (abs(high) - abs(low))/2;
		}
		score = evaluateScore(solution);

		//consider the first solution as the best for now
		best_solution = solution;
		best_score = score;

		//Type of random number distribution
	        std::uniform_real_distribution<double> dist(0.0, 1.0);  //(min, max)
	        //Mersenne Twister: Good quality random number generator
	        std::mt19937 rng; 
	        //Initialize with non-deterministic seeds
	        rng.seed(std::random_device{}()); 

		while (T > 1){
			new_solution = getRandomNeighbor(steps, solution);
			new_score = evaluateScore(new_solution);

			if(new_score > score){
				solution = new_solution;
				score = new_score;
				//we keep the best
				if(new_score > best_score){
					best_solution = new_solution;
					best_score = new_score;
				}
			} else if(dist(rng) < exp((score-new_score)/T)){
				solution = new_solution;
				score = new_score;
			}

			//reduce temperature
			T *= alpha;

			//reduce each step/*
			for(int i = 0; i < m_numberOfParameters; i++){
				steps[i] *= alpha;
			}
		}		
		return best_solution;
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
			double max = (expectedMax > max) ? max : expectedMax;
			
			//Type of random number distribution
		        std::uniform_real_distribution<double> dist(min, max);  //(min, max)
		        //Mersenne Twister: Good quality random number generator
		        std::mt19937 rng; 
		        //Initialize with non-deterministic seeds
		        rng.seed(std::random_device{}()); 

			res[i] = dist(rng); //neighboor
		}
		return res;
	}
};

#endif //SimulatedAnnealing_SimulatedAnnealing_H

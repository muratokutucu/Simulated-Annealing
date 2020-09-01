#ifndef TSP_H
#define TSP_H

#include <utility> //pair, make_pair
#include <vector>
#include <random>
#include <ctime>
#include <algorithm> //for std::remove
#include <string>

struct city{
	std::string name;
	double x, y;
};

class TSP{
private:
	int m_numberOfCity;
	std::vector<city> m_cities; //current_max_number of supported city = 9

public:

	/**
	 * Constructor. Takes a vector containing the range of each parameter.
	 * @param constraints - Range of each parameter of the problem. Example : {[xmin, xmax], [ymin, ymax], ...}
	 */
	TSP(std::vector<city> cities){
		m_cities = cities;
		m_numberOfCity = cities.size();
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
	double getRandomNeighbor(double step, double from){
		double expectedMin = from-step;
		double expectedMax = from+step;

		double low = 0, high = m_numberOfCity-1;

		double min = (expectedMin < low) ? low : expectedMin;
		double max = (expectedMax > high) ? high : expectedMax;
	    //std::cout << "from " << from[i] << " min " << min << " max " << std::endl;
        std::uniform_real_distribution<double> dist(min, max); //Type of random number distribution
        std::mt19937 rng; //Mersenne Twister: Good quality random number generator
        rng.seed(std::random_device{}()); //Initialize with non-deterministic seeds

		return dist(rng); //neighboor
	}

	/**
	 * Apply the simulated annealing algorithm in order to find the global optimum of the function 'evaluateScore' (passed as argument).
	 * @param initialTemperature - initial temperature. Determined empirically. 
	 * @param alpha - Temperature's decrease speed. Recommended value : 0.99
	 */
	template<typename Function>
	std::pair<std::vector<city>, std::vector<city>> run(Function energy, uint64_t iteration, double alpha, double initialTemperature){
		double T = initialTemperature;
		double T_MAX = initialTemperature;
		// definition of step for each parameter
		double step = m_numberOfCity/4.0;

		double best_solution, new_solution, solution;
		double best_e, new_e, e;
		
		//generate a first random solution
		std::uniform_real_distribution<double> distSolution(0, m_numberOfCity-1);  
   		std::mt19937 rngSolution; 
   		rngSolution.seed(std::random_device{}()); //Initialize with non-deterministic seeds
		solution = distSolution(rngSolution);

		e = energy(transform(solution));
		best_solution = solution;
		best_e = e;
		
		//generator of 'double' between 0 and 1
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        std::mt19937 rng; 
        rng.seed(std::random_device{}()); 

		while (iteration--){
			new_solution = getRandomNeighbor(step, solution); 
			new_e = energy(transform(new_solution));
			//std::cout << "x " << new_solution[0] << " step " << steps[0] << std::endl;
			if(new_e < e){
				solution = new_solution;
				e = new_e;
				if(new_e < best_e){ best_solution = new_solution; best_e = new_e;}
			} else if(dist(rng) < sigmoideFunction(T, T_MAX)){ //exp(-abs(score-new_score)/T);
				solution = new_solution;
				e = new_e;
			}

			T *= alpha; //reduce temperature
		}
		//we have to transform a vector<int> into a vector<city>	
		return {transform(best_solution), transform(solution)};
	}

	std::vector<city> transform(double sol){
		//std::cout << "sol" << sol << std::endl;
		int total_nb_of_city = m_cities.size();
		std::vector<int> path(total_nb_of_city+1);
		std::vector<int> cities(total_nb_of_city-1);

		int nbOfPossibleCity = total_nb_of_city-1;

		//we consider the departure as the city 0
		path[0] = 0;
		for(int i = 0; i < nbOfPossibleCity; i++){
			cities[i] = i+1;
			//std::cout << cities[i] << " * ";
		}
		//std::cout << std::endl;

		int i = 1;
		double x = sol;
		int iteration = 0;

		while(nbOfPossibleCity > 0){
			int city_i = (int)x;
			nbOfPossibleCity--;
			
			path[i] = cities[city_i];
			//std::cout << "path" << path[i] << std::endl;
			double remaining = x - city_i;
			x = remaining * nbOfPossibleCity;
			//x = 0 + remaining * nbOfPossibleCity; //remaining
			//std::cout << " x " << x << " " << city_i << " " << cities.size() << std::endl;
 			cities.erase(std::remove(cities.begin(), cities.end(), cities[city_i]), cities.end()); //remove added city from possible ones

			i++;
		}
		//std::cout << std::endl;
		path[path.size()-1] = 0;

		//transform
		std::vector<city> res(path.size());
		for(int i = 0; i < path.size(); i++){
			res[i] = m_cities[path[i]];
			//std::cout << res[i].name << " ";
		}
		//std::cout << std::endl;
		return res;
	}

};

#endif //SIMULATEDANNEALING_H

#include <iostream>
#include <tuple>
#include <vector>
#include <utility> //pair, make_pair
#include <random>
#include <cmath>
#include <ctime>
#include "SimulatedAnnealing.hpp"

using namespace std;

/**
 * Constructor. Takes a vector containing the range of each parameter.
 * @param constraints - Range of each parameter of the problem. Example : {[xmin, xmax], [ymin, ymax], ...}
 */
SimulatedAnnealing::SimulatedAnnealing(vector<pair<double, double>> constraints){
	m_constraints = constraints;
	m_numberOfParameters = constraints.size();
}

/**
 * Returns a neighbor solution.
 * @param steps - half length of search zone's. 
 * @param from - value of each parameter of current solution.
 */
vector<double> SimulatedAnnealing::getRandomNeighbor(vector<double> steps, vector<double> from){
	vector<double> res(m_numberOfParameters);
	for(int i = 0; i < m_numberOfParameters; i++){
		double step = steps[i];
		double expectedMin = from[i]-step;
		double expectedMax = from[i]+step;
		pair<double, double> range = m_constraints[i];
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

/**
 * Apply the simulated annealing algorithm in order to find the global optimum of the function 'evaluateScore' (passed as argument).
 * @param initialTemperature - initial temperature. Determined empirically. 
 * @param alpha - Temperature's decrease speed. Recommended value : 0.99
 */
template<typename Function>
vector<double> SimulatedAnnealing::run(double initialTemperature, double alpha, Function evaluateScore){
	double T = initialTemperature;
	// definition of step for each parameter
	vector<double> steps(m_numberOfParameters);
	for(int i = 0; i < steps.size(); i++){
		pair<double, double> range = m_constraints[i];
		double low = range.first, high = range.second;
		steps[i] = (high - low);
	}


	vector<double> best_solution(m_numberOfParameters), new_solution(m_numberOfParameters), solution(m_numberOfParameters);
	double best_score, new_score, score;
	
	//generate a first random solution
	for(int i = 0; i < m_numberOfParameters; i++){
		pair<double, double> range = m_constraints[i];
		double low = range.first, high = range.second;
		solution[i] = (abs(high) - abs(low))/2;
	}
	score = evaluateScore(solution);

	//consider the first solution as the best for now
	best_solution = solution;
	best_score = score;

	//Type of random number distribution
    uniform_real_distribution<double> dist(0.0, 1.0);  //(min, max)
    //Mersenne Twister: Good quality random number generator
    mt19937 rng; 
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
 * Evaluates the score of the problem that we want to find the global optimum.
 * @param values - values of each parameter of the problem.
 */
double evaluateScore(vector<double> values){
	double x = values[0];
	double y = values[1];	
	return sin(x)*cos(y)*x;
}

//g++ -std=c++11 -O3 SimulatedAnnealing.cpp -o exe
//./exe

//vector<double> result = s.run(100000000000, 0.99995, evaluateScore); -> 4,277 s
//après amélioration : -> 
int main(void){
	pair<double, double> x_range = make_pair(-7, 7);
	pair<double, double> y_range = make_pair(-7, 7);
	vector<pair<double, double>> constraints {x_range, y_range};
	
	SimulatedAnnealing s(constraints);
	vector<double> result = s.run(100000000000, 0.99995, evaluateScore);
	
	cout << "Best result : ";
	for(int i = 0; i < result.size(); i++){
		cout << result[i] << " ";
	}
	cout << endl;
	return 0;
}
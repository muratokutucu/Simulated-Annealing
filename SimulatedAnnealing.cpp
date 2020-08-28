#include <iostream>
#include <tuple>
#include <vector>
#include <utility> //pair, make_pair
#include <random>
#include <cmath>
#include <ctime>
#include "sa.hpp"

using namespace std;

SimulatedAnnealing::SimulatedAnnealing(double temperature, double alpha, vector< pair<double, double> > constraints){
	m_temperature = temperature;
	m_alpha = alpha;
	m_constraints = constraints;
}

//	std::vector<double> m_parametersWithOptimizedValues; 
//	std::vector< std::pair<double, double> > m_parametersWithConstraints;
vector<double> SimulatedAnnealing::getRandomNeighbor(vector<double> steps, vector<double> from){
	vector<double> res(steps.size());
	for(int i = 0; i < res.size(); i++){
		double step = steps[i];
		pair<double, double> range = m_constraints[i];
		double low = range.first, high = range.second;

		double min = (from[i]-step < low) ? low : from[i]-step;
		double max = (from[i]+step > max) ? max : from[i]+step;
		
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

double SimulatedAnnealing::evaluateScore(std::vector<double> values){
	//here, in this example, the vector is only composed of the value of the parameter x
	//let's say that the function that we want to optimize iz f(x) = x^2 and that we want the value of x that provides the greater f(x).
	//So, the higher will be f(x) the higher must be the score. An obvious solution is simply to return the value of x^2.
	double x = values[0];
	double y = values[1];	
	return sin(x)*cos(y)*x;
}

vector<double> SimulatedAnnealing::run(){
	double T = m_temperature;
	// definition of step for each parameter
	vector<double> steps(m_constraints.size());
	for(int i = 0; i < steps.size(); i++){
		pair<double, double> range = m_constraints[i];
		double low = range.first, high = range.second;
		steps[i] = (high - low);
	}


	vector<double> best_solution(m_constraints.size()), new_solution(m_constraints.size()), solution(m_constraints.size());
	double best_score, new_score, score;
	
	//generate a first random solution
	for(int i = 0; i < m_constraints.size(); i++){
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
		T *= 0.99;

		//reduce each step/*
		for(int i = 0; i < steps.size(); i++){
			steps[i] *= 0.99;
		}
	}		
	return best_solution;
}

int main(void){
	pair<double, double> x_range = make_pair(-7, 7);
	pair<double, double> y_range = make_pair(-7, 7);
	vector< pair<double, double> > constraints {x_range, y_range};
	
	SimulatedAnnealing s(100.0, 0.99, constraints);
	vector<double> result = s.run();
	
	cout << "Best result : ";
	for(int i = 0; i < result.size(); i++){
		cout << result[i] << " ";
	}
	cout << endl;
	return 0;
}
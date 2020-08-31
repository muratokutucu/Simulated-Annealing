#include <iostream>
#include <utility>
#include <vector>
#include <cmath>

#include "SimulatedAnnealing.hpp"

using namespace std;

/**
 * Evaluates the score of the problem that we want to find the global optimum.
 * @param values - values of each parameter of the problem.
 */
double evaluateScore(vector<double> values){
	double x = values[0];
	double y = values[1];	
	//ackley function : https://en.wikipedia.org/wiki/Ackley_function
	double pi = 3.14, e = 2.718;
	double f = -20*exp(-0.2*sqrt(0.5*(x*x+y*y)))-exp(0.5*(cos(2*pi*x)+cos(2*pi*y)))+e+20;
	return 1.0/f;
}

typedef pair<double, double> pdd;

int main(void){
	pdd x_range = make_pair(-2, 2);
	pdd y_range = make_pair(-2, 2);
	vector<pdd> constraints {x_range, y_range};
	
	SimulatedAnnealing s(constraints);
	pair<vector<double>, vector<double>> result = s.run(evaluateScore, 100000, 0.9995, 10);
	
	cout << "Best result  : ";
	for(int i = 0; i < result.first.size(); i++){
		cout << result.first[i] << " ";
	}
	cout << endl;


	cout << "Final result : ";
	for(int i = 0; i < result.second.size(); i++){
		cout << result.second[i] << " ";
	}
	cout << endl;

	return 0;
}
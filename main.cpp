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
	return sin(x)*cos(y)*x;
}

//g++ -std=c++11 -O3 main.cpp SimulatedAnnealing.cpp -o exe
//./exe
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
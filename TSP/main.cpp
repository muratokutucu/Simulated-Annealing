#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <string>
#include "TSP.hpp"

using namespace std;

//g++ -std=c++11 -O3 main.cpp SimulatedAnnealing.cpp -o exe
//./exe

/**
 * Evaluates the score of the problem that we want to find the global optimum.
 * @param values - values of each parameter of the problem.
 */
double energy(vector<city> path){
	double score = 0;
	double total_distance = 0;
	for(int i = 0; i < path.size()-1; i++){
		//calculate distance between two consecutive cities
		double x1 = path[i].x, x2 = path[i+1].x;
		double y1 = path[i].y, y2 = path[i+1].y;
		double a = abs(x2-x1), b = abs(y2-y1); 
		double distance = sqrt(a*a + b*b);
		total_distance += distance;
	}
	return total_distance;
}


int main(void){
	//ci corresponds to coordinates of cities in a 2d plan
	city c0 = {"Tokyo", 0, 0};
	city c1 = {"Berlin", 0, 1};
	city c2 = {"Paris", 0, 2};
	city c3 = {"New York", 1, 0};
	city c4 = {"Tunis", 1, 1};
	city c5 = {"Ankara", 1, 2};
	city c6 = {"Mexico City", 2, 0};
	city c7 = {"Jakarta", 2, 1};
	city c8 = {"Lima", 2, 2};
	city c9 = {"Prague", 3, 0};
	city c10 = {"Rome", 3, 1};
	city c11 = {"Vienna", 3, 2};
	vector<city> cities {c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11};
	
	TSP s(cities);
	pair<vector<city>, vector<city>> result = s.run(energy, 1000, 0.99995, 10);
	
	cout << "Best result  : ";
	for(int i = 0; i < result.first.size(); i++){
		cout << result.first[i].name << " > ";
	}
	cout << endl;


	cout << "Final result : ";
	for(int i = 0; i < result.second.size(); i++){
		cout << result.second[i].name << " > ";
	}
	cout << endl;
	return 0;
} 
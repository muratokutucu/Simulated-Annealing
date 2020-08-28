#ifndef SimulatedAnnealing_SimulatedAnnealing_H
#define SimulatedAnnealing_SimulatedAnnealing_H

#include <tuple>
#include <utility> //pair, make_pair
#include <string>
#include <vector>

class SimulatedAnnealing{
private:
	int m_numberOfParameters;
	std::vector<std::pair<double, double>> m_constraints;

public:
	SimulatedAnnealing(std::vector<std::pair<double, double>> constraints);
	std::vector<double> getRandomNeighbor(std::vector<double> steps, std::vector<double> from);
	template<typename Function>
	std::vector<double> run(double initialTemperature, double alpha, Function evaluateScore); 
};

#endif 
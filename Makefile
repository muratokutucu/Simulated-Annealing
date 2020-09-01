e: main.cpp SimulatedAnnealing.hpp
	g++ -std=c++11 -O3 -g main.cpp -o e

clean:
	rm e
	rm -r e.dSYM
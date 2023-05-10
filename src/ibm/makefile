EXE=social_learning_network.exe
CPP=individual.cpp network_main.cpp simulation.cpp patch.cpp
HPP=individual.hpp simulation.hpp patch.hpp 

CXX=g++
CXXFLAGS=-O3 -Wall -std=c++17

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

clean :
	rm -rf $(EXE)


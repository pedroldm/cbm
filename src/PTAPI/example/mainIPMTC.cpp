/*
 * main.cpp
 * 
 * Copyright 2021 André Luis <André Luis@DESKTOP-HDL2CBS>
 *
 * Compilation Instructions 
 * g++ mainIPMTC.cpp ../include/*.h ./IPMTC.cpp -std=c++20 -Wshadow -Wall -o mainIPMTC -O3 -Wno-unused-result -lpthread
 *
 * Example Execution
 * ./mainIPMTC ../instances/IPMTC/instanceLarge1_m\=3_n\=50_l\=30_c\=10_s\=0.PMTC --TEMP_INIT 0.2 --TEMP_FIM 1 --N_REPLICAS 16 --MCL 400 --PTL 1000 --TEMP_DIST 1 --TYPE_UPDATE 2 --TEMP_UPDATE 35000
 *
 */
#include <cstdlib>
#include "../include/ExecTime.h"
#include "./IPMTC.h"
#include "../include/PT.h"


using namespace std;

int main(int argc, char* argv[])
{
	//varibles
	float tempIni = 0.01;
	float tempfim = 2.0;
	int tempN = 10;
	int MCL = 0;
	int PTL = 2;	
	int tempUp = 50;
	int tempD = 1;
	int uType = 0;
	int thN = thread::hardware_concurrency();	
	vector<std::string> arguments(argv + 1, argv + argc);	
	
	// Instance file name
	std::string fn = arguments[0];
	
	// Read arguments
	for(unsigned int i=1; i<arguments.size(); i+=2)
	{
		            
        if(arguments[i]== "--TEMP_FIM")
            tempfim = std::stof(arguments[i+1]);
        else if(arguments[i]== "--TEMP_INIT")
            tempIni = std::stof(arguments[i+1]);
        else if(arguments[i]== "--N_REPLICAS")
            tempN = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--MCL")
            MCL  = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--PTL")
            PTL = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--TEMP_DIST")
            tempD = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--TYPE_UPDATE")
            uType = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--TEMP_UPDATE")
            tempUp = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--THREAD_USED")
            thN = std::stoi(arguments[i+1]);
    }
	
	// Create SSP object
	IPMTC* prob = new IPMTC(fn);
	
	// Create and start PT 
	PT<solIPMTC> algo(tempIni,tempfim,tempN,MCL,PTL,tempD,uType,tempUp);
	ExecTime et;
	solIPMTC sol = algo.start(thN, prob);
	// Divides between machines
	sol = prob->setCut(sol); 
	et.printTime();
	std::cout<<sol.evalSol<<"\n";

	return 0;
}


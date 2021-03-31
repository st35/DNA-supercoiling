#include <iostream>
#include <fstream>
#include <ctime>
#include <random>
#include <boost/numeric/odeint.hpp>
#include <mpi.h>
#include "../../code/Constants.hpp"
#include "../../code/Model_Dynamics.hpp"
#include "Simulation_Setups.hpp"

void Run_Simulation(double promoterarg, int argid, int world_rank)
{
	double force = 1.0;
	double gene_length0 = 5300.0, gene_length1 = 5300.0, spacer = 2500.0;
	int config = 1, clamp0_flag = 1, clamp1_flag = 1;
	double T = 100.0;
	int controlgene = 0;
	double promoter0 = 10.0, promoter1 = 0.0, topoisomerase = 1.0;
	int stepwiseflag = 1, runid = 0;
	std::string outputfolder = "outputfiles";
	int fileflag = 0;

	clamp0_flag = 1;
	clamp1_flag = 1;
	topoisomerase = 1.0;

	gene_length0 = 5300.0;
	gene_length1 = 0.0;
	spacer = 0.0;
	config = 1;
	promoter0 = promoterarg;
	promoter1 = 0.0;
	controlgene = 0;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, stepwiseflag, runid, outputfolder, fileflag);
	}

//	gene_length0 = 0.0;
//	gene_length1 = 5300.0;
//	spacer = 0.0;
//	config = 1;
//	promoter0 = 0.0;
//	promoter1 = promoterarg;
//	controlgene = 1;
//
//	for(int i = 0; i < 16; i++)
//	{
//		outputfolder = "outputfiles/Rev_RUN_" + std::to_string(argid);
//		runid = 16*world_rank + i;
//		Gillespie_Simulation(gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, stepwiseflag, runid, outputfolder, fileflag);
//	}

	clamp0_flag = 0;
	clamp1_flag = 0;

	gene_length0 = 5300.0;
	gene_length1 = 0.0;
	spacer = 0.0;
	config = 1;
	promoter0 = promoterarg;
	promoter1 = 0.0;
	controlgene = 0;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/Linear_RUN_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, stepwiseflag, runid, outputfolder, fileflag);
	}

//	gene_length0 = 0.0;
//	gene_length1 = 5300.0;
//	spacer = 0.0;
//	config = 1;
//	promoter0 = 0.0;
//	promoter1 = promoterarg;
//	controlgene = 1;
//
//	for(int i = 0; i < 16; i++)
//	{
//		outputfolder = "outputfiles/Linear_Rev_RUN_" + std::to_string(argid);
//		runid = 16*world_rank + i;
//		Gillespie_Simulation(gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, stepwiseflag, runid, outputfolder, fileflag);
//	}

	return;
}

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	generator = std::mt19937(std::time(NULL) + world_rank);

	double promoterarg = std::stod(argv[1]);
	int argid = std::stod(argv[2]);

	Run_Simulation(promoterarg, argid, world_rank);

	MPI_Finalize();

	return(0);
}

#include <iostream>
#include <fstream>
#include <ctime>
#include <random>
#include <boost/numeric/odeint.hpp>
#include <mpi.h>
#include "../../code/Constants.hpp"
#include "../../code/Model_Dynamics.hpp"
#include "Simulation_Setups.hpp"

void Run_Simulation(double spacerarg, int spacerid, int world_rank)
{
	double force = 1.0;
	double gene_length0 = 5300.0, gene_length1 = 5300.0, spacer = 2500.0;
	int config = 1, clamp0_flag = 1, clamp1_flag = 1;
	double T = 50.0;
	int controlgene = 0;
	double promoter0 = 10.0, promoter1 = 0.0, topoisomerase = 1.0;
	int stepwiseflag = 1, runid = 0, argid;
	std::string outputfolder = "outputfiles";
	int fileflag = 0;

	clamp0_flag = 1;
	clamp1_flag = 1;

	topoisomerase = 1.0;
	gene_length0 = 5300.0;
	gene_length1 = 5300.0;
	spacer = spacerarg;
	T = 10.0;

	config = 0;
	promoter1 = 1.0e-2;
	promoter0 = 0.0;
	controlgene = 1;
	argid = 0;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_" + std::to_string(spacerid) + "_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, stepwiseflag, runid, outputfolder, fileflag);
	}

	config = 0;
	promoter1 = 1.0e-2;
	promoter0 = 1.0;
	controlgene = 1;
	argid = 1;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_" + std::to_string(spacerid) + "_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, stepwiseflag, runid, outputfolder, fileflag);
	}

	config = 1;
	promoter1 = 1.0e-2;
	promoter0 = 0.0;
	controlgene = 1;
	argid = 2;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_" + std::to_string(spacerid) + "_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, stepwiseflag, runid, outputfolder, fileflag);
	}

	config = 1;
	promoter1 = 1.0e-2;
	promoter0 = 1.0;
	controlgene = 1;
	argid = 3;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_" + std::to_string(spacerid) + "_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, stepwiseflag, runid, outputfolder, fileflag);
	}

	config = 2;
	promoter1 = 1.0e-2;
	promoter0 = 0.0;
	controlgene = 1;
	argid = 4;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_" + std::to_string(spacerid) + "_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, stepwiseflag, runid, outputfolder, fileflag);
	}

	config = 2;
	promoter1 = 1.0e-2;
	promoter0 = 1.0;
	controlgene = 1;
	argid = 5;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_" + std::to_string(spacerid) + "_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, stepwiseflag, runid, outputfolder, fileflag);
	}

	return;
}

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	double spacerarg = std::stod(argv[1]);
	double spacerid = std::stod(argv[2]);

	generator = std::mt19937(std::time(NULL) + world_rank);

	Run_Simulation(spacerarg, spacerid, world_rank);

	MPI_Finalize();

	return(0);
}

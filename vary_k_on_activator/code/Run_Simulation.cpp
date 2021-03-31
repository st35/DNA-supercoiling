#include <iostream>
#include <fstream>
#include <ctime>
#include <random>
#include <boost/numeric/odeint.hpp>
#include <mpi.h>
#include "../../code/Constants.hpp"
#include "../../code/Model_Dynamics.hpp"
#include "Simulation_Setups.hpp"

void Run_Simulation(double TFarg, int argid, int world_rank)
{
	double force = 1.0;
	double gene_length0 = 5300.0, gene_length1 = 5300.0, spacer = 2500.0;
	int config = 1, clamp0_flag = 1, clamp1_flag = 1;
	double T = 100.0;
	int controlgene = 0;
	double promoter0 = 10.0, promoter1 = 0.0, topoisomerase = 1.0, TF = 0.0;
	int stepwiseflag = 1, runid = 0;
	std::string outputfolder = "outputfiles";
	int fileflag = 0;

	TF = TFarg;

	clamp0_flag = 1;
	clamp1_flag = 1;

	gene_length0 = 5300.0;
	gene_length1 = 0.0;
	spacer = 0.0;
	config = 1;
	promoter0 = 10.0;
	promoter1 = 0.0;
	controlgene = 0;

	topoisomerase = 1.0e-2;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_0_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, TF, stepwiseflag, runid, outputfolder, fileflag);
	}

	topoisomerase = 1.0e-1;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_1_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, TF, stepwiseflag, runid, outputfolder, fileflag);
	}

	topoisomerase = 1.0;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_2_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, TF, stepwiseflag, runid, outputfolder, fileflag);
	}

	topoisomerase = 10.0;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_3_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, TF, stepwiseflag, runid, outputfolder, fileflag);
	}

	topoisomerase = 100.0;

	for(int i = 0; i < 16; i++)
	{
		outputfolder = "outputfiles/RUN_4_" + std::to_string(argid);
		runid = 16*world_rank + i;
		Gillespie_Simulation(force, gene_length0, gene_length1, spacer, config, clamp0_flag, clamp1_flag, T, controlgene, promoter0, promoter1, topoisomerase, TF, stepwiseflag, runid, outputfolder, fileflag);
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

	generator = std::mt19937(std::time(NULL) + world_rank);

	double TFarg = std::stod(argv[1]);
	int argid = std::stod(argv[2]);

	Run_Simulation(TFarg, argid, world_rank);

	MPI_Finalize();

	return(0);
}

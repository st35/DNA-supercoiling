std::mt19937 generator;

void Gillespie_Simulation(double force, double gene_length0, double gene_length1, double spacer, int config, int clamp0_flag, int clamp1_flag, double T, int controlgene, double promoter0, double promoter1, double topoisomerase, int stepwiseflag, int argid, std::string outputfolder, int fileflag)
{
	//if(gene_length0 > 0.0 && !(promoter0 > 0.0))
	//{
		//std::cout << "This is utter malarkey." << "\n";
		//return;
	//}
	//if(gene_length1 > 0.0 && !(promoter1 > 0.0))
	//{
		//std::cout << "Malarkey 2: Electric Boogaloo." << "\n";
		//return;
	//}

	double k_on_0 = (0.5 / 60.0)*promoter0;
	double k_on_1 = (0.5 / 60.0)*promoter1;
	double k_off = 0.0;
	double k_topo = (0.5 / 60.0)*topoisomerase;

	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	std::vector<std::vector<double>> phi, x, in_times, out_times, in_pos, out_pos;
	std::vector<double> Gene_Lengths, TSSs, Distance_Traveled, phi_x, dphi_xdt, Segments, Sigma, Torques;
	std::vector<int> Gene_Direction, PolII_Direction, PolII_Count, Finish_Count;

	Gene_Lengths.push_back(0.34*gene_length0);
	Gene_Lengths.push_back(0.34*gene_length1);
	Finish_Count.push_back(0);
	Finish_Count.push_back(0);

	double barrier = 0.34*10000.0;
	spacer = 0.34*spacer;

	clamp0 = 0.0;
	clamp1 = clamp0 + barrier + Gene_Lengths[1] + spacer + Gene_Lengths[0] + barrier;

	if(config == 0) // Tandem genes
	{
		TSSs.push_back(clamp0 + barrier + Gene_Lengths[1] + spacer);
		Gene_Direction.push_back(1);
		TSSs.push_back(clamp0 + barrier);
		Gene_Direction.push_back(1);
	}
	else if(config == 1) // Divergent genes
	{
		TSSs.push_back(clamp0 + barrier + Gene_Lengths[1] + spacer);
		Gene_Direction.push_back(1);
		TSSs.push_back(clamp0 + barrier + Gene_Lengths[1]);
		Gene_Direction.push_back(-1);
	}
	else if(config == 2) // Convergent genes
	{
		TSSs.push_back(clamp0 + barrier + Gene_Lengths[1] + spacer + Gene_Lengths[0]);
		Gene_Direction.push_back(-1);
		TSSs.push_back(clamp0 + barrier);
		Gene_Direction.push_back(1);
	}
	else
	{
		std::cout << "We killed Yamamoto." << "\n";
	}

	for(int i = 0; i < Gene_Lengths.size(); i++)
	{
		x.push_back(std::vector<double>());
		phi.push_back(std::vector<double>());
		in_times.push_back(std::vector<double>());
		out_times.push_back(std::vector<double>());
		in_pos.push_back(std::vector<double>());
		out_pos.push_back(std::vector<double>());
		PolII_Count.push_back(0);
	}

	int Total_PolII_Count = 0, finish_count = 0, numevents = 4, event, check_flag = 0, orig_count = 0, zero_flag = 0, spacer_flag = 0;
	std::vector<double> R;
	for(int i = 0; i < numevents; i++)
	{
		R.push_back(0.0);
	}

	double t = 0.0, dt, a0, p0, p1, sum_new, sum_prev, ini_phi, sub_t, sub_dt, epsilon = 1.0e-12, last_finish_0 = 0.0, last_finish_1 = 0.0;
	int counts_sum = 0, index, switchflag = 0;

	std::ofstream segment_file, sigma_file, torque_file, phi_file, x_file, v_file, av_file, spacer_sigma_file, RNA_file_0, RNA_file_1, density_file_0, density_file_1;
	if(fileflag == 1)
	{
		segment_file.open(outputfolder + "/segments_" + std::to_string(argid) + ".log");
		sigma_file.open(outputfolder + "/sigma_" + std::to_string(argid) + ".log");
		torque_file.open(outputfolder + "/torque_" + std::to_string(argid) + ".log");
		phi_file.open(outputfolder + "/phi_" + std::to_string(argid) + ".log");
		x_file.open(outputfolder + "/x_" + std::to_string(argid) + ".log");
		v_file.open(outputfolder + "/v_" + std::to_string(argid) + ".log");
		av_file.open(outputfolder + "/angular_v_" + std::to_string(argid) + ".log");
	}
	spacer_sigma_file.open(outputfolder + "/spacer_sigma_" + std::to_string(argid) + ".log");
	RNA_file_0.open(outputfolder + "/RNA_file_0_" + std::to_string(argid) + ".log");
	RNA_file_1.open(outputfolder + "/RNA_file_1_" + std::to_string(argid) + ".log");
	density_file_0.open(outputfolder + "/density_file_0_" + std::to_string(argid) + ".log");
	density_file_1.open(outputfolder + "/density_file_1_" + std::to_string(argid) + ".log");

	while(out_times[controlgene].size() < T)
	{
		counts_sum = 0;
		for(int i = 0; i < Gene_Lengths.size(); i++)
		{
			counts_sum += PolII_Count[i];
			if(x[i].size() != phi[i].size())
			{
				std::cout << "...home in time for Ellen." << "\n";
			}
		}
		if(Total_PolII_Count != counts_sum)
		{
			std::cout << "Google those rat tumors!" << "\n";
		}

		R[0] = k_on_0;
		R[1] = k_on_1;
		R[2] = k_off*Total_PolII_Count;
		R[3] = k_topo;
		a0 = 0.0;
		for(int i = 0; i < numevents; i++)
		{
			a0 += R[i];
		}
		p0 = distribution(generator);
		dt = (1.0 / a0)*std::log(1.0 / p0);

		t += dt;

		if(Total_PolII_Count > 0)
		{
			spacer_flag = 1;
			sub_t = t - dt;
			if(stepwiseflag == 1)
			{
				sub_dt = 1.0;
			}
			else
			{
				sub_dt = dt;
			}
			while(sub_t < t && Total_PolII_Count > 0)
			{
				if(sub_t + sub_dt > t && stepwiseflag == 1)
				{
					sub_dt = t - sub_t;
				}
				sub_t += sub_dt;
				phi_x.clear();
				PolII_Direction.clear();
				Distance_Traveled.clear();
				Total_PolII_Count = 0;
				for(int i = 0; i < Gene_Lengths.size(); i++)
				{
					if(Gene_Direction[i] == 1)
					{
						for(int j = 0; j < PolII_Count[i]; j++)
						{
							PolII_Direction.push_back(Gene_Direction[i]);
							Distance_Traveled.push_back(TSSs[i]);
							phi_x.push_back(phi[i][j]);
							Total_PolII_Count += 1;
						}
					}
					else if(Gene_Direction[i] == -1)
					{
						for(int j = PolII_Count[i] - 1; j >= 0; j--)
						{
							PolII_Direction.push_back(Gene_Direction[i]);
							Distance_Traveled.push_back(TSSs[i]);
							phi_x.push_back(phi[i][j]);
							Total_PolII_Count += 1;
						}
					}
				}
				for(int i = 0; i < Gene_Lengths.size(); i++)
				{
					if(Gene_Direction[i] == 1)
					{
						for(int j = 0; j < PolII_Count[i]; j++)
						{
							phi_x.push_back(x[i][j]);
						}
					}
					else if(Gene_Direction[i] == -1)
					{
						for(int j = PolII_Count[i] - 1; j >= 0; j--)
						{
							phi_x.push_back(x[i][j]);
						}
					}
				}

				counts_sum = 0;
				for(int i = 0; i < Gene_Lengths.size(); i++)
				{
					counts_sum += PolII_Count[i];
					if(x[i].size() != phi[i].size())
					{
						std::cout << "...home in time for Ellen." << "\n";
					}
				}
				if(Total_PolII_Count != counts_sum)
				{
					std::cout << "Google those rat tumors!" << "\n";
				}

				for(int i = 0; i < Total_PolII_Count - 1; i++)
				{
					if(phi_x[i + Total_PolII_Count] < phi_x[i + Total_PolII_Count + 1])
					{
						std::cout << "Hans Bethe." << "\n";
					}
				}

				if(phi_x.size() != Total_PolII_Count*2)
				{
					std::cout << "Treat yo'self!" << "\n";
				}

				if(check_flag == 1)
				{
					Segments = Get_Segment_Lengths(phi_x);
        				Sigma = Get_Supercoiling_Density(clamp0_flag, clamp1_flag, phi_x, Segments);
					zero_flag = 0;
					for(int i = 0; i < Sigma.size() - 1; i++)
					{
						if(std::abs(Sigma[i] - Sigma[i + 1]) < epsilon)
						{
							zero_flag = 1;
						}
					}
					if(zero_flag == 0)
					{
						std::cout << "Our knight is a witch." << "\n";
					}
					check_flag = 0;
				}

				boost::numeric::odeint::integrate_const(boost::numeric::odeint::runge_kutta4<std::vector<double>>(), odesystem{clamp0_flag, clamp1_flag, (int) Gene_Lengths.size(), PolII_Direction, Distance_Traveled, force}, phi_x, 0.0, sub_dt, 1e-3);

				dphi_xdt.clear();
				for(int i = 0; i < phi_x.size(); i++)
				{
					dphi_xdt.push_back(0.0);
				}
				odesystem{clamp0_flag, clamp1_flag, (int) Gene_Lengths.size(), PolII_Direction, Distance_Traveled, force}(phi_x, dphi_xdt, 0.0);
				for(int i = 0; i < PolII_Direction.size(); i++)
				{
					if(PolII_Direction[i]*dphi_xdt[i + Total_PolII_Count] < 0.0)
					{
						std::cout << "Saturday night is laundry night." << "\n";
					}
				}

				Segments = Get_Segment_Lengths(phi_x);
        			Sigma = Get_Supercoiling_Density(clamp0_flag, clamp1_flag, phi_x, Segments);
        			Torques = Get_Torques(Sigma, Segments, force);
				if(spacer_flag == 1)
				{
					spacer_sigma_file << Sigma[PolII_Count[0]] << " " << Torques[PolII_Count[0]] << "\n";
					spacer_flag = 0;
				}

				for(int i = 0; i < Sigma.size(); i++)
				{
					if(Torques[i]*Sigma[i] < 0.0)
					{
						std::cout << "National Constitution Center." << "\n";
					}
					if(Segments[i] < 0.0)
					{
						std::cout << "Supreme Disorder. By Ilya Shapiro." << "\n";
					}
				}

				if((Segments.size() != Sigma.size()) || (Sigma.size() != Total_PolII_Count + 1))
				{
					std::cout << "It's the way of the world." << "\n";
				}

				if(fileflag == 1)
				{
        				segment_file << sub_t << " ";
        				sigma_file << sub_t << " ";
					torque_file << sub_t << " ";
        				for(int i = 0; i < Segments.size(); i++)
        				{
            					segment_file << Segments[i] << " ";
            					sigma_file << Sigma[i] << " ";
						torque_file << Torques[i] << " ";
        				}
        				segment_file << "\n";
        				sigma_file << "\n";
					torque_file << "\n";

					phi_file << sub_t << " ";
        				x_file << sub_t << " ";
					v_file << sub_t << " ";
					av_file << sub_t << " ";
        				for(int i = 0; i < Total_PolII_Count; i++)
        				{
						phi_file << phi_x[i] << " ";
            					x_file << phi_x[i + Total_PolII_Count] << " ";
						v_file << dphi_xdt[i + Total_PolII_Count] << " ";
						av_file << dphi_xdt[i] << " ";
        				}
					phi_file << "\n";
        				x_file << "\n";
					v_file << "\n";
					av_file << "\n";
				}

				index = 0;
				finish_count = 0;
				for(int i = 0; i < Gene_Lengths.size(); i++)
				{
					phi[i].clear();
					x[i].clear();
					Finish_Count[i] = 0;
					for(int j = 0; j < PolII_Count[i]; j++)
					{
						if(std::abs(phi_x[index + Total_PolII_Count] - TSSs[i]) > Gene_Lengths[i])
						{
							out_pos[i].push_back(phi_x[index + Total_PolII_Count]);
							index += 1;
							finish_count += 1;
							Finish_Count[i] += 1;
							out_times[i].push_back(sub_t);
							continue;
						}
						phi[i].push_back(phi_x[index]);
						x[i].push_back(phi_x[index + Total_PolII_Count]);
						index += 1;
					}
					if(Gene_Direction[i] == -1)
					{
						std::reverse(phi[i].begin(), phi[i].end());
						std::reverse(x[i].begin(), x[i].end());
					}
					PolII_Count[i] = phi[i].size();
				}
				Total_PolII_Count -= finish_count;
				if(Finish_Count[0] > 0)
				{
					RNA_file_0 << ((double) Finish_Count[0]) / (sub_t - last_finish_0) << "\n";
					last_finish_0 = sub_t;
				}
				if(Finish_Count[1] > 0)
				{
					RNA_file_1 << ((double) Finish_Count[1]) / (sub_t - last_finish_1) << "\n";
					last_finish_1 = sub_t;
				}

				for(int i = 0; i < Gene_Lengths.size(); i++)
				{
					if(Gene_Direction[i] == 1)
					{
						for(int j = 0; j < PolII_Count[i] - 1; j++)
						{
							if(x[i][j] < x[i][j + 1])
							{
								std::cout << "Someone will die... of fun!" << "\n";
							}
						}
					}
					else if(Gene_Direction[i] == -1)
					{
						for(int j = 0; j < PolII_Count[i] - 1; j++)
						{
							if(x[i][j] > x[i][j + 1])
							{
								std::cout << "Someone will die... of fun!" << "\n";
							}
						}
					}
				}

				for(int i = 0; i < Gene_Lengths.size(); i++)
				{
					for(int j = 0; j < PolII_Count[i]; j++)
					{
						if(std::isnan(phi[i][j]))
						{
							std::cout << "Her fingernail is blue." << "\n";
						}
					}
				}

				if(stepwiseflag == 0)
				{
					break;
				}
			}
		}

		p1 = distribution(generator);
		sum_prev = 0.0;
		sum_new = 0.0;
		event = -1;
		for(int i = 0; i < numevents; i++)
		{
			sum_new = sum_prev + (R[i] / a0);
			if(p1 >= sum_prev && p1 < sum_new)
			{
				event = i;
				break;
			}
			sum_prev = sum_new;
		}

		if(event != 0 && event != 1 && event != 3)
		{
			std::cout << "Allergic to chestnuts... and good haircuts." << "\n";
		}

		orig_count = Total_PolII_Count;
		if(config == 0 && event < Gene_Lengths.size()) // Tandem genes
		{
			for(int i = 0; i < Gene_Lengths.size(); i++)
			{
				if(switchflag == 1 && event == i)
				{
					break;
				}
				if(event != i)
				{
					continue;
				}
				if(i == 0)
				{
					if(PolII_Count[0] == 0 && PolII_Count[1] == 0)
					{
						phi[i].push_back(0.0);
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] == 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							if(clamp0_flag == 0)
							{
								phi[i].push_back(phi[i][PolII_Count[i] - 1]);
							}
							else if(clamp0_flag == 1)
							{
								phi[i].push_back((TSSs[i] / x[i][PolII_Count[i] - 1])*phi[i][PolII_Count[i] - 1]);
							}
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
					else if(PolII_Count[0] == 0 && PolII_Count[1] > 0)
					{
						if(clamp1_flag == 0)
						{
							phi[i].push_back(phi[1][0]);
						}
						else if(clamp1_flag == 1)
						{
							phi[i].push_back(((TSSs[i] - clamp1) / (x[1][0] - clamp1))*phi[1][0]);
						}
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] > 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							phi[i].push_back(((TSSs[i] - x[i][PolII_Count[i] - 1]) / (x[1][0] - x[i][PolII_Count[i] - 1]))*(phi[1][0] - phi[i][PolII_Count[i] - 1]) + phi[i][PolII_Count[i] - 1]);
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
				}
				if(i == 1)
				{
					if(PolII_Count[0] == 0 && PolII_Count[1] == 0)
					{
						phi[i].push_back(0.0);
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] == 0)
					{
						if(clamp0_flag == 0)
						{
							phi[i].push_back(phi[0][PolII_Count[0] - 1]);
						}
						else if(clamp0_flag == 1)
						{
							phi[i].push_back((TSSs[i] / x[0][PolII_Count[0] - 1])*phi[0][PolII_Count[0] - 1]);
						}
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] == 0 && PolII_Count[1] > 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							if(clamp0_flag == 0)
							{
								phi[i].push_back(phi[i][PolII_Count[i] - 1]);
							}
							else if(clamp0_flag == 1)
							{
								phi[i].push_back((TSSs[i] / x[i][PolII_Count[i] - 1])*phi[i][PolII_Count[i] - 1]);
							}
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] > 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							if(clamp0_flag == 0)
							{
								phi[i].push_back(phi[i][PolII_Count[i] - 1]);
							}
							else if(clamp0_flag == 1)
							{
								phi[i].push_back((TSSs[i] / x[i][PolII_Count[i] - 1])*phi[i][PolII_Count[i] - 1]);
							}
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
				}
				for(int j = 0; j < PolII_Count[i]; j++)
				{
					if(std::isnan(phi[i][j]))
					{
						std::cout << "Lorenzo Wibberly." << "\n";
					}
				}
			}
		}
		else if(config == 1 && event < Gene_Lengths.size()) // Divergent genes
		{
			for(int i = 0; i < Gene_Lengths.size(); i++)
			{
				if(switchflag == 1 && event == i)
				{
					break;
				}
				if(event != i)
				{
					continue;
				}
				if(i == 0)
				{
					if(PolII_Count[0] == 0 && PolII_Count[1] == 0)
					{
						phi[i].push_back(0.0);
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] == 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							if(clamp0_flag == 0)
							{
								phi[i].push_back(phi[i][PolII_Count[i] - 1]);
							}
							else if(clamp0_flag == 1)
							{
								phi[i].push_back((TSSs[i] / x[i][PolII_Count[i] - 1])*phi[i][PolII_Count[i] - 1]);
							}
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
					else if(PolII_Count[0] == 0 && PolII_Count[1] > 0)
					{
						if(clamp1_flag == 0)
						{
							phi[i].push_back(phi[1][PolII_Count[1] - 1]);
						}
						else if(clamp1_flag == 1)
						{
							phi[i].push_back(((TSSs[i] - clamp1) / (x[1][PolII_Count[1] - 1] - clamp1))*phi[1][PolII_Count[1] - 1]);
						}
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] > 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							phi[i].push_back(((TSSs[i] - x[i][PolII_Count[i] - 1]) / (x[1][PolII_Count[1] - 1] - x[i][PolII_Count[i] - 1]))*(phi[1][PolII_Count[1] - 1] - phi[i][PolII_Count[i] - 1]) + phi[i][PolII_Count[i] - 1]);
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
				}
				if(i == 1)
				{
					if(PolII_Count[0] == 0 && PolII_Count[1] == 0)
					{
						phi[i].push_back(0.0);
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] == 0)
					{
						if(clamp0_flag == 0)
						{
							phi[i].push_back(phi[0][PolII_Count[0] - 1]);
						}
						else if(clamp0_flag == 1)
						{
							phi[i].push_back((TSSs[i] / x[0][PolII_Count[0] - 1])*phi[0][PolII_Count[0] - 1]);
						}
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] == 0 && PolII_Count[1] > 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							if(clamp1_flag == 0)
							{
								phi[i].push_back(phi[i][PolII_Count[i] - 1]);
							}
							else if(clamp1_flag == 1)
							{
								phi[i].push_back(((clamp1 - TSSs[i]) / (clamp1 - x[i][PolII_Count[i] - 1]))*phi[i][PolII_Count[i] - 1]);
							}
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] > 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							phi[i].push_back(((TSSs[i] - x[0][PolII_Count[0] - 1]) / (x[1][PolII_Count[1] - 1] - x[0][PolII_Count[0] - 1]))*(phi[1][PolII_Count[1] - 1] - phi[0][PolII_Count[0] - 1]) + phi[0][PolII_Count[0] - 1]);
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
				}
				for(int j = 0; j < PolII_Count[i]; j++)
				{
					if(std::isnan(phi[i][j]))
					{
						std::cout << "Lorenzo Wibberly." << "\n";
					}
				}
			}
		}
		else if(config == 2 && event < Gene_Lengths.size()) // Convergent genes
		{
			for(int i = 0; i < Gene_Lengths.size(); i++)
			{
				if(switchflag == 1 && event == i)
				{
					break;
				}
				if(event != i)
				{
					continue;
				}
				if(i == 0)
				{
					if(PolII_Count[0] == 0 && PolII_Count[1] == 0)
					{
						phi[i].push_back(0.0);
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] == 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							if(clamp1_flag == 0)
							{
								phi[i].push_back(phi[i][PolII_Count[i] - 1]);
							}
							else if(clamp1_flag == 1)
							{
								phi[i].push_back(((clamp1 - TSSs[i]) / (clamp1 - x[i][PolII_Count[i] - 1]))*phi[i][PolII_Count[i] - 1]);
							}
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
					else if(PolII_Count[0] == 0 && PolII_Count[1] > 0)
					{
						if(clamp1_flag == 0)
						{
							phi[i].push_back(phi[1][0]);
						}
						else if(clamp1_flag == 1)
						{
							phi[i].push_back(((TSSs[i] - clamp1) / (x[1][0] - clamp1))*phi[1][0]);
						}
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] > 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							if(clamp1_flag == 0)
							{
								phi[i].push_back(phi[0][PolII_Count[0] - 1]);
							}
							else if(clamp1_flag == 1)
							{
								phi[i].push_back(((TSSs[i] - clamp1) / (x[0][PolII_Count[0] - 1] - clamp1))*phi[0][PolII_Count[0] - 1]);
							}
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
				}
				if(i == 1)
				{
					if(PolII_Count[0] == 0 && PolII_Count[1] == 0)
					{
						phi[i].push_back(0.0);
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] == 0)
					{
						if(clamp0_flag == 0)
						{
							phi[i].push_back(phi[0][0]);
						}
						else if(clamp0_flag == 1)
						{
							phi[i].push_back((TSSs[i] / x[0][0])*phi[0][0]);
						}
						x[i].push_back(TSSs[i]);
						PolII_Count[i] += 1;
						in_times[i].push_back(t);
						in_pos[i].push_back(TSSs[i]);
						Total_PolII_Count += 1;
					}
					else if(PolII_Count[0] == 0 && PolII_Count[1] > 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							if(clamp0_flag == 0)
							{
								phi[i].push_back(phi[i][PolII_Count[i] - 1]);
							}
							else if(clamp0_flag == 1)
							{
								phi[i].push_back((TSSs[i] / x[i][PolII_Count[i] - 1])*phi[i][PolII_Count[i] - 1]);
							}
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
					else if(PolII_Count[0] > 0 && PolII_Count[1] > 0)
					{
						if(std::abs(x[i][PolII_Count[i] - 1] - TSSs[i]) > delta)
						{
							if(clamp0_flag == 0)
							{
								phi[i].push_back(phi[i][PolII_Count[i] - 1]);
							}
							else if(clamp0_flag == 1)
							{
								phi[i].push_back((TSSs[i] / x[i][PolII_Count[i] - 1])*phi[i][PolII_Count[i] - 1]);
							}
							x[i].push_back(TSSs[i]);
							PolII_Count[i] += 1;
							in_times[i].push_back(t);
							in_pos[i].push_back(TSSs[i]);
							Total_PolII_Count += 1;
						}
					}
				}
				for(int j = 0; j < PolII_Count[i]; j++)
				{
					if(std::isnan(phi[i][j]))
					{
						std::cout << "Lorenzo Wibberly." << "\n";
					}
				}
			}
		}
		else if(event < Gene_Lengths.size())
		{
			std::cout << "It's never lupus." << "\n";
		}

		if(event == Gene_Lengths.size() + 1)
		{
			for(int i = 0; i < Gene_Lengths.size(); i++)
			{
				for(int j = 0; j < PolII_Count[i]; j++)
				{
					phi[i][j] = 0.0;
				}
			}
		}
		density_file_0 << PolII_Count[0] << "\n";
		density_file_1 << PolII_Count[1] << "\n";

		check_flag = 0;
		if(Total_PolII_Count > orig_count)
		{
			check_flag = 1;
		}
		if(Total_PolII_Count - orig_count > 1)
		{
			std::cout << "May be they just bought their parasites at the mall." << "\n";
		}
//		if(in_times[0].size() > 0)
//		{
//			switchflag = 1;
//		}
		if(switchflag == 1 && out_times[0].size() > 0)
		{
			break;
		}
	}

	for(int i = 0; i < 2; i++)
	{
		if(in_times[i].size() != in_pos[i].size())
		{
			std::cout << "Dr. Boots Hofstadter." << "\n";
		}
		if(out_times[i].size() != out_pos[i].size())
		{
			std::cout << "Dr. Boots Hofstadter." << "\n";
		}
	}

	std::ofstream transcription_rate_0, transcription_rate_1;
	transcription_rate_0.open(outputfolder + "/transcription_rates_" + std::to_string(argid) + "_0.log");
	transcription_rate_1.open(outputfolder + "/transcription_rates_" + std::to_string(argid) + "_1.log");

	for(int i = 0; i < out_times[0].size() ; i++)
	{
		transcription_rate_0 << in_times[0][i] << " " << (std::abs(out_pos[0][i] - in_pos[0][i]) / 0.34) / (out_times[0][i] - in_times[0][i]) << "\n";
	}

	for(int i = 0; i < out_times[1].size() ; i++)
	{
		transcription_rate_1 << in_times[1][i] << " " << (std::abs(out_pos[1][i] - in_pos[1][i]) / 0.34) / (out_times[1][i] - in_times[1][i]) << "\n";
	}

	if(fileflag == 1)
	{
		segment_file.close();
		sigma_file.close();
		torque_file.close();
		phi_file.close();
		x_file.close();
		v_file.close();
		av_file.close();
	}
	spacer_sigma_file.close();
	transcription_rate_0.close();
	transcription_rate_1.close();
	RNA_file_0.close();
	RNA_file_1.close();
	density_file_0.close();
	density_file_1.close();

	return;
}

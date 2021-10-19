double Get_Torque(double sigma, double f, double segment)
{
	double cutoff = 0.34*1000.0;

	double A = 50.0, C = 95.0, P = 24.0;
	double kBT = 4.1; // T = 300K.
	double A_m = 4.0, C_m = 1.75, e_m = 6.0*kBT, sigma_0 = -1.0;

	double c = kBT*C*w0*w0;
	double p = kBT*P*w0*w0;

	double cs = c*(1.0 - ((C / (4.0*A))*std::sqrt(kBT / (A*f))));
	double g = f - std::sqrt((kBT*f / A));

	double factor = std::sqrt((2.0*p*g) / (1.0 - (p / cs)));

	double sigma_s = factor / cs;
	sigma_s = sigma_s*(1.0 + std::pow((cutoff / segment), 2.0));
	double sigma_p = factor / p;

	double g_m = 1.2*(f - std::sqrt(kBT*f / A_m));
	double c_m = kBT*C_m*w0*w0;

	double sigma_sm = (c_m / (cs - c_m))*(-sigma_0 - std::sqrt(sigma_0*sigma_0 + (2.0*(cs - c_m) / (cs*c_m))*(g + e_m - g_m)));
	double sigma_m = sigma_0 + (cs / (cs - c_m))*(-sigma_0 - std::sqrt(sigma_0*sigma_0 + (2.0*(cs - c_m) / (cs*c_m))*(g + e_m - g_m)));

	double torque = 0.0;

	if(sigma < -2.0)
	{
		if(sigma_error_flag == 0)
		{
			std::cout << "Going to cscQ." << "\n";
			sigma_error_flag = 1;
		}
	}
	if(sigma > 2.6)
	{
		if(sigma_error_flag == 0)
		{
			std::cout << "Going to Honululu." << "\n";
			sigma_error_flag = 1;
		}
	}

	if(sigma <= sigma_m)
	{
		torque = (c_m / w0)*(sigma - sigma_0);
	}
	else if(sigma > sigma_m && sigma <= sigma_sm)
	{
		torque = (c_m / w0)*(sigma_m - sigma_0);
	}
	else if(sigma > sigma_sm && sigma <= sigma_s)
	{
		torque = (cs / w0)*sigma;
	}
	else if(sigma > sigma_s && sigma <= sigma_p)
	{
		torque = factor / w0;
	}
	else if(sigma > sigma_p)
	{
		torque = (p / w0)*sigma;
		if(torque > 40.0)
		{
			torque = 40.0;
		}
	}

	return(torque);
}

std::vector<double> Get_Segment_Lengths(std::vector<double> phi_x)
{
	std::vector<double> Segments;
	int PolII_Count = (phi_x.size() - 2) / 2;
	double seg_left, seg_right;
	for(int i = 0; i < PolII_Count; i++)
	{
		if(PolII_Count == 1)
		{
			seg_left = phi_x[i + PolII_Count] - clamp0;
			seg_right = clamp1 - phi_x[i + PolII_Count];
		}
		else
		{
			if(i == 0)
			{
				seg_left = phi_x[i + PolII_Count] - phi_x[i + 1 + PolII_Count];
				seg_right = clamp1 - phi_x[i + PolII_Count];
			}
			else if(i == PolII_Count - 1)
			{
				seg_left = phi_x[i + PolII_Count] - clamp0;
				seg_right = phi_x[i - 1 + PolII_Count] - phi_x[i + PolII_Count];
			}
			else
			{
				seg_left = phi_x[i + PolII_Count] - phi_x[i + 1 + PolII_Count];
				seg_right = phi_x[i - 1 + PolII_Count] - phi_x[i + PolII_Count];
			}
		}
		Segments.push_back(seg_right);
		if(i == PolII_Count - 1)
		{
			Segments.push_back(seg_left);
		}
	}

	return(Segments);
}

std::vector<double> Get_Supercoiling_Density(int clamp0_flag, int clamp1_flag, std::vector<double> phi_x, std::vector<double> Segments)
{
	std::vector<double> Sigma;
	int PolII_Count = (phi_x.size() - 2) / 2;
	double sigma_left, sigma_right, seg_left, seg_right;
	for(int i = 0; i < PolII_Count; i++)
	{
		seg_left = Segments[i + 1];
		seg_right = Segments[i];
		if(PolII_Count == 1)
		{
			if(clamp0_flag == 0)
			{
				sigma_left = 0.0;
			}
			else
			{
				sigma_left = (phi_x[2*PolII_Count] - phi_x[i]) / (w0*seg_left);
			}
			if(clamp1_flag == 0)
			{
				sigma_right = 0.0;
			}
			else
			{
				sigma_right = (phi_x[i] - phi_x[2*PolII_Count + 1]) / (w0*seg_right);
			}
		}
		else
		{
			if(i == 0)
			{
				sigma_left = (phi_x[i + 1] - phi_x[i]) / (w0*seg_left);
				if(clamp1_flag == 0)
				{
					sigma_right = 0.0;
				}
				else
				{
					sigma_right = (phi_x[i] - phi_x[2*PolII_Count + 1]) / (w0*seg_right);
				}
			}
			else if(i == PolII_Count - 1)
			{
				if(clamp0_flag == 0)
				{
					sigma_left = 0.0;
				}
				else
				{
					sigma_left = (phi_x[2*PolII_Count] - phi_x[i]) / (w0*seg_left);
				}
				sigma_right = (phi_x[i] - phi_x[i - 1]) / (w0*seg_right);
			}
			else
			{
				sigma_left = (phi_x[i + 1] - phi_x[i]) / (w0*seg_left);
				sigma_right = (phi_x[i] - phi_x[i - 1]) / (w0*seg_right);
			}
		}
//		Sigma.push_back(sigma_right / (2.0*3.14));
		Sigma.push_back(sigma_right);
		if(i == PolII_Count - 1)
		{
//			Sigma.push_back(sigma_left / (2.0*3.14));
			Sigma.push_back(sigma_left);
		}
	}

	return(Sigma);
}

std::vector<double> Get_Torques(std::vector<double> Sigma, std::vector<double> Segments, double force)
{
	std::vector<double> Torques;
	for(int i = 0; i < Sigma.size(); i++)
	{
		Torques.push_back(Get_Torque(Sigma[i], force, Segments[i]));
	}

	return(Torques);
}

std::vector<double> Get_Velocities(std::vector<int> Direction, std::vector<double> Torques, std::vector<double> Segments)
{
	std::vector<double> Velocities;
	double tau_left, tau_right, seg_left, seg_right;
	if(Torques.size() == 0)
	{
		return(Velocities);
	}
	for(int i = 0; i < Torques.size() - 1; i++)
	{
		tau_left = Torques[i + 1];
		tau_right = Torques[i];
		seg_left = Segments[i + 1];
		seg_right = Segments[i];
		if(Direction[i] == 1)
		{
			if(seg_right > delta)
			{
				Velocities.push_back((v0 / 2.0)*(1.0 - std::tanh((tau_right - tau_left) / tau_c)));
			}
			else
			{
				Velocities.push_back(0.0);
			}
		}
		else if(Direction[i] == -1)
		{
			if(seg_left > delta)
			{
				Velocities.push_back(-(v0 / 2.0)*(1.0 - std::tanh(-(tau_right - tau_left) / tau_c)));
			}
			else
			{
				Velocities.push_back(0.0);
			}
		}
	}

	return(Velocities);
}

struct odesystem
{
	int clamp0_flag, clamp1_flag, num_genes;
	double force, motor0, motor1;
	std::vector<int> Counts, Direction;
	std::vector<double> TSS;
	odesystem(int param1, int param2, int param3, std::vector<int> param4, std::vector<double> param5, double param6, double param7, double param8)
	{
		clamp0_flag = param1;
		clamp1_flag = param2;
		num_genes = param3;
		for(int i = 0; i < param4.size(); i++)
		{
			Direction.push_back(param4[i]);
		}
		for(int i = 0; i < param5.size(); i++)
		{
			TSS.push_back(param5[i]);
		}
		force = param6;
		motor0 = param7;
		motor1 = param8;
	}

	void operator()(const std::vector<double> &phi_x, std::vector<double> &dphi_xdt, double t)
	{
		std::vector<double> Segments = Get_Segment_Lengths(phi_x);
		std::vector<double> Sigma = Get_Supercoiling_Density(clamp0_flag, clamp1_flag, phi_x, Segments);
		std::vector<double> Torques = Get_Torques(Sigma, Segments, force);
		std::vector<double> Velocities = Get_Velocities(Direction, Torques, Segments);

		int PolII_Count = (phi_x.size() - 2) / 2;
		double denom;
		for(int i = 0; i < PolII_Count; i++)
		{
			denom = chi + eta*std::pow(std::abs(phi_x[i + PolII_Count] - TSS[i]), alpha);
			dphi_xdt[i] = 0.0;
			dphi_xdt[i + PolII_Count] = 0.0;
			if(std::abs(Velocities[i]) > 0.0)
			{
				dphi_xdt[i] = ((w0*Velocities[i])*(eta*std::pow(std::abs(phi_x[i + PolII_Count] - TSS[i]), alpha)) / denom) - ((Torques[i] - Torques[i + 1]) / denom);
			}
			dphi_xdt[i + PolII_Count] = Velocities[i];
		}
		dphi_xdt[2*PolII_Count] = motor0;
		dphi_xdt[2*PolII_Count + 1] = motor1;
	}
};

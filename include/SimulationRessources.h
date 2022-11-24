#ifndef SIMULATIONRESSOURCES_H
#define SIMULATIONRESSOURCES_H
#include<qstringlist.h>
#include<math.h>
#include<vector>
#include<array>
#pragma once

struct PennesEquationParameter
{




	double rho = 1068.0;	//tissue density[kg / m3]		///	1000.0		  ->darf nicht null sein!! 1068.0
	double c = 3600.0; //tissue specific heat capacity [J/kg °C]   //4200.0      ->darf nicht null sein!! 3600.0
	double T =21.0;	//tissue temperature [◦C]
	
	double k_a = 0.06;
	double k_b = 0.5;
	double k = 2.0; 	// tissue thermal conductivity[W / (m °C)]         -> darf nicht null sein!! //0.59;
	
	double getK()
	{
		return k_a * T + k_b;
	}
	double M = rho * c;//volumetric heat capacity {J/⋅K⋅m3]


	//double td = k / M;//thermal diffusivity (m^2/s)
	double td = 1.5;//thermal diffusivity (mm^2/s) 1.4
	double td_a = 0.155;
	double td_b = 4.95;
	double td_c = 2.01;

	double getTd()
	{
		return td * 1e-6;
		//return (td_a + T * td_b) * 1e-6;
		/*double td = (td_a + td_b * 1e-10 * exp(td_c * 1e-1 * T)) * 1e-6;
		if (td > 0)
			return td;
		else
			return 0.1;*/
	}

	//Volumetric heat capacity

	double w_b = 0.0;	//blood perfusion rate [kg/(m3 s] 500.0
	double c_b = 0.0;  //blood specific capacity[J / kg  °C] 4200.0
	double T_a = 0.0;	//arterial temperature[°C] 37.0
	double Q_m = 0.0;	//metabolic heat generation rate[W / m3]  //33800.0
	double Q_r = 0.0;	//regional heat source[W / m3]
	double Q_r_rel = 0.0;
	double T_ind_max = 0.0; //Q_r/(rho * c)	max induzierte Temperatur

	//aus bu lin
	bool isNecrotic = false;
	bool keepTconstant = false;
	bool isVessel = false;
	//PennesEquationParameter(const PennesEquationParameter& pep) : rho(pep.rho), c(pep.c) {}
};
struct defaultValues
{

};
struct HasToBeOptimized
{
	bool rho = false;
	bool c = false;
	bool T = false;
	bool k = false;
	bool w_b = false;
	bool c_b = false;
	bool T_a = false;
	bool Q_m = false;
	bool Q_r = false;
};
enum ParameterFlags
{
	td = 0,
	k = 1,

	k_a = 2,
	k_b = 3,

	td_a = 4,
	td_b = 5,
	td_c = 6,

	T_ind_max = 7,

	c_b = 8,
};
//namespace {
	class GlobalVariables {
	public:
		//float GlobalVariables::needleDir[3] = { 0 };
		
		static const int numberOfSlices = 8;
		static const int numberOfOptiParam = 2;
		static const float dt;
		static const float timeBetweenMeas;

		static int slice_index;
		static std::vector<std::array<int,3>> needleInVolume;
		static float needleDir[3];
		static float pointOnNeedle_world[3];
		static int pointOnNeedle_volume[3];
		static int boundaryCondition; //1 == ditrlichet, 0==von neumann

		static std::string angleList[numberOfSlices];//"112_5"
		static QStringList slices_time;
		static float slices_time_f[GlobalVariables::numberOfSlices];
		static const int baseLineT = 17;
		static ParameterFlags paramToOpti [numberOfOptiParam];

	};

	
	
//}
#endif // !SIMULATIONRESSOURCES_H
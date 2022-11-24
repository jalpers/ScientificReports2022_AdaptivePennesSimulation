//#include "..\include\ThermalDose.h"
#include<ThermalDose.h>
ThermalDose::ThermalDose()
{
}
ThermalDose::ThermalDose(int _N)
{
	N = _N;
	CEM43 = new float[N];
	ThermalDamage = new float[N];
	for (int i = 0; i < N; i++)
		CEM43[i] = 0;
}
ThermalDose::~ThermalDose()
{
	delete CEM43;
	delete ThermalDamage;
}

void ThermalDose::update(float* _T_i, float _t_i)
{
	float R = 0;
	for (int i = 0; i < N; i++)
	{

		if (_T_i[i] < threshold43)
			R = R_lesser43;
		else
			R = R_greater43;

		
		//std::cout <<_T_i[i] <<" "<< powf(R, (threshold43 - _T_i[i])) << " "<<R <<std::endl;
		CEM43[i] += _t_i * powf(R, threshold43 - _T_i[i]);

		
	}
}
void ThermalDose::computeNecrosisZones()
{
	for (int i = 0; i < N; i++)
	{

		if (CEM43[i] < 250)
			ThermalDamage[i] = CEM43[i];
		else
			ThermalDamage[i] = 250;
		//if (CEM43[i] > necroticZone_thresh)
		//	ThermalDamage[i] = necroticZone_flag;
		//else if(CEM43[i] > transitionZone_thresh)
		//	ThermalDamage[i] = transitionZone_flag;
		//else 
		//	ThermalDamage[i] = healthyZone_flag;
	}
}

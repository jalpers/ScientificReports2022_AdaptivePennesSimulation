#ifndef THERMALDOSE_H
#define THERMALDOSE_H
#include <math.h>
#include<iostream>
class ThermalDose
{
public:
	ThermalDose();
	ThermalDose(int _N);
	~ThermalDose();
	void update(float* _T_i, float t_i);
	void computeNecrosisZones();
	float* getNecrosiszones() { return ThermalDamage; }
	float* getCEM43() { return CEM43; }
	
private:
	const float threshold43 = 43.0;
	const float R_greater43 = 0.5;
	const float R_lesser43 = 0.25;

	const int necroticZone_flag = 2;
	const float necroticZone_thresh = 240;
	const int transitionZone_flag = 1;
	const float transitionZone_thresh = 30;
	const int healthyZone_flag = 0;


	int N;
	float* CEM43;
	float* ThermalDamage;
};


#endif // !THERMALDOSE_H
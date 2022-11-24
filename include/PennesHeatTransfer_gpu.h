#ifndef PENNES_HEAT_TRANSFER_GPU_H
#define PENNES_HEAT_TRANSFER_GPU_H
#include <cuda_runtime.h>
#include <cuda.h>
#include<iostream>
#include<Matrix3D.h>
#include<SimulationRessources.h>
#include<omp.h>
#define BUFFER_HEATSOURCES 128
#define GPU 1
#define DEBUG 0
__device__ __host__
struct Parameters
{
	float td_a;
	float td_b;
	float td_c;
	float cb;
	float Ta;
	float k;

	float defaultT;
	int boundCond;
	float dx_2[3];
	float dt;
	dim3 dim;

	// buffer[32 - 14];
	int N_heat;
	float T_heat[BUFFER_HEATSOURCES];
	int heat_x[BUFFER_HEATSOURCES];
	int heat_y[BUFFER_HEATSOURCES];
	int heat_z[BUFFER_HEATSOURCES];

	Parameters(
		float _td_a,
		float _td_b,
		float _td_c,
		float _cb,
		float _Ta,
		float _k,
		float _defaultT,
		int _boundCond,
		float* _dx_2,
		float _dt,
		dim3 _dim,
		int _N_heat,
		float* _T_heat,
		float* _heat_x,
		float* _heat_y,
		float* _heat_z)
	{
		td_a = _td_a;
		td_b = _td_b;
		td_c = _td_c;
		cb = _cb;
		Ta = _Ta;
		k = _k;

		defaultT = _defaultT;
		boundCond = _boundCond;
		for (int i = 0; i < 3; i++)
		{
			dx_2[i] = _dx_2[i];
		}
		dim = _dim;
		N_heat = _N_heat;
		for (int i = 0; i < BUFFER_HEATSOURCES; i++)
		{
			T_heat[i] = _T_heat[i];
			heat_x[i] = _heat_x[i];
			heat_y[i] = _heat_y[i];
			heat_z[i] = _heat_z[i];
		}
		
	}
		
	Parameters() {
		td_a = 1.4;
		td_b = 0.1;
		td_c = 0.1;
		cb = 4182;
		Ta = 25;
		k = 2;
		defaultT = 21;
		boundCond = 1;
		for (int i = 0; i < 3; i++)
		{
			dx_2[i] = 4.0 * 1e-6;
		}
		dim = { 60,60,60 };
		N_heat = 60;
		for (int i = 0; i < BUFFER_HEATSOURCES; i++)
		{
			T_heat[i] = 100;
			heat_x[i] = 30;
			heat_y[i] = i;
			heat_z[i] =30;
		}
	}
		
	


};
__device__ __host__
struct test_str
{
	float u;
	float x;
};
enum Parameters_ind
{
	/*td_a = 0,
	td_b,
	cb,
	T_a,
	defaultT,
	boundCond*/
};
class PennesHeatTransfer_Gpu 
{
public:
	PennesHeatTransfer_Gpu();
	PennesHeatTransfer_Gpu(unsigned int *_dim);
	PennesHeatTransfer_Gpu(const PennesHeatTransfer_Gpu& pht);
	~PennesHeatTransfer_Gpu();

	void setT_h(Matrix3D<float>* _T_h);
	void setT_d(float *_Tsaved_d);
	void setAll2(float _T);
	void  updateTnew_h();
	float* getTsaved_d() { return Tsaved_d; }
	Matrix3D<float>* getTsaved_h() { return Tsaved_h; }
	Matrix3D<float>* getTnew_h() { updateTnew_h(); return Tnew_h; }
	
	void setDefaultParameters();
	void setPerfusion(Matrix3D<bool>& _map);
	void finiteStep(int _NrOfTimeSteps = 1 );
	void saveTimestep();
	
	Parameters* par_h;
	static float* Tsaved_d;
	static Matrix3D<float>* Tsaved_h;
private:
#if GPU
	float* T_d, * Ts_d, * Tss_d, * Tnew_d, *wb_d;
	Parameters* par_d;
#endif

	
	//dim3 dim;
	int n;
	int handle;
#if GPU	
	void finiteStep_GPU(int _NrOfTimeSteps);
#endif
	void finiteStep_CPU(int _NrOfTimeSteps);
	Matrix3D<float>* T_h;
	Matrix3D<float>* Tnew_h;
	
	Matrix3D<float>* wb_h;
	
	void updateParameters();
	void printParameters();
	
};

#endif // !PENNES_HEAT_TRANSFER_GPU_H

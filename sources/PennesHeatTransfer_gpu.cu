#include <PennesHeatTransfer_gpu.h>
#include<SimulationRessources.h>
#include<chrono>
#define NUMBER_OF_ARRAYS 4
#include<cuda_util.h>
Matrix3D<float>* PennesHeatTransfer_Gpu::Tsaved_h = nullptr;
#if GPU
float *PennesHeatTransfer_Gpu::Tsaved_d = nullptr;

__global__ void cuda_hello(Parameters *_par_d) {//float* _T_d, float* _Ts_d
	printf("Hello World from GPU!\n");
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	//printf("x %d, y %d, z %d \n", dim.x, dim.y, dim.z);
	//test_str* _par = (test_str *) _par_d;
	//printf("%f ", (*_par).u);

}
__device__ __host__
float getTd(Parameters* _p, float _T) 
{ 
	//(td_a + td_b * 1e-10 * exp(td_c * 1e-1 * T)) * 1e-6;
	if (_T >= 21 && _T < 200)
	{
		//printf("%f ", exp(_p->td_c * _T));
		return float((_p->td_a + _p->td_b * exp(_p->td_c * _T)) * 1e-6);//
		//return float((1 + _T * _p->td_b) * _p->td_a * 1e-6);
	}

	else
		return float(_p->td_a*1e-6);

}
__global__ void copyT(float* _src, float* _dest, int n)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= n)
		return;

	_dest[tid] = _src[tid];

}
__global__ void setHeatSources(float* _T_d, Parameters* _par_d)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	//printf("tid %d ",tid);
	Matrix3D<float> T_m = Matrix3D<float>(_par_d->dim.x, _par_d->dim.y, _par_d->dim.z, _T_d);

	if((int)_par_d->T_heat[tid] != (int) _par_d->defaultT)
		T_m.setValue({ _par_d->heat_x[tid],_par_d->heat_y[tid] ,_par_d->heat_z[tid] }, _par_d->T_heat[tid]);
}
__device__ __host__ void thomas(float* a, float* b, float* c, float* d, int n) {
	/*
	// n is the number of unknowns

	|b0 c0 0 ||x0| |d0|
	|a1 b1 c1||x1|=|d1|
	|0  a2 b2||x2| |d2|

	1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

		x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

	2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
		from 1st it.: -| a1x0 + a1g0x1        = a1r0
					-----------------------------
						  (b1 - a1g0)x1 + c1x2 = d1 - a1r0

		x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

	3rd iteration:      | a2x1 + b2x2   = d2
		from 2nd it. : -| a2x1 + a2g1x2 = a2r2
					   -----------------------
					   (b2 - a2g1)x2 = d2 - a2r2
		x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
	Finally we have a triangular matrix:
	|1  g0 0 ||x0| |r0|
	|0  1  g1||x1|=|r1|
	|0  0  1 ||x2| |r2|

	Condition: ||bi|| > ||ai|| + ||ci||

	in this version the c matrix reused instead of g
	and             the d matrix reused instead of r and x matrices to report results
	Written by Keivan Moradi, 2014
	*/
	n--; // since we start from x0 (not x1)
	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++) {
		c[i] /= b[i] - a[i] * c[i - 1];
		d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
	}

	d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

	for (int i = n; i-- > 0;) {

		d[i] -= c[i] * d[i + 1];
	}
}
__global__ void zSweep(float* _T_d, float* _Ts_d, float* _Tss_d, float * _Tnew_d, float* _wb_d, Parameters* _par_d )
{
	extern __shared__ float arrZ[];


	dim3 dim = { _par_d->dim.x, _par_d->dim.y, _par_d->dim.z };
	float defaultValue = _par_d->defaultT;
	float r1, r2, r3;
	float td;
	float dt = _par_d->dt;
	float dx_2[3] = {_par_d->dx_2[0], _par_d->dx_2[1], _par_d->dx_2[2] };


	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= (dim.x * dim.y))
		return;




	float* a = &arrZ[threadIdx.x * NUMBER_OF_ARRAYS * dim.z + 0];
	float* b = &arrZ[threadIdx.x * NUMBER_OF_ARRAYS * dim.z + dim.z];
	float* c = &arrZ[threadIdx.x * NUMBER_OF_ARRAYS * dim.z + dim.z * 2];
	float* d = &arrZ[threadIdx.x * NUMBER_OF_ARRAYS * dim.z + dim.z * 3];

	Matrix3D<float> T_m = Matrix3D<float>(dim.x, dim.y, dim.z, _T_d);
	Matrix3D<float> Ts_m = Matrix3D<float>(dim.x, dim.y, dim.z, _Ts_d);
	Matrix3D<float> Tss_m = Matrix3D<float>(dim.x, dim.y, dim.z, _Tss_d);
	Matrix3D<float> Tnew_m = Matrix3D<float>(dim.x, dim.y, dim.z, _Tnew_d);
	Matrix3D<float> wb_m = Matrix3D<float>(dim.x, dim.y, dim.z, _wb_d);


	

	T_m.setAmbientValue(defaultValue);
	Ts_m.setAmbientValue(defaultValue);
	Tss_m.setAmbientValue(defaultValue);


	int x = tid % dim.x;
	int y = (int)(tid / dim.x);
	
	for (int i = 0; i < dim.z; i++)
	{
		int z = i;
		/*if (x >= _dim.x || x < 0 ||
			y >= _dim.y || y < 0 ||
			z >= _dim.z || z < 0)
			printf("%d , %d , %d \n", x, y, z);*/
		float T_act = (T_m)[{x, y, z}];
		td = getTd(_par_d, T_act);
		
		r1 = r2 = r3 = td * dt / dx_2[2];
		if (r1 < 0)
			printf("r1 %f \n", r1);
		//if ((T_m)[{x, y, z}] > 100 || (T_m)[{x, y, z}] < 0)
		//	printf("Tm %f ", (T_m)[{x, y, z}]);
		//if ((Ts_m)[{x, y, z}] > 100 || (Ts_m)[{x, y, z}] < 0)
		//	printf("Tsm %f ", (Ts_m)[{x, y, z}]);
		//if ((Tss_m)[{x, y, z}] > 100 || (Tss_m)[{x, y, z}] < 0)
		//	printf("Tssm %f \n", (Tss_m)[{x, y, z}]);


		a[i] = -r3 / 2;
		b[i] = 1 + r3;
		c[i] = -r3 / 2;

		d[i] = r1 / 2 * (T_m)[{x - 1, y, z}] + r1 / 2 * (T_m)[{x + 1, y, z}]
			+ r1 / 2 * (Ts_m)[{x - 1, y, z}] + r1 / 2 * (Ts_m)[{x + 1, y, z}]
			+ r2 / 2 * (T_m)[{x, y - 1, z}] + r2 / 2 * (T_m)[{x, y + 1, z}]
			+ r2 / 2 * (Tss_m)[{x, y - 1, z}] + r2 / 2 * (Tss_m)[{x, y + 1, z}]
			+ r3 / 2 * (T_m)[{x, y, z - 1}] + r3 / 2 * (T_m)[{x, y, z + 1}]
			- r1 * (Ts_m)[{x, y, z}] - r2 * (Tss_m)[{x, y, z}]
			+ (1 - r1 - r2 - r3) * (T_m)[{x, y, z}];


		//if (d[i] > 22 || d[i] < 20)
		//{
		//	//printf("falsch  %d  % 0.3f \n" , tid, d[i]);
		//	printf("falsch  %d  % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, \n", tid, (T_m)[{x - 1, y, z}], (T_m)[{x + 1, y, z}], (Ts_m)[{x - 1, y, z}], (Ts_m)[{x + 1, y, z}], (T_m)[{x, y - 1, z}], (T_m)[{x, y + 1, z}], (T_m)[{x, y, z - 1}], (T_m)[{x, y, z + 1}], (Ts_m)[{x, y, z}], (T_m)[{x, y, z}], d[i]);
		//}
		////printf(" %.3f ", d[i]);

	}

	a[0] = 0;
	c[dim.z - 1] = 0;



	if (_par_d->boundCond == 0)//neumann -gleiche t_-1 = t_0
	{
		b[0] -= r3 / 2;
		b[dim.z - 1] -= r3 / 2;
	}
	else
	{
		d[0] += r3 / 2 * defaultValue;
		d[dim.z - 1] += r3 / 2 * defaultValue;
	}

	thomas(a, b, c, d, dim.z);
	float P_t;
	for (int i = 0; i < dim.z; i++)
	{

		int z = i;
		//_Tss_d[(tid * _dim.y) + i] = d[i];
		//(Ts_m)[{x, y, z}] = (T_m)[{x, y, z}];
		P_t = (wb_m)[{x, y, z}] * _par_d->cb * (_par_d->Ta - d[z]) * (getTd(_par_d, (T_m)[{x, y, z}]) / _par_d->k);
	//	printf("Pn %.3f wb %.3f \n", P_t, (wb_m)[{x, y, z}]);
		Tnew_m.setValue({ x,y,z }, d[z] + dt * P_t);
		
	}

}
__global__ void ySweep(float* _T_d, float* _Ts_d, float* _Tss_d, Parameters* _par_d)
{
	extern __shared__ float arrY[];

	dim3 dim = { _par_d->dim.x, _par_d->dim.y, _par_d->dim.z };
	float defaultValue = _par_d->defaultT;
	float r1, r2, r3;
	float td;
	float dt = _par_d->dt;
	float dx_2[3] = { _par_d->dx_2[0], _par_d->dx_2[1], _par_d->dx_2[2] };

	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= (dim.x * dim.z))
		return;

	
	float* a = &arrY[threadIdx.x * NUMBER_OF_ARRAYS * dim.y + 0];
	float* b = &arrY[threadIdx.x * NUMBER_OF_ARRAYS * dim.y + dim.y];
	float* c = &arrY[threadIdx.x * NUMBER_OF_ARRAYS * dim.y + dim.y * 2];
	float* d = &arrY[threadIdx.x * NUMBER_OF_ARRAYS * dim.y + dim.y * 3];

	Matrix3D<float> T_m = Matrix3D<float>(dim.x, dim.y, dim.z, _T_d);
	Matrix3D<float> Ts_m = Matrix3D<float>(dim.x, dim.y, dim.z, _Ts_d);
	Matrix3D<float> Tss_m = Matrix3D<float>(dim.x, dim.y, dim.z, _Tss_d);
	T_m.setAmbientValue(defaultValue);
	Ts_m.setAmbientValue(defaultValue);



	
	int x = tid % dim.x;
	int z = (int) (tid / dim.x);
	
	for (int i = 0; i < dim.y; i++)
	{
		int y = i;
		/*if (x >= dim.x || x < 0 ||
			y >= dim.y || y < 0 ||
			z >= dim.z || z < 0)
			printf("%d , %d , %d \n", x, y, z);*/

		td = getTd(_par_d, (T_m)[{x, y, z}]);

		r1 = r2 = r3 = td * dt / dx_2[1];

		a[i] = -r2 / 2;
		b[i] = 1 + r2;
		c[i] = -r2 / 2;

		d[i] = r1 / 2 * (T_m)[{x - 1, y, z}] + r1 / 2 * (T_m)[{x + 1, y, z}]
			+ r1 / 2 * (Ts_m)[{x - 1, y, z}] + r1 / 2 * (Ts_m)[{x + 1, y, z}]
			+ r2 / 2 * (T_m)[{x, y - 1, z}] + r2 / 2 * (T_m)[{x, y + 1, z}]
			+ r3 * (T_m)[{x, y, z - 1}] + r3 * (T_m)[{x, y, z + 1}]
			- r1 * (Ts_m)[{x, y, z}]
			+ (1 - r1 - r2 - 2 * r3) * (T_m)[{x, y, z}];


		//if (d[i] > 22 || d[i] < 20)
		//{
		//	//printf("falsch  %d  % 0.3f \n" , tid, d[i]);
		//	printf("falsch  %d  % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, \n", tid, (T_m)[{x - 1, y, z}], (T_m)[{x + 1, y, z}], (Ts_m)[{x - 1, y, z}], (Ts_m)[{x + 1, y, z}], (T_m)[{x, y - 1, z}], (T_m)[{x, y + 1, z}], (T_m)[{x, y, z - 1}], (T_m)[{x, y, z + 1}], (Ts_m)[{x, y, z}], (T_m)[{x, y, z}], d[i]);
		//}
		////printf(" %.3f ", d[i]);

	}

	a[0] = 0;
	c[dim.y - 1] = 0;



	if (_par_d->boundCond == 0)//neumann -gleiche t_-1 = t_0
	{
		b[0] -= r2 / 2;
		b[dim.y - 1] -= r2 / 2;
	}
	else
	{
		d[0] += r2 / 2 * defaultValue;
		d[dim.y - 1] += r2 / 2 * defaultValue;
	}

	thomas(a, b, c, d, dim.y);

	for (int i = 0; i < dim.y; i++)
	{

		int y = i;
		//_Tss_d[(tid * dim.y) + i] = d[i];
		//(Tss_m)[{x, y, z}] = (T_m)[{x, y, z}];
		Tss_m.setValue({ x,y,z }, d[i]);
	/*	if (d[i] > 100 || d[i] < 0)
		{
			printf("falsch Thomas y % 0.3f \n", d[i]);
		}*/
	}

}
__global__ void xSweep(float * _T_d, float* _Ts_d, Parameters* _par_d)
{
	extern __shared__ float arr[];

	dim3 dim = { _par_d->dim.x, _par_d->dim.y, _par_d->dim.z };
	//printf("Hallo %d ,%d , %d ", dim.x, dim.y, dim.z);
	
	float defaultValue = _par_d->defaultT;
	float r1, r2, r3;
	float td;
	float dt = _par_d->dt;
	float dx_2[3] = { _par_d->dx_2[0], _par_d->dx_2[1], _par_d->dx_2[2] };

	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= (dim.y * dim.z))
		return;



	float* a = &arr[threadIdx.x * NUMBER_OF_ARRAYS * dim.x + 0];
	float* b = &arr[threadIdx.x * NUMBER_OF_ARRAYS * dim.x + dim.x];
	float* c = &arr[threadIdx.x * NUMBER_OF_ARRAYS * dim.x + dim.x * 2];
	float* d = &arr[threadIdx.x * NUMBER_OF_ARRAYS * dim.x + dim.x * 3];
	
	Matrix3D<float> T_m = Matrix3D<float>(dim.x, dim.y, dim.z,_T_d);
	Matrix3D<float> Ts_m = Matrix3D<float>(dim.x, dim.y, dim.z, _Ts_d);
	T_m.setAmbientValue(defaultValue);
	

	
	//int y = (int)(tid / dim.z);
	//int z = tid % dim.z;
	/*if (tid >= 3600)
		printf("TID %d \n ", tid);*/
	int y = tid % dim.y;
	int z = (int)(tid / dim.y);

	for (int i = 0; i < dim.x; i++)
	{
		
		int x = i;
		/*if (x >= dim.x || x < 0 ||
			y >= dim.y || y < 0 ||
			z >= dim.z || z < 0)
			printf("%d , %d , %d \n", x, y, z);*/
		td = getTd(_par_d, (T_m)[{x, y, z}]);
		
		
		r1 = r2 = r3 = td * dt / dx_2[0];

		a[i] = -r1 / 2;
		b[i] = 1 + r1;
		c[i] = -r1 / 2;
		//printf("td %f, dt %f, dx %f ,r1 %f \n", td, dt, dx_2[0], r1);

		float T_im, T_ip, T_jm, T_jp,T_km, T_kp;

		//d[i] = (T_m)[{x, y, z}];
		d[i] = r1 / 2 * (T_m)[{x - 1, y, z}] + r1 / 2 * (T_m)[{x, y, z}]
		+ r2 * (T_m)[{x, y - 1, z}] + r2 * (T_m)[{x, y + 1, z}]
		+ r3 * (T_m)[{x, y, z - 1}] + r3 * (T_m)[{x, y, z + 1}]
		+ (1 - r1 - 2 * r2 - 2 * r3) * (T_m)[{x, y, z}];
		/*if (d[i] != 21.0)
			printf("ahhhh %.6f \n", d[i]);*/ //das noch rausfinden warum nicht gneau 21
		//if (d[i] > 22 || d[i] < 20)
		//{
		//	printf("falsch  %d  % 0.3f \n" , tid, d[i]);
		//	//printf("falsch  %d  % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f, % 0.3f \n", tid, _T_d[im], _T_d[ip], _T_d[jm], _T_d[jp], _T_d[km], _T_d[kp], _T_d[tid * dim.x + i], d[i]);
		//}
		//printf(" %.3f ", d[i]);
			
	}
	
	a[0] = 0;
	c[dim.x - 1] = 0;

	//printf("snow %0.3f \n ", (Ts_m)[{1, 1, 1}]);
	////(Ts_m).setValue({1, 1, 1},15);
	//Ts_m.setAll2(15.0);
	//printf("set %0.3f \n", (Ts_m)[{1, 1, 1}]);
	//Ts_m.setAll2(0.0);
	if (_par_d->boundCond == 0)//neumann -gleiche t_-1 = t_0
	{
		b[0] -= r1 / 2;
		b[dim.x - 1] -= r1 / 2;
	}
	else
	{
		d[0] += r1 / 2 * defaultValue;
		d[dim.x - 1] += r1 / 2 * defaultValue;
	}

	thomas(a, b, c, d, dim.x);
	
	for (int i = 0; i < dim.x; i++)
	{
		
		int x = i;
		//_Ts_d[(tid * dim.x) + i] = d[i];
		//(Ts_m)[{x, y, z}] = d[i];
		Ts_m.setValue({ x,y,z }, d[i]);
		//if (d[i] > 100 || d[i] < 0)
		//{
		//	printf("falsch Thomas z  %0.3f \n", d[i]);
		//}
	}

	

}
#else
float getTd(Parameters* _p, float _T)
{
	//(td_a + td_b * 1e-10 * exp(td_c * 1e-1 * T)) * 1e-6;
	//(td_a + T * td_b) * 1e-6

	return float(_p->td_a * 1e-6);

}
void thomas(float* a, float* b, float* c, float* d, int n) {
	/*
	// n is the number of unknowns

	|b0 c0 0 ||x0| |d0|
	|a1 b1 c1||x1|=|d1|
	|0  a2 b2||x2| |d2|

	1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

		x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

	2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
		from 1st it.: -| a1x0 + a1g0x1        = a1r0
					-----------------------------
						  (b1 - a1g0)x1 + c1x2 = d1 - a1r0

		x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

	3rd iteration:      | a2x1 + b2x2   = d2
		from 2nd it. : -| a2x1 + a2g1x2 = a2r2
					   -----------------------
					   (b2 - a2g1)x2 = d2 - a2r2
		x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
	Finally we have a triangular matrix:
	|1  g0 0 ||x0| |r0|
	|0  1  g1||x1|=|r1|
	|0  0  1 ||x2| |r2|

	Condition: ||bi|| > ||ai|| + ||ci||

	in this version the c matrix reused instead of g
	and             the d matrix reused instead of r and x matrices to report results
	Written by Keivan Moradi, 2014
	*/
	n--; // since we start from x0 (not x1)
	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++) {
		c[i] /= b[i] - a[i] * c[i - 1];
		d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
	}

	d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

	for (int i = n; i-- > 0;) {

		d[i] -= c[i] * d[i + 1];
	}
}
#endif

PennesHeatTransfer_Gpu::PennesHeatTransfer_Gpu()
{
#if GPU
	handle = 0;
	bool ft = InitDevice(handle);
	std::cout <<"initDevice: " << (ft ? "True" : "False") << std::endl;
#endif

}
PennesHeatTransfer_Gpu::PennesHeatTransfer_Gpu(unsigned int *_dim)
{
	par_h = (Parameters*)malloc(sizeof(Parameters));
	par_h->dim = { _dim[0],_dim[1],_dim[2] };
	n = par_h->dim.x * par_h->dim.y * par_h->dim.z;

	std::cout << "param0 " << par_h->dim.x << " " << par_h->dim.y << " " << par_h->dim.z << std::endl;

	
	handle = 0;


	

	T_h = new Matrix3D<float>(par_h->dim.x, par_h->dim.y, par_h->dim.z);
	Tnew_h = new Matrix3D<float>(par_h->dim.x, par_h->dim.y, par_h->dim.z);
	wb_h = new Matrix3D<float>(par_h->dim.x, par_h->dim.y, par_h->dim.z);
	//T_h->setAll2(GlobalVariables::boundaryCondition);
	//setDefaultParameters();
#if GPU
	bool ft = InitDevice(handle);
	std::cout << "initDevice: " << (ft ? "True" : "False") << std::endl;

	checkErrorsCuda(cudaMalloc((void**)&T_d, sizeof(float) * n));
	checkErrorsCuda(cudaMalloc((void**)&Ts_d, sizeof(float) * n));
	checkErrorsCuda(cudaMalloc((void**)&Tss_d, sizeof(float) * n));
	checkErrorsCuda(cudaMalloc((void**)&Tnew_d, sizeof(float) * n));

	checkErrorsCuda(cudaMalloc((void**)&wb_d, sizeof(float) * n));
	std::cout << "sizeofPara" << sizeof(Parameters) << std::endl;
	checkErrorsCuda(cudaMalloc((void**)&par_d, sizeof(Parameters)));


	checkLastCudaError("Kernel launch failed.");
#endif



}
PennesHeatTransfer_Gpu::PennesHeatTransfer_Gpu(const PennesHeatTransfer_Gpu& pht)
{
	par_h = (Parameters*)malloc(sizeof(Parameters));
	//par_h->dim = {pht.par_h->dim.x,pht.par_h->dim.y,pht.par_h->dim.z };
	(*par_h) = Parameters((*pht.par_h));
	n = par_h->dim.x * par_h->dim.y * par_h->dim.z;
	

	std::cout << "param0 " << par_h->dim.x << " " << par_h->dim.y << " " << par_h->dim.z << std::endl;


	handle = 0;


	

	T_h = new Matrix3D<float>(par_h->dim.x, par_h->dim.y, par_h->dim.z);
	Tnew_h = new Matrix3D<float>(par_h->dim.x, par_h->dim.y, par_h->dim.z);
	Tsaved_h = pht.Tsaved_h;
	wb_h = new Matrix3D<float>(*pht.wb_h);

#if GPU
	bool ft = InitDevice(handle);
	std::cout << "initDevice: " << (ft ? "True" : "False") << std::endl;
	checkErrorsCuda(cudaMalloc((void**)&T_d, sizeof(float) * n));
	checkErrorsCuda(cudaMalloc((void**)&Ts_d, sizeof(float) * n));
	checkErrorsCuda(cudaMalloc((void**)&Tss_d, sizeof(float) * n));
	checkErrorsCuda(cudaMalloc((void**)&Tnew_d, sizeof(float) * n));

	checkErrorsCuda(cudaMalloc((void**)&par_d, sizeof(Parameters)));
	checkErrorsCuda(cudaMalloc((void**)&wb_d, sizeof(float) * n));

	const int NUM_THREADS_PER_BLOCK = 1024;
	int num_of_threads = n;
	int num_threads_per_block = NUM_THREADS_PER_BLOCK;
	int num_blocks = num_of_threads / num_threads_per_block;
	if (n % num_threads_per_block != 0)
		num_blocks++;

	
	copyT << < num_blocks, num_threads_per_block >> > (pht.wb_d, wb_d, n);
	cudaDeviceSynchronize();
	checkLastCudaError("Kernel launch failed.");
#endif
}
PennesHeatTransfer_Gpu::~PennesHeatTransfer_Gpu()
{
#if GPU
	checkErrorsCuda(cudaFree(T_d));
	checkErrorsCuda(cudaFree(Tnew_d));
	
	checkErrorsCuda(cudaFree(Ts_d));
	checkErrorsCuda(cudaFree(Tss_d));
	checkErrorsCuda(cudaFree(par_d));
	checkErrorsCuda(cudaFree(wb_d));
#endif
	delete T_h;
	delete Tnew_h;
	delete Tsaved_h;
	delete wb_h;

}
void PennesHeatTransfer_Gpu::setT_h(Matrix3D<float>* _T_h)
{

	delete T_h;
	T_h = new Matrix3D<float>(*_T_h);
#if	DEBUG
	long float sum = 0;
	for (int i = 0; i < n; i++)
	{
		sum += T_h->getData()[i];
		//std::cout << _Tnew_h[i] << " ";

	}
	printf("T_h: %f avag %f \n", sum, sum / (n));
#endif
#if GPU
	std::chrono::steady_clock::time_point begin_memcpy_malloc = std::chrono::steady_clock::now();
	checkErrorsCuda(cudaMemcpy(T_d, T_h->getData(), n * sizeof(float), cudaMemcpyHostToDevice));
	std::chrono::steady_clock::time_point end_memcpy_malloc = std::chrono::steady_clock::now();
	std::cout << "memcpy_malloc: " << std::chrono::duration_cast<std::chrono::microseconds>(end_memcpy_malloc - begin_memcpy_malloc).count() << std::endl;
#endif

}
#if GPU
void PennesHeatTransfer_Gpu::setT_d(float* _Tsaved_d)
{
	const int NUM_THREADS_PER_BLOCK = 1024;
	int num_of_threads = n;
	int num_threads_per_block = NUM_THREADS_PER_BLOCK;
	int num_blocks = num_of_threads / num_threads_per_block;
	if (n % num_threads_per_block != 0)
		num_blocks++;

	/*std::cout << "coypTnew2T: num_blocks = " << num_blocks << " :: "
		<< "num_threads_per_block = " << num_threads_per_block << " :: num_of_threads: " << num_of_threads << std::endl;*/

	copyT << < num_blocks, num_threads_per_block >> > (_Tsaved_d, T_d, n);
	cudaDeviceSynchronize();
	checkLastCudaError("Kernel launch failed.");

	//Eigentlich nur für debug zwecke
	checkErrorsCuda(cudaMemcpy(T_h->getData(), T_d, sizeof(float) * n, cudaMemcpyDeviceToHost));
	long float sum = 0;

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		sum += T_h->getData()[i];
		//std::cout << _Tnew_h[i] << " ";

	}
	//std::cout <<std::endl<< "SumIs: " << sum << std::endl;// (dim.x * dim.y * dim.z) 
	printf("T_h: %f avag %f \n", sum, sum / (n));
#endif
}
#endif
void PennesHeatTransfer_Gpu::setAll2(float _T) 
{ 
	T_h->setAll2(_T); 
	Tnew_h->setAll2(_T);
#if GPU
	std::chrono::steady_clock::time_point begin_memcpy_malloc = std::chrono::steady_clock::now();
	checkErrorsCuda(cudaMemcpy(T_d, T_h->getData(), n * sizeof(float), cudaMemcpyHostToDevice));
	checkErrorsCuda(cudaMemcpy(Tnew_d, Tnew_h->getData(), n * sizeof(float), cudaMemcpyHostToDevice));
	std::chrono::steady_clock::time_point end_memcpy_malloc = std::chrono::steady_clock::now();
	std::cout << "memcpy_malloc: " << std::chrono::duration_cast<std::chrono::microseconds>(end_memcpy_malloc - begin_memcpy_malloc).count() << std::endl;
#endif
}
void PennesHeatTransfer_Gpu::updateTnew_h()
{
	//std::chrono::steady_clock::time_point begin_memcpy = std::chrono::steady_clock::now();
#if GPU	
	checkErrorsCuda(cudaMemcpy(Tnew_h->getData(), Tnew_d, sizeof(float) * n, cudaMemcpyDeviceToHost));
#endif
	/*std::chrono::steady_clock::time_point end_memcpy = std::chrono::steady_clock::now();
	std::cout << "memcpy: " << std::chrono::duration_cast<std::chrono::microseconds>(end_memcpy - begin_memcpy).count() << std::endl;*/
#if DEBUG
	long float sum = 0;

	for (int i = 0; i < n; i++)
	{
		sum += Tnew_h->getData()[i];
		//std::cout << _Tnew_h[i] << " ";

	}
	//std::cout <<std::endl<< "SumIs: " << sum << std::endl;// (dim.x * dim.y * dim.z) 

	printf("Tnew_h: %f avag %f \n", sum, sum/(n));
#endif

	//return Tnew_h;
}
//void PennesHeatTransfer_Gpu::configTimeSteps(float _dt, int _N = 1)
//{
//
//}
void PennesHeatTransfer_Gpu::saveTimestep()
{
#if GPU
	if (Tsaved_d == nullptr)
	{
		checkErrorsCuda(cudaFree(PennesHeatTransfer_Gpu::Tsaved_d));
		checkErrorsCuda(cudaMalloc((void**)&Tsaved_d, sizeof(float) * n));
		Tsaved_h = new Matrix3D<float>(par_h->dim.x, par_h->dim.y, par_h->dim.z);
	}


	const int NUM_THREADS_PER_BLOCK = 1024;
	int num_of_threads = n;
	int num_threads_per_block = NUM_THREADS_PER_BLOCK;
	int num_blocks = num_of_threads / num_threads_per_block;
	if (n % num_threads_per_block != 0)
		num_blocks++;


	copyT << < num_blocks, num_threads_per_block >> > (Tnew_d,Tsaved_d, n);
	cudaDeviceSynchronize();
	checkLastCudaError("Kernel launch failed.");

	//DEbug
	
	/*std::chrono::steady_clock::time_point end_memcpy = std::chrono::steady_clock::now();
	std::cout << "memcpy: " << std::chrono::duration_cast<std::chrono::microseconds>(end_memcpy - begin_memcpy).count() << std::endl;*/
#else
	if (Tsaved_h == nullptr)
	{
		Tsaved_h = new Matrix3D<float>(par_h->dim.x, par_h->dim.y, par_h->dim.z);
	}
	for (int i = 0; i < n; i++)
		Tsaved_h->setValue(i, Tnew_h->getData()[i]);
#endif

#if DEBUG
#if GPU
	checkErrorsCuda(cudaMemcpy(Tsaved_h->getData(), Tsaved_d, sizeof(float) * n, cudaMemcpyDeviceToHost));
#endif
	long float sum = 0;

	for (int i = 0; i < n; i++)
	{
		sum += Tsaved_h->getData()[i];
		//std::cout << _Tnew_h[i] << " ";

	}
	//std::cout <<std::endl<< "SumIs: " << sum << std::endl;// (dim.x * dim.y * dim.z) 
	printf("Tsaved_h: %f avag %f \n", sum, sum / (n));
#endif
}
#if GPU
void  PennesHeatTransfer_Gpu::updateParameters()
{
	//std::cout << "param " << par_h->dim.x << " "<< par_h->dim.y << " " << par_h->dim.z << std::endl;

	checkErrorsCuda(cudaMemcpy(par_d, par_h, sizeof(Parameters), cudaMemcpyHostToDevice));
	//printParameters();
	/*cuda_hello << <1,1>> > (par_d);
	cudaDeviceSynchronize();
	checkLastCudaError("Kernel launch failed.");
	std::cout << "hallo" << std::endl;*/

}
#endif
void PennesHeatTransfer_Gpu::setDefaultParameters()
{
	par_h->td_a = 1.5;//0.155;
	par_h->td_b = 0.0;
	par_h->td_c = 0.0;
	//par_h->td_b = 4.95;
	//par_h->td_c = 2.01;
	par_h->boundCond = GlobalVariables::boundaryCondition;
	par_h->cb = 4182;
	par_h->k = 2;

	par_h->Ta = 25;
	par_h->defaultT = GlobalVariables::baseLineT;
	par_h->dt = 2.0;

	par_h->dx_2[0] = 0.0025;// 4.0 * 1e-6;
	par_h->dx_2[1] = 0.0025; //4.0 * 1e-6;
	par_h->dx_2[2] = 0.0025; //4.0 * 1e-6;

	for (int i = 0; i < BUFFER_HEATSOURCES; i++)
	{
		par_h->T_heat[i] = (float)GlobalVariables::baseLineT;
		par_h->heat_x[i] = 0;
		par_h->heat_y[i] = 0;
		par_h->heat_z[i] = 0;
	}
	par_h->N_heat = 0;
#if GPU
	updateParameters();
#endif
}
void PennesHeatTransfer_Gpu::setPerfusion(Matrix3D<bool>& _map)
{
	std::cout << "SetPerfusion " <<n<< std::endl;
	double w_b = 800 / (1e6 * 60) * 997 / (0.005 * 0.005 * 1 * 3.1415);
	
	for (int i = 0; i < n; i++)
	{
		if (_map.getData()[i] == false)
		{
			//wb_h->setValue(i, 0);
			wb_h->getData()[i] = 0;
		}
		else
		{
			//wb_h->setValue(i, w_b);
			wb_h->getData()[i] = w_b;
			//std::cout << "WB" << std::endl;
		}
		//if(wb_h->getData()[i] != 0)
		//	std::cout << " " << wb_h->getData()[i];
	}
#if GPU
	std::cout << "ersterErlfog" << std::endl;
	std::chrono::steady_clock::time_point begin_memcpy_malloc = std::chrono::steady_clock::now();
	checkErrorsCuda(cudaMemcpy(wb_d, wb_h->getData(), sizeof(float)* n, cudaMemcpyHostToDevice));
	std::chrono::steady_clock::time_point end_memcpy_malloc = std::chrono::steady_clock::now();
	std::cout << "memcpy: " << std::chrono::duration_cast<std::chrono::microseconds>(end_memcpy_malloc - begin_memcpy_malloc).count() << std::endl;
	cudaDeviceSynchronize();
	checkLastCudaError("Kernel launch failed.");
#endif
}
void PennesHeatTransfer_Gpu::finiteStep(int _NrOfTimeSteps)
{
	//printParameters();
#if GPU
	finiteStep_GPU(_NrOfTimeSteps);
#else 
	finiteStep_CPU(_NrOfTimeSteps);
#endif 
}
//
//
#if GPU
void  PennesHeatTransfer_Gpu::finiteStep_GPU(int _NrOfTimeSteps)
{
	//std::cout << "finiteStep " << _NrOfTimeSteps << std::endl;
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	
	updateParameters();
	//checkErrorsCuda(cudaMemcpy(par_d, (void*)&par_h, sizeof(PennesEquationParameter), cudaMemcpyHostToDevice));


	const int NUM_THREADS_PER_BLOCK = 32;//optimerung damit chache voll ausgefüllt (ist aber glaube ich lagnsamer)
	dim3 num_blocks, num_threads_per_block, number_of_threads;

	//xSweep
	number_of_threads.x = par_h->dim.y * par_h->dim.z;
	num_threads_per_block.x = NUM_THREADS_PER_BLOCK;
	num_blocks.x = number_of_threads.x / num_threads_per_block.x;
	if (0 != number_of_threads.x % num_threads_per_block.x) {
		num_blocks.x++;
	}
	/*std::cout << "X: num_blocks = " << num_blocks.x << " :: "
		<< "num_threads_per_block = " << num_threads_per_block.x << " :: num_of_threads: " << number_of_threads.x<< std::endl;*/

		//ySweep
	number_of_threads.y = par_h->dim.x * par_h->dim.z;
	num_threads_per_block.y = NUM_THREADS_PER_BLOCK;
	num_blocks.y = number_of_threads.y / num_threads_per_block.y;
	if (0 != number_of_threads.y % num_threads_per_block.y) {
		num_blocks.y++;
	}
	/*std::cout << "Y: num_blocks = " << num_blocks.y << " :: "
		<< "num_threads_per_block = " << num_threads_per_block.y << " :: num_of_threads: " << number_of_threads.y << std::endl;*/
		//zSweep
	number_of_threads.z = par_h->dim.x * par_h->dim.y;
	num_threads_per_block.z = NUM_THREADS_PER_BLOCK;
	num_blocks.z = number_of_threads.z / num_threads_per_block.z;
	if (0 != number_of_threads.z % num_threads_per_block.z) {
		num_blocks.z++;
	}
	/*std::cout << "Z: num_blocks = " << num_blocks.z << " :: "
		<< "num_threads_per_block = " << num_threads_per_block.z << " :: num_of_threads: " << number_of_threads.z << std::endl;*/

	for (int i = 0; i < _NrOfTimeSteps; i++)
	{


		//std::cout << "setHeatSources" << std::endl;
		setHeatSources << <1, par_h->N_heat >> > (T_d, par_d);
		cudaDeviceSynchronize();
		checkLastCudaError("Kernel launch failed.");
		// MBTotal shared memory per block: 49152

		//Execute XSweep
		//std::cout << "Execute XSweep" << std::endl;
		//std::chrono::steady_clock::time_point begin_xSweep = std::chrono::steady_clock::now();

		int reservedSharedMemoryX = (num_threads_per_block.x * NUMBER_OF_ARRAYS * par_h->dim.x * sizeof(float));
		//std::cout << "sharedMemoryX: " << reservedSharedMemoryX << std::endl;;
		xSweep << <num_blocks.x, num_threads_per_block.x, reservedSharedMemoryX >> > (T_d, Ts_d, par_d);
		cudaDeviceSynchronize();
		checkLastCudaError("Kernel launch failed.");
		/*std::chrono::steady_clock::time_point end_xSweep = std::chrono::steady_clock::now();
		std::cout << "xSweep_duration: " << std::chrono::duration_cast<std::chrono::microseconds>(end_xSweep - begin_xSweep).count() << std::endl;*/

		//Execute YSweep
		//std::cout << "Execute YSweep" << std::endl;
		//std::chrono::steady_clock::time_point begin_ySweep = std::chrono::steady_clock::now();
		int reservedSharedMemoryY = (num_threads_per_block.y * NUMBER_OF_ARRAYS * par_h->dim.y * sizeof(float));
		//std::cout << "sharedMemoryY: " << reservedSharedMemoryY << std::endl;;
		ySweep << <num_blocks.y, num_threads_per_block.y, reservedSharedMemoryY >> > (T_d, Ts_d, Tss_d, par_d);
		cudaDeviceSynchronize();
		checkLastCudaError("Kernel launch failed.");
		/*std::chrono::steady_clock::time_point end_ySweep = std::chrono::steady_clock::now();
		std::cout << "ySweep_duration: " << std::chrono::duration_cast<std::chrono::microseconds>(end_ySweep - begin_ySweep).count() << std::endl;*/

		////Execute ZSweep
		//std::cout << "Execute ZSweep" << std::endl;
		//std::chrono::steady_clock::time_point begin_zSweep = std::chrono::steady_clock::now();
		int reservedSharedMemoryZ = (num_threads_per_block.z * NUMBER_OF_ARRAYS * par_h->dim.z * sizeof(float));
		//std::cout << "sharedMemoryZ: " << reservedSharedMemoryZ<< std::endl;;
		zSweep << <num_blocks.z, num_threads_per_block.z, reservedSharedMemoryZ >> > (T_d, Ts_d, Tss_d, Tnew_d, wb_d, par_d);
		cudaDeviceSynchronize();
		checkLastCudaError("Kernel launch failed.");
		/*std::chrono::steady_clock::time_point end_zSweep = std::chrono::steady_clock::now();
		std::cout << "zSweep_duration: " << std::chrono::duration_cast<std::chrono::microseconds>(end_zSweep - begin_zSweep).count() << std::endl;*/

		setHeatSources << <1, par_h->N_heat >> > (T_d, par_d);
		cudaDeviceSynchronize();
		checkLastCudaError("Kernel launch failed.");

		if (_NrOfTimeSteps > 1)
		{
			const int NUM_THREADS_PER_BLOCK = 1024;
			int num_of_threads = n;
			int num_threads_per_block = NUM_THREADS_PER_BLOCK;
			int num_blocks = num_of_threads / num_threads_per_block;
			if (n % num_threads_per_block != 0)
				num_blocks++;

			/*std::cout << "coypTnew2T: num_blocks = " << num_blocks << " :: "
				<< "num_threads_per_block = " << num_threads_per_block << " :: num_of_threads: " << num_of_threads << std::endl;*/

			copyT << < num_blocks, num_threads_per_block >> > (Tnew_d, T_d, n);
			cudaDeviceSynchronize();
			checkLastCudaError("Kernel launch failed.");

		}



	}

	//DEBUG



	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	//	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << std::endl;// "finiteStep_duration: " <<
}
#endif
void  PennesHeatTransfer_Gpu::finiteStep_CPU(int _NrOfTimeSteps)
{
	//std::cout << "StepImlicitADI-Start" << std::endl;

	//crank-Nicolson - von nuemann adiabatisch randbedingung

	dim3 dim = { par_h->dim.x, par_h->dim.y, par_h->dim.z };
	float defaultValue = par_h->defaultT;
	float dt = par_h->dt;
	float dx_2[3] = { par_h->dx_2[0], par_h->dx_2[1], par_h->dx_2[2] };

	//std::cout << "T:mitte: " << (*P_n)[{dim.x / 2, dim[1] / 2, dim.z / 2}].T << std::endl;

	//Matrix3D<double>* T = new Matrix3D<double>(dim.x, dim[1], dim.z);

	
	Matrix3D<float> Ts_m = Matrix3D<float>(par_h->dim.x, par_h->dim.y, par_h->dim.z);
	Matrix3D<float> Tss_m = Matrix3D<float>(par_h->dim.x, par_h->dim.y, par_h->dim.z);

	
	if (GlobalVariables::boundaryCondition == 1)
	{
		T_h->setAmbientValue(defaultValue);
		Ts_m.setAmbientValue(defaultValue);
		Tss_m.setAmbientValue(defaultValue);
		
	}
	for (int t = 0; t < _NrOfTimeSteps; t++)
	{
		//Heating
	//std::cout << "Heating Temperatur: " << std::endl;
#pragma omp parallel for 
		for (int i = 0; i < par_h->N_heat; i++)
		{

			T_h->setValue({ par_h->heat_x[i],par_h->heat_y[i] ,par_h->heat_z[i] }, par_h->T_heat[i]);
		}
		//x-direction	

//#pragma omp parallel for collapse(2)
		for (int y = 0; y < dim.y; y++)
		{
			for (int z = 0; z < dim.z; z++)
			{


				float* a = new float[dim.x]();
				float* b = new float[dim.x]();
				float* c = new float[dim.x]();
				float* d = new float[dim.x]();

				float r1, r2, r3;
				float td; //thermal diffusivity (m^2/s)

				for (int x = 0; x < dim.x; x++)
				{



					//TODO k flux conservative 
						//TODO boundary conditition
					//td = (*P_n)[{i, j, k}].getK() / ((*P_n)[{i, j, k}].c* (*P_n)[{i, j, k}].rho);
					//td = (*P_n)[{i, j, k}].k / (*P_n)[{i, j, k}].M;
					td = getTd(par_h, (*T_h)[{x, y, z}]);
					r1 = r2 = r3 = td * dt / dx_2[0];

					a[x] = -r1 / 2.0;
					b[x] = 1.0 + r1;
					c[x] = -r1 / 2.0;


					d[x] = r1 / 2.0 * (*T_h)[{x - 1, y, z}] + r1 / 2.0 * (*T_h)[{x, y, z}]
						+ r2 * (*T_h)[{x, y - 1, z}] + r2 * (*T_h)[{x, y + 1, z}]
						+ r3 * (*T_h)[{x, y, z - 1}] + r3 * (*T_h)[{x, y, z + 1}]
						+ (1.0 - r1 - 2.0 * r2 - 2.0 * r3) * (*T_h)[{x, y, z}];

				}
				a[0] = 0;
				c[dim.x - 1] = 0;
				if (GlobalVariables::boundaryCondition == 0)//neumann -gleiche t_-1 = t_0
				{
					b[0] -= r1 / 2.0;
					b[dim.x - 1] -= r1 / 2.0;
				}
				else
				{

					d[0] += r1 / 2.0 * defaultValue;

					d[dim.x - 1] += r1 / 2.0 * defaultValue;
					/*d[0] += r3 / 2 * defaultValue + r2 / 2 * defaultValue +  r1 / 2 * defaultValue;
					d[par_h->dim.x - 1] += r3 / 2 * defaultValue + r2 / 2 * defaultValue +  r1 / 2 * defaultValue;*/
				}
				//d[0] = d[0] + a[0] * (*P_n)[{0, j, k}].T;
				//d[par_h->dim.x - 1] = d[par_h->dim.x - 1] + c[0] * (*P_n)[{par_h->dim.x - 1, j, k}].T;
				thomas(a, b, c, d, par_h->dim.x);

				for (int x = 0; x < par_h->dim.x; x++)
				{
				/*	if (d[x] < 20)
						std::cout << "dx :" <<x << " " << y << " " << z << " : " << d[x] << std::endl;*/
					(Ts_m)[{x, y, z}] = d[x];
				}
				delete[] a;
				delete[] b;
				delete[] c;
				delete[] d;
			}
		}





		//y-Direction

//#pragma omp parallel for collapse(2)
		for (int x = 0; x < par_h->dim.x; x++)
		{
			for (int z = 0; z < par_h->dim.z; z++)
			{
				float* a = new float[par_h->dim.y]();
				float* b = new float[par_h->dim.y]();
				float* c = new float[par_h->dim.y]();
				float* d = new float[par_h->dim.y]();
				float r1, r2, r3;
				float td; //thermal diffusivity (m^2/s)

				for (int y = 0; y < par_h->dim.y; y++)
				{

					//td = (*P_n)[{i, j, k}].getK() / ((*P_n)[{i, j, k}].c* (*P_n)[{i, j, k}].rho);
					//td = (*P_n)[{i, j, k}].k / (*P_n)[{i, j, k}].M;
					td = getTd(par_h, (*T_h)[{x, y, z}]);
					r1 = r2 = r3 = td * dt / dx_2[1];

					a[y] = -r2 / 2.0;
					b[y] = 1.0 + r2;
					c[y] = -r2 / 2.0;

					d[y] = r1 / 2.0 * (*T_h)[{x - 1, y, z}] + r1 / 2.0 * (*T_h)[{x + 1, y, z}]
						+ r1 / 2.0 * (Ts_m)[{x - 1, y, z}] + r1 / 2.0 * (Ts_m)[{x + 1, y, z}]
						+ r2 / 2.0 * (*T_h)[{x, y - 1, z}] + r2 / 2.0 * (*T_h)[{x, y + 1, z}]
						+ r3 * (*T_h)[{x, y, z - 1}] + r3 * (*T_h)[{x, y, z + 1}]
						- r1 * (Ts_m)[{x, y, z}]
						+ (1.0 - r1 - r2 - 2.0 * r3) * (*T_h)[{x, y, z}];
				}
				a[0] = 0;
				c[par_h->dim.y - 1] = 0;
				//cout << "d[0]_y " << d[0] << std::endl << std::endl;
				if (GlobalVariables::boundaryCondition == 0)//neumann -gleiche t_-1 = t_0
				{
					b[0] -= r2 / 2.0;
					b[par_h->dim.y - 1] -= r2 / 2.0;
				}
				else
				{
					d[0] += r2 / 2.0 * defaultValue;
					d[par_h->dim.y - 1] += r2 / 2.0 * defaultValue;
					/*d[0] += r3 / 2 * defaultValue + r2 / 2 * defaultValue + 2 * r1 / 2 * defaultValue;
					d[par_h->dim.y - 1] += r3 / 2 * defaultValue + r2 / 2 * defaultValue + 2 * r1 / 2 * defaultValue;*/
				}


				thomas(a, b, c, d, par_h->dim.y);
				//cout << "d[0]_y " << d[0] << std::endl << std::endl;

				for (int y = 0; y < par_h->dim.y; y++)
				{
					/*if (d[l] < 20)
						std::cout << "dy :" << i << " " << l << " " << k << " : " << d[l] << std::endl;*/
						/*if ((*P_n)[{i, l, k}].keepTconstant)
							(*T_ss)[{i, l, k}] = (*T_s)[{i, l, k}];
						else*/
					(Tss_m)[{x, y, z}] = d[y];
				}
				delete[] a;
				delete[] b;
				delete[] c;
				delete[] d;
			}
		}



		//z.direction
//#pragma omp parallel for collapse(2)
		for (int x = 0; x < par_h->dim.x; x++)
		{
			for (int y = 0; y < par_h->dim.y; y++)
			{
				float* a = new float[par_h->dim.z]();
				float* b = new float[par_h->dim.z]();
				float* c = new float[par_h->dim.z]();
				float* d = new float[par_h->dim.z]();
				float r1, r2, r3;
				float td; //thermal diffusivity (m^2/s)
				float res;

				for (int z = 0; z < par_h->dim.z; z++)
				{

					//td = (*P_n)[{i, j, k}].getK() / ((*P_n)[{i, j, k}].c* (*P_n)[{i, j, k}].rho);
					//td = (*P_n)[{i, j, k}].k / (*P_n)[{i, j, k}].M;
					td = getTd(par_h, (*T_h)[{x, y, z}]);
					r1 = r2 = r3 = td * dt / dx_2[2];

					a[z] = -r3 / 2.0;
					b[z] = 1.0 + r3;
					c[z] = -r3 / 2.0;

					d[z] = r1 / 2.0 * (*T_h)[{x - 1, y, z}] + r1 / 2.0 * (*T_h)[{x + 1, y, z}]
						+ r1 / 2.0 * (Ts_m)[{x - 1, y, z}] + r1 / 2.0 * (Ts_m)[{x + 1, y, z}]
						+ r2 / 2.0 * (*T_h)[{x, y - 1, z}] + r2 / 2.0 * (*T_h)[{x, y + 1, z}]
						+ r2 / 2.0 * (Tss_m)[{x, y - 1, z}] + r2 / 2.0 * (Tss_m)[{x, y + 1, z}]
						+ r3 / 2.0 * (*T_h)[{x, y, z - 1}] + r3 / 2.0 * (*T_h)[{x, y, z + 1}]
						- r1 * (Ts_m)[{x, y, z}] - r2 * (Tss_m)[{x, y, z}]
						+ (1.0 - r1 - r2 - r3) * (*T_h)[{x, y, z}];
					/*		d[k] = r1 / 2 * (*T)[{i - 1, j, k}] + r1 / 2 * (*T)[{i + 1, j, k}]
								+ r1 / 2 * (*T_s)[{i - 1, j, k}] + r1 / 2 * (*T_s)[{i + 1, j, k}]
								+ r2 / 2 * (*T)[{i, j - 1, k}] + r2 / 2 * (*T)[{i, j + 1, k}]
								+ r2 / 2 * (*T_ss)[{i, j - 1, k}] + r2 / 2 * (*T_ss)[{i, j + 1, k}]
								+ r3 / 2 * (*T)[{i, j, k - 1}] + r3 / 2 * (*T)[{i, j, k + 1}]
								- r1 * (*T_s)[{i, j, k}] - r2 * (*T_ss)[{i, j, k}]
								+ (1 - r1 - r2 - r3) * (*T)[{i, j, k}];*/

				}
				//cout <<"d[0]_z "<< d[0]<<std::endl;
				a[0] = 0;
				c[par_h->dim.z - 1] = 0;

				if (GlobalVariables::boundaryCondition == 0)//neumann -gleiche t_-1 = t_0
				{

					b[0] -= r3 / 2.0;
					b[par_h->dim.z - 1] -= r3 / 2.0;
				}
				else
				{
					d[0] += r3 / 2.0 * defaultValue;
					//cout << "d[0] " << d[0] << std::endl;
					d[par_h->dim.z - 1] += r3 / 2.0 * defaultValue;
					/*d[0] += r3 / 2 * defaultValue + 2 * r2 / 2 * defaultValue + 2 * r1 / 2 * defaultValue;
					d[par_h->dim.z-1] += r3 / 2 * defaultValue + 2 * r2 / 2 * defaultValue + 2 * r1 / 2 * defaultValue;*/
				}


				//cout << "d[0] " << d[0] << std::endl;
				thomas(a, b, c, d, par_h->dim.z);
				//cout << "d[0]_z " << d[0] << std::endl << std::endl;
				float P_t;
				for (int z = 0; z < par_h->dim.z; z++)
				{

					//if (d[z] > 100 || d[z] < 10)
					//{
					//	std::cout << "Falsch d pen " <<x<<" "<<y<<" "<<z<<" " << d[z]<<std::endl;
					//}
					/*if(hf_t > 0.1|| hf_t < -0.1)
						std::cout  <<"hft_t: " <<hf_t << std::endl;*/
						/*if (d[l] < 20)
							std::cout<<"d :" << i << " " << j << " " << l <<" : "<<d[l]<<  std::endl;*/
							/*if (hf_t> 0)
								std::cout <<"hft: "<<hf_t<<" " << i << " " << j << " " << l << std::endl;*/
								/*if (d[l] > 100)
									std::cout << i<< " " << j << " " <<l << std::endl;*/
									//int v;
									//if (v_n.isVessel)
									//{
									//	v = 1;
									//	//std::cout << "isVessel" << std::endl;
									//}
									//	
									//else
									//	v = 0;
									//Pt = rho_b * w * Cb / (rho * Ct) * (T_amb - T_new)
					//float P_t = v_n.w_b * v_n.c_b * (v_n.T_a - d[l]) * (v_n.getTd() / v_n.k);
					/*if (P_t > 0)
						std::cout <<"P_t" << P_t << " Diff "<< (v_n.T_a - d[l])<<std::endl;*/
						//res = d[l] + dt * P_t;//+ hf_t   /((v_n.getTd()/v_n.k))
						//+dt * (v_n.w_b * v_n.c_b * (v_n.T_a - v_n.T))
						/*if (v_n.keepTconstant == true)
						{
							(*P_new)[{i, j, l}].T = v_n.T;
							(*T_new)[{i, j, l}] = v_n.T;

						}
						else*/
					
					P_t = (*wb_h)[{x, y, z}] * par_h->cb * (par_h->Ta - d[z]) *(getTd(par_h, (*T_h)[{x, y, z}]) / par_h->k);
						(*Tnew_h)[{x, y, z}] = d[z] + dt * P_t;
					

				}
				delete[] a;
				delete[] b;
				delete[] c;
				delete[] d;
			}
		}


		long float Heatsum = 0;

		if (_NrOfTimeSteps > 1)
		{
		//	std::cout << "multi" << std::endl;

			for (int i = 0; i < par_h->dim.x; i++)
			{
				for (int j = 0; j < par_h->dim.y; j++)
				{
					for (int k = 0; k < par_h->dim.z; k++)
					{

						(*T_h)[{i, j, k}] = (*Tnew_h)[{i, j, k}];

						Heatsum += (*T_h)[{i, j, k}];
						//if ((*P_new)[{i, j, k}].T <= 20 || (*P_new)[{i, j, k}].T >= 22)
						//{
						//	std::cout <<"null: " << i << " " << " " << j << " " << k <<" "<< (*P_new)[{i, j, k}].T << std::endl;
						//}

					}
				}
			}
		}

		//debugViewer(T_s);
		//ParameterOptimization::debugViewer(T_ss);
		//ParameterOptimization::debugViewer(P_new);


		//delete T;
		/*delete T_s;
		delete T_ss;
		delete T_new;*/
		//std::cout << "stepImplicitADI - Heatsum: "<<Heatsum<<" " << (double)(Heatsum / (P_n->getpar_h->dimensions()[0] * P_n->getpar_h->dimensions()[1] * P_n->getDimensions()[2])) << std::endl;
		//debugViewer(P_new);
		//std::cout << "stepImplicitADI-End" << std::endl;
	}
	



}

void PennesHeatTransfer_Gpu::printParameters()
{
	std::cout<< "dx_2 " << par_h->dx_2[0] << " " << par_h->dx_2[1] << " " << par_h->dx_2[2] << std::endl;
	std::cout << "dt " << par_h->dt << std::endl;
	std::cout << "td_a " << par_h->td_a << std::endl;
	std::cout << "defT " << par_h->defaultT << std::endl;
	std::cout << "bound " << par_h->boundCond << std::endl;
	std::cout << "dim " << par_h->dim.x<<" "<< par_h->dim.y<<" "<<par_h->dim.z << std::endl;
	std::cout << "N_heat " << par_h->N_heat << std::endl;
	for (int i = 0; i < par_h->N_heat; i++)
	{
		std::cout << par_h->T_heat[i] << " "<< par_h->heat_x[i]<<" "<< par_h->heat_y[i] << " " << par_h->heat_z[i] << std::endl;
	}
	
}



//wenn man mehere zeitschrritte will
//fkt(--int tsteps)
//atomic int currentstep
// while(current step < tsptes)
// {
// syncthreads
// currentstep++
// }
//
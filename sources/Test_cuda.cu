//cudaWorker.cu




#include "Test_cuda.h"
#include<cuda_util.h>
#include<TestClass.h>
//__global__ void kernel1(float deviceData[])
//
//{
//
//}
__global__ void cuda_hello(A a) {//float* _T_d, float* _Ts_d
	printf("Hello World from GPU!\n");
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	printf("tid %d: \n", tid);
	a.increment();
	a.print_data();
}
void Test_Cuda::doCuda()

{
	std::cout << "testCuda" << std::endl;
	int handle = 0;
	bool ft = InitDevice(handle);
	std::cout << (ft ? "True" : "False") << std::endl;
	const int n = 1000;
	float a[n];
	float b[n];
	float* a_d, * b_d;
	for (int i = 0; i < n; i++)
	{
		a[i] = 10;
		b[i] = 1;
	}
	//checkErrorsCuda(cudaMalloc(((void**)&a_d, n * sizeof(float)));

	/*checkErrorsCuda(cudaMalloc(((void**)&b_d, n * sizeof(float)));

	checkErrorsCuda(cudaMemcpy(a_d, a, cudaMemcpyHostToDevice));*/

	A h_a;
	h_a.increment();
	h_a.print_data();
	cuda_hello << <1, 1 >> > (h_a);
	cudaDeviceSynchronize();
	// 
	//unpack the vector into an array

	//use CUDA functions to to memcpy data and to launch a kernel

	//kernel1 << <64, 64 >> > (deviceData); //launching kernel legal, it's a .cu file

}
Test_Cuda::Test_Cuda()
{

}
Test_Cuda::~Test_Cuda()
{

}

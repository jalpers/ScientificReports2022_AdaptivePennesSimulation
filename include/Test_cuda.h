#ifndef CUDA_TEST_H
#define CUDA_TEST_H
#include<iostream>
#include <cuda_runtime.h>
#include <cuda.h>


class Test_Cuda
{

public:
	Test_Cuda();
	~Test_Cuda();
	void doCuda();

};

#endif // !CUDA_TEST_H
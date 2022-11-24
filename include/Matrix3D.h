#pragma once
//https://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
//https://codereview.stackexchange.com/questions/237962/template-matrix-class-implemented-some-basic-functionalities
//https://www.bestprog.net/en/2019/08/23/c-an-example-of-creating-a-template-class-matrix-dynamic-memory-allocation/
//https://stackoverflow.com/questions/30448976/c-opencv-creating-a-3d-matrix-and-access-its-elements
//https://stackoverflow.com/questions/9498384/how-delete-a-pointer-of-classes-which-has-pointer-members
//Implemnteriung so das ein groﬂer chunk im speicher inialitsiert ist und nicht random
//alle paramter liegen beianander 

#include<iostream>
#include<SimulationRessources.h>
#include <stdexcept>
#include<utility>
#include<cmath>
//#include <omp.h>
#include<array>
#include <cuda_runtime.h>
#include <cuda.h>
//#include <stdexcept>
template<class T>
class Matrix3D
{
public:

	/*size_t index(int x, int y, int z)
	{
		const return(data + x * m_height * m_depth + y * m_depth + z);
	};*/
	__host__ __device__
	struct indices
	{
		int x, y, z;
	
		__host__ __device__
		indices(int _x, int _y, int _z) : x(_x), y(_y), z(_z) {}
		__host__ __device__
		indices(int ind[3]) : x(ind[0]), y(ind[1]), z(ind[2]) {}
		__host__ __device__
		indices(std::array<int,3> ind) : x(ind[0]), y(ind[1]), z(ind[2]) {}
	};

	__host__ __device__ 
	Matrix3D();
	__host__ __device__
	Matrix3D(int m, int n, int o, T* _data);
	__host__ __device__
	Matrix3D(int m, int n, int o, bool angles = false);
	__host__ __device__
	Matrix3D(const Matrix3D& m);
	__host__ __device__
	~Matrix3D();

	__host__ __device__
	float getAngle(indices _idx);
	__host__ __device__
	T& operator[](indices _idx);
	//void set(unsigned _x, unsigned _x,  T _v);
	//T get(size_t _i);
	__host__ __device__
	void setAll2(T _value);
	__host__ __device__
	int* getDimensions();
	__host__ __device__
	void setAmbientValue(T _d);
	__host__ __device__
	T getAmbientValue();
	__host__ __device__
	void setSpacing(double* _sp);
	__host__ __device__
	double* getSpacing();
	
	__device__ __host__
	void setValue(indices _idx, T _value);
	__device__ __host__
	void setValue(int _idx, T _value);
	//void setBoundaryCondition(bool _bound)
	T* getData() { return data; };
private: 

	int m_width;
	int m_height;
	int m_depth;
	T* data;
	float* angles;
	T defaultValue;
	double spacing[3];
	void setCylinderCoordinates();
	bool deleteData;
};


template<class T>
Matrix3D<T>::Matrix3D()
{
	m_width = m_height = m_depth = 1;
	defaultValue = T();
	spacing[0] = 0;
	data = new T[m_width * m_height * m_depth];
	angles = new float[m_width * m_height * m_depth];
}
template<class T>
Matrix3D<T>::Matrix3D(int m, int n, int o, T* _data)
{
	m_width = m;
	m_height = n;
	m_depth = o;
	defaultValue = T();
	spacing[0] = 0;
	data = _data;
	angles = nullptr;
	deleteData = false;
	
}
template<class T>
Matrix3D<T>::Matrix3D(int m, int n,  int o, bool _angles)
{
	m_width = m;
	m_height = n;
	m_depth = o;
	defaultValue = T();
	spacing[0] = 0;
	data = new T[m_width * m_height * m_depth];
	angles = new float[m_width * m_height * m_depth];
	if(_angles)
		setCylinderCoordinates();
}
template<class T>
Matrix3D<T>::Matrix3D(const Matrix3D& m)
{
	m_width = m.m_width;
	m_height = m.m_height;
	m_depth = m.m_depth;
	defaultValue = m.defaultValue;
	for (int d = 0; d < 3; d++)
	{
		spacing[d] = m.spacing[d];
	}
	data = new T[m_width * m_height * m_depth];
	angles = new float[m_width * m_height * m_depth];
	//#pragma omp parallel for collapse(3)
	for (int i = 0; i < m_width; i++)
	{
		for (int j = 0; j < m_height; j++)
		{
			for (int k = 0; k < m_depth; k++)
			{
				*(data + i * m_height * m_depth + j * m_depth + k) = *(m.data + i * m_height * m_depth + j * m_depth + k);
				*(angles + i * m_height * m_depth + j * m_depth + k) = *(m.angles + i * m_height * m_depth + j * m_depth + k);
			}
		}
	}
}
template<class T>
Matrix3D<T>::~Matrix3D()
{
	//std::cout << "DeleteDAta-start" << std::endl;
	if(data != nullptr && deleteData)
		delete[] data;
	if(angles != nullptr)
		delete[] angles;
	//std::cout << "DeleteDAta-end" << std::endl;
}
template<class T>
float Matrix3D<T>::getAngle(indices _idx)
{
	 return angles[_idx.x * m_height * m_depth + _idx.y * m_depth + _idx.z];
}
template<class T>
void Matrix3D<T>::setCylinderCoordinates()
{
	//std::cout << "setCylinderCoordinates" << std::endl;
	double pi = 3.14159;


	
	std::pair<int, int> center = { GlobalVariables::pointOnNeedle_volume[0], GlobalVariables::pointOnNeedle_volume[2] };
	std::cout << m_width << " " <<" "<< m_height<<" <"<< m_depth << std::endl;
	std::pair<int, int> ref = { 0,1 };
//	#pragma omp parallel for collapse(3)
	for (int i = 0; i < m_width; i++)
	{
		for (int j = 0; j < m_height; j++)
		{
			for (int k = 0; k < m_depth; k++)
			{
				std::pair<int, int> vx = { i,k };
				std::pair<int, int> vx_dir = { vx.first - center.first, vx.second - center.second };
				
				
				*(angles + i * m_height * m_depth + j * m_depth + k) = std::atan2(vx_dir.first, vx_dir.second) * (180 / pi);
			}
		}
	}	
}

template<class T>
T& Matrix3D<T>::operator[](indices _idx)
{
	if (_idx.x < -1 || _idx.x > m_width ||
		_idx.y < -1 || _idx.y > m_height ||
		_idx.z < -1 || _idx.z > m_depth)
		printf("BOUNDARY: %d , %d , %d \n", _idx.x, _idx.y, _idx.z);
		//throw std::invalid_argument("OUt of Boundary");
	//adibatisch - von neumann 
	if (false)//GlobalVariables::boundaryCondition == 0
	{
		int x = _idx.x;
		int y = _idx.y;
		int z = _idx.z;

		if (_idx.x < 0)
			x = 0;	
		if (_idx.x > (m_width - 1))
			x = (m_width - 1);

		if (_idx.y < 0)
			y = 0;
		if (_idx.y > (m_height - 1))
			y = (m_height - 1);
	
		if (_idx.z < 0)
			z = 0;	
		if (_idx.z > (m_depth - 1))
			z = (m_depth - 1);

		/*if (_idx.x == -1)
			x = 0;	
		else if (_idx.x == (m_width))
			x = (m_width - 1);
		else
			throw std::out_of_range("raus");


		if (_idx.y < 0)
			y = 0;
		if (_idx.y > (m_height - 1))
			y = (m_height - 1);
	
		if (_idx.z < 0)
			z = 0;	
		if (_idx.z > (m_depth - 1))
			z = (m_depth - 1);*/

		return data[x * m_height * m_depth + y * m_depth + z];
	}
	else
	{
		if (_idx.x < 0 || _idx.x >(m_width - 1)
			|| _idx.y < 0 || _idx.y >(m_height - 1)
			|| _idx.z < 0 || _idx.z >(m_depth - 1))
		{
			//std::cout << _idx.x << " " << _idx.y << " " << _idx.z << std::endl;
			return defaultValue;
		}

		
		return data[_idx.x * m_height * m_depth + _idx.y * m_depth + _idx.z];
	}

}

template<class T>
void Matrix3D<T>::setAll2(T _value)
{
//	#pragma omp parallel for collapse(3)
	for (int i = 0; i < m_width; i++)
	{
		for (int j = 0; j < m_height; j++)
		{
			for (int k = 0; k < m_depth; k++)
			{
				*(data + i * m_height * m_depth + j * m_depth + k) = _value;
			}
		}
	}
}

template<class T>
int* Matrix3D<T>::getDimensions()
{
	static int dim []= { m_width, m_height, m_depth };
	//std::cout << "correct" << m_width << m_height << m_depth << std::endl;
	return dim;
}
template<class T>
void Matrix3D<T>::setAmbientValue(T _d)
{
	defaultValue = _d;
}
template<class T>
T Matrix3D<T>::getAmbientValue()
{
	return defaultValue;
}


template<class T>
void Matrix3D<T>::setSpacing(double* _sp)
{
	for (int i = 0; i < 3; i++)
	{
		spacing[i] = _sp[i];
	}
}

template<class T>
double * Matrix3D<T>::getSpacing()
{
	return spacing;
}
template<class T>
void Matrix3D<T>::setValue(indices _idx, T _value)
{
	*(data + _idx.x * m_height * m_depth + _idx.y * m_depth + _idx.z) = _value;
}
template<class T>
void Matrix3D<T>::setValue(int _idx, T _value)
{
	*(data + _idx) = _value;
}
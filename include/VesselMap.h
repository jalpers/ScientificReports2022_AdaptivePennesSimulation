#pragma once
#ifndef VESSELMAP_H
#define VESSELMAP_H
#include<Matrix3D.h>
#include<CoordinatesConverter.h>


#include<vtkImageData.h>
#include<vtkSmartPointer.h>
#include<vtkDICOMReader.h>
#include<vtkIntArray.h>
#include<vtkImageExtractComponents.h>
#include<vtkNIFTIImageReader.h>
class VesselMap
{
public:
	VesselMap();
	VesselMap(CoordinatesConverter _world2VolumeConverter, int * _size);
	~VesselMap();
	void load(std::string _str);
	//void getMap(Matrix3D<bool>& _map) {
	//	_map = Matrix3D<bool>(*map);
	//}
	Matrix3D<bool>* getMap()
	{
		return map;
	}
private:
	Matrix3D<bool>* map;
	CoordinatesConverter world2VolumeConverter;
	int volumeSize[3];
};




#endif // !VESSELMAP_H

#pragma once

#ifndef DATAVOLUME_H
#define DATAVOLUME_H

#include<vtkImageData.h>
#include<CoordinatesConverter.h>
#include<DicomHandler.h>
#include<SimulationRessources.h>
#include<Pathfinder.h>
//vtk
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkRenderer.h>
#include<vtkNamedColors.h>
//#include <vtkOpenGLGPUVolumeRayCastMapper.h>
#include <vtkVolumeProperty.h>

#include <vtkLookupTable.h>
#include <vtkPolyData.h>

#include <vtkContourValues.h>
#include<qdebug.h>
#include<vector>
#include <deque>
#include<array>

#include <vtkImageMathematics.h>

#include<vtkOBJExporter.h>
#include<vtkDataSetMapper.h>
#include<vtkXMLImageDataWriter.h>
#include<vtkImageDataGeometryFilter.h>
#include<vtkPolyDataWriter.h>
#include<vtkImageDataGeometryFilter.h>
#include<vtkDICOMWriter.h>
#include <vtkDICOMMRGenerator.h>
#include <vtkImageData.h>
#include<vtkMatrix3x3.h>
#include<vtkImageResize.h>
#include<vtkImageMapper.h>
#include<vtkImageHybridMedian2D.h>

#include <vtkFixedPointVolumeRayCastMapper.h>
#define VOLUMESIZE 60
class DataVolume
{
public:

	DataVolume();
	/*
	* Constructor: 
	* Initialisates vtkSmartPointer<vtkImageData> volume 
	* calculates CoordinatesConverter world2VolumeConverter
	* @param _filename: directory of dicom files
	*/
	DataVolume(QString _filenName);

	DataVolume(vtkSmartPointer<vtkImageData> _imageData);
	~DataVolume();


	/*
	* writes the created volume to a stack of dicom files
	* @_targetfileName directory where files are saved
	*/
	void writeToFile(QString _targetfileName);

	/*
	* load and add 2D data to vtkSmartPointer<vtkImageData> volume
	* 
	* 
	*/
	void addSlice(QString _angle = "0", QString _timeStep = "0");
	
	
	void display(vtkSmartPointer<vtkRenderWindow> _q_vtk_renderWindow, bool _isCropped = true, bool _cameraAbove = false);

	
	
	
	vtkSmartPointer<vtkImageData> getImageData();
	
	
	

	void crop(vtkSmartPointer<vtkImageData> volume, vtkSmartPointer<vtkImageData> croppedvolume, float percentage);
	void crop(float percentage);

	CoordinatesConverter getWorld2VoxelConverter() { return world2VolumeConverter; };
	vtkSmartPointer<vtkImageData> getcurrentSlice();
	std::vector<double> *currentMeasValues() { return &(currentSlice->measpoints); };
	std::vector<std::array<int, 3>>* getcurrentMeasPos() { return &(currentSlice->measPointsInVolume); };
	std::vector<std::array<int, 3>>* getcurrentSlicePos() { return &(currentSlice->slicePositionInVolume); };

	std::vector<double>* getMeasValues() { 
		int nrOfMeasVal = 0;
		measpoints.clear();
		for (int i = 0; i < GlobalVariables::numberOfSlices; i++)
		{
			nrOfMeasVal += slicesBuffer[i].measpoints.size();
		}
		measpoints.reserve(nrOfMeasVal);
		for (int i = 0; i < GlobalVariables::numberOfSlices; i++)
		{
			measpoints.insert(measpoints.end(), slicesBuffer[i].measpoints.begin(), slicesBuffer[i].measpoints.end());
		}
		return &measpoints;
	
	};
	std::vector<std::array<int, 3>>* getMeasPos() { 
		int nrOfMeasPos = 0;
		measPointsInVolume.clear();
		for (int i = 0; i < GlobalVariables::numberOfSlices; i++)
		{
			nrOfMeasPos += slicesBuffer[i].measPointsInVolume.size();
		}
		measpoints.reserve(nrOfMeasPos);
		for (int i = 0; i < GlobalVariables::numberOfSlices; i++)
		{
			measPointsInVolume.insert(measPointsInVolume.end(), slicesBuffer[i].measPointsInVolume.begin(), slicesBuffer[i].measPointsInVolume.end());
		}
		return &measPointsInVolume;
	};

	static unsigned char reverse(unsigned char b) {
		b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
		b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
		b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
		return b;
	};
	static int decode_vlq(unsigned char* input)
	{
		int result = 0;
		do
		{
			result = (result << 7) | (*input & 0x7F);
		} while (*input++ & 0x80);
		return result;
	}
	static int encode(unsigned char a)
	{
		return a; //>> 1;
	}
	vector<int> getTempAlongNeedle()
	{
		return tempAlongNeedle;
	}
	vector<double> getQAlongNeedle()
	{
		return q_relAlongNeedle;
	}
	void defaultTempAlongNeedle();
	std::vector<std::array<int, 3>>* posNeedleInVolume()
	{
		return &needleInVolume;
	}
	
private:
	std::vector<std::array<int, 3>> needleInVolume;
	void init(QString _fileName);
	vtkSmartPointer<vtkImageData> volume;
	CoordinatesConverter world2VolumeConverter;
	int volumeSize = VOLUMESIZE;
	QString fileName;
	
	struct Slice
	{
		vtkSmartPointer<vtkImageData> data;
		std::vector<double> measpoints;
		CoordinatesConverter voxel2WorldConverter;
		int timestep = -1;
		std::vector<std::array<int,3>> slicePositionInVolume;
		std::vector<std::array<int,3>> measPointsInVolume;
		std::vector<int> needleAxisInSlice;
		int mainNeedleDirSlice;
		vector<double> distance_avag;
		int maxHeat;
	
		Slice() : data(vtkSmartPointer<vtkImageData>::New()), measpoints(1, 0),
			voxel2WorldConverter(CoordinatesConverter()), timestep(-1), slicePositionInVolume(1, std::array<int, 3>({ 0,0,0 })), measPointsInVolume(1, std::array<int, 3>({ 0,0,0 })), needleAxisInSlice(0), mainNeedleDirSlice(0), distance_avag(0), maxHeat(0) {}
		Slice(vtkSmartPointer<vtkImageData> _data,
			std::vector<double> _measpoints,
			CoordinatesConverter _voxel2WorldConverter,
			int _timestep ,
			std::vector<std::array<int, 3>> _slicePositionInVolume,
			std::vector<std::array<int, 3>>  _measPointsInVolume,
			std::vector<int> _needleAxisInSlice,
			int _mainNeedleDirSlice,
			vector<double> _distance_avag,
			int _maxHeat):
			data(_data),
			measpoints(_measpoints),
			voxel2WorldConverter(_voxel2WorldConverter),
			timestep(_timestep),
			slicePositionInVolume(_slicePositionInVolume),
			measPointsInVolume(_measPointsInVolume),
			needleAxisInSlice(_needleAxisInSlice),
			distance_avag(_distance_avag),
			mainNeedleDirSlice(_mainNeedleDirSlice),
			maxHeat(_maxHeat) {}
			

	};
	Slice slicesBuffer[GlobalVariables::numberOfSlices];
    Slice* currentSlice;
	
	std::vector<double> measpoints;
	std::vector<std::array<int, 3>> measPointsInVolume;

	void createCoordNeedleAxis(DicomHandler::dicomDataProperties _p0, DicomHandler::dicomDataProperties _p90);
	//std::vector<std::array<double,3>()> coordNeedleAxis;
	void crossProduct(float v_A[], float v_B[], float c_P[]);
	float dotProduct(float v_A[], float v_B[]);


	void updateTempNeedle();
	vector<int> tempAlongNeedle;
	vector<double> q_relAlongNeedle;
	int timestep;

	float l2Norm(float* a, float* b, int N);
	void insertionSort(int arr[], int n)
	{
		int i, key, j;
		for (i = 1; i < n; i++)
		{
			key = arr[i];
			j = i - 1;

			/* Move elements of arr[0..i-1], that are
			greater than key, to one position ahead
			of their current position */
			while (j >= 0 && arr[j] > key)
			{
				arr[j + 1] = arr[j];
				j = j - 1;
			}
			arr[j + 1] = key;
		}
	}
	
	

	

	
	

	

	
	
	QStringList slices_angle = { "0","22_5","45","67_5","90","112_5","135","157_5" };
	//QStringList slices_time = { "0","90","45","135_5","22_5","112_5","67_5","157_5" };
	QStringList slices_time = { "0","90","45","135_5","22_5","112_5","67_5","157_5" };
	//egal ob volumeSize ungerade oder gerade -> bei ungerade wird Layer 0 nicht besetzt
	//Alle gesetzten Punkte auf einem Layer

	DicomHandler* dicomHandler;
	int mainneedledir;
	int maindirVol;

	void setAllVoxel2Value(int _value);

	
	
};

#endif // DATAVOLUME_H
#include "DataVolume.h"
#include <math.h> 
#include<algorithm>
#include<string>
#include<Debugger.h>
#include<vector>
//#pragma optimize( "", off )
#define mytype uchar
DataVolume::DataVolume()
{
	volume = vtkSmartPointer<vtkImageData>::New();
	
	dicomHandler = new DicomHandler();
	//slicesBuffer = new Slice[6];
	//qDebug() << sizeof(slicesBuffer);
	//volumeSize = VOLUMESIZE;
}

DataVolume::DataVolume(QString _fileName)
{
	
	DataVolume();
	timestep = 0;
	currentSlice = new Slice();
	init(_fileName);
}

DataVolume::DataVolume(vtkSmartPointer<vtkImageData> _imageData)
{
	volume = _imageData;
	
}






DataVolume::~DataVolume()
{
	
}



void DataVolume::init(QString _fileName)
{

	qDebug() << "DataVolume::CreateVolume";
	fileName = _fileName;


	

	//get properties of orthogonal images to define volume.
	QString filename0 = _fileName + "/0/" + "0.dcm";
	QString filename90 = _fileName + "/90/" + "0.dcm";
	//qDebug() << _fileName + "/0.dcm";
	//QString filename0 = _fileName + "/0.dcm";
	//QString filename90 = _fileName +  "/90.dcm";
	//QString filename0 = _fileName + "/0_14.IMA";
	//QString filename90 = _fileName +  "/90_14.IMA";

	vtkSmartPointer<vtkDICOMReader> reader = vtkSmartPointer<vtkDICOMReader>::New();
	
	reader->SetFileName(filename0.toStdString().c_str());
	reader->Update();

	DicomHandler dHandler;
	DicomHandler::dicomDataProperties properties0;
	dHandler.getImageDataProperties(reader, &properties0);
	
	CoordinatesConverter converter(properties0);

	reader->SetFileName(filename90.toStdString().c_str());
	reader->Update();
	DicomHandler::dicomDataProperties properties90;
	dHandler.getImageDataProperties(reader, &properties90);

	volumeSize = properties90.dimension[0];

	CoordinatesConverter::CoordinateSystemProperties worldCoordSystem;
	//-----------------------------------------------------------------------------------------------------
	std::cout << "STAAAAAAAAArt" << std::endl;
	std::vector<float> coords[3];

	for (int i = 0; i < 8; i++)
	{
		QString filename = _fileName + "/"+ slices_angle[i] +"/"	+ "0"	 + ".dcm";
		qDebug() << filename;
		vtkSmartPointer<vtkDICOMReader> reader = vtkSmartPointer<vtkDICOMReader>::New();

		reader->SetFileName(filename.toStdString().c_str());
		reader->Update();

		DicomHandler dHandler;
		DicomHandler::dicomDataProperties properties;
		dHandler.getImageDataProperties(reader, &properties);



		CoordinatesConverter converter(properties);
		float* worldCoord1 = new float[3];
		float* worldCoord2 = new float[3];

		//get values of corners _properties90.
		worldCoord1 = converter.transform(0, 0);
		worldCoord2 = converter.transform(properties.dimension[0] - 1, properties.dimension[1] - 1);

		for(int i = 0; i < 3; i++)
		{
			coords[i].push_back(worldCoord1[i]);
			coords[i].push_back(worldCoord2[i]);
		}

		delete[] worldCoord1;
		delete[] worldCoord2;


	}
	for (int i = 0; i < 3; i++)
	{
		worldCoordSystem.maxExtend[i] = *std::max_element(coords[i].begin(), coords[i].end()) + properties0.pixelSpacing[1];
		worldCoordSystem.minExtend[i] = *std::min_element(coords[i].begin(), coords[i].end()) - properties0.pixelSpacing[1];
		worldCoordSystem.origin[i] = worldCoordSystem.minExtend[i];
		qDebug() << "mine:	" << worldCoordSystem.minExtend[i];
		qDebug() << "maxe:	" << worldCoordSystem.maxExtend[i];
		qDebug() << "origincal:	" << worldCoordSystem.origin[i];

	}

	//------------------------------------------------------------------------------------------

	
	//get extend of World Coordinate System
	//converter.getCoordinateSystemProperties(properties90, &worldCoordSystem);

	CoordinatesConverter::CoordinateSystemProperties volumeCoordSystem;
	int vs[3];
	std::cout << "asdfasdfdsf";
	std::cout << "GALLLLO " << properties0.pixelSpacing[1]<<std::endl;
	for (int i = 0; i < 3; i++)
	{
		vs[i] = std::ceil((1 / properties0.pixelSpacing[1]) * (worldCoordSystem.maxExtend[i] - worldCoordSystem.minExtend[i]))+1;//
		std::cout << "vs " << vs[i] << std::endl;
		volumeCoordSystem.maxExtend[i] = vs[i]-1;//volumesize- 1
		volumeCoordSystem.minExtend[i] = 0;
		volumeCoordSystem.origin[i] = 0;
	}
	std::cout << "geht" << std::endl;
	volume = vtkSmartPointer<vtkImageData>::New();

	//volume->SetDimensions(volumeSize, volumeSize, volumeSize);
	volume->SetDimensions(vs[0], vs[1], vs[2]);
	volume->SetExtent(0, vs[0] - 1, 0, vs[1] - 1, 0, vs[2] - 1);
	
		//
	volume->AllocateScalars(VTK_UNSIGNED_CHAR, 1);//
	volume->SetSpacing(properties0.pixelSpacing[0], properties0.pixelSpacing[1], properties0.pixelSpacing[0]);
	volume->Modified();

	CoordinatesConverter::TransformationMatrixFeatures world2VolumeTransform;
	converter.getTransMatrixFeatures(volumeCoordSystem, worldCoordSystem, &world2VolumeTransform);


	


	world2VolumeConverter = CoordinatesConverter(world2VolumeTransform);
	volume->SetOrigin(-world2VolumeConverter.getTransformMatrix()->GetElement(0, 3), world2VolumeConverter.getTransformMatrix()->GetElement(1, 3), world2VolumeConverter.getTransformMatrix()->GetElement(2, 3));
	volume->Modified();

	qDebug ()<< "_________________________________________________________________________________________";
	qDebug() << "DataVolume: ";
	world2VolumeConverter.getTransformMatrix()->Print(std::cout);
	std::cout << "Dimension: " << volume->GetDimensions()[0] << " " << volume->GetDimensions()[1] << " " << volume->GetDimensions()[2] << std::endl;
	std::cout << "Spacing: " << volume->GetSpacing()[0] << " " << volume->GetSpacing()[1] << " " << volume->GetSpacing()[2] << std::endl;
	qDebug() << "_________________________________________________________________________________________";

	setAllVoxel2Value(0);
//tEst################################################################
	//vtkSmartPointer<vtkImageData> a = dHandler.loadDicom(_fileName + "/" + "0" + "/0.IMA");
	//vtkSmartPointer<vtkImageData> b = dHandler.loadDicom(_fileName + "/" + "90" + "/0.IMA");

	//double* worldCoord1 = new double[3];
	//double* worldCoord2 = new double[3];
	//a->TransformContinuousIndexToPhysicalPoint(0, 0, 0, worldCoord1);
	//a->TransformContinuousIndexToPhysicalPoint(59, 59, 0, worldCoord2);


	//double* worldCoord3 = new double[3];
	//double* worldCoord4 = new double[3];


	//for (int i = 0; i < 3; i++)
	//{
	//}

	//double direction[9] = { 0 };
	//direction[0] = properties90.imageOrientationX[0];
	//direction[1] = properties90.imageOrientationX[1];
	//direction[2] = properties90.imageOrientationX[2];
	//direction[3] = properties90.imageOrientationY[0];
	//direction[4] = properties90.imageOrientationY[1];
	//direction[5] = properties90.imageOrientationY[2];
	//double result[16] = { 0 };
	//b->ComputeIndexToPhysicalMatrix((double*)properties90.imagePosition, (double*)properties90.pixelSpacing, direction, result);
	//std::cout << "result" << std::endl;
	//for (int i = 0; i < 16; i++)
	//	std::cout << result[i] << "" << std::endl;
	//b->SetOrigin((double*)properties90.imagePosition);
	//b->SetSpacing(properties90.pixelSpacing);
	//b->SetDirectionMatrix(direction);
	//b->Modified();

	//
	//b->TransformIndexToPhysicalPoint(0, 0, 0, worldCoord3);
	//b->TransformIndexToPhysicalPoint(59, 59, 0, worldCoord4);
	//for (int i = 0; i < 3; i++)
	//{
	//	qDebug() << "worlfCoord_new: " << worldCoord1[i] << worldCoord2[i] << worldCoord3[i] << worldCoord4[i];
	//}

	//b->GetIndexToPhysicalMatrix()->Print(std::cout);



	//Berechnung der Nadelachse
	DicomHandler::dicomDataProperties _ref;
	dHandler.getImageDataProperties(_fileName + "/" + "0" + "/0.dcm",&_ref);


	double wC[3];

	float* vC = new float[3];

	for (int i = 0; i < slices_angle.size(); i++)
	{

		DicomHandler::dicomDataProperties _in;
		dHandler.getImageDataProperties(_fileName + "/" + slices_angle.at(i) + "/0.dcm", &_in);
		createCoordNeedleAxis(_ref,_in);

		for (int j = 0; j < 3; j++)
			wC[j] = GlobalVariables::pointOnNeedle_world[j] + (30)*GlobalVariables::needleDir[j] * volume->GetSpacing()[j];

		vC = world2VolumeConverter.transformTranslationFirst(wC[0], wC[1], wC[2]);
		GlobalVariables::pointOnNeedle_volume[0] = (int)round(vC[0]);
		GlobalVariables::pointOnNeedle_volume[1] = (int)round(vC[1]);
		GlobalVariables::pointOnNeedle_volume[2] = (int)round(vC[2]);


		std::cout << "GlobalVariables::pointOnNeedle_volume: " << GlobalVariables::pointOnNeedle_volume[0] << " " << GlobalVariables::pointOnNeedle_volume[1] << " " << GlobalVariables::pointOnNeedle_volume[2] << std::endl;
		std::cout << "GlobalVariables::pointOnNeedle_world: " << GlobalVariables::pointOnNeedle_world[0] << " " << GlobalVariables::pointOnNeedle_world[1] << " " << GlobalVariables::pointOnNeedle_world[2] << std::endl;
		std::cout << "GlobalVariables::needleDir[0]: " << GlobalVariables::needleDir[0] << " " << GlobalVariables::needleDir[1] << " " << GlobalVariables::needleDir[2] << std::endl;
		
	}
	//volume
	
	float ptr[3];
	ptr[0] = abs((float)GlobalVariables::needleDir[0]);
	ptr[1] = abs((float)GlobalVariables::needleDir[1]);
	ptr[2] = abs((float)GlobalVariables::needleDir[2]);
	maindirVol = std::distance(ptr, std::max_element(ptr, ptr + 3));
	std::cout << "maindirVol " << maindirVol << std::endl;
	for (int i = 0; i < volume->GetDimensions()[maindirVol]; i++)
	{
		float z = (i - GlobalVariables::pointOnNeedle_volume[maindirVol]) / (GlobalVariables::needleDir[maindirVol]);
		//std::array<int, 2> h = { i, (int)needlePosSlice[1 - currentSlice->mainNeedleDirSlice] + z * inSliceDirection[1 - currentSlice->mainNeedleDirSlice] };
		//std::cout << h[0] << " " << h[1] << std::endl;
		int h[3];
		for (int j = 0; j < 3; j++)
		{
			h[j] = GlobalVariables::pointOnNeedle_volume[j] + z * (GlobalVariables::needleDir[j]);
		}
		std::cout << h[0] << " " << h[1] << " "<<h[2] << std::endl;

		needleInVolume.push_back(std::array<int, 3>({ h[0],h[1],h[2] }));
	}

	int sz = volume->GetDimensions()[maindirVol];
	std::cout << "size " << sz<< std::endl;
	//entlang der Nadelachse
	//tempAlongNeedle.resize(vs[maindirVol]);
	for (int i = 0; i < sz; i++)
	{
		//std::cout << i << std::endl;
		tempAlongNeedle.push_back(GlobalVariables::baseLineT);
	}
		

	//q_relAlongNeedle.resize(vs[maindirVol]);
	for (int i = 0; i < sz; i++)
	{
		//std::cout << i << std::endl;
		q_relAlongNeedle.push_back(0.0);
	}



	//Slice
	for (int i = 0; i < slices_angle.size(); i++)
	{
		DicomHandler::dicomDataProperties dicomHeader;
		dHandler.getImageDataProperties(_fileName + "/" + slices_angle.at(i) + "/0.dcm", &dicomHeader);
		CoordinatesConverter world2VoxelConverter(dicomHeader);
		world2VoxelConverter.invert();

		float* needlePosSlice;
		needlePosSlice = world2VoxelConverter.transform(GlobalVariables::pointOnNeedle_world[0] + 2 * GlobalVariables::needleDir[0], GlobalVariables::pointOnNeedle_world[1] + 2 * GlobalVariables::needleDir[1], GlobalVariables::pointOnNeedle_world[2] + 2 * GlobalVariables::needleDir[2]);
		//std::cout << "Hallo " << needlePosSlice[0] << " " << needlePosSlice[1] << " " << needlePosSlice[2] << std::endl;


		float* needlePosSlice2;
		needlePosSlice2 = world2VoxelConverter.transform(GlobalVariables::pointOnNeedle_world[0], GlobalVariables::pointOnNeedle_world[1], GlobalVariables::pointOnNeedle_world[2]);
		//std::cout << "Hallo " << needlePosSlice2[0] << " " << needlePosSlice2[1] << " " << needlePosSlice2[2] << std::endl;

		float inSliceDirection[2];
		for (int j = 0; j < 2; j++)
		{
			inSliceDirection[j] = needlePosSlice[j] - needlePosSlice2[j];
		}
		//std::cout << "Hallo " << inSliceDirection[0] << " " << inSliceDirection[1] << std::endl;

		currentSlice = &slicesBuffer[i];
		currentSlice->mainNeedleDirSlice = inSliceDirection[0] > inSliceDirection[1] ? 0 : 1;
		//std::cout << "Mainneedle " << currentSlice->mainNeedleDirSlice << "dim " << dicomHeader.dimension[0] << " " << dicomHeader.dimension[1] << std::endl;
		for (int k = 0; k < dicomHeader.dimension[currentSlice->mainNeedleDirSlice]; k++)
		{
			float z = (k - needlePosSlice[currentSlice->mainNeedleDirSlice]) / inSliceDirection[currentSlice->mainNeedleDirSlice];
			std::array<int,2> h = { k, (int)needlePosSlice[1 - currentSlice->mainNeedleDirSlice] + z * inSliceDirection[1 - currentSlice->mainNeedleDirSlice] };
			//std::cout << h[0] << " " << h[1] << std::endl;
			currentSlice->needleAxisInSlice.push_back(h[1]);
			
		}
		delete[] needlePosSlice;
		delete[] needlePosSlice2;
	}






	

	delete[] vC;
	qDebug() << "DataVolume::CreateVolume - Done";




}

void DataVolume::writeToFile(QString _targetfileName)
{
	DicomHandler dhandler;

	vtkSmartPointer<vtkDICOMWriter> writer_d =
		vtkSmartPointer<vtkDICOMWriter>::New();
	
	// Create a generator for MR images.
	vtkNew<vtkDICOMMRGenerator> generator;
	
	vtkSmartPointer<vtkMatrix4x4> p = vtkSmartPointer<vtkMatrix4x4>::New();
	p->GlobalWarningDisplayOff();
	p->Identity();
	p->Modified();

	generator->Modified();
	writer_d->SetGenerator(generator.GetPointer());

	// Set the output filename format as a printf-style string.
	writer_d->SetFilePattern("%s/IM-0001-%04.3d.dcm");
	
	std::string l = _targetfileName.toStdString() +"/output";
	writer_d->SetFilePrefix(l.c_str());

	writer_d->SetInputData(volume);

	writer_d->Modified();
	writer_d->Write();

	qDebug() << "write Dicom -Done";


}

void DataVolume::addSlice(QString _angle, QString _timeStep)
{
	qDebug() << "DataVolume::addSlice";

	int index_slice = slices_angle.lastIndexOf(_angle);
	currentSlice = &slicesBuffer[index_slice];
	QString fileDataName = fileName + "/" + _angle + "/"  + _timeStep + ".dcm";
	qDebug() << fileDataName;
	//QString fileDataName = fileName + "/" + _angle  +".dcm";//.dcm _14.IMA
	currentSlice->data = dicomHandler->loadDicom(fileDataName);
	
	DicomHandler::dicomDataProperties dicomHeader;
	dicomHandler->getImageDataProperties(fileDataName, &dicomHeader);

	CoordinatesConverter voxel2WorldConverter(dicomHeader);
	


	timestep++;
	currentSlice->timestep = timestep;
	currentSlice->voxel2WorldConverter = voxel2WorldConverter;
	
	float* voxelVolumePosition = new float[3];
	float* voxelWorldPosition = new float[3];

	//double* voxelVolumePositio = new double[3];
	int* dim = dicomHeader.dimension;
	qDebug() << dim[0]<< " "<< dim[1];
	int mD = currentSlice->mainNeedleDirSlice;
	std::cout << "md " << mD << std::endl;
	int dimNA = currentSlice->needleAxisInSlice.size();
	bool isNeedleYAxisAligned = mD;

	int index = 0;
	int Heatsum = 0;

	const int sizeWindow = 5;
	
	const int nrIsoL = 8;
	int isolinesValues[nrIsoL] = { 28,31,34,37 , 40, 43, 46, 49 };
//,45,50



	int N_mp = dimNA * nrIsoL *2;
	//int N_mp = dim[0] * dim[1];
	
	//std::vector<std::array<int,3>> positionInVolume;
	currentSlice->slicePositionInVolume.clear();
	currentSlice->slicePositionInVolume.reserve(dim[0]* dim[1]);

	currentSlice->measPointsInVolume.clear();
	currentSlice->measPointsInVolume.reserve(N_mp);

	currentSlice->measpoints.clear();
	currentSlice->measpoints.reserve(N_mp);
	
	std::cout << "Hallo" << std::endl;



	vector<pair<int, int>> border_left;
	vector<pair<int, int>> border_right;

	border_left.resize(dim[0]);
	border_right.resize(dim[0]);

	std::vector<std::vector<int>> grid_left[nrIsoL];
	std::vector<std::vector<int>> grid_right[nrIsoL];
	std::vector<std::vector<int>> grid_iso;
	std::vector<std::vector<int>> grid[nrIsoL];
	std::vector<std::vector<double>> thermomap;

	for (int i = 0; i < nrIsoL; i++)
	{
		grid_left[i].resize(dim[0]);
		grid_right[i].resize(dim[0]);
		grid[i].resize(dim[0]);
		grid_iso.resize(dim[0]);
		thermomap.resize(dim[0]);

		for (size_t r = 0; r < dim[0]; r++)
		{
			grid[i][r].resize(dim[1]);
			grid_left[i][r].resize(dim[1]);
			grid_right[i][r].resize(dim[1]);
			grid_iso[r].resize(dim[1]);
			thermomap[r].resize(dim[1]);
		}
	}

	Pathfinder::Pair start;
	Pathfinder::Pair end;
	Pathfinder::Pair back;
	std::cout << "Hallo" << std::endl;
	//bool isNeedleYAxisAligned = voxel2WorldConverter.isNeedleYAxisAligned();



	
	Pathfinder pf;

	std::cout << "Hallo " << std::endl;
	/*for (int i = 0; i < dim[0]; i++)
	{
		
		if (isNeedleYAxisAligned)
		{
			int y_ = dim[1] - i - 1;
			border_left[i] = std::pair<int, int>(GlobalVariables::pointOnNeedle_volume[0] - 1, y_ );
			border_right[i] = std::pair<int, int>(GlobalVariables::pointOnNeedle_volume[0] + 1, y_ );

		}
		else
		{
			border_left[i] = std::pair<int, int>(i , dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] - 1 );
			border_right[i] = std::pair<int,int>(i , dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] + 1 );
		}
	}*/

	//if (isNeedleYAxisAligned)
	//{
	//	grid_left[i][GlobalVariables::pointOnNeedle_volume[0] - 1][y_] = UINT8_MAX;
	//	grid_right[i][GlobalVariables::pointOnNeedle_volume[0] + 1][y_] = UINT8_MAX;

	//}
	//else
	//{
	//	grid_left[i][x][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] - 1] = UINT8_MAX;
	//	grid_right[i][x][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] + 1] = UINT8_MAX;
	//}
	//std::cout << "Hall1o" << volume->GetDimensions()[0]<< " "<< volume->GetDimensions()[1] <<" "<< volume->GetDimensions()[2] <<std::endl;
	for (unsigned int x = 0; x < dim[0] ; x++)
	{
		
		for (unsigned int y = 0; y < dim[1]; y++)
		{
			int y_ = dim[1] - y - 1;
			voxelWorldPosition = voxel2WorldConverter.transform(x, y_);
			voxelVolumePosition = world2VolumeConverter.transformTranslationFirst(voxelWorldPosition[0], voxelWorldPosition[1], voxelWorldPosition[2]);
			//std::cout <<"voxVolalt: "<< voxelVolumePosition[0] << " " << voxelVolumePosition[1] << " " << voxelVolumePosition[2] << std::endl;
			//volume->TransformPhysicalPointToContinuousIndex((double)voxelWorldPosition[0], (double)voxelWorldPosition[1], (double)voxelWorldPosition[2],voxelVolumePositio);
			
			
			//std::cout << "voxVolstart: " << voxelVolumePosition[0] << " " << voxelVolumePosition[1] << " " << voxelVolumePosition[2] << std::endl;

			mytype* pixelImage = static_cast<mytype*>(currentSlice->data->GetScalarPointer(x, y, 0));
			if ((int)round(voxelVolumePosition[0]) > volume->GetDimensions()[0] - 1 ||

				(int)round(voxelVolumePosition[2]) > volume->GetDimensions()[2] - 1)
			{
				std::cout << "grö0er " << (voxelVolumePosition[0]) << " " << (voxelVolumePosition[1]) << " " << (voxelVolumePosition[2]) << std::endl;


			}
			if ((int)round(voxelVolumePosition[1]) > volume->GetDimensions()[1] - 1)
			{

			}




			if (round(voxelVolumePosition[0]) < 0)
			{
				std::cout << "kleiner " << (voxelVolumePosition[0]) << " " << (voxelVolumePosition[1]) << " " << (voxelVolumePosition[2]) << std::endl;
				voxelVolumePosition[2] = 0;
			}
			if (round(voxelVolumePosition[1]) < 0)
			{
				std::cout << "kleiner " << (voxelVolumePosition[0]) << " " << (voxelVolumePosition[1]) << " " << (voxelVolumePosition[2]) << std::endl;
				voxelVolumePosition[2] = 0;
			}
			if (round(voxelVolumePosition[2]) < 0)
			{
				std::cout << "kleiner "<<(voxelVolumePosition[0]) << " " << (voxelVolumePosition[1]) << " " << (voxelVolumePosition[2]) << std::endl;
				voxelVolumePosition[2] =0;
			}
			mytype* voxelVolume = static_cast<mytype*>(volume->GetScalarPointer(
				(int)round(voxelVolumePosition[0]), 
				(int)round(voxelVolumePosition[1]), 
				(int)round(voxelVolumePosition[2])));
		//	qDebug() << pixelImage[0];
			// 
			// 
			// 
			//TODO richtige Messpunkte aussuchen
			//alle punkte in der Slice 
			currentSlice->slicePositionInVolume[index][0] = (int)round(voxelVolumePosition[0]);
			currentSlice->slicePositionInVolume[index][1]= (int)round(voxelVolumePosition[1]);
			currentSlice->slicePositionInVolume[index][2] = (int)round(voxelVolumePosition[2]);
			////std::cout << (int)static_cast<unsigned char*>(volume->GetScalarPointer(0, x, y))[0];
			////if((int)round(voxelVolumePosition[1])==30)

			//alle Punkte
			/*currentSlice->measPointsInVolume[index][0] = (int)round(voxelVolumePosition[0]);
			currentSlice->measPointsInVolume[index][1] = (int)round(voxelVolumePosition[1]);
			currentSlice->measPointsInVolume[index][2] = (int)round(voxelVolumePosition[2]);
			currentSlice->measpoints[index] = voxelVolume[0];*/
			/*currentSlice->measPointsInVolume.push_back({ (int)round(voxelVolumePosition[0]) ,(int)round(voxelVolumePosition[1]),(int)round(voxelVolumePosition[2]) });

			currentSlice->measpoints.push_back(voxelVolume[0]);*/


			index++;



			voxelVolume[0] = encode(pixelImage[0]);
			//std::cout << (int)voxelVolume[0] << " ";
			Heatsum += voxelVolume[0];
			thermomap[x][y_] = abs(voxelVolume[0]);
			for (int i = 0; i < nrIsoL; i++)
			{
				
				grid[i][x][y_] = abs(voxelVolume[0] - isolinesValues[i]);
				grid_left[i][x][y_] = abs(voxelVolume[0] - isolinesValues[i]);
				grid_right[i][x][y_] = abs(voxelVolume[0] - isolinesValues[i]);



				/*if (isNeedleYAxisAligned)
				{
					grid_iso[GlobalVariables::pointOnNeedle_volume[0]][y_] = UINT8_MAX;
					grid_left[i][GlobalVariables::pointOnNeedle_volume[0] - 1][y_] = UINT8_MAX;
					grid_right[i][GlobalVariables::pointOnNeedle_volume[0] + 1][y_] = UINT8_MAX;
				}
				else
				{
					grid_iso[x][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2]] = UINT8_MAX;
					grid_left[i][x][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] - 1] = UINT8_MAX;
					grid_right[i][x][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] + 1] = UINT8_MAX;
				}*/


			}



		}
		//std::cout << std::endl;
	}
	int window[9];
	std::cout << "Hall2o" << std::endl;
	
	//start variable
	double* needlesumStore = new double[dimNA]();
	for (int d = sizeWindow; d < dimNA - sizeWindow; d++)
	{
		double sumneedle = 0;

		for (int w = -sizeWindow; w < sizeWindow; w++)
		{
			sumneedle = sumneedle + thermomap[currentSlice->needleAxisInSlice[d]][d + w];

		}
		needlesumStore[d] = sumneedle;
	}
	int maxnr = std::distance(needlesumStore, std::max_element(needlesumStore, needlesumStore + dimNA));
	delete[] needlesumStore;
	if (isNeedleYAxisAligned)
	{
		std::cout << "sNeedleYAxisAligned" << std::endl;
		start = { currentSlice->needleAxisInSlice[0], 0 };
		end = {   currentSlice->needleAxisInSlice[dimNA-1], dimNA-1 };
		grid_iso[start.first][start.second] = UINT8_MAX;
		//std::cout << "isNeedleYAxisAligned == true" << std::endl;
		//start pixel verbreitern
		for (int i = 0; i < nrIsoL; i++)
		{
			//median filter
			for (int x = 1; x < dim[0]-1; ++x)
			{
				for (int y = 1; y < dim[1]-1; ++y)
				{
					window[0] = grid[i][x - 1][y - 1];
					window[1] = grid[i][x + 1][y - 1];
					window[2] = grid[i][x][y - 1];
					window[3] = grid[i][x - 1][y];
					window[4] = grid[i][x + 1][y];
					window[5] = grid[i][x][y];
					window[6] = grid[i][x - 1][y + 1];
					window[7] = grid[i][x + 1][y + 1];
					window[8] = grid[i][x][y + 1];

					//sort window array
					insertionSort(window, 9);
					//put the median to the new array
					grid_left[i][x][y] = window[4];
					grid_right[i][x][y] = window[4];
					

				}
			}
			//Debugger::debugViewer(grid_left[i], "gridleft");
			/*double* grad = new double[dimNA - 1]();
			double* grad_acc = new double[dimNA - 1]();
			double sum = 0;
			for (int d = 0; d < dimNA - 1; d++)
			{
				grad[d] = (double)abs(grid_left[i][currentSlice->needleAxisInSlice[d]][d] - grid_left[i][currentSlice->needleAxisInSlice[d]][d + 1]);
				sum = sum + grad[d];
				
				grad_acc[d] = sum;
				std::cout << grad[d] << " ";
			}
			std::cout << std::endl;
			delete[] grad;
			delete[] grad_acc;*/
			
			
			double* sumStore = new double[dimNA]();
			
			for (int d = 0; d < dimNA; d++)
			{
				double sumquer = 0;
				
				for (int w = -sizeWindow; w < sizeWindow; w++)
				{
					sumquer = sumquer + grid_left[i][currentSlice->needleAxisInSlice[d] + w][d];
					
				}
				sumStore[d] = sumquer;
			}
			delete[] sumStore;




			
			int nr1 = std::distance(sumStore, std::min_element(sumStore, sumStore + maxnr));
		
			int nr2 = std::distance(sumStore, std::min_element(sumStore + maxnr, sumStore + dimNA));

			std::cout <<"manr " << maxnr<<" nr " << nr1 <<" nr2 " << nr2<< std::endl;
			
			//ende besser machen-> das ist urpsrünglich
			int while_y = dimNA-1;
			//while(while_y >  (dimNA - nr1))// && (grid_left[i][currentSlice->needleAxisInSlice[while_y]][while_y+1] - grid_left[i][currentSlice->needleAxisInSlice[while_y]][while_y]) < 2)
			//{
			//	grid_left[i][currentSlice->needleAxisInSlice[while_y]][while_y] = 0;
			//	grid_right[i][currentSlice->needleAxisInSlice[while_y]][while_y] = 0;
			//	while_y--;
			//}

			for (int f = 0; f < dimNA; f++)
			{
				//std::cout << f << std::endl;
				if (f < nr1 || f > nr2)
				{
					grid_left[i][currentSlice->needleAxisInSlice[f]][f] = 0;
					grid_right[i][currentSlice->needleAxisInSlice[f]][f] = 0;
				}

			}

//Alternative
			//while_y = 0;
			//while (while_y < nr2-2)// && (grid_left[i][currentSlice->needleAxisInSlice[while_y]][while_y+1] - grid_left[i][currentSlice->needleAxisInSlice[while_y]][while_y]) < 2)
			//{
			//	grid_left[i][currentSlice->needleAxisInSlice[while_y]][while_y] = 0;
			//	grid_right[i][currentSlice->needleAxisInSlice[while_y]][while_y] = 0;
			//	while_y++;
			//}

			//while_y = dimNA - 1;
			//while (while_y > dimNA- nr1)// && (grid_left[i][currentSlice->needleAxisInSlice[while_y]][while_y+1] - grid_left[i][currentSlice->needleAxisInSlice[while_y]][while_y]) < 2)
			//{
			//	grid_left[i][currentSlice->needleAxisInSlice[while_y]][while_y] = 0;
			//	grid_right[i][currentSlice->needleAxisInSlice[while_y]][while_y] = 0;
			//	while_y--;
			//}



			//mittelinie
			std::cout << "mitelllinie" << std::endl;
			for (int y = 0; y < dimNA; y++)
			{
			/*	grid_iso[GlobalVariables::pointOnNeedle_volume[0]][y] = UINT8_MAX;
				grid_iso[GlobalVariables::pointOnNeedle_volume[0] - 1][y] = UINT8_MAX/4;
				grid_iso[GlobalVariables::pointOnNeedle_volume[0] + 1][y] = UINT8_MAX/4;*/
			//std::cout << "a " << y << " " << currentSlice->needleAxisInSlice[y] << std::endl;
				grid_left[i][currentSlice->needleAxisInSlice[y] - 1][y] = UINT8_MAX;
				grid_right[i][currentSlice->needleAxisInSlice[y] + 1][y] = UINT8_MAX;
			}
			//anfang verschoben
			for (int d = 0; d < 5; d++)
			{
				/*grid_left[i][GlobalVariables::pointOnNeedle_volume[0] - d][0] = 0;
				grid_right[i][GlobalVariables::pointOnNeedle_volume[0] + d][0] = 0;*/

				//grid_iso[GlobalVariables::pointOnNeedle_volume[0] - d][0] = UINT8_MAX / 2;
				//grid_iso[GlobalVariables::pointOnNeedle_volume[0] + d][0] = UINT8_MAX / 2;
			}

		}
	}
	
	else
	{
		start = { 0, dimNA - 1 - currentSlice->needleAxisInSlice[0] };
		end = { dimNA-1, dimNA - 1 - currentSlice->needleAxisInSlice[dimNA-1] };
	
		std::cout << "isNeedleYAxisAligned == false" << std::endl;

		for (int i = 0; i < nrIsoL; i++)
		{
			//median filter
			for (int x = 1; x < dim[0]-1; x++)
			{
				for (int y = 1; y < dim[1]-1; y++)
				{
					window[0] = grid[i][x - 1][y - 1];
					window[1] = grid[i][x + 1][y - 1];
					window[2] = grid[i][x][y - 1];
					window[3] = grid[i][x - 1][y];
					window[4] = grid[i][x + 1][y];
					window[5] = grid[i][x][y];
					window[6] = grid[i][x - 1][y + 1];
					window[7] = grid[i][x + 1][y + 1];
					window[8] = grid[i][x][y + 1];

					//sort window array
					insertionSort(window, 9);
					//put the median to the new array
					grid_left[i][x][y] = window[4];
					grid_right[i][x][y] = window[4];
				}
			}
			

			double* grad = new double[dimNA-1]();
			double* grad_acc = new double[dimNA-1]();
			double sum = 0;
			for (int d = 0; d < dimNA-1; d++)
			{
				grad[d] = (double)abs(grid_left[i][d][dimNA - 1 - currentSlice->needleAxisInSlice[d]] - grid_left[i][d][dimNA - 1 - currentSlice->needleAxisInSlice[d] + 1]);
				sum = sum + grad[d];
				std::cout << grad[d] << " ";
				grad_acc[d] = sum;
			}
			std::cout <<std::endl;
			delete[] grad;
			delete[] grad_acc;
			


			int while_x = dimNA-1;
			while (while_x > dimNA/2 )//&& (grid_left[i][while_x][dimNA - 1 - currentSlice->needleAxisInSlice[while_x] +1 ] - grid_left[i][while_x][dimNA - 1 - currentSlice->needleAxisInSlice[while_x]]) < 2)
			{
				grid_left[i][while_x][dimNA - 1 - currentSlice->needleAxisInSlice[while_x]]= 0;
				grid_right[i][while_x][dimNA - 1 - currentSlice->needleAxisInSlice[while_x]] = 0;
				while_x--;
			}
			//Debugger::debugViewer(grid[i], "grid");
			//Debugger::debugViewer(grid_left[i], "gridleft");


			//mittelinie
		//	std::cout << "mitelllinie" << std::endl;
			for (int x = 0; x < dimNA; x++) 
			{
				/*grid_iso[x][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2]] = UINT8_MAX;
				grid_iso[x][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2]- 1] = UINT8_MAX/4;
				grid_iso[x][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] + 1] = UINT8_MAX / 4;*/

				//std::cout << "a " <<x<<" "<< currentSlice->needleAxisInSlice[x] << std::endl;
				grid_left[i][x][currentSlice->needleAxisInSlice[x] - 1] = UINT8_MAX;
				grid_right[i][x][currentSlice->needleAxisInSlice[x] + 1] = UINT8_MAX;
			}

		//	std::cout << "mitelllinieende" << std::endl;
			//start pixel verbreitern
			for (int d = 0; d < 5; d++)
			{
			/*	grid_left[i][0][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] - d] = 0;
				grid_right[i][0][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] + d] = 0;*/

				//grid_iso[0][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] + d] = UINT8_MAX / 2;
				//grid_iso[0][dim[1] - 1 - GlobalVariables::pointOnNeedle_volume[2] - d] = UINT8_MAX / 2;
			}

		}
	}

	//for (int i = 0; i < nrIsoL; i++)
	//{

	//	grid_left[i][start.first][start.second] = UINT8_MAX / 2;
	//	grid_left[i][end.first][end.second] = UINT8_MAX / 2;

	//	grid_right[i][start.first][start.second] = UINT8_MAX / 2;
	//	grid_right[i][end.first][end.second] = UINT8_MAX / 2;

	//}
	vector<std::pair<int, int>> res_left[nrIsoL];
	vector<std::pair<int, int>> res_right[nrIsoL];

	vector<pair<double, double>> distances[nrIsoL];
	//vector<double> distance_avag; 

	pair<double, double> max_d[nrIsoL];
	double max_d_avag = 0;

	int path_size;
	index = 0;
	std::cout << "start search" << std::endl;
	for (int i = 0; i < nrIsoL; i++)
	{
		//Debugger::debugViewer(grid[i], "grid");
		//Debugger::debugViewer(grid_left[i],"gridleft");
		//Debugger::debugViewer(grid_right[i], "gridright");
		max_d[i] = { 0,0 };

		//std::cout << start.first << " " << start.second << " " << end.first << " " << end.second << std::endl;
		pf.aStarSearch(grid_left[i], start, end ,res_left[i]);
		pf.aStarSearch(grid_right[i], start, end, res_right[i]);
		//std::cout <<"res.size"<< res_left[i].size() << " " << res_right[i].size() << std::endl;
		
		//readjust border line 
		/*for (int o = 0; o < dim[0]; o++)
		{

			if (isNeedleYAxisAligned)
			{
				int y_ = dim[1] - o - 1;
				border_left[o] = std::pair<int, int>(res_left[i][o].first, y_);
				border_right[o] = std::pair<int, int>(res_right[i][o].first, y_);
				grid_left[i][border_left[o].first][border_left[o].second] = UINT8_MAX;
				grid_right[i][border_right[o].first][border_right[o].second] = UINT8_MAX;

			}
			else
			{
				border_left[o] = std::pair<int, int>(o, res_left[i][o].second);
				border_right[o] = std::pair<int, int>(o, res_left[i][o].second);
				grid_left[i][border_left[o].first][border_left[o].second] = UINT8_MAX;
				grid_right[i][border_right[o].first][border_right[o].second] = UINT8_MAX;
			}
		}*/



		path_size = res_left[i].size();
		
		std::cout << "path size " << path_size << std::endl;

		


		distances[i].resize(path_size);
		currentSlice->distance_avag.resize(path_size);

		
		//std::cout << "Distances" << std::endl;
		for (int o = 0; o < path_size; o++)
		{
			if (isNeedleYAxisAligned)
			{
				distances[i][o].first = abs(currentSlice->needleAxisInSlice[o] - res_left[i][o].first);
				distances[i][o].second = abs(currentSlice->needleAxisInSlice[o] - res_right[i][o].first);
				
			}
			else
			{
				distances[i][o].first = abs(dimNA - 1 - currentSlice->needleAxisInSlice[o] - res_left[i][o].second);
				distances[i][o].second = abs(dimNA - 1 - currentSlice->needleAxisInSlice[o] - res_right[i][o].second);
			}
			currentSlice->distance_avag[o] += (distances[i][o].first + distances[i][o].second);
			//std::cout <<"( " << distances[i][o].first << " "<<distances[i][o].second <<" ) , ";
			if (distances[i][o].first > max_d[i].first)
			{
				max_d[i].first = distances[i][o].first;
			}
			if (distances[i][o].second > max_d[i].second)
			{
				max_d[i].second = distances[i][o].second;
			}
			if (currentSlice->distance_avag[o] > max_d_avag)
			{
				max_d_avag = currentSlice->distance_avag[o];
			}
		}
		//std::cout <<"maxDist: "<< std::endl << max_d[i].first << " " << max_d[i].second << std::endl;
		//std::cout << "normal Distances" << std::endl;
		for (int o = 0; o < path_size; o++)
		{
			if (max_d[i].first == 0)
				distances[i][o].first = 0;
			else
				distances[i][o].first /= max_d[i].first;

			if(max_d[i].second == 0)
				distances[i][o].second = 0;
			else
				distances[i][o].second /= max_d[i].second;
			

			//std::cout << "( " << distances[i][o].first << " , " << distances[i][o].second << " ) ,";
		}
		//std::cout<<std::endl;







		for (int o = 0; o < path_size ;o++)
		{
			/*if (grid_iso[res_left.at(o).first][res_left.at(o).second] == 0)
			{*/
			if (distances[i][o].first != 0)
			{
				grid_iso[res_left[i].at(o).first][res_left[i].at(o).second] = isolinesValues[i];

				int y_ = res_left[i].at(o).second; //dim[1] - res_left[i].at(o).second - 1;
				voxelWorldPosition = voxel2WorldConverter.transform(res_left[i].at(o).first, y_);
				voxelVolumePosition = world2VolumeConverter.transformTranslationFirst(voxelWorldPosition[0], voxelWorldPosition[1], voxelWorldPosition[2]);


				//Isothermpunkte
				currentSlice->measPointsInVolume.push_back({ (int)round(voxelVolumePosition[0]) ,(int)round(voxelVolumePosition[1]),(int)round(voxelVolumePosition[2]) });
			
				currentSlice->measpoints.push_back(isolinesValues[i]);




				//currentSlice->measPointsInVolume[index][0] = (int)round(voxelVolumePosition[0]);
				//currentSlice->measPointsInVolume[index][1] = (int)round(voxelVolumePosition[1]);
				//currentSlice->measPointsInVolume[index][2] = (int)round(voxelVolumePosition[2]);

				//currentSlice->measpoints[index] = isolinesValues[i];
				//index++;
			}


			//}

			


		}
		for (int o = 0; o < path_size; o++)
		{
			if (distances[i][o].second != 0)
			{
				grid_iso[res_right[i].at(o).first][res_right[i].at(o).second] = isolinesValues[i];


				int y_ = res_right[i].at(o).second;//dim[1] - res_right[i].at(o).second - 1;
				voxelWorldPosition = voxel2WorldConverter.transform(res_right[i].at(o).first, y_);
				voxelVolumePosition = world2VolumeConverter.transformTranslationFirst(voxelWorldPosition[0], voxelWorldPosition[1], voxelWorldPosition[2]);

				//Isothermpunkte
				currentSlice->measPointsInVolume.push_back({ (int)round(voxelVolumePosition[0]) ,(int)round(voxelVolumePosition[1]),(int)round(voxelVolumePosition[2]) });

				currentSlice->measpoints.push_back(isolinesValues[i]);

				/*currentSlice->measPointsInVolume[index][0] = (int)round(voxelVolumePosition[0]);
				currentSlice->measPointsInVolume[index][1] = (int)round(voxelVolumePosition[1]);
				currentSlice->measPointsInVolume[index][2] = (int)round(voxelVolumePosition[2]);

				currentSlice->measpoints[index] = isolinesValues[i];
				index++;*/

			}
			
			
		}


	}
	std::cout << "maindirVol " << GlobalVariables::needleDir[maindirVol] << std::endl;
	voxelWorldPosition = voxel2WorldConverter.transform(start.first, start.second);
	voxelVolumePosition = world2VolumeConverter.transformTranslationFirst(voxelWorldPosition[0], voxelWorldPosition[1], voxelWorldPosition[2]);
	std::cout << "voxVolstart: " << voxelVolumePosition[0] << " " << voxelVolumePosition[1] << " " << voxelVolumePosition[2] << std::endl;
	
	int offset_start = static_cast<int>(round(voxelVolumePosition[maindirVol]));

	voxelWorldPosition = voxel2WorldConverter.transform(end.first, end.second);
	voxelVolumePosition = world2VolumeConverter.transformTranslationFirst(voxelWorldPosition[0], voxelWorldPosition[1], voxelWorldPosition[2]);
	std::cout << "voxVolend: " << voxelVolumePosition[0] << " " << voxelVolumePosition[1] << " " << voxelVolumePosition[2] << std::endl;
	int offset_end = static_cast<int>(round(voxelVolumePosition[maindirVol]));

	int offset = std::min(offset_start,offset_end);

	int maxHeatInSlice = encode(currentSlice->data->GetScalarRange()[1]);

	int maxHeatNeedle = std::distance(currentSlice->distance_avag.begin(), std::max_element(currentSlice->distance_avag.begin(), currentSlice->distance_avag.end()));
	std::cout << "maxNeedle " << maxHeatNeedle << std::endl;


	//std::cout << "maxHeat " << maxHeatInSlice << std::endl;
	
	//if (offset_start > offset_end)
	//{
	//	for (int o = 0; o < (path_size); o++)
	//	{
	//		if (max_d_avag == 0)
	//			currentSlice->distance_avag[o] = 0;
	//		else
	//			currentSlice->distance_avag[o] /= max_d_avag;
	//		//tempAlongNeedle[offset + o] =  GlobalVariables::baseLineT + (int) (currentSlice->distance_avag[o] * ((maxHeatInSlice > 90.0 ? 90.0 :maxHeatInSlice)- GlobalVariables::baseLineT));
	//		q_relAlongNeedle[offset_start - o] = currentSlice->distance_avag[o];
	//		tempAlongNeedle[offset_start - o] = GlobalVariables::baseLineT + (int)(currentSlice->distance_avag[o] * (120 - GlobalVariables::baseLineT));
	//		std::cout << "( " << currentSlice->distance_avag[o] <<" "<<tempAlongNeedle[offset_start - o] << " ) ,";
	//	}
	//}

	if (GlobalVariables::needleDir[maindirVol] < 0)
	{
		std::cout << "reverse" << std::endl;
		std::reverse(currentSlice->distance_avag.begin(), currentSlice->distance_avag.end());
		offset = path_size - offset_start;//oder logisch wäre q_relAlongNeedle.size()
	}
	std::cout <<"pathsize:  "<<path_size<< "offset: " << offset << "start " << offset_start << "end " << offset_end << "tempTsize : " << tempAlongNeedle.size() << std::endl;
	/*tempAlongNeedle.clear();
	q_relAlongNeedle.clear();*/
	for (int o = 0; o < (q_relAlongNeedle.size()); o++)
		q_relAlongNeedle[o] = 0;
	for (int o = 0; o < (path_size); o++)
	{
		//std::cout << "o " << offset + o << std::endl;;
		if (max_d_avag == 0)
			currentSlice->distance_avag[o] = 0;
		else
			currentSlice->distance_avag[o] /= max_d_avag;
		//tempAlongNeedle[offset + o] =  GlobalVariables::baseLineT + (int) (currentSlice->distance_avag[o] * ((maxHeatInSlice > 90.0 ? 90.0 :maxHeatInSlice)- GlobalVariables::baseLineT));
		q_relAlongNeedle[offset + o ] = currentSlice->distance_avag[o];
		tempAlongNeedle[offset + o] = GlobalVariables::baseLineT + (int)(currentSlice->distance_avag[o] * (120 - GlobalVariables::baseLineT));
		//std::cout << "( " << currentSlice->distance_avag[o] << " " << tempAlongNeedle[o] << " ) ,";
	
	}
	for (int o = 0; o < tempAlongNeedle.size(); o++)
		std::cout << "( " << tempAlongNeedle[o] << " ) ,";

	std::cout << std::endl;
	std::map<int, int> m;// = { {30,50},{35,100},{40,150},{45,200},{50,250} };
	for (int i = 0; i < nrIsoL; i++)
	{
		m.insert(std::pair<int,int>(isolinesValues[i], 25 * i));
	}
	cv::String str = "gridiso_" + cv::String(_angle.toStdString()) + "_" + cv::String(_timeStep.toStdString());
	
	grid_iso[start.first][start.second] = UINT8_MAX;
	grid_iso[end.first][end.second] = UINT8_MAX;
	//Debugger::debugViewer(grid_iso, str);//1,&m

	// 
	//std::cout << "Heatsum Slice: " << Heatsum / (60 * 60) << std::endl;
	//slicesBuffer[index] =  Slice(imageDicom, voxel2WorldConverter, timestep, positionInVolume);
	/*for (int x = 0; x < 60; x++)
	{
		for (int y = 0; y < 60; y++)
		{
			for (int z = 0; z < 60; z++)
			{
				mytype* voxelVolume = static_cast<mytype*>(volume->GetScalarPointer(
					x,
					y,
					z));
				if (x == round(GlobalVariables::pointOnNeedle_volume[0]) && z == round(GlobalVariables::pointOnNeedle_volume[2]))
				{
					voxelVolume[0] = 50;
					index_++;
				}
			}
		}
	}*/

	//std::cout << "index: " << index_ << std::endl;
	volume->Modified();
	//updateTempNeedle();
	delete[] voxelVolumePosition;
	delete[] voxelWorldPosition;
	//imageDicom->Delete();
	//Matrix3D<double>* m = new Matrix3D<double>(60, 60, 1);
	//for (int i = 0; i < 60; i++)
	//{
	//	for (int j = 0; j < 60; j++)
	//	{
	//		(*m)[{i, j, 0}] = grid_left[0][i][j];
	//	}
	//}
	//Debugger::debugViewer(m);
	std::cout<<std::endl << "DataVolume::addSlice - Done" << std::endl;
}

void DataVolume::display(vtkSmartPointer<vtkRenderWindow> _q_vtk_renderWindow, bool _isCropped, bool _isCameraAbove)
{
	qDebug() << "DataVolume::display";
	int low = 10;//10
	//int middle = 1;
	int high = 200;//220
	/*
	*	vtkNamedColors bietet eine einfache Farbdatenbank mit entsprechenden Namen, die man sich einfach ausgeben lassen kann.
	*	Beim Setzen der Hintergrundfarbe habe ich zum Beispiel die Farbe "Wheat" genommen.
	*/
	vtkSmartPointer<vtkNamedColors> namedColors = vtkSmartPointer<vtkNamedColors>::New();
	/*
	*	Hier wird das entsprechende RenderWindow und der dazu gehrige Renderer geholt.
	*/
	vtkSmartPointer<vtkRenderer> q_vtk_renderer = vtkSmartPointer<vtkRenderer>::New();
	_q_vtk_renderWindow->AddRenderer(q_vtk_renderer);
	/*
	*	Hier gibts den Interactor und den 3D InteractorStyle.
	*/
	vtkSmartPointer<vtkRenderWindowInteractor> q_vtk_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> q_vtk_interactorStyle = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	q_vtk_interactor->SetInteractorStyle(q_vtk_interactorStyle);
	q_vtk_interactor->SetRenderWindow(_q_vtk_renderWindow);
	q_vtk_interactor->Modified();
	/*
	*	Als Mapper hab ich mich in der Flle von Mglichkeiten fr den OpenGL GPU Volume Cast Mapper entschieden, der in der Lage ist
	*	IsoSurfaces vernnftig darzustellen weil die entsprechenden Shader bereits integriert sind.
	*/
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volumeMapper =vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	//vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper> volumeMapper = vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper>::New();
	//vtkFixedPointVolumeRayCastMapper
	volumeMapper->SetInputData(volume);

	/*
	*	Die automatische Sample Distanz knnte je nach Auflsung des Datensatzes zu Problemen fhren. Ich hab sie jetzt auf 0.5 gesetzt aber
	*	da knnt ihr einfach mal ein bisschen rumspielen mit den Einstellungen. Wann sieht es gut aus und ist trotzdem performant?
	*	Anschlieend wird der BlendMode auf IsoSurface gesetzt. Da knnt ihr aber auch nochmal rumspielen, vielleicht findet ihr was cooleres.
	*/
	volumeMapper->AutoAdjustSampleDistancesOff();
	volumeMapper->SetSampleDistance(0.5);
	//volumeMapper->SetBlendModeToIsoSurface();
	volumeMapper->Modified();
	/*
	*	Jetzt wird es ekelhaft. Um vernnftig zu visualisieren msst ihr eine richtige Transferfunktion definieren fr die Isowerte. Das geschieht mit
	*	zwei sogenannten LookupTables. In der vtkColorTransferFunction definiert ihr, welcher Voxelwert welche Farbe bekommt. In diesem Fall ist alles
	*	mit dem Wert 10 grn und alles mit dem Wert 220 rot.
	*/
	vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
	color->RemoveAllPoints();
	for (int i = 0; i < low; i++)
		color->AddRGBPoint(i, 0.0, 1.0, 0.0);
	//color->AddRGBPoint(middle, 0.0, 0.0, 1.0);
	for (int j = low; j < high; j++)
		color->AddRGBPoint(j, 1.0, 0.0, 0.0);
	color->Modified();
	//color->AddRGBPoint(1, 0.0, 1.0, 0.0);
	//color->AddRGBPoint(0, 1.0, 0.0, 0.0);
	/*
	*	Mit der vtkPiecewiseFunction definiert ihr die Opacitt eines entsprechenden Wertes. Die Punkte, die ihr hier eintragt mssen identisch sein mit den
	*	Werten in der vtkColorTransferFunction. Hier haben Voxel mit dem Wert 10 eine Opazitt von 30% und Voxel mit einem Wert von 220 eine Opazitt von 60%.
	*/
	vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = vtkSmartPointer<vtkPiecewiseFunction>::New();
	compositeOpacity->RemoveAllPoints();

	for (int i = 0; i < low; i++)
		compositeOpacity->AddPoint(i, 0.0);
	//compositeOpacity->AddPoint(middle, 0.3);
	for (int j = low; j < high; j++)
		compositeOpacity->AddPoint(j, 1.0);
	compositeOpacity->Modified();
	//compositeOpacity->AddPoint(1, 0.3);
	//compositeOpacity->AddPoint(220, 0.6);
	/*
	*	Jetzt mssen wir die VolumeProperty definieren, dazu stellen wir das Shading an und eine einfache lineare interpolation. Bei der Interpolation gibt
	*	es bestimmt auch noch andere Mglichkeiten. SetColor bekommt die vtkColorTransferFunction und setScalarOpacity bekommt die
	*	vtkPiecewiseFunction.
	*/
	vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->ShadeOn();
	volumeProperty->SetInterpolationTypeToLinear();
	volumeProperty->SetColor(color);
	volumeProperty->SetScalarOpacity(compositeOpacity);
	volumeProperty->Modified();
	/*
	*	Dann mssen wir das Volumen definieren mit dem entsprechenden Mapper und den vorab eingestellenten Properties.
	*/
	vtkSmartPointer<vtkVolume> vtk_volume = vtkSmartPointer<vtkVolume>::New();
	vtk_volume->SetMapper(volumeMapper);
	vtk_volume->SetProperty(volumeProperty);
	vtk_volume->Modified();
	/*
	*	Jetzt schmeien wir alles in den Renderer und setzen am Ende das DepthPeeling auf True, damit der Renderer wei
	*	in welcher Reihenfolge er die Strukturen visualisieren muss.
	*/
	q_vtk_renderer->AddVolume(vtk_volume);
	q_vtk_renderer->SetBackground(namedColors->GetColor3d("Wheat").GetData());//
	//q_vtk_renderer->GetActiveCamera()->UpdateViewport(q_vtk_renderer);
	//q_vtk_renderer->SetViewPoint(0, 1, 0);

	if (false)
	{
		q_vtk_renderer->GetActiveCamera()->SetFocalPoint(vtk_volume->GetCenter());
		q_vtk_renderer->GetActiveCamera()->SetPosition(0, 200, 0);
	}


	q_vtk_renderer->SetUseDepthPeeling(true);
	/*
	*	Ganz am Ende weisen wir unseren entsprechenden Isowerten auch noch die Isovalues zu fr die interne Verarbeitung.
	*	Hier reprsentiert unsere IsoSurface 0 den Wert 10 und unsere IsoSurface 1 den Wert 220. Auch diese Werte mssen
	*	identisch sein mit den Werten in der Transfer Funktion.
	*/
	volumeProperty->GetIsoSurfaceValues()->SetValue(0, low);
	volumeProperty->GetIsoSurfaceValues()->SetValue(2, high);

	_q_vtk_renderWindow->Render();
	qDebug() << "DataVolume::display - Done";

}


 


void DataVolume::crop(vtkSmartPointer<vtkImageData> volume, vtkSmartPointer<vtkImageData> croppedvolume, float percentage)
{
	//qDebug() << "Constructor Data Volume";
	

	int newSize = (int) volumeSize * percentage;
	
	int SizeDiff = (int) (volumeSize - newSize)/2;

	croppedvolume->SetDimensions(newSize, volumeSize, newSize);
	croppedvolume->SetOrigin(0, 0, 0);

#if VTK_MAJOR_VERSION <= 5
#else
	croppedvolume->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
#endif

	for (int x = 0; x < newSize; x++)
	{
		for (int y = 0; y < volumeSize; y++)
		{
			for (int z = 0; z < newSize; z++)
			{
			
					mytype* voxelCroppedVolume = static_cast<mytype*>(croppedvolume->GetScalarPointer(x, y, z));
					mytype* voxelVolume = static_cast<mytype*>(volume->GetScalarPointer(x + SizeDiff, y , z + SizeDiff));
					voxelCroppedVolume[0] = voxelVolume[0];
			}
		}
	}
	croppedvolume->Modified();
}

void DataVolume::crop(float percentage)
{
	int newSize = (int)volumeSize * percentage;
	//qDebug() << "newSize" << newSize;
	int SizeDiff = (int)(volumeSize - newSize) / 2;
	//qDebug() << "SizeDiff" << SizeDiff;
	vtkSmartPointer<vtkImageData> tempVolume = vtkSmartPointer<vtkImageData>::New();
	tempVolume->DeepCopy(volume);
	volume->Delete();

	volume = vtkSmartPointer<vtkImageData>::New();
	volume->SetDimensions(newSize, newSize, newSize);
	volume->SetSpacing(tempVolume->GetSpacing());
	volume->SetOrigin(tempVolume->GetOrigin());
	
	volume->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

	for (int x = 0; x < newSize; x++)
	{
		for (int y = 0; y < newSize; y++)
		{
			for (int z = 0; z < newSize; z++)
			{

				mytype* voxelCroppedVolume = static_cast<mytype*>(volume->GetScalarPointer(x, y, z));
				mytype* voxelVolume = static_cast<mytype*>(tempVolume->GetScalarPointer(x + SizeDiff, y + SizeDiff, z + SizeDiff));
				voxelCroppedVolume[0] = voxelVolume[0];
			}
		}
	}
	volume->Modified();
}

vtkSmartPointer<vtkImageData> DataVolume::getcurrentSlice()
{
	return currentSlice->data;
}





vtkSmartPointer<vtkImageData> DataVolume::getImageData()
{
	return volume;
}









void DataVolume::setAllVoxel2Value(int _value)
{
	//qDebug() << "DataVolume::setAllVoxel2Value";
	int* dim = volume->GetDimensions();
	for (unsigned int x = 0; x < dim[0]; x++)
	{
		for (unsigned int y = 0; y < dim[1]; y++)
		{
			for (unsigned int z = 0; z < dim[2]; z++)
			{
				mytype* pixel = static_cast<mytype*>(volume->GetScalarPointer(x, y, z));

				pixel[0] = _value;
				
			}
		}
	}
	volume->Modified();
	//qDebug() << "DataVolume::setAllVoxel2Value - Done";
}

void DataVolume::createCoordNeedleAxis(DicomHandler::dicomDataProperties _ref, DicomHandler::dicomDataProperties _in)
{
	//(0,1,0) needle direction
	float eta_ref[3];
	float eta_in[3];
	crossProduct(_ref.imageOrientationX, _ref.imageOrientationY, eta_ref);
	crossProduct(_in.imageOrientationX, _in.imageOrientationY, eta_in);

	float needleDir[3];
	crossProduct(eta_ref, eta_in, GlobalVariables::needleDir);




	float d_ref = dotProduct(eta_ref, _ref.imagePosition);
	float d_in = dotProduct(eta_in, _in.imagePosition);

	float q[3];
	float a = (d_ref * dotProduct(eta_in,eta_in) - d_in * dotProduct(eta_ref, eta_in)) / (dotProduct(eta_ref,eta_ref) *dotProduct(eta_in, eta_in) - (dotProduct(eta_ref, eta_in) * dotProduct(eta_ref, eta_in)));
	float b = (d_in *dotProduct(eta_ref,eta_ref) - d_ref * dotProduct(eta_ref, eta_in)) / (dotProduct(eta_ref, eta_ref) * dotProduct(eta_in, eta_in) - (dotProduct(eta_ref, eta_in) * dotProduct(eta_ref, eta_in)));
	for (int i = 0; i < 3; i++)
	{
		q[i] = a * eta_ref[i] + b * eta_in[i];
 		GlobalVariables::pointOnNeedle_world[i] = a * eta_ref[i] + b * eta_in[i];
	}
	

}
void DataVolume::defaultTempAlongNeedle()
{
	for(int i=0; i< tempAlongNeedle.size();i++)
	{
		tempAlongNeedle[i] = GlobalVariables::baseLineT;
		q_relAlongNeedle[i] = 1.0;
	}
}
//#pragma optimize( "", on )
void DataVolume::updateTempNeedle()
{
	std::cout << "updateTempNeedle" << std::endl;
	int index;
	std::array<int, 9> heatNeighbors;
	std::cout << std::endl;
	for (int y = 0; y < volume->GetDimensions()[1]; y++)
	{
		index = 0;

		
		for (int x = GlobalVariables::pointOnNeedle_volume[0]-1; x <= GlobalVariables::pointOnNeedle_volume[0] +1; x++)
		{
			for (int z= GlobalVariables::pointOnNeedle_volume[2] -1; z <= GlobalVariables::pointOnNeedle_volume[2]+1; z++)
			{
				//std::cout << "xyz: " << x <<" " << y <<" " << z;
				mytype* voxelVolume = static_cast<mytype*>(volume->GetScalarPointer(
					x,
					y,
					z));
				
				int u =(int) voxelVolume[0];
				//std::cout <<" " <<u << " " << std::endl;
				heatNeighbors.at(index) = u;
				index++;
			}
		}

		//std::cout << std::endl;
		int maxH = *std::max_element(heatNeighbors.begin(), heatNeighbors.end());
		std::cout << maxH << " ";
		//array brauch man theoretisch nicht mehr 
		tempAlongNeedle[y] = maxH;
	}
	std::cout << std::endl;
}
//int* DataVolume::getTempAlongNeedle()
//{
//	return tempAlongNeedle;
//}
void DataVolume::crossProduct(float v_A[], float v_B[], float c_P[]) {
	c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
	c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
	c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}
float DataVolume::dotProduct(float v_A[], float v_B[]) {
	return v_A[0] * v_B[0] + v_A[1] * v_B[1] + v_A[2] * v_B[2];
}
float DataVolume::l2Norm(float* a, float* b, int N)
{
	float accum = 0;
	for (int i = 0; i < N; i++)
	{
		accum += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return sqrt(accum);
}
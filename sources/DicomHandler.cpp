#include "DicomHandler.h"
#include "vtkImageReslice.h"
#include<CoordinatesConverter.h>
DicomHandler::DicomHandler()
{
	//reader = vtkSmartPointer<vtkDICOMImageReader>::New();
}


DicomHandler::~DicomHandler()
{
}

vtkSmartPointer<vtkImageData> DicomHandler::loadDicom(QString q_c_dicomFilePath)
{
	//qDebug() << "DicomHandler::loadDicom";
	//vtkDICOMReader
	//vtkDICOMImageReader
	vtkSmartPointer<vtkDICOMReader> reader = vtkSmartPointer<vtkDICOMReader>::New();
	
	reader->SetFileName(q_c_dicomFilePath.toStdString().c_str());
	reader->Update();
	//qDebug() << "DicomHandler::loadDicom  - Done";
	return  reader->GetOutput();
}

void DicomHandler::getImageDataProperties(vtkSmartPointer<vtkDICOMReader> _reader, dicomDataProperties* _properties)
{
	//qDebug() << "DicomHandler::getImageDataProperties -Reader";
	vtkSmartPointer<vtkDICOMMetaData> metaData = _reader->GetMetaData();

	if (metaData->HasAttribute(DC::EchoTime))
		_properties->echoTime = metaData->Get(DC::EchoTime).AsFloat();
	if (metaData->HasAttribute(DC::MagneticFieldStrength))
		_properties->magneticFieldStrength = metaData->Get(DC::MagneticFieldStrength).AsFloat();

	metaData->Get(DC::ImagePositionPatient).GetValues(_properties->imagePosition, 3);
	//auto temp = metaData->Get(DC::PixelSpacing);
	if(metaData->HasAttribute(DC::PixelSpacing))
		metaData->Get(DC::PixelSpacing).GetValues(_properties->pixelSpacing, 2);
	if (metaData->HasAttribute(DC::SliceThickness))
		metaData->Get(DC::SliceThickness).GetValues(&_properties->pixelSpacing[2], 1);
	if (metaData->HasAttribute(DC::Columns))
		metaData->Get(DC::Columns).GetValues(&_properties->dimension[0], 1);
	if (metaData->HasAttribute(DC::Rows))
		metaData->Get(DC::Rows).GetValues(&_properties->dimension[1], 1);
	//testValue = metaData->Get(DC::PixelSpacing);
	float* orientation = new float[6];
	if (metaData->HasAttribute(DC::ImageOrientationPatient))
		metaData->Get(DC::ImageOrientationPatient).GetValues(orientation, 6);
	_properties->imageOrientationX[0] = orientation[0];
	_properties->imageOrientationX[1] = orientation[1];
	_properties->imageOrientationX[2] = orientation[2];

	_properties->imageOrientationY[0] = orientation[3];
	_properties->imageOrientationY[1] = orientation[4];
	_properties->imageOrientationY[2] = orientation[5]; 
	_properties->pixelSpacing[2] = 5;
	//qDebug() << "Position:		" << _properties->imagePosition[0] << _properties->imagePosition[1] << _properties->imagePosition[2];
	//qDebug() << "Orientation X:		" << _properties->imageOrientationX[0] << _properties->imageOrientationX[1] << _properties->imageOrientationX[2];
	//qDebug() << "Orientation Y:		" << _properties->imageOrientationY[0] << _properties->imageOrientationY[1] << _properties->imageOrientationY[2];
	//qDebug() << "Spacing:		" << _properties->pixelSpacing[0] << _properties->pixelSpacing[1] << _properties->pixelSpacing[2];
	//qDebug() << "Echotime:		" << _properties->echoTime;
	//qDebug() << "fStrength:		" << _properties->magneticFieldStrength;
	//qDebug() << "fStrength:		" << _properties->magneticFieldStrength;
	//qDebug() << "Columns:		" << _properties->dimension[0] << "Rows:	" << _properties->dimension[1] ;

}

void DicomHandler::getImageDataProperties(QString _fileDataName, dicomDataProperties * _properties)
{
	//qDebug() << "DicomHandler::getImageDataProperties - FileName";
	vtkSmartPointer<vtkDICOMReader> reader = vtkSmartPointer<vtkDICOMReader>::New();

	reader->SetFileName(_fileDataName.toStdString().c_str());
	reader->Update();

	DicomHandler dHandler;
	dHandler.getImageDataProperties(reader, _properties);
	
	
	//qDebug() << "DicomHandler::getImageDataProperties - Done";

}


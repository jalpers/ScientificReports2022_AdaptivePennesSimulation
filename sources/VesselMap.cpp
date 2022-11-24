#include<VesselMap.h>

VesselMap::VesselMap()
{
	map = new Matrix3D<bool>();
	
}

VesselMap::~VesselMap()
{
	delete map;
}
VesselMap::VesselMap(CoordinatesConverter _world2VolumeConverter, int* _size)
{
	
	map = new Matrix3D<bool>(_size[0], _size[1], _size[2]);
	
	map->setAll2(false);
	world2VolumeConverter = _world2VolumeConverter;
	for (int i = 0; i < 3; i++)
		volumeSize[i] = map->getDimensions()[i];
}
void VesselMap::load(std::string _filename)
{
	std::cout << "LoadVesselMap" << std::endl;
	/*vtkSmartPointer<vtkNIFTIImageReader> reader_nii =
		vtkSmartPointer<vtkNIFTIImageReader>::New();

	reader_nii->SetFileName(_filename.c_str());
	reader_nii->Update();

	reader_nii->GetImageDataInput(0);*/


	vtkSmartPointer <vtkDICOMReader> reader =
		vtkSmartPointer <vtkDICOMReader>::New();
	// Provide a multi-frame, multi-stack file
	reader->SetFileName(_filename.c_str());
	// Read the meta data, get a list of stacks
	reader->UpdateInformation();

	DicomHandler dHandler;
	DicomHandler::dicomDataProperties properties;
	dHandler.getImageDataProperties(reader, &properties);

	properties.pixelSpacing[2] = 1.1;
	CoordinatesConverter vesselMapVolume2World = CoordinatesConverter(properties);

	//https://dgobbi.github.io/vtk-dicom/doc/api/attributes.html
	// Get the arrays that map slice to file and frame.
	vtkIntArray* fileMap = reader->GetFileIndexArray();
	vtkIntArray* frameMap = reader->GetFrameIndexArray();
	// Get the image data and meta data.
	vtkImageData* image = reader->GetOutput();
	vtkDICOMMetaData* meta = reader->GetMetaData();
	// Get the number of components in the data.
	int numComponents = image->GetNumberOfScalarComponents();

	// Get the full vector dimension for the DICOM data.
	int vectorDimension = fileMap->GetNumberOfComponents();

	// Compute the samples per pixel in original files.
	int samplesPerPixel = numComponents / vectorDimension;
	samplesPerPixel = 1;

	// Check for time dimension
	int timeDimension = reader->GetTimeDimension();
	//if (timeDimension == 0)
	{
		timeDimension = 1;
	}

	//Extract an image at the desired time slot (e.g. for display).
	int componentIndex = 1 * vectorDimension / timeDimension * samplesPerPixel;

	for (int z = 0; z < 100; z++)
	{
		componentIndex = z;
		//qDebug() << "componentIndex" << componentIndex;
		vtkNew<vtkImageExtractComponents> extractor;
		extractor->SetInputConnection(reader->GetOutputPort());
		if (samplesPerPixel == 1)
		{
			extractor->SetComponents(componentIndex);
		}
		extractor->Update();

	

		vtkSmartPointer<vtkImageData> tmp = extractor->GetOutput();
		#pragma omp parallel for collapse(2)
		for (int x = 0; x < 256; x++)
		{
			for (int y = 0; y < 256; y++)
			{
				float* voxelVolumePosition = new float[3];
				float* voxelWorldPosition = new float[3];
				uchar* tmpPxl = static_cast<uchar*>(tmp->GetScalarPointer(x, y, 0));
				voxelWorldPosition = vesselMapVolume2World.transform(x, 256 - 1 - y, z);
				voxelVolumePosition = world2VolumeConverter.transformTranslationFirst(voxelWorldPosition[0], voxelWorldPosition[1], voxelWorldPosition[2]);

				if (voxelVolumePosition[0] >= 0 && voxelVolumePosition[0] < volumeSize[0]
					&& voxelVolumePosition[1] >= 0 && voxelVolumePosition[1] < volumeSize[1]
					&& voxelVolumePosition[2] >= 0 && voxelVolumePosition[2] < volumeSize[2])
				{

					(*map)[{(int)round(voxelVolumePosition[0]), (int)round(voxelVolumePosition[1]), (int)round(voxelVolumePosition[2])}]  = tmpPxl[0] > 0 ? true : false;
				}
				delete[] voxelVolumePosition;
				delete[] voxelWorldPosition;
			}
		}
		

	}
	
}
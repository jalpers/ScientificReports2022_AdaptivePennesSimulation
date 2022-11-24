#pragma once
#ifndef DEBUGGER_H
#define DEBUGGER_H
#include<opencv2/opencv.hpp>
#include<Matrix3D.h>
#include<iostream>
#include<vector>
#include <map>
#include<SimulationRessources.h>
//using namespace cv;
class Debugger
{
public:
	/*Debugger();
	~Debugger();*/
	static void  debugViewer(Matrix3D<double>* _b, cv::String _str, int _f = 1)
	{
		//m.lock();
		//##########DEBUG###############
		int* dim = _b->getDimensions();

		int sz[2] = { dim[0],dim[1] };
		cv::Mat L(2, sz, CV_8UC1, cv::Scalar::all(0));

		double heatsum = 0;

		for (int i = 0; i < dim[0]; i++)
		{
			for (int j = 0; j < dim[1]; j++)
			{
				for (int k = 0; k < dim[2]; k++)
				{
					{

						double l = (*_b)[{i, j, k}];
						//std::cout << l << " ";
						heatsum += l;

						L.at<uchar>(j, i) = (int)(l * _f);
					}


				}
			}

		//	std::cout << std::endl;

		}
		std::cout << "debugViewer- heatsum: " << heatsum / (dim[0] * dim[1] * dim[2]) << std::endl;




		cv::namedWindow(_str, cv::WINDOW_NORMAL);
		// Display the image.
		imshow(_str, L);

		// Wait for a keystroke.   
		//waitKey(0);
		while (cv::waitKey(1) != 27); // 27 = ascii value of ESC
		// Destroys all the windows created                         
		cv::destroyAllWindows();

		// Write the image in the same directory

		imwrite(_str, L);
		//m.unlock();
	}
	static void debugViewer(const std::vector<std::vector<int>>& _b, cv::String _str, int _f = 1, const std::map<int,int>* _map = nullptr )

	{
		int dim[2];
		dim[0] = _b.size();
		dim[1] = _b[0].size();

		int sz[2] = { dim[0],dim[1] };
		cv::Mat L(2, sz, CV_8UC1, cv::Scalar::all(0));

		double heatsum = 0;

		for (int i = 0; i < dim[0]; i++)
		{
			for (int j = 0; j < dim[1]; j++)
			{

				int l = _b[i][j];
				std::cout << l << " ";
				heatsum += l;

	
				if (_map != nullptr)
				{
					//std::map<int, int>::iterator it;
					auto it = _map->find(l);
					if (it == _map->end())
						L.at<uchar>(j, i) = l;
					else
						L.at<uchar>(j, i) = it->second;
				}
				else
				{
					L.at<uchar>(j, i) = l * _f;
				}



			}
			std::cout << std::endl;
		}
		std::cout << "debugViewer - heatsum: " << heatsum / (dim[0] * dim[1]) << std::endl;





		cv::namedWindow(_str,cv::WINDOW_NORMAL);

		cv::resizeWindow(_str, 600, 600);
		// Display the image.
		imshow(_str, L);

		// Wait for a keystroke.   
		cv::waitKey(0);

		// Destroys all the windows created                         
		cv::destroyAllWindows();

		// Write the image in the same directory

		//"C:/Users/maxro/Documents/MasterInformatik/Masterarbeit/Ergebnisse/DebugViewer_Bilder/"+
		cv::resize(L, L, cv::Size(600, 600), cv::INTERSECT_NONE);
		//imwrite(_str + ".jpg", L);


	}
	static void debugViewer(Matrix3D<PennesEquationParameter>* _b)
	{
		//##########DEBUG###############
		int* dim = _b->getDimensions();

		int sz[2] = { dim[0],dim[2] };
		cv::Mat L(2, sz, CV_8UC1, cv::Scalar::all(0));

		double heatsum = 0;
		double Qsum = 0;
		for (int i = 0; i < dim[0]; i++)
		{

			for (int j = 0; j < dim[1]; j++)
			{
				double l = (*_b)[{i, j, 0}].T;
				std::cout << l << " ";
				heatsum += l;
				L.at<uchar>(i, j) = (int)l;

			}

			std::cout << std::endl;
		}

		std::cout << "debugViewer-heatsum: " << heatsum / (dim[0]*  dim[1]) << std::endl;




		//Mat img_color;
		// Apply the colormap:
		//applyColorMap(L, img_color, COLORMAP_JET);

		cv::namedWindow("Display frame", cv::WINDOW_NORMAL);

		// Show the result:
		//imshow("colorMap", img_color);
		// Display the image.
		imshow("Display frame", L);

		// Wait for a keystroke.   
		cv::waitKey(0);

		// Destroys all the windows created                         
		cv::destroyAllWindows();

		// Write the image in the same directory
		//imwrite("grayscale.jpg", L);
	}

private:

};

//Debugger::Debugger()
//{
//}
//
//Debugger::~Debugger()
//{
//}

#endif // !DEBUGGER_H


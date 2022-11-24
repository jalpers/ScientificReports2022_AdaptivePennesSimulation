#include<SimulationRessources.h>

const float GlobalVariables::dt = .5;
const float GlobalVariables::timeBetweenMeas = 3.0;
float GlobalVariables::needleDir[3] = { 0 };
float GlobalVariables::pointOnNeedle_world[3] = { 0 };
int GlobalVariables::pointOnNeedle_volume[3] = { 0 };
int GlobalVariables::slice_index = 0;
int GlobalVariables::boundaryCondition = 1;//bei matrix.h muss man das auch noch ändern 
std::string GlobalVariables::angleList[GlobalVariables::numberOfSlices] = { "0","22_5","45","67_5","90","112_5","135","157_5" };
//QStringList GlobalVariables::slices_time = { "0","90","45","135_5","22_5","112_5","67_5","157_5" };
float GlobalVariables::slices_time_f[GlobalVariables::numberOfSlices] = { 0,90,45,135.5,22.5,112.5,67.5,157.5 };
ParameterFlags GlobalVariables::paramToOpti[GlobalVariables::numberOfOptiParam]{ ParameterFlags::T_ind_max, ParameterFlags::td_a};//, ParameterFlags::c_b, ParameterFlags::td_c ParameterFlags GlobalVariables::paramToOpti[GlobalVariables::numberOfOptiParam]{ ParameterFlags::T_ind_max, ParameterFlags::td_a };//, ParameterFlags::c_b, ParameterFlags::td_c 
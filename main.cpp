#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Windows.h>
#include "euler_krylov.h"
#include "ostreams.h"


int main() 
{
	std::ifstream nav_path("Nav_1.txt");
	std::ifstream init_path("init_1.txt");
	std::ifstream imu_path("IMU_1.txt");
	
	elem init = {};

	if (init_path.is_open()) {
		std::string rubbish;
		std::getline(init_path, rubbish);
		double dt, lat, lon, h, vx, vy, vz, pitch, roll, thdg;
		init_path >> dt >> lat >> lon >> h >> vx >> vy >> vz >> pitch >> roll >> thdg;
		init = elem{ dt, lat, lon, h, vx, vy, vz, pitch, roll, thdg };
	}

	std::vector<elem> elements_check = { init };
	std::vector<elem> elements = { init };
	std::vector<imu> imus;

	if (nav_path.is_open()) {
		std::string rubbish;
		std::getline(nav_path, rubbish);
		double dt, lat, lon, h, vx, vy, vz, pitch, roll, thdg;
		while (nav_path >> dt >> lat >> lon >> h >> vx >> vy >> vz >> pitch >> roll >> thdg)
		{
			elements_check.push_back(elem{ dt, lat, lon, h, vx, vy, vz, pitch, roll, thdg });
		}
	}

	if (imu_path.is_open()) {
		std::string rubbish;
		std::getline(imu_path, rubbish);
		double t, ax, ay, az, wx, wy, wz;
		while (imu_path >> t >> ax >> ay >> az >> wx >> wy >> wz)
		{
			imus.push_back(imu{ t, ax, ay, az, wx, wy, wz });
		}
	}

	for (int i = 0; i < imus.size(); i++) {
		// std::vector<double> cords = count_cords(elements[i].t, elements[i].vxg, elements[i].vyg, elements[i].vzg, elements[i].lat, elements[i].lon, elements[i].h);
		// std::vector<double> wgs = count_wgs(elements[i].vxg, elements[i].vzg, cords[1], cords[2]);
		// std::vector<double> speeds = count_speeds(elements[i].t, elements[i].vxg, elements[i].vyg, elements[i].vzg, imus[i].ax, imus[i].ay, imus[i].az, cords[1], cords[2]);
		// std::vector<double> rot_angles = count_rot_angles(elements[i].t, imus[i].wx, imus[i].wy, imus[i].wz, wgs[0], wgs[1], wgs[2], 360 - elements[i].thdg, elements[i].pitch, elements[0].roll);

		std::vector<double> wgs = count_wgs(elements[i].vxg, elements[i].vzg, elements[i].lat, elements[i].h); // {wxg, wyg, wzg}
		std::vector<double> speeds = count_speeds(elements[i].t, elements[i].vxg, elements[i].vyg, elements[i].vzg, imus[i].ax, imus[i].ay, imus[i].az, elements[i].lat, elements[i].h); //{vxg, vyg, vzg};
		std::vector<double> cords = count_cords(elements[i].t, elements[i].vxg, elements[i].vyg, elements[i].vzg, elements[i].lon, elements[i].lat, elements[i].h); // {lat, lon, h}
		std::vector<double> rot_angles = count_rot_angles(elements[i].t, imus[i].wx, imus[i].wy, imus[i].wz, wgs[0], wgs[1], wgs[2], elements[i].pitch, elements[i].roll, elements[i].thdg); //
		elements.push_back(elem{ imus[i].t, cords[0], cords[1], cords[2], speeds[0], speeds[1], speeds[2], rot_angles[0], rot_angles[1], 360 - rot_angles[2] }); 
		if (i > 6794){
			std::cout << "File info       " << elements_check[i + 1] << "\n";
			std::cout << "Calculated info " << elements[i + 1] << "\n";
		}
	}
}
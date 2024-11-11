#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Windows.h>
#include "euler_krylov.h"
#include "ostreams.h"


int main() 
{
	std::ifstream nav_path("D:/практикум/Nav_1.txt");
	std::ifstream init_path("D:/практикум/init_1.txt");
	std::ifstream imu_path("D:/практикум/IMU_1.txt");
	
	elem init = {};
	std::vector<elem> elements_check = { init };
	std::vector<elem> elements = { init };
	std::vector<imu> imus;

	if (init_path.is_open()) {
		std::string rubbish;
		std::getline(init_path, rubbish);
		double dt, lat, lon, height, vx, vy, vz, pitch, roll, thdg;
		init_path >> dt >> lat >> lon >> height >> vx >> vy >> vz >> pitch >> roll >> thdg;
		init = elem{ dt, lat, lon, height, vx, vy, vz, pitch, roll, thdg };
	}

	if (nav_path.is_open()) {
		std::string rubbish;
		std::getline(nav_path, rubbish);
		double dt, lat, lon, height, vx, vy, vz, pitch, roll, thdg;
		while (nav_path >> dt >> lat >> lon >> height >> vx >> vy >> vz >> pitch >> roll >> thdg)
		{
			elements_check.push_back(elem{ dt, lat, lon, height, vx, vy, vz, pitch, roll, thdg });
		}
	}

	if (imu_path.is_open()) {
		std::string rubbish;
		std::getline(imu_path, rubbish);
		double t, ax, ay, az, wx, wy, wz;
		while (imu_path >> t >> ax >> ay >> az >> wx >> wy >> wz)
		{
			imus.push_back(imu{t, ax, ay, az, wx, wy, wz});
		}
	}

	for (int i = 0; i < imus.size(); i++) {
		std::vector<double> cords = count_cords(elements[i].t, elements[i].vxg, elements[i].vyg, elements[i].vzg, elements[i].lon, elements[i].lat, elements[i].height);
		std::vector<double> wgs = count_wgs(elements[i].vxg, elements[i].vyg, elements[i].vzg, cords[2], cords[0]);
		std::vector<double> speeds = count_speeds(elements[i].t, elements[i].vxg, elements[i].vyg, elements[i].vzg, imus[i].ax, imus[i].ay, imus[i].az, cords[0], cords[2]);
		std::vector<double> rot_angles = count_rot_angles(elements[i].t, imus[i].wx, imus[i].wy, imus[i].wz, wgs[0], wgs[1], wgs[2], 360 - elements[i].thdg, elements[i].pitch, elements[0].roll);
		elements.push_back(elem{ imus[i].t, cords[1], cords[0], cords[2], speeds[0], speeds[1], speeds[2], rot_angles[1], rot_angles[2], 360 - rot_angles[0] });

		std::cout << "File info       " << elements_check[i + 1] << "\n";
		std::cout << "Calculated info " << elements[i + 1] << "\n";
	}
}
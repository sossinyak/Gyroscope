#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <string>
#include <vector>
#include <Windows.h>
#include "euler_krylov.h"

inline std::ostream& operator << (std::ostream& o, const elem element) {
    o << std::fixed << std::setprecision(7) <<  element.t << " " << element.lat << " " << element.lon << " " << element.h << " " << element.vxg << " " 
	<< element.vyg << " " << element.vzg << " " << element.pitch << " " << element.roll << element.thdg;
    return o;
}

int main() 
{
	std::ifstream nav_path("Nav_1.txt");
	std::ifstream init_path("init_1.txt");
	std::ifstream imu_path("IMU_1.txt");
	std::ofstream fout("Conclusion.txt", std::ios_base::trunc);
	
	elem init = {};

	if (init_path.is_open()) {
		std::string rubbish;
		std::getline(init_path, rubbish);
		double dt, lat, lon, h, vx, vy, vz, pitch, roll, thdg;
		init_path >> dt >> lat >> lon >> h >> vx >> vy >> vz >> pitch >> roll >> thdg;
		init = elem{ dt, lat, lon, h, vx, vy, vz, pitch, roll, -thdg };
	}

	std::vector<elem> elements = { init };
	std::vector<imu> imus;
	double dt = elements[0].t;
	std::vector<elem> elements_check = { init };

	if (nav_path.is_open())
	{
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
			imus.push_back(imu{ t, ax, ay, az, wx, wy, wz });
		}
	}

	for (int i = 0; i < imus.size(); i++) {
		std::vector<double> corners = radians(elements[i].pitch, elements[i].roll, elements[i].thdg, elements[i].lat);
		std::vector<double> wgs = count_wgs(elements[i].vxg, elements[i].vzg, corners[3], elements[i].h); // {wxg, wyg, wzg}
		std::vector<double> rot_angles = count_rot_angles(dt, imus[i].wx, imus[i].wy, imus[i].wz, wgs[0], wgs[1], wgs[2], corners[0], corners[1], corners[2]); // {pitch, roll, thdg}
		corners = radians(rot_angles[0], rot_angles[1], rot_angles[2], elements[i].lat);
		std::vector<double> speeds = count_speeds(dt, elements[i].vxg, elements[i].vyg, elements[i].vzg, imus[i].ax, imus[i].ay, imus[i].az, corners[3], elements[i].h, corners[0], corners[1], corners[2]); //{vxg, vyg, vzg};
		std::vector<double> cords = count_cords(dt, speeds[0], speeds[1], speeds[2], elements[i].lat, elements[i].lon, elements[i].h); // {lat, lon, h}
		elements.push_back(elem{ imus[i].t, cords[0], cords[1], cords[2], speeds[0], speeds[1], speeds[2], rot_angles[0], rot_angles[1], rot_angles[2] }); 
		fout << elements[i + 1] << "\n";
	}
}
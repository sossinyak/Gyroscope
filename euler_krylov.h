#pragma once
#ifndef EULER_KRYLOV_H
#define EULER_KRYLOV_H

#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#define _USE_MATH_DEFINES
#define PI 3.14159265358979323846
#define R_EARTH 6400000.0
#define U_EARTH 7.2921158553 * std::pow(10, -5)
#define G_EARTH 9.8

struct elem
{
	double t, lat, lon, h, vxg, vyg, vzg, pitch, roll, thdg; // Nav
};


struct imu
{
	double t, ax, ay, az, wx, wy, wz; // IMU
};


std::vector<double> radians(double p, double r, double t, double l)
{
	return std::vector<double> {PI * p / 180, PI * r / 180, PI * t / 180, PI * l / 180}; 
};

// double radians(double degrees)
// {
// 	return PI *  degrees / 180; 
// };

double degrees(double radians)
{
	return 180 * radians / PI;
}

// абсолютная угловая скорость 
std::vector<double> count_wgs(double vxg, double vzg, double phi, double h)
{
	double R = R_EARTH + h;
	// проекции вектора абсолютной угловой скорости (3.6)
	double wxg = U_EARTH * std::cos(phi) + vzg / R;
	double wyg = U_EARTH * std::sin(phi) + vzg / R * std::tan(phi);
	double wzg = -vzg / R; 

	return std::vector<double> {wxg, wyg, wzg};
}

std::vector<double> count_rot_angles(double dt, double wx, double wy, double wz, double wxg, double wyg, double wzg, double theta, double gamma, double psi)
{
	// проекция вектора угловой скорости относительно географической системы координат (3.24)
	double wx_rel = wx - (wxg * std::cos(theta) * std::cos(psi) + wyg * std::sin(theta) - wzg * std::cos(theta) * std::sin(psi));
	double wy_rel = wy - (wxg * (-std::cos(gamma) * std::cos(psi) * std::sin(theta) + std::sin(gamma) * std::sin(psi)) + 
		wyg * std::cos(gamma) * std::cos(theta) + wzg * (std::cos(gamma) * std::sin(psi) * std::sin(theta) + std::sin(gamma) * std::cos(psi)));
	double wz_rel = wz - (wxg * (std::sin(gamma) * std::cos(psi) * std::sin(theta) - std::cos(gamma) * std::sin(psi)) - 
		wyg * std::sin(gamma) * std::cos(theta) + wzg * (-std::sin(gamma) * std::sin(psi) * std::sin(theta) + std::cos(gamma) * std::cos(psi)));

	// кинематические уравнения в углах Эйлера-Крылова (3.30)
	double psi_der = 1 / std::cos(theta) * (wy_rel * std::cos(gamma) - wz_rel * std::sin(gamma));
	double theta_der = wy_rel * std::sin(gamma) + wz_rel * std::cos(gamma);
	double gamma_der = wx_rel - std::tan(theta) * (wy_rel * std::cos(gamma) - wz_rel * std::sin(gamma));

	double psi_new = degrees(psi) + psi_der * dt;
	double theta_new = degrees(theta) + theta_der * dt;
	double gamma_new = degrees(gamma) + gamma_der * dt;

	return std::vector<double> {theta_new, gamma_new, psi_new};
}

// относительная скорость 
std::vector<double> count_speeds(double dt, double vxg, double vyg, double vzg, double ax, double ay, double az, 
								double phi, double h, double theta, double gamma, double psi)
{
	double R = R_EARTH + h;

	double nx = ax * std::cos(theta) * std::cos(psi) + ay * (-std::cos(gamma) * std::cos(psi) * std::sin(theta) + std::sin(gamma) * std::sin(psi)) 
	            + az * (std::sin(gamma) * std::cos(psi) * std::sin(theta) + std::cos(gamma) * std::sin(psi));
	double ny = ax * std::sin(theta) + ay * std::cos(gamma) * std::cos(theta) - az * std::sin(gamma) * std::cos(theta); 
	double nz = - ax * std::cos(theta) * std::sin(psi) + ay * (std::cos(gamma) * std::sin(psi) * std::sin(theta) + std::sin(gamma) * std::cos(psi)) 
	             + az * (-std::sin(gamma) * std::sin(psi) * std::sin(theta) + std::cos(gamma) * std::cos(psi));

	// компенсирующие составляющие ускорения (3.15)
	double akx = std::pow(vzg, 2) * std::tan(phi) / R + vxg * vyg / R + 2 * U_EARTH * vzg * std::sin(phi);
	double aky = - std::pow(vzg, 2)  / R - std::pow(vxg, 2)/ R - 2 * U_EARTH * vzg * std::cos(phi) + G_EARTH;
	double akz = vzg * vyg / R - vxg * vzg * std::tan(phi) / R + 2 * (vyg * U_EARTH * std::cos(phi) - U_EARTH * vxg * std::sin(phi));
	
	// составляющие относительной скорости (3.17)
	double vxg_new = vxg + (nx - akx) * dt;
	double vyg_new = vyg + (ny - aky) * dt;
	double vzg_new = vzg + (nz - akz) * dt;

	return std::vector<double>{vxg_new, vyg_new, vzg_new};
}

// высота, широта, долгота
std::vector<double> count_cords(double dt, double vxg, double vyg, double vzg, double phi, double lambda, double h)
{
	// координаты местоположения (3.18) 
	double h_new = h + vyg * dt;
	double R = R_EARTH + h_new;
	double phi_new = phi + vxg / R * dt;
	double lambda_new = lambda + vzg / (R * std::cos(phi))  * dt;
	return std::vector<double> {phi_new, lambda_new, h_new};
}

#endif
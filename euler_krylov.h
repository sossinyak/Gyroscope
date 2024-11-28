#pragma once
#ifndef EULER_KRYLOV_H
#define EULER_KRYLOV_H

#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#define _USE_MATH_DEFINES
#define PI 3.14159265358979323846
#define R 6380000.
#define U_EARTH 15.04 / 180. * PI / 3600.
#define G_EARTH 9.81

const double D2r = PI / 180.;
const double R2d = 180. / PI;


struct elem
{
	double t, lat, lon, h, vxg, vyg, vzg, pitch, roll, thdg;
};

struct imu
{
	double t, ax, ay, az, wx, wy, wz; 
};


std::vector<double> count_speeds(double dt, double vxg, double vyg, double vzg, double ax, double ay, double az, 
								double phi, double h, double theta, double gamma, double psi)
{
	double nx = ax * std::cos(theta) * std::cos(psi) + ay * (-std::cos(gamma) * std::cos(psi) * std::sin(theta) + std::sin(gamma) * std::sin(psi)) 
	            + az * (std::sin(gamma) * std::cos(psi) * std::sin(theta) + std::cos(gamma) * std::sin(psi));
	double ny = ax * std::sin(theta) + ay * std::cos(gamma) * std::cos(theta) - az * std::sin(gamma) * std::cos(theta); 
	double nz = - ax * std::cos(theta) * std::sin(psi) + ay * (std::cos(gamma) * std::sin(psi) * std::sin(theta) + std::sin(gamma) * std::cos(psi)) 
	             + az * (-std::sin(gamma) * std::sin(psi) * std::sin(theta) + std::cos(gamma) * std::cos(psi));

	double akx = std::pow(vzg, 2) * std::tan(phi) / R + vxg * vyg / R + 2 * U_EARTH * vzg * std::sin(phi);
	double aky = -std::pow(vzg, 2) / R - std::pow(vxg, 2)/ R - 2. * U_EARTH * vzg * std::cos(phi) + G_EARTH;
	double akz = vzg * vyg / R - vxg * vzg * std::tan(phi) / R + 2. * (vyg * U_EARTH * std::cos(phi) - U_EARTH * vxg * std::sin(phi));
	
	double vxg_new = vxg + (nx - akx) * dt;
	double vyg_new = vyg + (ny - aky) * dt;
	double vzg_new = vzg + (nz - akz) * dt;

	return std::vector<double>{vxg_new, vyg_new, vzg_new};
}

std::vector<double> count_cords(double dt, double vxg, double vyg, double vzg, double phi, double lambda, double h)
{
	double h_new = h + vyg * dt;
	double phi_new = phi * R2d + vxg / R * dt;
	double lambda_new = lambda * R2d + vzg / (R * std::cos(phi))  * dt;

	return std::vector<double> {phi_new, lambda_new, h_new};
}

std::vector<double> count_wgs(double vxg, double vzg, double phi)
{
	double wxg = U_EARTH * std::cos(phi) + vzg / R;
	double wyg = U_EARTH * std::sin(phi) + vzg / R * std::tan(phi);
	double wzg = -vxg / R; 

	return std::vector<double> {wxg, wyg, wzg};
}

std::vector<double> count_rot_angles(double dt, double wx, double wy, double wz, double wxg, double wyg, double wzg, double theta, double gamma, double psi)
{
	double wx_rel = wx - (wxg * std::cos(theta) * std::cos(psi) + wyg * std::sin(theta) - wzg * std::cos(theta) * std::sin(psi));
	double wy_rel = wy - (wxg * (-std::cos(gamma) * std::cos(psi) * std::sin(theta) + std::sin(gamma) * std::sin(psi)) + 
		wyg * std::cos(gamma) * std::cos(theta) + wzg * (std::cos(gamma) * std::sin(psi) * std::sin(theta) + std::sin(gamma) * std::cos(psi)));
	double wz_rel = wz - (wxg * (std::sin(gamma) * std::cos(psi) * std::sin(theta) - std::cos(gamma) * std::sin(psi)) - 
		wyg * std::sin(gamma) * std::cos(theta) + wzg * (-std::sin(gamma) * std::sin(psi) * std::sin(theta) + std::cos(gamma) * std::cos(psi)));

	double psi_der = 1 / std::cos(theta) * (wy_rel * std::cos(gamma) - wz_rel * std::sin(gamma));
	double theta_der = wy_rel * std::sin(gamma) + wz_rel * std::cos(gamma);
	double gamma_der = wx_rel - std::tan(theta) * (wy_rel * std::cos(gamma) - wz_rel * std::sin(gamma));

	double theta_new = theta * R2d + theta_der * dt;
	double gamma_new = gamma * R2d + gamma_der * dt;
	double psi_new = psi * R2d + psi_der * dt;

	return std::vector<double> {theta_new, gamma_new, psi_new};
}

#endif
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


double deg2rad(double degrees)
{
	return PI * degrees / 180;
};

// абсолютная угловая скорость 
std::vector<double> count_wgs(double vxg, double vzg, double phi0, double h)
{
	double phi = deg2rad(phi0), R = R_EARTH + h;
	// проекции вектора абсолютной угловой скорости (3.6)
	double wxg = U_EARTH * std::cos(phi) + vzg / R;
	double wyg = U_EARTH * std::sin(phi) + vzg / R * std::tan(phi);
	double wzg = -vxg / R;
	return std::vector<double> {wxg, wyg, wzg};
}

// относительная скорость (vxg, vyg, vzg)
std::vector<double> count_speeds(double t, double vxg, double vyg, double vzg, double ax, double ay, double az, double phi0, double h)
{
	double R = R_EARTH + h, phi = deg2rad(phi0);
	// проекции вектора кажещегося ускорения вершин трегранника на его оси (3.14)
	double nx = ax + std::pow(vzg, 2) * std::tan(phi)/R + vxg * vyg / R + 2 * U_EARTH * vzg * std::sin(phi);
	double ny = ay - std::pow(vzg, 2) / R - std::pow(vxg, 2) / R - 2 * U_EARTH * vzg * std::cos(phi) + G_EARTH;
	double nz = az + vzg * vyg / R - vxg * vzg * std::tan(phi) / R + 2 * (vyg * U_EARTH * std::cos(phi) - U_EARTH * vxg * std::sin(phi));

	// компенсирующие составляющие ускорения (3.15)
	double akx = std::pow(vzg, 2) * std::tan(phi) / R + vxg * vyg / R + 2 * U_EARTH * vzg * std::sin(phi);
	double aky = -std::pow(vzg, 2) / R - std::pow(vxg, 2) / R - 2 * U_EARTH * vzg * std::cos(phi) + G_EARTH;
	double akz = vzg * vyg / R - vxg * vzg * std::tan(phi) / R + 2 * (vyg * U_EARTH * std::cos(phi) - U_EARTH * vxg * std::sin(phi));
	
	// составляющие относительной скорости (3.17)
	double dt = 0.02;
	double vxg_new = vxg + (nx - akx) * dt;
	double vyg_new = vyg + (ny - aky) * dt;
	double vzg_new = vzg + (nz - akz) * dt;

	return std::vector<double>{vxg_new, vyg_new, vzg_new};
}

// высота, широта, долгота (phi, lambda, height)
std::vector<double> count_cords(double t, double vxg, double vyg, double vzg, double phi0, double lambda0, double h0)
{
	// координаты местоположения (3.18) 
	double dt = 0.02;
	double h = h0 + vyg * dt;
	double R = R_EARTH + h;
	double phi = phi0 + vxg / R * dt;
	double lambda = lambda0 + vzg / (R * std::cos(phi))  * dt;

	return std::vector<double> {phi, lambda, h};
}

std::vector<double> count_rot_angles(double t, double wx, double wy, double wz, double wxg, double wyg, double wzg, double teta0, double gamma0, double psi0)
{
	double gamma = deg2rad(gamma0), psi = deg2rad(psi0), teta = deg2rad(teta0);
	// проекция вектора угловой скорости относительно географической системы координат (3.24)
	double wx_rel = wx - (wxg * std::cos(teta) * std::cos(psi) + wyg * std::sin(teta) - wzg * std::cos(teta) * std::sin(psi));
	double wy_rel = wy - (wxg * (-std::cos(gamma) * std::cos(psi) * std::sin(teta) + std::sin(gamma) * std::sin(psi)) + 
		wyg * std::cos(gamma) * std::cos(teta) + wzg * (std::cos(gamma) * std::sin(psi) * std::sin(teta) + std::sin(gamma) * std::cos(psi)));
	double wz_rel = wz - (wxg * std::sin(gamma) * std::cos(psi) * std::sin(teta) - std::cos(gamma) * std::sin(psi) - 
		wyg * std::sin(gamma) * std::cos(teta) + wzg * (-std::sin(gamma) * std::sin(psi) * std::sin(teta) + std::cos(gamma) * std::cos(psi)));

	// кинематические уравнения в углах Эйлера-Крылова (3.30)
	double psi_der = 1 / std::cos(teta) * (wy_rel * std::cos(gamma) - wz_rel * std::sin(gamma)) + psi0;
	double teta_der = wy_rel * std::sin(gamma) + wz_rel * std::cos(gamma) + teta0;
	double gamma_der = wx_rel - std::tan(teta) * (wy_rel * std::cos(gamma) - wz_rel * std::sin(gamma)) + gamma0;
	return std::vector<double> {teta_der, gamma_der, psi_der};
}

#endif
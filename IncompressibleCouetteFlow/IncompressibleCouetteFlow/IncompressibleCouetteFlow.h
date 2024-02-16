#pragma once
#include <iostream>
#include "Data.h"

class IncompressibleCouetteFlow {
	double L;
	double D;
	double rho;
	double u_e;
	double mu;
	int Nx;
	int Ny;
	double dx;
	double dy;
	double dt;
	double alpha_p;
	double tol;
	void initialize(Data& data, std::vector<std::vector<double>>& u_star, std::vector<std::vector<double>>& v_star, std::vector<std::vector<double>>& p_star, std::vector<std::vector<double>>& rho_u_star, std::vector<std::vector<double>>& rho_v_star, std::vector<std::vector<double>>& p_prime_old, std::vector<std::vector<double>>& p_prime_new) {
		data.x_u.resize(Nx + 2);
		for (int i = 0; i < Nx + 2; i++) {
			data.x_u[i] = (i - 0.5) * dx;
		}
		data.y_u.resize(Ny + 1);
		for (int i = 0; i < Ny + 1; i++) {
			data.y_u[i] = i * dy;
		}
		data.x_v.resize(Nx + 3);
		for (int i = 0; i < Nx + 3; i++) {
			data.x_v[i] = (i - 1) * dx;
		}
		data.y_v.resize(Ny + 2);
		for (int i = 0; i < Ny + 2; i++) {
			data.y_v[i] = (i - 0.5) * dy;
		}
		data.x_p.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			data.x_p[i] = i * dx;
		}
		data.y_p.resize(Ny + 1);
		for (int i = 0; i < Ny + 1; i++) {
			data.y_p[i] = i * dy;
		}
		u_star.resize(Nx + 2);
		for (int i = 0; i < Nx + 2; i++) {
			u_star[i].resize(Ny + 1);
		}
		for (int i = 0; i < Nx + 2; i++) {
			u_star[i][Ny] = u_e;
		}
		v_star.resize(Nx + 3);
		for (int i = 0; i < Nx + 3; i++) {
			v_star[i].resize(Ny + 2);
		}
		v_star[14][4] = u_e / 2;
		p_star.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			p_star[i].resize(Ny + 1);
		}
		rho_u_star.resize(Nx + 2);
		for (int i = 0; i < Nx + 2; i++) {
			rho_u_star[i].resize(Ny + 1);
		}
		rho_v_star.resize(Nx + 3);
		for (int i = 0; i < Nx + 3; i++) {
			rho_v_star[i].resize(Ny + 2);
		}
		p_prime_old.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			p_prime_old[i].resize(Ny + 1);
		}
		p_prime_new.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			p_prime_new[i].resize(Ny + 1);
		}
		data.t.push_back(0);
		data.u.push_back(u_star);
		data.v.push_back(v_star);
		data.p.push_back(p_star);
	}
	double residual(std::vector<std::vector<double>>& x1, std::vector<std::vector<double>>& x2) {
		double res = 0;
		for (int i = 0; i < x1.size(); i++) {
			for (int j = 0; j < x1[i].size(); j++) {
				res = res < abs(x1[i][j] - x2[i][j]) ? abs(x1[i][j] - x2[i][j]) : res;
			}
		}
		return res;
	}
	double residual(std::vector<std::vector<double>>& x) {
		double res = 0;
		for (int i = 0; i < x.size(); i++) {
			for (int j = 0; j < x[i].size(); j++) {
				res = res < abs(x[i][j]) ? abs(x[i][j]) : res;
			}
		}
		return res;
	}
	double residual(Data& data) {
		double res_u = residual(data.u[data.u.size() - 1], data.u[data.u.size() - 2]);
		double res_v = residual(data.v[data.v.size() - 1], data.v[data.v.size() - 2]);
		double res_p = residual(data.p[data.p.size() - 1], data.p[data.p.size() - 2]);
		double res = std::max(std::max(res_u, res_v), res_p);
		return res;
	}
public:
	IncompressibleCouetteFlow() :L(0.5), D(0.01), rho(0.002377), u_e(1), mu(rho * u_e * D / 63.6), Nx(20), Ny(10), dx(L / Nx), dy(D / Ny), dt(0.001), alpha_p(0.1), tol(1e-6) {}
	void pressureCorrection(Data& data) {
		double res;
		double a = 2 * (dt / pow(dx, 2) + dt / pow(dy, 2));
		double b = -dt / pow(dx, 2);
		double c = -dt / pow(dy, 2);
		std::vector<std::vector<double>> u_star, v_star, p_star, rho_u_star, rho_v_star, p_prime_old, p_prime_new;
		initialize(data, u_star, v_star, p_star, rho_u_star, rho_v_star, p_prime_old, p_prime_new);
		do {
			do {
				for (int i = 1; i < Nx + 1; i++) {
					for (int j = 1; j < Ny; j++) {
						double v_bar = (v_star[i][j + 1] + v_star[i + 1][j + 1]) / 2;
						double v = (v_star[i][j] + v_star[i + 1][j]) / 2;
						double A_star = -((rho * pow(u_star[i + 1][j], 2) - rho * pow(u_star[i - 1][j], 2)) / (2 * dx) + (rho * u_star[i][j + 1] * v_bar - rho * u_star[i][j - 1] * v) / (2 * dy)) + mu * ((u_star[i + 1][j] - 2 * u_star[i][j] + u_star[i - 1][j]) / pow(dx, 2) + (u_star[i][j + 1] - 2 * u_star[i][j] + u_star[i][j - 1]) / pow(dy, 2));
						rho_u_star[i][j] = rho * u_star[i][j] + A_star * dt - dt / dx * (p_star[i][j] - p_star[i - 1][j]);
					}
				}
				for (int i = 1; i < Nx + 2; i++) {
					for (int j = 1; j < Ny + 1; j++) {
						double u_bar = (u_star[i][j - 1] + u_star[i][j]) / 2;
						double u = (u_star[i - 1][j - 1] + u_star[i - 1][j]) / 2;
						double B_star = -((rho * v_star[i + 1][j] * u_bar - rho * v_star[i - 1][j] * u) / (2 * dx) + (rho * pow(v_star[i][j + 1], 2) - rho * pow(v_star[i][j - 1], 2)) / (2 * dy)) + mu * ((v_star[i + 1][j] - 2 * v_star[i][j] + v_star[i - 1][j]) / pow(dx, 2) + (v_star[i][j + 1] - 2 * v_star[i][j] + v_star[i][j - 1]) / pow(dy, 2));
						rho_v_star[i][j] = rho * v_star[i][j] + B_star * dt - dt / dy * (p_star[i - 1][j] - p_star[i - 1][j - 1]);
					}
				}
				for (int i = 0; i < Nx + 2; i++) {
					rho_u_star[i][Ny] = rho * u_e;
				}
				for (int i = 0; i < Ny + 1; i++) {
					rho_u_star[0][i] = rho_u_star[1][i];
				}
				for (int i = 0; i < Ny + 1; i++) {
					rho_u_star[Nx + 1][i] = rho_u_star[Nx][i];
				}
				for (int i = 0; i < Ny + 2; i++) {
					rho_v_star[Nx + 2][i] = rho_v_star[Nx + 1][i];
				}
				for (int i = 0; i < Nx + 1; i++) {
					for (int j = 0; j < Ny + 1; j++) {
						p_prime_old[i][j] = 0;
					}
				}
				do {
					for (int i = 1; i < Nx; i++) {
						for (int j = 1; j < Ny; j++) {
							double d = (rho_u_star[i + 1][j] - rho_u_star[i][j]) / dx + (rho_v_star[i + 1][j + 1] - rho_v_star[i + 1][j]) / dy;
							p_prime_new[i][j] = -(b * p_prime_old[i + 1][j] + b * p_prime_old[i - 1][j] + c * p_prime_old[i][j + 1] + c * p_prime_old[i][j - 1] + d) / a;
						}
					}
					res = residual(p_prime_old, p_prime_new);
					p_prime_old = p_prime_new;
				} while (res > tol);
				for (int i = 0; i < Nx + 1; i++) {
					for (int j = 0; j < Ny + 1; j++) {
						p_star[i][j] += alpha_p * p_prime_new[i][j];
					}
				}
				res = residual(p_prime_new);
			} while (res > tol);
			for (int i = 1; i < Nx + 1; i++) {
				for (int j = 1; j < Ny; j++) {
					double v_bar = (v_star[i][j + 1] + v_star[i + 1][j + 1]) / 2;
					double v = (v_star[i][j] + v_star[i + 1][j]) / 2;
					double A_star = -((rho * pow(u_star[i + 1][j], 2) - rho * pow(u_star[i - 1][j], 2)) / (2 * dx) + (rho * u_star[i][j + 1] * v_bar - rho * u_star[i][j - 1] * v) / (2 * dy)) + mu * ((u_star[i + 1][j] - 2 * u_star[i][j] + u_star[i - 1][j]) / pow(dx, 2) + (u_star[i][j + 1] - 2 * u_star[i][j] + u_star[i][j - 1]) / pow(dy, 2));
					rho_u_star[i][j] = rho * u_star[i][j] + A_star * dt - dt / dx * (p_star[i][j] - p_star[i - 1][j]);
				}
			}
			for (int i = 1; i < Nx + 2; i++) {
				for (int j = 1; j < Ny + 1; j++) {
					double u_bar = (u_star[i][j - 1] + u_star[i][j]) / 2;
					double u = (u_star[i - 1][j - 1] + u_star[i - 1][j]) / 2;
					double B_star = -((rho * v_star[i + 1][j] * u_bar - rho * v_star[i - 1][j] * u) / (2 * dx) + (rho * pow(v_star[i][j + 1], 2) - rho * pow(v_star[i][j - 1], 2)) / (2 * dy)) + mu * ((v_star[i + 1][j] - 2 * v_star[i][j] + v_star[i - 1][j]) / pow(dx, 2) + (v_star[i][j + 1] - 2 * v_star[i][j] + v_star[i][j - 1]) / pow(dy, 2));
					rho_v_star[i][j] = rho * v_star[i][j] + B_star * dt - dt / dy * (p_star[i - 1][j] - p_star[i - 1][j - 1]);
				}
			}
			for (int i = 0; i < Nx + 2; i++) {
				rho_u_star[i][Ny] = rho * u_e;
			}
			for (int i = 0; i < Ny + 1; i++) {
				rho_u_star[0][i] = rho_u_star[1][i];
			}
			for (int i = 0; i < Ny + 1; i++) {
				rho_u_star[Nx + 1][i] = rho_u_star[Nx][i];
			}
			for (int i = 0; i < Ny + 2; i++) {
				rho_v_star[Nx + 2][i] = rho_v_star[Nx + 1][i];
			}
			for (int i = 0; i < Nx + 2; i++) {
				for (int j = 0; j < Ny + 1; j++) {
					u_star[i][j] = rho_u_star[i][j] / rho;
				}
			}
			for (int i = 0; i < Nx + 3; i++) {
				for (int j = 0; j < Ny + 2; j++) {
					v_star[i][j] = rho_v_star[i][j] / rho;
				}
			}
			data.t.push_back(data.t[data.t.size() - 1] + dt);
			data.u.push_back(u_star);
			data.v.push_back(v_star);
			data.p.push_back(p_star);
			res = residual(data);
			printf("time=%f res=%f\n", data.t[data.t.size() - 1], res);
		} while (res > tol);
	}
};

/*
 * ModelParam.h
 *
 *  Created on: Mar 13, 2018
 *      Author: chuckjia
 */

#ifndef INCLUDE_MODELPARAM_H_
#define INCLUDE_MODELPARAM_H_

class ModelParam {
public:

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Constructors
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	ModelParam() {
		x0_ = xf_ = 0;
		pA_ = 0;
		pB_ = 0;
		pB_x_ = 0;
		Nx_ = Np_ = 0;
		Dx_ = 0;
		num_cell_x_ = num_cell_p_ = num_gridpt_x_ = num_gridpt_p_ = 0;
		init_T_ = init_q_ = init_u_ = 0;
	}

	ModelParam(double _x0, double _xf, double _pA, double (*_pB)(double), double (*_pB_x)(double), int _Nx, int _Np,
			double (*_T)(double, double, double), double (*_q)(double, double, double), double (*_u)(double, double, double)) {
		x0_ = _x0;
		xf_ = _xf;
		pA_ = _pA;
		pB_ = _pB;
		pB_x_ = _pB_x;

		Nx_ = _Nx;
		Np_ = _Np;

		Dx_ = (_xf - _x0) / _Nx;

		num_cell_x_ = Nx_ + 2;
		num_cell_p_ = Np_ + 2;

		num_gridpt_x_ = num_cell_x_ + 1;
		num_gridpt_p_ = num_cell_p_ + 1;

		init_T_ = _T;
		init_q_ = _q;
		init_u_ = _u;
	}

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * I/O
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	void write_to_file() {
		FILE *f = fopen("Output/param.csv", "wb");
		fprintf(f, "%1.20e,%1.20e,%1.20e,%1.20e,%1.20e,%1.20e,%d,%d",
				x0_, xf_, pA_, (*pB_)(x0_), (*pB_)(0.5 * (x0_ + xf_)), (*pB_)(xf_), Nx_, Np_);
		fclose(f);
	}

	/* ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
	 * Informational
	 * ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== */

	inline double x0() const { return x0_; }

	inline double xf() const { return xf_; }

	inline double pA() const { return pA_; }

	inline double (*pB())(double) { return pB_; }

	inline double (*pB_x())(double) { return pB_x_; }

	inline int Nx() const { return Nx_; }

	inline int Np() const { return Np_; }

	inline double Dx() const { return Dx_; }

	inline int numCellX() const { return num_cell_x_; }

	inline int numCellP() const { return num_cell_p_; }

	inline int numGridPtX() const { return num_gridpt_x_; }

	inline int numGridPtP() const { return num_gridpt_p_; }

protected:

	double x0_, xf_, pA_, Dx_;
	double (*pB_)(double), (*pB_x_)(double);

	int Nx_, Np_, num_cell_x_, num_cell_p_, num_gridpt_x_, num_gridpt_p_;

	// Initial conditions. Assume each function has (x, p, t) as input parameters
	double (*init_T_)(double, double, double), (*init_q_)(double, double, double), (*init_u_)(double, double, double);
};

#endif /* INCLUDE_MODELPARAM_H_ */

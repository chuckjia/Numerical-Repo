/*
 * MeshBuild.h
 *
 *  Created on: Mar 12, 2018
 *      Author: chuckjia
 */

#ifndef SRC_MESHBUILD_H_
#define SRC_MESHBUILD_H_

#include "../include/MeshGrid.h"
#include "../include/ModelParam.h"
#include "../include/Vector1D.h"

class MeshBuild {
protected:
	// Calculate the coordinates of the intersection of 2 lines. For each line, the coordinates of two points are given.
	// Assume (x1, p1) and (x3, p3) are on line 1 with slope k1; (x2, p2) and (x4, p4) are on line 2 with slope k2. See details in notes
	// Calculates the coordinates of the intersection point (x, y) = (ans[0], ans[1])
	static void calcLineIntersection(double ans[2], double x1, double x2, double x3, double x4, double p1, double p2, double p3, double p4) {
		double k1 = (p3 - p1) / (x3 - x1), k2 = (p4 - p2) / (x4 - x2);

		double x = (k1 * x1 - p1 - k2 * x2 + p2) / (k1 - k2);
		ans[0] = x;
		ans[1] = k1 * (x - x1) + p1;
	}

	// Assume that (x1, p1), (x2, p2), (x3, p3), and (x4, p4) are arranged in counter-clockwise order
	static void calcTrapezoidCenter(double center[2], double x1, double x2, double x3, double x4, double p1, double p2, double p3, double p4) {
		double x124 = (x1 + x2 + x4) / 3, x134 = (x1 + x3 + x4) / 3, x123 = (x1 + x2 + x3) / 3, x234 = (x2 + x3 + x4) / 3;
		double p124 = (p1 + p2 + p4) / 3, p134 = (p1 + p3 + p4) / 3, p123 = (p1 + p2 + p3) / 3, p234 = (p2 + p3 + p4) / 3;

		calcLineIntersection(center, x124, x123, x234, x134, p124, p123, p234, p134);
	}

public:

	static bool sameSize(Matrix2D &a, Matrix2D &b) {
		return a.nrow() == b.nrow() && a.ncol() == b.ncol();
	}

	// Given model parameters, build the grid points
	static void buildMeshGrid(ModelParam &model, MeshGrid &grid_x, MeshGrid &grid_p, Vector1D &cell_left_Dp) {
		int num_gridpt_x = model.numGridPtX(), num_gridpt_p = model.numGridPtP();
		assert(num_gridpt_x == grid_x.nrow() && num_gridpt_p == grid_x.ncol() &&
				num_gridpt_x == grid_p.nrow() && num_gridpt_p == grid_p.ncol());
		assert(num_gridpt_x = cell_left_Dp.len());

		int Nx = model.Nx(), Np = model.Np();
		double x0 = model.x0(), xf = model.xf(), pA = model.pA();
		double (*pB)(double) = model.pB();

		double Dx = (xf - x0) / Nx;

		for (int i = 0; i < num_gridpt_x; ++i) {
			double x = x0 + (i - 1) * Dx,
					Dp = ((*pB)(x) - pA) / Np;
			cell_left_Dp[i] = Dp;
			for (int j = 0; j < num_gridpt_p; ++j) {
				double p = pA + (j - 1) * Dp;
				grid_x[i][j] = x;
				grid_p[i][j] = p;
			}
		}
	}

	// Given grid points, build the cell center coordinates
	static void calcCellCenters(MeshGrid &grid_x, MeshGrid &grid_p, Matrix2D &cell_centers_x, Matrix2D &cell_centers_p) {
		assert(sameSize(grid_x, grid_p));
		double num_cell_x = grid_x.nrow() - 1, num_cell_p = grid_x.ncol() - 1;
		assert(cell_centers_x.nrow() == num_cell_x && cell_centers_p.nrow() == num_cell_x &&
				cell_centers_x.ncol() == num_cell_p && cell_centers_p.ncol() == num_cell_p);

		for (int i = 0; i < num_cell_x; ++i)
			for (int j = 0; j < num_cell_p; ++j) {
				double center[2];
				double x1 = grid_x.bottLeft(i, j), x2 = grid_x.bottRight(i, j), x3 = grid_x.topRight(i, j), x4 = grid_x.topLeft(i, j);
				double p1 = grid_p.bottLeft(i, j), p2 = grid_p.bottRight(i, j), p3 = grid_p.topRight(i, j), p4 = grid_p.topLeft(i, j);
				calcTrapezoidCenter(center, x1, x2, x3, x4, p1, p2, p3, p4);
				cell_centers_x[i][j] = center[0];
				cell_centers_p[i][j] = center[1];
			}
	}

	static void buildMesh(ModelParam &model, MeshGrid &grid_x, MeshGrid &grid_p,
			Matrix2D &cell_centers_x, Matrix2D &cell_centers_p, Vector1D &cell_left_Dp) {
		buildMeshGrid(model, grid_x, grid_p, cell_left_Dp);
		calcCellCenters(grid_x, grid_p, cell_centers_x, cell_centers_p);

		grid_x.printToFile("grid_x");
		grid_p.printToFile("grid_p");
		cell_centers_x.printToFile("cell_centers_x");
		cell_centers_p.printToFile("cell_centers_p");
	}

	static void calcCellVol(ModelParam &model, Vector1D &cell_left_Dp, Vector1D &cell_vol) {
		int num_cell_x = model.numCellX();
		assert(num_cell_x == cell_vol.len() && num_cell_x + 1 == cell_left_Dp.len());
		double Dx = model.Dx();
		for (int i = 0; i < num_cell_x; ++i)
			cell_vol[i] = 0.5 * (cell_left_Dp[i] + cell_left_Dp[i + 1]) * Dx;
	}

	static Matrix2D calcCellBottLen(MeshGrid &grid_x, MeshGrid &grid_p) {
		assert(sameSize(grid_x, grid_p));
		int num_cell_x = grid_x.nrow() - 1, num_cell_p = grid_x.ncol() - 1;

		Matrix2D cell_bott_len(num_cell_x, num_cell_p);
		for (int i = 0; i < num_cell_x; ++i)
			for (int j = 0; j < num_cell_p; ++j) {
				double x_diff = grid_x.bottLeft(i, j) - grid_x.bottRight(i, j),
						p_diff = grid_p.bottLeft(i, j) - grid_p.bottRight(i, j);
				cell_bott_len[i][j] = sqrt(x_diff * x_diff + p_diff * p_diff);
			}
		return cell_bott_len;
	}

	static void calcCellBottNorm(MeshGrid &grid_x, MeshGrid &grid_p, Matrix2D &cell_bott_norm_x, Matrix2D &cell_bott_norm_p) {
		assert(sameSize(grid_x, grid_p) && sameSize(cell_bott_norm_x, cell_bott_norm_p));
		int num_cell_x = grid_x.nrow() - 1, num_cell_p = grid_x.ncol() - 1;
		assert(num_cell_x == cell_bott_norm_x.nrow() && num_cell_p == cell_bott_norm_x.ncol());

		Matrix2D cell_bott_norm(num_cell_x, num_cell_p);

		for (int i = 0; i < num_cell_x; ++i) {

		}
	}

	static void fillcache_mesh(MeshGrid &grid_x, MeshGrid &grid_p, Matrix2D &cell_vol, Matrix2D &cell_bott_len) {
		cell_bott_len = calcCellBottLen(grid_x, grid_p);
	}
};


#endif /* SRC_MESHBUILD_H_ */

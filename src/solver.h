#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "grid.h"

class SmokeSolver {
public:
    SmokeSolver(int grid_num_x, int grid_num_y, int grid_num_z);
    void run_dt(int);
	void dump_density(int i);

private:
    void init_source();
    void advect_velocity();
    void advect_density_temp();
    void apply_force();
    void project_pressure();
	void apply_boundary();

    Vec3f trace(const Vec3f& orig, const Vec3f& u, Float t);

    // pressure matrix solve
    Eigen::VectorXd b,x;
    Eigen::SparseMatrix<Float, Eigen::RowMajor> A;
	std::vector<Eigen::Triplet<Float> > triple;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<Float>, Eigen::Lower | Eigen::Upper> CG;

    Float dt;
    Float total_time;

    // MAC grid data structure.
    MACGrid grid;
};

#endif
#include "solver.h"
#include "common.h"
#include <iostream>
#include <fstream>
#include <random>

SmokeSolver::SmokeSolver(int grid_num_x, int grid_num_y, int grid_num_z)
	:grid(grid_num_x, grid_num_y, grid_num_z) {
	dt = 0.02;
	total_time = 0.0;
	int size = grid_num_x * grid_num_y * grid_num_z;
	A = Eigen::SparseMatrix<Float, Eigen::RowMajor>(size, size);
	b = Eigen::VectorXd(size);
	x = Eigen::VectorXd(size);

	triple.reserve(7 * size);
	CG.setTolerance(eps);
	init_source();              
	dump_density(0);
}


void SmokeSolver::run_dt(int step) {
	project_pressure();
	advect_velocity();
	advect_density_temp();
    apply_force();               
    total_time += dt;	
	init_source();
}

Vec3f SmokeSolver::trace(const Vec3f& orig, const Vec3f& u, Float t) {
   // return orig - u * t;
	Vec3f dir = u * t;
	Vec3f target = orig - dir;
	while (grid.occupy.get((int)target.x, (int)target.y, (int)target.z) != 0) {
		dir = dir * 0.9;
		target = orig - dir;
	}
	return target;
}

void SmokeSolver::init_source() {
	// init temp to ambient.
	for (int i = 0; i < grid.grid_num_x; i++) {
		for (int j = 0; j < grid.grid_num_y; j++) {
			for (int k = 0; k < grid.grid_num_z; k++) {
				if (grid.occupy.get(i, j, k) == 0) {
					grid.temp0.get(i, j, k) = T_ambient;
				}
			}
		}
	}
	// single source
	for (int k = 12; k < 20; ++k) {
		for (int i = 12; i < 20; ++i) {
			for (int j = 58; j < 61; ++j) {
				grid.density0.get(i, j, k) = 1.0;
				grid.temp0.get(i, j, k) = T_ambient + 100;
			}
		}
	}
	//// double source
	//for (int k = 12; k < 20; ++k) {
	//	for (int i = 12; i < 20; ++i) {
	//		for (int j = 58; j < 61; ++j) {
	//			grid.density0.get(i, j, k) = 1.0;
	//			grid.temp0.get(i, j, k) = T_ambient + 100;
	//		}
	//		for (int j = 3; j < 6; ++j) {
	//			grid.density0.get(i, j, k) = 1.0;
	//			grid.F[1].get(i, j, k) = 2000;
	//			grid.temp0.get(i, j, k) = T_ambient + 100;
	//		}
	//	}
	//}

	//// rising 
	//for (int k = 26; k < 30; ++k) {
	//	for (int i = 2; i < 8; ++i) {
	//		for (int j = grid.grid_num_y - 3 - 3; j < grid.grid_num_y - 3; ++j) {
	//			grid.density0.get(i, j, k) = 1.0;
	//			grid.temp0.get(i, j, k) = T_ambient + 100;
	//			grid.F[2].get(i, j, k) = -2000;
	//		}
	//	}
	//}

	// apply force in initial state
	for (int i = 0; i< grid.grid_num_x; i++) {
		for (int j = 0; j< grid.grid_num_y; j++) {
			for (int k = 0; k< grid.grid_num_z; k++) {
				if (grid.occupy.get(i, j, k) != 0)
					continue;

				if (i < grid.grid_num_x - 1) {
					grid.v0[0].get(i + 1, j, k) += (grid.F[0].get(i, j, k) + grid.F[0].get(i + 1, j, k)) * 0.5 * dt;
				}
				if (j < grid.grid_num_y - 1) {
					grid.v0[1].get(i, j + 1, k) += (grid.F[1].get(i, j, k) + grid.F[1].get(i, j + 1, k)) * 0.5 * dt;
				}
				if (k < grid.grid_num_z - 1) {
					grid.v0[2].get(i, j, k + 1) += (grid.F[2].get(i, j, k) + grid.F[2].get(i, j, k + 1)) * 0.5 * dt;
				}
			}
		}
	}
	apply_boundary();
}


void SmokeSolver::apply_boundary()
{
	int center_x = grid.grid_num_x / 2, center_y = grid.grid_num_y / 2, center_z = grid.grid_num_z / 2;
	int rad = 8;
	/*for (int i = center_x - 10; i < center_x + 10; ++i) {
		for (int j = center_y - 1; j < center_y + 1; ++j) {
			for (int k = 20; k < grid.grid_num_z; ++k) {*/
	for (int i = 0; i < grid.grid_num_x; i++) {
		for (int j = 0; j < grid.grid_num_y; j++) {
			for (int k = 0; k < grid.grid_num_z; k++) {
				Vec3f vec(i - center_x, j - center_y, k - center_z);
				if (vec.length() < rad) {
					grid.occupy.get(i, j, k) = 1.0;
					grid.density0.get(i, j, k) = 0;
					grid.v0[0].get(i, j, k) = 0;
					grid.v0[1].get(i, j, k) = 0;
					grid.v0[2].get(i, j, k) = 0;
					grid.v1[0].get(i, j, k) = 0;
					grid.v1[1].get(i, j, k) = 0;
					grid.v1[2].get(i, j, k) = 0;
					grid.v0[0].get(i + 1, j, k) = 0;
					grid.v0[1].get(i, j + 1, k) = 0;
					grid.v0[2].get(i, j, k + 1) = 0;
					grid.v1[0].get(i + 1, j, k) = 0;
					grid.v1[1].get(i, j + 1, k) = 0;
					grid.v1[2].get(i, j, k + 1) = 0;
				}
			}
		}
	}
}


void SmokeSolver::apply_force() {
    // buoyant
	for (int i = 0; i < grid.grid_num_x; i++) {
			for (int j = 0; j < grid.grid_num_y; j++) {
				for (int k = 0; k < grid.grid_num_z; k++) {
					if (grid.occupy.get(i, j, k) != 0)
						continue;

					grid.F[0].get(i, j, k) = 0.0;
					grid.F[1].get(i, j, k) = alpha * grid.density0.get(i, j, k) - beta * (grid.temp0.get(i, j, k) - T_ambient);
					grid.F[2].get(i, j, k) = 0.0;
				}
			}
		}
	
    // vorticity
    // 1. compute center velocities.
    for(int i=0;i< grid.grid_num_x; i++) {
        for(int j=0;j< grid.grid_num_y; j++) {
            for(int k=0;k< grid.grid_num_z; k++) {
				if (grid.occupy.get(i, j, k) != 0) continue;
                grid.avg_v[0].get(i,j,k) = (grid.v0[0].get(i,j,k) + grid.v0[0].get(i+1,j,k) )/2.0;
                grid.avg_v[1].get(i,j,k) = (grid.v0[1].get(i,j,k) + grid.v0[1].get(i,j+1,k) )/2.0;
                grid.avg_v[2].get(i,j,k) = (grid.v0[2].get(i,j,k) + grid.v0[2].get(i,j,k+1) )/2.0;              
            }
        }
    }
	

    // 2. compute vorticity.
    for(int i=0;i< grid.grid_num_x; i++) {
        for(int j=0;j< grid.grid_num_y; j++) {
            for(int k=0;k< grid.grid_num_z; k++) {
                if(i == 0 || j == 0 || k == 0 || i == grid.grid_num_x - 1 
                    || j == grid.grid_num_y - 1  || k == grid.grid_num_z - 1 
                    ||  grid.occupy.get(i,j,k) != 0) {
                    continue;
                }

                Float vort_x,vort_y, vort_z;
                vort_x = (grid.avg_v[2].get(i,j+1,k) - grid.avg_v[2].get(i,j-1,k) - grid.avg_v[1].get(i,j,k+1) + grid.avg_v[1].get(i,j,k-1)) / 2.0;
                vort_y = (grid.avg_v[0].get(i,j,k+1) - grid.avg_v[0].get(i,j,k-1) - grid.avg_v[2].get(i+1,j,k) + grid.avg_v[2].get(i-1,j,k)) / 2.0;
                vort_z = (grid.avg_v[1].get(i+1,j,k) - grid.avg_v[1].get(i-1,j,k) - grid.avg_v[0].get(i,j+1,k) + grid.avg_v[0].get(i,j-1,k)) / 2.0;           
                
                grid.vorticity_norm.get(i,j,k) = Vec3f(vort_x, vort_y, vort_z).length();
                grid.vorticity[0].get(i,j,k) = vort_x;
                grid.vorticity[1].get(i,j,k) = vort_y;
                grid.vorticity[2].get(i,j,k) = vort_z;
            }
        }
    }

    // 3. compute vorticity force
    for(int i=0;i< grid.grid_num_x; i++) {
        for(int j=0;j< grid.grid_num_y; j++) {
            for(int k=0;k< grid.grid_num_z; k++) {
				if (i == 0 || j == 0 || k == 0 || i == grid.grid_num_x - 1
					|| j == grid.grid_num_y - 1 || k == grid.grid_num_z - 1
					|| grid.occupy.get(i, j, k) != 0) {
					grid.vorticity[0].get(i, j, k) = 0;
					grid.vorticity[1].get(i, j, k) = 0;
					grid.vorticity[2].get(i, j, k) = 0;
					continue;
				}

                Float Nx = grid.vorticity_norm.get(i+1,j,k) - grid.vorticity_norm.get(i-1,j,k);
                Float Ny = grid.vorticity_norm.get(i,j+1,k) - grid.vorticity_norm.get(i,j-1,k);
                Float Nz = grid.vorticity_norm.get(i,j,k+1) - grid.vorticity_norm.get(i,j,k-1);

                Vec3f N(Nx, Ny, Nz);
                Float length = N.length();
                if(length != 0)
					N = N * (1.0/length);
				else N = Vec3f(0, 0, 0);
                
                Vec3f vort(grid.vorticity[0].get(i,j,k), grid.vorticity[1].get(i,j,k), grid.vorticity[2].get(i,j,k));
				Vec3f f_vort = cross(vort, N) * epsilon;

                grid.F[0].get(i,j,k) += f_vort.x;
				grid.F[1].get(i, j, k) += f_vort.y;
                grid.F[2].get(i,j,k) += f_vort.z;
            }
        }
    }
	
    // apply force
     for(int i=0;i< grid.grid_num_x; i++) {
        for(int j=0;j< grid.grid_num_y; j++) {
            for(int k=0;k< grid.grid_num_z; k++) {
				if (grid.occupy.get(i, j, k) != 0 && grid.occupy.get(i+1, j, k) != 0 && grid.occupy.get(i, j + 1, k) && grid.occupy.get(i, j, k + 1))
					continue;

                if(i < grid.grid_num_x - 1) {
                    grid.v0[0].get(i+1,j,k) += (grid.F[0].get(i,j,k) + grid.F[0].get(i+1,j,k)) * 0.5 * dt;
                } 
                if(j < grid.grid_num_y - 1) {
                    grid.v0[1].get(i,j+1,k) += (grid.F[1].get(i,j,k) + grid.F[1].get(i,j+1,k)) * 0.5 * dt;
                } 
                if(k < grid.grid_num_z - 1) {
                    grid.v0[2].get(i,j,k+1) += (grid.F[2].get(i,j,k) + grid.F[2].get(i,j,k+1)) * 0.5 * dt;
                } 
            }
        }
     }
	
}
void SmokeSolver::advect_velocity() {
     //swap v0, v1
    std::copy(grid.v0[0].data.begin(), grid.v0[0].data.end(), grid.v1[0].data.begin());
    std::copy(grid.v0[1].data.begin(), grid.v0[1].data.end(), grid.v1[1].data.begin());
    std::copy(grid.v0[2].data.begin(), grid.v0[2].data.end(), grid.v1[2].data.begin());
	
    // x 
    for(int i=0;i< grid.grid_num_x; i++) {
        for(int j=0;j< grid.grid_num_y; j++) {
            for(int k=0;k< grid.grid_num_z; k++) {
                // ignore solid object
               if(grid.occupy.get(i,j,k) != 0) {
                    continue;
                } 
			   
			   if (i + 1 < grid.grid_num_x && grid.occupy.get(i + 1, j, k) == 0) {
				   Vec3f orig(i + 1.0, 0.5 + j, 0.5 + k);
				   Float vx = get_face_value(grid.v1[0], orig, AXIS::X);
				   Float vy = get_face_value(grid.v1[1], orig, AXIS::Y);
				   Float vz = get_face_value(grid.v1[2], orig, AXIS::Z);
				   Vec3f p = trace(orig, Vec3f(vx,vy,vz), dt);
				   grid.v0[0].get(i+1,j,k) = get_face_value(grid.v1[0], p, AXIS::X);
				}
			   if (j + 1 < grid.grid_num_y && grid.occupy.get(i, j + 1, k) == 0) {
				   Vec3f orig(i + 0.5, j + 1.0, 0.5 + k);
				   Float vx = get_face_value(grid.v1[0], orig, AXIS::X);
				   Float vy = get_face_value(grid.v1[1], orig, AXIS::Y);
				   Float vz = get_face_value(grid.v1[2], orig, AXIS::Z);
				   Vec3f p = trace(orig, Vec3f(vx, vy, vz), dt);
				   grid.v0[1].get(i, j + 1, k) = get_face_value(grid.v1[1], p, AXIS::Y);
			   }
			   if (k + 1 < grid.grid_num_z && grid.occupy.get(i, j, k + 1) == 0) {
				   Vec3f orig(i + 0.5, 0.5 + j, 1.0 + k);
				   Float vx = get_face_value(grid.v1[0], orig, AXIS::X);
				   Float vy = get_face_value(grid.v1[1], orig, AXIS::Y);
				   Float vz = get_face_value(grid.v1[2], orig, AXIS::Z);
				   Vec3f p = trace(orig, Vec3f(vx, vy, vz), dt);
				   grid.v0[2].get(i, j, k + 1) = get_face_value(grid.v1[2], p, AXIS::Z);
			   }
            }
        }
    }
   
	//apply_boundary();
}

void SmokeSolver::advect_density_temp() {
     // swap density0,1 , temp0,1
    std::copy(grid.density0.data.begin(), grid.density0.data.end(), grid.density1.data.begin());
    std::copy(grid.temp0.data.begin(), grid.temp0.data.end(), grid.temp1.data.begin());
    for(int i=0;i< grid.grid_num_x; i++) {
        for(int j=0;j< grid.grid_num_y; j++) {
            for(int k=0;k< grid.grid_num_z; k++) {
                // ignore solid object
                if(grid.occupy.get(i,j,k) != 0) {
                    continue;
                } 
                Vec3f center = 	Vec3f(0.5 + i, 0.5 + j ,0.5 + k);
                Vec3f v_orig = Vec3f(get_face_value(grid.v1[0], center, AXIS::X), 
                                get_face_value(grid.v1[1],center, AXIS::Y),
                                get_face_value(grid.v1[2],center, AXIS::Z));
                Vec3f p = trace(center, v_orig, dt);
                grid.density0.get(i,j,k) =  get_center_value(grid.density1,p);
				grid.temp0.get(i, j, k) = get_center_value(grid.temp1, p);
            }
        }
    }   

}

void SmokeSolver::project_pressure() {
    A.setZero();
    b.setZero();
    x.setZero();
	triple.clear();

    // 1. compute divergence and matrix A
    Float inv_dt = 1.0 / dt;
	for (int k = 0; k < grid.grid_num_z; k++) {
		for (int j = 0; j < grid.grid_num_y; j++) {
			for (int i = 0; i < grid.grid_num_x; i++) {
				int idx = ijk2idx(i, j, k);

				// solid boundary condition
				if (grid.occupy.get(i, j, k) != 0) {
					b(idx) = 0;
					continue;
				}
				Float vx1 = grid.v0[0].get(i + 1, j, k);
				Float vx0 = grid.v0[0].get(i, j, k);

				Float vy1 = grid.v0[1].get(i, j + 1, k);
				Float vy0 = grid.v0[1].get(i, j, k);

				Float vz1 = grid.v0[2].get(i, j, k + 1);
				Float vz0 = grid.v0[2].get(i, j, k);

				b(idx) = (vx1 + vy1 + vz1 - vx0 - vy0 - vz0) * inv_dt;
			}
		}
	}

	Grid3DInt diag(grid.grid_num_x, grid.grid_num_y, grid.grid_num_z);
	for (int k = 0; k < grid.grid_num_z; k++) {
		for (int j = 0; j < grid.grid_num_y; j++) {
			for (int i = 0; i < grid.grid_num_x; i++) {
				bool fluid = grid.occupy.get(i, j, k) == 0;
				bool fluid_x = ( (i + 1) < grid.grid_num_x) && (grid.occupy.get(i + 1, j, k) == 0);
				bool fluid_y = ( (j + 1) < grid.grid_num_y) && (grid.occupy.get(i, j + 1, k) == 0);
				bool fluid_z = ( (k + 1) < grid.grid_num_z) && (grid.occupy.get(i, j, k + 1) == 0);


				if (fluid && fluid_x) {
					diag.get(i, j, k) += 1;
					diag.get(i + 1, j, k) += 1;
					triple.push_back(Eigen::Triplet<Float>(ijk2idx(i, j, k), ijk2idx(i + 1, j, k), 1));
					triple.push_back(Eigen::Triplet<Float>(ijk2idx(i+1, j, k), ijk2idx(i, j, k), 1));
				}

				if (fluid && fluid_y) {
					diag.get(i, j, k) += 1;
					diag.get(i, j+1, k) += 1;
					triple.push_back(Eigen::Triplet<Float>(ijk2idx(i, j, k), ijk2idx(i, j + 1, k), 1));
					triple.push_back(Eigen::Triplet<Float>(ijk2idx(i, j + 1, k), ijk2idx(i, j, k), 1));
				}

				if (fluid && fluid_z) {
					diag.get(i, j, k) += 1;
					diag.get(i, j, k + 1) += 1;
					triple.push_back(Eigen::Triplet<Float>(ijk2idx(i, j, k), ijk2idx(i, j, k + 1), 1));
					triple.push_back(Eigen::Triplet<Float>(ijk2idx(i, j, k + 1), ijk2idx(i, j, k), 1));
				}

            }
        }
    }
	for (int k = 0; k < grid.grid_num_z; k++) {
		for (int j = 0; j < grid.grid_num_y; j++) {
			for (int i = 0; i < grid.grid_num_x; i++) {
				int cnt = diag.get(i, j, k);
				if(cnt != 0)
					triple.push_back(Eigen::Triplet<Float>(ijk2idx(i, j, k), ijk2idx(i, j, k), -Float(cnt)));
			}
		}
	}
    // 2. solve linear system
    A.setFromTriplets(triple.begin(), triple.end());
    CG.compute(A);
	if (CG.info() != Eigen::Success) {
		std::cerr << "error occured, CG failed" << std::endl;
	}
    x = CG.solve(b);

    for(int i=0;i< grid.grid_num_x; i++) {
        for(int j=0;j< grid.grid_num_y; j++) {
            for(int k=0;k< grid.grid_num_z; k++) {
				grid.pressure.get(i,j,k) = x(ijk2idx(i,j,k));
            }
        }
     }

     // 3. update v
    for(int i=0;i< grid.grid_num_x; i++) {
        for(int j=0;j< grid.grid_num_y; j++) {
            for(int k=0;k< grid.grid_num_z; k++) {
				if (grid.occupy.get(i, j, k) != 0) continue;
				if (i + 1 < grid.grid_num_x && grid.occupy.get(i+1,j,k)==0 ) grid.v0[0].get(i + 1, j, k) -= (grid.pressure.get(i + 1, j, k) - grid.pressure.get(i, j, k)) * dt;
				if (j + 1 < grid.grid_num_y && grid.occupy.get(i, j + 1, k) == 0) grid.v0[1].get(i, j + 1, k) -= (grid.pressure.get(i, j + 1, k) - grid.pressure.get(i, j, k)) * dt;
				if (k + 1 < grid.grid_num_z && grid.occupy.get(i, j, k + 1) == 0) grid.v0[2].get(i, j, k + 1) -= (grid.pressure.get(i, j, k + 1) - grid.pressure.get(i, j, k)) * dt;
            }
        }
    }
	apply_boundary();
}

void SmokeSolver::dump_density(int frameCnt) {
	char basename[256];
	snprintf(basename, sizeof(basename), "solid/frame_%06d.vol", frameCnt);
	std::string filename = std::string(basename);
	std::ofstream file(filename.c_str(), std::ofstream::binary);

	int xres = grid.grid_num_x,
		yres = grid.grid_num_y,
		zres = grid.grid_num_z;
	int cnt = 0;
	float scale = 1.0f / std::max(std::max(xres, yres), zres);
	if (file) {
		printf("Writing %s...\n", filename.c_str());
		char header[48];
		memset(header, 0, sizeof(header));

		header[0] = 'V';
		header[1] = 'O';
		header[2] = 'L';
		header[3] = 3;
		int32_t* encoding = reinterpret_cast<int32_t*>(header + 4);
		encoding[0] = 1; 
		encoding[1] = static_cast<int32_t>(xres);
		encoding[2] = static_cast<int32_t>(yres);
		encoding[3] = static_cast<int32_t>(zres);
		encoding[4] = 1;  
		float* bbox = reinterpret_cast<float*>(encoding + 5);

		float minX = -0.5;
		float minY = -1.0;
		float minZ = -0.5;
		float maxX = 0.5;
		float maxY = 1.0;
		float maxZ = 0.5;
		bbox[0] = static_cast<float>(minX);
		bbox[1] = static_cast<float>(minY);
		bbox[2] = static_cast<float>(minZ);
		bbox[3] = static_cast<float>(maxX);
		bbox[4] = static_cast<float>(maxY);
		bbox[5] = static_cast<float>(maxZ);

		file.write(header, sizeof(header));

		for (int i = 0; i < grid.grid_num_x; i++) {
			for (int j = 0; j < grid.grid_num_y; j++) {
				for (int k = 0; k < grid.grid_num_z; k++) {

					double value = grid.density0.get(i, j, k);

					if (grid.occupy.get(i, j, k) != 0) value = -1.0;
					float tmp = value;
					file.write(reinterpret_cast<const char*>(&tmp), sizeof(float));
					if (tmp != 0) cnt++;
				}
			}
		}
		file.close();
	}
	printf("%d\n", cnt);
}
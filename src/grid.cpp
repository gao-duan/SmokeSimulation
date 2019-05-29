#include "grid.h"
#include <fstream>
#include <assert.h>
Grid3D::Grid3D(int _num_x, int _num_y, int _num_z)
	: num_x(_num_x), num_y(_num_y), num_z(_num_z){
		size = _num_x * _num_y * _num_z;
		data.resize(size,0);
}

Float Grid3D::get(int x, int y, int z) const {
	if (x < 0) x = 0;
	if (x > num_x - 1) x = num_x - 1;
	if (y < 0) y = 0;
	if (y > num_y - 1) y = num_y - 1;
	if (z < 0) z = 0;
	if (z > num_z - 1) z = num_z - 1;
	
	return data[x + num_x * y + num_x * num_y * z];
}

Float& Grid3D::get(int x, int y, int z) {
	if (x < 0) x = 0;
	if (x > num_x - 1) x = num_x - 1;
	if (y < 0) y = 0;
	if (y > num_y - 1) y = num_y - 1;
	if (z < 0) z = 0;
	if (z > num_z - 1) z = num_z - 1;
	
	return data[x + num_x * y + num_x * num_y * z];
}

Grid3DInt::Grid3DInt(int _num_x, int _num_y, int _num_z)
	: num_x(_num_x), num_y(_num_y), num_z(_num_z){
		size = _num_x * _num_y * num_z;
		data.resize(size, 0);
}

int Grid3DInt::get(int x, int y, int z) const {
	/*if (x < 0 || x >= num_x || y < 0 || y >= num_y || z < 0 || z >= num_z) {
		cout << "ERROR in Grid3DInt::get " << x << "," << y << "," << z << "; " << num_x << "," << num_y << "," << num_z << endl;
	}*/
	if (x < 0) x = 0;
	if (x > num_x - 1) x = num_x - 1;
	if (y < 0) y = 0;
	if (y > num_y - 1) y = num_y - 1;
	if (z < 0) z = 0;
	if (z > num_z - 1) z = num_z - 1;
	return data[x + num_x * y + num_x * num_y * z];
}

int& Grid3DInt::get(int x, int y, int z) {
	/*if (x < 0 || x >= num_x || y < 0 || y >= num_y || z < 0 || z >= num_z) {
		cout << "ERROR in Grid3DInt::get " << x << "," << y << "," << z << "; " << num_x << "," << num_y << "," << num_z << endl;
	}*/
	if (x < 0) x = 0;
	if (x > num_x - 1) x = num_x - 1;
	if (y < 0) y = 0;
	if (y > num_y - 1) y = num_y - 1;
	if (z < 0) z = 0;
	if (z > num_z - 1) z = num_z - 1;
	return data[x + num_x * y + num_x * num_y * z];
}


Float Grid3D::get(const Vec3f& p) const {
	return trilinear_interp(p);
	//return monotonic_cubic_interp(p);
}

void Grid3D::dump(const std::string & name) const
{
	std::ofstream os(name);
	for (int i = 0; i < size; ++i) {
		os << data[i] << std::endl;
	}
	os.close();
}

// tri-linear: https://en.wikipedia.org/wiki/Trilinear_interpolation
Float Grid3D::trilinear_interp(const Vec3f & p) const
{
	Float x = p.x, y = p.y, z = p.z;
	x = std::min(x, ((Float)num_x - 1 - 1e-6));
	y = std::min(y, ((Float)num_y - 1 - 1e-6));
	z = std::min(z, ((Float)num_z - 1 - 1e-6));

	x = std::max(0.0, x);
	y = std::max(0.0, y);
	z = std::max(0.0, z);

	int x0 = int(x);
	int y0 = int(y);
	int z0 = int(z);

	int x1 = x0 + 1;
	int y1 = y0 + 1;
	int z1 = z0 + 1;
	Float xd = (x - x0); // (x-x0) / (x1-x0)
	Float yd = (y - y0); // (y-z0) / (y1-y0)
	Float zd = (z - z0); // (z-z0) / (z1-z0)
	Float c000 = get(x0, y0, z0);
	Float c001 = get(x0, y0, z1);
	Float c010 = get(x0, y1, z0);
	Float c011 = get(x0, y1, z1);
	Float c100 = get(x1, y0, z0);
	Float c101 = get(x1, y0, z1);
	Float c110 = get(x1, y1, z0);
	Float c111 = get(x1, y1, z1);

	Float c00 = c000 * (1.0 - xd) + c100 * xd;
	Float c01 = c001 * (1.0 - xd) + c101 * xd;
	Float c10 = c010 * (1.0 - xd) + c110 * xd;
	Float c11 = c011 * (1.0 - xd) + c111 * xd;
	Float c0 = c00 * (1.0 - yd) + c10 * yd;
	Float c1 = c01 * (1.0 - yd) + c11 * yd;
	Float c = c0 * (1.0 - zd) + c1 * zd;
	return c;
}

// arr: f_k-1, f_k, f_k+1, f_k+2
Float monntonic_cubic_1D(const std::vector<Float>& arr, Float t) {
	Float dk = (arr[2] - arr[0]) / 2.0;
	Float dk1 = (arr[3] - arr[1]) / 2.0;
	Float delta_k = (arr[2] - arr[1]);

	if (delta_k == 0) {
		dk = 0;
		dk1 = 0;
	}
	else {
		if (delta_k > 0) {
			dk = std::abs(dk);
			dk1 = std::abs(dk1);
		}
		else {
			dk = -std::abs(dk);
			dk1 = -std::abs(dk1);
		}
	}

	Float a3 = dk + dk1 - 2 * delta_k;
	Float a2 = 3.0 * delta_k - 2 * dk - dk1;
	Float a1 = dk;
	Float a0 = arr[1];

	return a3 * t * t * t + a2 * t * t + a1 * t + a0;
}

Float Grid3D::monotonic_cubic_interp(const Vec3f & p) const
{
	Float x = p.x, y = p.y, z = p.z;
	x = std::min(x, ((Float)num_x - 1 - 1e-6));
	y = std::min(y, ((Float)num_y - 1 - 1e-6));
	z = std::min(z, ((Float)num_z - 1 - 1e-6));

	x = std::max(0.0, x);
	y = std::max(0.0, y);
	z = std::max(0.0, z);

	int x0 = int(x);
	int y0 = int(y);
	int z0 = int(z);

	Float xd = (x - x0); // (x-x0) / (x1-x0)
	Float yd = (y - y0); // (y-z0) / (y1-y0)
	Float zd = (z - z0); // (z-z0) / (z1-z0)

	std::vector<Float> Z(4,0);
	for (int k = 0; k < 4; ++k) {
		std::vector<Float> X(4, 0);
		for (int i = 0; i < 4; ++i) {
			int x1 = x0 + i - 1;
			int z1 = z0 + k - 1;
			X[i] = monntonic_cubic_1D({ get(x1, y0-1, z1),
									get(x1, y0, z1), get(x1, y0+1, z1), get(x1, y0+2, z1) }, yd);
		}
		Z[k] = monntonic_cubic_1D(X, xd);
	}

	return monntonic_cubic_1D(Z, zd);
}


Float get_center_value(const Grid3D& data, const Vec3f& p) {
    Vec3f offset(0.5,0.5,0.5);
    return data.get(p - offset);
}

Float get_face_value(const Grid3D& data, const Vec3f& p, AXIS axis) {
    Vec3f offset;
    switch(axis){
        case AXIS::X: offset = Vec3f(0,0.5,0.5); break;
        case AXIS::Y: offset = Vec3f(0.5,0,0.5); break;
        case AXIS::Z: offset = Vec3f(0.5,0.5,0); break;
    }
    return data.get(p - offset);
}


MACGrid::MACGrid(int _grid_num_x, int _grid_num_y, int _grid_num_z)
: grid_num_x(_grid_num_x), grid_num_y(_grid_num_y), grid_num_z(_grid_num_z) {
    v0[0] = Grid3D(grid_num_x + 1, grid_num_y, grid_num_z);
    v1[0] = Grid3D(grid_num_x + 1, grid_num_y, grid_num_z);
    v0[1] = Grid3D(grid_num_x, grid_num_y + 1, grid_num_z);
    v1[1] = Grid3D(grid_num_x, grid_num_y + 1, grid_num_z);
    v0[2] = Grid3D(grid_num_x, grid_num_y, grid_num_z + 1);
    v1[2] = Grid3D(grid_num_x, grid_num_y, grid_num_z + 1);

    F[0] = Grid3D(grid_num_x, grid_num_y, grid_num_z);
    F[1] = Grid3D(grid_num_x, grid_num_y, grid_num_z);
    F[2] = Grid3D(grid_num_x, grid_num_y, grid_num_z);
    
    density0 = Grid3D(grid_num_x, grid_num_y, grid_num_z);
    density1 = Grid3D(grid_num_x, grid_num_y, grid_num_z);
    
    temp0 = Grid3D(grid_num_x, grid_num_y, grid_num_z);
    temp1 = Grid3D(grid_num_x, grid_num_y, grid_num_z);

    pressure = Grid3D(grid_num_x, grid_num_y, grid_num_z);

    occupy = Grid3DInt(grid_num_x, grid_num_y, grid_num_z);


	avg_v[0] = Grid3D(grid_num_x, grid_num_y, grid_num_z);
    avg_v[1] = Grid3D(grid_num_x, grid_num_y, grid_num_z);
    avg_v[2] = Grid3D(grid_num_x, grid_num_y, grid_num_z);

	vorticity[0] = Grid3D(grid_num_x, grid_num_y, grid_num_z);
    vorticity[1] = Grid3D(grid_num_x, grid_num_y, grid_num_z);
    vorticity[2] = Grid3D(grid_num_x, grid_num_y, grid_num_z);

	vorticity_norm = Grid3D(grid_num_x, grid_num_y, grid_num_z);

}
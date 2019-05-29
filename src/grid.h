#ifndef __GRID_H__
#define __GRID_H__

#include "common.h"
#include <vector>
#include <cmath>
class SmokeSolver;

template<typename T>
struct Vector3 {
	T x, y, z;
	Vector3(T _x=0, T _y=0, T _z=0):x(_x), y(_y), z(_z) {}
	Vector3<T> operator+(const Vector3<T>& rhs) const {
		return Vector3<T>(x + rhs.x, y + rhs.y, z + rhs.z);
	}
	Vector3<T> operator-(const Vector3<T>& rhs) const {
		return Vector3<T>(x - rhs.x, y - rhs.y, z - rhs.z);
	}
	Vector3<T> operator*(const Float t) const {
		return Vector3<T>(x * t, y * t ,z * t);
	}
	Float length() const {
		return std::sqrt(x * x + y * y + z * z);
	}
};

template<typename T>
inline T dot(const Vector3<T> &v1, const Vector3<T> &v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

template<typename T>
inline Vector3<T> cross(const Vector3<T> &v1, const Vector3<T> &v2) {
	return Vector3<T>(
		(v1.y * v2.z) - (v1.z * v2.y), 
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x)
	);
}


typedef Vector3<Float> Vec3f;
typedef Vector3<int> Vec3i;

struct Grid3D {
	Grid3D(){}
	~Grid3D(){}
	Grid3D(int _num_x, int _num_y, int _num_z);
	Float get(int x, int y, int z) const;
	Float& get(int x, int y, int z); 

	Float get(const Vec3f& p) const;
	void dump(const std::string& name) const;
	std::vector<Float> data;
	int num_x, num_y, num_z;
	int size;
private:
	// tri-linear: https://en.wikipedia.org/wiki/Trilinear_interpolation
	Float trilinear_interp(const Vec3f &p) const;
	Float monotonic_cubic_interp(const Vec3f& p) const;
};


struct Grid3DInt {
	Grid3DInt() {}
	Grid3DInt(int _num_x, int _num_y, int _num_z);
	int get(int x, int y, int z) const;
	int& get(int x, int y, int z); 
	
	
	std::vector<int> data;
	int num_x, num_y, num_z;
	int size;
};

enum class AXIS {X,Y,Z};
Float get_center_value(const Grid3D& data, const Vec3f& p);
Float get_face_value(const Grid3D& data, const Vec3f& p, AXIS axis);

class MACGrid {
public:
	MACGrid(int _x, int _y, int _z);
	~MACGrid() {}
	friend class SmokeSolver;
private:
	Grid3D v0[3],v1[3];   	 		  // velocity
	Grid3D density0, density1;		 // density
	Grid3D temp0, temp1;			// temperature
	Grid3D F[3]; 		  			// force
	Grid3D pressure; 			   // pressure
	Grid3DInt occupy;			   // voxel is smoke or solid
	int grid_num_x, grid_num_y, grid_num_z;
	Float dt;

	Grid3D avg_v[3];
	Grid3D vorticity[3];
	Grid3D vorticity_norm;
};
#endif
#ifndef VECTOR3D_H
#define VECTOR3D_H

#include "vector3d.h"

class vector3d{
public:
    double x,y,z;

    vector3d();

    vector3d(double X,double Y,double Z);

    bool operator==(vector3d& a);

    double length();

    vector3d operator/(double a);

    vector3d operator+(vector3d a);

    vector3d operator-(vector3d a);

    vector3d cross(vector3d a);

    vector3d norm();

    double operator*(vector3d a);

    vector3d operator*(double a);



};


#endif // VECTOR3D_H

#include "vector3d.h"
#include <math.h>

vector3d::vector3d(){
    x=0;
    y=0;
    z=0;
}

vector3d::vector3d(double X, double Y, double Z){
    x=X;
    y=Y;
    z=Z;
}

bool vector3d::operator==(vector3d &a){
    if(x==a.x&& y==a.y && z==a.z)
        return true;
    else
        return false;
}

double vector3d::length(){
    double L=0;
    L=sqrt(x*x+y*y+z*z);
    return L;
}

vector3d vector3d::operator/(double a){
    if(a!=0){
        vector3d b(x/a,y/a,z/a);
        return b;
    }
    else{
        vector3d b(0,0,0);
        return b;
    }
}

vector3d vector3d::operator+(vector3d a){
    vector3d b(x+a.x,y+a.y,z+a.z);
    return b;
}

vector3d vector3d::operator-(vector3d a){
    vector3d b(x-a.x,y-a.y,z-a.z);
    return b;
}

vector3d vector3d::cross(vector3d a){
    vector3d b(0,0,0);
    b.x=y*a.z-z*a.y;
    b.y=-x*a.z+z*a.x;
    b.z=x*a.y-y*a.x;
    return b;

}

vector3d vector3d::norm(){
    vector3d a=*this;
    if(this->length()!=0){
        vector3d b=a/this->length();
        return b;
    }
    return a;

}

vector3d vector3d::operator*(double a){
    vector3d b(a*x,a*y,a*z);
    return b;
}

double vector3d::operator*(vector3d a){
    return x*a.x+y*a.y+z*a.z;
}



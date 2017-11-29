#ifndef ANGULARSPRING_H
#define ANGULARSPRING_H
#include <iostream>
using namespace std;
class vector3d;
class mass;

class AngularSpring
{
public:
    AngularSpring();
    AngularSpring(double theta0,double k);

    double Theta0,K;
    double Theta;

    vector3d *r0,*r1;
    mass *m[3];

    double getTheta();

    void AssignMass(mass *Mi,mass *Mmid,mass *Mf);
    void stretch(double dL);

    vector3d getRM() ;  //middle
    vector3d getR1() ;
    vector3d getR2() ;

    vector3d getF1();
    vector3d getF2();
    vector3d getFMid();
};

#endif // ANGULARSPRING_H

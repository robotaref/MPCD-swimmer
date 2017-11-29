#include "angularspring.h"
#include "vector3d.h"
#include "mass.h"
#include <math.h>

AngularSpring::AngularSpring()
{

}

AngularSpring::AngularSpring(double theta0, double k)
{
    this->Theta0=theta0;
    this->K=k;
}

double AngularSpring::getTheta()
{
    double a=acos(getR1().norm()*getR2().norm());

    return a;
}

void AngularSpring::AssignMass(mass *Mi, mass *Mmid, mass *Mf)
{

    Mmid->AddAngularSpring(this,0);
    Mi->AddAngularSpring(this,-1);
    Mf->AddAngularSpring(this,1);


}

vector3d AngularSpring::getRM()
{
    return (this->m[1]->r);
}

vector3d AngularSpring::getR1()
{
    vector3d r0=(m[0]->r-m[1]->r);
    return r0;
}

vector3d AngularSpring::getR2()
{
    vector3d r1=(m[2]->r-m[1]->r);
    return r1;
}

vector3d AngularSpring::getF1()
{
    vector3d r0=getR1();
    vector3d r1=getR2();
    vector3d T=r1.norm().cross(r0.norm());
    return T.norm().cross(r0.norm())/(r0.length())*K*(-Theta0+this->getTheta());

}

vector3d AngularSpring::getF2()
{
    vector3d r0=getR1();
    vector3d r1=getR2();
    vector3d T=r1.norm().cross(r0.norm());
    return T.norm().cross(r1.norm())/(r1.length())*K*(-Theta0+this->getTheta());

}

vector3d AngularSpring::getFMid()
{
    return getF2()*(-1) - getF1();
}

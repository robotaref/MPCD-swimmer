#include "spring.h"
#include "vector3d.h"
#include "mass.h"
#include <iostream>
#include <math.h>
using namespace std;
spring::spring()
{

}

spring::spring(double l0,double k)
{
    L0=l0;
    K=k;

}

vector3d *spring::getR1() const
{
    return &(this->m[1]->r);
}

vector3d *spring::getR0() const
{
    return &(this->m[0]->r);
}


double spring::getL()
{
    vector3d c;
    c.x=(this->getR1()->x)-(this->getR0()->x);
    c.y=(this->getR1()->y)-(this->getR0()->y);
    c.z=(this->getR1()->z)-(this->getR0()->z);


//    if(fabs((this->getR1()->x)-(this->getR0()->x))<10)
//        c.x= fabs();
//    else
//        c.x= fabs((this->getR1()->x)-(this->getR0()->x))-20;

//    if(fabs((this->getR1()->y)-(this->getR0()->y))<15)
//        c.y= fabs((this->getR1()->y)-(this->getR0()->y));
//    else
//        c.y= fabs((this->getR1()->y)-(this->getR0()->y))+30;

//    if(fabs((this->getR1()->z)-(this->getR0()->z))<15)
//        c.z= fabs((this->getR1()->z)-(this->getR0()->z));
//    else
//        c.z= fabs((this->getR1()->z)-(this->getR0()->z))+30;

    return c.length();
    
}

vector3d spring::getF()
{
    vector3d F;
    F=*getR0()-*getR1();
    if(F.x<-10)
        F.x+=20;

    F=F.norm();
    F=F*(1*K*(getL()-L0));
//    cout<<getL()<<"\t"<<this->L0<<"\t"<<this->K<<"\t"<<F.x<<endl;
    return F;

}

void spring::stretch(double dL)
{
    L0+=dL;

}

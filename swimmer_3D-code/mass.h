#ifndef MASS_H
#define MASS_H
#include <vector>
#include "vector3d.h"
using namespace std;
class environment;
class vector3d;
class spring;
class AngularSpring;

class mass{

public:
    vector<spring*> SpringList;
    vector<int> SpringType;
    vector<AngularSpring*> ASpringList;
    vector<int> ASpringType;

    vector3d r;
    vector3d v;
    vector3d F;

    double m;
    int ID;
    environment *e;

    mass(int Id,environment *E);
    mass(double M, vector3d R,vector3d V);


    void AddSpring(spring *ss,int st);
    void AddAngularSpring(AngularSpring *ss,int st);

    void FindF();
    void Move(double dt);
};

#endif // MASS_H

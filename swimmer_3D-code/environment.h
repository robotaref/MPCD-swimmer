#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h>
#include <cmath>
#include <iomanip>
#include <limits>
#include <iomanip>
#include <cstdlib>
//#include <random>
//int TRYNUM =10;
//double len=2;


////-----------------general
//#define step  8000// number of MD stps steps
#define savingstep   1//saving trajectory with distance of savingstep
#define  h     0.001 // time length of steps in MD
//#define KT     0.5  // KT the quanta of energy
#define pi     3.141592 // clear !
//#define save_condition 0
//#define rand_init 0
//#define active_swimmer 0
//#define shiftforce   0.1
//double NodeMass=  100000;
////-----------------general



////-------------------------------Sovlent-----------------------------
//// periodic in x direction and in  walls are in y and z direction
#define showsolvent   1 //  1=show   0=dont show they are too many!
#define BoundaryType  2 //(1=slip and no force) (2=no slip and no force) (3=slip  and no force  and initial flow)
////(4=noslip and no force and  initial velocity)  (5=noslip and constant force) (6=shear flow)
//#define Fsolvent   0 // if BoundaryType==5
//#define Vinitialflow 10.0 // initial flow if  BoundaryType==3,4
#define Vshear   20.0 // if BoundaryType==6
#define collisiontimestep  30  /// collision/inreaction == must be integer
#define interactionstep    1  /// collision/inreaction == must be integer
//#define  msolvent 0.2  //mass of solvent particle
//#define  lbs  0.4 // thikness of membrane film in solvebt-triangle interaction
#define wallthickness 0.2 // thickness of wall in slip or no-slip boundary condition
int const Lx = 20.0;/// must be even number!!!         size of simulation box
int const Ly = 30.0;/// must be even number!!!         size of simulation box
int const Lz = 30.0;/// must be even number!!!         size of simulation box
//#define VmaxSolvent 1.0 // range of randomly distributed velocities in initialization
int const solventdensityinxox=10;// number of solvent particles in eachb box in initialization
const int nb = Lx*Ly*Lz ; //number of boxes
////const int nbody = 0;
const int nsolvent = Lx*Ly*Lz*solventdensityinxox ; // number of solvent particles
 #define spliteachboxfordiagram   1

using namespace std;
class mass;
class spring;
class AngularSpring;

class environment
{
public:

    int eachSolventinwichBox [nb][1000*solventdensityinxox]; //  wich particle is in each box!
    int howmanySolventsineachBox[nb];                       // how many particle in each box
    double vdiagram[Lz*spliteachboxfordiagram][2];

    vector<mass*> M;
    vector<spring*> Ss;
    vector<AngularSpring*> AS;

    ofstream trajectory;
    ofstream com;
    ofstream vd;

    string address;

    int rand_init;
    int HeavyNum;
    double VmaxSolvent;
    double len;
    double NodeMass;
    double LX,LY,LZ;

    void initialzesolvent();
    void solventBoundaries();
    void SaveLastCondition();
    void MDStep(int time);
    void saveTrajectory(int time);
    void collision(int time);
    void updateHowmanySolventsineachBox();
    void SwimmimgSenarioNG(int time);

    environment();
    environment(int TryNum);

};

#endif // ENVIRONMENT_H

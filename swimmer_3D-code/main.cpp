
#include<time.h>
//-----------------------HEADERS
#include <iostream>
#include<fstream>
#include<cmath>
#include <iomanip>
#include <limits>
#include <iomanip>
#include <cstdlib>
//#include <random>
#include "vector3d.h"
#include "spring.h"
#include <math.h>
#include "environment.h"
int TRYNUM=1;
int step=10000;
using namespace std;

void SPRING(int t){
//    double dl=0;
//    if(active_swimmer==1){
//        if(t%1200<100)
//            Ss[0]->stretch(-dl);
//        else if(t%1200<200)
//            Ss[2]->stretch(-dl);
//        else if(t%1200<300)
//            Ss[3]->stretch(-dl);
//        else if(t%1200<400)
//            Ss[5]->stretch(-dl);
//        else if(t%1200<500)
//            Ss[4]->stretch(-dl);
//        else if(t%1200<600)
//            Ss[1]->stretch(-dl);
//        else if(t%1200<700)
//            Ss[1]->stretch( dl);
//        else if(t%1200<800)
//            Ss[4]->stretch( dl);
//        else if(t%1200<900)
//            Ss[5]->stretch( dl);
//        else if(t%1200<1000)
//            Ss[3]->stretch( dl);
//        else if(t%1200<1100)
//            Ss[2]->stretch( dl);
//        else if(t%1200<1200)
//            Ss[0]->stretch( dl);
//    }

//    for(int i=0;i<4;i++){
//        M[i]->FindF();
//    }

//    for(int i=0;i<4;i++)
//        M[i]->Move(interactionstep*h);


}


int main()
{


    for(int TryNum=0;TryNum<TRYNUM;TryNum++)
    {
        environment *e= new environment(TryNum);
        cout<<"Simulation "<<TryNum+1<<" out of "<<TRYNUM<<" starting: "<<endl;
        for(int i=1 ;i<=step ; i++) // main MD loop
        {
            e->MDStep(i);
            e->saveTrajectory(i);
        } /// END of MD

        e->SaveLastCondition();
    }

}













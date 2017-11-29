#include "mass.h"
#include "spring.h"
#include <iostream>
#include "environment.h"
#include "angularspring.h"

using namespace std;
mass::mass(int Id,environment *E){
    m=0;
    this->e=E;
    ID=Id;
}

mass::mass(double M, vector3d R, vector3d V){
    m=M;
    r=R;
    v=V;
}

void mass::AddSpring(spring *ss,int st)
{
    SpringList.push_back(ss);
    SpringType.push_back(st);
    if(st==-1)
        ss->m[0]=this;
    if(st==1)
        ss->m[1]=this;

}

void mass::AddAngularSpring(AngularSpring *ss, int st)
{
    ASpringList.push_back(ss);
    ASpringType.push_back(st);
    if(st==-1)
        ss->m[0]=this;
    if(st==1)
        ss->m[2]=this;
    if(st==0)
        ss->m[1]=this;
}

void mass::FindF()
{
    vector3d zero(0,0,0);
    //    F=zero;

    for(int i=0;i<SpringList.size();i++){
        F=F+SpringList[i]->getF()*SpringType[i];
    }

    for(int i=0;i<ASpringList.size();i++){
        if(SpringType[i]==0)
            F=F+ASpringList[i]->getFMid();
        else if(SpringType[i]==-1)
            F=F+ASpringList[i]->getF1();
        else if(SpringType[i]==1)
            F=F+ASpringList[i]->getF2();

    }


//    int counter=0;
    int ini=e->HeavyNum;
    if(ID+1 > e->HeavyNum)
        ini=ID+1;
    for(int i=ini;i<this->e->M.size();i++){
        if(this->ID!=e->M[i]->ID){
            vector3d g=(this->r - e->M[i]->r);
            while(g.x>e->LX/2)
                g.x=g.x-e->LX/2;
            while(g.x<= -e->LX/2)
                g.x=g.x+e->LX/2;
            if(g.length()<1)
                g=g.norm();
            if(g.length()<2){
//                counter++;
                double lw=pow(g.length()*2,-6);
                g=g.norm()*1000*(-lw+ lw*lw);

                this->F=F-g;
                e->M[i]->F=e->M[i]->F+g;
            }

        }

    }
//    cout<<this->F.x<<"\t";

}

void mass::Move(double dt)
{
    v=v+F*dt/this->m;
    r=r+v*dt;

}

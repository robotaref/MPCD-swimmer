#include "environment.h"
#include "mass.h"
#include "spring.h"

void environment::initialzesolvent()
{
    for(int i=0;i<nsolvent;i++){
        M.push_back( new mass(i,this));
    }
    if(rand_init==1){
        cout<<"Random initial condition"<<endl;
        if(BoundaryType ==1 || BoundaryType==2  || BoundaryType==5 || BoundaryType==6)
        {
            for (int i = 0; i < nsolvent; i++) // random place and velocity
            {
                (M[i]->r.x)= (double) (rand() % (Lx * 1000)) / 1000 - (Lx/2);
                (M[i]->r.y) = (double) (rand() % (Ly * 1000)) / 1000 - (Ly/2);
                (M[i]->r.z) = (double) (rand() % (Lz * 1000)) / 1000 - (Lz/2);

                (M[i]->v.z) =(double) (rand() % (int) (VmaxSolvent * 1000)) / 1000 - ((double)VmaxSolvent/2);
                (M[i]->v.y) =(double) (rand() % (int) (VmaxSolvent * 1000)) / 1000 - ((double) VmaxSolvent/2);
                (M[i]->v.x) =(double) (rand() % (int) (VmaxSolvent * 1000)) / 1000 - ((double)VmaxSolvent/2);

                M[i]->m=1;
            }
        }
        if(BoundaryType ==3 || BoundaryType==4)
        {
            for (int i = 0; i < nsolvent; i++) // random place and velocity
            {
                M[i]->r.x = (double) (rand() % (Lx * 1000)) / 1000 - (Lx/2);
                M[i]->r.y = (double) (rand() % (Ly * 1000)) / 1000 - (Ly/2);
                M[i]->r.z = (double) (rand() % (Lz * 1000)) / 1000 - (Lz/2);

                M[i]->v.z = 0;
                M[i]->v.y = 0;
                M[i]->v.x = 0;

            }
        }
    }
    else{

        cout<<"Loading from file:"<<endl;
        ifstream init;
        init.open("init.txt");
        if(init.is_open())
            cout<<endl<<"Readable File Opened!"<<endl;
        else{
            cout<<"File corrupted!"<<endl;

        }
        double line;
        int i=0;
        while(init>>line){
            M[i/6]->m=1;
            //        double line1=atof(line.c_str());
            if(i%6==0)
                M[i/6]->r.x=line;
            else if(i%6==1)
                M[i/6]->r.y=line;
            else if(i%6==2)
                M[i/6]->r.z=line;
            else if(i%6==3)
                M[i/6]->v.x=line;
            else if(i%6==4)
                M[i/6]->v.y=line;
            else if(i%6==5)
                M[i/6]->v.z=line;
            i++;
        }
        vector3d zero(0,0,0);
//        vector3d Zero(0.5,0,0);
        vector3d initL(-len,0,0);
//        M[2]->r=Zero;
//        M[3]->r=Zero;
//        M[4]->r=Zero;
//        M[5]->r=Zero;

        M[0]->v=zero;
        M[0]->r=zero;
        M[0]->m=NodeMass;
        M[1]->v=zero;
        M[1]->r=initL;
        M[1]->m=NodeMass;
        Ss.push_back(new spring(len,100000));
        M[0]->AddSpring(Ss[0],-1);
        M[1]->AddSpring(Ss[0], 1);
        /*
        double a=4;
        double RRR=a/sqrt(3.0);
        double HHH=a*sqrt(2.0/3.0);

        (M[0]->r.x=RRR);
        (M[1]->r.x=RRR*cos(2.0*pi/3.0));
        (M[2]->r.x=RRR*cos(4.0*pi/3.0));
        (M[3]->r.x=0);

        (M[0]->r.y=0);
        (M[1]->r.y=RRR*sin(2.0*pi/3.0));
        (M[2]->r.y=RRR*sin(4.0*pi/3.0));
        (M[3]->r.y=0);

        (M[0]->r.z= -HHH/4.0);
        (M[1]->r.z= -HHH/4.0);
        (M[2]->r.z= -HHH/4.0);
        (M[3]->r.z=  HHH*3.0/4.0);


        M[0]->v.x=M[1]->v.x=M[2]->v.x=M[3]->v.x=0;
        M[0]->v.y=M[1]->v.y=M[2]->v.y=M[3]->v.y=0;
        M[0]->v.z=M[1]->v.z=M[2]->v.z=M[3]->v.z=0;

        M[0]->m=M[1]->m=M[2]->m=M[3]->m=1000;

        Ss[0]=new spring(a,100000);
        Ss[1]=new spring(a,100000);
        Ss[2]=new spring(a,100000);
        Ss[3]=new spring(a,100000);
        Ss[4]=new spring(a,100000);
        Ss[5]=new spring(a,100000);



        M[0]->AddSpring(Ss[0],-1);
        M[0]->AddSpring(Ss[2],1);
        M[0]->AddSpring(Ss[3],-1);

        M[1]->AddSpring(Ss[0],1);
        M[1]->AddSpring(Ss[1],-1);
        M[1]->AddSpring(Ss[5],1);

        M[2]->AddSpring(Ss[1],1);
        M[2]->AddSpring(Ss[2],-1);
        M[2]->AddSpring(Ss[4],1);

        M[3]->AddSpring(Ss[3],1);
        M[3]->AddSpring(Ss[4],-1);
        M[3]->AddSpring(Ss[5],-1);

*/

    }
}

void environment::solventBoundaries()
{

    if(BoundaryType==1  || BoundaryType==3)
    {
        for (int j = 0; j < nsolvent; j++)
        {
            // ___________________ Boundary Condition Solvent
            //x:
//            if ((M[j]->r.x) > Lx/2.0)
//                (M[j]->r.x) = (M[j]->r.x) - Lx ;

//            if ((M[j]->r.x) < -Lx/2.0)
//                (M[j]->r.x) = (M[j]->r.x) + Lx ;

            // slip condition
            //y and z :
            double ly =Ly/2.0-wallthickness;
            double lz =Lz/2.0-wallthickness;
            if (   (((M[j]->r.y)>ly) && ((M[j]->v.y) >0.0))   ||   (((M[j]->r.y)<-ly) && (((M[j]->v.y) <0.0)))  )
            {
                (M[j]->v.y) = -(M[j]->v.y);
            }

            else  if (   (((M[j]->r.z)>lz) && ((M[j]->v.z) >0.0))   ||   (((M[j]->r.z)<-lz) && ((M[j]->v.z) <0.0))  )
            {
                (M[j]->v.z) = -(M[j]->v.z);
            }
        }
    }

    else if(BoundaryType==2 ||  BoundaryType==4  || BoundaryType==5)
    {
        for (int j = 0; j < nsolvent; j++)
        {
            // ___________________ Boundary Condition Solvent
            //x:
//            if ((M[j]->r.x) > Lx/2.0)
//                (M[j]->r.x) =( M[j]->r.x) - Lx ;

//            if ((M[j]->r.x )< -Lx/2.0)
//                (M[j]->r.x) = (M[j]->r.x) + Lx ;

            // no-slip condition
            //y and z :
            double ly =Ly/2.0-wallthickness;
            double lz =Lz/2.0-wallthickness;
            if (   (((M[j]->r.y)>ly) &&( (M[j]->v.y) >0.0))   ||   (((M[j]->r.y)<-ly) && ((M[j]->v.y) <0.0))  )
            {
                (M[j]->v.x) = -(M[j]->v.x);
                (M[j]->v.y) = -(M[j]->v.y);
                (M[j]->v.z) = -(M[j]->v.z);
            }

            else    if (   (((M[j]->r.z)>lz) && ((M[j]->v.z) >0.0))   ||   (((M[j]->r.z)<-lz) && ((M[j]->v.z) <0.0))  )
            {
                (M[j]->v.x) = -(M[j]->v.x);
                (M[j]->v.y) = -(M[j]->v.y);
                (M[j]->v.z) = -(M[j]->v.z);
            }
        }
    }


    else if(BoundaryType==6)
    {
        for (int j = 0; j < nsolvent; j++)
        {
            // ___________________ Boundary Condition Solvent
            //x:
//            if ((M[j]->r.x) > Lx/2.0)
//                (M[j]->r.x) = (M[j]->r.x) - Lx ;

//            if ((M[j]->r.x) < -Lx/2.0)
//                (M[j]->r.x) = (M[j]->r.x) + Lx ;

            // shear flow in z plane and slip in y plane
            //y and z :
            double ly =Ly/2.0-wallthickness;
            double lz =Lz/2.0-wallthickness;

            if (   (((M[j]->r.y)>ly) && ((M[j]->v.y) >0.0))   ||   (((M[j]->r.y)<-ly) && ((M[j]->v.y) <0.0))  )
            {
                (M[j]->v.y) = -(M[j]->v.y);
            }
            else   if (   (((M[j]->r.z)>lz) & ((M[j]->v.z) >0.0))   ||   (((M[j]->r.z)<-lz) & ((M[j]->v.z) <0.0))  )
            {
                (M[j]->v.z) = -(M[j]->v.z);
            }


            if (  (M[j]->r.z)>lz     )
            {
                (M[j]->v.x) = +Vshear;
            }
            else   if (  (M[j]->r.z)<(-lz )    )
            {
                (M[j]->v.x) = -Vshear;
            }

        }

    }


}

void environment::SaveLastCondition()
{
    cout<<"Writing to file:"<<endl;
    ofstream init;
    init.open("init.txt");
    init << std:: fixed;

    for(int i=0;i<nsolvent;i++){
        init<<M[i]->r.x<<"\t"<<M[i]->r.y<<"\t"<<M[i]->r.z<<"\t"<<M[i]->v.x<<"\t"<<M[i]->v.y<<"\t"<<M[i]->v.z<<endl;
    }

}

void environment::MDStep(int i)
{

    if(i%400==0){
        cout<<M[0]->r.x <<"\t\t"<< M[1]->r.x<<"\t\t";
        cout<<"now at step:\t"<<i<<endl;
    }



    if (i % interactionstep == 0)
    {
        // update position:

        if(BoundaryType==1 || BoundaryType==2 ||BoundaryType==3 ||BoundaryType==4 || BoundaryType==6 ) // no force
        {
//            for (int j = 0; j < nsolvent; j++) // update place of solvent particles + Randomshift
//            {

//                M[j]->r.x = (M[j]->r.x) + h*interactionstep * (M[j]->v.x); // Streaming
//                M[j]->r.y = (M[j]->r.y) + h*interactionstep * (M[j]->v.y);
//                M[j]->r.z = (M[j]->r.z) + h*interactionstep * (M[j]->v.z);
//            }

            //                    SPRING(i);
            //            force();
            vector3d zero(0,0,0);
            for(int k=0;k<M.size();k++){
                M[k]->F=zero;
            }
            for(int k=0;k<HeavyNum;k++){
                M[k]->FindF();
            }
//            cout<<M[0]->F.x <<"\t"<< M[1]->F.x<<"\t";
            vector3d externalF=vector3d(NodeMass*1,0,0);
            M[0]->F=M[0]->F+externalF;

            for(int k=0;k<M.size();k++){
                M[k]->Move(h*interactionstep);
            }
            solventBoundaries();
        }

        else if(BoundaryType==5  ) // with force(Velocity verlet)
        {
            for (int j = 0; j < nsolvent; j++) // update place of solvent particles + Randomshift
            {
                (M[j]->r.x) = (M[j]->r.x) + h*interactionstep * (M[j]->v.x); // Streaming
                (M[j]->r.y) = (M[j]->r.y) + h*interactionstep * (M[j]->v.y);
                (M[j]->r.z) = (M[j]->r.z) + h*interactionstep * (M[j]->v.z);
            }


            solventBoundaries();

        }

    }



    if (i % collisiontimestep == 0)// collsion step
    {
        double randomshftvector[3]; // stores the random shift vector in order to inverse it after collision
        // choose random vector

        double TheTa=0;
        double Phi=0;
        TheTa = (double) (rand() % 20000) / 10000-1;  // create random vector  Random coordinates -0.5 +0.5
        Phi = (double) (rand() % 10000) / 10000 ;

        Phi=Phi*2*pi;
        TheTa=asin(TheTa)+pi/2;
        double n1,n2,n3;
        n1=randomshftvector[0] = sin(TheTa)*cos(Phi); // each compenent between -0.5 and 0.5
        n2=randomshftvector[1] = sin(TheTa)*sin(Phi);
        n3=randomshftvector[2] = cos(TheTa);

        if(i==7642||i==5921||i==7652||i==7973){
            if(n1==0||n2==0||n3==0){
                cout<<i<<"        random shift"<<"\t"<<n1<<"\t"<<n2<<"\t"<<n3<<endl;
            }
        }
        for (int j = 0; j < nsolvent; j++) // update place of solvent particles + Randomshift
        {


            //_______________ Randomshift BLOCK :
            (M[j]->r.x)=(M[j]->r.x) +  randomshftvector[0] ;
            (M[j]->r.y)=(M[j]->r.y) +  randomshftvector[1] ;
            (M[j]->r.z)=(M[j]->r.z) +  randomshftvector[2] ;
            //keep solent in BOX
//            if ((M[j]->r.x) > Lx/2.0)
//                (M[j]->r.x) = (M[j]->r.x) - Lx ;
//            if ((M[j]->r.x) < -Lx/2.0)
//                (M[j]->r.x) = (M[j]->r.x) + Lx ;

            if ((M[j]->r.y) > Ly/2.0)
                (M[j]->r.y) = (M[j]->r.y) - Ly ;
            if (M[j]->r.y < -Ly/2.0)
                (M[j]->r.y) = (M[j]->r.y) + Ly ;

            if ((M[j]->r.z) > Lz/2.0)
                ( M[j]->r.z) =( M[j]->r.z) - Lz ;
            if ((M[j]->r.z) < -Lz/2.0)
                (M[j]->r.z) = (M[j]->r.z) + Lz ;
            //_______________ Randomshift BLOCK :
        }


        collision(i); // update velocity in solvent


        for (int j = 0; j < nsolvent; j++) // Randomshift - INVERSE
        {

            //_______________ Randomshift BLOCK  : inversing
            (M[j]->r.x)=(M[j]->r.x) -  randomshftvector[0] ;
            (M[j]->r.y)=(M[j]->r.y) -  randomshftvector[1] ;
            (M[j]->r.z)=(M[j]->r.z) -  randomshftvector[2] ;
            //keep solent in BOX
//            if ((M[j]->r.x) > Lx/2.0)
//                (M[j]->r.x) = (M[j]->r.x) - Lx ;
//            if ((M[j]->r.x) < -Lx/2.0)
//                (M[j]->r.x) = (M[j]->r.x) + Lx ;

            if ((M[j]->r.y) > Ly/2.0)
                (M[j]->r.y) = (M[j]->r.y) - Ly ;
            if ((M[j]->r.y) < -Ly/2.0)
                (M[j]->r.y) = (M[j]->r.y) + Ly ;

            if ((M[j]->r.z) > Lz/2.0)
                (M[j]->r.z) = (M[j]->r.z) - Lz ;
            if ((M[j]->r.z) < -Lz/2.0)
                (M[j]->r.z) = (M[j]->r.z) + Lz ;
            //_______________ Randomshift BLOCK :
        }


    }


    //--------------------------------------------------------------------saving and friends!
}

void environment::saveTrajectory(int time)
{
    com.open(address,fstream::app);
    trajectory.open("trajectory.xyz",fstream::app);
    int i=time;
    int counter= time;
    if (counter%savingstep==0  )  // SAVING AND CALCULATION SECTION
    {

        if (showsolvent==1  )  // SAVING AND CALCULATION SECTION
        {

            trajectory <<2<<endl;    // saving trajectories
            trajectory << " nodes  "<<endl;


            for(int j=0; j<2;j++) // saving trajectory
            {
                trajectory << "s"<<  "\t"<<M[j]->r.x<< "\t" <<M[j]->r.y<< "\t" <<M[j]->r.z<<endl;

            }
            vector3d rcom(0,0,0);
            for(int kkk=0;kkk<2;kkk++){
                rcom=rcom+M[kkk]->r;
            }
//            com<<"test";
            com<<i<<"\t"<<rcom.x/2<<"\t"<<M[0]->v.x<<endl;

        }


        counter=0;


        for(int i=0;i<Lz*spliteachboxfordiagram;i++)
        {
            vdiagram[i][0]=0.0;
            vdiagram[i][1]=0.0;
        }

        for(int i=0;i<nsolvent;i++)
        {
            int index=  ( floor(M[i]->r.z*spliteachboxfordiagram) + spliteachboxfordiagram*Lz/2 );
            vdiagram[index][0]= vdiagram[index][0]+(M[i]->v.x);
            vdiagram[index][1]= vdiagram[index][1]+1;
        }
        for(int i=0;i<Lz*spliteachboxfordiagram;i++)
        {
            if(vdiagram[i][1]!=0)
                vd<<vdiagram[i][0]/vdiagram[i][1]<<endl;
        }
        for(int i=0;i<5;i++)
        {
            vd<<endl;
        }







    }
    counter=counter+1;
    trajectory.close();
    com.close();
}

environment::environment()
{

}

environment::environment(int TryNum)
{
//    string address;
    HeavyNum=2;
    len=1*pow(1.3,TryNum) +TryNum;
    NodeMass=500;
    string num=to_string(TryNum);
    address= "swimmerDualMMD"+num+".txt";
    com.open(address);
//    com<<"tetst";
    LX=Lx;LY=Ly;LZ=Lz;
    com.close();
    trajectory.open("trajectory.xyz");
    trajectory.close();
    initialzesolvent();




}

void environment::collision(int time)
{


    double u[3][nb];

    double n1, n2, n3; // arbitrary vector in space around which particles rotate

    double x, y, z;


    for (int i = 0; i < nb; i++)
    {
        u[0][i] = 0;
        u[1][i] = 0;
        u[2][i] = 0;
    }



    this->updateHowmanySolventsineachBox();



    // calculate velocity of center of mass for each box

    double Mcom[nb]={0};
    for (int i = 0; i < nb; i++)
    {
        for (int j = 0; j < howmanySolventsineachBox[i]; j++)
        {
            int jjj=eachSolventinwichBox[i][j];
            u[0][i] += (M[jjj]->v.x)*(M[jjj]->m);
            u[1][i] += (M[jjj]->v.y)*(M[jjj]->m);
            u[2][i] += (M[jjj]->v.z)*(M[jjj]->m);
            Mcom[i] += (M[jjj]->m);
        }

        if (howmanySolventsineachBox[i] != 0)
        {
            u[0][i] = (double) u[0][i] / Mcom[i];
            u[1][i] = (double) u[1][i] / Mcom[i];
            u[2][i] = (double) u[2][i] / Mcom[i];
        }
        if (howmanySolventsineachBox[i] == 0)
        {

            u[0][i] = 0.0;
            u[1][i] = 0.0;
            u[2][i] = 0.0;
        }
        //        if (howmanySolventsineachBox[i] == 1 )
        //        {
        //            cout<<time;
        //        }


    }



    // rotate velocities in boxes
    for (int i = 0; i < nb; i++)
    {
        //        double VC1=0,VC2=0;



        double TheTa=0;
        double Phi=0;
        TheTa = (double) (rand() % 20000) / 10000-1;  // create random vector  Random coordinates -0.5 +0.5
        Phi = (double) (rand() % 10000) / 10000 ;

        Phi=Phi*2*pi;
        TheTa=asin(TheTa)+pi/2;

        n1 = sin(TheTa)*cos(Phi);
        n2 = sin(TheTa)*sin(Phi);
        n3 = cos(TheTa);

        //        if(time==7642||time==5921||time==7652||time==7973){
        //            if(n1==0||n2==0||n3==0){
        //                cout<<time<<"        Davaran"<<"\t"<<n1<<"\t"<<n2<<"\t"<<n3<<endl;
        //            }
        //        }

        for (int j = 0; j < howmanySolventsineachBox[i]; j++)
        {
            //            int k=j;

            if (eachSolventinwichBox[i][j] >= 0 && eachSolventinwichBox[i][j] < nsolvent)
            {

                x =-u[0][i]+ (M[eachSolventinwichBox[i][j]]->v.x);
                y =-u[1][i]+ (M[eachSolventinwichBox[i][j]]->v.y);
                z =-u[2][i]+ (M[eachSolventinwichBox[i][j]]->v.z);

                //                if(x==0&&y==0&&z==0)
                //                    cout<<time<<endl;

                (M[eachSolventinwichBox[i][j]]->v.x)=u[0][i] +  n1*(n1*x+n2*y+n3*z) + (-n3*y + n2*z );
                (M[eachSolventinwichBox[i][j]]->v.y)=u[1][i] +  n2*(n1*x+n2*y+n3*z) + ( n3*x - n1*z );
                (M[eachSolventinwichBox[i][j]]->v.z)=u[2][i] +  n3*(n1*x+n2*y+n3*z) + (-n2*x + n1*y );


            }


        }


    }

    double Vcom=0;
    for(int j= 0;j<nsolvent;j++){
        Vcom+=M[j]->v.x;
    }
    Vcom/=(nsolvent);

    //    if(fabs(Vcom-Vcom2)>.001)
    //        cout<<time<<"\t"<<Vcom<<"\t"<<Vcom2<<endl;

}

//--------------------------------+UPDATE  howmanySolventsineachBox-----------------

void environment::updateHowmanySolventsineachBox()
{   int x, y, z;
    int index;
    for (int i = 0; i < nb; i++)
    {
        howmanySolventsineachBox[i] = 0;
    }
    // specify the box of each particle
    for (int i = HeavyNum; i < nsolvent; i++)
    {
        if ( ((M[i]->r.x) > (-Lx/2)) && ((M[i]->r.x) < (Lx/2)) && ((M[i]->r.y) > (-Ly/2)) && ((M[i]->r.y) < (Ly/2)) && ((M[i]->r.z) > (-Lz/2)) && ((M[i]->r.z) < Lz/2))
        {
            x =floor ((M[i]->r.x) + (Lx/2));// x y z will lable the corresponding box of each particle
            while(x>Lx)
                x-=Lx;
            while(x<0)
                x+=Lx;
            y =floor ((M[i]->r.y) + (Ly/2));
            z =floor ((M[i]->r.z) + (Lz/2));

            index = (x *  Ly * Lz) + (z * Ly) + y;  // the lable of each box in another way
            eachSolventinwichBox[index][howmanySolventsineachBox[index]] = i;                // says that i th particle is here!
            howmanySolventsineachBox[index]++;                            //says that  this box have how many particles!
        }

    }


}






































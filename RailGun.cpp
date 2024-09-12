#include <iostream>
#include <cmath>
#include <conio.h>
#include "progressbar.hpp"

int main() {
    
    //Mathematical Constants
    const double e = 2.7182;
    const double pi = 3.141592;

    //Time
    double t = 0;

    //Time step
    const long double dt = 0.000000001;

    //Total simulated time
    const long double TargetT = 1;

    //Armature size
    const double ArmWidth = 0.05;
    const double ArmHeight = 0.05;
    const double ArmLength = 0.1;

    //Armature density
    const double Density = 2.710 * pow(10, 3);

    //Calculate armature mass and Frontal Area
    double mass = ArmWidth*ArmHeight*ArmLength*Density;
    double FrontArea = ArmHeight * ArmWidth;

    //Rail size
    const double RailLen = 6;
    const double RW_H = 0.2;


    //Resistance

    //Electrical Resistivity
    const double ResMultip = 0.0265 * pow(10, -6);

    //Total circuit length
    double TotalLen = ArmWidth + 2*RailLen;

    /*
    Area of wires

    Resistance needs more precise calculations
    */
    double Area = ArmHeight*ArmLength;
    double Res;


    //Capacitors

    //Rated voltage
    const double v0 = 400;

    //Capacitance
    const double capa = 16000*pow(10,-6);

    //Current Voltage
    double v = 0;

    //Constants for drag calculation
    const double AirDens = 1.293;
    const double Cd = 1.1;

    //Magnetic Field
    //https://quickfield.com/advanced/biot-savart_law.htm
    const double mu0 = 4 * pi * pow(10, -7);
    double MagDist = ArmWidth / 2;
    double magField;

    //Inductance gradient for railgun force
    double IndGrad = 1;

    //Armature variables
    double acc;
    double vel;
    double dist;

    double I;
    double F_l;
    double F_d;

    //Friction
    double F_r = mass*9.81*0.4;

    double F;

    //Variables for results
    double MaxVel;
    double Maxt;
    double MaxFl;
    double MaxFd;
    double MaxI;
    double MaxMag;

    //Calcultate amount of iterations
    long long int Iterations = (long long int)(TargetT / dt);
    std::cout << "Iterations to go: " << Iterations << '\n';

    //Progressbar
    progressbar bar(100);

    //Loop code with set timesteps until target time is reached
    while(t < TargetT) {

            //Check if projectile is in barrel
            if(dist < RailLen) {
                
                //Calculate electrical resistivity
                Res = (ResMultip*(ArmWidth+(2*dist)))/Area;

                //Calculate voltage
                v = v0 / (pow(e, (t / (Res*capa))));

                //calculate current
                I = v / Res;

                //Calculate magnetic field strength
                magField = 2 * ((mu0 * I) / (2 * pi * MagDist));

                //Calculate Lorentz force
                F_l = magField*I*ArmWidth;

                /*
                With railgun force equation
                F_l = 0.5 * IndGrad * pow(I, 2);

                */
                
                //Save highest current for results, which is at t=0
                if(t == 0) {

                    MaxI = I;
                    MaxMag = magField;

                }



            }
            //if projectile is outside of barrel Lorentz force is 0
            else {
                F_l = 0;
            }
            
            //Calculate drag
            F_d = 0.5*AirDens*Cd*FrontArea*pow(vel, 2);

            //Calculate total force
            F = F_l - F_d - F_r;

            //Calculate acceleration
            acc = F / mass;

            //Calculate new velocity
            vel += acc*dt;

            //Save highest vel for results 
            if(vel > MaxVel) {
                MaxVel = vel;
                Maxt = t;
                MaxFl = F_l;
                MaxFd = F_d;

            }


            //Calculate new distance
            dist += vel*dt;

            //Increase time
            t += dt;

            //check if progressbar should update
            long int CurrentIt = t/dt;

            if((CurrentIt) % (Iterations/100) == 0) {
                bar.update();
            }

    }

    //Output results
    std::cout << '\n';
    std::cout << "Dist: " << dist << '\n';
    std::cout << "Vel: " << vel << '\n';
    std::cout << "Acc: " << acc << '\n';
    std::cout << "F_d: " << F_d << '\n';
    std::cout << "F: " << F << '\n';
    std::cout << "MaxVel: " << MaxVel << '\n';
    std::cout << "Max t: " << Maxt << '\n';
    std::cout << "Max Fl: " << MaxFl << '\n';
    std::cout << "Max Fd: " << MaxFd << '\n';
    std::cout << "Max I: " << MaxI << '\n';
    std::cout << "Max Mag: " << MaxMag << '\n';
    std::cout << "Resis: " << Res << '\n';

    return 0;
}
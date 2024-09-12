#include <iostream>
#include <cmath>
#include <conio.h>
#include "progressbar.hpp"

int main() {

    const double e = 2.7182;
    const double pi = 3.141592;

    //Time
    double t = 0;
    const long double dt = 0.000000001;
    const long double TargetT = 1;

    //Armature
    const double ArmWidth = 0.05;
    const double ArmHeight = 0.05;
    const double ArmLength = 0.1;
    const double Density = 2.710 * pow(10, 3);
    double mass = ArmWidth*ArmHeight*ArmLength*Density;
    double FrontArea = ArmHeight * ArmWidth;

    //Rails
    const double RailLen = 6;
    const double RW_H = 0.2;

    //Resistance
    const double ResMultip = 0.0265 * pow(10, -6);
    double TotalLen = ArmWidth + 2*RailLen;
    double Area = ArmHeight*ArmLength;
    double Res;

    //Capacitors
    const double v0 = 400;
    const double capa = 16000*pow(10,-6);
    double v = 0;

    //Drag
    const double AirDens = 1.293;
    const double Cd = 1.1;

    //Magnetic Field
    //https://quickfield.com/advanced/biot-savart_law.htm
    const double mu0 = 4 * pi * pow(10, -7);
    double MagDist = ArmWidth / 2;
    double magField;

    //Projectile
    double acc;
    double vel;
    double dist;

    double I;
    double F_l;
    double F_d;
    double F_r = mass*9.81*0.4;
    double F;

    double MaxVel;
    double Maxt;
    double MaxFl;
    double MaxFd;
    double MaxI;
    double MaxMag;

    long long int Iterations = (long long int)(TargetT / dt);
    std::cout << "Iterations to go: " << Iterations << '\n';

    //Progressbar
    progressbar bar(100);

    while(t < TargetT) {

            if(dist < RailLen) {

                Res = (ResMultip*(ArmWidth+(2*dist)))/Area;

                v = v0 / (pow(e, (t / (Res*capa))));

                I = v / Res;

                magField = 2 * ((mu0 * I) / (2 * pi * MagDist));

                F_l = magField*I*ArmWidth;

                if(t == 0) {

                    MaxI = I;
                    MaxMag = magField;

                }

            }
            else {
                F_l = 0;
            }
            
            F_d = 0.5*AirDens*Cd*FrontArea*pow(vel, 2);

            F = F_l - F_d - F_r;

            acc = F / mass;

            vel += acc*dt;

            if(vel > MaxVel) {
                MaxVel = vel;
                Maxt = t;
                MaxFl = F_l;
                MaxFd = F_d;

            }

            dist += vel*dt;

            t += dt;

            long int CurrentIt = t/dt;

            if((CurrentIt) % (Iterations/100) == 0) {
                bar.update();
            }

    }

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
#include <iostream>
#include <cmath>
#include <conio.h>
#include "progressbar.hpp"

int main() {

    const double e = 2.7182;
    const double pi = 3.141592;

    //Time
    double t = 0;
    const double dt = 0.00000001;
    const double TargetT = 1;

    //Armature
    const double ArmLen = 0.1;
    const double AW_H = 0.2;
    const double Density = 8.94 * pow(10, 3);
    double mass = ArmLen*pow(AW_H, 2)*Density;
    double FrontArea = AW_H * ArmLen;

    //Rails
    const double RailLen = 3;
    const double RW_H = 0.2;

    //Resistance
    const double ResMultip = 0.0000000175;
    double TotalLen = ArmLen + 2*RailLen;
    double Area = pow(AW_H, 2);
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
    double MagDist = ArmLen / 2;
    double magField;

    //Projectile
    double acc;
    double vel;
    double dist;

    double I;
    double F_l;
    double F_d;
    double F_r = mass*9.81*0.2;
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

                Res = (ResMultip*(ArmLen+(2*dist)))/Area;

                v = v0 / (pow(e, (t / (Res*capa))));

                I = v / Res;

                magField = 2 * ((mu0 * I) / (2 * pi * MagDist));

                F_l = magField*I*ArmLen;

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

    return 0;
}
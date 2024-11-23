#include <iostream>
#include <cmath>
#include <conio.h>
#include <vector>
#include "progressbar.hpp"

//Quickly change standard datatype
using t_simfloat = double;

/*
        Changeable variables
*/

//General simulation settings
const t_simfloat dt_in = 0;         //Timestep inside railgun
const t_simfloat dt_out = 0;        //Timestep outside railgun

//Standard units
const t_simfloat AirDens = 0;
const t_simfloat g = 0;
const t_simfloat RoomTemp = 273;    //In Kelvin

//Friction coefficients between armature, and rails and plates
const t_simfloat mu_s = 0;
const t_simfloat mu_k = 0;

//Properties of individual rails
const t_simfloat l_r = 0;
const t_simfloat w_r = 0;
const t_simfloat h_r = 0;
const t_simfloat dens_r = 0;
const t_simfloat resistiv_r = 0;
const t_simfloat SpecHeat_r = 0;

//Properties of armature
const t_simfloat l_a = 0;
const t_simfloat w_a = 0;
const t_simfloat h_a = 0;
const t_simfloat dens_a = 0;
const t_simfloat resistiv_a = 0;
const t_simfloat SpecHeat_a = 0;
const t_simfloat alpha_V = 0;
const t_simfloat c_w_in = 0;    //c_w changes because the air can't go around as easily inside the railgun
const t_simfloat c_w_out = 0;

//Properties of plates
const t_simfloat h_pl = 0;
const t_simfloat dens_pl = 0;
const t_simfloat SpecHeat_pl = 0;

//Properties of power wires
const t_simfloat l_pw = 0;
const t_simfloat A_pw = 0;
const t_simfloat resistiv_pw = 0;

//Properties of powersupply
const bool ConstPower = false;
const t_simfloat C = 0;
const t_simfloat U0 = 0;




/*
    Variables during simulation
*/

//Current time
t_simfloat t = 0;
t_simfloat dt = 0;

//Positional variables
t_simfloat acc = 0;
t_simfloat speed = 0;
t_simfloat dist = 0;
t_simfloat ddist = 0;

//Vector positional variables
std::vector<t_simfloat> acc_v = {0,0};
std::vector<t_simfloat> vel_v = {0,0};
std::vector<t_simfloat> loc_v = {0,0};

//Electrical variables
t_simfloat R = 0;
t_simfloat U = 0;
t_simfloat I = 0;
t_simfloat E_tot = 0;
t_simfloat E_use = 0;

//Forces
t_simfloat F_d = 0;
t_simfloat F_l = 0;
t_simfloat F_f_pl_u = 0;
t_simfloat F_f_pl_d = 0;
t_simfloat F_f_r = 0;
t_simfloat F_f_tot = 0;
t_simfloat F = 0;

//Coefficients
t_simfloat mu_curr;

//Thermals
t_simfloat Q_pl_u = 0;
t_simfloat Q_pl_d = 0;
t_simfloat Q_r = 0;
t_simfloat Q_a = 0;
t_simfloat Qtot = 0;
t_simfloat T_a = RoomTemp;
t_simfloat T_r = RoomTemp;
t_simfloat T_pl_u = RoomTemp;
t_simfloat T_pl_d = RoomTemp;
t_simfloat dV_a = 0;
t_simfloat P = 0;



/*
    Formulas
*/

//Drag (Air Resistance)
inline t_simfloat calc_F_d(t_simfloat AirDensity, t_simfloat c_w, t_simfloat A, t_simfloat ProjectileSpeed) {
    return ((t_simfloat)0.5) * AirDens * c_w * A * pow(ProjectileSpeed, 2);
}

//Current
inline t_simfloat calc_I(t_simfloat CurrentZero, t_simfloat time, t_simfloat Resistance, t_simfloat Capacity) {
    return CurrentZero * exp(-t / (Resistance*Capacity));
}

// Lorentz force
inline t_simfloat calc_F_l(t_simfloat Current, t_simfloat InductionGradient) {
    return ((t_simfloat)0.5) * InductionGradient * pow(Current, 2);
}

//Friction forces
inline t_simfloat calc_F_f_pl_d(t_simfloat mu, t_simfloat Pressure, t_simfloat Area, t_simfloat Gravity) {
    return mu * ((Pressure*Area) + Gravity);
}
inline t_simfloat calc_F_f_other(t_simfloat mu, t_simfloat Pressure, t_simfloat Area) {
    return mu * (Pressure*Area);
}

//Heat absorbed by object
inline t_simfloat calc_Q_obj(t_simfloat ObjectFriction, t_simfloat deltaDistance) {
    return ((t_simfloat)0.5) * ObjectFriction * deltaDistance;
}

//New object temperature
inline t_simfloat calc_dT_obj(t_simfloat Heat, t_simfloat SpecificHeat, t_simfloat mass) {
    return (Heat / (SpecificHeat * mass));
}

//Armature Volume Increase
inline t_simfloat calc_dV(t_simfloat VolumeZero, t_simfloat VolumetricConstant, t_simfloat Temp) {
    return VolumeZero * VolumetricConstant * (Temp -RoomTemp);
}

//Armature Pressure
inline t_simfloat calc_P(t_simfloat Heat, t_simfloat deltaVolume) {
    return Heat / deltaVolume;
}

//Resistance
inline t_simfloat calc_R_obj(t_simfloat Temp, t_simfloat length, t_simfloat Area) {
    return (((0.0077*Temp) - 0.7175) * pow(10, -8) * length) / Area;
}


int main() {

    /*
        Calculations based on properties
    */

   // Inductance gradient, this should only be computed once, and then cached. Can also be changed to a constant if known
    t_simfloat numerator = pow(10, -6);
    t_simfloat denominator =
        (((t_simfloat)0.5986) * (h_r / w_a)) +
        (((t_simfloat)0.9683) * (h_r / (w_a + ((t_simfloat)2) * w_r))) +
        (((t_simfloat)4.3157) * (((t_simfloat)1) / log(
                (((t_simfloat)4) * (w_a + w_r)) /
                w_r
        ))) -
        ((t_simfloat)0.7831);
    t_simfloat IndGrad = t_simfloat(numerator / denominator);

    //Armature
    t_simfloat Afront_a = w_a * h_a;
    t_simfloat Aside_a = h_a * l_a;
    t_simfloat Atop_a = w_a * l_a;
    t_simfloat V_a = l_a * w_a * h_a;
    t_simfloat V0_a = V_a;
    t_simfloat m_a = V_a * dens_a;
    t_simfloat R_a = (resistiv_a * l_a) / Aside_a;      //For R0, changes during simulating
    t_simfloat F_g = m_a * g;                           //Gravity doesn't change since mass doesn't change

    //Rails
    t_simfloat Afront_r = w_r * h_r;
    t_simfloat V_r = l_r * w_r * h_r;
    t_simfloat m_r = V_r * dens_r;
    t_simfloat R_r = (resistiv_a * l_a) / Afront_a;     //For R0. Starts with armature at front of rails

    //Plates
    t_simfloat V_pl = h_pl * l_r * (w_a + (2 * w_r));
    t_simfloat m_pl = V_pl * dens_pl;

    //Power wires
    t_simfloat R_pw = (resistiv_pw * l_pw) / A_pw;

    //Electrical
    t_simfloat R0 = R_a + R_pw + (2*R_r);       //Current flows through 2 rails
    t_simfloat I0 = U0 / R0;
    R = R0;
    U = U0;
    I = I0;



    /*
        Loops
    */

    //Loop while armature is inside railgun
    while (dist < l_r) {

        //Use current dt
        dt = dt_in;

        //Calculate drag
        F_d = calc_F_d(AirDens, c_w_in, Afront_a, speed);

        //Calculate Lorentzforce
        I = calc_I(I0, t, R, C);
        F_l = calc_F_l(I, IndGrad);

        //First check if static or moving for friction coefficient
        if (speed == 0) {
            mu_curr = mu_s;
        } else {
            mu_curr = mu_k;
        }
        
        //Calculate friction forces
        F_f_pl_d = calc_F_f_pl_d(mu_curr, P, Atop_a, F_g);
        F_f_pl_u = calc_F_f_other(mu_curr, P, Atop_a);
        F_f_r = calc_F_f_other(mu_curr, P, Aside_a);
        F_f_tot = F_f_pl_d + F_f_pl_u + (2*F_f_r);

        //Calculate total force
        F = F_l - F_d - F_f_tot;

        //Calculate movement
        acc = F / m_a;
        speed += acc*dt;
        ddist = speed*dt;
        dist += ddist;

        //Calculate released heat
        E_tot = pow(I, 2) * R * dt;
        E_use = F_l * ddist;
        Qtot += (E_tot - E_use);

        //Calculate heat of objects
        Q_pl_d = calc_Q_obj(F_f_pl_d, ddist);
        Q_pl_u = calc_Q_obj(F_f_pl_u, ddist);
        Q_r = calc_Q_obj(F_f_r, ddist);
        Q_a = Q_pl_d + Q_pl_u + Q_r;

        //Calculate new temperatures
        T_pl_d += calc_dT_obj(Q_pl_d, SpecHeat_pl, m_pl);
        T_pl_u += calc_dT_obj(Q_pl_u, SpecHeat_pl, m_pl);
        T_r += calc_dT_obj(Q_r, SpecHeat_r, m_r);
        T_a += calc_dT_obj(Q_a, SpecHeat_a, m_a);

        //Calculate pressure of armature
        dV_a = calc_dV(V0_a, alpha_V, T_a);
        P = calc_P(Q_a, dV_a);

        //Calculate new resistances
        R_r = calc_R_obj(T_r, (dist + l_a), Afront_r);
        R_a = calc_R_obj(T_a, w_a, Aside_a);

    }

    //Convert to vectors
    vel_v[0] = speed;
    loc_v[0] = dist;

    //Loop for armature outside railgun
    while (speed >= 0) {

        break;

    }


    std::cout << "Done! \n";
    std::cout << "Time: " << t << "\n";
    return 0;
}
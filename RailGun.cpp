#include <iostream>
#include <cmath>
#include <conio.h>
#include <vector>

// Quickly change the datatype used by the simulation
using t_simfloat = long double;

// -----------------------------
// Changeable variables
// -----------------------------
#pragma region constvars

// General simulation settings
const t_simfloat dt_in = 1 * pow(10, -5);        // Timestep inside railgun
const t_simfloat dt_out = 0;                        // Timestep outside railgun

// Standard units
const t_simfloat AirDens = 1.293;   // kg/m3
const t_simfloat g = 9.81;         // TODO g can be more precise than 9.81
const t_simfloat RoomTemp = 293;  // In Kelvin

// Friction coefficients between armature, and rails and plates
const t_simfloat mu_s = 1.5;
const t_simfloat mu_k = 1.1;

// Properties of individual rails
const t_simfloat l_r = 2;
const t_simfloat w_r = 0.05;
const t_simfloat h_r = 0.05;
const t_simfloat dens_r = 8.933 * pow(10, 3);
const t_simfloat resistiv_r = 1.678 * pow(10, -8);
const t_simfloat SpecHeat_r = 384;

// Properties of armature
const t_simfloat l_a = 0.08;
const t_simfloat w_a = 0.05;
const t_simfloat h_a = 0.05;
const t_simfloat dens_a = 8.933 * pow(10, 3);
const t_simfloat resistiv_a = 1.678 * pow(10, -8);
const t_simfloat SpecHeat_a = 384;
const t_simfloat alpha_V = 49.5 * pow(10, -6);
const t_simfloat c_w_in = 1.05;    //c_w changes because the air can't go around as easily inside the railgun
const t_simfloat c_w_out = 1.05;

// Properties of plates
const t_simfloat h_pl = 0.03;
const t_simfloat dens_pl = 8.933 * pow(10, 3);
const t_simfloat SpecHeat_pl =384;

// Properties of power wires
const t_simfloat l_pw = 1;
const t_simfloat A_pw = 0.000314159265359;
const t_simfloat resistiv_pw = 1.678 * pow(10, -8);

// Properties of powersupply
const bool ConstPower = false;
const t_simfloat C = 160000 * pow(10, -6);
const t_simfloat U0 = 400;

#pragma endregion constvars
// -----------------------------
// Variables during simulation
// -----------------------------
#pragma region simvars

// Iterations
t_simfloat it = 0;

// Current time
t_simfloat t = 0;
t_simfloat dt = 0;

// Positional variables
t_simfloat acc = 0;
t_simfloat speed = 0;
t_simfloat dist = 0;
t_simfloat ddist = 0;

// Vector positional variables
std::vector<t_simfloat> acc_v = {0,0};
std::vector<t_simfloat> vel_v = {0,0};
std::vector<t_simfloat> loc_v = {0,0};

// Electrical variables
t_simfloat R = 0;
t_simfloat U = 0;
t_simfloat I = 0;
t_simfloat E_tot = 0;
t_simfloat E_use = 0;

// Forces
t_simfloat F_d = 0;
t_simfloat F_l = 0;
t_simfloat F_f_pl_u = 0;
t_simfloat F_f_pl_d = 0;
t_simfloat F_f_r = 0;
t_simfloat F_f_tot = 0;
t_simfloat F = 0;

// Coefficients
t_simfloat mu_curr;

// Thermals
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

// Debug
t_simfloat F_l0;
t_simfloat speed0;
t_simfloat speedmax = 0;

#pragma endregion simvars
// -----------------------------
// Formulas
// -----------------------------
#pragma region formulas

// Drag (Air Resistance)
inline t_simfloat calc_F_d(t_simfloat c_w, t_simfloat A, t_simfloat ProjectileSpeed) {
    return ((t_simfloat)0.5) * AirDens * c_w * A * pow(ProjectileSpeed, 2);
}

// Current
inline t_simfloat calc_I(t_simfloat CurrentZero, t_simfloat Resistance, t_simfloat Capacity) {
    return CurrentZero * exp(((t_simfloat)-1)* t / (Resistance*Capacity));
}

// Lorentz force
inline t_simfloat calc_F_l(t_simfloat Current, t_simfloat InductionGradient) {
    return ((t_simfloat)0.5) * InductionGradient * pow(Current, 2);
}

// Friction forces
inline t_simfloat calc_F_f_pl_d(t_simfloat mu, t_simfloat Pressure, t_simfloat Area, t_simfloat Gravity) {
    return mu * ((Pressure*Area) + Gravity);
}
inline t_simfloat calc_F_f_other(t_simfloat mu, t_simfloat Pressure, t_simfloat Area) {
    return mu * (Pressure*Area);
}

// Heat absorbed by object
inline t_simfloat calc_Q_obj(t_simfloat ObjectFriction, t_simfloat deltaDistance) {
    return ((t_simfloat)0.5) * ObjectFriction * deltaDistance;
}

// New object temperature
inline t_simfloat calc_dT_obj(t_simfloat Heat, t_simfloat SpecificHeat, t_simfloat mass) {
    return (Heat / (SpecificHeat * mass));
}

// Armature Volume Increase
inline t_simfloat calc_dV(t_simfloat VolumeZero, t_simfloat VolumetricConstant, t_simfloat Temp) {
    return VolumeZero * VolumetricConstant * (Temp - RoomTemp);
}

// Armature Pressure
inline t_simfloat calc_P(t_simfloat Heat, t_simfloat deltaVolume) {
    return Heat / deltaVolume;
}

// Resistance
inline t_simfloat calc_R_obj(t_simfloat Temp, t_simfloat length, t_simfloat Area) {
    return (((0.0077*Temp) - 0.7175) * pow(10, -8) * length) / Area;
}

#pragma endregion formulas
// -----------------------------
// The simulation
// -----------------------------

int main() {

    // -----------------------------
    // Calculations based on properties
    // -----------------------------
    #pragma region constcals

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
    t_simfloat IndGrad = numerator / denominator;
    std::cout << "IndGrad: " << IndGrad << "\n";

    // Armature
    t_simfloat Afront_a = w_a * h_a;
    t_simfloat Aside_a = h_a * l_a;
    t_simfloat Atop_a = w_a * l_a;
    t_simfloat V_a = l_a * w_a * h_a;
    t_simfloat V0_a = V_a;
    t_simfloat m_a = V_a * dens_a;
    t_simfloat R_a = (resistiv_a * l_a) / Aside_a;      //For R0, changes during simulating
    t_simfloat F_g = m_a * g;                           //Gravity doesn't change since mass doesn't change

    // Rails
    t_simfloat Afront_r = w_r * h_r;
    t_simfloat V_r = l_r * w_r * h_r;
    t_simfloat m_r = V_r * dens_r;
    t_simfloat R_r = (resistiv_a * l_a) / Afront_a;     //For R0. Starts with armature at front of rails

    // Plates
    t_simfloat V_pl = h_pl * l_r * (w_a + (2 * w_r));
    t_simfloat m_pl = V_pl * dens_pl;

    // Power wires
    t_simfloat R_pw = (resistiv_pw * l_pw) / A_pw;

    // Electrical
    t_simfloat R0 = R_a + R_pw + (2*R_r);       //Current flows through 2 rails
    t_simfloat I0 = U0 / R0;
    R = R0;
    U = U0;
    I = I0;

    #pragma endregion constcals
    // -----------------------------
    // Simulation loop
    // -----------------------------

    // Loop whilst armature is inside the railgun
    dt = dt_in; // Set the right dt
    while (dist < l_r) {
        //Calculate drag
        F_d = calc_F_d(c_w_in, Afront_a, speed);

        //Calculate Lorentzforce
        I = calc_I(I0, R, C);
        F_l = calc_F_l(I, IndGrad);

        //First check if static or moving for friction coefficient
        if (speed == 0) {
            mu_curr = mu_s;
        } else {
            mu_curr = mu_k;
        }

        //Check if friction forces apply
        if (speed < 0) {
            std::cout << "Speed negative aborting. Iteration: ";
            std::cout << it << "\n";
            std::cout << "Speed: " << speed << "\n";
            return 1;
        } else {
            //Calculate friction forces
            F_f_pl_d = calc_F_f_pl_d(mu_curr, P, Atop_a, F_g);
            F_f_pl_u = calc_F_f_other(mu_curr, P, Atop_a);
            F_f_r = calc_F_f_other(mu_curr, P, Aside_a);
            F_f_tot = F_f_pl_d + F_f_pl_u + (2*F_f_r);
        }

        //Calculate total force
        F = F_l - F_d - F_f_tot;
        //std::cout << it << "\n";
        std::cout << "speed: " << speed << "\n";
        //std::cout << "dist: " << dist << "\n";
        //std::cout << "P: " << P << "\n";
        //std::cout << "F_l: " << F_l << "\n";
        //std::cout << "F_d: " << F_d << "\n";
        //std::cout << "F_f_tot: " << F_f_tot << "\n";
        //std::cout << F << "\n";

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
        P = 0; //calc_P(Q_a, dV_a); // TODO Pressure is way too high

        //Calculate new resistances
        R_r = calc_R_obj(T_r, (dist + l_a), Afront_r);
        R_a = calc_R_obj(T_a, w_a, Aside_a);

        //Increment time
        t += dt;
        it += 1;


        /*
            Debug
        */

        /*if (it == 1) {
            F_l0 = F_l;
            speed0 = speed;
            std::cout << "Rendement: " << E_use/E_tot << "\n";
            std::cout << "F_f0: " << F_f_tot << "\n";
        }
        if (it == 2) {

            std::cout << "Almost last itertation \n";
            std::cout << "Time: " << t << "\n";
            std::cout << "Iterations: " << it << "\n";
            std::cout << "Current: " << I << "\n";
            std::cout << "T_a: " << T_a << "\n";
            std::cout << "acc: " << acc << "\n";
            std::cout << "F_l: " << F_l << "\n";
            std::cout << "F_l0: " << F_l0 << "\n";
            std::cout << "R, U, I" << R0 << " " << U0 << " " << I0 << " " << "\n";
            std::cout << "Dist: " << dist << "\n";
            std::cout << "F: " << F << "\n";
            std::cout << "F_d: " << F_d << "\n";
            std::cout << "F_f_tot: " << F_f_tot << "\n";
            std::cout << "speed: " << speed << "\n";
            std::cout << "speed0: " << speed0 << "\n";
            std::cout << "F_f_r: " << F_f_r << "\n";
            std::cout << "F_f_pl_u: " << F_f_pl_u << "\n";
            std::cout << "F_f_pl_d: " << F_f_pl_d << "\n";
            std::cout << "Pressure: " << P << "\n";
            std::cout << "dV: " << dV_a << "\n";
            std::cout << "Heat: " << Q_a << "\n";
            std::cout << "\n";

        }*/
        //std::cout << it;

        if (speed > speedmax) {
            speedmax = speed;
        }
    }
    return 0;

    t_simfloat E_cap = 0.5*C*pow(U0, 2);
    t_simfloat E_kin = 0.5*m_a*pow(speed, 2);

    t_simfloat AfvuurRendement = E_kin / E_cap;
    std::cout << "Totaal rendement: " << AfvuurRendement << "\n";

    //Convert to vectors
    vel_v[0] = speed;
    loc_v[0] = dist;

    //Loop for armature outside railgun
    while (speed >= 0) {

        break;

    }


    std::cout << "Done! \n";
    std::cout << "Time: " << t << "\n";
    std::cout << "Iterations: " << it << "\n";
    std::cout << "Current: " << I << "\n";
    std::cout << "T_a: " << T_a << "\n";
    std::cout << "acc: " << acc << "\n";
    std::cout << "F_l: " << F_l << "\n";
    std::cout << "F_l0: " << F_l0 << "\n";
    std::cout << "R, U, I" << R0 << " " << U0 << " " << I0 << " " << "\n";
    std::cout << "Dist: " << dist << "\n";
    std::cout << "F: " << F << "\n";
    std::cout << "F_d: " << F_d << "\n";
    std::cout << "F_f_tot: " << F_f_tot << "\n";
    std::cout << "speed: " << speed << "\n";
    std::cout << "speed0: " << speed0 << "\n";
    std::cout << "SpeedMax: " << speedmax << "\n";
    return 0;
}
#include <iostream>
#include <cmath>
#include <conio.h>
#include <vector>
#include <fstream>

// Quickly change the datatype used by the simulation
using t_simfloat = long double;
// Use a slightly less precise float type for the exported data to save disk space
using t_exportfloat = double;

// -----------------------------
// Changeable variables
// -----------------------------
#pragma region constvars

// General simulation settings
constexpr t_simfloat dt_in = 1 * pow(10, -6);        // Timestep inside railgun
constexpr t_simfloat dt_out = 1 * pow(10, -3);       // Timestep outside railgun
constexpr t_simfloat height = 1;

// Standard units
constexpr t_simfloat AirDens = 1.293;  // kg/m3
constexpr t_simfloat g = 9.81;         // TODO g can be more precise than 9.81
constexpr t_simfloat RoomTemp = 293;   // In Kelvin

// Friction coefficients between armature, and rails and plates
constexpr t_simfloat mu_s = 1.5;
constexpr t_simfloat mu_k = 1.1;

// Properties of individual rails
constexpr t_simfloat l_r = 6;
constexpr t_simfloat w_r = 0.04;
constexpr t_simfloat h_r = 0.06;
constexpr t_simfloat dens_r = 8.933 * pow(10, 3);
constexpr t_simfloat resistiv_r = 1.678 * pow(10, -8);
constexpr t_simfloat SpecHeat_r = 383.9;

// Properties of armature
constexpr t_simfloat l_a = 0.0114;
constexpr t_simfloat w_a = 0.0114;
constexpr t_simfloat h_a = 0.056;
constexpr t_simfloat dens_a = 8.933 * pow(10, 3);
constexpr t_simfloat resistiv_a = 1.678 * pow(10, -8);
constexpr t_simfloat SpecHeat_a = 383.9;
constexpr t_simfloat k_T = 140 * pow(10,9);             //Bulk modulus
constexpr t_simfloat alpha_V = 49.5 * pow(10, -6);
constexpr t_simfloat c_w_in = 2;    //c_w changes because the air can't go around as easily inside the railgun
constexpr t_simfloat c_w_out = 1.05;

// Properties of plates
constexpr t_simfloat h_pl = 0.03;
constexpr t_simfloat dens_pl = 8.933 * pow(10, 3);
constexpr t_simfloat SpecHeat_pl =384;

// Properties of power wires
constexpr t_simfloat l_pw = 3;
constexpr t_simfloat A_pw = 0.000314159265359;
constexpr t_simfloat resistiv_pw = 1.678 * pow(10, -8);

// Properties of powersupply
constexpr bool ConstPower = false;
constexpr t_simfloat C = 0.9363885766;
constexpr t_simfloat U0 = 4000;
constexpr t_simfloat ConstR0 = 0.00468 * pow(10, 0) * .5666666666666666666666667;
constexpr t_simfloat ConstRho = 1.678 * pow(10, -8);
constexpr t_simfloat dUdt = -1172400;

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
std::vector<t_simfloat> F_d_v = {0,0};
std::vector<t_simfloat> ddist_v = {0,0};

// Electrical variables
t_simfloat R = 0;
t_simfloat U = 0;
t_simfloat I = 0;
t_simfloat E_tot = 0;
t_simfloat E_use = 0;
t_simfloat ConstR = ConstR0;

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
t_simfloat Q_a_tot = 0;
t_simfloat Qtot = 0;
t_simfloat T_a = RoomTemp;
t_simfloat T_r = RoomTemp;
t_simfloat T_pl_u = RoomTemp;
t_simfloat T_pl_d = RoomTemp;
t_simfloat dV_a = 0;
t_simfloat P = 0;

// Misc
t_simfloat F_l0;
t_simfloat speed0;
t_simfloat speedmax;
t_simfloat ExitVelocity;
t_simfloat CapEnergy;
t_simfloat ExitI;
t_simfloat Discharged = false;
t_simfloat DischargeTime;
t_simfloat KineticEnergy;
t_simfloat RestEnergy;

#pragma endregion simvars
// -----------------------------
// Formulas
// -----------------------------
#pragma region formulas

// Drag (Air Resistance)
inline t_simfloat calc_F_d(t_simfloat c_w, t_simfloat A, t_simfloat ProjectileSpeed) {
    return 0.5 * AirDens * c_w * A * (ProjectileSpeed * ProjectileSpeed);
}

// Current
inline t_simfloat calc_I(t_simfloat CurrentZero, t_simfloat Resistance, t_simfloat Capacity) {
    return CurrentZero * expl((-1 * t) / (Resistance*Capacity));
}

inline t_simfloat calc_I_custom(t_simfloat Cap, t_simfloat dU) {
    return Cap * dU;
}

// Lorentz force
inline t_simfloat calc_F_l(t_simfloat Current, t_simfloat InductionGradient) {
    return 0.5 * InductionGradient * (Current * Current);
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
    return 0.5 * ObjectFriction * deltaDistance;
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
inline t_simfloat calc_P_eqs(t_simfloat deltaV, t_simfloat Volume0) {
    return k_T * (deltaV/Volume0);
}

// Resistance
inline t_simfloat calc_R_obj(t_simfloat Temp, t_simfloat length, t_simfloat Area) {
    return (((0.0077*Temp) - 0.7175)            // Formula for resistivity at higher temperatures 
            * pow(10, -8) * length) / Area;    
}

// Const R Temp Dependency
inline t_simfloat calc_R_const(t_simfloat Res0, t_simfloat Temp) {
    return Res0 * ((((0.0077*Temp) - 0.7175) * pow(10, -8))/(ConstRho));
}

#pragma endregion formulas
// -----------------------------
// Internal defs
// -----------------------------
#pragma region internaldefs
// Define the struct of a single export data point, fields are written to the export buffer in the same order as how they are defined here
struct Datapoint {
    t_exportfloat t;
    t_exportfloat speed_x;
    t_exportfloat speed_y;
    t_exportfloat loc_x;
    t_exportfloat loc_y;
    t_exportfloat F_l;
    t_exportfloat I;
    t_exportfloat T_a;
    t_exportfloat R_a;
    t_exportfloat P;
    t_exportfloat R_r;
    t_exportfloat R;
    t_exportfloat ConstR;
    t_exportfloat Q_tot;
};
#pragma endregion internaldefs

// We keep a buffer with datapoints that we flush to the disk once it's full. This allows us to step more than once without a disk write, reducing IO overhead.
constexpr size_t DATAPOINTS_BUF_COUNT = 1 << 20; // 2 to the power 20

// -----------------------------
// The simulation
// -----------------------------

int main() {
    // Python requires IEEE floats, on most systems t_exportfloat should be IEEE. We log if it is just to make sure.
    std::cout << "Is IEEE float? " << std::numeric_limits<t_exportfloat>::is_iec559 << "\n";
    std::cout << "Float size: " << sizeof(t_exportfloat) << "\n";
    std::cout << "Buffer size: " << (sizeof(Datapoint) * DATAPOINTS_BUF_COUNT) / 1024 / 1024 << " MiB\n";

    // -----------------------------
    // Calculations based on properties
    // -----------------------------
    #pragma region constcals

    // Inductance gradient, this should only be computed once, and then cached. Can also be changed to a constant if known
    t_simfloat numerator = pow(10, -6);
    t_simfloat denominator =
        (static_cast<t_simfloat>(0.5986) * (h_r / w_a)) +
        (static_cast<t_simfloat>(0.9683) * (h_r / (w_a + ((t_simfloat)2) * w_r))) +
        (static_cast<t_simfloat>(4.3157) * (static_cast<t_simfloat>(1) / log(
            (static_cast<t_simfloat>(4) * (w_a + w_r)) /
            w_r
        ))) -
        static_cast<t_simfloat>(0.7831);
    t_simfloat IndGrad = numerator / denominator;
    std::cout << "IndGrad: " << IndGrad << "\n";

    // Armature
    t_simfloat Afront_a = w_a * h_a;
    t_simfloat Aside_a = h_a * l_a;
    t_simfloat Atop_a = w_a * l_a;
    t_simfloat V_a = l_a * w_a * h_a;
    t_simfloat V0_a = V_a;
    t_simfloat m_a = 0.740;//V_a * dens_a;
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
    t_simfloat R0 = R_a + R_pw + (2*R_r) + ConstR;       //Current flows through 2 rails
    t_simfloat I0 = U0 / R0;
    R = R0;
    U = U0;
    I = I0;

    #pragma endregion constcals
    // -----------------------------
    // Initialize export buffer
    // -----------------------------
    // Allocate our buffer, this will be flushed to disk every time it is full
    auto* datapoints = static_cast<Datapoint*>(malloc(
        sizeof(Datapoint) * DATAPOINTS_BUF_COUNT
    ));
    int datapoint_idx = 0;
    std::ofstream out_file;
    // Open out.bin, which is where we will write our simulation results to
    out_file.open("out.bin", std::ios::out | std::ios::trunc | std::ios::binary);

    // Debug
    CapEnergy = 0.5 * C * pow(U0, 2);

    // -----------------------------
    // Simulation loop 1: Inside the railgun
    // -----------------------------
    dt = dt_in; // Set the right dt
    while (dist < l_r) {

        //Calculate drag
        F_d = calc_F_d(c_w_in, Afront_a, speed);

        //Calculate Lorentzforce
        I = calc_I(I0, R, C);
        //I = calc_I_custom(C, dUdt);
        F_l = calc_F_l(I, IndGrad);

        //First check if static or moving for friction coefficient
        if (speed == 0) {
            mu_curr = mu_s;
        } else {
            mu_curr = mu_k;
        }

        if (speed < 0) {
            break;
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
        E_tot = (I * I) * R * dt;
        E_use = F_l * ddist;
        Qtot += (E_tot - E_use);

        //Calculate heat of objects
        Q_r = (m_r / (m_a + (2*m_r))) * (E_tot - E_use);
        Q_a = (m_a / (m_a + (2*m_r))) * (E_tot - E_use); 
        Q_a_tot += Q_a;

        //Calculate new temperatures
        T_pl_d += calc_dT_obj(Q_pl_d, SpecHeat_pl, m_pl);
        T_pl_u += calc_dT_obj(Q_pl_u, SpecHeat_pl, m_pl);
        T_r += calc_dT_obj(Q_r, SpecHeat_r, m_r);
        T_a += calc_dT_obj(Q_a, SpecHeat_a, m_a);

        //Calculate pressure of armature
        dV_a = calc_dV(V0_a, alpha_V, T_a);
        V_a += dV_a;
        P = 0 * calc_P_eqs(dV_a, V0_a);

        //Calculate new resistances
        R_r = calc_R_obj(T_r, (dist + l_a), Afront_r);
        R_a = calc_R_obj(T_a, w_a, Aside_a);
        ConstR = calc_R_const(ConstR0, T_a);

        R = R_a + R_pw + (2*R_r) + ConstR;

        //Increment time
        t += dt;
        it += 1;

        // Store datapoint
        Datapoint* datapoint = &datapoints[datapoint_idx]; // Get a reference to the first datapoint in the buffer that has not yet been written
        datapoint_idx++; // Increment our buffer index, so that the next iteration of the loop doesn't overwrite this datapoint
        // Fill in the datapoint fields with the values from the current step
        datapoint->t = t;
        datapoint->speed_x = speed;
        datapoint->speed_y = 0;
        datapoint->loc_x = dist;
        datapoint->loc_y = height;
        datapoint->F_l = F_l;
        datapoint->I = I;
        datapoint->T_a = T_a;
        datapoint->R_a = R_a;
        datapoint->P = P;
        datapoint->R_r = R_r;
        datapoint->R = R;
        datapoint->ConstR = ConstR;
        datapoint->Q_tot = Qtot;
        // Flush the buffer to the disk if it is full
        if (datapoint_idx == DATAPOINTS_BUF_COUNT) {
            // Flush the buffer
            out_file.write(reinterpret_cast<char*>(datapoints), sizeof(Datapoint) * DATAPOINTS_BUF_COUNT);
            datapoint_idx = 0; // Reset the buffer index back to the start
        }

        // Store topspeed and discharge time (Time until 1% of initial current)
        if (speed > speedmax) {
            speedmax = speed;
        }
        if (I < (I0*0.01) && Discharged == false) {
            DischargeTime = t;
            Discharged = true;
        }
    }

    //Calculate efficiency
    t_simfloat E_cap = 0.5*C*powl(U0, 2);
    t_simfloat E_kin = 0.5*m_a*powl(speed, 2);

    t_simfloat Efficiency = E_kin / E_cap;
    std::cout << "Total Efficiency: " << Efficiency << "\n";

    //Save exit velocity
    ExitVelocity = speed;
    ExitI = I;
    KineticEnergy = 0.5 * m_a * powl(speed, 2);
    Qtot = SpecHeat_a * (m_a + (2*m_r)) * (T_a - RoomTemp);

    RestEnergy = CapEnergy - (KineticEnergy + Qtot);


    //Convert to vectors
    vel_v[0] = speed;
    loc_v[0] = dist;
    loc_v[1] = height;

///*
    // -----------------------------
    // Simulation loop 2: Outside the railgun
    // -----------------------------
    dt = dt_out;
    while (false && speed >= 1 && loc_v[1] > 0) {
        // Calculate drag in both directions
        F_d_v[0] = calc_F_d(c_w_out, Afront_a, vel_v[0]);
        F_d_v[1] = calc_F_d(c_w_out, Atop_a, vel_v[1]);

        // Calculate acceleration in both directions
        acc_v[0] = -F_d_v[0] / m_a;
        acc_v[1] = (-F_g + F_d_v[1]) / m_a;

        // Calculate new velocity
        vel_v[0] += acc_v.at(0)*dt;
        vel_v[1] += acc_v.at(1)*dt;

        // Calculate delta distance
        ddist_v[0] = vel_v[0]*dt;
        ddist_v[1] = vel_v[1]*dt;

        // Calculate new location
        loc_v[0] += ddist_v[0];
        loc_v[1] += ddist_v[1];

        // New distance
        dist += sqrtl(powl(ddist_v[0], 2L) + powl(ddist_v[1], 2));

        // Speed
        speed = sqrtl(powl(vel_v[0], 2) + powl(vel_v[1], 2));

        t += dt;
        it += 1;

        // Store datapoint
        Datapoint* datapoint = &datapoints[datapoint_idx]; // Get a reference to the first datapoint in the buffer that has not yet been written
        datapoint_idx++; // Increment our buffer index, so that the next iteration of the loop doesn't overwrite this datapoint
        // Fill in the datapoint fields with the values from the current step
        datapoint->t = t;
        datapoint->speed_x = vel_v[0];
        datapoint->speed_y = vel_v[1];
        datapoint->loc_x = loc_v[0];
        datapoint->loc_y = loc_v[1];
        datapoint->F_l = 0;
        datapoint->I = I;
        datapoint->T_a = T_a;
        datapoint->R_a = R_a;
        datapoint->P = P;
        datapoint->R_r = R_r;
        // Flush the buffer to the disk if it is full
        if (datapoint_idx == DATAPOINTS_BUF_COUNT) {
            // Flush the buffer
            out_file.write(reinterpret_cast<char*>(datapoints), sizeof(Datapoint) * DATAPOINTS_BUF_COUNT);
            datapoint_idx = 0;
        }
    }
//*/

    // -----------------------------
    // Finalize and cleanup
    // -----------------------------
    // Flush and free the export buffer
    if (datapoint_idx > 0) { // If the buffer is not empty, make sure to write the bits of the buffer that have been filled to the disk.
        // Flush unflushed buffer
        out_file.write(reinterpret_cast<char*>(datapoints), sizeof(Datapoint) * (datapoint_idx - 1));
    }
    // Close out.bin
    out_file.close();
    free(datapoints); // De-allocate the buffer in memory, to prevent memory leaks.

    // Print stats
    std::cout << "Done! \n";
    std::cout << "Time: " << t << "\n";
    std::cout << "Iterations: " << it << "\n";
    std::cout << "Current: " << I << "\n";
    std::cout << "T_a: " << T_a << "\n";
    std::cout << "T_r: " << T_r << "\n";
    std::cout << "m_a: " << m_a << "\n";
    std::cout << "acc: " << acc << "\n";
    std::cout << "F_l: " << F_l << "\n";
    std::cout << "F_l0: " << F_l0 << "\n";
    std::cout << "R, U, I" << R0 << " " << U0 << " " << I0 << " " << "\n";
    std::cout << "CapEnergy: " << CapEnergy << "\n";
    std::cout << "Dist: " << dist << "\n";
    std::cout << "F: " << F << "\n";
    std::cout << "F_d: " << F_d << "\n";
    std::cout << "F_f_tot: " << F_f_tot << "\n";
    std::cout << "speed: " << speed << "\n";
    std::cout << "X Velocity: " << vel_v[0] << "\n";
    std::cout << "Y Velocity: " << vel_v[1] << "\n";
    std::cout << "y drag: " << F_d_v[1] << "\n";
    std::cout << "speed0: " << speed0 << "\n";
    std::cout << "SpeedMax: " << speedmax << "\n";
    std::cout << "Exit Velocity: " << ExitVelocity << "\n";
    std::cout << "Q_tot: " << Qtot << "\n";
    std::cout << "m_tot: " << m_a + (2*m_r) << "\n";
    std::cout << "V0_a: " << V0_a << "\n";
    std::cout << "V_a: " << V_a << "\n";
    std::cout << "ExitI: " << ExitI << "\n";
    std::cout << "Q_a_tot: " << Q_a_tot << "\n"; 
    std::cout << "m_r: " << m_r << "\n";
    std::cout << "DischargeTime: " << DischargeTime << "\n";
    std::cout << "Kinetic Energy: " << KineticEnergy << "\n";
    std::cout << "RestEnergy: " << RestEnergy << "\n";
    return 0;
}

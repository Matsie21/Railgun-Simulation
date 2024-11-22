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
const t_simfloat dt_in = 0;     //Timestep inside railgun
const t_simfloat dt_out = 0;    //Timestep outside railgun

//Standard units
const t_simfloat AirDens = 0;
const t_simfloat g = 0;

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
t_simfloat F_z = 0;
t_simfloat F = 0;

//Thermals
t_simfloat Q_pl_u = 0;
t_simfloat Q_pl_d = 0;
t_simfloat Q_r = 0;
t_simfloat Qtot = 0;
t_simfloat T_a = 0;
t_simfloat T_r = 0;
t_simfloat T_pl_u = 0;
t_simfloat T_pl_d = 0;
t_simfloat dV_a = 0;
t_simfloat P = 0;



/*
    Formulas
*/

//Drag (Air Resistance)
inline t_simfloat F_drag(t_simfloat AirDensity, t_simfloat c_w, t_simfloat A, t_simfloat ProjectileSpeed) {
    return ((t_simfloat)0.5) * AirDens * c_w * A * pow(ProjectileSpeed, 2);
}

// Inductance gradient, this should only be computed once, and then cached
t_simfloat numerator = pow(10, -6);
t_simfloat denominator =
    (((t_simfloat)0.5986) * (h_r / w_a)) +
    (((t_simfloat)0.9683) * (h_r / (w_a + ((t_simfloat)2) * w_r))) +
    (((t_simfloat)4.3157) * (((t_simfloat)1) / logf(
            (((t_simfloat)4) * (w_a + w_r)) /
            w_r
    ))) -
    ((t_simfloat)0.7831);
t_simfloat IndGrad = t_simfloat(numerator / denominator);

int main() {

    /*
        Calculations based on properties
    */

    //Armature
    t_simfloat Afront_a = w_a * h_a;
    t_simfloat Aside_a = h_a * l_a;
    t_simfloat Atop_a = w_a * l_a;
    t_simfloat V_a = l_a * w_a * h_a;
    t_simfloat m_a = V_a * dens_a;
    t_simfloat R_a = (resistiv_a * l_a) / Aside_a;      //For R0, changes during simulating

    //Rails
    t_simfloat Afront_r = w_r * h_r;
    t_simfloat V_r = l_r * w_r * h_r;
    t_simfloat m_r = V_r * dens_r;
    t_simfloat R_r = (resistiv_r * l_r) / Afront_r;     //For R0, changes during simulating

    //Plates
    t_simfloat V_pl = h_pl * l_r * (w_a + (2 * w_r));
    t_simfloat m_pl = V_pl * dens_pl;

    //Power wires
    t_simfloat R_pw = (resistiv_pw * l_pw) / A_pw;

    //Electrical
    t_simfloat R0 = R_a + R_pw + R_r;
    t_simfloat I0 = U0 / R0;
    R = R0;
    U = U0;
    I = I0;


    //Loop while armature is inside railgun
    while (dist < l_r) {

        //

    }

    //Convert to vectors
    vel_v[0] = speed;
    loc_v[0] = dist;

    //Loop for armature outside railgun
    while (speed >= 0) {

        //

    }

    return 0;
}
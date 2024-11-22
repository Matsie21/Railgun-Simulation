#include <iostream>
#include <cmath>
#include <conio.h>
#include "progressbar.hpp"

//Quickly change standard datatype
using t_simfloat = double;

/*
        Changeable variables
*/

//General simulation settings
const t_simfloat t = 0;
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
const t_simfloat c_w_in = 0;
const t_simfloat c_w_out = 0;

//Properties of plates
const t_simfloat h_pl = 0;
const t_simfloat SpecHeat_pl = 0;

//Properties of power wires
const t_simfloat l_pw = 0;
const t_simfloat A_pw = 0;
const t_simfloat resistiv_pw = 0;

//Properties of powersupply
bool ConstPower = false;
const t_simfloat C = 0;
const t_simfloat U0 = 0;



int main() {

    /*
        Calculations based on properties
    */

   //Armature
   

    return 0;
}
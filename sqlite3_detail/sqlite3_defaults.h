//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_DEFAULTS_H
#define LSSEFT_SQLITE3_DEFAULTS_H


constexpr auto SQLITE3_DEFAULT_FRW_TABLE_NAME                 = "models";
constexpr auto SQLITE3_DEFAULT_REDSHIFT_CONFIGURATION_TABLE   = "z_config";
constexpr auto SQLITE3_DEFAULT_WAVENUMBER_CONFIGURATION_TABLE = "k_config";
constexpr auto SQLITE3_DEFAULT_IR_CUTOFF_CONFIGURATION_TABLE  = "IR_cutoff";
constexpr auto SQLITE3_DEFAULT_UV_CUTOFF_CONFIGURATION_TABLE  = "UV_cutoff";
constexpr auto SQLITE3_DEFAULT_IR_RESUM_CONFIGURATION_TABLE   = "IR_resum";
constexpr auto SQLITE3_DEFAULT_TRANSFER_TABLE                 = "transfer";
constexpr auto SQLITE3_DEFAULT_GROWTH_G_TABLE                 = "g_factors";
constexpr auto SQLITE3_DEFAULT_GROWTH_F_TABLE                 = "f_factors";
constexpr auto SQLITE3_DEFAULT_LOOP_FG_META_TABLE             = "fg_metadata";
constexpr auto SQLITE3_DEFAULT_LOOP_AA_TABLE                  = "delta22_AA";
constexpr auto SQLITE3_DEFAULT_LOOP_AB_TABLE                  = "delta22_AB";
constexpr auto SQLITE3_DEFAULT_LOOP_BB_TABLE                  = "delta22_BB";
constexpr auto SQLITE3_DEFAULT_LOOP_D_TABLE                   = "delta13_D";
constexpr auto SQLITE3_DEFAULT_LOOP_E_TABLE                   = "delta13_E";
constexpr auto SQLITE3_DEFAULT_LOOP_F_TABLE                   = "delta13_F";
constexpr auto SQLITE3_DEFAULT_LOOP_G_TABLE                   = "delta13_G";
constexpr auto SQLITE3_DEFAULT_LOOP_J1_TABLE                  = "delta13_J1";
constexpr auto SQLITE3_DEFAULT_LOOP_J2_TABLE                  = "delta13_J2";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD13_A_TABLE             = "rsd13_a";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD13_B_TABLE             = "rsd13_b";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD13_C_TABLE             = "rsd13_c";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD13_D_TABLE             = "rsd13_d";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD13_E_TABLE             = "rsd13_e";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD13_F_TABLE             = "rsd13_f";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD13_G_TABLE             = "rsd13_g";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_A1_TABLE            = "rsd22_A1";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_A2_TABLE            = "rsd22_A2";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_A3_TABLE            = "rsd22_A3";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_A4_TABLE            = "rsd22_A4";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_A5_TABLE            = "rsd22_A5";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_B2_TABLE            = "rsd22_B2";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_B3_TABLE            = "rsd22_B3";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_B6_TABLE            = "rsd22_B6";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_B8_TABLE            = "rsd22_B8";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_B9_TABLE            = "rsd22_B9";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_C1_TABLE            = "rsd22_C1";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_C2_TABLE            = "rsd22_C2";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_C4_TABLE            = "rsd22_C4";
constexpr auto SQLITE3_DEFAULT_LOOP_RSD22_D1_TABLE            = "rsd22_D1";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_TABLE                  = "dd_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_RSD_MU0_TABLE          = "dd_rsd_mu0_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_RSD_MU2_TABLE          = "dd_rsd_mu2_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_RSD_MU4_TABLE          = "dd_rsd_mu4_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_RSD_MU6_TABLE          = "dd_rsd_mu6_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_RSD_MU8_TABLE          = "dd_rsd_mu8_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_P0_TABLE                  = "dd_P0";
constexpr auto SQLITE3_DEFAULT_LOOP_P2_TABLE                  = "dd_P2";
constexpr auto SQLITE3_DEFAULT_LOOP_P4_TABLE                  = "dd_P4";

constexpr auto SQLITE3_DEFAULT_TEMPORARY_TABLE                = "temp";


#endif //LSSEFT_SQLITE3_DEFAULTS_H

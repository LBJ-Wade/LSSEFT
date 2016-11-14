//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_DEFAULTS_H
#define LSSEFT_SQLITE3_DEFAULTS_H


constexpr auto SQLITE3_DEFAULT_FRW_TABLE_NAME                 = "models";
constexpr auto SQLITE3_DEFAULT_REDSHIFT_CONFIGURATION_TABLE   = "z_config";
constexpr auto SQLITE3_DEFAULT_WAVENUMBER_CONFIGURATION_TABLE = "k_config";
constexpr auto SQLITE3_DEFAULT_IR_CONFIGURATION_TABLE         = "IR_cutoff";
constexpr auto SQLITE3_DEFAULT_UV_CONFIGURATION_TABLE         = "UV_cutoff";
constexpr auto SQLITE3_DEFAULT_TRANSFER_TABLE                 = "transfer";
constexpr auto SQLITE3_DEFAULT_GROWTH_G_TABLE                 = "g_factors";
constexpr auto SQLITE3_DEFAULT_GROWTH_F_TABLE                 = "f_factors";
constexpr auto SQLITE3_DEFAULT_LOOP_FG_META_TABLE             = "fg_metadata";
constexpr auto SQLITE3_DEFAULT_LOOP_AA_TABLE                  = "AA";
constexpr auto SQLITE3_DEFAULT_LOOP_AB_TABLE                  = "AB";
constexpr auto SQLITE3_DEFAULT_LOOP_BB_TABLE                  = "BB";
constexpr auto SQLITE3_DEFAULT_LOOP_D_TABLE                   = "D";
constexpr auto SQLITE3_DEFAULT_LOOP_E_TABLE                   = "E";
constexpr auto SQLITE3_DEFAULT_LOOP_F_TABLE                   = "F";
constexpr auto SQLITE3_DEFAULT_LOOP_G_TABLE                   = "G";
constexpr auto SQLITE3_DEFAULT_LOOP_J1_TABLE                  = "J1";
constexpr auto SQLITE3_DEFAULT_LOOP_J2_TABLE                  = "J2";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_TABLE                  = "dd_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_RSD_MU0_TABLE          = "dd_rsd_mu0_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_RSD_MU2_TABLE          = "dd_rsd_mu2_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_RSD_MU4_TABLE          = "dd_rsd_mu4_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_RSD_MU6_TABLE          = "dd_rsd_mu6_Pk";
constexpr auto SQLITE3_DEFAULT_LOOP_DD_RSD_MU8_TABLE          = "dd_rsd_mu8_Pk";

constexpr auto SQLITE3_DEFAULT_TEMPORARY_TABLE                = "temp";


#endif //LSSEFT_SQLITE3_DEFAULTS_H

#include "embedana.h"

void embedana::InitParams()
{
    BBC_VERTEX_CUT = 10;
    TRIGGERBIT = 4; // 4 - BBCLL1(>1 tubes) narrowvtx, 5 - BBCLL1(>1 tubes), 6 - BBCLL1(>1 tubes) novertex
    E_PT = 0.2;     // 0.15 in RDanalyzer for simulation to study smearing effect
    MAX_PT = 5;
    QUALITY[0] = 31;
    QUALITY[1] = 51;
    QUALITY[2] = 63;
    Z_GLOBAL = 90;
    DC_DEADMAP = 1; // by default, DC deadmap is on
    N0 = 2;
    DISP = 5;
    CHI2_NPE0 = 10;
    EOVERP = -999;
    DEP[0] = -2; // -2, can be applied offline
    DEP[1] =  5;  // 5
    PROB = 0.01; // 0.01
    TEMC = 5;
    EMCDPHI = 0.05; // ~5 sigma 0.05
    EMCDZ = 25;     // ~5 sigma 25
    EMCSDPHI = 5;   // remove sigmalized cut for Run14AuAu because the recalibrator doesn't work for simulation
    EMCSDZ = 5;
    RICH_GHOST = 10; // by default, RICH ghost cut is off

    std::cout << "Run14AuAuCut loaded... " << std::endl;
    std::cout << "***Event selection*** " << std::endl;
    std::cout <<  "BBC VERTEX CUT = "<<BBC_VERTEX_CUT<<std::endl;
    if (TRIGGERBIT > -99)
        std::cout << "trigger bit = " << TRIGGERBIT << std::endl;
    std::cout << "ZDC conincidence required" << std::endl;

    std::cout << "***Single track selection***" << std::endl;
    if (QUALITY[0] > -99)
        std::cout << "quality = " << QUALITY[0] << std::endl;
    if (QUALITY[1] > -99)
        std::cout << "quality = " << QUALITY[1] << std::endl;
    if (QUALITY[2] > -99)
        std::cout << "quality = " << QUALITY[2] << std::endl;
    if (Z_GLOBAL > -99)
        std::cout << "|z_global| < " << Z_GLOBAL << std::endl;
    if (DC_DEADMAP)
    {
        std::cout << "dc dead map in" << std::endl;
        init_rungroup();
        init_dead_area();
    }

    std::cout << "***e+/- track selection***" << std::endl;
    if (E_PT > -99)
        std::cout << "pT > " << E_PT << std::endl;
    if (MAX_PT > -99)
        std::cout << "pT < " << MAX_PT << std::endl;
    if (N0 > -99)
        std::cout << "n0 > " << N0 << std::endl;
    if (DISP > -99)
        std::cout << "disp < " << DISP << std::endl;
    if (CHI2_NPE0 > -99)
        std::cout << "chi2/npe0 < " << CHI2_NPE0 << std::endl;
    if (EMCDPHI > -99)
        std::cout << "emc matching cut in: emcdphi < " << EMCDPHI << "rad" << std::endl;
    if (EMCDZ > -99)
        std::cout << "emc matching cut in: emcdz < " << EMCDZ << "cm" << std::endl;
    if (EMCSDPHI > -99)
        std::cout << "emc matching cut (sigmalized) in: emcsdphi_e < " << EMCSDPHI << "sigma" << std::endl;
    if (EMCSDZ > -99)
        std::cout << "emc matching cut (sigmalized) in: emcsdz_e < " << EMCSDZ << "sigma" << std::endl;
    if (PROB > -99)
        std::cout << "prob > " << PROB << std::endl;
    if (TEMC > -99)
        std::cout << "TEMC < " << TEMC << std::endl;
    if (EOVERP > -99)
        std::cout << "E/p > " << EOVERP << std::endl;
    if (DEP[0] > -99)
        std::cout << "dep > " << DEP[0] << std::endl;
    if (DEP[1] > -99)
        std::cout << "dep < " << DEP[1] << std::endl;
    if (RICH_GHOST > -99)
    {
        std::cout << "RICH ghost cut in: decenter < " << RICH_GHOST << "cm" << std::endl;
    }

    if(!DC_DEADMAP)
    {
        std::cout << "setting init_dead_area to -9999" << std::endl;
        for (int i = 0; i < N_RUN_GRP; ++i)
        {
            for (int j = 0; j < N_SIDE; ++j)
            {
                for (int k = 0; k < N_ARM; ++k)
                {
                    for (int l = 0; l < MAX_DEAD_AREA; ++l)
                    {
                        dcmap_xx1[i][j][k][l] = -9999;
                        dcmap_yy1[i][j][k][l] = -9999;
                        dcmap_xx2[i][j][k][l] = -9999;
                        dcmap_yy2[i][j][k][l] = -9999;
                        dcmap_xx3[i][j][k][l] = -9999;
                        dcmap_yy3[i][j][k][l] = -9999;
                        dcmap_xx4[i][j][k][l] = -9999;
                        dcmap_yy4[i][j][k][l] = -9999;
                    }
                }
            }
        }
    }
}

float embedana::get_board(float phi_dc)
{
    if (phi_dc > 1.57)
        return ((3.72402 - phi_dc + 0.008047 * cos(phi_dc + 0.87851)) / 0.01963496);
    return ((0.573231 + phi_dc - 0.0046 * cos(phi_dc + 0.05721)) / 0.01963496);
}

void embedana::init_rungroup()
{
    for (int igrp = 0; igrp < N_RUN_GRP; ++igrp)
    {
        for (int irun = 0; irun < MAX_RUN; ++irun)
        {
            dcmap_runs[igrp][irun] = -9999;
        }
    }

    // run group 1
    int igrp = 1;
    for (int irun = 0; irun < N_dcmap_runs_rg1; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg1[irun];
    }

    // run group 2
    igrp = 2;
    for (int irun = 0; irun < N_dcmap_runs_rg2; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg2[irun];
    }

    // run group 3
    igrp = 3;
    for (int irun = 0; irun < N_dcmap_runs_rg3; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg3[irun];
    }

    // run group 4
    igrp = 4;
    for (int irun = 0; irun < N_dcmap_runs_rg4; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg4[irun];
    }

    // run group 5
    igrp = 5;
    for (int irun = 0; irun < N_dcmap_runs_rg5; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg5[irun];
    }

    // run group 6
    igrp = 6;
    for (int irun = 0; irun < N_dcmap_runs_rg6; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg6[irun];
    }

    // run group 7
    igrp = 7;
    for (int irun = 0; irun < N_dcmap_runs_rg7; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg7[irun];
    }

    // run group 8
    igrp = 8;
    for (int irun = 0; irun < N_dcmap_runs_rg8; ++irun)
    {
        dcmap_runs[igrp][irun] = dcmap_runs_rg8[irun];
    }
}

int embedana::get_rungroup(int run_num)
{
    for (int igrp = 0; igrp < N_RUN_GRP; ++igrp)
    {
        for (int irun = 0; irun < MAX_RUN; ++irun)
        {
            if (run_num == dcmap_runs[igrp][irun])
                return igrp;
        }
    }
    return 0; // no match (bad runs)
}

void embedana::init_dead_area()
{

    // run group 1
    int igrp = 1;
    for (int iarea = 0; iarea < (int)(N_ES_rg1 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg1[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg1[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg1[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg1[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg1[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg1[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg1[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg1[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg1 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg1[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg1[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg1[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg1[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg1[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg1[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg1[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg1[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg1 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg1[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg1[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg1[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg1[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg1[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg1[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg1[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg1[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg1 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg1[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg1[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg1[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg1[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg1[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg1[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg1[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg1[4 * iarea + 3];
    }

    // run group 2
    igrp = 2;
    for (int iarea = 0; iarea < (int)(N_ES_rg2 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg2[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg2[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg2[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg2[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg2[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg2[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg2[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg2[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg2 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg2[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg2[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg2[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg2[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg2[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg2[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg2[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg2[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg2 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg2[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg2[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg2[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg2[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg2[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg2[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg2[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg2[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg2 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg2[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg2[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg2[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg2[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg2[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg2[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg2[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg2[4 * iarea + 3];
    }

    // run group 3
    igrp = 3;
    for (int iarea = 0; iarea < (int)(N_ES_rg3 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg3[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg3[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg3[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg3[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg3[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg3[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg3[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg3[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg3 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg3[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg3[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg3[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg3[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg3[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg3[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg3[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg3[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg3 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg3[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg3[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg3[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg3[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg3[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg3[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg3[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg3[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg3 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg3[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg3[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg3[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg3[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg3[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg3[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg3[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg3[4 * iarea + 3];
    }

    // run group 4
    igrp = 4;
    for (int iarea = 0; iarea < (int)(N_ES_rg4 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg4[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg4[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg4[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg4[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg4[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg4[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg4[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg4[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg4 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg4[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg4[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg4[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg4[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg4[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg4[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg4[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg4[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg4 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg4[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg4[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg4[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg4[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg4[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg4[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg4[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg4[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg4 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg4[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg4[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg4[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg4[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg4[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg4[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg4[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg4[4 * iarea + 3];
    }

    // run group 5
    igrp = 5;
    for (int iarea = 0; iarea < (int)(N_ES_rg5 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg5[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg5[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg5[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg5[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg5[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg5[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg5[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg5[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg5 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg5[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg5[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg5[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg5[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg5[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg5[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg5[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg5[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg5 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg5[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg5[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg5[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg5[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg5[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg5[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg5[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg5[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg5 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg5[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg5[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg5[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg5[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg5[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg5[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg5[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg5[4 * iarea + 3];
    }

    // run group 6
    igrp = 6;
    for (int iarea = 0; iarea < (int)(N_ES_rg6 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg6[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg6[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg6[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg6[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg6[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg6[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg6[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg6[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg6 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg6[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg6[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg6[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg6[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg6[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg6[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg6[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg6[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg6 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg6[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg6[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg6[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg6[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg6[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg6[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg6[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg6[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg6 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg6[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg6[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg6[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg6[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg6[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg6[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg6[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg6[4 * iarea + 3];
    }

    // run group 7
    igrp = 7;
    for (int iarea = 0; iarea < (int)(N_ES_rg7 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg7[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg7[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg7[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg7[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg7[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg7[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg7[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg7[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg7 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg7[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg7[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg7[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg7[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg7[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg7[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg7[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg7[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg7 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg7[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg7[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg7[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg7[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg7[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg7[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg7[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg7[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg7 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg7[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg7[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg7[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg7[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg7[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg7[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg7[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg7[4 * iarea + 3];
    }

    // run group 8
    igrp = 8;
    for (int iarea = 0; iarea < (int)(N_ES_rg8 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][0][iarea] = dcmap_x_ES_rg8[4 * iarea];
        dcmap_yy1[igrp][0][0][iarea] = dcmap_y_ES_rg8[4 * iarea];
        dcmap_xx2[igrp][0][0][iarea] = dcmap_x_ES_rg8[4 * iarea + 1];
        dcmap_yy2[igrp][0][0][iarea] = dcmap_y_ES_rg8[4 * iarea + 1];
        dcmap_xx3[igrp][0][0][iarea] = dcmap_x_ES_rg8[4 * iarea + 2];
        dcmap_yy3[igrp][0][0][iarea] = dcmap_y_ES_rg8[4 * iarea + 2];
        dcmap_xx4[igrp][0][0][iarea] = dcmap_x_ES_rg8[4 * iarea + 3];
        dcmap_yy4[igrp][0][0][iarea] = dcmap_y_ES_rg8[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WS_rg8 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][0][1][iarea] = dcmap_x_WS_rg8[4 * iarea];
        dcmap_yy1[igrp][0][1][iarea] = dcmap_y_WS_rg8[4 * iarea];
        dcmap_xx2[igrp][0][1][iarea] = dcmap_x_WS_rg8[4 * iarea + 1];
        dcmap_yy2[igrp][0][1][iarea] = dcmap_y_WS_rg8[4 * iarea + 1];
        dcmap_xx3[igrp][0][1][iarea] = dcmap_x_WS_rg8[4 * iarea + 2];
        dcmap_yy3[igrp][0][1][iarea] = dcmap_y_WS_rg8[4 * iarea + 2];
        dcmap_xx4[igrp][0][1][iarea] = dcmap_x_WS_rg8[4 * iarea + 3];
        dcmap_yy4[igrp][0][1][iarea] = dcmap_y_WS_rg8[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_EN_rg8 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][0][iarea] = dcmap_x_EN_rg8[4 * iarea];
        dcmap_yy1[igrp][1][0][iarea] = dcmap_y_EN_rg8[4 * iarea];
        dcmap_xx2[igrp][1][0][iarea] = dcmap_x_EN_rg8[4 * iarea + 1];
        dcmap_yy2[igrp][1][0][iarea] = dcmap_y_EN_rg8[4 * iarea + 1];
        dcmap_xx3[igrp][1][0][iarea] = dcmap_x_EN_rg8[4 * iarea + 2];
        dcmap_yy3[igrp][1][0][iarea] = dcmap_y_EN_rg8[4 * iarea + 2];
        dcmap_xx4[igrp][1][0][iarea] = dcmap_x_EN_rg8[4 * iarea + 3];
        dcmap_yy4[igrp][1][0][iarea] = dcmap_y_EN_rg8[4 * iarea + 3];
    }
    for (int iarea = 0; iarea < (int)(N_WN_rg8 / 4.); ++iarea)
    {
        dcmap_xx1[igrp][1][1][iarea] = dcmap_x_WN_rg8[4 * iarea];
        dcmap_yy1[igrp][1][1][iarea] = dcmap_y_WN_rg8[4 * iarea];
        dcmap_xx2[igrp][1][1][iarea] = dcmap_x_WN_rg8[4 * iarea + 1];
        dcmap_yy2[igrp][1][1][iarea] = dcmap_y_WN_rg8[4 * iarea + 1];
        dcmap_xx3[igrp][1][1][iarea] = dcmap_x_WN_rg8[4 * iarea + 2];
        dcmap_yy3[igrp][1][1][iarea] = dcmap_y_WN_rg8[4 * iarea + 2];
        dcmap_xx4[igrp][1][1][iarea] = dcmap_x_WN_rg8[4 * iarea + 3];
        dcmap_yy4[igrp][1][1][iarea] = dcmap_y_WN_rg8[4 * iarea + 3];
    }
}

bool embedana::left_bound(float x, float y, float xx1, float yy1, float xx2, float yy2)
{ // bound line defined by (xx1, yy1), (xx2, yy2), check if (x, y) is on the right-hand side of the bound line
    if (xx1 == xx2 && x > xx1)
        return true;
    else
    {
        const float k = (yy2 - yy1) / (xx2 - xx1); // slope
        if (k > 0 && y < k * (x - xx2) + yy2)
            return true;
        if (k < 0 && y > k * (x - xx2) + yy2)
            return true;
    }
    return false;
}

bool embedana::right_bound(float x, float y, float xx1, float yy1, float xx2, float yy2)
{ // need to also consider slope infinity case (xx1==xx2)
    if (xx1 == xx2 && x < xx1)
        return true;
    else
    {
        const float k = (yy2 - yy1) / (xx2 - xx1); // slope
        if (k > 0 && y > k * (x - xx2) + yy2)
            return true;
        if (k < 0 && y < k * (x - xx2) + yy2)
            return true;
    }
    return false;
}

bool embedana::dead_region(float x, float y, float xx1, float yy1, float xx2, float yy2, float xx3, float yy3, float xx4, float yy4)
{ // if in dead region, return true
    // alpha cut is equivalent to pt, so the middle region better be kept to keep the large pt
    if (xx1 < -999)
        return false; // not a dead region
    if (xx1 > -999 && xx1 < 0 && right_bound(x, y, xx3, yy3, xx4, yy4))
        return true; // only right bound
    if (xx3 > -999 && xx3 < 0 && left_bound(x, y, xx1, yy1, xx2, yy2))
        return true; // only left bound
    if (left_bound(x, y, xx1, yy1, xx2, yy2) && right_bound(x, y, xx3, yy3, xx4, yy4))
        return true;
    return false;
}

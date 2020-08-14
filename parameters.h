#pragma once

// parameters that will be dynamically initialized within readArguments function
extern double gjj, gjk, E, B, A;
extern double quenched_gjj, quenched_gjk, quenched_E, quenched_B, quenched_A;
extern int Jtot, M_J, Mtot;
extern int START_IT, END_IT;
extern double step;
extern int QNmode, vsWhat;
extern double PMode;

extern bool ONLY_ENERGIES;
extern bool MAGNETIZATION;
extern bool J2;
extern bool PRINT_EIGENVECTORS;

extern bool QUENCH_DYNAMICS;
extern bool QUENCH_R2, QUENCH_MG, QUENCH_J, QUENCH_JTOT;
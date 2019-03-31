#ifndef TESTS_H
#define TESTS_H

#include <chrono>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TRandom3.h>
#include <TRandom1.h>
#include "PolynomialFit.h"
#include "FunctionTable.h"
#include "argon_cross.h"

void test_polynomial_fit (void);
void test_2_dim_table (void);
void test_phase_shift_fit (ArDataTables *ArTables);
void test_factor_helper (void);
void test_legendre_polynomial(void);
void test_legendre_intregral (void);
void test_colored_interval (void);
void test_diff_tot_cross (ArDataTables *ArTables);
void test_J_L_probabilities (ArDataTables *ArTables);
void test_time_delay(ArDataTables *ArTables);
void test_backward_scatter_prob (ArDataTables *ArTables);
void test_TM_forward (ArDataTables *ArTables);
void test_TM_backward (ArDataTables *ArTables);
void test_total_cross_all (ArDataTables *ArTables);
void test_resonance_NBrS (ArDataTables *ArTables); //TODO:
void test_all(ArDataTables *ArTables);

#endif

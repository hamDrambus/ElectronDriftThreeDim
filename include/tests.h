#ifndef TESTS_H
#define TESTS_H

#include <chrono>
#include <boost/random.hpp>
#ifndef _NO_CERN_ROOT
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TRandom3.h>
#include <TRandom1.h>
#endif //_NO_CERN_ROOT
#include "PolynomialFit.h"
#include "FunctionTable.h"
#include "ParticleTable.h"

void test_polynomial_fit (void);
void test_2_dim_table (void);
void test_phase_shift_fit (void);
void test_factor_helper (void);
void test_legendre_polynomial(void);
void test_legendre_intregral (void);
void test_colored_interval (void);
void test_diff_tot_cross (void);
void test_J_L_probabilities (void);
void test_backward_scatter_prob (void);
void test_TM_forward (void);
void test_TM_backward (void);
void test_total_cross_all (void);
void test_time_delay(void);
//void test_resonance_NBrS (ArDataTables *ArTables); //TODO:
void test_all(void);

#endif

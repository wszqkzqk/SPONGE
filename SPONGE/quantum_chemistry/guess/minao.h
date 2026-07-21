#pragma once

#include "../structure/molecule.h"
#include "../structure/scf_workspace.h"

// Build a MINAO-like diagonal density from neutral-atom angular populations.
// ECP core electrons are removed from the real atom configuration, then a
// bounded global projection enforces the molecule's charge and spin targets.
void QC_Build_Minao_Guess(const QC_MOLECULE& mol,
                          const QC_SCF_Runtime_State& runtime, float* d_P,
                          float* d_P_beta);

#pragma once

#include "../Lennard_Jones_force/LJ_soft_core.h"
#include "../Lennard_Jones_force/Lennard_Jones_force.h"
#include "../common.h"
#include "../control.h"

enum SITS_MODE
{
    SITS_MODE_OBSERVATION = 0,
    SITS_MODE_ITERATION = 1,
    SITS_MODE_PRODUCTION = 2,
    SITS_MODE_EMPIRICAL = 4,
    SITS_MODE_AMD = 5,
    SITS_MODE_GAMD = 6
};

struct SELECT
{
    std::vector<float*> select_atom_energy;
    std::vector<float*> select_energy;
    std::vector<VECTOR*> select_force;
    std::vector<LTMatrix3*> select_atom_virial_tensor;
    std::vector<LTMatrix3*> select_virial_tensor;

    void Initial();

    int Add_One_Energy(int atom_numbers);

    int Add_One_Force(int atom_numbers);

    int Add_One_Virial(int atom_numbers);
};

struct CLASSIC_SITS_INFORMATION;
struct SITS_INFORMATION;

struct SITS_STATE_SNAPSHOT
{
    float enhancing_energy = 0.0f;
    float bias = 0.0f;
    float factor = 1.0f;
    float fb_reference_energy = 0.0f;
    float fb_reference_bias = 0.0f;
    float effective_potential = 0.0f;
};

struct CLASSIC_SITS_INFORMATION
{
    int is_initialized = 0;

    SITS_INFORMATION* sits_controller;
    int k_numbers;
    float* beta_k;
    float* NkExpBetakU;
    float* Nk;
    float* sum_a;
    float* sum_b;
    float* factor;
    float* fb_reference_energy;
    float* fb_reference_bias;

    int record_count;
    int record_interval;
    int update_interval;
    int fb_interval;
    int reset = 1;

    // ene_recorded - vshift - ene
    // gf - gf - log( n_k * exp(-beta_k * ene) )
    // gfsum - gfsum - log( Sum_(k=1)^N ( log( n_k * exp(-beta_k * ene) ) ) )
    // log_weight - rb - log of the weighting function
    // log_mk_inverse - ratio - log(m_k^-1)
    // log_norm_old - normlold - W(j-1)
    // log_norm - norml - W(j)
    // log_pk - rbfb - log(p_k)
    // log_nk_inverse - pratio - log(n_k^-1)
    // log_nk - fb - log(n_k)
    float* d_bias;
    float* ene_recorded;
    float* gf;
    float* gfsum;
    float* log_weight;
    float* log_mk_inverse;
    float* log_norm_old;
    float* log_norm;
    float* log_pk;
    float* log_nk_inverse;
    float* log_nk;

    float pe_a = 1.0;
    float pe_b = 0.0;
    float fb_bias = 0.0;
    float T_high = 600;
    float T_low = 250;
    int nk_fix = 0;

    // record
    FILE* nk_traj_file = NULL;
    FILE* nk_rest_file = NULL;

    std::string nk_rest_file_name;
    std::string nk_traj_file_name;

    float* nk_record_cpu;
    float* log_norm_record_cpu;

    void Initial(CONTROLLER* controller, SITS_INFORMATION* sits);

    void Memory_Allocate();

    void SITS_Record_Ene();

    void SITS_Update_gf();

    void SITS_Update_gfsum();

    void SITS_Update_log_pk();

    void SITS_Update_log_mk_inverse();

    void SITS_Update_log_nk_inverse();

    void SITS_Update_nk();

    void SITS_Update_Common(const float beta);

    void SITS_Update_Nk();

    void SITS_Write_Nk_Norm();

    void SITS_Update_Fb(float beta_0, int step);
};

struct SITS_INFORMATION
{
    char module_name[CHAR_LENGTH_MAX];
    char print_ab_name[CHAR_LENGTH_MAX];
    char print_aa_name[CHAR_LENGTH_MAX];
    char print_aa_kab_name[CHAR_LENGTH_MAX];
    char print_bias_name[CHAR_LENGTH_MAX];
    char print_fb_name[CHAR_LENGTH_MAX];
    int is_initialized = 0;
    CONTROLLER* controller = NULL;
    int sits_mode = 0;
    int atom_numbers;
    int* atom_sys_mark;  // 标记溶剂，0为体系，1为溶剂
    int* atom_sys_mark_local = NULL;
    int local_atom_numbers = 0;
    int ghost_numbers = 0;
    int* d_local_metadata_error = NULL;
    bool local_metadata_is_ready = false;

    float pwwp_enhance_factor = 0.5;
    float h_factor = 1.0;

    float h_enhancing_energy;

    SELECT pw_select;

    bool selectively_applied;
    CLASSIC_SITS_INFORMATION classic_sits;

    void Initial(CONTROLLER* controller, int atom_numbers_,
                 const char* module_name = NULL);
    void Memory_Allocate();

    void Reset_Force_Energy(int* md_need_potential);
    void Update_And_Enhance(const int step, float* d_total_potential,
                            const int need_pressure, LTMatrix3* d_total_virial,
                            VECTOR* frc, float beta0,
                            bool update_statistics = true);
    void Save_State(SITS_STATE_SNAPSHOT* snapshot,
                    const float* d_effective_potential);
    void Restore_State(const SITS_STATE_SNAPSHOT& snapshot,
                       float* d_effective_potential);

    void SITS_LJ_Direct_CF_Force_With_Atom_Energy_And_Virial(
        const int atom_numbers, const int local_atom_numbers,
        const int solvent_numbers, const int ghost_numbers, const VECTOR* crd,
        const float* charge, LENNARD_JONES_INFORMATION* lj_info, VECTOR* md_frc,
        const LTMatrix3 cell, const LTMatrix3 rcell, const ATOM_GROUP* nl,
        const float cutoff, const float pme_beta, const int need_energy,
        float* atom_energy_ww, const int need_pressure,
        LTMatrix3* atom_virial_ww, float* elect_atom_ene);

    void SITS_LJ_Soft_Core_Direct_CF_Force_With_Atom_Energy_And_Virial(
        const int atom_numbers, const int local_atom_numbers,
        const int solvent_numbers, const int ghost_numbers, const VECTOR* crd,
        const float* charge, LJ_SOFT_CORE* lj_info, VECTOR* frc,
        const LTMatrix3 cell, const LTMatrix3 rcell, const ATOM_GROUP* nl,
        const float cutoff, const float pme_beta, const int need_energy,
        float* atom_energy_ww, const int need_pressure,
        LTMatrix3* atom_virial_ww, float* elect_atom_ene);

    void Step_Print(CONTROLLER* controller, const float beta0);

    void Check_Solvent(CONTROLLER* controller, int atom_numbers,
                       int solvent_numbers);

    /*
        以下用于区域分解
    */
    void Get_Local(int* atom_local, int local_atom_numbers, int ghost_numbers);
    bool Validate_Local_State(const char* error_by, int global_atom_numbers,
                              int local_atom_numbers, int ghost_numbers);
};

#pragma once

struct MINIMIZATION_iteration
{
    MD_INFORMATION* md_info = NULL;
    CONTROLLER* controller = NULL;
    int* d_invalid_atom = NULL;
    float max_move = 0.02f;
    int dynamic_dt = 1;
    int last_decrease_step = 0;
    float momentum_keep = 0;
    float beta1 = 0.9f;
    float beta2 = 0.99f;
    float epsilon = 1e-4f;
    float learning_rate = 3e-4f;
    void Gradient_Descent(int atom_numbers, const int* atom_local, VECTOR* crd,
                          const VECTOR* frc, VECTOR* vel,
                          const float* d_mass_inverse, const VECTOR* adam_move);
    void Scale_Force_For_Dynamic_Dt(int atom_numbers, const int* atom_local,
                                    const float* d_mass_inverse,
                                    const VECTOR* frc, VECTOR* first_moment,
                                    VECTOR* root_second_moment,
                                    VECTOR* adam_move);
    void Initial(CONTROLLER* controller, MD_INFORMATION* md_info);
    void Validate_Final_Time_Step();
    void Clear();
};

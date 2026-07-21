#include "tersoff.h"

enum TersoffParam
{
    p_m = 0,
    p_gamma,
    p_lam3,
    p_c,
    p_d,
    p_h,
    p_n,
    p_beta,
    p_lam2,
    p_B,
    p_R,
    p_D,
    p_lam1,
    p_A,
    p_c1,
    p_c2,
    p_c3,
    p_c4,
    PARAM_STRIDE = 18
};

static __device__ __forceinline__ float ters_fc(float r, const float* param)
{
    float R = param[p_R];
    float D = param[p_D];
    if (r < R - D) return 1.0f;
    if (r > R + D) return 0.0f;
    return 0.5f * (1.0f - sinf(CONSTANT_Pi * 0.5f * (r - R) / D));
}

template <int N>
static __device__ __forceinline__ SADfloat<N> ters_fc_sad(SADfloat<N> r,
                                                          const float* param)
{
    float R = param[p_R];
    float D = param[p_D];
    if (r.val < R - D) return SADfloat<N>(1.0f);
    if (r.val > R + D) return SADfloat<N>(0.0f);
    return 0.5f * (1.0f - sinf(CONSTANT_Pi * 0.5f * (r - R) / D));
}

template <int N>
static __device__ __forceinline__ SADfloat<N> ters_gijk_sad(
    SADfloat<N> costheta, const float* param)
{
    float c = param[p_c];
    float d = param[p_d];
    float h = param[p_h];
    float gamma = param[p_gamma];
    SADfloat<N> diff = h - costheta;
    SADfloat<N> term = d * d + diff * diff;
    return gamma * (1.0f + c * c / (d * d) - c * c / term);
}

static __device__ __forceinline__ float ters_bij(float zeta, const float* param)
{
    float beta = param[p_beta];
    float n = param[p_n];
    float c1 = param[p_c1];
    float c2 = param[p_c2];
    float c3 = param[p_c3];
    float c4 = param[p_c4];
    float tmp = beta * zeta;
    if (tmp > c1) return 1.0f / sqrtf(tmp);
    if (tmp > c2) return (1.0f - powf(tmp, -n) / (2.0f * n)) / sqrtf(tmp);
    if (tmp < c4) return 1.0f;
    if (tmp < c3) return 1.0f - powf(tmp, n) / (2.0f * n);
    return powf(1.0f + powf(tmp, n), -1.0f / (2.0f * n));
}

static __device__ __forceinline__ float ters_bij_d(float zeta,
                                                   const float* param)
{
    float beta = param[p_beta];
    float n = param[p_n];
    float c1 = param[p_c1];
    float c2 = param[p_c2];
    float c3 = param[p_c3];
    float c4 = param[p_c4];
    float tmp = beta * zeta;
    if (tmp > c1) return beta * -0.5f * powf(tmp, -1.5f);
    if (tmp > c2)
        return beta * (-0.5f * powf(tmp, -1.5f) *
                       (1.0f - (1.0f + 1.0f / (2.0f * n)) * powf(tmp, -n)));
    if (tmp < c4) return 0.0f;
    if (tmp < c3) return -0.5f * beta * powf(tmp, n - 1.0f);
    float tmp_n = powf(tmp, n);
    return -0.5f * powf(1.0f + tmp_n, -1.0f - (1.0f / (2.0f * n))) * tmp_n /
           zeta;
}

template <bool need_energy, bool need_virial>
static __global__ void Tersoff_Force_CUDA(
    const int atom_numbers, const VECTOR* crd, VECTOR* frc,
    const LTMatrix3 cell, const LTMatrix3 rcell, const ATOM_GROUP* nl,
    int* atom_types, float* params, int* map, int ntypes, float* atom_energy,
    LTMatrix3* atom_virial, float* d_energy_sum)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        int type_i = atom_types[i];
        VECTOR ri = crd[i];
        ATOM_GROUP nl_i = nl[i];
        VECTOR fi = {0, 0, 0};
        float en_i = 0;
        LTMatrix3 vi = {0, 0, 0, 0, 0, 0};

        for (int jj = 0; jj < nl_i.atom_numbers; jj++)
        {
            int j = nl_i.atom_serial[jj];
            int type_j = atom_types[j];
            VECTOR rj = crd[j];
            VECTOR drij = Get_Periodic_Displacement(ri, rj, cell, rcell);
            float rij = norm3df(drij.x, drij.y, drij.z);
            int param_idx_ij =
                map[type_i * ntypes * ntypes + type_j * ntypes + type_j];
            const float* param_ij = params + param_idx_ij * PARAM_STRIDE;
            if (rij > param_ij[p_R] + param_ij[p_D]) continue;

            SADfloat<1> rij_sad(rij, 0);
            SADfloat<1> fr = param_ij[p_A] * expf(-param_ij[p_lam1] * rij_sad) *
                             ters_fc_sad(rij_sad, param_ij);

            float f_rep = -0.5f * fr.dval[0] / rij;
            fi.x += f_rep * drij.x;
            fi.y += f_rep * drij.y;
            fi.z += f_rep * drij.z;
            atomicAdd(&frc[j].x, -f_rep * drij.x);
            atomicAdd(&frc[j].y, -f_rep * drij.y);
            atomicAdd(&frc[j].z, -f_rep * drij.z);
            if (need_energy)
            {
                float ev = 0.5f * fr.val;
                en_i += ev;
                atomicAdd(d_energy_sum, ev);
            }
            if (need_virial)
            {
                vi = vi + Get_Virial_From_Force_Dis(f_rep * drij, drij);
            }

            float zeta = 0;
            for (int kk = 0; kk < nl_i.atom_numbers; kk++)
            {
                if (jj == kk) continue;
                int k = nl_i.atom_serial[kk];
                int type_k = atom_types[k];
                int param_idx_ijk =
                    map[type_i * ntypes * ntypes + type_j * ntypes + type_k];
                const float* param_ijk = params + param_idx_ijk * PARAM_STRIDE;
                VECTOR rk = crd[k];
                VECTOR drik = Get_Periodic_Displacement(ri, rk, cell, rcell);
                float rik = norm3df(drik.x, drik.y, drik.z);
                if (rik > param_ijk[p_R] + param_ijk[p_D]) continue;
                float costheta =
                    (drij.x * drik.x + drij.y * drik.y + drij.z * drik.z) /
                    (rij * rik);
                float fc_ik = ters_fc(rik, param_ijk);
                float g = param_ijk[p_gamma] *
                          (1.0f +
                           param_ijk[p_c] * param_ijk[p_c] /
                               (param_ijk[p_d] * param_ijk[p_d]) -
                           param_ijk[p_c] * param_ijk[p_c] /
                               (param_ijk[p_d] * param_ijk[p_d] +
                                (param_ijk[p_h] - costheta) *
                                    (param_ijk[p_h] - costheta)));
                float diff = rij - rik;
                float arg = param_ijk[p_lam3] * diff;
                if (fabsf(param_ijk[p_m] - 3.0f) < 1e-5) arg = arg * arg * arg;
                zeta += fc_ik * g * expf(arg);
            }

            float bij = ters_bij(zeta, param_ij);
            float fa_val = -param_ij[p_B] * expf(-param_ij[p_lam2] * rij) *
                           ters_fc(rij, param_ij);
            SADfloat<1> fa_sad = -param_ij[p_B] *
                                 expf(-param_ij[p_lam2] * rij_sad) *
                                 ters_fc_sad(rij_sad, param_ij);
            float f_attr_direct = -0.5f * bij * fa_sad.dval[0] / rij;
            fi.x += f_attr_direct * drij.x;
            fi.y += f_attr_direct * drij.y;
            fi.z += f_attr_direct * drij.z;
            atomicAdd(&frc[j].x, -f_attr_direct * drij.x);
            atomicAdd(&frc[j].y, -f_attr_direct * drij.y);
            atomicAdd(&frc[j].z, -f_attr_direct * drij.z);
            if (need_energy)
            {
                float ev = 0.5f * bij * fa_val;
                en_i += ev;
                atomicAdd(d_energy_sum, ev);
            }
            if (need_virial)
            {
                vi = vi + Get_Virial_From_Force_Dis(f_attr_direct * drij, drij);
            }

            float pre = -0.5f * fa_val * ters_bij_d(zeta, param_ij);
            for (int kk = 0; kk < nl_i.atom_numbers; kk++)
            {
                if (jj == kk) continue;
                int k = nl_i.atom_serial[kk];
                int param_idx_ijk = map[type_i * ntypes * ntypes +
                                        type_j * ntypes + atom_types[k]];
                const float* param_ijk = params + param_idx_ijk * PARAM_STRIDE;
                VECTOR rk = crd[k];
                VECTOR drik = Get_Periodic_Displacement(ri, rk, cell, rcell);
                float rik = norm3df(drik.x, drik.y, drik.z);
                if (rik > param_ijk[p_R] + param_ijk[p_D]) continue;
                float costheta =
                    (drij.x * drik.x + drij.y * drik.y + drij.z * drik.z) /
                    (rij * rik);
                SADfloat<3> s_rij(rij, 0);
                SADfloat<3> s_rik(rik, 1);
                SADfloat<3> s_cos(costheta, 2);
                SADfloat<3> s_arg = param_ijk[p_lam3] * (s_rij - s_rik);
                if (fabsf(param_ijk[p_m] - 3.0f) < 1e-5)
                    s_arg = s_arg * s_arg * s_arg;
                SADfloat<3> s_zeta_ijk = ters_fc_sad(s_rik, param_ijk) *
                                         ters_gijk_sad(s_cos, param_ijk) *
                                         expf(s_arg);

                VECTOR rj_hat = (1.0f / rij) * drij;
                VECTOR rk_hat = (1.0f / rik) * drik;
                VECTOR dcos_dxj = (-1.0f / rij) * (rk_hat - costheta * rj_hat);
                VECTOR dcos_dxk = (-1.0f / rik) * (rj_hat - costheta * rk_hat);

                VECTOR fj_tri;
                fj_tri.x = pre * (-s_zeta_ijk.dval[0] * rj_hat.x +
                                  s_zeta_ijk.dval[2] * dcos_dxj.x);
                fj_tri.y = pre * (-s_zeta_ijk.dval[0] * rj_hat.y +
                                  s_zeta_ijk.dval[2] * dcos_dxj.y);
                fj_tri.z = pre * (-s_zeta_ijk.dval[0] * rj_hat.z +
                                  s_zeta_ijk.dval[2] * dcos_dxj.z);

                VECTOR fk_tri;
                fk_tri.x = pre * (-s_zeta_ijk.dval[1] * rk_hat.x +
                                  s_zeta_ijk.dval[2] * dcos_dxk.x);
                fk_tri.y = pre * (-s_zeta_ijk.dval[1] * rk_hat.y +
                                  s_zeta_ijk.dval[2] * dcos_dxk.y);
                fk_tri.z = pre * (-s_zeta_ijk.dval[1] * rk_hat.z +
                                  s_zeta_ijk.dval[2] * dcos_dxk.z);

                fi.x -= (fj_tri.x + fk_tri.x);
                fi.y -= (fj_tri.y + fk_tri.y);
                fi.z -= (fj_tri.z + fk_tri.z);
                atomicAdd(&frc[j].x, fj_tri.x);
                atomicAdd(&frc[j].y, fj_tri.y);
                atomicAdd(&frc[j].z, fj_tri.z);
                atomicAdd(&frc[k].x, fk_tri.x);
                atomicAdd(&frc[k].y, fk_tri.y);
                atomicAdd(&frc[k].z, fk_tri.z);
                if (need_virial)
                {
                    vi = vi - Get_Virial_From_Force_Dis(drij, fj_tri);
                    vi = vi - Get_Virial_From_Force_Dis(drik, fk_tri);
                }
            }
        }
        atomicAdd(&frc[i].x, fi.x);
        atomicAdd(&frc[i].y, fi.y);
        atomicAdd(&frc[i].z, fi.z);
        if (need_energy) atom_energy[i] += en_i;
        if (need_virial) atomicAdd(atom_virial + i, vi);
    }
}

void TERSOFF_INFORMATION::Initial(CONTROLLER* controller, int atom_numbers,
                                  const char* module_name,
                                  bool* need_full_nl_flag)
{
    if (module_name == NULL)
        strcpy(this->module_name, "TERSOFF");
    else
        strcpy(this->module_name, module_name);
    this->atom_numbers = atom_numbers;
    if (!controller->Command_Exist(this->module_name, "in_file")) return;
    controller->printf("START INITIALIZING TERSOFF FORCE\n");
    FILE* fp;
    Open_File_Safely(
        &fp, controller->Original_Command(this->module_name, "in_file"), "r");
    if (fscanf(fp, "%d %d", &this->atom_numbers, &this->atom_type_numbers) != 2)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, "TERSOFF_INFORMATION::Initial",
            "Reason:\n\tThe number of atoms and types can not be found\n");
    }
    std::vector<std::string> type_names(this->atom_type_numbers);
    for (int i = 0; i < this->atom_type_numbers; ++i)
    {
        char name[10];
        if (fscanf(fp, "%s", name) != 1) break;
        type_names[i] = name;
    }

    struct ParamEntry
    {
        std::string e1, e2, e3;
        float params[PARAM_STRIDE];
    };
    std::vector<ParamEntry> entries;
    while (true)
    {
        char e1[10], e2[10], e3[10];
        long current_pos = ftell(fp);
        if (fscanf(fp, "%s", e1) != 1) break;
        if (e1[0] == '#')
        {
            char tmp[1024];
            fgets(tmp, 1024, fp);
            continue;
        }
        // If we hit a number, it's probably the Atom types section
        if (isdigit(e1[0]))
        {
            fseek(fp, current_pos, SEEK_SET);
            break;
        }

        if (fscanf(fp, "%s %s", e2, e3) != 2) break;
        float p[14];
        int n_read = 0;
        for (int k = 0; k < 14; ++k) n_read += fscanf(fp, "%f", &p[k]);
        if (n_read != 14) break;
        ParamEntry entry;
        entry.e1 = e1;
        entry.e2 = e2;
        entry.e3 = e3;
        for (int i = 0; i < 14; ++i) entry.params[i] = p[i];
        entry.params[p_A] *= CONSTANT_EV_TO_KCAL_MOL;
        entry.params[p_B] *= CONSTANT_EV_TO_KCAL_MOL;
        float n = p[p_n];
        if (n > 0)
        {
            entry.params[p_c1] = powf(2.0f * n * 1.0e-16f, -1.0f / n);
            entry.params[p_c2] = powf(2.0f * n * 1.0e-8f, -1.0f / n);
            entry.params[p_c3] = 1.0f / entry.params[p_c2];
            entry.params[p_c4] = 1.0f / entry.params[p_c1];
        }
        else
            entry.params[p_c1] = entry.params[p_c2] = entry.params[p_c3] =
                entry.params[p_c4] = 0;
        entries.push_back(entry);
    }

    this->n_unique_params = entries.size();
    Malloc_Safely((void**)&h_params,
                  sizeof(float) * n_unique_params * PARAM_STRIDE);
    for (int i = 0; i < n_unique_params; ++i)
        for (int j = 0; j < PARAM_STRIDE; ++j)
            h_params[i * PARAM_STRIDE + j] = entries[i].params[j];

    int map_size = atom_type_numbers * atom_type_numbers * atom_type_numbers;
    Malloc_Safely((void**)&h_map, sizeof(int) * map_size);
    for (int i = 0; i < map_size; ++i) h_map[i] = -1;

    for (int i = 0; i < atom_type_numbers; ++i)
    {
        for (int j = 0; j < atom_type_numbers; ++j)
        {
            for (int k = 0; k < atom_type_numbers; ++k)
            {
                std::string e1 = type_names[i];
                std::string e2 = type_names[j];
                std::string e3 = type_names[k];
                for (int p = 0; p < n_unique_params; ++p)
                {
                    if (entries[p].e1 == e1 && entries[p].e2 == e2 &&
                        entries[p].e3 == e3)
                    {
                        h_map[i * atom_type_numbers * atom_type_numbers +
                              j * atom_type_numbers + k] = p;
                        break;
                    }
                }
            }
        }
    }

    Device_Malloc_And_Copy_Safely(
        (void**)&d_params, h_params,
        sizeof(float) * n_unique_params * PARAM_STRIDE);
    Device_Malloc_And_Copy_Safely((void**)&d_map, h_map,
                                  sizeof(int) * map_size);
    Device_Malloc_Safely((void**)&d_energy_sum, sizeof(float));

    Malloc_Safely((void**)&h_atom_type, sizeof(int) * atom_numbers);
    for (int i = 0; i < atom_numbers; i++)
    {
        if (fscanf(fp, "%d", h_atom_type + i) != 1)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, "TERSOFF_INFORMATION::Initial",
                "Reason:\n\tSome atom types can not be found\n");
        }
    }
    Device_Malloc_And_Copy_Safely((void**)&d_atom_type, h_atom_type,
                                  sizeof(int) * atom_numbers);

    fclose(fp);
    if (need_full_nl_flag != NULL) *need_full_nl_flag = true;
    is_initialized = true;
    controller->printf("END INITIALIZING TERSOFF FORCE\n\n");
}

void TERSOFF_INFORMATION::TERSOFF_Force_With_Atom_Energy_And_Virial(
    const int atom_numbers, const VECTOR* crd, VECTOR* frc,
    const LTMatrix3 cell, const LTMatrix3 rcell, const ATOM_GROUP* nl,
    const int need_atom_energy, float* atom_energy, const int need_virial,
    LTMatrix3* atom_virial)
{
    if (!is_initialized) return;
    if (need_atom_energy) deviceMemset(d_energy_sum, 0, sizeof(float));
    dim3 blockSize(128);
    dim3 gridSize((atom_numbers + blockSize.x - 1) / blockSize.x);

    auto force_kernel = Tersoff_Force_CUDA<false, false>;
    if (need_atom_energy && need_virial)
    {
        force_kernel = Tersoff_Force_CUDA<true, true>;
    }
    else if (need_atom_energy)
    {
        force_kernel = Tersoff_Force_CUDA<true, false>;
    }
    else if (need_virial)
    {
        force_kernel = Tersoff_Force_CUDA<false, true>;
    }

    Launch_Device_Kernel(force_kernel, gridSize, blockSize, 0, NULL,
                         atom_numbers, crd, frc, cell, rcell, nl, d_atom_type,
                         d_params, d_map, atom_type_numbers, atom_energy,
                         atom_virial, d_energy_sum);
}

void TERSOFF_INFORMATION::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized) return;
    deviceMemcpy(&h_energy_sum, d_energy_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    controller->Step_Print(this->module_name, h_energy_sum, true);
}

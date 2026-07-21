#ifdef _WIN32
#define PLUGIN_API extern "C" __declspec(dllexport)
#include <windows.h>
#ifndef RTLD_NOW
#define RTLD_NOW 0
#endif
#ifndef RTLD_GLOBAL
#define RTLD_GLOBAL 0
#endif
#define dlopen(filename, mode) LoadLibrary(filename)
#else
#define PLUGIN_API extern "C"
#include <dlfcn.h>
#endif

#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <sstream>

#include "../include/sponge_plugin_api.h"
#include "Python.h"
#include "dlpack.h"

static const SPONGE_PLUGIN_API* sponge_api = NULL;
static int is_initialized = 0;
static std::string py_script_path;
static DLDeviceType dlpack_device_type = kDLCPU;
static int dlpack_device_id = 0;

[[noreturn]] static void Python_Plugin_Contract_Error(
    const std::string& message);

static bool controller_command_exist(const char* key)
{
    return sponge_api != NULL && sponge_api->get_command != NULL &&
           sponge_api->get_command(key) != NULL;
}

static const char* controller_command(const char* key)
{
    if (!controller_command_exist(key)) return NULL;
    return sponge_api->get_command(key);
}

static void controller_printf(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    char buffer[4096];
    vsnprintf(buffer, sizeof(buffer), fmt, args);
    va_end(args);
    if (sponge_api != NULL && sponge_api->log_message != NULL)
    {
        sponge_api->log_message(buffer);
    }
}

[[noreturn]] static void PRIPS_Fatal_Error(const char* source,
                                           const std::string& message)
{
    if (sponge_api != NULL && sponge_api->report_fatal_error != NULL)
    {
        sponge_api->report_fatal_error(source, message.c_str());
    }
    fprintf(stderr, "%s: %s", source, message.c_str());
    abort();
}

PLUGIN_API std::string Name()
{
    return std::string("Python Runtime Interface Plugin");
}

#ifndef PRIPS_VERSION
#define PRIPS_VERSION "0.0.0"
#endif

PLUGIN_API std::string Version() { return std::string(PRIPS_VERSION); }

PLUGIN_API std::string Version_Check(int i)
{
    if (i != SPONGE_PRIPS_API_VERSION)
    {
        return std::string("Reason:\n\tPRIPS v" PRIPS_VERSION
                           " expects SPONGE plugin API version " +
                           std::to_string(SPONGE_PRIPS_API_VERSION) +
                           ", but got " + std::to_string(i));
    }
    return std::string();
}

static void delete_dltensor(DLManagedTensor* managed_tensor)
{
    DLTensor* tensor = &managed_tensor->dl_tensor;
    free(tensor->shape);
    if (tensor->strides != NULL) free(tensor->strides);
    free(managed_tensor);
}

static void delete_dltensor_capsule(PyObject* capsule)
{
    if (!PyCapsule_IsValid(capsule, "dltensor")) return;
    DLManagedTensor* managed_tensor = static_cast<DLManagedTensor*>(
        PyCapsule_GetPointer(capsule, "dltensor"));
    if (managed_tensor != NULL && managed_tensor->deleter != NULL)
    {
        managed_tensor->deleter(managed_tensor);
    }
}

// strides虽然dlpack支持，但是很多引擎不支持，所以尽量不给
static PyObject* create_dltensor(void* data, int N_dim, int64_t* shape,
                                 int64_t* strides, DLDataTypeCode type_code,
                                 uint8_t bits = 32, uint16_t lanes = 1)
{
    DLManagedTensor* managed_tensor =
        (DLManagedTensor*)malloc(sizeof(DLManagedTensor));
    memset(managed_tensor, 0, sizeof(DLManagedTensor));
    managed_tensor->deleter = delete_dltensor;
    DLTensor* tensor = &managed_tensor->dl_tensor;
    tensor->data = data;
    tensor->device = {dlpack_device_type, dlpack_device_id};
    tensor->ndim = N_dim;
    tensor->dtype = {(uint8_t)type_code, bits, lanes};
    tensor->shape = (int64_t*)malloc(sizeof(int64_t) * N_dim);
    memcpy(tensor->shape, shape, sizeof(int64_t) * N_dim);
    if (strides != NULL)
    {
        tensor->strides = (int64_t*)malloc(sizeof(int64_t) * N_dim);
        memcpy(tensor->strides, strides, sizeof(int64_t) * N_dim);
    }
    PyObject* a = PyCapsule_New((void*)managed_tensor, "dltensor",
                                delete_dltensor_capsule);
    return a;
}

// MD information
static PyObject* Atom_Numbers(PyObject* self, PyObject* args)
{
    return Py_BuildValue(
        "i", sponge_api == NULL || sponge_api->get_atom_numbers == NULL
                 ? 0
                 : sponge_api->get_atom_numbers());
}

static PyObject* Steps(PyObject* self, PyObject* args)
{
    return Py_BuildValue("i",
                         sponge_api == NULL || sponge_api->get_steps == NULL
                             ? 0
                             : sponge_api->get_steps());
}

static PyObject* Coordinate(PyObject* self, PyObject* args)
{
    const int atom_numbers =
        sponge_api == NULL || sponge_api->get_atom_numbers == NULL
            ? 0
            : sponge_api->get_atom_numbers();
    void* crd = sponge_api == NULL || sponge_api->get_coordinate_ptr == NULL
                    ? NULL
                    : sponge_api->get_coordinate_ptr();
    int64_t shape[2] = {atom_numbers, 3};
    return create_dltensor(crd, 2, shape, NULL, kDLFloat);
}

static PyObject* Force(PyObject* self, PyObject* args)
{
    const int atom_numbers =
        sponge_api == NULL || sponge_api->get_atom_numbers == NULL
            ? 0
            : sponge_api->get_atom_numbers();
    void* frc = sponge_api == NULL || sponge_api->get_force_ptr == NULL
                    ? NULL
                    : sponge_api->get_force_ptr();
    int64_t shape[2] = {atom_numbers, 3};
    return create_dltensor(frc, 2, shape, NULL, kDLFloat);
}

static PyObject* Force_Evaluation_Commits_Sampling_State(PyObject* self,
                                                         PyObject* args)
{
    return PyBool_FromLong(
        sponge_api != NULL &&
        sponge_api->get_force_evaluation_commits_sampling_state != NULL &&
        sponge_api->get_force_evaluation_commits_sampling_state() != 0);
}

static PyObject* Force_Evaluation_Is_Exact(PyObject* self, PyObject* args)
{
    return PyBool_FromLong(sponge_api != NULL &&
                           sponge_api->get_force_evaluation_is_exact != NULL &&
                           sponge_api->get_force_evaluation_is_exact() != 0);
}

static PyObject* Force_Evaluation_Needs_Energy(PyObject* self, PyObject* args)
{
    return PyBool_FromLong(
        sponge_api != NULL &&
        sponge_api->get_force_evaluation_needs_energy != NULL &&
        sponge_api->get_force_evaluation_needs_energy() != 0);
}

static PyObject* Force_Evaluation_Needs_Virial(PyObject* self, PyObject* args)
{
    return PyBool_FromLong(
        sponge_api != NULL &&
        sponge_api->get_force_evaluation_needs_virial != NULL &&
        sponge_api->get_force_evaluation_needs_virial() != 0);
}

// Domain decomposition
static PyObject* Local_Atom_Numbers(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_local_atom_numbers == NULL)
        return Py_BuildValue("");
    return Py_BuildValue("i", sponge_api->get_local_atom_numbers());
}

static PyObject* Local_Ghost_Numbers(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_local_ghost_numbers == NULL)
        return Py_BuildValue("");
    return Py_BuildValue("i", sponge_api->get_local_ghost_numbers());
}

static PyObject* Local_PP_Rank(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_local_pp_rank == NULL)
        return Py_BuildValue("");
    return Py_BuildValue("i", sponge_api->get_local_pp_rank());
}

static PyObject* Atom_Local(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_atom_local_ptr == NULL ||
        sponge_api->get_local_atom_numbers == NULL ||
        sponge_api->get_local_ghost_numbers == NULL)
    {
        return Py_BuildValue("");
    }
    void* atom_local = sponge_api->get_atom_local_ptr();
    if (atom_local == NULL) return Py_BuildValue("");
    int64_t shape[1] = {sponge_api->get_local_atom_numbers() +
                        sponge_api->get_local_ghost_numbers()};
    return create_dltensor(atom_local, 1, shape, NULL, kDLInt);
}

static PyObject* Atom_Local_Label(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_atom_local_label_ptr == NULL ||
        sponge_api->get_local_max_atom_numbers == NULL)
    {
        return Py_BuildValue("");
    }
    void* atom_local_label = sponge_api->get_atom_local_label_ptr();
    if (atom_local_label == NULL) return Py_BuildValue("");
    int64_t shape[1] = {sponge_api->get_local_max_atom_numbers()};
    return create_dltensor(atom_local_label, 1, shape, NULL, kDLUInt, 8);
}

static PyObject* Atom_Local_Id(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_atom_local_id_ptr == NULL ||
        sponge_api->get_local_max_atom_numbers == NULL)
    {
        return Py_BuildValue("");
    }
    void* atom_local_id = sponge_api->get_atom_local_id_ptr();
    if (atom_local_id == NULL) return Py_BuildValue("");
    int64_t shape[1] = {sponge_api->get_local_max_atom_numbers()};
    return create_dltensor(atom_local_id, 1, shape, NULL, kDLInt);
}

static PyObject* Local_Coordinate(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_local_coordinate_ptr == NULL ||
        sponge_api->get_local_atom_numbers == NULL ||
        sponge_api->get_local_ghost_numbers == NULL)
    {
        return Py_BuildValue("");
    }
    void* crd = sponge_api->get_local_coordinate_ptr();
    if (crd == NULL) return Py_BuildValue("");
    int64_t shape[2] = {sponge_api->get_local_atom_numbers() +
                            sponge_api->get_local_ghost_numbers(),
                        3};
    return create_dltensor(crd, 2, shape, NULL, kDLFloat);
}

static PyObject* Local_Force(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_local_force_ptr == NULL ||
        sponge_api->get_local_atom_numbers == NULL ||
        sponge_api->get_local_ghost_numbers == NULL)
    {
        return Py_BuildValue("");
    }
    void* frc = sponge_api->get_local_force_ptr();
    if (frc == NULL) return Py_BuildValue("");
    int64_t shape[2] = {sponge_api->get_local_atom_numbers() +
                            sponge_api->get_local_ghost_numbers(),
                        3};
    return create_dltensor(frc, 2, shape, NULL, kDLFloat);
}

static PyObject* Local_Energy(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_local_energy_ptr == NULL ||
        sponge_api->get_local_atom_numbers == NULL)
    {
        return Py_BuildValue("");
    }
    void* energy = sponge_api->get_local_energy_ptr();
    if (energy == NULL) return Py_BuildValue("");
    int64_t shape[1] = {sponge_api->get_local_atom_numbers()};
    return create_dltensor(energy, 1, shape, NULL, kDLFloat);
}

static PyObject* Local_Virial(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_local_virial_ptr == NULL ||
        sponge_api->get_local_atom_numbers == NULL)
    {
        return Py_BuildValue("");
    }
    void* virial = sponge_api->get_local_virial_ptr();
    if (virial == NULL) return Py_BuildValue("");
    int64_t shape[2] = {sponge_api->get_local_atom_numbers(), 6};
    return create_dltensor(virial, 2, shape, NULL, kDLFloat);
}

static bool Validate_Replacement_Tensor(const DLTensor* tensor,
                                        const int64_t* expected_shape,
                                        int expected_ndim,
                                        const char* buffer_name,
                                        uint64_t* byte_count)
{
    if (tensor == NULL || tensor->ndim != expected_ndim ||
        (expected_ndim > 0 && tensor->shape == NULL))
    {
        PyErr_Format(PyExc_ValueError,
                     "%s replacement must have exactly %d dimensions",
                     buffer_name, expected_ndim);
        return false;
    }
    if (tensor->dtype.code != kDLFloat || tensor->dtype.bits != 32 ||
        tensor->dtype.lanes != 1)
    {
        PyErr_Format(PyExc_TypeError, "%s replacement must have dtype float32",
                     buffer_name);
        return false;
    }
    if (tensor->device.device_type != dlpack_device_type ||
        tensor->device.device_id != dlpack_device_id)
    {
        PyErr_Format(
            PyExc_ValueError,
            "%s replacement is on DLPack device (%d, %d), expected (%d, %d)",
            buffer_name, static_cast<int>(tensor->device.device_type),
            tensor->device.device_id, static_cast<int>(dlpack_device_type),
            dlpack_device_id);
        return false;
    }

    uint64_t elements = 1;
    int64_t contiguous_stride = 1;
    for (int axis = expected_ndim - 1; axis >= 0; --axis)
    {
        if (tensor->shape[axis] != expected_shape[axis])
        {
            PyErr_Format(PyExc_ValueError,
                         "%s replacement shape mismatch at axis %d: got %lld, "
                         "expected %lld",
                         buffer_name, axis,
                         static_cast<long long>(tensor->shape[axis]),
                         static_cast<long long>(expected_shape[axis]));
            return false;
        }
        if (tensor->strides != NULL &&
            tensor->strides[axis] != contiguous_stride)
        {
            PyErr_Format(PyExc_ValueError,
                         "%s replacement must be C-contiguous", buffer_name);
            return false;
        }
        const uint64_t extent = static_cast<uint64_t>(expected_shape[axis]);
        if (extent != 0 &&
            elements > std::numeric_limits<uint64_t>::max() / extent)
        {
            PyErr_Format(PyExc_OverflowError,
                         "%s replacement size overflows uint64_t", buffer_name);
            return false;
        }
        elements *= extent;
        if (expected_shape[axis] != 0 &&
            contiguous_stride >
                std::numeric_limits<int64_t>::max() / expected_shape[axis])
        {
            PyErr_Format(PyExc_OverflowError,
                         "%s replacement strides overflow int64_t",
                         buffer_name);
            return false;
        }
        contiguous_stride *= expected_shape[axis];
    }
    if (elements > std::numeric_limits<uint64_t>::max() / sizeof(float))
    {
        PyErr_Format(PyExc_OverflowError,
                     "%s replacement byte size overflows uint64_t",
                     buffer_name);
        return false;
    }
    *byte_count = elements * sizeof(float);
    if (*byte_count != 0 && tensor->data == NULL)
    {
        PyErr_Format(PyExc_ValueError, "%s replacement has a null data pointer",
                     buffer_name);
        return false;
    }
    const uintptr_t base = reinterpret_cast<uintptr_t>(tensor->data);
    if (tensor->byte_offset >
        static_cast<uint64_t>(std::numeric_limits<uintptr_t>::max() - base))
    {
        PyErr_Format(PyExc_OverflowError,
                     "%s replacement byte offset overflows the address "
                     "space",
                     buffer_name);
        return false;
    }
    const uintptr_t source_address =
        base + static_cast<uintptr_t>(tensor->byte_offset);
    if (*byte_count >
        static_cast<uint64_t>(std::numeric_limits<uintptr_t>::max() -
                              source_address))
    {
        PyErr_Format(PyExc_OverflowError,
                     "%s replacement address range overflows the address "
                     "space",
                     buffer_name);
        return false;
    }
    return true;
}

static PyObject* Replace_Device_Buffer(PyObject* source, void* destination,
                                       const int64_t* expected_shape,
                                       int expected_ndim,
                                       const char* buffer_name)
{
    if (sponge_api == NULL || sponge_api->copy_device_buffer == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "SPONGE does not provide synchronous device-buffer "
                        "writeback support");
        return NULL;
    }
    PyObject* capsule = PyObject_CallMethod(source, "__dlpack__", NULL);
    if (capsule == NULL) return NULL;

    const bool is_versioned =
        PyCapsule_IsValid(capsule, "dltensor_versioned") != 0;
    const bool is_legacy = PyCapsule_IsValid(capsule, "dltensor") != 0;
    if (!is_versioned && !is_legacy)
    {
        Py_DECREF(capsule);
        PyErr_Format(PyExc_TypeError,
                     "%s replacement __dlpack__() did not return a valid "
                     "DLPack capsule",
                     buffer_name);
        return NULL;
    }

    DLManagedTensorVersioned* versioned = NULL;
    DLManagedTensor* legacy = NULL;
    DLTensor* tensor = NULL;
    if (is_versioned)
    {
        versioned = static_cast<DLManagedTensorVersioned*>(
            PyCapsule_GetPointer(capsule, "dltensor_versioned"));
        if (versioned == NULL)
        {
            Py_DECREF(capsule);
            return NULL;
        }
        if (versioned->version.major != DLPACK_MAJOR_VERSION)
        {
            const unsigned int producer_major = versioned->version.major;
            Py_DECREF(capsule);
            PyErr_Format(PyExc_RuntimeError,
                         "%s replacement uses unsupported DLPack major "
                         "version %u (expected %u)",
                         buffer_name, producer_major, DLPACK_MAJOR_VERSION);
            return NULL;
        }
        tensor = &versioned->dl_tensor;
    }
    else
    {
        legacy = static_cast<DLManagedTensor*>(
            PyCapsule_GetPointer(capsule, "dltensor"));
        if (legacy == NULL)
        {
            Py_DECREF(capsule);
            return NULL;
        }
        tensor = &legacy->dl_tensor;
    }

    uint64_t byte_count = 0;
    if (!Validate_Replacement_Tensor(tensor, expected_shape, expected_ndim,
                                     buffer_name, &byte_count))
    {
        Py_DECREF(capsule);
        return NULL;
    }
    const uintptr_t source_address =
        reinterpret_cast<uintptr_t>(tensor->data) + tensor->byte_offset;

    const char* used_capsule_name =
        is_versioned ? "used_dltensor_versioned" : "used_dltensor";
    if (PyCapsule_SetName(capsule, used_capsule_name) != 0)
    {
        Py_DECREF(capsule);
        return NULL;
    }
    const int copy_status = sponge_api->copy_device_buffer(
        destination, reinterpret_cast<const void*>(source_address), byte_count);
    if (is_versioned)
    {
        if (versioned->deleter != NULL) versioned->deleter(versioned);
    }
    else if (legacy->deleter != NULL)
    {
        legacy->deleter(legacy);
    }
    Py_DECREF(capsule);
    if (copy_status != 0)
    {
        PyErr_Format(PyExc_RuntimeError,
                     "synchronous writeback to %s failed with device status "
                     "%d",
                     buffer_name, copy_status);
        return NULL;
    }
    Py_RETURN_NONE;
}

static bool Local_Buffer_Extents(int* local_atoms, int* local_with_ghosts)
{
    if (sponge_api == NULL || sponge_api->get_local_atom_numbers == NULL ||
        sponge_api->get_local_ghost_numbers == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "SPONGE local-domain metadata is unavailable");
        return false;
    }
    const int owned = sponge_api->get_local_atom_numbers();
    const int ghosts = sponge_api->get_local_ghost_numbers();
    if (owned < 0 || ghosts < 0 ||
        owned > std::numeric_limits<int>::max() - ghosts)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "SPONGE reported invalid local-domain extents");
        return false;
    }
    *local_atoms = owned;
    *local_with_ghosts = owned + ghosts;
    return true;
}

static PyObject* Replace_Local_Force(PyObject* self, PyObject* args)
{
    PyObject* source = NULL;
    if (!PyArg_ParseTuple(args, "O", &source)) return NULL;
    int local_atoms = 0;
    int local_with_ghosts = 0;
    if (!Local_Buffer_Extents(&local_atoms, &local_with_ghosts)) return NULL;
    if (sponge_api->get_local_force_ptr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "SPONGE local force buffer is unavailable");
        return NULL;
    }
    int64_t shape[2] = {local_with_ghosts, 3};
    return Replace_Device_Buffer(source, sponge_api->get_local_force_ptr(),
                                 shape, 2, "local force");
}

static PyObject* Replace_Local_Energy(PyObject* self, PyObject* args)
{
    PyObject* source = NULL;
    if (!PyArg_ParseTuple(args, "O", &source)) return NULL;
    int local_atoms = 0;
    int local_with_ghosts = 0;
    if (!Local_Buffer_Extents(&local_atoms, &local_with_ghosts)) return NULL;
    if (sponge_api->get_local_energy_ptr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "SPONGE local energy buffer is unavailable");
        return NULL;
    }
    int64_t shape[1] = {local_atoms};
    return Replace_Device_Buffer(source, sponge_api->get_local_energy_ptr(),
                                 shape, 1, "local energy");
}

static PyObject* Replace_Local_Virial(PyObject* self, PyObject* args)
{
    PyObject* source = NULL;
    if (!PyArg_ParseTuple(args, "O", &source)) return NULL;
    int local_atoms = 0;
    int local_with_ghosts = 0;
    if (!Local_Buffer_Extents(&local_atoms, &local_with_ghosts)) return NULL;
    if (sponge_api->get_local_virial_ptr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "SPONGE local virial buffer is unavailable");
        return NULL;
    }
    int64_t shape[2] = {local_atoms, 6};
    return Replace_Device_Buffer(source, sponge_api->get_local_virial_ptr(),
                                 shape, 2, "local virial");
}

// Neighbor List
static PyObject* Neighbor_List_Index(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_neighbor_list_index_ptr == NULL ||
        sponge_api->get_atom_numbers == NULL ||
        sponge_api->get_neighbor_list_max_numbers == NULL)
    {
        return Py_BuildValue("");
    }
    void* nl_index = sponge_api->get_neighbor_list_index_ptr();
    if (nl_index == NULL) return Py_BuildValue("");
    int64_t shape[2] = {sponge_api->get_atom_numbers(),
                        sponge_api->get_neighbor_list_max_numbers()};
    return create_dltensor(nl_index, 2, shape, NULL, kDLInt);
}

static PyObject* Neighbor_List_Numbers(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_neighbor_list_count == NULL ||
        sponge_api->get_atom_numbers == NULL)
    {
        return Py_BuildValue("");
    }
    const int atom_numbers = sponge_api->get_atom_numbers();
    PyObject* numbers = PyList_New(atom_numbers);
    for (int i = 0; i < atom_numbers; i++)
    {
        PyList_SET_ITEM(
            numbers, i,
            PyLong_FromLong(sponge_api->get_neighbor_list_count(i)));
    }
    return numbers;
}

static PyObject* Neighbor_List_Max_Numbers(PyObject* self, PyObject* args)
{
    if (sponge_api == NULL || sponge_api->get_neighbor_list_max_numbers == NULL)
    {
        return Py_BuildValue("");
    }
    return Py_BuildValue("i", sponge_api->get_neighbor_list_max_numbers());
}

// Domain information
// CONTROLLER
static PyObject* Control_Printf(PyObject* self, PyObject* args, PyObject* kw)
{
    static char* kwlist[] = {(char*)"toprint", NULL};
    char* buffer;
    if (!PyArg_ParseTupleAndKeywords(args, kw, "s", kwlist, &buffer))
    {
        return NULL;
    }
    controller_printf("%s", buffer);
    return Py_BuildValue("");
}

static PyObject* Control_MPI_Rank(PyObject* self, PyObject* args)
{
    return Py_BuildValue("i",
                         sponge_api == NULL || sponge_api->get_mpi_rank == NULL
                             ? 0
                             : sponge_api->get_mpi_rank());
}

static PyObject* Fatal_Error(PyObject* self, PyObject* args)
{
    const char* source = NULL;
    const char* message = NULL;
    if (!PyArg_ParseTuple(args, "ss", &source, &message)) return NULL;
    if (sponge_api != NULL && sponge_api->report_fatal_error != NULL)
    {
        sponge_api->report_fatal_error(source, message);
    }
    fprintf(stderr, "%s: %s\n", source, message);
    abort();
}

static PyMethodDef SpongeMethods[] = {
    {"_atom_numbers", (PyCFunction)Atom_Numbers, METH_VARARGS, ""},
    {"_steps", (PyCFunction)Steps, METH_VARARGS, ""},
    {"_coordinate", (PyCFunction)Coordinate, METH_VARARGS, ""},
    {"_force", (PyCFunction)Force, METH_VARARGS, ""},
    {"_force_evaluation_commits_sampling_state",
     (PyCFunction)Force_Evaluation_Commits_Sampling_State, METH_VARARGS, ""},
    {"_force_evaluation_is_exact", (PyCFunction)Force_Evaluation_Is_Exact,
     METH_VARARGS, ""},
    {"_force_evaluation_needs_energy",
     (PyCFunction)Force_Evaluation_Needs_Energy, METH_VARARGS, ""},
    {"_force_evaluation_needs_virial",
     (PyCFunction)Force_Evaluation_Needs_Virial, METH_VARARGS, ""},
    {"_local_atom_numbers", (PyCFunction)Local_Atom_Numbers, METH_VARARGS, ""},
    {"_local_ghost_numbers", (PyCFunction)Local_Ghost_Numbers, METH_VARARGS,
     ""},
    {"_local_pp_rank", (PyCFunction)Local_PP_Rank, METH_VARARGS, ""},
    {"_atom_local", (PyCFunction)Atom_Local, METH_VARARGS, ""},
    {"_atom_local_label", (PyCFunction)Atom_Local_Label, METH_VARARGS, ""},
    {"_atom_local_id", (PyCFunction)Atom_Local_Id, METH_VARARGS, ""},
    {"_local_crd", (PyCFunction)Local_Coordinate, METH_VARARGS, ""},
    {"_local_frc", (PyCFunction)Local_Force, METH_VARARGS, ""},
    {"_local_energy", (PyCFunction)Local_Energy, METH_VARARGS, ""},
    {"_local_virial", (PyCFunction)Local_Virial, METH_VARARGS, ""},
    {"_replace_local_force", (PyCFunction)Replace_Local_Force, METH_VARARGS,
     ""},
    {"_replace_local_energy", (PyCFunction)Replace_Local_Energy, METH_VARARGS,
     ""},
    {"_replace_local_virial", (PyCFunction)Replace_Local_Virial, METH_VARARGS,
     ""},
    {"_neighbor_list_number", (PyCFunction)Neighbor_List_Numbers, METH_VARARGS,
     ""},
    {"_neighbor_list_index", (PyCFunction)Neighbor_List_Index, METH_VARARGS,
     ""},
    {"_neighbor_list_max_numbers", (PyCFunction)Neighbor_List_Max_Numbers,
     METH_VARARGS, ""},
    {"_printf", (PyCFunction)Control_Printf, METH_VARARGS | METH_KEYWORDS, ""},
    {"_MPI_rank", (PyCFunction)Control_MPI_Rank, METH_VARARGS, ""},
    {"_fatal_error", (PyCFunction)Fatal_Error, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}};

static PyModuleDef SpongeModule = {PyModuleDef_HEAD_INIT,
                                   "Sponge",
                                   NULL,
                                   -1,
                                   SpongeMethods,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL};

PyMODINIT_FUNC PyInit_sponge(void) { return PyModule_Create(&SpongeModule); }

static PyMethodDef PripsMethods[] = {{NULL, NULL, 0, NULL}};

static PyModuleDef pripsModule = {PyModuleDef_HEAD_INIT,
                                  "prips",
                                  NULL,
                                  -1,
                                  PripsMethods,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL};

PyMODINIT_FUNC PyInit_prips(void)
{
    PyObject* m = PyModule_Create(&pripsModule);
    PyModule_AddStringConstant(m, "__version__", Version().c_str());
    PyObject* m0 = PyInit_sponge();
    PyModule_AddIntConstant(m0, "_backend",
                            static_cast<int>(dlpack_device_type));
    PyModule_AddIntConstant(m0, "_device_id", dlpack_device_id);
    PyModule_AddIntConstant(m0, "FORCE_ENERGY_COMPLETE",
                            SPONGE_PLUGIN_FORCE_ENERGY_COMPLETE);
    PyModule_AddIntConstant(m0, "FORCE_VIRIAL_COMPLETE",
                            SPONGE_PLUGIN_FORCE_VIRIAL_COMPLETE);
    PyModule_AddIntConstant(m0, "FORCE_PURE", SPONGE_PLUGIN_FORCE_PURE);
    PyModule_AddIntConstant(m0, "FORCE_TRANSACTIONAL",
                            SPONGE_PLUGIN_FORCE_TRANSACTIONAL);
    PyModule_AddObject(m, "Sponge", m0);
    return m;
}

PLUGIN_API void Set_Backend_Device_Type(int device_type)
{
    dlpack_device_type = static_cast<DLDeviceType>(device_type);
}

PLUGIN_API void Initial_Stable(const SPONGE_PLUGIN_API* api)
{
    // Publish no Python module or capabilities until every initialization
    // step, including the required force entry point, has succeeded.
    is_initialized = 0;
    py_script_path.clear();
    sponge_api = api;
    if (sponge_api == NULL)
    {
        fprintf(stderr,
                "        PRIPS initialization failed: null stable API.\n");
        abort();
    }
    if (sponge_api->api_version != SPONGE_PRIPS_API_VERSION ||
        sponge_api->report_fatal_error == NULL ||
        sponge_api->copy_device_buffer == NULL ||
        sponge_api->get_device_id == NULL)
    {
        const std::string reason =
            "Reason:\n\tPRIPS requires stable plugin API v" +
            std::to_string(SPONGE_PRIPS_API_VERSION) +
            " with synchronous device-buffer writeback\n";
        // report_fatal_error was added before v5; negotiated loaders always
        // provide it.  Avoid reading beyond the v2/v3 prefix if an obsolete
        // core bypasses Version_Check and calls this entry point directly.
        if (sponge_api->api_version >= 4 &&
            sponge_api->report_fatal_error != NULL)
        {
            sponge_api->report_fatal_error("PRIPS::Initial_Stable",
                                           reason.c_str());
        }
        fprintf(stderr, "%s", reason.c_str());
        abort();
    }
    dlpack_device_type = static_cast<DLDeviceType>(sponge_api->device_type);
    dlpack_device_id = sponge_api->get_device_id();
    if (dlpack_device_id < 0)
    {
        sponge_api->report_fatal_error(
            "PRIPS::Initial_Stable",
            "Reason:\n\tSPONGE reported a negative active device id\n");
        abort();
    }
    controller_printf("    initializing pyplugin\n");
    if (controller_command_exist("py"))
    {
        py_script_path = controller_command("py");
    }
    else
    {
        const char* py_env = std::getenv("SPONGE_PRIPS_PY");
        if (py_env != NULL && py_env[0] != '\0')
        {
            py_script_path = py_env;
            controller_printf(
                "        No 'py' command found. Falling back to "
                "SPONGE_PRIPS_PY.\n");
        }
        else
        {
            PRIPS_Fatal_Error(
                "PRIPS::Initial_Stable",
                "Reason:\n\tPRIPS requires a Python force module; configure "
                "the 'py' command or SPONGE_PRIPS_PY\n");
        }
    }
    if (PyImport_AppendInittab("prips", &PyInit_prips) != 0)
    {
        PRIPS_Fatal_Error(
            "PRIPS::Initial_Stable",
            "Reason:\n\tfailed to register the embedded prips module\n");
    }
    Py_Initialize();
    if (!Py_IsInitialized())
    {
        if (sponge_api->report_fatal_error != NULL)
        {
            sponge_api->report_fatal_error(
                "PRIPS::Initial_Stable",
                "Reason:\n\tPython initialization failed\n");
        }
        abort();
    }
    else
    {
        controller_printf("        Python Initialized\n");
    }

    wchar_t* temp_args[1] = {(wchar_t*)L"SPONGE"};
    PySys_SetArgv(1, temp_args);
    if (PyRun_SimpleString(R"XYJ(
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import importlib.util as iu
import traceback as traceback_module

from prips import Sponge

def new_hook(exctype, value, traceback):
    detail = ''.join(
        traceback_module.format_exception(exctype, value, traceback)
    )
    Sponge._fatal_error("PRIPS Python", detail)
sys.excepthook = new_hook

class SpongeDLPackTensor:
    def __init__(self, capsule):
        self.capsule = capsule

    def __dlpack__(self, *args, **kwargs):
        return self.capsule

    def __dlpack_device__(self):
        return (Sponge._backend, Sponge._device_id)

Sponge.SpongeDLPackTensor = SpongeDLPackTensor

class MD_INFORMATION:
    def __init__(self):
        self._crd = Sponge.SpongeDLPackTensor(Sponge._coordinate())
        self._frc = Sponge.SpongeDLPackTensor(Sponge._force())

    class system_information:
        @property
        def steps(self):
            return Sponge._steps()

    @property
    def atom_numbers(self):
        return Sponge._atom_numbers()

    @property
    def crd(self):
        return self._crd

    @property
    def frc(self):
        return self._frc

Sponge.MD_INFORMATION = MD_INFORMATION
Sponge.md_info = MD_INFORMATION()
Sponge.md_info.sys = MD_INFORMATION.system_information()

class FORCE_EVALUATION_CONTEXT:
    @property
    def commits_sampling_state(self):
        return Sponge._force_evaluation_commits_sampling_state()

    @property
    def is_exact(self):
        return Sponge._force_evaluation_is_exact()

    @property
    def needs_energy(self):
        return Sponge._force_evaluation_needs_energy()

    @property
    def needs_virial(self):
        return Sponge._force_evaluation_needs_virial()

Sponge.FORCE_EVALUATION_CONTEXT = FORCE_EVALUATION_CONTEXT
Sponge.force_evaluation = FORCE_EVALUATION_CONTEXT()

class FORCE_RESULT:
    def __init__(self, force, *, energy=None, virial=None):
        if force is None:
            raise ValueError("force_result requires a complete updated force buffer")
        self.force = force
        self.energy = energy
        self.virial = virial

def force_result(force, *, energy=None, virial=None):
    """Return complete updated local buffers from a functional force hook."""
    return FORCE_RESULT(force, energy=energy, virial=virial)

Sponge.FORCE_RESULT = FORCE_RESULT
Sponge.force_result = force_result

class NEIGHBOR_LIST:
    def __init__(self):
        index_capsule = Sponge._neighbor_list_index()
        self._index = None if index_capsule is None else Sponge.backend(Sponge.SpongeDLPackTensor(index_capsule))

    @property
    def index(self):
        return self._index

    @property
    def number(self):
        return Sponge._neighbor_list_number()

    @property
    def max_neighbor_numbers(self):
        return Sponge._neighbor_list_max_numbers()

Sponge.NEIGHBOR_LIST = NEIGHBOR_LIST
Sponge.neighbor_list = None

class DOMAIN_INFORMATION:
    @property
    def atom_numbers(self):
        return Sponge._local_atom_numbers()

    @property
    def ghost_numbers(self):
        return Sponge._local_ghost_numbers()

    @property
    def pp_rank(self):
        return Sponge._local_pp_rank()

    def __init__(self):
        atom_local_capsule = Sponge._atom_local()
        atom_local_label_capsule = Sponge._atom_local_label()
        atom_local_id_capsule = Sponge._atom_local_id()
        crd_capsule = Sponge._local_crd()
        frc_capsule = Sponge._local_frc()
        energy_capsule = Sponge._local_energy()
        virial_capsule = Sponge._local_virial()
        self._atom_local = None if atom_local_capsule is None else Sponge.backend(Sponge.SpongeDLPackTensor(atom_local_capsule))
        self._atom_local_label = None if atom_local_label_capsule is None else Sponge.backend(Sponge.SpongeDLPackTensor(atom_local_label_capsule))
        self._atom_local_id = None if atom_local_id_capsule is None else Sponge.backend(Sponge.SpongeDLPackTensor(atom_local_id_capsule))
        self._crd = None if crd_capsule is None else Sponge.backend(Sponge.SpongeDLPackTensor(crd_capsule))
        self._frc = None if frc_capsule is None else Sponge.backend(Sponge.SpongeDLPackTensor(frc_capsule))
        self._energy = None if energy_capsule is None else Sponge.backend(Sponge.SpongeDLPackTensor(energy_capsule))
        self._virial = None if virial_capsule is None else Sponge.backend(Sponge.SpongeDLPackTensor(virial_capsule))

    @property
    def atom_local(self):
        return self._atom_local

    @property
    def atom_local_label(self):
        return self._atom_local_label

    @property
    def atom_local_id(self):
        return self._atom_local_id

    @property
    def crd(self):
        return self._crd

    @property
    def frc(self):
        return self._frc

    @property
    def energy(self):
        return self._energy

    @property
    def virial(self):
        return self._virial

Sponge.DOMAIN_INFORMATION = DOMAIN_INFORMATION
Sponge.dd = None

class CONTROLLER:
    @property
    def MPI_rank(self):
        return Sponge._MPI_rank()

    def printf(self, *values, sep=" ", end="\n"):
        return Sponge._printf(sep.join([f"{i}" for i in values]) + end)

Sponge.CONTROLLER = CONTROLLER
Sponge.controller = CONTROLLER()

def _numpy_backend(dlpack_tensor):
    import ctypes
    import numpy as np

    class DLDevice(ctypes.Structure):
        _fields_ = [('device_type', ctypes.c_int), ('device_id', ctypes.c_int)]

    class DLDataType(ctypes.Structure):
        _fields_ = [('code', ctypes.c_uint8), ('bits', ctypes.c_uint8), ('lanes', ctypes.c_uint16)]

    class DLTensor(ctypes.Structure):
        _fields_ = [
            ('data', ctypes.c_void_p),
            ('device', DLDevice),
            ('ndim', ctypes.c_int),
            ('dtype', DLDataType),
            ('shape', ctypes.POINTER(ctypes.c_int64)),
            ('strides', ctypes.POINTER(ctypes.c_int64)),
            ('byte_offset', ctypes.c_uint64),
        ]

    class DLManagedTensor(ctypes.Structure):
        _fields_ = [('dl_tensor', DLTensor), ('manager_ctx', ctypes.c_void_p), ('deleter', ctypes.c_void_p)]

    pycapsule_get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    pycapsule_get_pointer.restype = ctypes.c_void_p
    pycapsule_get_pointer.argtypes = [ctypes.py_object, ctypes.c_char_p]

    capsule = dlpack_tensor.capsule
    ptr = pycapsule_get_pointer(capsule, b'dltensor')
    managed = ctypes.cast(ptr, ctypes.POINTER(DLManagedTensor)).contents
    tensor = managed.dl_tensor
    shape = tuple(tensor.shape[i] for i in range(tensor.ndim))
    if tensor.dtype.lanes != 1:
        raise TypeError(f'unsupported dlpack lanes: {tensor.dtype.lanes}')
    dtype_map = {
        (0, 32): (ctypes.c_int32, np.int32),
        (1, 8): (ctypes.c_uint8, np.uint8),
        (2, 32): (ctypes.c_float, np.float32),
    }
    key = (int(tensor.dtype.code), int(tensor.dtype.bits))
    if key not in dtype_map:
        raise TypeError(
            f'unsupported dlpack dtype: code={tensor.dtype.code}, '
            f'bits={tensor.dtype.bits}, lanes={tensor.dtype.lanes}'
        )
    ctype, np_dtype = dtype_map[key]
    size = int(np.prod(shape, dtype=np.int64))
    data_ptr = tensor.data + tensor.byte_offset
    buffer = (ctype * size).from_address(data_ptr)
    return np.ctypeslib.as_array(buffer).view(np_dtype).reshape(shape)

def _resolve_backend(name):
    backend = name.lower()
    if backend == "numpy":
        return _numpy_backend
    if backend == "jax":
        import jax.dlpack
        return jax.dlpack.from_dlpack
    if backend == "cupy":
        import cupy
        return cupy.from_dlpack
    if backend in ("pytorch", "torch"):
        import torch
        return torch.from_dlpack
    raise ValueError(
        "backend must be one of: 'numpy', 'jax', 'cupy', 'pytorch'"
    )

def _device_name(device_type):
    names = {
        1: "cpu",
        2: "cuda",
        10: "rocm",
    }
    return names.get(device_type, f"device_type={device_type}")

def set_backend(backend):
    if not isinstance(backend, str):
        raise TypeError("backend must be a string")
    if Sponge.backend_name is not None:
        raise RuntimeError("backend has already been set")
    backend_name = backend.lower()
    device_type = Sponge._backend
    if backend_name == "numpy" and device_type != 1:
        raise ValueError(
            f"backend 'numpy' requires CPU tensors, but current device is "
            f"{_device_name(device_type)}"
        )
    if backend_name == "cupy" and device_type == 1:
        raise ValueError(
            "backend 'cupy' requires GPU tensors, but current device is "
            f"{_device_name(device_type)}"
        )
    Sponge.backend = _resolve_backend(backend_name)
    Sponge.md_info._crd = Sponge.backend(Sponge.md_info._crd)
    Sponge.md_info._frc = Sponge.backend(Sponge.md_info._frc)
    Sponge.backend_name = backend_name

Sponge.set_backend = set_backend
Sponge.backend_name = None

def _refresh_backend_views():
    if Sponge.backend_name != "jax":
        return
    Sponge.md_info._crd = Sponge.backend(
        Sponge.SpongeDLPackTensor(Sponge._coordinate())
    )
    Sponge.md_info._frc = Sponge.backend(
        Sponge.SpongeDLPackTensor(Sponge._force())
    )
    if Sponge.neighbor_list is not None:
        Sponge.neighbor_list = Sponge.NEIGHBOR_LIST()
    if Sponge.dd is not None:
        Sponge.dd = Sponge.DOMAIN_INFORMATION()

def _block_until_ready(value):
    block = getattr(value, "block_until_ready", None)
    if callable(block):
        ready_value = block()
        return value if ready_value is None else ready_value
    return value

def _apply_force_result(result, *, require_result):
    if result is None:
        if require_result:
            raise RuntimeError(
                "a JAX Calculate_Force hook must return Sponge.force_result("
                "force, energy=..., virial=...) with complete updated local "
                "buffers; JAX arrays cannot update SPONGE in place"
            )
        return
    if not isinstance(result, Sponge.FORCE_RESULT):
        raise TypeError(
            "Calculate_Force return value must be Sponge.force_result(...)"
        )
    if Sponge.dd is None:
        raise RuntimeError(
            "functional force writeback requires local domain buffers"
        )

    capabilities = getattr(
        sponge_pyplugin, "SPONGE_FORCE_CAPABILITIES", 0
    )
    if (
        Sponge.force_evaluation.needs_energy
        and capabilities & Sponge.FORCE_ENERGY_COMPLETE
        and result.energy is None
    ):
        raise RuntimeError(
            "ENERGY_COMPLETE Calculate_Force must return the complete updated "
            "local energy buffer when energy is requested"
        )
    if (
        Sponge.force_evaluation.needs_virial
        and capabilities & Sponge.FORCE_VIRIAL_COMPLETE
        and result.virial is None
    ):
        raise RuntimeError(
            "VIRIAL_COMPLETE Calculate_Force must return the complete updated "
            "local virial buffer when virial is requested"
        )

    # The result object and each producer-owned array remain strongly
    # referenced through the synchronous DLPack copy.  block_until_ready()
    # closes JAX's asynchronous producer side before the core consumes the
    # device pointer.
    Sponge._replace_local_force(_block_until_ready(result.force))
    if result.energy is not None:
        Sponge._replace_local_energy(_block_until_ready(result.energy))
    if result.virial is not None:
        Sponge._replace_local_virial(_block_until_ready(result.virial))

Sponge._refresh_backend_views = _refresh_backend_views
Sponge._apply_force_result = _apply_force_result

)XYJ") != 0)
    {
        Python_Plugin_Contract_Error(
            "failed to initialize the embedded Python API");
    }

    PyObject* main_module = PyImport_AddModule("__main__");
    PyObject* main_globals =
        main_module == NULL ? NULL : PyModule_GetDict(main_module);
    PyObject* sponge_module =
        main_globals == NULL ? NULL
                             : PyDict_GetItemString(main_globals, "Sponge");
    PyObject* python_script_path =
        PyUnicode_DecodeFSDefault(py_script_path.c_str());
    if (sponge_module == NULL || python_script_path == NULL ||
        PyObject_SetAttrString(sponge_module, "fname", python_script_path) != 0)
    {
        Py_XDECREF(python_script_path);
        Python_Plugin_Contract_Error("failed to expose the Python plugin path");
    }
    Py_DECREF(python_script_path);

    if (PyRun_SimpleString(R"XYJ(sponge_pyplugin = None
sponge_pyplugin_path = Path(Sponge.fname)
spec = iu.spec_from_file_location('sponge_pyplugin', sponge_pyplugin_path)
if spec is None or spec.loader is None:
    raise ImportError(
        f"cannot create an import specification for PRIPS module "
        f"'{sponge_pyplugin_path}'"
    )
candidate_module = iu.module_from_spec(spec)
spec.loader.exec_module(candidate_module)
if not callable(getattr(candidate_module, 'Calculate_Force', None)):
    raise TypeError(
        f"PRIPS module '{sponge_pyplugin_path}' must define callable "
        "Calculate_Force"
    )
sponge_pyplugin = candidate_module
Sponge.controller.printf(
    "        module '%s' imported." % sponge_pyplugin_path.stem
)
)XYJ") != 0)
    {
        Python_Plugin_Contract_Error("failed to load the Python plugin");
    }

    is_initialized = 1;
    controller_printf("    end initializing pyplugin\n");
}

PLUGIN_API void Initial(void* md, void* ctrl, void* nl, void* cv, void* cv_map,
                        void* cv_instance_map)
{
    (void)md;
    (void)ctrl;
    (void)nl;
    (void)cv;
    (void)cv_map;
    (void)cv_instance_map;
    fprintf(stderr,
            "PRIPS requires SPONGE stable plugin API support. "
            "Rebuild SPONGE with Initial_Stable loader support.\n");
    exit(1);
}

static PyObject* Python_Plugin_Module()
{
    if (!is_initialized || !Py_IsInitialized()) return NULL;
    PyObject* main_module = PyImport_AddModule("__main__");
    if (main_module == NULL) return NULL;
    PyObject* main_globals = PyModule_GetDict(main_module);
    if (main_globals == NULL) return NULL;
    PyObject* plugin = PyDict_GetItemString(main_globals, "sponge_pyplugin");
    return plugin == NULL || plugin == Py_None || !PyModule_Check(plugin)
               ? NULL
               : plugin;
}

[[noreturn]] static void Python_Plugin_Contract_Error(
    const std::string& message)
{
    std::string python_detail;
    if (PyErr_Occurred())
    {
        // PyErr_Print invokes sys.excepthook, which enters the core fatal path
        // recursively.  Extract the exception text directly so this function
        // reports exactly one controller error.
        PyObject* exception_type = NULL;
        PyObject* exception_value = NULL;
        PyObject* exception_traceback = NULL;
        PyErr_Fetch(&exception_type, &exception_value, &exception_traceback);
        PyErr_NormalizeException(&exception_type, &exception_value,
                                 &exception_traceback);
        if (exception_value != NULL)
        {
            PyObject* text = PyObject_Str(exception_value);
            if (text != NULL)
            {
                const char* utf8 = PyUnicode_AsUTF8(text);
                if (utf8 != NULL) python_detail = std::string(": ") + utf8;
                Py_DECREF(text);
            }
        }
        Py_XDECREF(exception_type);
        Py_XDECREF(exception_value);
        Py_XDECREF(exception_traceback);
        PyErr_Clear();
    }
    const std::string reason =
        "Reason:\n\tPRIPS force contract error: " + message + python_detail +
        "\n";
    if (sponge_api != NULL && sponge_api->report_fatal_error != NULL)
    {
        sponge_api->report_fatal_error("PRIPS", reason.c_str());
    }
    // This fallback is reachable only when PRIPS was not given a valid core
    // API; normal plugin contract failures terminate through the controller.
    fprintf(stderr, "%s", reason.c_str());
    abort();
}

static bool Python_Plugin_Has_Callable(PyObject* plugin, const char* name)
{
    PyObject* attribute = PyObject_GetAttrString(plugin, name);
    if (attribute == NULL)
    {
        PyErr_Clear();
        return false;
    }
    const bool is_callable = PyCallable_Check(attribute) != 0;
    Py_DECREF(attribute);
    return is_callable;
}

PLUGIN_API uint32_t Get_Force_Capabilities()
{
    PyObject* plugin = Python_Plugin_Module();
    if (plugin == NULL)
    {
        if (!is_initialized) return 0;
        Python_Plugin_Contract_Error(
            "initialized Python force module is unavailable");
    }

    PyObject* value =
        PyObject_GetAttrString(plugin, "SPONGE_FORCE_CAPABILITIES");
    if (value == NULL)
    {
        if (PyErr_ExceptionMatches(PyExc_AttributeError))
        {
            PyErr_Clear();
            return 0;
        }
        Python_Plugin_Contract_Error(
            "failed to read SPONGE_FORCE_CAPABILITIES");
    }
    if (PyBool_Check(value) || !PyLong_Check(value))
    {
        Py_DECREF(value);
        Python_Plugin_Contract_Error(
            "SPONGE_FORCE_CAPABILITIES must be an integer bit mask");
    }
    const unsigned long long raw_capabilities =
        PyLong_AsUnsignedLongLong(value);
    Py_DECREF(value);
    if (PyErr_Occurred() ||
        raw_capabilities > std::numeric_limits<uint32_t>::max())
    {
        Python_Plugin_Contract_Error(
            "SPONGE_FORCE_CAPABILITIES must fit in uint32_t");
    }
    const uint32_t capabilities = static_cast<uint32_t>(raw_capabilities);

    if ((capabilities & SPONGE_PLUGIN_FORCE_TRANSACTIONAL) != 0)
    {
        const char* required_hooks[] = {"Begin_Force_Transaction",
                                        "Commit_Force_Transaction",
                                        "Rollback_Force_Transaction"};
        for (const char* hook : required_hooks)
        {
            if (!Python_Plugin_Has_Callable(plugin, hook))
            {
                Python_Plugin_Contract_Error(
                    std::string("transactional Python plugin must define ") +
                    hook);
            }
        }
    }
    return capabilities;
}

static void Call_Python_Force_Transaction_Hook(const char* name)
{
    PyObject* plugin = Python_Plugin_Module();
    if (plugin == NULL)
    {
        Python_Plugin_Contract_Error(
            std::string("cannot call ") + name +
            " because the Python plugin module is unavailable");
    }
    PyObject* hook = PyObject_GetAttrString(plugin, name);
    if (hook == NULL || !PyCallable_Check(hook))
    {
        Py_XDECREF(hook);
        Python_Plugin_Contract_Error(std::string(name) +
                                     " is missing or is not callable");
    }
    PyObject* result = PyObject_CallObject(hook, NULL);
    Py_DECREF(hook);
    if (result == NULL)
    {
        Python_Plugin_Contract_Error(std::string(name) + " failed");
    }
    Py_DECREF(result);
}

PLUGIN_API void Begin_Force_Transaction()
{
    Call_Python_Force_Transaction_Hook("Begin_Force_Transaction");
}

PLUGIN_API void Commit_Force_Transaction()
{
    Call_Python_Force_Transaction_Hook("Commit_Force_Transaction");
}

PLUGIN_API void Rollback_Force_Transaction()
{
    Call_Python_Force_Transaction_Hook("Rollback_Force_Transaction");
}

PLUGIN_API void After_Initial()
{
    if (!is_initialized) return;
    if (PyRun_SimpleString(R"XYJ(
Sponge._refresh_backend_views()
if Sponge.neighbor_list is None:
    Sponge.neighbor_list = Sponge.NEIGHBOR_LIST()
if hasattr(sponge_pyplugin, "After_Initial"):
    sponge_pyplugin.After_Initial()
    )XYJ") != 0)
    {
        Python_Plugin_Contract_Error("After_Initial failed");
    }
}

PLUGIN_API void Set_Domain_Information(void* domain_info)
{
    (void)domain_info;
    if (!is_initialized) return;
    if (PyRun_SimpleString(R"XYJ(
Sponge.neighbor_list = Sponge.NEIGHBOR_LIST()
Sponge.dd = Sponge.DOMAIN_INFORMATION()
)XYJ") != 0)
    {
        Python_Plugin_Contract_Error("Set_Domain_Information failed");
    }
}

PLUGIN_API void Calculate_Force()
{
    if (!is_initialized) return;
    if (PyRun_SimpleString(R"XYJ(
Sponge._refresh_backend_views()
_force_result = sponge_pyplugin.Calculate_Force()
Sponge._apply_force_result(
    _force_result, require_result=Sponge.backend_name == "jax"
)
del _force_result
Sponge._refresh_backend_views()
    )XYJ") != 0)
    {
        Python_Plugin_Contract_Error("Calculate_Force failed");
    }
}

PLUGIN_API void Mdout_Print()
{
    if (!is_initialized) return;
    if (PyRun_SimpleString(R"XYJ(
Sponge._refresh_backend_views()
if hasattr(sponge_pyplugin, "Mdout_Print"):
    sponge_pyplugin.Mdout_Print()
    )XYJ") != 0)
    {
        Python_Plugin_Contract_Error("Mdout_Print failed");
    }
}

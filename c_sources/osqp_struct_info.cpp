#include <osqp.h>
#include "osqp_struct.h"


/*
 * Specialization of the struct wrapper for the OSQPInfo struct.
 */
template<>
void OSQPStructWrapper<OSQPInfo>::registerFields() {
    m_struct = static_cast<OSQPInfo*>(c_calloc(1, sizeof(OSQPInfo)));

    if(!m_struct)
        mexErrMsgTxt("Failed to allocate a OSQPInfo object.");

    /*
     * Register the mapping between struct field name and the info struct memory location
     */
    // Solver status
    m_structFields.push_back(new OSQPStructFieldCharArray(m_struct->status, 32,     "status"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->status_val,    "status_val"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->status_polish, "status_polish"));

    // Solution quality
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->obj_val,  "obj_val"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->prim_res, "prim_res"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->dual_res, "dual_res"));

    // Algorithm information
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->iter,           "iter"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->rho_updates,    "rho_updates"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->rho_estimate, "rho_estimate"));

    // Timing information
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->setup_time,  "setup_time"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->solve_time,  "solve_time"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->update_time, "update_time"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->polish_time, "polish_time"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->run_time,    "run_time"));
}


// Instantiate the OSQPInfo wrapper class
template class OSQPStructWrapper<OSQPInfo>;
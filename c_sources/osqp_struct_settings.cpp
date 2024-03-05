#include <osqp.h>
#include "osqp_struct.h"

/*
 * Specialization for the settings struct
 */
template<>
void OSQPStructWrapper<OSQPSettings>::registerFields() {
    m_struct = static_cast<OSQPSettings*>(c_calloc(1, sizeof(OSQPSettings)));

    if(!m_struct)
        mexErrMsgTxt("Failed to allocate a OSQPSettings object.");

    osqp_set_default_settings(m_struct);

    /*
     * Register the mapping between struct field name and the settings memory location
     */
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->device,        "device"));
    m_structFields.push_back(new OSQPStructField<enum osqp_linsys_solver_type>(&m_struct->linsys_solver, "linsys_solver"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->verbose,       "verbose"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->warm_starting, "warm_starting"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->scaling,       "scaling"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->polishing,     "polishing"));

    // ADMM parameters
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->rho,      "rho"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->rho_is_vec, "rho_is_vec"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->sigma,    "sigma"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->alpha,    "alpha"));

    // CG settings
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->cg_max_iter,          "cg_max_iter"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->cg_tol_reduction,     "cg_tol_reduction"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->cg_tol_fraction,    "cg_tol_fraction"));
    m_structFields.push_back(new OSQPStructField<osqp_precond_type>(&m_struct->cg_precond, "cg_precond"));

    // adaptive rho logic
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->adaptive_rho,             "adaptive_rho"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->adaptive_rho_interval,    "adaptive_rho_interval"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->adaptive_rho_fraction,  "adaptive_rho_fraction"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->adaptive_rho_tolerance, "adaptive_rho_tolerance"));

    // termination parameters
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->max_iter,           "max_iter"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->eps_abs,          "eps_abs"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->eps_rel,          "eps_rel"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->eps_prim_inf,     "eps_prim_inf"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->eps_dual_inf,     "eps_dual_inf"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->scaled_termination, "scaled_termination"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->check_termination,  "check_termination"));
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->time_limit,       "time_limit"));

    // polishing parameters
    m_structFields.push_back(new OSQPStructField<OSQPFloat>(&m_struct->delta,            "delta"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->polish_refine_iter, "polish_refine_iter"));
}


// Instantiate the OSQPSettings wrapper class
template class OSQPStructWrapper<OSQPSettings>;
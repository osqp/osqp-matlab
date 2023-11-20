#include <osqp.h>

#include "memory_matlab.h"
#include "settings_matlab.h"

#include <cstring>


void OSQPSettingsWrapper::registerFields() {
    m_settings = static_cast<OSQPSettings*>(c_calloc(1, sizeof(OSQPSettings)));

    if(!m_settings)
        mexErrMsgTxt("Failed to allocate a OSQPSettings object.");

    osqp_set_default_settings(m_settings);

    /*
     * Register the mapping between struct field name and the settings memory location
     */
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->device,        "device"));
    m_settingsFields.push_back(new OSQPSettingsField<enum osqp_linsys_solver_type>(&m_settings->linsys_solver, "linsys_solver"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->verbose,       "verbose"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->warm_starting, "warm_starting"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->scaling,       "scaling"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->polishing,     "polishing"));

    // ADMM parameters
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->rho,      "rho"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->rho_is_vec, "rho_is_vec"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->sigma,    "sigma"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->alpha,    "alpha"));

    // CG settings
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->cg_max_iter,          "cg_max_iter"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->cg_tol_reduction,     "cg_tol_reduction"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->cg_tol_fraction,    "cg_tol_fraction"));
    m_settingsFields.push_back(new OSQPSettingsField<osqp_precond_type>(&m_settings->cg_precond, "cg_precond"));

    // adaptive rho logic
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->adaptive_rho,             "adaptive_rho"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->adaptive_rho_interval,    "adaptive_rho_interval"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->adaptive_rho_fraction,  "adaptive_rho_fraction"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->adaptive_rho_tolerance, "adaptive_rho_tolerance"));

    // termination parameters
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->max_iter,           "max_iter"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->eps_abs,          "eps_abs"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->eps_rel,          "eps_rel"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->eps_prim_inf,     "eps_prim_inf"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->eps_dual_inf,     "eps_dual_inf"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->scaled_termination, "scaled_termination"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->check_termination,  "check_termination"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->time_limit,       "time_limit"));

    // polishing parameters
    m_settingsFields.push_back(new OSQPSettingsField<OSQPFloat>(&m_settings->delta,            "delta"));
    m_settingsFields.push_back(new OSQPSettingsField<OSQPInt>(&m_settings->polish_refine_iter, "polish_refine_iter"));
}


OSQPSettingsWrapper::~OSQPSettingsWrapper() {
    for(auto& s : m_settingsFields) {
        delete s;
    }
    
    c_free(m_settings);
}


mxArray* OSQPSettingsWrapper::GetMxStruct() {
    // No fields are added right now, they are added in the for loop when they are set
    mxArray* mxSettings = mxCreateStructMatrix(1, 1, 0, NULL);

    // Copy the current settings into the struct to return
    for(const auto& s : m_settingsFields) {
        s->ToMxStruct(mxSettings);
    }

    return mxSettings;
}


void OSQPSettingsWrapper::ParseMxStruct(const mxArray* aStruct) {
    for(const auto& s : m_settingsFields) {
        s->ToOSQPSettings(aStruct);
    }
}


OSQPSettings* OSQPSettingsWrapper::GetOSQPSettingsCopy() {
    // Allocate the default settings
    OSQPSettings* ret = static_cast<OSQPSettings*>(c_calloc(1, sizeof(OSQPSettings)));
    
    // Copy the current settings for their return
    std::memcpy(ret, m_settings, sizeof(OSQPSettings));

    return ret;
}


void OSQPSettingsWrapper::ParseOSQPSettings(const OSQPSettings* aSettings) {
    std::memcpy(m_settings, aSettings, sizeof(OSQPSettings));
}

#include <osqp.h>
#include "osqp_struct.h"

/*
 * Specialization for the codegen_defines struct
 */
template<>
void OSQPStructWrapper<OSQPCodegenDefines>::registerFields() {
    m_struct = static_cast<OSQPCodegenDefines*>(c_calloc(1, sizeof(OSQPCodegenDefines)));
    if(!m_struct)
        mexErrMsgTxt("Failed to allocate a OSQPCodegenDefines object.");

    osqp_set_default_codegen_defines(m_struct);

    /*
     * Register the mapping between struct field name and the settings memory location
     */
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->embedded_mode,     "embedded_mode"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->float_type,        "float_type"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->printing_enable,   "printing_enable"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->profiling_enable,  "profiling_enable"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->interrupt_enable,  "interrupt_enable"));
    m_structFields.push_back(new OSQPStructField<OSQPInt>(&m_struct->derivatives_enable,"derivatives_enable"));
}


// Instantiate the OSQPCodegenDefines wrapper class
template class OSQPStructWrapper<OSQPCodegenDefines>;
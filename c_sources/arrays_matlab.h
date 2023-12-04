#ifndef ARRAYS_MATLAB_H_
#define ARRAYS_MATLAB_H_

#include <osqp.h>

/**
 * Copy the data from one array to another provided array.
 */
template <typename outArr, typename inArr>
void copyVector(outArr* out, inArr* in, OSQPInt numel) {
    // Don't bother doing anything if there is no input data
    if(!in || !out || (numel == 0))
        return;
    
    // Copy the data
    for(OSQPInt i=0; i < numel; i++){
        out[i] = static_cast<outArr>(in[i]);
    }
}


/**
 * Copy the data from one array to another newly allocated array.
 * The caller gains ownership of the returned array.
 */
template <typename outArr, typename inArr>
outArr* cloneVector(inArr* in, OSQPInt numel) {
    // Don't bother doing anything if there is no input data
    if(!in || (numel == 0))
        return NULL;

    // Allocate new array
    outArr* out = static_cast<outArr*>(c_malloc(numel * sizeof(outArr)));

    if(!out)
        mexErrMsgTxt("Failed to allocate a vector object.");

    // Copy the data
    copyVector<outArr, inArr>(out, in, numel);
    return out;
}

#endif
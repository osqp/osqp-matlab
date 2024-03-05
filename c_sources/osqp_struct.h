#ifndef OSQP_STRUCT_H_
#define OSQP_STRUCT_H_

#include <cstring>
#include <functional>
#include <string>
#include <vector>

#include <mex.h>
#include <matrix.h>

#include "memory_matlab.h"
#include <osqp.h>

/**
 * Base class used to store the field types for a struct.
 */
class OSQPStructFieldBase {
public:
    OSQPStructFieldBase() {}

    /**
     * Set the field in the given Matlab struct to the value of this field
     */
    virtual void ToMxStruct(mxArray* aStruct) = 0;

    /**
     * Set the field in the internal struct with the data from aStruct
     */
    virtual void ToOSQPStruct(const mxArray* aStruct) = 0;
};

/**
 * Class to hold a numeric struct field (e.g. float/double/int/enum, etc.).
 */
template<class T>
class OSQPStructField : public OSQPStructFieldBase {
public:
    OSQPStructField(T* aStructPtr, std::string aName) :
        m_structPtr(aStructPtr),
        m_name(aName) {
    }

    void ToMxStruct(mxArray* aStruct) override {
        mxAddField(aStruct, m_name.data());
        mxSetField(aStruct, 0, m_name.data(), mxCreateDoubleScalar(*m_structPtr));
    }

    void ToOSQPStruct(const mxArray* aStruct) override {
        *(m_structPtr) = static_cast<T>(mxGetScalar(mxGetField(aStruct, 0, m_name.data())));
    }

private:
    T* m_structPtr;
    std::string m_name;
};

/**
 * Class to hold a character array (actual array, not char* array) field in a struct.
 */
class OSQPStructFieldCharArray : public OSQPStructFieldBase {
public:
    OSQPStructFieldCharArray(char* aStructPtr, size_t aLength, std::string aName) :
        m_structPtr(aStructPtr),
        m_name(aName),
        m_length(aLength) {
    }

    void ToMxStruct(mxArray* aStruct) override {
        mxAddField(aStruct, m_name.data());
        mxSetField(aStruct, 0, m_name.data(), mxCreateString(m_structPtr));
    }

    void ToOSQPStruct(const mxArray* aStruct) override {
        mxArray* tmp = mxGetField(aStruct, 0, m_name.data());
        mxGetString(tmp, m_structPtr, m_length);
    }

private:
    char* m_structPtr;
    std::string m_name;
    size_t m_length;
};

/**
 * Wrap a struct from OSQP to automatically transfer the data between OSQP and Matlab.
 */
template<class T>
class OSQPStructWrapper {
public:
    /**
     * Initialize the wrapper using the default values.
     */
    OSQPStructWrapper() {
        // Allocate the default struct and register field handlers
        registerFields();
    }

    /**
     * Initialize the wrapper using the values from the OSQP struct.
     */
    OSQPStructWrapper(const T* aStruct) {
        // Allocate the default struct and register field handlers
        registerFields();
        ParseOSQPStruct(aStruct);
    }

    /**
     * Initialize the wrapper using the values from the Matlab struct
     */
    OSQPStructWrapper(const mxArray* aStruct) {
        // Allocate the default struct and register field handlers
        registerFields();
        ParseMxStruct(aStruct);
    }

    ~OSQPStructWrapper() {
        for(auto& s : m_structFields) {
            delete s;
        }
        
        c_free(m_struct);
    }

    /**
     * Return a Matlab struct populated with the values of the current struct
     * contained in this wrapper.
     * 
     * @return a Matlab struct with a copy of the struct (caller owns this copy and must free it)
     */
    mxArray* GetMxStruct() {
        // No fields are added right now, they are added in the for loop when they are set
        mxArray* matStruct = mxCreateStructMatrix(1, 1, 0, NULL);

        // Copy the current struct into the struct to return
        for(const auto& s : m_structFields) {
            s->ToMxStruct(matStruct);
        }

        return matStruct;
    }

    /**
     * Read a Matlab struct and populate the wrapper with its values.
     */
    void ParseMxStruct(const mxArray* aStruct) {
        for(const auto& s : m_structFields) {
            s->ToOSQPStruct(aStruct);
        }
    }

    /**
     * Get a copy of the struct contained inside this wrapper.
     *
     * @return a copy of the struct (caller owns this copy and must free it)
     */
    T* GetOSQPStructCopy() {
        // Allocate the default struct
        T* ret = static_cast<T*>(c_calloc(1, sizeof(T)));
        
        // Copy the current values for their return
        std::memcpy(ret, m_struct, sizeof(T));
        return ret;
    }

    /**
     * Get the pointer to the internal struct object.
     */
    T* GetOSQPStruct() {
        return m_struct;
    }

    /*
     * Read an existing OSQP struct object into this wrapper.
     * The struct elements are copied, so no ownership of the aStruct pointer is transferred.
     */
    void ParseOSQPStruct(const T* aStruct) {
        std::memcpy(m_struct, aStruct, sizeof(T));
    }

private:
    /**
     * Register all the fields for the wrapper.
     * This function should be specialized for each struct type to map the fields appropriately.
     */
    void registerFields();

    // All struct fields
    std::vector<OSQPStructFieldBase*> m_structFields;

    // Base OSQP struct object. Owned by this wrapper.
    T* m_struct;
};

/**
 * Wrapper around the OSQPSettings struct
 */
typedef OSQPStructWrapper<OSQPSettings> OSQPSettingsWrapper;

/**
 * Wrapper around the OSQPInfo struct
 */
typedef OSQPStructWrapper<OSQPInfo> OSQPInfoWrapper;

#endif
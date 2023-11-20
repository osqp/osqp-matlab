#ifndef SETTINGS_MATLAB_H_
#define SETTINGS_MATLAB_H_

#include <functional>
#include <string>
#include <vector>

#include <mex.h>
#include <matrix.h>
#include <osqp.h>

/*
 * Base class used to store the templated settings field types.
 */
class OSQPSettingsFieldBase {
public:
    OSQPSettingsFieldBase() {}

    virtual void ToMxStruct(mxArray* aStruct) = 0;
    virtual void ToOSQPSettings(const mxArray* aStruct) = 0;
};

template<class T>
class OSQPSettingsField : public OSQPSettingsFieldBase {
public:
    OSQPSettingsField(T* aSettingPtr, std::string aName) :
        m_settingsPtr(aSettingPtr),
        m_name(aName) {
    }

    /*
     * Set the field in the given Matlab struct to the value of this settings field
     */
    void ToMxStruct(mxArray* aStruct) override {
        mxAddField(aStruct, m_name.data());
        mxSetField(aStruct, 0, m_name.data(), mxCreateDoubleScalar(*m_settingsPtr));
    }

    /*
     * Set the field in the internal OSQPSettings struct with the data from aStruct
     */
    void ToOSQPSettings(const mxArray* aStruct) override {
        *(m_settingsPtr) = static_cast<T>(mxGetScalar(mxGetField(aStruct, 0, m_name.data())));
    }

private:
    T* m_settingsPtr;
    std::string m_name;
};

class OSQPSettingsWrapper {
public:
    /*
     * Initialize the settings wrapper using the default settings.
     */
    OSQPSettingsWrapper() {
        // Allocate the default settings and register field handlers
        registerFields();
    }

    /*
     * Initialize the settings wrapper using the values from aSettings.
     */
    OSQPSettingsWrapper(const OSQPSettings* aSettings) {
        // Allocate the default settings and register field handlers
        registerFields();
        ParseOSQPSettings(aSettings);
    }

    /*
     * Initialize the settings wrapper using the values from aStruct
     */
    OSQPSettingsWrapper(const mxArray* aStruct) {
        // Allocate the default settings and register field handlers
        registerFields();
        ParseMxStruct(aStruct);
    }

    ~OSQPSettingsWrapper();

    /*
     * Return a Matlab structu populated with the values of the current settings
     * contained in this wrapper.
     * 
     * @return a Matlab struct with a copy of the settings (caller owns this copy and must free it)
     */
    mxArray* GetMxStruct();

    /*
     * Read a Matlab struct and populate the wrapper with its values.
     */
    void ParseMxStruct(const mxArray* aStruct);

    /*
     * Get a copy of the settings contained inside this wrapper.
     *
     * @return a copy of the settings (caller owns this copy and must free it)
     */
    OSQPSettings* GetOSQPSettingsCopy();

    /*
     * Get the pointer to the internal settings object.
     */
    OSQPSettings* GetOSQPSettings() {
        return m_settings;
    }

    /*
     * Read an existing OSQPSettings object into this wrapper.
     * The settings are copied, so no ownership of the aSettings pointer is transferred.
     */
    void ParseOSQPSettings(const OSQPSettings* aSettings);

private:
    // Register all the fields
    void registerFields();

    // All settings fields
    std::vector<OSQPSettingsFieldBase*> m_settingsFields;

    // Base OSQP settings object. Owned by this wrapper.
    OSQPSettings* m_settings;
};

#endif
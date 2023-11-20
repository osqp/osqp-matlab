#include <map>

#include "mex.h"
#include "matrix.h"
#include "osqp.h"

// Mex-specific functionality
#include "osqp_mex.hpp"
#include "memory_matlab.h"
#include "settings_matlab.h"

//c_int is replaced with OSQPInt
//c_float is replaced with OSQPFloat

//TODO: Check if this definition is required, and maybe replace it with:
//   enum linsys_solver_type { QDLDL_SOLVER, MKL_PARDISO_SOLVER };
#define QDLDL_SOLVER 0 //Based on the previous API

// all of the OSQP_INFO fieldnames as strings
const char* OSQP_INFO_FIELDS[] = {"status",         //char*     
                                  "status_val",     //OSQPInt 
                                  "status_polish",  //OSQPInt
                                  "obj_val",        //OSQPFloat
                                  "prim_res",       //OSQPFloat
                                  "dual_res",       //OSQPFloat
                                  "iter",           //OSQPInt
                                  "rho_updates",    //OSQPInt
                                  "rho_estimate",   //OSQPFloat
                                  "setup_time",     //OSQPFloat
                                  "solve_time",     //OSQPFloat
                                  "update_time",    //OSQPFloat
                                  "polish_time",    //OSQPFloat
                                  "run_time",       //OSQPFloat
                                  };      

#define NEW_SETTINGS_TOL (1e-10)

// wrapper class for all osqp data and settings
class OsqpData
{
public:
  OsqpData() : solver(NULL){}
  OSQPSolver     * solver;
};

// internal utility functions
OSQPSolver*    initializeOSQPSolver();
void           castToDoubleArr(OSQPFloat *arr, double* arr_out, OSQPInt len);
void           setToNaN(double* arr_out, OSQPInt len);
//void           castCintToDoubleArr(OSQPInt *arr, double* arr_out, OSQPInt len); //DELETE HERE
void           freeCscMatrix(OSQPCscMatrix* M);
OSQPInt*       copyToOSQPIntVector(mwIndex * vecData, OSQPInt numel);
OSQPInt*       copyDoubleToOSQPIntVector(double* vecData, OSQPInt numel);
OSQPFloat*     copyToOSQPFloatVector(double * vecData, OSQPInt numel);
mxArray*       copyInfoToMxStruct(OSQPInfo* info);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    OsqpData* osqpData;
    //OSQPSolver* osqpSolver = NULL;
    // Exitflag
    OSQPInt exitflag = 0;
    // Static string for static methods
    char stat_string[64];
    // Get the command string
    char cmd[64];
	  if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    // new object
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1){
            mexErrMsgTxt("New: One output expected.");
          }
        // Return a handle to a new C++ wrapper instance
        osqpData = new OsqpData;
        //osqpData->solver = initializeOSQPSolver();
        osqpData->solver = NULL;
        plhs[0] = convertPtr2Mat<OsqpData>(osqpData);
        return;
    }

    // Check for a second input, which should be the class instance handle or string 'static'
    if (nrhs < 2)
    		mexErrMsgTxt("Second input should be a class instance handle or the string 'static'.");

    if(mxGetString(prhs[1], stat_string, sizeof(stat_string))){
        // If we are dealing with non-static methods, get the class instance pointer from the second input
        osqpData = convertMat2Ptr<OsqpData>(prhs[1]);
    } else {
        if (strcmp("static", stat_string)){
            mexErrMsgTxt("Second argument for static functions is string 'static'");
        }
    }
    // delete the object and its data
    if (!strcmp("delete", cmd)) {
        
        osqp_cleanup(osqpData->solver);
        destroyObject<OsqpData>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }

    // report the current settings
    if (!strcmp("current_settings", cmd)) {
        // Throw an error if this is called before solver is configured
        if(!osqpData->solver) {
            mexErrMsgTxt("Solver is uninitialized.  No settings have been configured.");
        }
        if(!osqpData->solver->settings) {
            mexErrMsgTxt("Solver settings is uninitialized.  No settings have been configured.");
        }

        // Report the current settings
        OSQPSettingsWrapper settings(osqpData->solver->settings);
        plhs[0] = settings.GetMxStruct();
        return;
    }

    // write_settings
    if (!strcmp("update_settings", cmd)) {
        // Overwrite the current settings.  Mex function is responsible
        // for disallowing overwrite of selected settings after initialization,
        // and for all error checking
        // throw an error if this is called before solver is configured
        if(!osqpData->solver){
            mexErrMsgTxt("Solver is uninitialized.  No settings have been configured.");
        }

        OSQPSettingsWrapper settings(prhs[2]);
        osqp_update_settings(osqpData->solver, settings.GetOSQPSettings());
        return;
    }

    // Update rho value
    if (!strcmp("update_rho", cmd)) {
        //throw an error if this is called before solver is configured
        if(!osqpData->solver){
            mexErrMsgTxt("Solver is uninitialized.  No settings have been configured.");
        }
    
        OSQPFloat rho = (OSQPFloat)mxGetScalar(prhs[2]);

        osqp_update_rho(osqpData->solver, rho);
        return;
    }

    // report the default settings
    if (!strcmp("default_settings", cmd)) {
        // Warn if other commands were ignored
        if (nrhs > 2)
            mexWarnMsgTxt("Default settings: unexpected number of arguments.");


        // Create a Settings structure in default form and report the results
        // Useful for external solver packages (e.g. Yalmip) that want to
        // know which solver settings are supported
        OSQPSettingsWrapper settings;
        plhs[0] = settings.GetMxStruct();
        return;
    }

    // setup
    if (!strcmp("setup", cmd)) {
        //throw an error if this is called more than once
        if(osqpData->solver){
          mexErrMsgTxt("Solver is already initialized with problem data.");
        }

        // handle the problem data first.  Matlab-side
        // class wrapper is responsible for ensuring that
        // P and A are sparse matrices,  everything
        // else is a dense vector and all inputs are
        // of compatible dimension

        // Get mxArrays
        const mxArray* P  = prhs[4];
        const mxArray* q  = prhs[5];
        const mxArray* A  = prhs[6];
        const mxArray* l = prhs[7];
        const mxArray* u = prhs[8];


        OSQPInt dataN      = (OSQPInt)mxGetScalar(prhs[2]);
        OSQPInt dataM      = (OSQPInt)mxGetScalar(prhs[3]);
        OSQPFloat* dataQ   = copyToOSQPFloatVector(mxGetPr(q), dataN);
        OSQPFloat* dataL   = copyToOSQPFloatVector(mxGetPr(l), dataM);
        OSQPFloat* dataU   = copyToOSQPFloatVector(mxGetPr(u), dataM);

        // Matrix P:  nnz = P->p[n]
        OSQPInt * Pp = (OSQPInt*)copyToOSQPIntVector(mxGetJc(P), dataN + 1);
        OSQPInt * Pi = (OSQPInt*)copyToOSQPIntVector(mxGetIr(P), Pp[dataN]);
        OSQPFloat * Px = copyToOSQPFloatVector(mxGetPr(P), Pp[dataN]);
        OSQPCscMatrix* dataP = (OSQPCscMatrix*)c_calloc(1,sizeof(OSQPCscMatrix));
        csc_set_data(dataP, dataN, dataN, Pp[dataN], Px, Pi, Pp);

        // Matrix A: nnz = A->p[n]
        OSQPInt* Ap = (OSQPInt*)copyToOSQPIntVector(mxGetJc(A), dataN + 1);
        OSQPInt* Ai = (OSQPInt*)copyToOSQPIntVector(mxGetIr(A), Ap[dataN]);
        OSQPFloat * Ax = copyToOSQPFloatVector(mxGetPr(A), Ap[dataN]);
        OSQPCscMatrix* dataA = (OSQPCscMatrix*)c_calloc(1,sizeof(OSQPCscMatrix));
        csc_set_data(dataA, dataM, dataN, Ap[dataN], Ax, Ai, Ap);

        // Create Settings
        OSQPSettingsWrapper settings;

        if(!mxIsEmpty(prhs[9])){
          // Populate settings structure from mxArray input, otherwise the default settings are used
          settings.ParseMxStruct(prhs[9]);
        }

        // Setup workspace
        //exitflag = osqp_setup(&(osqpData->work), data, settings);
        exitflag = osqp_setup(&(osqpData->solver), dataP, dataQ, dataA, dataL, dataU, dataM, dataN, settings.GetOSQPSettings());
        //cleanup temporary structures
        // Data
        if (Px)       c_free(Px);
        if (Pi)       c_free(Pi);
        if (Pp)       c_free(Pp);
        if (Ax)       c_free(Ax);
        if (Ai)       c_free(Ai);
        if (Ap)       c_free(Ap);
        if (dataQ)    c_free(dataQ);
        if (dataL)    c_free(dataL);
        if (dataU)    c_free(dataU);
        if (dataP)    c_free(dataP);
        if (dataA)    c_free(dataA);

        // Report error (if any)
        if(exitflag){
            mexErrMsgTxt("Invalid problem setup");
        }

        return;

    }

    // get #constraints and variables
    if (!strcmp("get_dimensions", cmd)) {

        //throw an error if this is called before solver is configured
        if(!osqpData->solver){
          mexErrMsgTxt("Solver has not been initialized.");
        }
        OSQPInt n, m;
        osqp_get_dimensions(osqpData->solver, &m, &n);
        plhs[0] = mxCreateDoubleScalar(n);
        plhs[1] = mxCreateDoubleScalar(m);

        return;
    }

    // report the version
    if (!strcmp("version", cmd)) {

        plhs[0] = mxCreateString(osqp_version());

        return;
    }

    // update linear cost and bounds
    if (!strcmp("update", cmd)) {

        //throw an error if this is called before solver is configured
        if(!osqpData->solver){
          mexErrMsgTxt("Solver has not been initialized.");
        }

        // Fill q, l, u
        const mxArray *q = prhs[2];
        const mxArray *l = prhs[3];
        const mxArray *u = prhs[4];
        const mxArray *Px = prhs[5];
        const mxArray *Px_idx = prhs[6];
        const mxArray *Ax = prhs[8];
        const mxArray *Ax_idx = prhs[9];

        int Px_n, Ax_n;
        Px_n = mxGetScalar(prhs[7]);
        Ax_n = mxGetScalar(prhs[10]);

        // Copy vectors to ensure they are cast as OSQPFloat
        OSQPFloat *q_vec      = NULL;
        OSQPFloat *l_vec      = NULL;
        OSQPFloat *u_vec      = NULL;
        OSQPFloat *Px_vec     = NULL;
        OSQPFloat *Ax_vec     = NULL;
        OSQPInt *Px_idx_vec   = NULL;
        OSQPInt *Ax_idx_vec   = NULL;

        OSQPInt n, m;
        osqp_get_dimensions(osqpData->solver, &m, &n);
        if(!mxIsEmpty(q)){
            q_vec = copyToOSQPFloatVector(mxGetPr(q), n);
        }
        if(!mxIsEmpty(l)){
            l_vec = copyToOSQPFloatVector(mxGetPr(l), m);
        }
        if(!mxIsEmpty(u)){
            u_vec = copyToOSQPFloatVector(mxGetPr(u), m);
        }
        if(!mxIsEmpty(Px)){
            Px_vec = copyToOSQPFloatVector(mxGetPr(Px), Px_n);
        }
        if(!mxIsEmpty(Ax)){
            Ax_vec = copyToOSQPFloatVector(mxGetPr(Ax), Ax_n);
        }
        if(!mxIsEmpty(Px_idx)){
            Px_idx_vec = copyDoubleToOSQPIntVector(mxGetPr(Px_idx), Px_n);
        }
        if(!mxIsEmpty(Ax_idx)){
            Ax_idx_vec = copyDoubleToOSQPIntVector(mxGetPr(Ax_idx), Ax_n);
        }

        if (!exitflag && (!mxIsEmpty(q) || !mxIsEmpty(l) || !mxIsEmpty(u))) {
          exitflag = osqp_update_data_vec(osqpData->solver, q_vec, l_vec, u_vec);
          if (exitflag) exitflag=1; 
        }
   
        if (!exitflag && (!mxIsEmpty(Px) || !mxIsEmpty(Ax))) {
          exitflag = osqp_update_data_mat(osqpData->solver, Px_vec, Px_idx_vec, Px_n, Ax_vec, Ax_idx_vec, Ax_n);
          if (exitflag) exitflag=2;
        }
                                      

        // Free vectors
        if(!mxIsEmpty(q))  c_free(q_vec);
        if(!mxIsEmpty(l))  c_free(l_vec);
        if(!mxIsEmpty(u))  c_free(u_vec);
        if(!mxIsEmpty(Px)) c_free(Px_vec);
        if(!mxIsEmpty(Ax)) c_free(Ax_vec);
        if(!mxIsEmpty(Px_idx)) c_free(Px_idx_vec);
        if(!mxIsEmpty(Ax_idx)) c_free(Ax_idx_vec);

        // Report errors (if any)
        switch (exitflag) {
            case 1:
                mexErrMsgTxt("Data update error!");
            case 2:
                mexErrMsgTxt("Matrix update error!");
        }

        return;
    }

    if (!strcmp("warm_start", cmd) || !strcmp("warm_start_x", cmd) || !strcmp("warm_start_y", cmd)) {
      
      //throw an error if this is called before solver is configured
      if(!osqpData->solver){
          mexErrMsgTxt("Solver has not been initialized.");
        }

      // Fill x and y
      const mxArray *x = NULL;
      const mxArray *y = NULL;
      if (!strcmp("warm_start", cmd)) {  
        x = prhs[2];
        y = prhs[3];
      }
      else if (!strcmp("warm_start_x", cmd)) {
        x = prhs[2];
        y = NULL;
      }

      else if (!strcmp("warm_start_y", cmd)) {
        x = NULL;
        y = prhs[2];
      }

      // Copy vectors to ensure they are cast as OSQPFloat
      OSQPFloat *x_vec = NULL;
      OSQPFloat *y_vec = NULL;

      OSQPInt n, m;
      osqp_get_dimensions(osqpData->solver, &m, &n);
      if(!mxIsEmpty(x)){
          x_vec = copyToOSQPFloatVector(mxGetPr(x),n);
      }
      if(!mxIsEmpty(y)){
          y_vec = copyToOSQPFloatVector(mxGetPr(y),m);
      }

      // Warm start x and y
      osqp_warm_start(osqpData->solver, x_vec, y_vec);

      // Free vectors
      if(!x_vec) c_free(x_vec);
      if(!y_vec) c_free(y_vec);

      return;
    }

    // SOLVE
    if (!strcmp("solve", cmd)) {
        if (nlhs != 5 || nrhs != 2){
          mexErrMsgTxt("Solve : wrong number of inputs / outputs");
        }
        if(!osqpData->solver){
            mexErrMsgTxt("No problem data has been given.");
        }
        // solve the problem
        osqp_solve(osqpData->solver);

        OSQPInt n, m;
        osqp_get_dimensions(osqpData->solver, &m, &n);
        // Allocate space for solution
        // primal variables
        plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
        // dual variables
        plhs[1] = mxCreateDoubleMatrix(m,1,mxREAL);
        // primal infeasibility certificate
        plhs[2] = mxCreateDoubleMatrix(m,1,mxREAL);
        // dual infeasibility certificate
        plhs[3] = mxCreateDoubleMatrix(n,1,mxREAL);

        //copy results to mxArray outputs
        //assume that five outputs will always
        //be returned to matlab-side class wrapper
        if ((osqpData->solver->info->status_val != OSQP_PRIMAL_INFEASIBLE) &&
            (osqpData->solver->info->status_val != OSQP_DUAL_INFEASIBLE)){

            //primal and dual solutions
            castToDoubleArr(osqpData->solver->solution->x, mxGetPr(plhs[0]), n);
            castToDoubleArr(osqpData->solver->solution->y, mxGetPr(plhs[1]), m);

            //infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), m);
            setToNaN(mxGetPr(plhs[3]), n);

        } else if (osqpData->solver->info->status_val == OSQP_PRIMAL_INFEASIBLE ||
        osqpData->solver->info->status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE){ //primal infeasible

            //primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), n);
            setToNaN(mxGetPr(plhs[1]), m);

            //primal infeasibility certificates
            castToDoubleArr(osqpData->solver->solution->prim_inf_cert, mxGetPr(plhs[2]), m);

            //dual infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[3]), n);

            // Set objective value to infinity
            osqpData->solver->info->obj_val = mxGetInf();

        } else { //dual infeasible

            //primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), n);
            setToNaN(mxGetPr(plhs[1]), m);

            //primal infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), m);

            //dual infeasibility certificates
            castToDoubleArr(osqpData->solver->solution->dual_inf_cert, mxGetPr(plhs[3]), n);

            // Set objective value to -infinity
            osqpData->solver->info->obj_val = -mxGetInf();
        }

        if (osqpData->solver->info->status_val == OSQP_NON_CVX) {
            osqpData->solver->info->obj_val = mxGetNaN();
        }

        plhs[4] = copyInfoToMxStruct(osqpData->solver->info); // Info structure

        return;
    }

    if (!strcmp("constant", cmd)) { // Return solver constants
        static std::map<std::string, OSQPFloat> floatConstants{
            // Numerical constants
            {"OSQP_INFTY", OSQP_INFTY}
        };

        static std::map<std::string, OSQPInt> intConstants{
            // Return codes
            {"OSQP_SOLVED",                       OSQP_SOLVED},
            {"OSQP_SOLVED_INACCURATE",            OSQP_SOLVED_INACCURATE},
            {"OSQP_UNSOLVED",                     OSQP_UNSOLVED},
            {"OSQP_PRIMAL_INFEASIBLE",            OSQP_PRIMAL_INFEASIBLE},
            {"OSQP_PRIMAL_INFEASIBLE_INACCURATE", OSQP_PRIMAL_INFEASIBLE_INACCURATE},
            {"OSQP_DUAL_INFEASIBLE",              OSQP_DUAL_INFEASIBLE},
            {"OSQP_DUAL_INFEASIBLE_INACCURATE",   OSQP_DUAL_INFEASIBLE_INACCURATE},
            {"OSQP_MAX_ITER_REACHED",             OSQP_MAX_ITER_REACHED},
            {"OSQP_NON_CVX",                      OSQP_NON_CVX},
            {"OSQP_TIME_LIMIT_REACHED",           OSQP_TIME_LIMIT_REACHED},

            // Linear system solvers
            {"QDLDL_SOLVER",         QDLDL_SOLVER},
            {"OSQP_UNKNOWN_SOLVER",  OSQP_UNKNOWN_SOLVER},
            {"OSQP_DIRECT_SOLVER",   OSQP_DIRECT_SOLVER},
            {"OSQP_INDIRECT_SOLVER", OSQP_INDIRECT_SOLVER}
        };

        char constant[64];
        int  constantLength = mxGetN(prhs[2]) + 1;
        mxGetString(prhs[2], constant, constantLength);

        auto ci = intConstants.find(constant);

        if(ci != intConstants.end()) {
            plhs[0] = mxCreateDoubleScalar(ci->second);
            return;
        }

        auto cf = floatConstants.find(constant);

        if(cf != floatConstants.end()) {
            plhs[0] = mxCreateDoubleScalar(cf->second);
            return;
        }

        // NaN is special because we need the Matlab version
        if (!strcmp("OSQP_NAN", constant)){
            plhs[0] = mxCreateDoubleScalar(mxGetNaN());
            return;
        }

        mexErrMsgTxt("Constant not recognized.");

        return;
    }

    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}

/**
 * This function dynamically allocates OSQPSovler and sets all the properties of OSQPSolver to NULL.
 * WARNING: The memory allocated here (OSQPSolver*) needs to be freed.
 * WARNING: Any dynamically allocated pointers must be freed before calling this function.
*/
OSQPSolver* initializeOSQPSolver() {
  OSQPSolver* osqpSolver = new OSQPSolver;
  osqpSolver->info       = NULL;
  osqpSolver->settings   = NULL;
  osqpSolver->solution   = NULL;
  osqpSolver->work       = NULL;
  //osqp_set_default_settings(osqpSolver->settings);
  return osqpSolver;
}

//Dynamically creates a OSQPFloat vector copy of the input.
//Returns an empty pointer if vecData is NULL
OSQPFloat*    copyToOSQPFloatVector(double * vecData, OSQPInt numel){
  if (!vecData) return NULL;

  //This needs to be freed!
  OSQPFloat* out = (OSQPFloat*)c_malloc(numel * sizeof(OSQPFloat));

  //copy data
  for(OSQPInt i=0; i < numel; i++){
      out[i] = (OSQPFloat)vecData[i];
  }
  return out;
}

//Dynamically creates a OSQPInt vector copy of the input.
OSQPInt* copyToOSQPIntVector(mwIndex* vecData, OSQPInt numel){
  // This memory needs to be freed!
  OSQPInt* out = (OSQPInt*)c_malloc(numel * sizeof(OSQPInt));

  //copy data
  for(OSQPInt i=0; i < numel; i++){
      out[i] = (OSQPInt)vecData[i];
  }
  return out;

}

//Dynamically copies a double vector to OSQPInt.
OSQPInt* copyDoubleToOSQPIntVector(double* vecData, OSQPInt numel){
  // This memory needs to be freed!
  OSQPInt* out = (OSQPInt*)c_malloc(numel * sizeof(OSQPInt));

  //copy data
  for(OSQPInt i=0; i < numel; i++){
      out[i] = (OSQPInt)vecData[i];
  }
  return out;

}

/* DELETE HERE
void castCintToDoubleArr(OSQPInt *arr, double* arr_out, OSQPInt len) {
    for (OSQPInt i = 0; i < len; i++) {
        arr_out[i] = (double)arr[i];
    }
}*/

//This function frees the memory allocated in an OSQPCscMatrix M
void freeCscMatrix(OSQPCscMatrix* M) {
    if (!M) return;
    if (M->p) c_free(M->p);
    if (M->i) c_free(M->i);
    if (M->x) c_free(M->x);
    c_free(M);
}

void castToDoubleArr(OSQPFloat *arr, double* arr_out, OSQPInt len) {
    for (OSQPInt i = 0; i < len; i++) {
        arr_out[i] = (double)arr[i];
    }
}

void setToNaN(double* arr_out, OSQPInt len){
    OSQPInt i;
    for (i = 0; i < len; i++) {
        arr_out[i] = mxGetNaN();
    }
}

mxArray* copyInfoToMxStruct(OSQPInfo* info){

  //create mxArray with the right number of fields
  int nfields  = sizeof(OSQP_INFO_FIELDS) / sizeof(OSQP_INFO_FIELDS[0]);
  mxArray* mxPtr = mxCreateStructMatrix(1,1,nfields,OSQP_INFO_FIELDS);

  //map the OSQP_INFO fields one at a time into mxArrays
  //matlab all numeric values as doubles
  mxSetField(mxPtr, 0, "iter",          mxCreateDoubleScalar(info->iter));
  mxSetField(mxPtr, 0, "status",        mxCreateString(info->status));
  mxSetField(mxPtr, 0, "status_val",    mxCreateDoubleScalar(info->status_val));
  mxSetField(mxPtr, 0, "status_polish", mxCreateDoubleScalar(info->status_polish));
  mxSetField(mxPtr, 0, "obj_val",       mxCreateDoubleScalar(info->obj_val));
  mxSetField(mxPtr, 0, "prim_res",      mxCreateDoubleScalar(info->prim_res));
  mxSetField(mxPtr, 0, "dual_res",      mxCreateDoubleScalar(info->dual_res));

  mxSetField(mxPtr, 0, "setup_time",  mxCreateDoubleScalar(info->setup_time));
  mxSetField(mxPtr, 0, "solve_time",  mxCreateDoubleScalar(info->solve_time));
  mxSetField(mxPtr, 0, "update_time", mxCreateDoubleScalar(info->update_time));
  mxSetField(mxPtr, 0, "polish_time", mxCreateDoubleScalar(info->polish_time));
  mxSetField(mxPtr, 0, "run_time",    mxCreateDoubleScalar(info->run_time));


  mxSetField(mxPtr, 0, "rho_updates",    mxCreateDoubleScalar(info->rho_updates));
  mxSetField(mxPtr, 0, "rho_estimate",   mxCreateDoubleScalar(info->rho_estimate));


  return mxPtr;

}

#include "mex.h"
#include "matrix.h"
#include "osqp_mex.hpp"
#include "osqp.h"
#include "memory_matlab.h"

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

const char* OSQP_SETTINGS_FIELDS[] = {"device",                 //OSQPInt
                                      "linsys_solver",          //enum osqp_linsys_solver_type
                                      "verbose",                //OSQPInt    
                                      "warm_starting",          //OSQPInt
                                      "scaling",                //OSQPInt
                                      "polishing",              //OSQPInt
                                      "rho",                    //OSQPFloat
                                      "rho_is_vec",             //OSQPInt
                                      "sigma",                  //OSQPFloat
                                      "alpha",                  //OSQPFloat
                                      "cg_max_iter",            //OSQPInt
                                      "cg_tol_reduction",       //OSQPInt
                                      "cg_tol_fraction",        //OSQPFloat
                                      "cg_precond",             //osqp_precond_type
                                      "adaptive_rho",           //OSQPInt
                                      "adaptive_rho_interval",  //OSQPInt
                                      "adaptive_rho_fraction",  //OSQPFloat
                                      "adaptive_rho_tolerance", //OSQPFloat
                                      "max_iter",               //OSQPInt
                                      "eps_abs",                //OSQPFloat
                                      "eps_rel",                //OSQPFloat
                                      "eps_prim_inf",           //OSQPFloat
                                      "eps_dual_inf",           //OSQPFloat
                                      "scaled_termination",     //OSQPInt
                                      "check_termination",      //OSQPInt
                                      "time_limit",             //OSQPFloat
                                      "delta",                  //OSQPFloat
                                      "polish_refine_iter",     //OSQPInt
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
void           copyMxStructToSettings(const mxArray*, OSQPSettings*);
void           copyUpdatedSettingsToWork(const mxArray*, OSQPSolver*);
//void           castCintToDoubleArr(OSQPInt *arr, double* arr_out, OSQPInt len); //DELETE HERE
void           freeCscMatrix(OSQPCscMatrix* M);
OSQPInt*       copyToOSQPIntVector(mwIndex * vecData, OSQPInt numel);
OSQPInt*       copyDoubleToOSQPIntVector(double* vecData, OSQPInt numel);
OSQPFloat*     copyToOSQPFloatVector(double * vecData, OSQPInt numel);
mxArray*       copyInfoToMxStruct(OSQPInfo* info);
mxArray*       copySettingsToMxStruct(OSQPSettings* settings);


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
      //throw an error if this is called before solver is configured
      if(!osqpData->solver) mexErrMsgTxt("Solver is uninitialized.  No settings have been configured.");
      if(!osqpData->solver->settings){
        mexErrMsgTxt("Solver settings is uninitialized.  No settings have been configured.");
      }
      //report the current settings
      plhs[0] = copySettingsToMxStruct(osqpData->solver->settings);
      return;
    }

    // write_settings
    if (!strcmp("update_settings", cmd)) {
      //overwrite the current settings.  Mex function is responsible
      //for disallowing overwrite of selected settings after initialization,
      //and for all error checking
      //throw an error if this is called before solver is configured
      if(!osqpData->solver){
        mexErrMsgTxt("Solver is uninitialized.  No settings have been configured.");
      }

      copyUpdatedSettingsToWork(prhs[2],osqpData->solver);
      return;
    }

    // report the default settings
    if (!strcmp("default_settings", cmd)) {
        // Warn if other commands were ignored
        if (nrhs > 2)
            mexWarnMsgTxt("Default settings: unexpected number of arguments.");


      //Create a Settings structure in default form and report the results
      //Useful for external solver packages (e.g. Yalmip) that want to
      //know which solver settings are supported
      OSQPSettings* defaults = (OSQPSettings *)mxCalloc(1,sizeof(OSQPSettings));
      osqp_set_default_settings(defaults);
      plhs[0] = copySettingsToMxStruct(defaults);
      mxFree(defaults);
      return;
    }

    // setup
    if (!strcmp("setup", cmd)) {
        //throw an error if this is called more than once
        if(osqpData->solver){
          mexErrMsgTxt("Solver is already initialized with problem data.");
        }
        //Create data and settings containers
        OSQPSettings* settings = (OSQPSettings *)mxCalloc(1,sizeof(OSQPSettings));

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
        OSQPCscMatrix* dataP = new OSQPCscMatrix;
        csc_set_data(dataP, dataN, dataN, Pp[dataN], Px, Pi, Pp);

        // Matrix A: nnz = A->p[n]
        OSQPInt* Ap = (OSQPInt*)copyToOSQPIntVector(mxGetJc(A), dataN + 1);
        OSQPInt* Ai = (OSQPInt*)copyToOSQPIntVector(mxGetIr(A), Ap[dataN]);
        OSQPFloat * Ax = copyToOSQPFloatVector(mxGetPr(A), Ap[dataN]);
        OSQPCscMatrix* dataA = new OSQPCscMatrix;
        csc_set_data(dataA, dataM, dataN, Ap[dataN], Ax, Ai, Ap);

        // Create Settings
        const mxArray* mxSettings = prhs[9];
        if(mxIsEmpty(mxSettings)){
          // use defaults
          osqp_set_default_settings(settings);
        } else {
          //populate settings structure from mxArray input
          copyMxStructToSettings(mxSettings, settings);
        }

        // Setup workspace
        //exitflag = osqp_setup(&(osqpData->work), data, settings);
        exitflag = osqp_setup(&(osqpData->solver), dataP, dataQ, dataA, dataL, dataU, dataM, dataN, settings);
        //cleanup temporary structures
        // Data
        if (Px)       free(Px);
        if (Pi)       free(Pi);
        if (Pp)       free(Pp);
        if (Ax)       free(Ax);
        if (Ai)       free(Ai);
        if (Ap)       free(Ap);
        if (dataQ)    free(dataQ);
        if (dataL)    free(dataL);
        if (dataU)    free(dataU);
        if (dataP)    free(dataP);
        if (dataA)    free(dataA);
        // Settings
        if (settings) mxFree(settings);

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
        if(!mxIsEmpty(q))  free(q_vec);
        if(!mxIsEmpty(l))  free(l_vec);
        if(!mxIsEmpty(u))  free(u_vec);
        if(!mxIsEmpty(Px)) free(Px_vec);
        if(!mxIsEmpty(Ax)) free(Ax_vec);
        if(!mxIsEmpty(Px_idx)) free(Px_idx_vec);
        if(!mxIsEmpty(Ax_idx)) free(Ax_idx_vec);

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
      if(!x_vec) free(x_vec);
      if(!y_vec) free(y_vec);

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

        char constant[32];
        mxGetString(prhs[2], constant, sizeof(constant));

        if (!strcmp("OSQP_INFTY", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_INFTY);
            return;
        }
        if (!strcmp("OSQP_NAN", constant)){
            plhs[0] = mxCreateDoubleScalar(mxGetNaN());
            return;
        }

        if (!strcmp("OSQP_SOLVED", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_SOLVED);
            return;
        }

        if (!strcmp("OSQP_SOLVED_INACCURATE", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_SOLVED_INACCURATE);
            return;
        }

        if (!strcmp("OSQP_UNSOLVED", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_UNSOLVED);
            return;
        }

        if (!strcmp("OSQP_PRIMAL_INFEASIBLE", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_PRIMAL_INFEASIBLE);
            return;
        }

        if (!strcmp("OSQP_PRIMAL_INFEASIBLE_INACCURATE", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_PRIMAL_INFEASIBLE_INACCURATE);
            return;
        }

        if (!strcmp("OSQP_DUAL_INFEASIBLE", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_DUAL_INFEASIBLE);
            return;
        }

        if (!strcmp("OSQP_DUAL_INFEASIBLE_INACCURATE", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_DUAL_INFEASIBLE_INACCURATE);
            return;
        }

        if (!strcmp("OSQP_MAX_ITER_REACHED", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_MAX_ITER_REACHED);
            return;
        }

        if (!strcmp("OSQP_NON_CVX", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_NON_CVX);
            return;
        }

        if (!strcmp("OSQP_TIME_LIMIT_REACHED", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_TIME_LIMIT_REACHED);
            return;
        }

        // Linear system solvers
        if (!strcmp("QDLDL_SOLVER", constant)){
            plhs[0] = mxCreateDoubleScalar(QDLDL_SOLVER);
            return;
        }

        if (!strcmp("OSQP_UNKNOWN_SOLVER", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_UNKNOWN_SOLVER);
            return;
        }

        if (!strcmp("OSQP_DIRECT_SOLVER", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_DIRECT_SOLVER);
            return;
        }
        
        if (!strcmp("OSQP_INDIRECT_SOLVER", constant)){
            plhs[0] = mxCreateDoubleScalar(OSQP_INDIRECT_SOLVER);
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
  OSQPFloat* out = (OSQPFloat*)malloc(numel * sizeof(OSQPFloat));

  //copy data
  for(OSQPInt i=0; i < numel; i++){
      out[i] = (OSQPFloat)vecData[i];
  }
  return out;
}

//Dynamically creates a OSQPInt vector copy of the input.
OSQPInt* copyToOSQPIntVector(mwIndex* vecData, OSQPInt numel){
  // This memory needs to be freed!
  OSQPInt* out = (OSQPInt*)malloc(numel * sizeof(OSQPInt));

  //copy data
  for(OSQPInt i=0; i < numel; i++){
      out[i] = (OSQPInt)vecData[i];
  }
  return out;

}

//Dynamically copies a double vector to OSQPInt.
OSQPInt* copyDoubleToOSQPIntVector(double* vecData, OSQPInt numel){
  // This memory needs to be freed!
  OSQPInt* out = (OSQPInt*)malloc(numel * sizeof(OSQPInt));

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
    if (M->p) free(M->p);
    if (M->i) free(M->i);
    if (M->x) free(M->x);
    free(M);
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

mxArray* copySettingsToMxStruct(OSQPSettings* settings){

  int nfields  = sizeof(OSQP_SETTINGS_FIELDS) / sizeof(OSQP_SETTINGS_FIELDS[0]);
  mxArray* mxPtr = mxCreateStructMatrix(1,1,nfields,OSQP_SETTINGS_FIELDS);

  //map the OSQP_SETTINGS fields one at a time into mxArrays
  //matlab handles everything as a double
  mxSetField(mxPtr, 0, "rho",                    mxCreateDoubleScalar(settings->rho));
  mxSetField(mxPtr, 0, "sigma",                  mxCreateDoubleScalar(settings->sigma));
  mxSetField(mxPtr, 0, "scaling",                mxCreateDoubleScalar(settings->scaling));
  mxSetField(mxPtr, 0, "adaptive_rho",           mxCreateDoubleScalar(settings->adaptive_rho));
  mxSetField(mxPtr, 0, "adaptive_rho_interval",  mxCreateDoubleScalar(settings->adaptive_rho_interval));
  mxSetField(mxPtr, 0, "adaptive_rho_tolerance", mxCreateDoubleScalar(settings->adaptive_rho_tolerance));
  mxSetField(mxPtr, 0, "adaptive_rho_fraction",  mxCreateDoubleScalar(settings->adaptive_rho_fraction));
  mxSetField(mxPtr, 0, "max_iter",               mxCreateDoubleScalar(settings->max_iter));
  mxSetField(mxPtr, 0, "eps_abs",                mxCreateDoubleScalar(settings->eps_abs));
  mxSetField(mxPtr, 0, "eps_rel",                mxCreateDoubleScalar(settings->eps_rel));
  mxSetField(mxPtr, 0, "eps_prim_inf",           mxCreateDoubleScalar(settings->eps_prim_inf));
  mxSetField(mxPtr, 0, "eps_dual_inf",           mxCreateDoubleScalar(settings->eps_dual_inf));
  mxSetField(mxPtr, 0, "alpha",                  mxCreateDoubleScalar(settings->alpha));
  mxSetField(mxPtr, 0, "linsys_solver",          mxCreateDoubleScalar(settings->linsys_solver));
  mxSetField(mxPtr, 0, "delta",                  mxCreateDoubleScalar(settings->delta));
  mxSetField(mxPtr, 0, "polish_refine_iter",     mxCreateDoubleScalar(settings->polish_refine_iter));
  mxSetField(mxPtr, 0, "verbose",                mxCreateDoubleScalar(settings->verbose));
  mxSetField(mxPtr, 0, "scaled_termination",     mxCreateDoubleScalar(settings->scaled_termination));
  mxSetField(mxPtr, 0, "check_termination",      mxCreateDoubleScalar(settings->check_termination));
  mxSetField(mxPtr, 0, "warm_starting",          mxCreateDoubleScalar(settings->warm_starting));
  mxSetField(mxPtr, 0, "time_limit",             mxCreateDoubleScalar(settings->time_limit));
  mxSetField(mxPtr, 0, "device",                 mxCreateDoubleScalar(settings->device));
  mxSetField(mxPtr, 0, "polishing",              mxCreateDoubleScalar(settings->polishing));
  mxSetField(mxPtr, 0, "rho_is_vec",             mxCreateDoubleScalar(settings->rho_is_vec));
  mxSetField(mxPtr, 0, "cg_max_iter",            mxCreateDoubleScalar(settings->cg_max_iter));
  mxSetField(mxPtr, 0, "cg_tol_reduction",       mxCreateDoubleScalar(settings->cg_tol_reduction));
  mxSetField(mxPtr, 0, "cg_tol_fraction",        mxCreateDoubleScalar(settings->cg_tol_fraction));
  mxSetField(mxPtr, 0, "time_limit",             mxCreateDoubleScalar(settings->time_limit));
  return mxPtr;
}


// ======================================================================

void copyMxStructToSettings(const mxArray* mxPtr, OSQPSettings* settings){

  //this function assumes that only a complete and validated structure
  //will be passed.  matlab mex-side function is responsible for checking
  //structure validity

  //map the OSQP_SETTINGS fields one at a time into mxArrays
  //matlab handles everything as a double
  settings->rho                    = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "rho"));
  settings->sigma                  = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "sigma"));
  settings->scaling                = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "scaling"));
  settings->adaptive_rho           = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho"));
  settings->adaptive_rho_interval  = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho_interval"));
  settings->adaptive_rho_tolerance = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho_tolerance"));
  settings->adaptive_rho_fraction  = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "adaptive_rho_fraction"));
  settings->max_iter               = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "max_iter"));
  settings->eps_abs                = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs"));
  settings->eps_rel                = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_rel"));
  settings->eps_prim_inf           = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_dual_inf"));
  settings->eps_dual_inf           = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_dual_inf"));
  settings->alpha                  = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "alpha"));
  settings->linsys_solver          = (enum osqp_linsys_solver_type) (OSQPInt) mxGetScalar(mxGetField(mxPtr, 0, "linsys_solver"));
  settings->delta                  = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "delta"));
  settings->polish_refine_iter     = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "polish_refine_iter"));
  settings->verbose                = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "verbose"));
  settings->scaled_termination     = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "scaled_termination"));
  settings->check_termination      = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "check_termination"));
  settings->warm_starting          = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "warm_starting"));
  settings->time_limit             = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "time_limit"));
  settings->device                 = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "device"));
  settings->polishing              = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "polishing"));
  settings->rho_is_vec             = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "rho_is_vec"));
  settings->cg_max_iter            = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "cg_max_iter"));
  settings->cg_tol_reduction       = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "cg_tol_reduction"));
  settings->cg_tol_fraction        = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "cg_tol_fraction"));
  settings->cg_precond             = (osqp_precond_type) (OSQPInt) (mxGetField(mxPtr, 0, "cg_precond"));
}

void copyUpdatedSettingsToWork(const mxArray* mxPtr ,OSQPSolver* osqpSolver){

  OSQPInt exitflag;
  //TODO (Amit): Update this
  OSQPSettings* update_template = (OSQPSettings *)mxCalloc(1,sizeof(OSQPSettings));
  if (!update_template) mexErrMsgTxt("Failed to allocate temporary OSQPSettings object.");

  update_template->max_iter = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "max_iter"));
  update_template->eps_abs = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs"));
  update_template->eps_rel = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_rel"));
  update_template->eps_prim_inf = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_prim_inf"));
  update_template->eps_dual_inf = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_dual_inf"));
  update_template->alpha = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "alpha"));
  update_template->delta = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "delta"));
  update_template->polish_refine_iter = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "polish_refine_iter"));
  update_template->verbose = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "verbose"));
  update_template->scaled_termination = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "scaled_termination"));
  update_template->check_termination = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "check_termination"));
  update_template->warm_starting = (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "warm_starting"));
  update_template->time_limit = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "time_limit"));
  update_template->rho = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "rho"));

  osqp_update_settings(osqpSolver, update_template);
  if (update_template) free(update_template);
}
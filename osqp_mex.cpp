#include "mex.h"
#include "matrix.h"
#include "osqp_mex.hpp"
#include "osqp.h"
//#include "ctrlc.h"             // Needed for interrupt
#include "qdldl_interface.h"   // To extract workspace for codegen
//#include <string> //DELETE HERE
//#include <iostream> //DELETE HERE
//using std::cout; //DELETE HERE
//using std::endl; //DELETE HERE

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
                                  #ifdef PROFILING
                                  "setup_time",     //OSQPFloat
                                  "solve_time",     //OSQPFloat
                                  "update_time",    //OSQPFloat
                                  "polish_time",    //OSQPFloat
                                  "run_time",       //OSQPFloat
                                  #endif /* ifdef PROFILING */
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
void           castCintToDoubleArr(OSQPInt *arr, double* arr_out, OSQPInt len);
void           freeCscMatrix(OSQPCscMatrix* M);
OSQPInt*       copyToCintVector(mwIndex * vecData, OSQPInt numel);
OSQPInt*       copyDoubleToCintVector(double* vecData, OSQPInt numel);
OSQPFloat*     copyToOSQPFloatVector(double * vecData, OSQPInt numel);
mxArray*       copyInfoToMxStruct(OSQPInfo* info);
mxArray*       copySettingsToMxStruct(OSQPSettings* settings);
OSQPInt        osqp_update_max_iter(OSQPSolver* osqpSolver, OSQPInt max_iter_new);
OSQPInt        osqp_update_eps_abs(OSQPSolver* osqpSolver, OSQPFloat eps_abs_new);
OSQPInt        osqp_update_eps_rel(OSQPSolver* osqpSolver, OSQPFloat eps_rel_new);
OSQPInt        osqp_update_eps_prim_inf(OSQPSolver* osqpSolver, OSQPFloat eps_prim_inf_new);
OSQPInt        osqp_update_eps_dual_inf(OSQPSolver* osqpSolver, OSQPFloat eps_dual_inf_new);
OSQPInt        osqp_update_alpha(OSQPSolver* osqpSolver, OSQPFloat alpha_new);
OSQPInt        osqp_update_delta(OSQPSolver* osqpSolver, OSQPFloat delta_new);
OSQPInt        osqp_update_polish_refine_iter(OSQPSolver* osqpSolver, OSQPInt polish_refine_iter_new);
OSQPInt        osqp_update_verbose(OSQPSolver* osqpSolver, OSQPInt verbose_new);
OSQPInt        osqp_update_scaled_termination(OSQPSolver* osqpSolver, OSQPInt scaled_termination_new);
OSQPInt        osqp_update_check_termination(OSQPSolver* osqpSolver, OSQPInt check_termination_new);
OSQPInt        osqp_update_warm_start(OSQPSolver* osqpSolver, OSQPInt warm_start_new);
#ifdef PROFILING
OSQPInt        osqp_update_time_limit(OSQPSolver* osqpSolver, OSQPFloat time_limit_new);
#endif /* ifdef PROFILING */


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
    //std::cout << "cmd = " << cmd << std::endl; //DELETE HERE
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
        //OsqpData* osqpData;  // OSQP data identifier //DELETE HERE
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
        //clean up the handle object
        //if (prhs[1]) destroyObject<OSQPSolver>(prhs[1]);
        //osqp_cleanup(osqpSolver);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }

    // report the current settings
    if (!strcmp("current_settings", cmd)) {
      //throw an error if this is called before solver is configured
      if(!osqpData->solver){
        mexErrMsgTxt("Solver is uninitialized.  No settings have been configured.");
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
        OSQPInt * Pp = (OSQPInt*)copyToCintVector(mxGetJc(P), dataN + 1);
        OSQPInt * Pi = (OSQPInt*)copyToCintVector(mxGetIr(P), Pp[dataN]);
        OSQPFloat * Px = copyToOSQPFloatVector(mxGetPr(P), Pp[dataN]);
        OSQPCscMatrix* dataP = new OSQPCscMatrix;
        csc_set_data(dataP, dataN, dataN, Pp[dataN], Px, Pi, Pp);

        // Matrix A: nnz = A->p[n]
        OSQPInt* Ap = (OSQPInt*)copyToCintVector(mxGetJc(A), dataN + 1);
        OSQPInt* Ai = (OSQPInt*)copyToCintVector(mxGetIr(A), Ap[dataN]);
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
        exitflag = osqp_setup(&osqpData->solver, dataP, dataQ, dataA, dataL, dataU, dataM, dataN, settings);

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

        plhs[0] = mxCreateDoubleScalar(osqpData->solver->work->data->n);
        plhs[1] = mxCreateDoubleScalar(osqpData->solver->work->data->m);

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
        if(!mxIsEmpty(q)){
            q_vec = copyToOSQPFloatVector(mxGetPr(q),
                                       osqpData->solver->work->data->n);
        }
        if(!mxIsEmpty(l)){
            l_vec = copyToOSQPFloatVector(mxGetPr(l),
                                       osqpData->solver->work->data->m);
        }
        if(!mxIsEmpty(u)){
            u_vec = copyToOSQPFloatVector(mxGetPr(u),
                                       osqpData->solver->work->data->m);
        }
        if(!mxIsEmpty(Px)){
            Px_vec = copyToOSQPFloatVector(mxGetPr(Px), Px_n);
        }
        if(!mxIsEmpty(Ax)){
            Ax_vec = copyToOSQPFloatVector(mxGetPr(Ax), Ax_n);
        }
        if(!mxIsEmpty(Px_idx)){
            Px_idx_vec = copyDoubleToCintVector(mxGetPr(Px_idx), Px_n);
        }
        if(!mxIsEmpty(Ax_idx)){
            Ax_idx_vec = copyDoubleToCintVector(mxGetPr(Ax_idx), Ax_n);
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

      if(!mxIsEmpty(x)){
          x_vec = copyToOSQPFloatVector(mxGetPr(x),
                                      osqpData->solver->work->data->n);
      }
      if(!mxIsEmpty(y)){
          y_vec = copyToOSQPFloatVector(mxGetPr(y),
                                      osqpData->solver->work->data->m);
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


        // Allocate space for solution
        // primal variables
        plhs[0] = mxCreateDoubleMatrix(osqpData->solver->work->data->n,1,mxREAL);
        // dual variables
        plhs[1] = mxCreateDoubleMatrix(osqpData->solver->work->data->m,1,mxREAL);
        // primal infeasibility certificate
        plhs[2] = mxCreateDoubleMatrix(osqpData->solver->work->data->m,1,mxREAL);
        // dual infeasibility certificate
        plhs[3] = mxCreateDoubleMatrix(osqpData->solver->work->data->n,1,mxREAL);

        //copy results to mxArray outputs
        //assume that five outputs will always
        //be returned to matlab-side class wrapper
        if ((osqpData->solver->info->status_val != OSQP_PRIMAL_INFEASIBLE) &&
            (osqpData->solver->info->status_val != OSQP_DUAL_INFEASIBLE)){

            //primal and dual solutions
            castToDoubleArr(osqpData->solver->solution->x, mxGetPr(plhs[0]), osqpData->solver->work->data->n);
            castToDoubleArr(osqpData->solver->solution->y, mxGetPr(plhs[1]), osqpData->solver->work->data->m);

            //infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), osqpData->solver->work->data->m);
            setToNaN(mxGetPr(plhs[3]), osqpData->solver->work->data->n);

        } else if (osqpData->solver->info->status_val == OSQP_PRIMAL_INFEASIBLE ||
        osqpData->solver->info->status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE){ //primal infeasible

            //primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), osqpData->solver->work->data->n);
            setToNaN(mxGetPr(plhs[1]), osqpData->solver->work->data->m);

            //primal infeasibility certificates
            castToDoubleArr(osqpData->solver->solution->prim_inf_cert, mxGetPr(plhs[2]), osqpData->solver->work->data->m);

            //dual infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[3]), osqpData->solver->work->data->n);

            // Set objective value to infinity
            osqpData->solver->info->obj_val = mxGetInf();

        } else { //dual infeasible

            //primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), osqpData->solver->work->data->n);
            setToNaN(mxGetPr(plhs[1]), osqpData->solver->work->data->m);

            //primal infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), osqpData->solver->work->data->m);

            //dual infeasibility certificates
            castToDoubleArr(osqpData->solver->solution->dual_inf_cert, mxGetPr(plhs[3]), osqpData->solver->work->data->n);

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
  OSQPFloat* out = (OSQPFloat*)c_malloc(numel * sizeof(OSQPFloat));

  //copy data
  for(OSQPInt i=0; i < numel; i++){
      out[i] = (OSQPFloat)vecData[i];
  }
  return out;
}

//Dynamically creates a OSQPInt vector copy of the input.
OSQPInt* copyToCintVector(mwIndex* vecData, OSQPInt numel){
  // This memory needs to be freed!
  OSQPInt* out = (OSQPInt*)c_malloc(numel * sizeof(OSQPInt));

  //copy data
  for(OSQPInt i=0; i < numel; i++){
      out[i] = (OSQPInt)vecData[i];
  }
  return out;

}

//Dynamically copies a double vector to OSQPInt.
OSQPInt* copyDoubleToCintVector(double* vecData, OSQPInt numel){
  // This memory needs to be freed!
  OSQPInt* out = (OSQPInt*)c_malloc(numel * sizeof(OSQPInt));

  //copy data
  for(OSQPInt i=0; i < numel; i++){
      out[i] = (OSQPInt)vecData[i];
  }
  return out;

}

void castCintToDoubleArr(OSQPInt *arr, double* arr_out, OSQPInt len) {
    for (OSQPInt i = 0; i < len; i++) {
        arr_out[i] = (double)arr[i];
    }
}

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

  #ifdef PROFILING
  //if not profiling, these fields will be empty
  mxSetField(mxPtr, 0, "setup_time",  mxCreateDoubleScalar(info->setup_time));
  mxSetField(mxPtr, 0, "solve_time",  mxCreateDoubleScalar(info->solve_time));
  mxSetField(mxPtr, 0, "update_time", mxCreateDoubleScalar(info->update_time));
  mxSetField(mxPtr, 0, "polish_time", mxCreateDoubleScalar(info->polish_time));
  mxSetField(mxPtr, 0, "run_time",    mxCreateDoubleScalar(info->run_time));
  #endif /* ifdef PROFILING */

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
  mxSetField(mxPtr, 0, "time_limit",             mxCreateDoubleScalar(settings->cg_precond));
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

  //This does basically the same job as copyMxStructToSettings which was used
  //during setup, but uses the provided update functions in osqp.h to update
  //settings in the osqp workspace.  Protects against bad parameter writes
  //or future modifications to updated settings handling

  osqp_update_max_iter(osqpSolver,
    (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "max_iter")));
  osqp_update_eps_abs(osqpSolver,
    (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs")));
  osqp_update_eps_rel(osqpSolver,
    (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_rel")));
  osqp_update_eps_prim_inf(osqpSolver,
    (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_prim_inf")));
  osqp_update_eps_dual_inf(osqpSolver,
    (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "eps_dual_inf")));
  osqp_update_alpha(osqpSolver,
    (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "alpha")));
  osqp_update_delta(osqpSolver,
    (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "delta")));
  osqp_update_polish_refine_iter(osqpSolver,
    (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "polish_refine_iter")));
  osqp_update_verbose(osqpSolver,
    (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "verbose")));
  osqp_update_scaled_termination(osqpSolver,
    (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "scaled_termination")));
  osqp_update_check_termination(osqpSolver,
    (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "check_termination")));
  osqp_update_warm_start(osqpSolver,
    (OSQPInt)mxGetScalar(mxGetField(mxPtr, 0, "warm_start")));
#ifdef PROFILING    
  osqp_update_time_limit(osqpSolver,
    (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "time_limit")));
#endif /* ifdef PROFILING */

  // Check for settings that need special update
  // Update them only if they are different than already set values
  OSQPFloat rho_new = (OSQPFloat)mxGetScalar(mxGetField(mxPtr, 0, "rho"));
  // Check if it has changed
  if (c_absval(rho_new - osqpSolver->settings->rho) > NEW_SETTINGS_TOL){
      exitflag = osqp_update_rho(osqpSolver, rho_new);
      if (exitflag){
          mexErrMsgTxt("rho update error!");
      }
  }


}

/****************************
* Update problem settings  *
****************************/
OSQPInt osqp_update_max_iter(OSQPSolver* osqpSolver, OSQPInt max_iter_new) {

  // Check if osqpSolver has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that max_iter is positive
  if (max_iter_new <= 0) {
#ifdef PRINTING
    c_eprint("max_iter must be positive");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update max_iter
  osqpSolver->settings->max_iter = max_iter_new;

  return 0;
}

OSQPInt osqp_update_eps_abs(OSQPSolver* osqpSolver, OSQPFloat eps_abs_new) {

  // Check if osqpSolver has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that eps_abs is positive
  if (eps_abs_new < 0.) {
#ifdef PRINTING
    c_eprint("eps_abs must be nonnegative");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update eps_abs
  osqpSolver->settings->eps_abs = eps_abs_new;

  return 0;
}

OSQPInt osqp_update_eps_rel(OSQPSolver* osqpSolver, OSQPFloat eps_rel_new) {

  // Check if osqpSolver has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that eps_rel is positive
  if (eps_rel_new < 0.) {
#ifdef PRINTING
    c_eprint("eps_rel must be nonnegative");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update eps_rel
  osqpSolver->settings->eps_rel = eps_rel_new;

  return 0;
}

OSQPInt osqp_update_eps_prim_inf(OSQPSolver* osqpSolver, OSQPFloat eps_prim_inf_new) {

  // Check if osqpSolver has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that eps_prim_inf is positive
  if (eps_prim_inf_new < 0.) {
#ifdef PRINTING
    c_eprint("eps_prim_inf must be nonnegative");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update eps_prim_inf
  osqpSolver->settings->eps_prim_inf = eps_prim_inf_new;

  return 0;
}

OSQPInt osqp_update_eps_dual_inf(OSQPSolver* osqpSolver, OSQPFloat eps_dual_inf_new) {

  // Check if osqpSolver has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that eps_dual_inf is positive
  if (eps_dual_inf_new < 0.) {
#ifdef PRINTING
    c_eprint("eps_dual_inf must be nonnegative");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update eps_dual_inf
  osqpSolver->settings->eps_dual_inf = eps_dual_inf_new;


  return 0;
}

OSQPInt osqp_update_alpha(OSQPSolver* osqpSolver, OSQPFloat alpha_new) {

  // Check if osqpSolver has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that alpha is between 0 and 2
  if ((alpha_new <= 0.) || (alpha_new >= 2.)) {
#ifdef PRINTING
    c_eprint("alpha must be between 0 and 2");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update alpha
  osqpSolver->settings->alpha = alpha_new;

  return 0;
}

OSQPInt osqp_update_warm_start(OSQPSolver* osqpSolver, OSQPInt warm_start_new) {

  // Check if osqpSolver has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that warm_start is either 0 or 1
  if ((warm_start_new != 0) && (warm_start_new != 1)) {
#ifdef PRINTING
    c_eprint("warm_start should be either 0 or 1");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update warm_start
  osqpSolver->settings->warm_starting = warm_start_new;

  return 0;
}

OSQPInt osqp_update_delta(OSQPSolver* osqpSolver, OSQPFloat delta_new) {

  // Check if workspace has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that delta is positive
  if (delta_new <= 0.) {
# ifdef PRINTING
    c_eprint("delta must be positive");
# endif /* ifdef PRINTING */
    return 1;
  }

  // Update delta
  osqpSolver->settings->delta = delta_new;

  return 0;
}

OSQPInt osqp_update_polish_refine_iter(OSQPSolver* osqpSolver, OSQPInt polish_refine_iter_new) {

  // Check if workspace has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that polish_refine_iter is nonnegative
  if (polish_refine_iter_new < 0) {
# ifdef PRINTING
    c_eprint("polish_refine_iter must be nonnegative");
# endif /* ifdef PRINTING */
    return 1;
  }

  // Update polish_refine_iter
  osqpSolver->settings->polish_refine_iter = polish_refine_iter_new;

  return 0;
}

OSQPInt osqp_update_verbose(OSQPSolver* osqpSolver, OSQPInt verbose_new){

  // Check if workspace has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that verbose is either 0 or 1
  if ((verbose_new != 0) && (verbose_new != 1)) {
# ifdef PRINTING
    c_eprint("verbose should be either 0 or 1");
# endif /* ifdef PRINTING */
    return 1;
  }

  // Update verbose
  osqpSolver->settings->verbose = verbose_new;

  return 0;
}

OSQPInt osqp_update_scaled_termination(OSQPSolver* osqpSolver, OSQPInt scaled_termination_new) {

  // Check if osqpSolver has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that scaled_termination is either 0 or 1
  if ((scaled_termination_new != 0) && (scaled_termination_new != 1)) {
#ifdef PRINTING
    c_eprint("scaled_termination should be either 0 or 1");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update scaled_termination
  osqpSolver->settings->scaled_termination = scaled_termination_new;

  return 0;
}

OSQPInt osqp_update_check_termination(OSQPSolver* osqpSolver, OSQPInt check_termination_new) {

  // Check if osqpSolver has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that check_termination is nonnegative
  if (check_termination_new < 0) {
#ifdef PRINTING
    c_eprint("check_termination should be nonnegative");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update check_termination
  osqpSolver->settings->check_termination = check_termination_new;

  return 0;
}

#ifdef PROFILING

OSQPInt osqp_update_time_limit(OSQPSolver* osqpSolver, OSQPFloat time_limit_new) {

  // Check if osqpSolver has been initialized
  if (!osqpSolver) return osqp_error(OSQP_WORKSPACE_NOT_INIT_ERROR);

  // Check that time_limit is nonnegative
  if (time_limit_new < 0.) {
# ifdef PRINTING
    c_print("time_limit must be nonnegative\n");
# endif /* ifdef PRINTING */
    return 1;
  }

  // Update time_limit
  osqpSolver->settings->time_limit = time_limit_new;

  return 0;
}
#endif /* ifdef PROFILING */

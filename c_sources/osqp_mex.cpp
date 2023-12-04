#include <map>

#include "mex.h"
#include "matrix.h"
#include "osqp.h"

// Mex-specific functionality
#include "osqp_mex.hpp"
#include "osqp_struct.h"
#include "arrays_matlab.h"
#include "memory_matlab.h"

//TODO: Check if this definition is required, and maybe replace it with:
//   enum linsys_solver_type { QDLDL_SOLVER, MKL_PARDISO_SOLVER };
#define QDLDL_SOLVER 0 //Based on the previous API

// Wrapper class to pass the OSQP solver back and forth with Matlab
class OsqpData
{
public:
  OsqpData() :
    solver(NULL)
    {}
  OSQPSolver* solver;
};


// Internal utility function
static void setToNaN(double* arr_out, OSQPInt len){
    OSQPInt i;
    for (i = 0; i < len; i++) {
        arr_out[i] = mxGetNaN();
    }
}


// Main mex function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // OSQP solver wrapper
    OsqpData* osqpData;

    // Exitflag
    OSQPInt exitflag = 0;

    // Get the command string
    char cmd[64];

    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    /*
     * First check to see if a new object was requested
     */
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1){
            mexErrMsgTxt("New: One output expected.");
          }
        // Return a handle to a new C++ wrapper instance
        osqpData = new OsqpData;
        plhs[0] = convertPtr2Mat<OsqpData>(osqpData);
        return;
    }

    /*
     * Next check to see if any of the static methods were called
     */
    // Report the version
    if (!strcmp("version", cmd)) {
        plhs[0] = mxCreateString(osqp_version());
        return;
    }

    // Report the default settings
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

    // Return solver constants
    if (!strcmp("constant", cmd)) {
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
        int  constantLength = mxGetN(prhs[1]) + 1;
        mxGetString(prhs[1], constant, constantLength);

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

    /*
     * Finally, check to see if this is a function operating on a solver instance
     */

    // Check for a second input, which should be the class instance handle
    if (nrhs < 2)
        mexErrMsgTxt("Second input should be a class instance handle.");


    // Get the class instance pointer from the second input
    osqpData = convertMat2Ptr<OsqpData>(prhs[1]);

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
        osqp_update_settings(osqpData->solver, settings.GetOSQPStruct());
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
        OSQPFloat* dataQ   = cloneVector<OSQPFloat>(mxGetPr(q), dataN);
        OSQPFloat* dataL   = cloneVector<OSQPFloat>(mxGetPr(l), dataM);
        OSQPFloat* dataU   = cloneVector<OSQPFloat>(mxGetPr(u), dataM);

        // Matrix P:  nnz = P->p[n]
        OSQPInt * Pp = cloneVector<OSQPInt>(mxGetJc(P), dataN + 1);
        OSQPInt * Pi = cloneVector<OSQPInt>(mxGetIr(P), Pp[dataN]);
        OSQPFloat * Px = cloneVector<OSQPFloat>(mxGetPr(P), Pp[dataN]);
        OSQPCscMatrix* dataP = (OSQPCscMatrix*)c_calloc(1,sizeof(OSQPCscMatrix));
        csc_set_data(dataP, dataN, dataN, Pp[dataN], Px, Pi, Pp);

        // Matrix A: nnz = A->p[n]
        OSQPInt* Ap = cloneVector<OSQPInt>(mxGetJc(A), dataN + 1);
        OSQPInt* Ai = cloneVector<OSQPInt>(mxGetIr(A), Ap[dataN]);
        OSQPFloat * Ax = cloneVector<OSQPFloat>(mxGetPr(A), Ap[dataN]);
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
        exitflag = osqp_setup(&(osqpData->solver), dataP, dataQ, dataA, dataL, dataU, dataM, dataN, settings.GetOSQPStruct());
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
            q_vec = cloneVector<OSQPFloat>(mxGetPr(q), n);
        }
        if(!mxIsEmpty(l)){
            l_vec = cloneVector<OSQPFloat>(mxGetPr(l), m);
        }
        if(!mxIsEmpty(u)){
            u_vec = cloneVector<OSQPFloat>(mxGetPr(u), m);
        }
        if(!mxIsEmpty(Px)){
            Px_vec = cloneVector<OSQPFloat>(mxGetPr(Px), Px_n);
        }
        if(!mxIsEmpty(Ax)){
            Ax_vec = cloneVector<OSQPFloat>(mxGetPr(Ax), Ax_n);
        }
        if(!mxIsEmpty(Px_idx)){
            Px_idx_vec = cloneVector<OSQPInt>(mxGetPr(Px_idx), Px_n);
        }
        if(!mxIsEmpty(Ax_idx)){
            Ax_idx_vec = cloneVector<OSQPInt>(mxGetPr(Ax_idx), Ax_n);
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

    if (!strcmp("warm_start", cmd)) {

      //throw an error if this is called before solver is configured
      if(!osqpData->solver){
          mexErrMsgTxt("Solver has not been initialized.");
      }

      // Fill x and y
      const mxArray *x = prhs[2];
      const mxArray *y = prhs[3];

      // Copy vectors to ensure they are cast as OSQPFloat
      OSQPFloat *x_vec = NULL;
      OSQPFloat *y_vec = NULL;

      OSQPInt n, m;
      osqp_get_dimensions(osqpData->solver, &m, &n);

      if(!mxIsEmpty(x)){
          x_vec = cloneVector<OSQPFloat>(mxGetPr(x),n);
      }
      if(!mxIsEmpty(y)){
          y_vec = cloneVector<OSQPFloat>(mxGetPr(y),m);
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
            copyVector<double>(mxGetPr(plhs[0]), osqpData->solver->solution->x, n);
            copyVector<double>(mxGetPr(plhs[1]), osqpData->solver->solution->y, m);

            //infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), m);
            setToNaN(mxGetPr(plhs[3]), n);

        } else if (osqpData->solver->info->status_val == OSQP_PRIMAL_INFEASIBLE ||
        osqpData->solver->info->status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE){ //primal infeasible

            //primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), n);
            setToNaN(mxGetPr(plhs[1]), m);

            //primal infeasibility certificates
            copyVector<double>(mxGetPr(plhs[2]), osqpData->solver->solution->prim_inf_cert, m);

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
            copyVector<double>(mxGetPr(plhs[3]), osqpData->solver->solution->dual_inf_cert, n);

            // Set objective value to -infinity
            osqpData->solver->info->obj_val = -mxGetInf();
        }

        if (osqpData->solver->info->status_val == OSQP_NON_CVX) {
            osqpData->solver->info->obj_val = mxGetNaN();
        }

        // Populate the info structure
        OSQPInfoWrapper info(osqpData->solver->info);
        plhs[4] = info.GetMxStruct();

        return;
    }

    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}

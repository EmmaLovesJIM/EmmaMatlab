/*  File    : KF.c
 *  Abstract:
 *
 *      Implementation of Kalman Filter for parameter estimation,
        in C language for
 *      Matlab S-function, for the SISO case
 */

/* Preprozessor:
 *
 *  #define <name> <substitution>
 *  substitutes <name> by <substitution> everywhere in the following C-code
 *
 *  #define <name>(parameterlist) <substitution>
 *  possible to use parameters. Eg: #define A(i,k) a[i + k*n] allows to 
 *  access a matrix (which here is implemented as an one dimensional array 
 *  of size m*n) using the Matlab syntax.
 */
#define S_FUNCTION_NAME KF
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#define N_U 2
#define N_Y 2
#define TS 1

/*#undefine VERBOSE /* Define 'VERBOSE' when debugging */

#define U(element) (*uPtrs[element])  /* Pointer to Input Port0 */
#define PARAM(element) (*LambdaPtrs[element])  /* Pointer to Input Port0 */
#define THETA(element) (x[element])  /* Pointer to Theta Port0 */
#define P(element) (rwork[element])  /* Pointer to Theta Port0 */
#define MTEMP(element) (rwork[element + N_THETA*N_THETA])  /* Pointer to Theta Port0 */
#define X1(element) (rwork[element + 2*N_THETA*N_THETA])  /* Pointer to Theta Port0 */
#define VTEMP(element) (rwork[element + N_THETA*(2*N_THETA+1)])  /* Pointer to Theta Port0 */
#define VTEMP1(element) (rwork[element + N_THETA*(2*N_THETA+2)])  /* Pointer to Theta Port0 */
#define TEMP (rwork[N_THETA*(2*N_THETA+3)])  /* Pointer to Theta Port0 */
#define TEMP1 (rwork[N_THETA*(2*N_THETA+3)+1])  /* Pointer to Theta Port0 */

static int_T N_THETA=N_Y+N_U;

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }
    
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, N_THETA);
    
    if (!ssSetNumInputPorts(S, 2)) return;
    ssSetInputPortWidth(S, 0, 2);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    
    // This is for Parameters
    ssSetInputPortWidth(S, 1, 5);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    
    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, N_THETA);
    
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, N_THETA*(2*N_THETA+3)+1+5); //X1,P,TEMP
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);
    
    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    //ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy that we inherit our sample time from the driving block.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, TS);   
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
}

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    Initialize both discrete states to one.
 */
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *x    = ssGetRealDiscStates(S);
    InputRealPtrsType LambdaPtrs = ssGetInputPortRealSignalPtrs(S,1);
    real_T *rwork = ssGetRWork(S);
    int  j;
    
    /* initialise with zeros */
    for (j = 0; j < N_THETA*N_THETA; j++) {
        P(j)=0;
    }
    
    for (j = 0; j < N_THETA; j++) {
        /* Initialise covariance matrix with large term in diagonal */
        P(j + N_THETA*j )=.0001;
        X1(j)=0;
        //THETA(j)=0;
    }
    THETA(0) = -.4;
    THETA(1) = -.5;
    THETA(2) = 0.002;
    THETA(3) = 0.001;
}


/* Function: mdlOutputs =======================================================
 * Abstract:
 *      y = Theta
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T            *y     = ssGetOutputPortRealSignal(S,0);
    real_T            *x     = ssGetRealDiscStates(S);
    int i;
    
    for (i = 0; i < N_THETA; i++) {
        y[i] = THETA(i);
    }
    
    
}

#define MDL_UPDATE
/* Function: mdlUpdate ======================================================
 * Abstract:
 *      P = (P - P*x*(inv(1 + x'*P*x))*x'*P)/lambda;
 *      thetaRLS = thetaRLS + P*x*(y(t) - x'*thetaRLS);
 */
static void mdlUpdate(SimStruct *S, int_T tid)
{
    InputRealPtrsType uPtrs      = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType LambdaPtrs = ssGetInputPortRealSignalPtrs(S,1);
    real_T            *x         = ssGetRealDiscStates(S);
    real_T            *rwork     = ssGetRWork(S);
    int               i;
    double rv;
    rv = PARAM(0); /* Noise variance*/
    
    // Phi = Phi +Rw
    for (i = 0; i < N_THETA; i++) {
        P(i+N_THETA*i) = P(i+N_THETA*i) + PARAM(1+i);
    }
   
    /* 1 Multiply P*x -> VTEMP */
    matmul(&P(0), N_THETA, N_THETA, &X1(0),1,&VTEMP(0));
    
    /* 2 x'VTEMP  -> TEMP*/
    matmul(&X1(0), 1, N_THETA, &VTEMP(0), 1, &TEMP);
    
    /* 3 (rv + x'*Phi*x)^(-1) */   
    TEMP = 1/(rv + TEMP);
    
    /* 4 P*x * (rv + x'*P*x)^(-1) */
    for(i = 0; i < N_THETA; i++){
        VTEMP1(i) = VTEMP(i)*TEMP;
    }
    
    // 5 x'*theta
    matmul(&X1(0), 1, N_THETA,&THETA(0),1,&TEMP1);
   
    // 6 y - x'theta 
    TEMP1 = U(0) - TEMP1; 
    
    // 7 phi*(y - x'*theta)
    for(i = 0; i < N_THETA; i++){
        VTEMP(i) = VTEMP1(i)*TEMP1;
    }
    
    /* 8  ThetaHat = ThetaHat(t-1) + phi*(y - x'*ThetaHat) */
    for (i = 0; i < N_THETA; i++) {
        THETA(i) = THETA(i) + VTEMP(i);
    }
    
    /* 9 x'*P  -> VTEMP*/
    matmul(&X1(0), 1, N_THETA, &P(0), N_THETA, &VTEMP(0));
    
    // 10 phi*x'*P -> MTEMP
    matmul(&VTEMP1(0), N_THETA, 1, &VTEMP(0), N_THETA, &MTEMP(0));
    
    // 11 P = (I - phi*x') P
    for (i = 0; i < N_THETA*N_THETA; i++) {
        P(i) = P(i)- MTEMP(i);
    }
    
    /* Update the observation vector
     */
    for (i = N_Y-1; i >0; i--) {
        X1(i)=X1(i-1);
    }
    X1(0)=-U(0); //ouput
    
    for (i = N_U-1; i >0; i--) {
        X1(N_Y+i)=X1(N_Y+i-1);
    }
    X1(2)=U(1); //input
}

/* Function:    Matrix multiplication C=A*B, (A,B,C stored columnwise) 
 *              A: (m x n) matrix
 *              B: (n x l) matrix
 *              C: (m x l) matrix 
 */

#define A(i,k) a[i + k*n]
#define B(k,j) b[k + j*m]
#define C(i,j) c[i + j*n]

void matmul(a,n,m,b,l,c)
double a[],b[],c[];
int n,m,l;
{
    int i,j,k; double s;
    for( i=0 ; i < n; i++)
    {
        for( j=0; j < l; j++)
        {
            s = 0.;
            for( k=0; k< m; k++)
            {
                s += A(i,k)*B(k,j);
            }
            C(i,j) = s;
        }
    }
}

/* Function: Display of Matrix A (A stored columnwise) */

void matdisp(a,n,m)
double a[];
int n,m;
{
    int i,j;
    for( i=0 ; i < n; i++)
    {
        if (i==0){
            printf("=[%f",A(0,0));
        } else {
            printf("   %f",A(i,0));
        }
        for( j=1; j < m; j++)
        {
            printf(" %f",A(i,j));
        }
        if (i==n-1){
            printf("];\n");
        }else{
            printf("\n");
        }
    }
}

/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
    //UNUSED_ARG(S); /* unused input argument */
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

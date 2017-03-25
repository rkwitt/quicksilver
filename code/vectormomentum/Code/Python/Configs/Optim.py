"""Generalized optimization configs and routines"""

import Config

def _valspec(_):
    """Hook to do extra validation"""
    pass

OptimLongitudinalConfigSpec = {
    'Niter':
    Config.Param(default=100,
                 comment="Number gradient descent iterations to run"),
    
    'stepSizeGroup':
    Config.Param(default=0.01,
                 comment="Gradient descent step size for initial momenta for group geodesic"),

    'stepSizeResidual':
    Config.Param(default=0.01,
                 comment="Gradient descent step size for initial momenta for residual geodesic"),

    'maxPert':
    Config.Param(default=None,
                 comment="Maximum perturbation. Rough method of automatic " +
                 "step-size selection."),
    'method':
    Config.Param(default="FIXEDGD",
                 comment="Optimization scheme.  Only FIXEDGD is supported"),
    'NIterForInverse':
    Config.Param(default=5,
                 comment="Iterations for computing fixed point iterative inverse of a diffeomorphism."),
    'integMethodGroup':
    Config.Param(default="EULER",
                 comment="Integration scheme for group geodesic.  EULER or RK4"),

    'nTimeStepsGroup':
    Config.Param(default=10,
                 comment="Number of time discretization steps for group geodesic"),

    'integMethodResidual':
    Config.Param(default="EULER",
                 comment="Integration scheme for residual geodesics.  EULER or RK4"),

    'nTimeStepsResidual':
    Config.Param(default=10,
                 comment="Number of time discretization steps for residual geodesic"),

    '_resource': "OptimLong",
    '_validation_hook': _valspec}

OptimConfigSpec = {
    'Niter':
    Config.Param(default=100,
                 comment="Number gradient descent iterations to run"),
    
    'stepSize':
    Config.Param(default=0.01,
                 comment="Gradient descent step size for initial momenta"),
    'maxPert':
    Config.Param(default=None,
                 comment="Maximum perturbation. Rough method of automatic " +
                 "step-size selection. Currently is used for autoreduce fixed gradient descent step size if energy increases."),
    'method':
    Config.Param(default="FIXEDGD",
                 comment="Optimization scheme.  Only FIXEDGD is supported"),
    'NIterForInverse':
    Config.Param(default=5,
                 comment="Iterations for computing fixed point iterative inverse of a diffeomorphism."),
    'integMethod':
    Config.Param(default="EULER",
                 comment="Integration scheme.  EULER or RK4"),

    'nTimeSteps':
    Config.Param(default=10,
                 comment="Number of time discretization steps for geodesic"),
    '_resource': "Optim",
    '_validation_hook': _valspec}

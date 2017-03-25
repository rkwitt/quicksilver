"""Core components of Longitudinal algorithms based on vector momentum"""
import Config
import PyCA.Core as ca

VMConfigSpec = {
    'sigma':
    Config.Param(default=0.5,
                 required=True,
                 comment="Regularization weight on image match term, lower values assign higher weight to matching (sigma)"),
    'diffOpParams':
    Config.Param(default=[0.01, 0.01, 0.001],
                 required=True,
                 comment="Differential operator parameters: alpha, beta and gamma")
    }


VMLongitudinalConfigSpec = {
    'varIntercept':
    Config.Param(default=0.1,
                 required=True,
                 comment="Variance for intercept term for group geodesic (sigma_I)"),
    'varSlope':
    Config.Param(default=0.1,
                 required=True,
                 comment="Variance for slope term for group geodesic (sigma_S)"),
    'varInterceptReg':
    Config.Param(default=0.1,
                 required=True,
                 comment="Regularization weight for intercept matching (sigma_i)."),
    'diffOpParamsGroup':
    Config.Param(default=[0.01, 0.01, 0.001],
                 required=True,
                 comment="Differential operator parameters for Group: alpha, beta and gamma"),
    'diffOpParamsResidual':
    Config.Param(default=[0.01, 0.01, 0.001],
                 required=True,
                 comment="Differential operator parameters for Residual: alpha, beta and gamma")
    }

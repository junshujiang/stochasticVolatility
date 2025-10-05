// file: cloglike-stochvol.c
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "cgeneric.h"

#define Calloc(n_, type_) (type_ *)calloc((n_), sizeof(type_))
#define LOG_NORMC_GAUSSIAN (-0.91893853320467274178032973640560) /* -0.5*log(2π) */

/*
 * 模型： y ~ N(0, var = exp(x) + 1/tau)
 * theta[0] = log(tau)  (tau > 0)
 * ny 至少 1 列：y[0] = 观测值
 * nx 是需要评估的线性预测数组 x[0..nx-1]
 */
double *inla_cloglike_stochvol(inla_cloglike_cmd_tp cmd, double *theta,
                               inla_cgeneric_data_tp *data,
                               int ny, double *y, int nx, double *x, double *result)
{
    double *ret = NULL;

    /* 读取超参数 tau */
    double ltau = NAN, tau = NAN;
    if (theta) { ltau = theta[0]; tau = exp(ltau); }

    switch (cmd) {
    case INLA_CLOGLIKE_INITIAL: {
        /* 一个超参数：log(tau)，给个优化友好的初值 4 */
        ret = Calloc(2, double);
        ret[0] = 1;      /* #thetas */
        ret[1] = 4.0;    /* init of log(tau) */
    } break;

    case INLA_CLOGLIKE_LOG_PRIOR: {
        /* 先验：tau ~ Gamma(1,1)，对数变量需要加 Jacobian -> -tau + log(tau) = -exp(ltau) + ltau */
        ret = Calloc(1, double);
        ret[0] = -tau + ltau;
    } break;

    case INLA_CLOGLIKE_LOGLIKE: {
        assert(ny >= 1);
        const double y0 = y[0];
        const double var_offset = (isinf(tau) || isnan(tau)) ? 0.0 : 1.0 / tau;

        for (int i = 0; i < nx; i++) {
            /* 原函数里的 PREDICTOR_INVERSE_LINK(x[i], off) 这里采用最常见的 exp 链接 */
            const double var = exp(x[i]) + var_offset;
            result[i] = LOG_NORMC_GAUSSIAN - 0.5*log(var) - 0.5*(y0*y0)/var;
        }
    } break;

    case INLA_CLOGLIKE_CDF: {
        /* INLA 的 CDF: Prob(Ỹ < y | x,theta)；零均值正态 => Φ(y/√var) */
        assert(ny >= 1);
        const double y0 = y[0];
        const double var_offset = (isinf(tau) || isnan(tau)) ? 0.0 : 1.0 / tau;

        for (int i = 0; i < nx; i++) {
            const double var = exp(x[i]) + var_offset;
            const double z = y0 / sqrt(var);
            /* Φ(z) = 0.5 * (1 + erf(z / sqrt(2))) */
            result[i] = 0.5 * (1.0 + erf(z * M_SQRT1_2));
        }
    } break;

    case INLA_CLOGLIKE_QUIT:
    default:
        break;
    }

    return ret;
}

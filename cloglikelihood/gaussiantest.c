#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include "cgeneric.h"
#include <errno.h>
#include <string.h>
#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define Malloc(n_, type_)  (type_ *)malloc((n_) * sizeof(type_))
#define SQR(x_) ((x_)*(x_))




// 放在文件顶部（包含 enum 定义之后）
static const char* cloglike_cmd_name(inla_cloglike_cmd_tp cmd) {
    switch (cmd) {
    case INLA_CLOGLIKE_INITIAL:   return "INITIAL";
    case INLA_CLOGLIKE_LOG_PRIOR: return "LOG_PRIOR";
    case INLA_CLOGLIKE_LOGLIKE:   return "LOGLIKE";
    case INLA_CLOGLIKE_CDF:       return "CDF";
    case INLA_CLOGLIKE_QUIT:      return "QUIT";
    default:                      return "UNKNOWN";
    }
}


static int write_string_to_file(const char *path, const char *s)
{
    FILE *fp = fopen(path, "a");
    if (!fp) { fprintf(stderr, "open %s failed: %s\n", path, strerror(errno)); return -1; }
    fputs(s, fp);
    if (s[0] == '\0' || s[strlen(s)-1] != '\n') fputc('\n', fp);
    fclose(fp);
    return 0;
}

double *inla_cloglike_gaussian(inla_cloglike_cmd_tp cmd, double *theta,
			       inla_cgeneric_data_tp *data, int ny, double *y, int nx, double *x, double *result)
{
#define LOG_NORMC_GAUSSIAN (-0.91893853320467274178032973640560)	/* -1/2 * log(2*pi) */
#define CDF(x_) (0.5 * (1.0 + erf(M_SQRT1_2 * (x_))))


	double *ret = NULL, prec, lprec;

	if (theta) {
		lprec = theta[0];

		prec = exp(lprec);
	} else {
		prec = lprec = NAN;
	}

	switch (cmd) {
	case INLA_CLOGLIKE_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters

		ret = Malloc(2, double);
		ret[0] = 1;
		ret[1] = 4.0;
		char buf[64];
		snprintf(buf, sizeof buf, "init %.17g", ret[1]);
		write_string_to_file("log.txt", buf);


	}
		break;

	case INLA_CLOGLIKE_LOG_PRIOR:
	{
		// return c(LOG_PRIOR). with a Gamma(1,1) for precision, this is the log prior for the log(precision).
		ret = Malloc(1, double);
		ret[0] = -prec + lprec;

		char buf[64];
		snprintf(buf, sizeof buf, "theta prior %.17g", ret[0]);
		write_string_to_file("log.txt", buf);
	}
		break;

	case INLA_CLOGLIKE_LOGLIKE:
	{
		// y[0] ~ N(x[i], prec)
#pragma omp simd
		for (int i = 0; i < nx; i++) {
			result[i] = LOG_NORMC_GAUSSIAN + 0.5 * (lprec - (SQR(y[0] - x[i]) * prec));
		}
		char buf[64];
		snprintf(buf, sizeof buf, "theta loglike %.17g", result[nx-1]);
		write_string_to_file("log.txt", buf);
	}
		break;

	case INLA_CLOGLIKE_CDF:
	{
		// Prob(y[0] < x[i]) when y[0] ~ N(x[i], prec)
		double sprec = sqrt(prec);
		for (int i = 0; i < nx; i++) {
			double z = (y[0] - x[i]) * sprec;
			result[i] = CDF(z);
		}
		char buf[64];
		snprintf(buf, sizeof buf, "theta loglike %.17g", result[nx-1]);
		write_string_to_file("log.txt", buf);
	}
		break;

	case INLA_CLOGLIKE_QUIT:
		break;
	}


	return (ret);
}




/*
 * 模型： y ~ N(0, var = exp(x) + 1/tau)
 * theta[0] = log(tau)  (tau > 0)
 * ny 至少 1 列：y[0] = 观测值
 * nx 是需要评估的线性预测数组 x[0..nx-1]
 */
 double *inla_cloglike_stochvol(inla_cloglike_cmd_tp cmd, double *theta,
	inla_cgeneric_data_tp *data,
	int ny, double *y, int nx, double *x, double *result){
	double *ret = NULL;


    double ltau = NAN, tau = NAN;



	switch (cmd) {
	case INLA_CLOGLIKE_INITIAL: {
	ret = Calloc(1, double);
	ret[0] = 0;      /* #thetas */
	char buf[64];
	snprintf(buf, sizeof buf, "theta number %.17g", ret[0]);
	write_string_to_file("log.txt", buf);
	} break;

	case INLA_CLOGLIKE_LOG_PRIOR: {
	/* 先验：tau ~ Gamma(1,1)，对数变量需要加 Jacobian -> -tau + log(tau) = -exp(ltau) + ltau */
	char buf[64];
	snprintf(buf, sizeof buf, "theta initial start ", 0);
	write_string_to_file("log.txt", buf);
	ret = Calloc(1, double);
	ret[0] = -tau + ltau;
	
	snprintf(buf, sizeof buf, "theta initial %.17g", ret[0]);
	write_string_to_file("log.txt", buf);
	} break;

	case INLA_CLOGLIKE_LOGLIKE: {


		assert(ny >= 1);
		const double y0 = y[0];
		const double var_offset = 0.0;

		for (int i = 0; i < nx; i++) {
		/* 原函数里的 PREDICTOR_INVERSE_LINK(x[i], off) 这里采用最常见的 exp 链接 */
		const double var = exp(x[i]) + var_offset;
		result[i] = LOG_NORMC_GAUSSIAN - 0.5*log(var) - 0.5*(y0*y0)/var;
		}

	} break;

	case INLA_CLOGLIKE_CDF: {
	/* INLA 的 CDF: Prob(Ỹ < y | x,theta)；零均值正态 => Φ(y/√var) */
	char buf[64];
	snprintf(buf, sizeof buf, "CDF start %.17g", 0);
	write_string_to_file("log.txt", buf);
	assert(ny >= 1);
	const double y0 = y[0];
	const double var_offset =  0.0 ;

	for (int i = 0; i < nx; i++) {
	const double var = exp(x[i]) + var_offset;
	const double z = y0 / sqrt(var);
	/* Φ(z) = 0.5 * (1 + erf(z / sqrt(2))) */
	result[i] = 0.5 * (1.0 + erf(z * M_SQRT1_2));
	}
	buf[64];
	snprintf(buf, sizeof buf, "CDF %.17g", result[nx-1]);
	write_string_to_file("log.txt", buf);
	} break;

	case INLA_CLOGLIKE_QUIT:
	char buf[64];
	snprintf(buf, sizeof buf, "quit start %.17g", 0);
	write_string_to_file("log.txt", buf);
	default:
	break;
	}

	return ret;
}


/* 映射与导数：nu = 2 + exp(ldof) */
static inline double map_dof_forward_from_ldof(double ldof)     { return 2.0 + exp(ldof); }
static inline double map_dof_dforward_from_ldof(double ldof)    { return exp(ldof); }




/* 学生t：stochvol 的 cloglike 版本
   y / sqrt(exp(x)) ~ t_nu，等价于 y ~ (sqrt(exp(x))/sd) * t_nu，sd = sqrt(nu/(nu-2))
   超参数：theta[0] = log(nu - 2)  =>  nu = 2 + exp(theta[0]) > 2
*/
double *inla_cloglike_stochvol_t(inla_cloglike_cmd_tp cmd, double *theta,
	inla_cgeneric_data_tp *data,
	int ny, double *y, int nx, double *x, double *result){
	double *ret = NULL;

	/* 读超参数：nu */
	double ldof = NAN, nu = NAN, sd = NAN;
	if (theta) {
		ldof = theta[0];
		nu   = 2.0 + exp(ldof);                /* 确保 > 2 */
		sd   = sqrt(nu / (nu - 2.0));          /* 标准 t 的标准差 */
	}


	char bufhere[64];
	snprintf(bufhere, sizeof bufhere, "cmd=%s(%d)", cloglike_cmd_name(cmd), (int)cmd);
	write_string_to_file("log.txt", bufhere);




	switch (cmd) {
	case INLA_CLOGLIKE_INITIAL: {
		/* 返回 [M, init1, ..., initM]。这里只估一个超参：ldof=log(nu-2) */
		ret    = Calloc(2, double);
		ret[0] = 1;                 /* # thetas */
		ret[1] = log(4.0 - 2.0);    /* 初始 nu=8（可按需调整） */
	} break;

	case INLA_CLOGLIKE_LOG_PRIOR: {
		/*
		若想用平先验：直接 ret[0] = 0.0 即可。 */
		const double u = 15, alpha = 0.5;     /* 可改 */
		ret    = Calloc(1, double);
		char buf[64];
		snprintf(buf, sizeof buf, "theta prior %.17g", theta[0]);
		write_string_to_file("log.txt", buf);
		ret[0] = pc_logprior_ldof(theta[0], u, alpha);  /* theta[0] = ldof = log(nu-2) */
		snprintf(buf, sizeof buf, "theta after %.17g", theta[0]);
		write_string_to_file("log.txt", buf);
	} break;

	case INLA_CLOGLIKE_LOGLIKE: {
		/* log 密度：log f(y | x,nu) = const - ((nu+1)/2) * log(1 + (obs^2)/nu) - log(f)
		其中 f = sqrt(exp(x)) / sd，obs = y / f */
		assert(ny >= 1);
		const double y0   = y[0];
		const double lg1  = lgamma(nu * 0.5);
		const double lg2  = lgamma((nu + 1.0) * 0.5);
		const double logC = lg2 - lg1 - 0.5 * log(M_PI * nu);

		for (int i = 0; i < nx; i++) {
		const double var_u = exp(x[i]);        /* PREDICTOR_INVERSE_LINK(x) */
		const double f     = sqrt(var_u) / sd; /* 缩放 */
		const double obs   = y0 / f;
		result[i] = logC - 0.5 * (nu + 1.0) * log1p( (obs*obs) / nu ) - log(f);
		}
	} break;

	case INLA_CLOGLIKE_CDF: {
		/* 如需 CDF（用于 PIT），可用 GSL 的 gsl_cdf_tdist_* 或自己实现。
		这里先返回 NAN，表示未提供。 */
		for (int i = 0; i < nx; i++) result[i] = NAN;
	} break;

	case INLA_CLOGLIKE_QUIT:
	default:
	break;
	}
	return ret;
}


#undef LOG_NORMC_GAUSSIAN
#undef CDF
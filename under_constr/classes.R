# library(Rcpp)
# devtools::clean_dll()
Rcpp::compileAttributes()
devtools::load_all()


make_n_values <- function(.par_name, .n, .sd, .mean = NULL) {
    par_programmed <- c('aphid_density_0', 'aphid_surv_juv', 'aphid_surv_adult', 
                         'aphid_repro', 'K', 'disp_aphid', 'attack_surv')
    if (!.par_name %in% par_programmed) {
        stop(paste("parameter name", .par_name, "isn't programmed to be varied.",
                   "You have to do it manually."))
    }
    par_funs <- list(aphid_density_0 = rnorm, 
                     aphid_surv_juv = rbeta, 
                     aphid_surv_adult = rbeta, 
                     aphid_repro = rnorm, 
                     K = rbeta, 
                     disp_aphid = rbeta, 
                     attack_surv = rbeta)
    par_defaults <- list(
        aphid_density_0 = sap::populations$aphids_0,
        aphid_surv_juv = (sap::populations$surv_juv$high + 
                              sap::populations$surv_juv$low) / 2,
        aphid_surv_adult = (sap::populations$surv_adult$high + 
                                sap::populations$surv_adult$low) / 2,
        aphid_repro = (sap::populations$repro$high + 
                           sap::populations$repro$low) / 2,
        K = sap::populations$K,
        disp_aphid = sap::environ$disp_aphid,
        attack_surv = sap::wasp_attack$attack_surv / 2
    )
    if (is.null(.mean)) .mean <- par_defaults[[.par_name]]
    if (.par_name %in% c('aphid_density_0', 'aphid_repro')) {
        n_values <- rnorm(.n, .mean, .sd)
    } else {
        # Shape 1 and 2 for rbeta derived from .mean and .sd
        shape1 <- .mean * ((.mean*(1-.mean))/.sd^2-1)
        shape2 <- (1-.mean) * {.mean*(1-.mean)/.sd^2 - 1}
        n_values <- rbeta(.n, shape1, shape2)
    }
    
    # Left off --> Change above to the following (paste in LaTeX to see it):
    # For the `rbeta` parameters:
    # s_{i,j} &= \text{logit}^{-1}\left( \alpha_{i,j} \right) \\
    # \alpha_{i,j} &= \bar\alpha_i + \varepsilon_{\alpha,j} \\
    # \varepsilon_{\alpha,j} &\sim N \left(0,\sigma_{\alpha}^2 \right) \\
    # \bar\alpha_i &= \text{logit}\left( \bar s_i \right)
    # For the `rnorm` parameters:
    # r_{i,j} &= \text{exp}\left( \beta_{i,j} \right) \\
    # \beta_{i,j} &= \bar\beta_i + \varepsilon_{\beta,j} \\
    # \varepsilon_{\beta,j} &\sim N \left(0,\sigma_{\beta}^2 \right) \\
    # \bar\beta_i &= \text{log}\left( \bar r_i \right)
    
    
    
    # as.list(n_values)
    
}






z <- rbeta(10000, 20, 50); hist(z); mean(z); var(z)
20 / (20 + 50) # expected mean
(20 * 50) / {(20 + 50)^2 * (20 + 50 + 1)}  # expected variance
x = mean(z)
v=sd(z)^2
x * ((x*(1-x))/v-1) # expected shape1
(1-x) * { x*(1-x)/v - 1 }  # expected shape2
# v must be < x*(1-x)


x = sap::populations$K
v = sap::populations$K / 2
x * ((x*(1-x))/v-1) # expected shape1
(1-x) * { x*(1-x)/v - 1 }  # expected shape2

mean(rbeta(1e5, 0.0004665638, 0.9985994) > x)


aphid_density_0
sap::populations$aphids_0

aphid_surv_juv <- (sap::populations$surv_juv$high + sap::populations$surv_juv$low) / 2
aphid_surv_adult <- (sap::populations$surv_adult$high + 
                         sap::populations$surv_adult$low) / 2
aphid_repro <- (sap::populations$repro$high + 
                    sap::populations$repro$low) / 2
K <- sap::populations$K
disp_aphid <- sap::environ$disp_aphid
attack_surv <- sap::wasp_attack$attack_surv / 2

x <- (sap::populations$surv_adult$high + sap::populations$surv_adult$low) / 2
# x <- x[x > 0 & x < 1]
x <- x[1:21]
z <- log(x / (1-x))
z[1] <- 10; z[length(z)] <- -10
plot(z, cex = 1)

nls_df <- data.frame(surv = z, t = 1:(length(z)))

lmod <- lm(surv ~ t + I(t^2) + I(t^3) + I(t^4) + I(t^5), data = nls_df)
summary(lmod)
plot(predict(lmod), type = 'l')
points(surv ~ t, data = nls_df)

nlsfit <- nls(surv ~ 1 / (a * t) - exp(b * t), 
    start = list(a = 5, b = 1),
    lower = list(a = 0.1, b = -0.32156),
    data = nls_df, algorithm = 'port')
coef(nlsfit)



a = 0.5; b = -0.2
f <- function(x) 1 / (a * x) - exp(b * x)

curve(f(x), 1, length(z))


curve(exp(f2(x)) / (1 + exp(f2(x))), 0, 18)
lines(0:18, x[x>0 & x<1], col = 'red')
# plot(0:18, f2(5.484918, 0.142735, 0.000001, 0:18), type = 'l')



#





n_pops = 2
n_patches = 2
pl <- all_pop_lists(n_pops, n_patches, 
                    aphid_density_0 = list(
                        (1 - sap::populations$prop_resist) * sap::populations$aphids_0,
                        sap::populations$prop_resist * sap::populations$aphids_0),
                    sigma_x = 0, sigma_y = 0, rho = 0, demog_mult = 0,
                    attack_surv = list("susceptible", "resistant"), 
                    aphid_surv_adult = list('high', 'low'),
                    aphid_repro = list('high', 'low'))
sp <- new(SimPatches, pl, 
          rep(sap::environ$cycle_length, n_patches), c(0,15))
out <- sp$simulate(630, rng_seed = sample.int(2^31-1,1))
out

par(mar = c(b = 2, l = 4, t = 1, r = 1))
sap::populations$surv_juv$high

plot(sap::populations$surv_adult$high[1:21], type = 'l', ylab = 'surv_adult')
lines(sap::populations$surv_adult$low[1:21], col = 'red')
lines((sap::populations$surv_adult$high[1:21] + sap::populations$surv_adult$low[1:21])/2, 
      col = 'green')

plot(sap::populations$repro$high[1:24], type = 'l', ylab = 'repro')
lines(sap::populations$repro$low[1:24], col ='red')
lines((sap::populations$repro$high[1:24] + sap::populations$repro$low[1:24])/2, 
      col ='green')


sap::wasp_attack$attack_surv

# mixed stock analysis uses the method of Bolker et al. 2003, Ecol. App., 13:763-775
# see http://www2.uaem.mx/r-mirror/web/packages/mixstock/vignettes/mixstock.pdf 
# and https://github.com/bbolker/mixstock for information on the package mixstock

#### IF not already installed, install package mixstock using this command:
#  devtools::install_version("mixstock","0.9.5.1")
####

library(dplyr)
library(mixstock)

load("data/hap.freqs.rda")

# combining GRO.OAX and CentAm data as one source population
hap.freqs$CentAm <- hap.freqs$CentAm + hap.freqs$GRO.OAX
hap.freqs <- select(hap.freqs, c(CentAm, MMex, MCP))

# This analysis requires me to remove the 8 haplotypes that are found in MCP, not the other strata
# They are uninformative and cause an error
hap.dat.mcmc <- as.mixstock.data(hap.freqs[-which((hap.freqs[,1]+hap.freqs[,2])==0),])

# Run MCMC
hap.mcmc <- tmcmc(hap.dat.mcmc, n.iter = 40000, verbose = TRUE)

# Check for convergence
diag.1 <- calc.mult.RL(hap.dat.mcmc)
diag.2 <- calc.GR(hap.dat.mcmc, tot = 40000)

#format results
mcmc.res <- cbind(hap.mcmc$fit$input.freq, confint(hap.mcmc))

write.csv(mcmc.res, file = "Mix-stock_mcmc.csv")
write.csv(diag.2[[1]], file = "Gelman-Rubin_convergence_check.csv")

library(dplyr)
library(tidyverse)
library(mixstock)
library(ggplot2)
library(rsoi)
library(grid)
library(gridExtra)
library(cowplot)


load("data/hap.freqs.rda")

MCP.by.year.hapfreq <- read.csv("data-raw/MCP-by-year-hapfreqs.csv", header = TRUE, check.names = FALSE)
hap.freqs <- left_join(hap.freqs, MCP.by.year.hapfreq)
hap.freqs[is.na(hap.freqs)] <- 0

hap.freqs$CentAm <- hap.freqs$CentAm + hap.freqs$GRO.OAX
yrs <- names(MCP.by.year.hapfreq)[-c(1)] 

res.list <- lapply(yrs, function(y) {

  hap.freqs <- select(hap.freqs, c(CentAm, MMex, y)) 

  # This analysis requires me to remove the haplotypes that are found in MCP, but not the other strata
  hap.dat.mcmc <- as.mixstock.data(hap.freqs[-which((hap.freqs[,1]+hap.freqs[,2])==0),1:3])

  # Run MCMC  
  hap.mcmc <- tmcmc(hap.dat.mcmc, n.iter = 40000, verbose = TRUE)
  
  # Calculate diagnostics
  diag.1 <- calc.mult.RL(hap.dat.mcmc)
  diag.Gelman.Rubin <- calc.GR(hap.dat.mcmc, tot = 40000)
  
  # Format results and return
  mcmc.res <- cbind(hap.mcmc$fit$input.freq, confint(hap.mcmc))
  mcmc.res <- c(mcmc.res[1,], mcmc.res[2,])
  names(mcmc.res) <- c("CentAm", "CentAm 2.5", "CentAm 97.5", "MMex", "MMex 2.5", "MMex 97.5")
  
  return(list(diag.1,diag.Gelman.Rubin,mcmc.res,hap.mcmc))
})
names(res.list) <- yrs

mcmc.densities <- data.frame((do.call('cbind', lapply(res.list, function(y){y[[4]]$resample$contrib.MMex}))))
names(mcmc.densities) <- yrs
mcmc.densities[,7] <- NA #excluding 2016 due to low sample size
write.csv(colMeans(mcmc.densities), file = "results/Mix-stock_by-year_mcmc.csv")

### PLOT RESULTS AS BOXPLOTS OF PROBABILITY DENSITY
mcmc.densities <- pivot_longer(mcmc.densities, cols = everything())
p <- ggplot(mcmc.densities, aes(x = name, y = value)) + geom_boxplot() +
  labs(x = NULL, y = "Contribution from BB") +
#  scale_y_continuous(expand = c(0,0)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size=26, face="bold"),
        axis.title.y = element_text(size=26, face="bold"),
        axis.text = element_text(size = 20),
        panel.grid = element_blank(),
        axis.line = element_line())
tiff(filename = "results/Boxplot_mix-stock_by-year_mcmc.tif", width = 1000, height = 600)
p
dev.off()

save.image(file = "results/Mix-stock_by-year.rda")

### DOWNLOADING ENSO DATA FOR PLOTS
enso <- rsoi::download_enso() %>% 
  filter(Date > "2009-05-01") %>% filter(Date < "2021-09-01")
enso$rescaled <- 0.5 + (enso$ONI/5.28)

### ENSO PLOT TO DISPLAY UNDER BOXPLOTS
enso.plot = ggplot(data = enso, aes(x = Date, y = ONI, colour = phase, fill = phase))+
  #  geom_col(colour = "lightgrey", fill = "lightgrey")+
  #  geom_col(colour = "#999999")+
  geom_col() + 
  #  xlim(min = "2009-07-01", max = "2021-06-01") +
  scale_colour_manual(values = c("#56B4E9", "#999999", "#D55E00"), name = "") +
  scale_fill_manual(values = c("#56B4E9", "#999999", "#D55E00"), name = "") +
  theme_minimal() +
  labs(x = "Year", y = "ONI") +
  scale_y_continuous(breaks = c(-2.0, -1.0, 0, 1.0, 2.0),
                     labels = c("-2.0", "-1.0", "0", "1.0", "2.0"),
                     expand = c(0,0)) +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y",
               expand = c(0,0)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=26, face="bold"),
        axis.text.y = element_text(size = 20), 
        legend.position = "none",
        panel.grid = element_blank()) 
enso.plot
tiff(filename = "results/5358.R2_Fig5.tif", width = 1000, height = 600)
ggdraw() + draw_plot(p, 0, 0.25, 1, 0.75) + 
  draw_plot(enso.plot, 0,0,1,0.25) + 
  draw_plot_label(c("(a)","(b)"), c(0, 0), c(1, 0.25), size = 16)
dev.off()


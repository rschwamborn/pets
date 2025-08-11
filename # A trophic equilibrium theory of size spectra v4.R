
# Simulations and calculations for a 
# mass-specific “predator-prey-efficiency
# theory of size spectra” (PETS) 

# Ralf Schwamborn
# version 4.0
# August 11, 2025

# R script for the paper "Towards a compleat theory of ecosystem size spectra" 


# Numeric examples and simulations

# 1. First conversions -------------

PPMR = 121 # from Figueiredo et al., 2020

log10PPMR = log10(PPMR)  
log10PPMR # 2

TE = 0.1
log10(TE)  # -1

log10TE = log10(TE)  

m = 0.75

c = log10TE/log10PPMR
c # -0.5

b = c - m # size spectrum slope "b"

b= -1.23 # size spectrum slope "b"

# "b" values in De Figueiredo et al., 2025: 
# "b" values between -0.82 to -0.95 (NBSS) and -1.74 to  -1.9 (NNSS)

# 2. understanding log(a)/log(b) ratios:               ------------

# log(A)/ log(B) = logA(B)

# examples:
a = 20
b = 10
c = 100
d =  2

log10(a) # 1.3
log10 (b) # 1
log10 (c) #  2

a/b # 2
log10(a)/ log10(b) # 1.3

log2(a)/ log2(b) # 1.3
log10(a) # 1.3 

log2(10)/ log2(b) # 1.3
log10(10)/ log10(b) # 1.3
log(10)/ log(b) # 1.3

log10(10)/ log10(2) # 3.2
log2(10) # 3.2

#log(x, base)
log(10, base = 2) # log og x with any possible base
log(100, base = 0.1) # log og x with any possible base
log(1000, base = 0.1) # log og x with any possible base
log(100, base = 0.5) # log og x with any possible base
log(10000, base = 0.01) # log og x with any possible base
log(5, base = 0.1) # log og x with any possible base
-0.5-0.75
-2-0.75

# 3. Example data -------------- 
# Median NBSS (carbon/carbon) biomass data across the Altantic
# From Schwamborn et al. (submitted)

xvec <- 
  c(5.68e-14, 1.14e-13, 2.27e-13, 4.55e-13, 9.09e-13, 1.82e-12, 
    3.64e-12, 7.28e-12, 1.46e-11, 2.91e-11, 5.82e-11, 1.16e-10, 2.33e-10, 
    4.66e-10, 9.31e-10, 1.86e-09, 3.73e-09, 7.45e-09, 1.49e-08, 2.98e-08, 
    5.96e-08, 1.19e-07, 2.38e-07, 4.77e-07, 9.54e-07, 1.91e-06, 3.81e-06, 
    7.63e-06, 1.53e-05, 3.05e-05, 6.1e-05, 0.000122, 0.000244, 0.000488, 
    0.000977, 0.001953, 0.003906, 0.007813, 0.015625, 0.03125, 0.0625, 
    0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024
  )
# dput(X_vector_gCind)

xvec_log <- log10(xvec)
yvecNBSS <- median_vec_PHYTOP<- c(NA, NA, NA, NA, NA, NA, 241500000, 250500000, 8.9e+07, 37500000, 
                                  22950000, 15150000, 7e+06, 2730000, 1575000, 545000, NA, NA, 
                                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                                  NA, NA, NA, NA, NA)
# dput(median_vec_PHYTOP)
yvecNBSS_log <- log10(yvecNBSS)
# dev.off()
plot(yvecNBSS_log ~ xvec_log) 
lm0 <- lm(yvecNBSS_log ~ xvec_log)
summary(lm0)# -0.99 # slope = -1,NBSS, Based on Biomass Data    
#NBSS: log(Biomass) vs log(individual carbon biomass)
 
abline(lm0)


# 2. Simulations and transformations 
# - comparing units and models 
# A. NBSS (log(biomass) / log(weight) )
# B. NNSS (log(abundance) /  log(weight) 
# C. Abundance-size spectra (log(abundance) / log(linear size) )


# converting frm NNSS slopes to NBSS slopes and vice versa

# convert biomass to abundance biomass   B (g * m-3) = A(ind. m-3) * w (g/ind-1)  
#   B =  A * w  ;    A = B / w
# Question: When normalizing (dividing by bin width) do you actually convert from biomass to abundance? NO, it is just a crecrtion for the effects of varying bin sizes!
# Preliminarily Ignore normalization... 
# just convert Biomass to ABUNDANCE , by dividing by bin mids (gC ind-1)


yAbund_data <- yvecNBSS / xvec
yAbund_data_log <- log10( yAbund_data)

plot(yAbund_data_log ~ xvec_log) 
lmAbund <- lm(yAbund_data_log ~ xvec_log)
summary(lmAbund)# -0.99 # slope = -1,NBSS    
abline(lmAbund)
# slope = -2 !!! #  Based on Abundance Data (NNSS) 
#NNSS : log(abundance) vs log(individual carbon biomass)

# Abundance-size spectrum
# comvert biomass to linear size

# biomass ~ size^3
# size = cubic root (biomass),   cube root : x^(1/3)
x_vecsize = xvec^(1/3)

x_vecsize_log <- log10( x_vecsize)

plot(yAbund_data_log ~ x_vecsize_log) 
lmAbund_size <- lm(yAbund_data_log ~ x_vecsize_log)
summary(lmAbund_size)# -6# slope = -6,Abundance-size spectrum    
abline(lmAbund_size)


# Biomass-size (cm) spectrum
# convert biomass to linear size

# biomass ~ size^3
# size = cubic root (biomass),   cube root : x^(1/3)
x_vecsize = xvec^(1/3)
x_vecsize_log <- log10( x_vecsize)

plot(yvecNBSS_log ~ x_vecsize_log) 
lmbiom_size <- lm(yvecNBSS_log ~ x_vecsize_log)
summary(lmbiom_size)# -3# slope = -3,Biomass-size(cm) spectrum    
abline(lmbiom_size)

beta = E * S

beta = E / log(PPMR)


PPMR = 10^4

E  = -1   * log10(PPMR)
#  E = -4


# PPMR  ….  canbe cpnctes into the more convenient trophic scaling parameter (TS)
# TS =    log(TLi+1) - log(TLi)  /  log(wi+1) - log(wi)
# 
# From Figueiredo et al.,  2020:
# Ordinary least squares linear regression was used to investigate the
# linear relationship between log10-transformed average body size (ESD,
#                    μm) and trophic level (TL, estimated with TEF = 3.2 and 2.3). Body
# size (ESD) was log10-transformed to obtain linear relationships for
# analysis and to improve homoscedasticity. The slope of this linear regression
# model was used to estimate the average predator/ prey size
# ratio (PPSR) and predator/ prey mass ratio (PPMR), using the following
# equations: PPSR = 10(1/slope), if log10 is being used in the linear model
# (Hunt et al., 2015), and PPMR = PPSR3, assuming isometry and sizeinvariant
# density (Lins Silva et al., 2019).

#MeanPPMRfor the ARCH pelagic foodweb was 1.48 *  10^5:1,
PPMR_HUnt <- 1.48 *  10^5 # 148000
# intosizegroupsitbecameapparentthathighPPMRsintheARCH
# wereprimarilydrivenbythenekton.Mesozooplanktonand
# macrozooplankton/micronektonPPMR were1–2 orders of magni-
#   tude lower than the nekton PPMR

# PPMR = m(TL1i+1) / m(Tli)

# PPMR = m(TL1i+1) / m(Tli)

# PPMR = (mi+1 / mi)   
# TWS = S =     TLi+1 - TLi / log(mi+1) - log(mi))
# S = 1/ log10(PPMR)
x =  1:100000; y =  1/ (log10(x))
plot( x , y )
# S vas from 3.3 (PPMT = 1) to S = 0.33 (PPMR = 100000) within reasonble limits 

1/ (log10(100000)) #PPMR = 100000, S = 0.25

# if  PPMR = 10000
  S = 1/ log10(PPMR) # S = 0.25

  beta = -1
  
  # beta = E * S
  
  E = beta / S

  E # E = -4
  # can E actually be negative? YES!
  
  E =  log10(1) - log10(10000) / 1
  # E = -4, OK
  # Yes, E can and should always be negatove, OK
  # E =  log10(10) - log10(10000) / 1      
  # E =  log(Mi+1) - log(Mi)  /  TLi+1 - TLi
  # E =  log10(1) - log10(10000) / 1      
  
  # Numerocal example
  
  PPMR_HUnt # 148000
  S_Hunt <-  1/ log10(PPMR_HUnt) # S =  0.19
    E_Hunt = beta / S_Hunt  #E_Hunt = -5.17
    E =  log10(1) - log10( 1.48 *  10^5) / 1
    
     
size <- 1:200
w <- size^3
size = w^(1/3)
#slope_TL_logsize <- d(TL) / d(log(size))
slope_TL_logsize <-  0.1603

PPSR = 10^(1/slope_TL_logsize)

PPMR = PPSR^3

# PPMR <- 10^(1/slope_TL_logsize)^3

PPMR <- 10^(1/slope_TL_logsize)^3

slope_TL_logsize= 3 / (log10 (PPMR))
  

PPMR = 120

slope_TL_logsize = 1/ (log10(PPMR))^1/3
TSS <- slope_TL_logsize # 0.1603
TWS <- (1/3)  *TSS # 0.053
log10(0.05) # -1.3


PPMR_vec <- 1:10000
slope_TL_logsizevec = 1/ (log10(PPMR_vec))^1/3
TSSvec <- slope_TL_logsizevec # 0.1603
TWSvec <- (1/3)  *TSSvec # 0.053

plot (TWSvec ~ PPMR_vec)


# PPMR  ….  can be convetted  into the more convenient trophic scaling parameter (TS)
# example 3
#TWS =    (TLi+1) - (TLi)  /  log(wi+1) - log(wi)
# TWS =    1  /  log(10) - log(1) 
# TWS = 0.4, realistic example, ten times higheer weight per TL 

# example 4
#TWS =    (TLi+1) - (TLi)  /  log(wi+1) - log(wi)
#TWS =    1  /  log(3) - log(1) 
# TWS = 0.9, realistic example, tthree times higheer weight per TL 

# example 5
#TWS =    (TLi+1) - (TLi)  /  log(wi+1) - log(wi)
#TWS =    1  /  log(100) - log(1) 
# TWS = 0.217, realistic example, 100 times higheer weight per TL 

# exampl 6
#TWS =    (TLi+1) - (TLi)  /  log(wi+1) - log(wi)
#TWS =    1  /  log10(1000) - log10(1) 
# TWS = 0.333, realistic example, 1000 times higher weight per TL 

# example 7
#TWS =    (TLi) - (TLi+1)  /  log(wi) - log(wi+1)
#TWS =    -1  /  (log10(1) - log10(10000) ) 
# TWS = 0.25, exteme  example, 10000 times higher weight per TL 
# TWS depods o the base of thr logarithm used! use another eqution with log(TLi)?
# Solution : use th sam base throughghout ! E.g. always base 10, or always base 2!

# example 7, 10000 times higher weight per TL
# mte =  c * tws 
mte = -0.25 * 0.109 # mte = -0.027 # negtive mte?

#beta = log(Mi+1) - log(Mi)  /  log(wi+1) - log(wi)
#beta = (log10(10) - log10(100))  /  (log10(10) - log10(1))
# beta = -1 , OK

#beta = (log(10) - log(100))  /  (log(10) - log(1))
# beta = -1 , OK
# beta is -1 regrdless og th bas of thr logaithm used

#TWS =    1  /  log(10) - log(1) # tws = 0.4,dTL / d(log(weights))
#mte = -1 / TWS # -2


# Empirical solution  to the trophic equilibrium problem

# PPMR =  wPred / wPrey
# TE = Mi / Mi+1
# w - TL relationship 

beta = -1
#TWS = 0.05

#TWS = (TLi+1) - (TLi)  /  log(wi+1) - log(wi)
#beta =    log(Mi+1) - log(Mi)  /  log(wi+1) - log(wi)

#E =  log(Mi+1) - log(Mi)  /  TLi+1 - TLi


# Simple hypothetic example:

# MTE
# MTE_raw =  10 / 1    
# MTE_raw =  (Mi+1) - (Mi)  /  TLi+1 - TLi, ten-fold increase in total biomass = one TL step 
MTE_log =  1  
# MTE_log =  log10(10) - log10(1) / 1        #MTE_log =  log(Mi+1) - log(Mi)  /  TLi+1 - TLi
# MTE_log =  log10(10) - log10(1000) / 1        #MTE_log =  log(Mi+1) - log(Mi)  /  TLi+1 - TLi

E =  log10(10) - log10(1000) / 1
# E = -2, OK

#    TWS
# TWS = (TLi+1) - (TLi)  /  log(wi+1) - log(wi)
TWS <-  1 # <-   1   /  (log10(10) - log10(1) ) # ten-fold increase in individall weoght with TL

c = MTE_log  / TWS # = 1/1 = 1

TWS = 0.05


beta = c - m

m = 0.75
TWS = 0.05

c = -0.25

 MTE = 2.12


# loglogTE  =  log(gi+1) - log(gi)  /  log(wi+1) - log(wi),   where g is grwoth (i.e,  g =  (dw/dt))
# logTE  =  (gi+1) - (gi)  /  log(wi+1) - log(wi),   where g is grwoth (i.e,  g =  (dw/dt))
# MTE' =  Mi+1 - Mi  /  TLi+1 - TLi


beta = -1; m = 0.75
beta <- c - m

c <-   (log10(MTE)) / (log10(TWS))  

c_4 <-   (log10(MTE)) / ((TWS))  # c = 6.5 

TWS = 0.05

#c = -0.25

MTE = 2.12

# if MTE = 2.12,  
#MTE_log =  log(Mi+1) - log(Mi)  /  TLi+1 - TLi

MTE_log =  log10(3) - log10(1) / 1
log(MTE_log)

# EXAMPLE 1:
MTE_log =  log10(132) - log10(1) / 1 #, MTE = 2.12,
# if MTE = 2.12,  each trophic level is 132 times less bimoass than the previos one
# very low efficiency!

# if TWS = 0.05, ... 
PPMR = 120 # predator is 120 times mpre massive than its prey, 


# anoter simpe example (example 2)

mte = -0.0125
beta = -1
tws = 0.05 
  beta = mte  / tws

  MTE_log =  log10(0.97) - log10(1) / 1 #, MTE = -0.01,
  # if MTE = 2.12,  each trophic level is 132 times less bimoass than the previos one
  # very low efficiency!
  

# Plot examples with low and high TE! 
# ea scenrio with three plots: NBSS (beta), TL_size (PPMR) and Biomass_TL (MTE)


#   PLOT  Scenario 1: high PPMR, high TE, low S, neta = -1, b = beta
  #   PLOT  Scenario 2: low PPMR, high TE, low SIGMA, b = beta
    #   PLOT  Scnerio 3: with s unequal zero
  # neg s   (stress throgh fisheries), positive stress (nutrients),positive size specific stress effect (nlionfish 
  # size specific stress effect ("s")
  
  
  
  #   PLOT  Scenario 1: high PPMR, high E, low S, beta = -1, b = beta
  
  # NBSS
x <- xvec_log; y <- yvecNBSS_log
yvecNBSS_log <- log10(yvecNBSS)
# dev.off()
plot(yvecNBSS_log ~ xvec_log, xlim = c(-12, -8))
lm0 <- lm(yvecNBSS_log ~ xvec_log)
abline(lm0)
summary(lm0)# -0.99 # slope = -1,NBSS, Based on Biomass Data in Schwamborn et al, submitted.    
#NBSS: log(Biomass) vs log(individual carbon biomass)
# beta = -1
  # b = -0.99,a = -2.78, R-squared:  0.99
# b = beta, s = 0



  #   PLOT  Scnenario 2: low PPMR, high E, low SIGMA, b = beta
  
  
  #   PLOT  Scnerio 3: with s uneuall zero
  # neg s   (stress throgh fisheries), positive stress (nutrients),positive size specific stress effect (nlionfish 
  # size specific stress effect ("s")
  
  

#  Empirical solution  to the trophic equilibrium problem, II 
# lookig at the prationship of EE (bt TL), PSTE 
# and ME (bt wigh class)
 


# 3.b convert  PPMR to FCL, wha tis the FCL - PPMR relationship
   #solve empirically


Example1 <- data.frame(
  TL = c(1, 2, 3),
  w = c(10, 1e7, 1e13)
)

Example2 <- data.frame(
  w = c(10, 1e5, 1e9, 1e13),
  TL = c(1, 2, 3, 4)
)

Example3 <- data.frame(
  w = c(10, 1e4, 1e7, 1e10, 1e13),
  TL = c(1, 2, 3, 4, 5)
)


Example4 <- data.frame(
  w = c(10, 1e3, 1e5, 1e7, 1e9, 1e11, 1e13),
  TL = c(1, 2, 3, 4, 5, 6, 7)
)

Summary_table <- data.frame(
  example = c(1, 2, 3, 4),
  FCL = c(3, 4, 5, 7),
  PPMR = c(1e6, 1e4, 1e3, 1e2),
  logPPMR = c(6, 4, 3, 2)
)



# Step 3: Combine all into one list
FCL_PPMR_tables <- list(
  Example1 = Example1,
  Example2 = Example2,
  Example3 = Example3,
  Example4 = Example4,
  Summary_table = Summary_table
)

# Save the list as an .RData file for GitHub 

save(FCL_PPMR_tables, file = "FCL_PPMR_tables.RData")


#Create a data.frame 

df <- data.frame(
  example = c(1, 2, 3, 4),
  FCL = c(3, 4, 5, 7),
  PPMR = c(1e6, 1e4, 1e3, 1e2),
  logPPMR = c(6, 4, 3, 2)
)


# Given data (excel file with 4 simulated food chains)
FCL <- c(3, 4, 5, 7)
logPPMR <- c(6, 4, 3, 2)

invlogPPMR = 1/ logPPMR


# log(wtop /wbase) = 12

FCL_calc2 <-     12 * invlogPPMR
# Perfect fit 

FCL_calc2 
FCL


# FCL  =  log(wtop /wbase) * ( 1 / logPPMR ) !!!


plot( FCL  ~      invlogPPMR)

summary(m55 <- lm( FCL  ~      invlogPPMR +0 ) )
# slpe = 14.7692
abline(m55)

+14.77/12 # 1.23


summary(m56 <- lm( FCL  ~      invlogPPMR) )


# FCL  =   1 + ( log(wtop /wbase)  * ( 1 / logPPMR ) ) !!!

FCL_calc_ok  =   1 + (12  * invlogPPMR) # OK !!!
FCL_calc_ok
FCL
# OK!

# from populations to size spectra (mizer package)

library(mizer)
params <- newMultispeciesParams(NS_species_params, NS_interaction)
sim <- project(params, t_max = 10, effort = 0)

plot(sim)

summary(sim)


# Simulate A flat, predation-free size spectrum slope (bottom-up control)


-0.75 # Flat, predation-free abundance-weight spcertrum slop

-2 # Observed abundance-weight spectrum slope


-1 # Observed Biomass-weight spectrum slope


# -1 *( -0.75 / -2) = -0.375   # eq 276.x
# -0.375  = slope of a flat, predation-free Biomass-weight spectrum 

# check eq 276.x with numerical examples :



# 4. Fom M to P and P/B spectra  ------------
#   
# Mass specrum
#  Production spctrum
# Turnover spectrum
#   


#  how aare b_M,  beta_P andbetaaPB relted. Cam   we prpct this relltionship from Zpred?
## are these ratio  commstant? 

#For the sake of simplicity, we will preliminarily assume a 
# constant Kleiber’s law (ref Kleiber ) scaling of μ = 0.75, as in West et al., (1997) across the ecosystem.


# Mass-weight spectrum
M_vec <- yvecNBSS
logM_vec <- log10(M_vec)
plot(M_vec ~ xvec_log, xlim = c(-12, -8))
lmM <- lm(M_vec ~ xvec_log)
summary(lmM)
abline(lmM)

# the ubiquitous -1 slope of the Mass-weight spectrum can be used to 
# make predictions for the  Production-weight spectrum (-1.75)
# and for the Turnover-weight spectrum (-0.75)
# using first principles and MTE. 



#scenario 0: P/M = w^0.25

#  Turnover-weight spectrum,  Turnover = Production/Mass ratio each size bin  
PM_vec = M_vec^0.25
logPM_vec <- log10(PM_vec)
plot(logPM_vec ~ xvec_log, xlim = c(-12, -8))
lmPM <- lm(logPM_vec ~ xvec_log)
summary(lmPM)
abline(lmPM)
# slope of -0.75 (Kleiber's law)
#  Turnover-weight spectrum:  slope of 0.25 (1 - Kleiber's law)
# OK


# Production-weight spectrum, Production in each size bin
P_vec <-  PM_vec * M_vec
logP_vec <- log10(P_vec)
plot(logP_vec ~ xvec_log, xlim = c(-12, -8))
lmP <- lm(logP_vec ~ xvec_log)
summary(lmP)  # apprx -1.25  
abline(lmP)
# Slope of the Production-weight spectrum: approx - 1.75 , (ie. - 1 * ( 1 + Kleiber's law))

-1 + -0.25

#scenario 1: P/M = w^0.75

#  Turnover-weight spectrum,  Turnover = Production/Mass ratio each size bin  
PM_vec = M_vec^0.25
logPM_vec <- log10(PM_vec)
plot(logPM_vec ~ xvec_log, xlim = c(-12, -8))
lmPM <- lm(logPM_vec ~ xvec_log)
summary(lmPM)
abline(lmPM)
# slope of -0.75 (Kleiber's law)
#  Turnover-weight spectrum:  slope of 0.25 (1 - Kleiber's law)


# Production-weight spectrum, Production in each size bin
P_vec <-  PM_vec * M_vec
logP_vec <- log10(P_vec)
plot(logP_vec ~ xvec_log, xlim = c(-12, -8))
lmP <- lm(logP_vec ~ xvec_log)
summary(lmP)  # -0.74260   
abline(lmP)
# Slope of the Production-weight spectrum: approx - 1.75 , (ie. - 1 * ( 1 + Kleiber's law))


# Scenario 2: P = w^0.75 - leads to an inverted Turnover-weight spectrum!
#           not realistic!

#  Turnover-weight spectrum,  Turnover = Production/Mass ratio each size bin  
PM_vec = P_vec / M_vec
logPM_vec <- log10(PM_vec)
plot(logPM_vec ~ xvec_log, xlim = c(-12, -8))
lmPM <- lm(logPM_vec ~ xvec_log)
summary(lmPM)
abline(lmPM)
# slope of 0.25 (1 - Kleiber's law)
#  Turnover-weight spectrum:  slope of 0.25 (1 - Kleiber's law)


# Production-weight spectrum, Production in each size bin
P_vec <- M_vec^-0.75
logP_vec <- log10(P_vec)
plot(logP_vec ~ xvec_log, xlim = c(-12, -8))
lmP <- lm(logP_vec ~ xvec_log)
summary(lmP)  # -0.74260   
abline(lmP)
# Slope of the Production-weight spectrum: approx -0.75 (Kleiber's law)


# Cleared volume vs body size (weight)  ----------

# L ind⁻¹ d⁻¹

# cleared volume (also known as clearance rate, typically in mL ind⁻¹ d⁻¹ or L ind⁻¹ d⁻¹) and body weight  in copepods 

# Huntley & Lopez (1992) 
# Title: Temperature-dependent production of marine copepods
# Journal: Limnology and Oceanography

#example 
a = 2
b = 0.79   # b = 0.79

w = 10^ (1:10)

Vol_cleared <-  a * w^b

plot(Vol_cleared ~ w)
plot(log10(Vol_cleared) ~ log10(w))

# calculate N per 10^9 liters ( 1 cubic kilometer (km³))
N_ind.per.cubic.km = 10^9/ Vol_cleared  

# N = 6.3 ind  to 81 million organisms per cubic kilometer

plot( log10(N_ind.per.cubic.km) ~ log10(w)  )
summary(lm( log10(N_ind.per.cubic.km) ~ log10(w)  )) # slope = -0.79
# Abundance-wepght spctum is Kleiber-scaled, when botttom-up control

# Biomass spcrum: B = N * w

B_mu_g.per.cubic.km <- N_ind.per.cubic.km * w

plot( log10(B_mu_g.per.cubic.km) ~ log10(w)  )
summary(lm( log10(B_mu_g.per.cubic.km) ~ log10(w)  )) # slope = -0.79
# Abundance-wepght spctum is Kleiber-scaled


# slope of abundance - weight-spectrum based on cleard volume = b = -0.79

# slope of abundance - weight-spectrum observed:  -2 (Figueiredo et al., 2025)


# biomass predicted from NBSS slope = -1
B_NBSS_miuns1slope = w^-1   
#plot( log10(BNBSSmiuns1slope) ~ log10(w)  )
# summary(lm( log10(BNBSSmiuns1slope) ~ log10(w)  )) # slope = -1, OK

# abundance  predicted from -1 slope, Abund = Biomass /w 
Abund_NBSS_miuns1slope  = BNBSSmiuns1slopoe / w
plot( log10(Abund_NBSS_miuns1slope) ~ log10(w)  )
summary(lm( log10(Abund_NBSS_miuns1slope) ~ log10(w)  )) # slope = -2
# Abundance-weight spectrum slope = -2, OK as in Figueiredo et al, 2025


# Comments and unsolved problems

# PPMR increases with predator mass
# TE increases with predator mass and with PPMR

# Why does TE increase with predator mass? 
# Why does TE decrese with PPMR? Predation is most efcient when PPMR = 1 (low abnce? low pprey gooteth ?
# that is efficnt? ), less efficnt whenPPMTR is low
# size-specific trophic transfer efficiency TE (e, the ratio of the production of a
                #      mass category to that of its prey)
# Ecotrophic efficiency EE (e, the ratio of the production of a
#     tr
# Labrador  data 
# Globally Consistent Quantitative Observations of Planktonic Ecosystems
# April 2019Frontiers in Marine Science 6:196
# DOI: 10.3389/fmars.2019.00196
# LicenseCC BY 4.0
# Labs: Department of Marine Snow & Plankton ResearchJan Schulz's Lab
# Lombard FabienLombard FabienEmmanuel S BossEmmanuel S BossAnya M. WaiteAnya M. WaiteShow all 41 authorsWard AppeltansWard Appeltans


# Labrador <- read.csv("~/Papers/000 - generalized size spectra model of marine ecosystems/Labrador.csv")
#View(Labrador)

# Mass-esd size spectrum
#plot( log10(Labrador$y) ~ log10(Labrador$x ) ) 
#summary(lm( log10(Labrador$y) ~ log10(Labrador$x ) ) )
# slope = -3.4135 

# Size spectra slopes:  
# Mass-weight = -1
# Mass-size = -3 (e.g. Labrador sea example in Lombard et al., 2019)

# Abundncd-weight = -2
# Abundncd-size =  -6 ?




#References:

# Figueiredo, G. G. A. A., Schwamborn, R., Bertrand, A., Munaron, J. M., Le Loc’h, F. (2020). Body size and stable isotope composition of zooplankton in the western tropical Atlantic. J. Mar. Syst. 212, 103449. https://doi.org/10.1016/j.jmarsys.2020.103449
# De Figueiredo, G. G. A. A., de Albuquerque Lira, S. M., Bertrand, A., Neumann-Leitão, S., & Schwamborn, R. (2025). Zooplankton abundance and biovolume size-spectra in the western tropical Atlantic-From the shelf towards complex oceanic current systems. Marine Environmental Research, 204, 106906.

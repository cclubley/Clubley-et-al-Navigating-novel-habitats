# The following script carries out the analysis and presentation of results (Figures 4 & 5) from the laboratory experiments described and presented in "Navigating novel habitats: assessing the impacts of artificial complexity on space-use by a key intertidal grazer" by Clubley et al. 

# Data used are contained in the .csv file 'Master_limpet_movement.csv'

{library(dplyr)
library(ggpubr)
library(ggplot2)
library(pracma)
library(nlme)
library(performance)
library(emmeans)
library(multcomp)
library(plotrix)
library(viridis)}

# -------------------------------------------------------------------------
# ---------------------------- Gross distance -----------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# Gross (total) distance was calculated by summing the distance moved by an individual limpet between each timestep (5-minute) of the 24-hour period over which they were recorded.

# Load the movement data
limpet <- read.csv("./Master_limpet_movement.csv")

# Data description
head(limpet)
# Trial_number: Aribitrary, a product of the Ethovision analysis
# Replicate: The number of the experimental run (note: Replicate starts from 3 as replicates 1 & 2 were trial runs during which errors occurred with the equipment making data unusable)
# Date: Date on which the replicate occurred
# Rec_time: timestep number
# Actual_time: Time of day at which the photograph was recorded
# Tidal_state: The tidal state within the experimental system at the time of the photograph (either 'Low' or 'High')
# Period: The photoperiod at the time of the photograph (either 'Day' or 'Night')
# Limpet_ID: Arbitrary identification number
# Limpet_col: The colour of the nail varnish identifier on the limpet's shell
# Limpet_size: The width of the limpet in cm
# Panel_type: The experimental panel onto which the limpet was placed (either: Flat, Ripple, cm_25 (2.5 cm high), or cm_5 (5 cm high))
# Temp_C: The temperature in the experimental system at the time of the photograph in degrees celsius (either water or air temperature depending on Tidal_state)
# X_coord: The x coordinate of the limpet at the time of the photograph
# Y_coord: The y coordinate of the limpet at the time of the photograph
# Distance_moved_cm: The distance moved, in cm, from the position at the previous timestep
# Moving: Binary response variable - If movement was registered above the threshold distance (0.01 cm), equal to 1, otherwise equal to 0
# Not_moving: Binary response variable - If movement was not registered or the distance moved was below the threshold distance (0.01 cm), equal to 1, otherwise equal to 0

# Remove the variables that are not necessary for analysis and re-classify the remaining variables
str(limpet)
limpet <- limpet[, c(2, 8, 10, 11, 15)] # Keep: Replciate, Limpet_ID, Limpet_size, Panel_type, Distance_moved_cm
limpet$Replicate <- as.factor(limpet$Replicate)
limpet$Limpet_ID <- as.character(limpet$Limpet_ID)
limpet$Panel_type <- as.factor(limpet$Panel_type)

# Calculate the gross distance moved by each limpet by summing across timesteps
limpet <- limpet %>% group_by(Replicate, Limpet_ID, Limpet_size, Panel_type) %>% 
  summarise(Gross_distance=sum(Distance_moved_cm))

# Check for normality in the data
ggdensity(limpet$Gross_distance, main="Density plot of Gross distance", xlab="Gross distance (cm)")
ggqqplot(limpet$Gross_distance, main="Q-Q plot of Gross distance")
shapiro.test(limpet$Gross_distance)

# Try data transformation (fourth root)
ggdensity(nthroot(limpet$Gross_distance, 4), main="Density plot of Gross distance", xlab="Gross distance (cm)")
ggqqplot(nthroot(limpet$Gross_distance, 4), main="Q-Q plot of Gross distance")
shapiro.test(nthroot(limpet$Gross_distance, 4))

# Keep the fourth root transformed gross distance data
limpet$root_gross <- nthroot(limpet$Gross_distance, 4)

# Construct a linear mixed effects (LME) model using REML

names(limpet)
# Replicate = factor = random
# Limpet_size = continuous = covariate
# Panel_type = factor = fixed

set.seed(1234)
mod1 <- lme(root_gross~Panel_type+Limpet_size, random=~1|Replicate, data=limpet, method="REML")
summary(mod1)
# AIC = 217.5 

# # Create a simpler model without the random factor to see if this improves the fit
# mod2 <- lm(root_gross~Panel_type+Limpet_size, data=limpet)
# # Compare
# anova(mod1, mod2) 
# # p = 3e-04 - Replicate has an effect on response, so keep the mixed effects model (mod1)

# # Diagnostic plots
# check_model(mod1)
# plot(mod1)
# qqnorm(resid(mod1))
# qqline(resid(mod1))
# hist(residuals(mod1))

summary(mod1)
# Fixed effects:
#                      Value   Std.Error  DF    t-value  p-value 
# (Intercept)      1.1199579  0.29584383  91   3.785639   0.0003 ***
# Tile_typecm_5   -0.3217242  0.16391273  91  -1.962777   0.0527 
# Tile_typeFlat   -0.3956888  0.16388827  91  -2.414381   0.0178 *
# Tile_typeRipple  0.3659555  0.17790737  91   2.057000   0.0425 *
# Limpet_size      0.3671122  0.07803337  91   4.704554   0.0000 ***

# Correlation: 
#                 (Intr) Tl_t5c Tl_tyF Tl_tyR
# Tile_typecm_5   -0.233                     
# Tile_typeFlat   -0.236  0.515              
# Tile_typeRipple -0.252  0.442  0.442       
# Limpet_size     -0.837 -0.044 -0.041  0.001

# Carry out a post-hoc test for the effects of panel type
set.seed(1234)
em <- emmeans(mod1, ~Panel_type)
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate    SE df t.ratio p.value
# Flat - Ripple     -0.762 0.181 91  -4.211  0.0003 ***
# Flat - cm_25       0.396 0.164 91   2.414  0.0815
# Flat - cm_5        0.074 0.161 91   0.458  0.9678
# Ripple - cm_25    -0.366 0.178 91  -2.057  0.1752
# Ripple - cm_5     -0.688 0.181 91  -3.802  0.0015 **
# cm_25 - cm_5       0.322 0.164 91   1.963  0.2098

# Generate a compact letter display for use in plots
cld(em)
# Flat     a
# Ripple   b
# 2.5 cm   ab
# 5 cm     a

# Summarise the data
limpet_avg <- limpet %>% group_by(Panel_type) %>% 
  summarise(mean=mean(Gross_distance), se=std.error(Gross_distance), median=median(Gross_distance))
limpet_avg
# Panel_type  mean    se  median
# Flat       31.0  11.8    8.33
# Ripple     74.1  13.1    64.0
# 2.5 cm     43.0  8.36    34.3
# 5 cm       31.4  8.13    18.9

# Re-order the panel types for the plot
limpet$Panel_type = factor(limpet$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

# Plot (Figure 4A)
ggplot(limpet, aes(x=Panel_type, y=Gross_distance))+
  geom_boxplot(fill="lightgrey")+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Gross distance (cm)")+
  xlab("Panel complexity")+
  theme_light()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), legend.position="none")

# Plot limpet size to show the effect of the covariate (Figure S5)
limpet$Panel_type = factor(limpet$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(limpet, aes(x=Limpet_size, y=Gross_distance))+
  geom_point()+
  facet_wrap(~Panel_type)+
  geom_smooth(method="lm", colour="black")+
  scale_x_continuous("Limpet shell size (cm)", breaks=seq(1.5, 5, 0.5))+
  scale_y_continuous("Gross distance (cm)")+
  stat_cor(label.y=220)+                      
  stat_regline_equation(label.y=200)+ 
  theme_light()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16))

# -------------------------------------------------------------------------
# ----------------------------- Net distance ------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# Net distance was calculated as the straight-line distance between the initial position of the limpet and its position after the 24-hour recording period.

# load the net distances data
limpet2 <- read.csv("./Net_movement.csv")

# Data description
head(limpet2)
# Replicate: The number of the experimental run (note: Replicate starts from 3 as replicates 1 & 2 were trial runs during which errors occurred with the equipment making data unusable)
# Date: Date on which the replicate occurred
# Limpet_ID: Arbitrary identification number
# Limpet_col: The colour of the nail varnish identifier on the limpet's shell
# Limpet_size: The width of the limpet in cm
# Panel_type: The experimental panel onto which the limpet was placed (either: Flat, Ripple, cm_25 (2.5 cm high), or cm_5 (5 cm high))
# Net_distance: The distance, in cm, between the initial position of the limpet and its position after 24-hours

# Remove the variables that are not necessary for analysis and re-classify the remaining variables
str(limpet2)
limpet2 <- limpet2[, c(1, 3, 5, 6, 7)] # Keep: Replciate, Limpet_ID, Limpet_size, Panel_type, Net_distance
limpet2$Replicate <- as.factor(limpet2$Replicate)
limpet2$Limpet_ID <- as.character(limpet2$Limpet_ID)
limpet2$Panel_type <- as.factor(limpet2$Panel_type)

# Check for normality in the data
ggdensity(limpet2$Net_distance, main="Density plot of Net distance", xlab="Net distance (cm)")
ggqqplot(limpet2$Net_distance, main="Q-Q plot of Net distance")
shapiro.test(limpet2$Net_distance)

# Try data transformation (fourth root)
ggdensity(nthroot(limpet2$Net_distance, 4), main="Density plot of Net distance", xlab="Net distance (cm)")
ggqqplot(nthroot(limpet2$Net_distance, 4), main="Q-Q plot of Net distance")
shapiro.test(nthroot(limpet2$Net_distance, 4))

# Keep the fourth-root transformed net distance data
limpet2$root_net <- nthroot(limpet2$Net_distance, 4)

# Construct a linear mixed effects (LME) model using REML

set.seed(1234)
mod1 <- lme(root_net~Panel_type+Limpet_size, random=~1|Replicate, data=limpet2, method="REML")
summary(mod1)
# AIC = 184.0

# There is no effect of limpet size (p = 0.24) so try simplifying the model by removing it
mod2 <- lme(root_net~Panel_type, random=~1|Replicate, data=limpet2, method="REML")
anova(mod1, mod2) # p = 0.14
# keep the simplified model

# # Model diagnostic plots
# check_model(mod2)
# plot(mod2)
# hist(residuals(mod2))

summary(mod2)
# Fixed effects:
#                      Value  Std.Error  DF    t-value  p-value
# (Intercept)      1.3146356  0.1002405  92  13.114818   0.0000 ***
# Tile_type5 cm   -0.1254045  0.1426516  92  -0.879096   0.3816
# Tile_typeFlat   -0.3461424  0.1426516  92  -2.426488   0.0172 *
# Tile_typeRipple  0.2647517  0.1527441  92   1.733302   0.0864

# Correlation: 
#                 (Intr) Tl_t5c Tl_tyF
# Tile_type5 cm   -0.699              
# Tile_typeFlat   -0.699  0.491       
# Tile_typeRipple -0.653  0.458  0.458

# Carry out a post-hoc test for the effects of panel type
set.seed(1234)
em <- emmeans(mod2, ~Panel_type)
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate     SE  df  t.ratio  p.value
# Flat - Ripple     -0.611  0.154  92   -3.968   0.0008 ***
# Flat - 2.5 cm      0.346  0.143  92    2.426   0.0792
# Flat - 5 cm        0.221  0.144  92    1.534   0.4214
# Ripple - 2.5 cm   -0.265  0.153  92   -1.733   0.3125
# Ripple - 5 cm     -0.390  0.154  92   -2.534   0.0613
# 2.5 cm - 5 cm      0.125  0.143  92    0.879   0.8157

# Generate a compact letter display for use in plots
cld(em)
# Flat    a
# Ripple  b
# 2.5 cm  ab
# 5 cm    ab

# Summarise the data
net_avg <- limpet2 %>% group_by(Panel_type) %>% 
  summarise(mean=mean(Net_distance), se=std.error(Net_distance), median=median(Net_distance))
net_avg
# Panel_type  mean    se  median
# Flat       2.33  0.69    1.25
# Ripple     9.64  1.65    8.68
# 2.5 cm     5.85  1.18    3.13
# 5 cm       4.65  1.15    2.28

# Re-order the panel types for the plot
limpet2$Panel_type = factor(limpet2$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

# Plot (Figure 4B)
ggplot(limpet2, aes(x=Panel_type, y=Net_distance))+
  geom_boxplot(fill="lightgrey")+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Net distance (cm)")+
  xlab("Panel complexity")+
  ylim(c(0, 25))+
  theme_light()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), legend.position="none")

# -------------------------------------------------------------------------
# ------------------------ Likelihood of movement -------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# As data describing the distance moved by limpets for each timestep of the 24-hour time lapse were highly zero-inflated, we examined the likelihood of movement under different conditions.

# Load the movement data
limpet <- read.csv("./Master_limpet_movement.csv")

# Remove the variables that are not necessary for analysis and re-classify the remaining variables
str(limpet)
limpet <- limpet[, c(1, 2, 5, 6, 8, 10, 11, 12, 15, 16, 17)] # Keep: Trial_number, replicate, Actual_time, Tidal_state, Limpet_ID, Limpet_size, Panel_type, Temp_C, Distance_moved_cm, Moving, Not_moving
limpet$Replicate <- as.factor(limpet$Replicate)
limpet$Tidal_state <- as.factor(limpet$Tidal_state)
limpet$Limpet_ID <- as.character(limpet$Limpet_ID)
limpet$Panel_type <- as.factor(limpet$Panel_type)

# Run a linear mixed effects (LME) model with temporal autocorrelation

# Create a binary response variable for movement/no movement
Y <- cbind(limpet$Moving, limpet$Not_moving)
# Logit transform the response variable
Y <- logit(Y)

# Replicate = Factor = Random
# Actual_time = Continuous = covariate ~ check for temporal autocorrelation
# Tidal_state = Factor = Fixed
# Limpet_ID = Character = Random, nested in Replicate
# Limpet_size = Continuous = covariate
# Panel_type = Factor = Fixed
# Temp_C = Continuous = covariate

# First run a model without any autocorrelation
set.seed(1234)
mod.lme <- lme(Y ~ Panel_type*Tidal_state + Limpet_size + Temp_C + Actual_time, 
               random=~1|Replicate/Limpet_ID, data=limpet, method='REML')
# Check for temporal autocorrelation by plotting an ACF and PACF plot
acf(residuals(mod.lme, type="normalized"), lag=40) 
# There is temporal autocorrelation - Movement has a high correlation with a lagged version of itself (geometric decline) - i.e., lags 1-26 are outside the boundary of the blue box.
pacf(residuals(mod.lme, type="normalized"), lag=40)
# PACF has significant spikes at lags 1-3 - for everything in the blue box, there is no evidence that it is different from zero 
summary(mod.lme)$AIC # AIC = 147431.6

# Try a model with ARMA (autoregressive moving average) for autocorrelation 
# Based on the PACF plot above, set p value as 3
mod.lme.ARMA <- update(mod.lme, correlation=corARMA(p=3))
# Check for temporal autocorrelation by plotting an ACF and PACF plot
acf(residuals(mod.lme.ARMA, type="normalized"), lag=40) 
pacf(residuals(mod.lme.ARMA, type="normalized"), lag=40)
# There is no longer temporal autocorrelation
summary(mod.lme.ARMA)$AIC # AIC = 107100.4   

# Test whether the interaction term is truly non-significant
em <- emmeans(mod.lme.ARMA, ~Panel_type:Tidal_state)
contrast(em, "pairwise", adjust="Tukey")
# There are some significant combinations, so best to keep the full model

# # Diagnostic plots: QQ plots for levels of random effects
# qqnorm(mod.lme.ARMA, ~ranef(., level=2))
# qqnorm(mod.lme.ARMA, ~ranef(., level=1))

# Summarise
summary(mod.lme.ARMA)
#                                     Value   Std.Error     DF    t-value  p-value
# (Intercept)                    -2.6483812   0.7116299  29037  -3.721571   0.0002
# Tile_typecm_5                  -0.1985450   0.2491155     91  -0.797000   0.4275
# Tile_typeFlat                   0.5852541   0.2506526     91   2.334922   0.0217 *
# Tile_typeRipple                 0.4444811   0.2699804     91   1.646346   0.1031
# Tidal_stateLow                 -0.6224752   0.1492735  29037  -4.170032   0.0000 ***
# Limpet_size                     0.4958748   0.1074756     91   4.613837   0.0000 ***
# Temp_C                         -0.0396245   0.0392605  29037  -1.009271   0.3129
# Actual_time                    -0.0000210   0.0000568  29037  -0.369049   0.7121
# Tile_typecm_5:Tidal_stateLow    0.3464536   0.2130287  29037   1.626323   0.1039
# Tile_typeFlat:Tidal_stateLow   -0.2604628   0.2175002  29037  -1.197529   0.2311
# Tile_typeRipple:Tidal_stateLow -0.0433842   0.2280257  29037  -0.190260   0.8491

# Carry out a post-hoc test for the interaction term
em <- emmeans(mod.lme.ARMA, ~Panel_type:Tidal_state)
contrast(em, "pairwise", adjust="Tukey")

# Generate a compact letter display for the plot
cld(em)
# Panel  Tide  
# Flat   Low   abc
# Ripple Low   abd
# 2.5 cm Low   a
# 5 cm   Low   ab
# Flat   High  de
# Ripple High  ce
# 2.5 cm High  bcde
# 5 cm   High  abc

rm(acf_cor, mod.lme, mod.lme.AR1)

# Construct input data over which to predict the probabilities
# Hold temperature at the mean value as it was non-significant
# Hold time constant as it was non-significant
newdata <- with(limpet, data.frame(Panel_type=limpet$Panel_type,
                                   Tidal_state=limpet$Tidal_state, Limpet_size=limpet$Limpet_size,
                                   Temp_C=mean(limpet$Temp_C), Actual_time=1200))
newdata <- newdata[!duplicated(newdata),] # Remove duplicates

# Expand data to include all possible combinations of factors, and expand the limpet size to include all possible values (increasing in increments of 0.1) because it was a significant covariate
newdata <- newdata %>% expand(Panel_type, Tidal_state,  Actual_time, Temp_C, full_seq(Limpet_size, 0.1))
colnames(newdata) <- c("Panel_type", "Tidal_state", "Actual_time", "Temp_C", "Limpet_size")

# Extract the predicted probabilities 
newdata$prob <- predict(mod.lme.ARMA, newdata=newdata, level=0)
# Back transform the logit predictions (inverse logit)
newdata$prob <- exp(newdata$prob)/(1+exp(newdata$prob))
# # Check values seem reasonable
# min(newdata$prob) # 0.04
# max(newdata$prob) # 0.44

# Remove temperature as there was no significant effect
newdata$Temp_C <- NULL
newdata$Actual_time <- NULL
newdata$Panel_type = factor(newdata$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE) #Re-order panel type

# Plot the effect of limpet size on likelihood of movement
limp_size <- newdata %>% group_by(Limpet_size, Panel_type, Tidal_state) %>% 
  summarise(n = mean(prob), se = std.error(prob))
ggplot(limp_size, aes(x=Limpet_size, y=n, col=Tidal_state))+
  geom_line(linewidth=1)+
  facet_wrap(~Panel_type)+
  ylab("Likelihood of movement")+
  xlab("Limpet shell length (cm)")+
  ylim(c(0, 0.5))+
  cc_theme()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16))

# The effect of tidal state on likelihood of movement
summary_tide <- newdata %>% group_by(Tidal_state) %>% 
  summarise(n=mean(prob), se=std.error(prob))
summary_tide
# Immersed  0.209  0.00795
# Emersed   0.125  0.00469
((0.209-0.125)/0.209)*100  # Movement 40% more likely during immersion than emersion

# The effect of panel type on likelihood of movement
summary_panel <- newdata %>% group_by(Panel_type) %>%  
  summarise(n=mean(prob), se=std.error(prob))
summary_panel
# Flat    0.203  0.0124
# Ripple  0.194  0.0109
# 2.5 cm  0.138  0.00819
# 5 cm    0.132  0.00691

# Panel type & tidal state
summary <- newdata %>% group_by(Panel_type, Tidal_state) %>% 
  summarise(n=mean(prob), se=std.error(prob), median=median(prob))
summary 
# Treatment    mean   se
# Flat High    0.269  0.0161
# Ripple High  0.244  0.0151
# 2.5 cm High  0.174  0.0118
# 5 cm High    0.148  0.0104
# Flat Low     0.136  0.0097
# Ripple Low   0.145  0.0102
# 2.5 cm Low   0.103  0.0076
# 5 cm Low     0.117  0.0085

# % difference for significant contrasts
((0.174-0.103)/0.174)*100 # 41% more likely on 2.5 cm at high tide than low tide
((0.269-0.148)/0.269)*100 # 45% more likely on flat than 5 cm at high tide
((0.269-0.136)/0.269)*100 # 49% more likely on flat panels at high tide than low tide
((0.244-0.145)/0.244)*100 # 41% more likely on ripple panels at high tide than low tide

# Likelihood of movement for each panel at each tidal state, across all limpet sizes
# (Figure 5)
ggplot(newdata, aes(x=Tidal_state, y=prob, fill=Panel_type))+
  geom_boxplot()+
  ylab("Likelihood of movement")+
  xlab("Immersed/Emersed")+
  ylim(c(0, 0.5))+
  theme_light()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16))

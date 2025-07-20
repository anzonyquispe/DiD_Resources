rm(list = ls())
# Load Libraries
library(dplyr)
library(fixest)
library(did)
library(fect)
library(panelView)
library(PanelMatch)
library(ggplot2)
library(bacondecomp)
# library(didimputation)
library(doParallel)
library(DIDmultiplegtDYN)
library(data.table)
library(ggfixest)

# Data
e <- new.env()
load(url("https://raw.githubusercontent.com/xuyiqing/fect/master/data/fect.RData"))
data <- hh2019

# Panel View
panelview(nat_rate_ord ~ indirect, data = data, index = c("bfs","year"),
          xlab = "Year", ylab = "Unit", display.all = T,
          gridOff = TRUE, by.timing = TRUE)



# TWFE
model.twfe.0 <- feols(nat_rate_ord~indirect|bfs+year,
                      data=data, cluster = "bfs")
print(model.twfe.0)

# Bacon Decomposition
data.complete <- data[which(!is.na(data$nat_rate_ord)),]



df_bacon <- bacon(nat_rate_ord~indirect,
                  data = data.complete,
                  id_var = "bfs",
                  time_var = "year")
ggplot(df_bacon) +
  aes(x = weight, y = estimate, shape = factor(type), color = factor(type)) +
  labs(x = "Weight", y = "Estimate", shape = "Type", color = 'Type') +
  geom_point()
print(aggregate(df_bacon$estimate * df_bacon$weight,
                list(df_bacon$type), FUN=sum))

# Remove always treated units
df <- as.data.frame(data %>%
                      group_by(bfs) %>%
                      mutate(treatment_mean = mean(indirect,na.rm = TRUE)))
df.use <- df[which(df$treatment_mean<1),]


# Panel View
panelview(nat_rate_ord ~ indirect, data = df.use, index = c("bfs","year"),
          xlab = "Year", ylab = "Unit", display.all = T,
          gridOff = TRUE, by.timing = TRUE)


model.twfe.1 <- feols(nat_rate_ord~indirect|bfs+year,
                      data=df.use, cluster = "bfs")
print(model.twfe.1)

# Cohort & Event Study
df.use <- fect::get.cohort(df.use, D = "indirect", index=c("bfs","year"))
head(df.use[,-5],19)

# Dynamic TWFE
df.twfe <- df.use
df.twfe$treat <- as.numeric(df.twfe$treatment_mean>0)
table(df.twfe$treat)
df.twfe <- as.data.table(df.twfe)
fyr <- df.twfe[indirect ==1 , .(yeartr = min(year , na.rm = TRUE )), by = bfs  ]
df.twfe2 <- df.twfe %>% 
  left_join(fyr, "bfs")
df.twfe2$Time_to_Treatment <- 0
df.twfe2[treat == 1, ]$Time_to_Treatment <- df.twfe2[treat == 1, ]$year - df.twfe2[treat == 1, ]$yeartr
table(df.twfe2$Time_to_Treatment)
twfe.est <- feols(nat_rate_ord ~ i(Time_to_Treatment, treat, ref = -1)| bfs + year,
                  data = df.twfe2, cluster = "bfs")
ggiplot(twfe.est)


# Interaction Weighted
df.sa <- df.twfe2
df.sa[which(is.na(df.sa$yeartr)),"yeartr"] <- 1000
model.sa.1 <- feols(nat_rate_ord~sunab(yeartr,year)|bfs+year,
                    data = df.sa, cluster = "bfs")
summary(model.sa.1,agg = "ATT")

ggiplot(model.sa.1)



# CSDID
df.cs <- df.twfe2
df.cs$FirstTreat<- df.cs$yeartr
df.cs[which(is.na(df.cs$FirstTreat)),"FirstTreat"] <- 0
cs.est.1 <- att_gt(yname = "nat_rate_ord",
                   gname = "FirstTreat",
                   idname = "bfs",
                   tname = "year",
                   xformla = ~1,
                   allow_unbalanced_panel = TRUE,
                   data = df.cs,
                   cluster = "bfs", 
                   est_method = "dr", 
                   base_period = "universal")

cs.est.att.1 <- aggte(cs.est.1, type = "simple", na.rm=T, bstrap = F)
print(cs.est.att.1)
cs.att.1 <- aggte(cs.est.1, type = "dynamic",
                  bstrap=FALSE, cband=FALSE, na.rm=T)
print(cs.att.1)
ggdid(cs.att.1)



# Stacked DID
df.st <- NULL
df.use <- df.cs
df.use$Cohort <- df.use$FirstTreat
target.cohorts <- setdiff(unique(df.use$Cohort),0)
sort(target.cohorts)
k <- 1
for(cohort in target.cohorts){
  df.sub <- df.use[which(df.use$Cohort%in%c(cohort,0)),]
  df.sub$stack <- k
  df.sub$event_time <- df.sub$year - cohort
  df.sub$cohort <-cohort
  df.st <- rbind(df.st,df.sub)
  k <- k + 1
}
model.st <- feols(nat_rate_ord~indirect|bfs^stack  +year^stack ,
                data=df.st, cluster = ~bfs^stack)
df.st[stack==1,]
st.est <- feols(nat_rate_ord ~
                  i(event_time, treat, ref = -1)| bfs^stack + year^stack ,
                data = df.st,cluster = ~bfs^stack)

ggiplot(st.est)


# DIDmultiple
didm.results <- did_multiplegt_dyn(
  df = df.use,
  outcome = "nat_rate_ord",
  group = "bfs",
  controls = NULL,
  time = "year",
  treatment = "indirect",
  effects = 12,
  placebo = 9,
  cluster = "bfs"
)
print(didm.results)





df.pm <- df.use
# we need to convert the unit and time indicator to integer
df.pm[,"bfs"] <- as.integer(as.factor(df.pm[,"bfs"]))
df.pm[,"year"] <- as.integer(as.factor(df.pm[,"year"]))
df.pm <- df.pm[,c("bfs","year","nat_rate_ord","indirect")]

# Pre-processes and balances panel data
df.pm <- PanelData(panel.data = df.pm,
                   unit.id = "bfs",
                   time.id = "year",
                   treatment = "indirect",
                   outcome = "nat_rate_ord")

PM.results <- PanelMatch(lag=3, 
                         refinement.method = "none", 
                         panel.data = df.pm, 
                         qoi = "att", 
                         lead = c(0:3), 
                         match.missing = TRUE)

## For pre-treatment dynamic effects
PM.results.placebo <- PanelMatch(lag=3, 
                                 refinement.method = "none", 
                                 panel.data = df.pm, 
                                 qoi = "att", 
                                 lead = c(0:3), 
                                 match.missing = TRUE,
                                 placebo.test = TRUE)

# ATT
PE.results.pool <- PanelEstimate(PM.results, panel.data = df.pm, pooled = TRUE)
summary(PE.results.pool)


# Dynamic Treatment Effects
PE.results <- PanelEstimate(PM.results, panel.data = df.pm)
PE.results.placebo <- placebo_test(PM.results.placebo, panel.data = df.pm, plot = F)

# obtain lead and lag (placebo) estimates
est_lead <- as.vector(PE.results$estimate)
est_lag <- as.vector(PE.results.placebo$estimates)
sd_lead <- apply(PE.results$bootstrapped.estimates,2,sd)
sd_lag <- apply(PE.results.placebo$bootstrapped.estimates,2,sd)
coef <- c(est_lag, 0, est_lead)
sd <- c(sd_lag, 0, sd_lead)
pm.output <- cbind.data.frame(ATT=coef, se=sd, t=c(-2:4)) %>% as.data.table()
pm.output[, `:=`(
  ci_lower = ATT - 1.96 * se,
  ci_upper = ATT + 1.96 * se
)]
# plot
ggplot(pm.output, aes(x = t, y = ATT)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray60") +
  labs(
    x = "Event Time",
    y = "ATT",
    title = "Event Study Plot with 95% Confidence Intervals"
  ) +
  theme_minimal()





out.fect <- fect(nat_rate_ord~indirect, data = df, 
                 index = c("bfs","year"),
                 method = 'fe', se = TRUE )
print(out.fect$est.avg)
plot(out.fect)


################################################################################
########################## Reversal Treatments #################################
################################################################################

library(fect)
data(fect)
data <- gs2020
data$cycle <- as.integer(as.numeric(data$cycle/2))


y <- "general_sharetotal_A_all"
d <- "cand_A_all"
unit <- "district_final"
time <- "cycle"
controls <- c("cand_H_all", "cand_B_all")
index <- c("district_final", "cycle")

panelview(Y=y, D=d, X=controls, index = index, data = data, 
          xlab = "Time Period", ylab = "Unit", gridOff = TRUE, 
          by.timing = TRUE, cex.legend=5, cex.axis= 5, 
          cex.main = 10, cex.lab = 5)

model.twfe <- feols(general_sharetotal_A_all ~ cand_A_all + 
                      cand_H_all + cand_B_all | district_final + cycle,
                    data=data, cluster = "district_final") 
summary(model.twfe)



data_cohort <- get.cohort(data, index = index, D=d,start0 = TRUE)
# Generate a dummy variable treat
data_cohort$treat <- 0
data_cohort[which(data_cohort$Cohort!='Control'),'treat'] <- 1
data_cohort[which(is.na(data_cohort$Time_to_Treatment)), "treat"] <- 0

# remove observations that starts with treated status
remove <- intersect(which(is.na(data_cohort$Time_to_Treatment)),
                    which(data_cohort[,d]==1)) 
if(length(remove)>0){data_cohort <- data_cohort[-remove,]}

# replace missingness in Time_to_Treatment with an arbitrary number
data_cohort[which(is.na(data_cohort$Time_to_Treatment)), "Time_to_Treatment"] <- 999 

twfe.est <- feols(general_sharetotal_A_all ~ 
                    i(Time_to_Treatment, treat, ref = -1) + 
                    cand_H_all +cand_B_all | district_final + cycle,  
                  data = data_cohort, cluster = "district_final")



twfe.output <- as.data.frame(twfe.est$coeftable[c(1:25),])
twfe.output$Time <- c(c(-16:-2),c(0:9)) + 1 

# plot
p.twfe <- esplot(twfe.output,Period = 'Time',Estimate = 'Estimate',
                 SE = 'Std. Error', xlim = c(-15,1))
p.twfe




df.pm <- data_cohort
# we need to convert the unit and time indicator to integer
df.pm[,"district_final"] <- as.integer(as.factor(df.pm[,"district_final"]))
df.pm[,"cycle"] <- as.integer(as.factor(df.pm[,"cycle"]))
df.pm <- df.pm[,c("district_final","cycle","cand_A_all", 
                  "general_sharetotal_A_all")]

# Pre-processes and balances panel data
df.pm <- PanelData(panel.data = df.pm,
                   unit.id = "district_final",
                   time.id = "cycle",
                   treatment = "cand_A_all",
                   outcome = "general_sharetotal_A_all")

PM.results <- PanelMatch(lag=4, 
                         refinement.method = "none", 
                         panel.data = df.pm, 
                         qoi = "att", 
                         lead = 0, 
                         match.missing = TRUE)

## For pre-treatment dynamic effects
PM.results.placebo <- PanelMatch(lag=4, 
                                 refinement.method = "none", 
                                 panel.data = df.pm, 
                                 qoi = "att", 
                                 lead = 0, 
                                 match.missing = TRUE,
                                 placebo.test = TRUE)

PE.results.pool <- PanelEstimate(PM.results, panel.data = df.pm, pooled = TRUE)
summary(PE.results.pool)
#>      estimate  std.error       2.5%     97.5%
#> [1,] 0.124243 0.02198589 0.08458918 0.1726429
#> 
#> 
#> # Dynamic Treatment Effects
PE.results <- PanelEstimate(PM.results, panel.data = df.pm)
PE.results.placebo <- placebo_test(PM.results.placebo, panel.data = df.pm,
                                   plot = FALSE)

est_lead <- as.vector(PE.results$estimate)
est_lag <- as.vector(PE.results.placebo$estimates)
sd_lead <- apply(PE.results$bootstrapped.estimates,2,sd)
sd_lag <- apply(PE.results.placebo$bootstrapped.estimates,2,sd)
coef <- c(est_lag, 0, est_lead)
sd <- c(sd_lag, 0, sd_lead)
pm.output <- cbind.data.frame(ATT=coef, se=sd, t=c(-3:1))

# plot
p.pm <- esplot(data = pm.output,Period = 't',
               Estimate = 'ATT',SE = 'se')
p.pm



didm.results <- did_multiplegt_dyn(
  df = df.pm,
  outcome = "general_sharetotal_A_all",
  group = "district_final",
  controls = NULL,
  time = "cycle",
  treatment = "cand_A_all",
  effects = 2,
  placebo = 3,
  cluster = "district_final"
)
print(didm.results)






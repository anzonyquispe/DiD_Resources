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






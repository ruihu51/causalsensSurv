rm(list=ls())

library(haven)

aarp_data1 <- read_sas("../AARP-data/Chuck Matthews 11Oct2023/chuck_matthews_09oct2023.sas7bdat")
head(aarp_data1)
dim(aarp_data1)[1]

##################
# data exclusion
##################
table(aarp_data1$FL_RF_PROXY) # 1 - 1660
aarp_data2 <- aarp_data1[aarp_data1$FL_RF_PROXY !=1,]
table(aarp_data2$RF_PHYS_MODVIG_CURR) # 9 - 2157
aarp_data <- aarp_data2[aarp_data2$RF_PHYS_MODVIG_CURR !=9,]
dim(aarp_data)[1] 
# TODO: didn't match with 159964

# 159937/163781
# 163781-1660-2157
# 163781-1660-159937
# 159964/163781

######################
# data pre processing
#####################

# use functions in data_application_pipeline.R to pre-process the variables
# table 1
table(aarp_data$BMI_CUR1)
summary(aarp_data$BMI_CUR)
# TODO: 2728 NA
table(aarp_data$RF_PHYS_MODVIG_CURR1)
# 0     1     2 
# 22526 56613 80825
table(aarp_data$HEI2015_TOTAL_SCORE1)
# 0     1     2     3 
# 39991 39991 39991 39991
table(aarp_data$MPED_A_BEV_NOFOOD1)
# 0      1 
# 27897 132067

# 
# generate BMI
summary(aarp_data$HT_CUR)
summary(aarp_data$WT_CUR)
gen_bmi <- function(x, y){
if (is.na(x) | (is.na(y))) {
    return(NA)
  } else {
    return(y/(x^2))
  }
}
aarp_data$bmi <- mapply(gen_bmi, aarp_data$HT_CUR, aarp_data$WT_CUR)
summary(aarp_data$bmi)


# table 2
# The distribution of all the other baseline confounders looks reasonable.
# Reasonable means the percentages falls in the range shown in Table 2.
table(aarp_data$SEX)/159964
table(aarp_data$ENTRY_AGE1)/159964
table(aarp_data$RACEI1)/159964
table(aarp_data$EDUCM1)/159964
table(aarp_data$HEALTH1)/159964
table(aarp_data$SMOKE_DOSE)/159964
table(aarp_data$SMOKE_QUIT1)/159964

# with a mean (SD) follow-up of 18.9 (6.3) years
summary(aarp_data$PERSONYRS)
mean(aarp_data$PERSONYRS) # 18.97
sd(aarp_data$PERSONYRS) # 6.24


package_name_lme <- "lme4"
package_name_car <- "car"
package_name_emmeans <- "emmeans"


if (!requireNamespace(package_name_lme, quietly = TRUE)) {
  install.packages(package_name_lme)  
} else {
  message(paste(package_name_lme, "is already installed."))
}

if (!requireNamespace(package_name_car, quietly = TRUE)) {
 install.packages(package_name_car)  
} else {
 message(paste(package_name_car, "is already installed."))
}

if (!requireNamespace(package_name_emmeans, quietly = TRUE)) {
  install.packages(package_name_emmeans) 
} else {
  message(paste(package_name_emmeans, "is already installed."))
}

library(lme4)
library(car)
library(emmeans)


path <- "Network_size_table.csv"
data <- read.csv(path)

data$Stage <- factor(data$Stage)
data$Gender <- factor(data$Gender)
data$participant_id <- factor(data$participant_id)
data$epoch <- factor(data$epoch)


# 
# 'NS': network size
# 'Stage': factor
# 'Age'
# 'Gender'
# 'Education'
# 'participant_id'

# fit
model <- lmer(NS ~ Stage + Age + Gender + Education + (1|participant_id) + (1|epoch), data = data)



car::Anova(LMM.fit)

emout<-emmeans(model, pairwise ~ Stage,pbkrtest.limit = 3046,adjust = "fdr")
emout

# Load 'lavaan' package
library("lavaan")

# Import data as .csv
data <- read.csv("data.csv", header=F)
colnames(data) <- c("student_ID", "school_ID", "enrollment", 
"read", "math", "science", "L1_sample_size", "sampling_ratio")

# Define doubly latent model. Note: NA*variable allows the variable
# to be freely estimated. Here, the first factor loading is not 
# constrained to 1 by default, but freely estimated. c1*variable 
# allows for a constraint to be placed on that variable. The factor 
# loadings for the 3 observed variables are constrained to be equal 
# across L1 and L2. The L1 factor variance is fixed to 1, while the 
# L2 factor variance is freely estimated. Residual variances of the 
# 3 observed variables are estimated at both L1 and L2 (var ~~ var). 
model <- '
level: 1
Factor_w =~ NA*read + c1*read + c2*math + c3*science

read ~~ read
math ~~ math
science ~~ science

Factor_w ~~ 1*Factor_w 

level: 2
Factor_b =~ NA*read + c1*read + c2*math + c3*science

read ~~ read
math ~~ math
science ~~ science

Factor_b ~~ Factor_b
'

# Estimate the model using maximum likelihood with robust standard errors
fit <- sem(model = model, data = data, cluster = "school_ID", se = "robust.huber.white")

# Examine fit statistics
fitMeasures(fit, fit.measures = "all")

# Examine parameter estimates
summary(fit)


#############
# Using Random Forest algorithm to train model
#############

# Set working directory
# path <- "/Users/baoyinyuan/Nutstore Files/LCN82/Script_Bsite_firstRevision_20231001"
# setwd(path)

#---------------
# Load required packages
#---------------
packages <- c("data.table",
              "mlr",      # Interface to a large number of classification and regression techniques
              "mlrMBO",
              "readxl",        #<-- Load data in excel
              "vip",
              "pdp",
              'magrittr', # pipe operator
              "ggplot2", # for plotting
              "rpart",
              "caret",   # for modeling
              "dplyr",   # for data manipulation and joining
              "gridExtra",
              "cowplot",
              "gbm",
              "h2o",
              "xgboost", # for building XGBoost model
              "cowplot",  # for combining multiple plots
              "Ckmeans.1d.dp", # for xgboost.ggplot.importance
              "parallel",
              "parallelMap"
)
for(pkg in packages){
  library(pkg, character.only = TRUE)
}


#---------------
# 1. Import data
#---------------

# import data of changeable experiment condition 
df_expe <- read_excel("experiment_data.xlsx", sheet = "PredictOH")%>%
  as.data.frame()
dim(df_expe)
df_expe[is.na(df_expe)] <- 0
df_expe %>% head(3)
# import data of element attribute
df_elem <- read.csv("elements_data.csv") %>% 
  as.data.frame()
df_elem %>% head(5)

# specify the input sample
# Can be used as the input interface
(df_expe_colname <- df_expe %>% colnames()) 
AA1BB1O3_sample <- data.frame(matrix(ncol = length(df_expe_colname), nrow = 0,
                                     dimnames = list(NULL, df_expe_colname)))
AA1BB1O3_sample <- df_expe
AA1BB1O3_sample %>% head(3)
dim(AA1BB1O3_sample)

#-----------------------------
# 2. Descriptor customization
#-----------------------------

# Customize the descriptors based on the oxide composition and element information
source("initial_descriptor_customize.R")
df_sample_descriptor <- descriptor_customize(AA1BB1O3_sample) 
df_sample <- cbind(AA1BB1O3_sample[c("Composition_AA1BB1O3", "HPC")], df_sample_descriptor)

df_sample %>% colnames()
dim(df_sample)
df_sample %>% head(3)

## Correlation analysis
df_sample_x <- df_sample[, 3:ncol(df_sample)] # predictor variables
df_sample_y <- select(df_sample, 2) # response variable
df_sample_y$HPC <- df_sample_y$HPC

# Load data
df_sample %>% head(3)
data_NoName <- df_sample[, -1]
data_NoName %>% head(3)
data_NoName<- within(data_NoName, rm("fraction_A1", "fraction_B1"))
data_NoName %>% head(3)
idx_zeroHPC <- which(data_NoName$HPC == 0)
data_all <- data_NoName[-idx_zeroHPC,]
# data_all <- data_NoName
dim(data_all)

##-------------
# 3.1 Hyper-parameter tuning using mlr package#
##-------------

# create tasks
traintask <- makeRegrTask(data = data_all, target = "HPC")

set.seed(2023)
# create learner
rf.lrn <- makeLearner("regr.randomForest", predict.type = "response")
rf.lrn$par.vals <- list(importance = TRUE) 
rf.lrn <- makePreprocWrapperCaret(rf.lrn, ppc.scale = TRUE, ppc.center = TRUE)

# set parameter space
# getParamSet("regr.randomForest")
getParamSet("regr.randomForest")
params <- makeParamSet(makeIntegerParam("ntree", lower = 1, upper = 600), # nbr of trees to grow
                       makeIntegerParam("mtry", lower = 1, upper = 15), # it refers to how many variables we should select at node split
                       makeIntegerParam("nodesize", lower = 1L, upper = 50L) # it refers to how many observations we want in the terminal nodes
                       )

# set resampling strategy
r_desc <- makeResampleDesc("RepCV", folds = 5L, reps = 5)
r_inst <- makeResampleInstance(r_desc, task = traintask)

# search strategy
mbo.ctrl <- mlrMBO::makeMBOControl()
mbo.ctrl <- mlrMBO::setMBOControlTermination(mbo.ctrl, iters = 50)
ctrl <- mlr:::makeTuneControlMBO(mbo.control = mbo.ctrl)

# set parallel backend
library(parallel)
library(parallelMap)
(numCores <- detectCores())
parallelStartSocket(cpus = numCores - 1)

# parameter tuning
start_time <- Sys.time()
rf.par_tuning <- tuneParams(learner = rf.lrn, 
                         task = traintask, 
                         resampling = r_inst,
                         measures = mlr::rmse,
                         par.set = params,
                         control = ctrl,
                         show.info = TRUE)
end_time <- Sys.time()
(time_elapsed <- end_time - start_time)
parallelStop()
rf.par_tuning

# Show up all test performance via RMSE
rf.lrn_testPerf <- makeLearner("regr.randomForest", par.vals = rf.par_tuning$x)
rf.lrn_testPerf <- makePreprocWrapperCaret(rf.lrn_testPerf, ppc.scale = TRUE, ppc.center = TRUE)

rf.error_resample = resample(
  learner = rf.lrn_testPerf, 
  task = traintask, 
  resampling = r_inst, 
  measures = list(mlr::rmse, mlr::mae, mlr::rsq),
  models = TRUE,
  show.info = TRUE,
  keep.pred = TRUE)
rf.error_resample
rf.error_resample$aggr
all.equal(as.numeric(rf.error_resample$aggr[1]), as.numeric(rf.par_tuning$y)) 


##-------------
# 3.2 Predictive performance using rmse, mae, r2
##-------------
## mean and sd for rmse, mae and R-squared
(rmse_mean <- mean(rf.error_resample$measures.test$rmse))
(rmse_sd <- sd(rf.error_resample$measures.test$rmse))

(mae_mean <- mean(rf.error_resample$measures.test$mae))
(mae_sd <- sd(rf.error_resample$measures.test$mae))

(R2_mean <- mean(rf.error_resample$measures.test$rsq))
(R2_sd <- sd(rf.error_resample$measures.test$rsq))

#
rf.error_table <- data.frame(matrix(ncol = 2, nrow = 3))
colnames(rf.error_table) <- c("Mean", "SD")
rownames(rf.error_table) <- c("RMSE", "MAE", "R-squared")
rf.error_table$Mean <- c(rmse_mean, mae_mean, R2_mean)
rf.error_table$SD <- c(rmse_sd, mae_sd, R2_sd)
rf.error_table %<>% mutate_if(is.numeric, round, digits=3)
rf.error_table

## Visualize the true values and predicted values
# true values & predicted values
rf.pred = getRRPredictions(rf.error_resample)
rf.true_pred <- data.frame(Truth = getPredictionTruth(rf.pred),
                           Pred = pmax(getPredictionResponse(rf.pred), 0),
                           Repeats = as.factor((rf.pred$data$iter-1) %/% 5 + 1))


# calculate the 95% confidence interval of predicted values versus true values
dim(rf.true_pred)
(sample_num <- length(unique(rf.true_pred$Truth)))
rf.true_pred_ci <- data.frame(matrix(ncol = 5, nrow = sample_num))
colnames(rf.true_pred_ci) <- c("True", "Mean", "Upper","Median", "Lower")
for(x in seq(sample_num)){
  # x = 1
  x_true = unique(rf.true_pred$Truth)[x]
  x_index = which(rf.true_pred$Truth == x_true)
  x_pred = rf.true_pred[x_index, ]$Pred
  x_mean = mean(x_pred)
  rf.true_pred_ci[x,] = c(x_true, x_mean, quantile(x_pred, probs = c(0.975, 0.5, 0.025)))
}


## ggplot with confidence
ggplot(rf.true_pred_ci, aes(x = True, y = Median)) +
  geom_point(alpha = 0.9) +
  geom_errorbar(aes(ymin = Lower,  ymax = Upper), width = 0, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.8, alpha = 0.9) +
  labs(x = "True values", y = "Predicted values by random forest (95% CI)") + 
  scale_x_continuous(limits = c(0, 0.4), expand = c(0.02, 0)) +
  theme_bw(base_size=12) +
  theme(plot.margin = unit(c(3, 1, 1, 1), "lines"),
        axis.text.x=element_text(size=10, angle = 0, vjust= 0.5, hjust = 0.5, colour = "black"),
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title=element_text(size=12, colour = "black", face = "bold"),
        plot.title = element_text(size = 11, face="bold"),
        legend.background = element_rect(fill="transparent",
                                         linewidth = 0.3, linetype = "solid",
                                         colour = "black"),
        legend.box.background = element_rect(colour = "black"),
        legend.title=element_text(size= 8, colour = "black"),
        legend.text=element_text(size= 8, colour = "black"),
        legend.position = c("none"),#"left",#c(.75, .95),
        legend.direction="horizontal",
        legend.justification = c("right", "top"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.8)
  ) +
  annotate(geom = "text", x = 0.32, y =0.38, size= 5, col = "red",
           label = paste0("1:1 line"), hjust = "left") +
  annotate(geom = "text", x = 0.05, y = 0.44, size= 4, col = "black", fontface = "bold",
           label = paste0("RMSE = ", rf.error_table[1,1], "\U00B1", rf.error_table[1,2],";")) +
  annotate(geom = "text", x = 0.185, y = 0.44, size= 4, col = "black", fontface = "bold",
           label = paste0("MAE = ", rf.error_table[2,1],"\U00B1", rf.error_table[2,2], ";")) +
  annotate(geom = "text", x = 0.33, y = 0.44, size= 4, col = "black", fontface = "bold",
           label = paste0("R-squared = ", rf.error_table[3,1], "\U00B1", rf.error_table[3,2])) +
  coord_cartesian(xlim = c(0, 0.4), ylim = c(0, 0.4), expand = TRUE, clip="off") -> rf.true_pred_ci_plt

rf.true_pred_ci_plt


#
ggsave(filename = "rf.true_pred_ci_plt_rep10cv10_Asite_20240907.tiff",
       plot = rf.true_pred_ci_plt,
       width = 16,  # <=19.05cm
       height = 12, # <=22.225cm
       units= "cm",
       dpi= 300,
       compression = "lzw")  







#















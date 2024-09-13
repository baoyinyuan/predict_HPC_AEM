#############
# Using XGBoost algorithm to train model
#############

# Set working directory
path <- "/Users/baoyinyuan/Nutstore Files/LCN82/Script_Asite"
setwd(path)
getwd()

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

# ######
# # import data of changable experiment condition 
# df_expe2 <- read_excel("experiment_data_2.xlsx", sheet = "PredictOH")%>%
#   as.data.frame()
# dim(df_expe2)
# df_expe2$eletrode <- as.factor(df_expe2$eletrode)
# df_expe2 %>% head(2)
# 
# ggplot(df_expe2, aes(eletrode, HPC)) +
#   geom_violin(aes(fill = eletrode)) +
#   scale_x_discrete(labels = c("Electrolyte","Eletrode")) +
#   labs(x = "Sample type", y = "HPC") +
#   theme(legend.position = "none")
#   
# 
# df_expe_eletrode <- df_expe2[df_expe2$eletrode == 1,]
# df_expe_eletrode
# dim(df_expe_eletrode)
# df_expe %>% head(3)
# specify the input sample using only eletrode
# (df_expe_colname <- df_expe_eletrode %>% colnames()) 
# 
# AA1BB1O3_sample <- data.frame(matrix(ncol = length(df_expe_colname), nrow = 0,
#                                      dimnames = list(NULL, df_expe_colname)))
# AA1BB1O3_sample <- df_expe_eletrode
# AA1BB1O3_sample %>% head(3)
# dim(AA1BB1O3_sample)

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
dim(data_NoName)
(idx_zeroHPC <- which(data_NoName$HPC == 0))
data_all <- data_NoName[-idx_zeroHPC, ]
# data_all <- data_NoName
dim(data_all)
#-----------------------------
#3. Regression model training 
#-----------------------------

##-------------
# 3.1 Hyper-parameter tuning using mlr package#
##-------------

# create tasks
traintask <- makeRegrTask(data = data_all, target = "HPC")
# testtask <- makeRegrTask(data = data_stdd_test, target = "HPC")

# getTaskData(testtask) # extract data in defined task
# getTaskTargets(testtask) # extract target variable in defined task

set.seed(2023)

# create learner
xgb.lrn_tune <- makeLearner("regr.xgboost", predict.type = "response")
xgb.lrn_tune$par.vals <- list(booster = "gbtree",
                     objective = "reg:squarederror")
                     # eval_metric = "rmse", 
                     # nrounds = 100L # it controls the maximum nbr of iterations
                     # eta = 0.3) # it controls the learning rate
xgb.lrn_tune <- makePreprocWrapperCaret(xgb.lrn_tune, ppc.scale = TRUE, ppc.center = TRUE)

# set parameter space
xgb.params <- makeParamSet(makeIntegerParam("nrounds", lower = 100, upper = 500),
                          makeNumericParam("eta", lower = 0.1, upper = 1L), # it controls regularization
                           makeNumericParam("gamma", lower = 0, upper = 1L), # it controls regularization
                           makeIntegerParam("max_depth", lower = 3L, upper = 10L), # it controls the depth of the tree
                           makeNumericParam("min_child_weight", lower = 1L, upper = 10L), # the minimum nbr of instances required in child node
                           makeNumericParam("lambda", lower = 0, upper = 1), # L2 regularization term on weights
                           makeNumericParam("alpha", lower = 0, upper = 1) # L1 regularization term on weights
                           # makeNumericParam("subsample", lower = 0, upper = 1), # It controls the nbr of samples supplied to a tree
                           # makeNumericParam("colsample_bytree", lower = 0, upper = 1) # It controls the nbr of features supplied to a tree
                           )

# set resampling strategy
r_desc <- makeResampleDesc("RepCV", folds = 10L, reps = 10L)
r_inst <- makeResampleInstance(r_desc, task = traintask)

## search strategy
# ctrl <- makeTuneControlRandom(maxit = 1000L) # random search
mbo.ctrl <- makeMBOControl()
mbo.ctrl <- setMBOControlTermination(mbo.ctrl, iters = 100)
ctrl <- mlr:::makeTuneControlMBO(mbo.control = mbo.ctrl)

# set parallel backend
library(parallel)
library(parallelMap)
(numCores <- detectCores())
parallelStartSocket(cpus = numCores - 1)

# parameter tuning
start_time <- Sys.time()
xgb.par_tuning <- tuneParams(learner = xgb.lrn_tune, 
                     task = traintask, 
                     resampling = r_inst,
                     measures = mlr::rmse,
                     par.set = xgb.params,
                     control = ctrl,
                     show.info = TRUE)
end_time <- Sys.time()
(time_elapsed <- end_time - start_time)
parallelStop()
xgb.par_tuning
xgb.par_tuning$x
xgb.par_tuning$y
xgb.par_tuning$resampling$test.inds
xgb.par_tuning$opt.path


# Show up all test performance via RMSE
xgb.lrn_testPerf <- makeLearner("regr.xgboost", par.vals = xgb.par_tuning$x)
xgb.lrn_testPerf <- makePreprocWrapperCaret(xgb.lrn_testPerf, ppc.scale = TRUE, ppc.center = TRUE)

xgb.error_resample = resample(
  learner = xgb.lrn_testPerf, 
  task = traintask, 
  resampling = r_inst, 
  # measures = list(rmse,mae,rsq),
  measures = list(mlr::rmse, mlr::mae, mlr::rsq),
  models = TRUE,
  show.info = TRUE,
  keep.pred = TRUE)
xgb.error_resample
all.equal(as.numeric(xgb.error_resample$aggr[1]), as.numeric(xgb.par_tuning$y)) 


##-------------
# 3.2 Predictive performance using rmse, mae, r2
##-------------
## mean and sd for rmse, mae and R-squared
(rmse_mean <- mean(xgb.error_resample$measures.test$rmse))
(rmse_sd <- sd(xgb.error_resample$measures.test$rmse))

(mae_mean <- mean(xgb.error_resample$measures.test$mae))
(mae_sd <- sd(xgb.error_resample$measures.test$mae))

(R2_mean <- mean(xgb.error_resample$measures.test$rsq))
(R2_sd <- sd(xgb.error_resample$measures.test$rsq))

#
xgb.error_table <- data.frame(matrix(ncol = 2, nrow = 3))
colnames(xgb.error_table) <- c("Mean", "SD")
rownames(xgb.error_table) <- c("RMSE", "MAE", "R-squared")
xgb.error_table$Mean <- c(rmse_mean, mae_mean, R2_mean)
xgb.error_table$SD <- c(rmse_sd, mae_sd, R2_sd)
xgb.error_table %<>% mutate_if(is.numeric, round, digits=3)
xgb.error_table

# Visualize the true values and predicted values
# true values & predicted values
xgb.pred = getRRPredictions(xgb.error_resample)
xgb.true_pred <- data.frame(Truth = getPredictionTruth(xgb.pred),
                           Pred = pmax(getPredictionResponse(xgb.pred), 0),
                           Repeats = as.factor((xgb.pred$data$iter-1) %/% 10 + 1))


# calculate the 95% confidence interval of predicted values versus true values
dim(xgb.true_pred)
(sample_num <- length(unique(xgb.true_pred$Truth)))
xgb.true_pred_ci <- data.frame(matrix(ncol = 5, nrow = sample_num))
colnames(xgb.true_pred_ci) <- c("True", "Mean", "Upper","Median", "Lower")
for(x in seq(sample_num)){
  # x = 1
  x_true = unique(xgb.true_pred$Truth)[x]
  x_index = which(xgb.true_pred$Truth == x_true)
  x_pred = xgb.true_pred[x_index, ]$Pred
  x_mean = mean(x_pred)
  xgb.true_pred_ci[x,] = c(x_true, x_mean, quantile(x_pred, probs = c(0.975, 0.5, 0.025)))
}


## ggplot with confidence
ggplot(xgb.true_pred_ci, aes(x = True, y = Median)) +
  geom_point(alpha = 0.9) +
  geom_errorbar(aes(ymin = Lower,  ymax = Upper), width = 0, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.8, alpha = 0.9) +
  labs(x = "True values", y = "Predicted values (95% confidence interval)") + 
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
           label = paste0("RMSE = ", xgb.error_table[1,1], "\U00B1", xgb.error_table[1,2],";")) +
  annotate(geom = "text", x = 0.185, y = 0.44, size= 4, col = "black", fontface = "bold",
           label = paste0("MAE = ", xgb.error_table[2,1],"\U00B1", xgb.error_table[2,2], ";")) +
  annotate(geom = "text", x = 0.33, y = 0.44, size= 4, col = "black", fontface = "bold",
           label = paste0("R-squared = ", xgb.error_table[3,1], "\U00B1", xgb.error_table[3,2])) +
  coord_cartesian(xlim = c(0, 0.4), ylim = c(0, 0.4), expand = TRUE, clip="off") -> xgb.true_pred_ci_plt

xgb.true_pred_ci_plt

#
ggsave(filename = "xgb.true_pred_ci_plt_rep10cv10_Asite_20240907.tiff",
       plot = xgb.true_pred_ci_plt,
       width = 16,  # <=19.05cm
       height = 12, # <=22.225cm
       units= "cm",
       dpi= 300,
       compression = "lzw")  


# ######
# # ggplot
# xgb.true_pred_plt <- ggplot(xgb.true_pred, aes(x = Truth, y = Pred, color = Repeats)) +
#   geom_point(alpha = 0.7) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.8) +
#   labs(x = "True values", y = "Predicted values") + 
#   scale_x_continuous(limits = c(0, 0.4), expand = c(0.02,0)) +
#   theme_bw(base_size=12) +
#   theme(plot.margin = unit(c(3, 1, 1, 1), "lines"),
#         axis.text.x=element_text(size=10, angle = 0, vjust=0.5, hjust = 1, colour = "black"),
#         axis.text.y=element_text(size=10, colour = "black"),
#         axis.title=element_text(size=12, colour = "black", face = "bold"),
#         plot.title = element_text(size = 11, face="bold"),
#         legend.background = element_rect(fill="transparent",
#                                          linewidth = 0.3, linetype = "solid",
#                                          colour = "black"),
#         legend.box.background = element_rect(colour = "black"),
#         legend.title=element_text(size= 8, colour = "black"),
#         legend.text=element_text(size= 8, colour = "black"),
#         legend.position = c("none"),#"left",#c(.75, .95),
#         legend.direction="horizontal",
#         legend.justification = c("right", "top"),
#         panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.8)
#   ) +
#   annotate(geom = "text", x = 0.32, y =0.38, size= 5, col = "red",
#            label = paste0("1:1 line"), hjust = "left") +
#   annotate(geom = "text", x = 0.05, y = 0.45, size= 4, col = "black", fontface = "bold",
#            label = paste0("RMSE = ", xgb.error_table[1,1], "\U00B1", xgb.error_table[1,2],";")) +
#   annotate(geom = "text", x = 0.185, y = 0.45, size= 4, col = "black", fontface = "bold",
#            label = paste0("MAE = ", xgb.error_table[2,1],"\U00B1", xgb.error_table[2,2], ";")) +
#   annotate(geom = "text", x = 0.33, y = 0.45, size= 4, col = "black", fontface = "bold",
#            label = paste0("R-squared = ", xgb.error_table[3,1], "\U00B1", xgb.error_table[3,2])) +
#   coord_cartesian(xlim = c(0, 0.4), ylim = c(0, 0.4), expand = TRUE, clip="off")
# 
# xgb.true_pred_plt  
# 
# #
# ggsave(filename = "xgb.true_pred_plt_rep10cv10_Asite.tiff",
#        plot = xgb.true_pred_plt,
#        width = 16,  # <=19.05cm
#        height = 12, # <=22.225cm
#        units= "cm",
#        dpi= 300,
#        compression = "lzw")  
#####


## residual plot: true values - predicted values

###
xgb.res_df <- data.frame(Truth = getPredictionTruth(xgb.pred), 
                        Error =  getPredictionTruth(xgb.pred) - getPredictionResponse(xgb.pred),
                        Repeats = as.factor((xgb.pred$data$iter-1) %/% 10 + 1))

xgb.res_df %>% head(4)


# calculate the 95% confidence interval of residuals versus true values
dim(xgb.res_df)
(error_num <- length(unique(xgb.res_df$Truth)))
xgb.true_error_ci <- data.frame(matrix(ncol = 5, nrow = error_num))
colnames(xgb.true_error_ci) <- c("Sample_true", "Error_mean", "Error_upper","Error_median", "Error_lower")
for(x in seq(error_num)){
  # x = 1
  x_true = unique(xgb.res_df$Truth)[x]
  x_index = which(xgb.res_df$Truth == x_true)
  x_error = xgb.res_df[x_index, ]$Error
  x_mean = mean(x_error)
  xgb.true_error_ci[x,] = c(x_true, x_mean, quantile(x_error, probs = c(0.975, 0.5, 0.025)))
}

xgb.true_error_ci %>% head(5)
ggplot(xgb.true_error_ci, aes(x = Sample_true, y = Error_median)) +
  geom_point(alpha = 0.9) +
  geom_abline(intercept = 0, slope = 0, color = "red", linewidth = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin = Error_lower,  ymax = Error_upper), width = 0, alpha = 0.7) +
  geom_smooth() +
  labs(x = "True values", y = "Residuals (True values - predicted values)") +
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
  ) -> xgb.res_plt_ci
xgb.res_plt_ci


ggsave(filename = "xgb.true_resi_ci_plt_rep10cv10_Asite_20240907.tiff",
       plot = xgb.res_plt_ci,
       width = 16,  # <=19.05cm
       height = 12, # <=22.225cm
       units= "cm",
       dpi= 300,
       compression = "lzw")   
 
  
#
xgb.res_plt_byRepeat <- ggplot(xgb.res_df, aes(x = Truth, y = Error, color = Repeats)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_smooth() +
  labs(y = "Residuals (True values - predicted values)", x = "True values") +
  theme_bw(base_size=12) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"),
        axis.text.x=element_text(size=8, angle = 0, vjust=0.5, hjust = 0.2, colour = "black"),
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title=element_text(size=12, colour = "black", face = "bold"),
        plot.title = element_text(size = 11, face="bold"),
        legend.background = element_rect(fill="transparent",
                                         linewidth = 0.3, 
                                         linetype = "solid",
                                         colour = "black"),
        legend.box.background = element_rect(colour = "black"),
        # panel.grid= element_blank(),
        legend.title=element_text(size= 8, colour = "black"),
        legend.text=element_text(size= 8, colour = "black"),
        # legend.position =c(.55, .95),#"left",#c(.75, .95),
        legend.position = "top",#"left",#c(.75, .95),
        legend.direction="horizontal",
        legend.justification = c("top"),
        panel.border = element_rect(colour = "black", fill=NA, size= 1.5)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  facet_wrap(~ Repeats, nrow = 2)

xgb.res_plt_byRepeat 

#
ggsave(filename = "xgb.true_resi_byRepeat_plt_rep10cv10_Asite_20240907.tiff",
       plot = xgb.res_plt_byRepeat,
       width = 24,  # <=19.05cm
       height = 16, # <=22.225cm
       units= "cm",
       dpi= 300,
       compression = "lzw")   


  
# ##### 
# library(cowplot)
# xgb.res_plt <- ggplot(xgb.res_df, aes(x = Truth, y = Error, color = Repeats)) +
#   geom_point(alpha = 0.7) +
#   geom_hline(yintercept = 0, linewidth = 0.6) +
#   labs(y = "Residuals (True values - predicted values)", x = "True values") +
#   theme_bw(base_size=12) +
#   theme(plot.margin = unit(c(1, 1, 1, 1), "lines"),
#         axis.text.x=element_text(size=10, angle = 0, vjust=0.5, hjust = 1, colour = "black"),
#         axis.text.y=element_text(size=10, colour = "black"),
#         axis.title=element_text(size=12, colour = "black", face = "bold"),
#         plot.title = element_text(size = 11, face="bold"),
#         legend.background = element_rect(fill="transparent",
#                                          size = 0.3, linetype = "solid",
#                                          colour = "black"),
#         legend.box.background = element_rect(colour = "black"),
#         # panel.grid= element_blank(),
#         legend.title=element_text(size= 8, colour = "black"),
#         legend.text=element_text(size= 8, colour = "black"),
#         # legend.position =c(.55, .95),#"left",#c(.75, .95),
#         legend.position = "top",#"left",#c(.75, .95),
#         legend.direction="horizontal",
#         legend.justification = c("top"),
#         panel.border = element_rect(colour = "black", fill=NA, size= 1.5)
#         ) +
#   guides(color = guide_legend(nrow = 1, byrow = TRUE))
# xgb.res_plt
# 
# #
# ggsave(filename = "xgb.res_plt_rep10cv10.tiff",
#        plot = xgb.res_plt,
#        width = 16,  # <=19.05cm
#        height = 12, # <=22.225cm
#        units= "cm",
#        dpi= 300,
#        compression = "lzw")
############

##-------------
# 3.3 train model
##-------------
## 
# show up hyper-parameters
xgb.par_tuned <- setHyperPars(learner = xgb.lrn_tune, par.vals = xgb.par_tuning$x)
xgb.par_tuned$par.vals

# train model
xgb.model <- mlr::train(xgb.par_tuned, task = traintask) 
xgb.model

# save trained model
save(file="saveXGBmodel_rep10cv10_20240907.rdata", list=c("xgb.model"))
# load binary model to R
load(file="saveXGBmodel_rep10cv10.rdata")


##-------------
# 3.4 feature importance
##-------------
xgb.var_imp = getFeatureImportance(xgb.model) # type2(default):the mean decrease in node impurity; type1:the mean decrease in accuracy calculated on OOB data
xgb.var_imp_data = as.data.frame(xgb.var_imp$res)
xgb.var_imp_data_desc <- dplyr::arrange(xgb.var_imp_data, desc(importance))

xgb.var_imp_data_desc$variable
variable_NoBB1 <- subset(xgb.var_imp_data_desc, !grepl("(B1|B|A)$", xgb.var_imp_data_desc$variable))
xgb.var_imp_data_desc 
variable_NoBB1

as.matrix(variable_NoBB1)
xgb.ggplot.shap.summary(as.matrix(variable_NoBB1))

#
ggplot(variable_NoBB1, aes(x=variable, y=importance)) + 
  geom_bar(stat = "identity") +
  coord_flip()




library("SHAPforxgboost")
xgb.model$learner$next.learner
getLearnerModel(xgb.model)
model_trained = xgboost::xgboost(data = as.matrix(data_all[, -1]),
                 label = data_all$HPC,
                 nrounds = 304,
                 eta = 0.259,
                 gamma=5.81e-06,
                 max_depth=9,
                 min_child_weight=5.15,
                 lambda=0.627,
                 alpha=0.00902
                 )

shap_values <- shap.values(xgb_model = model_trained, X_train = as.matrix(data_all[, -1]))
shap_long_hpc <- shap.prep(xgb_model = model_trained, X_train = as.matrix(data_all[, -1]))
shap_long_hpc
shap_values$mean_shap_score
shap_values$shap_score
shap.plot.summary(shap_long_hpc)
shap.plot.dependence(data_long = shap_long_hpc, "temperature")
shap.plot.dependence(data_long = shap_long_hpc, "fraction_A")
shap.plot.dependence(data_long = shap_long_hpc, "ir_A1")



(varib_noBB1 <- subset(shap_long_hpc, !grepl("(B1|B|A)$", shap_long_hpc$variable)))
shap.plot.summary(varib_noBB1,
                  scientific = TRUE) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(colour = "black", face = "bold"),
        axis.title.y = element_text(colour = "black", face = "bold"),
        # panel.background = element_rect(fill = "gray90"),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.6),
        legend.position = "none") -> xgb.VarImp_plt_shap
xgb.VarImp_plt_shap

#
ggsave(filename = paste0("xgb.VarImp_plt_shap.tiff"),
       plot = xgb.VarImp_plt_shap,
       width = 19,  # <=19.05cm
       height = 12, # <=22.225cm
       units= "cm",
       dpi= 300,
       compression = "lzw")




#
variable_NoBB1[!variable_NoBB1$variable %in% c("temperature","pH2O"),] %>% 
  ggplot(aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", color = "#F8766D", fill = "#F8766D", width = 0.6) +
  coord_polar("y", start = 0) +
  xlab("Features") + ylab("Relative importance") +
  theme(plot.title = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(colour = "black", face = "bold"),
        axis.title.y = element_text(colour = "black", face = "bold"),
        panel.background = element_rect(fill = "gray90"),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.6),
        legend.position = "none") -> xgb.VarImp_plt

xgb.VarImp_plt


#
ggsave(filename = paste0("xgb.VarImp_plt_rep10cv10.tiff"),
       plot = xgb.VarImp_plt,
       width = 19,  # <=19.05cm
       height = 12, # <=22.225cm
       units= "cm",
       dpi= 300,
       compression = "lzw")



#




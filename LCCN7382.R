## predict the pc with the temperature range and water pressure
# (Temper_range = seq(300, 700, by = 10))
# (WaterPre_range = seq(0, 0.2, by = 0.01))
# for the composite  La0.7Ca0.3Co0.8Ni0.2O3


### Plot 3D surf

getwd()

library(readxl)
library(tidyr)
library(plyr)

df_elem <- read.csv("elements_data.csv") %>% 
  as.data.frame()
source("AA1B_descriptor_customize.R")

df_elem_user = as.data.frame(matrix(nrow = 1, ncol = 12))
colnames(df_elem_user) = c( "A", "A_valence", "A_fraction",
                       "A1", "A1_valence", "A1_fraction",
                       "B", "B_valence", "B_fraction", 
                       "B1", "B1_valence", "B1_fraction")

df_elem_user$A[1] = "La"
df_elem_user$A_valence = 3
df_elem_user$A_fraction = 0.8
df_elem_user$A1 = "Ba"
df_elem_user$A1_valence = 2
df_elem_user$A1_fraction = 0.2
df_elem_user$B = "Co"
df_elem_user$B_valence = 3
df_elem_user$B_fraction = 1
# df_elem_user$B1 = "Ni"
# df_elem_user$B1_valence = 2
# df_elem_user$B1_fraction = 0.2


(Temper_range = seq(450, 700, by = 1))
# (Temper_range = seq(300, 700, by = 10))
(WaterPre_range = seq(0, 0.05, by = 0.001))
# (WaterPre_range = seq(0, 0.05, by = 0.001))
(num_temper = length(Temper_range) )
(num_water = length(WaterPre_range) )
(total_iter = num_temper * num_water)

start_time <- Sys.time()
predict_pc_mesh_LBC <- matrix(NA, nrow = length(Temper_range), ncol = length(WaterPre_range))
num_iter = 1
for(i in seq(length(Temper_range))){
  for (j in seq(length(WaterPre_range))){
    print(paste0("Among ", total_iter, " , current calculating: ", num_iter))
    num_iter = num_iter + 1
    predict_pc_mesh_LBC[i,j] =  predict_pc(Temper_range[i], WaterPre_range[j])
  }
}
end_time <- Sys.time()
(time_elapsed <- (end_time - start_time))

# save(predict_pc_mesh_LCCN7382, file = "predict_pc_mesh_LCCN7382.RData" )
write.csv(predict_pc_mesh_LBC,"predict_pc_mesh_LBC.csv", row.names=FALSE) 
predict_pc_mesh_LBC <- read.csv("predict_pc_mesh_LBC.csv") %>% 
  as.matrix()

#rf_fit <- readRDS("rf_fit_0714_0719tune.rds") # tune parameters for rf


load(file="saveXGBmodel_rep10cv10_20231024.rdata")


predict_pc = function(x, y){
  #x = 600
  #y = 0.03
  df_expe_user = as.data.frame(matrix(nrow = 1, ncol = 2))
  colnames(df_expe_user) = c("temperature", "pH2O")
  df_expe_user$temperature = x
  df_expe_user$pH2O = y
  df_expe_elem_user = cbind(df_expe_user, df_elem_user)
  newdata = descriptor_customize(df_expe_elem_user)
  newdata <- within(newdata, rm("fraction_A1", "fraction_B1"))
  (new_pred = predict(xgb.model, newdata = newdata))
  #predicted_pc_rf = predict(rf_fit, newdata = newdata)
  # new_pred$data$response
  return(new_pred$data$response)
}


# predict_pc_composites <- numeric(num_composite)
# df_expe_elem_tosave <- data.frame()
# is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
# for(i in seq(num_composite)){
#   i = 1
#   (newdata = descriptor_customize(df_expe_elem_user[i,]))
#   newdata <- within(newdata, rm("fraction_A1", "fraction_B1"))
#   
#   df_expe_elem_tosave = rbind(df_expe_elem_tosave, newdata)
#   # new_pred = predict(xgb.model, newdata = newdata)
#   new_pred = predict(xgb.model, newdata = newdata)
#   predict_pc_composites[i] <- max(as.numeric(new_pred$data), 0)
#   
#   print(i)
# }





# using plot3D
library(plot3D)
expe_mesh = mesh(Temper_range, WaterPre_range)
# with(expe_mesh, predict_pc_mesh(Temper_range, WaterPre_range))
u <- expe_mesh$x
v <- expe_mesh$y
par(mar = c(1, 1, 1, 1))
tiff(file="predict_pc_mesh_LBC_xgb_360.tiff",
     width=10, height=8, units="in", res=300)
surf3D(u, v, predict_pc_mesh_LBC, 
       colvar = predict_pc_mesh_LBC, 
       colkey = list(shift = 0, las = 1, col.box = "white", side = 1, 
                     length = 0.5, width = 0.6, dist = -0.05, cex.axis = 0.8),
       box = TRUE, bty = "b2", 
       expand = 0.8,
       space = 0.1,
       breaks = ,
       ticktype = "detailed",
       d = 2, 
       # labelpad = 20,
       xlim = range(u), 
       ylim = range(v), 
       # zlim = range(predict_pc_mesh),
       zlim = c(0.02, 0.08),
       # breaks = range(predict_pc_mesh),
       xlab = "", ylab = "", zlab = "", cex.lab = 1,
       cex.axis = 1.5, 
       yaxt = "n",
       phi = 15, theta = 60) # adjust viewing direction:theta = 30o,80o and 360o were used

dev.off()











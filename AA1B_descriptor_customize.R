#######
# Based on the oxide-composition and the element information to customize the descriptors,
# which would be used for the following model training
#######
descriptor_customize <- function(AA1BB1O3_sample){
  
  # AA1BB1O3_sample = df_expe_elem
 #(AA1BB1O3_sample = df_expe_elem_user)

  # nbr of samples
  (N = dim(AA1BB1O3_sample)[1]) 
  
  # fraction of composition
  (fraction_A = AA1BB1O3_sample['A_fraction'])
  (fraction_A1 = AA1BB1O3_sample['A1_fraction'])
  (fraction_B = AA1BB1O3_sample['B_fraction'])
  (fraction_B1 = AA1BB1O3_sample['B1_fraction'])
  
  # extract element attribute of A-site host element in composition
  aw_A = numeric()
  ad_A = numeric()
  mp_A = numeric()
  fie_A = numeric()
  eln_A = numeric()
  ir_A = numeric()
  for (n in seq(N)){
    if (fraction_A[n,] > 0){
      a = AA1BB1O3_sample['A'][n,]
      v = AA1BB1O3_sample['A_valence'][n,]
      i = which((df_elem['Atom'] == a) & (df_elem['Valence'] == v))
      if(length(i) == 0){
        print("Invalid element!")
      }
      aw_A[n] = df_elem['Atomic.weight'][i,]
      ad_A[n] = df_elem['Atomic.density'][i,]
      mp_A[n] = df_elem['Melting.point'][i,] 
      fie_A[n] = df_elem['First.ionization.energy'][i,]
      eln_A[n] = df_elem['Electronegativity'][i,]
      ir_A[n] = df_elem['Ionic.radius.VI'][i,]
    }
  }
  
  
  # extract element attribute of A-site dopant element in composition
  aw_A1 = numeric()
  ad_A1 = numeric()
  mp_A1 = numeric()
  fie_A1 = numeric()
  eln_A1 = numeric()
  ir_A1 = numeric()
  for (n in seq(N)){
    if ((fraction_A1[n,] == 0) | (is.na(fraction_A1[n,]))){
      aw_A1[n] = 0
      ad_A1[n] = 0
      mp_A1[n] = 0
      fie_A1[n] = 0
      eln_A1[n] = 0
      ir_A1[n] = 0
    } else if (fraction_A1[n,] > 0){
      a = AA1BB1O3_sample["A1"][n,]
      v = AA1BB1O3_sample['A1_valence'][n,]
      i = which((df_elem['Atom'] == a) & (df_elem['Valence'] == v))
      if(length(i) == 0){
        print("Invalid element!")
      }
      aw_A1[n] = df_elem['Atomic.weight'][i,]
      ad_A1[n] = df_elem['Atomic.density'][i,]
      mp_A1[n] = df_elem['Melting.point'][i,] 
      fie_A1[n] = df_elem['First.ionization.energy'][i,]
      eln_A1[n] = df_elem['Electronegativity'][i,]
      ir_A1[n] = df_elem['Ionic.radius.VI'][i,]
    }
  }
  
  
  # extract element attribute of B-site host element in composition
  aw_B = numeric()
  ad_B = numeric()
  mp_B = numeric()
  fie_B = numeric()
  eln_B = numeric()
  ir_B = numeric()
  for (n in seq(N)){
    if ((fraction_B[n,] == 0) | (is.na(fraction_B[n,]))){
      next
    } else if ((fraction_B[n,] > 0)){
      a = AA1BB1O3_sample['B'][n,]
      v = AA1BB1O3_sample['B_valence'][n,]
      i = which((df_elem['Atom'] == a) & (df_elem['Valence'] == v))
      if(length(i) == 0){
        print("Invalid element!")
      }
      aw_B[n] = df_elem['Atomic.weight'][i,]
      ad_B[n] = df_elem['Atomic.density'][i,]
      mp_B[n] = df_elem['Melting.point'][i,] 
      fie_B[n] = df_elem['First.ionization.energy'][i,]
      eln_B[n] = df_elem['Electronegativity'][i,]
      ir_B[n] = df_elem['Ionic.radius.VI'][i,]
    }
  }
  
  # extract element attribute of B-site dopant element in composition
  aw_B1 = numeric()
  ad_B1 = numeric()
  mp_B1 = numeric()
  fie_B1 = numeric()
  eln_B1 = numeric()
  ir_B1 = numeric()
  for (n in seq(N)){
    if(is.na(AA1BB1O3_sample['B1'])){
      aw_B1[n] = 0
      ad_B1[n] = 0
      mp_B1[n] = 0
      fie_B1[n] = 0
      eln_B1[n] = 0
      ir_B1[n] = 0
      AA1BB1O3_sample['B1_fraction'] = 0
      AA1BB1O3_sample['B_fraction'] = 1
    } else if (AA1BB1O3_sample['B1'] == AA1BB1O3_sample['B']){
      aw_B1[n] = 0
      ad_B1[n] = 0
      mp_B1[n] = 0
      fie_B1[n] = 0
      eln_B1[n] = 0
      ir_B1[n] = 0
      AA1BB1O3_sample['B1_fraction'] = 0
      AA1BB1O3_sample['B_fraction'] = 1
    } else if ((fraction_B1[n,] == 0) | (is.na(fraction_B1[n,]))){
      aw_B1[n] = 0
      ad_B1[n] = 0
      mp_B1[n] = 0
      fie_B1[n] = 0
      eln_B1[n] = 0
      ir_B1[n] = 0
      } else if ((fraction_B1[n,] > 0)){
      a = AA1BB1O3_sample['B1'][n,]
      v = AA1BB1O3_sample['B1_valence'][n,]
      i = which((df_elem['Atom'] == a) & (df_elem['Valence'] == v))
      if(length(i) == 0){
        print("Invalid element!")
      }
      aw_B1[n] = df_elem['Atomic.weight'][i,]
      ad_B1[n] = df_elem['Atomic.density'][i,]
      mp_B1[n] = df_elem['Melting.point'][i,] 
      fie_B1[n] = df_elem['First.ionization.energy'][i,]
      eln_B1[n] = df_elem['Electronegativity'][i,]
      ir_B1[n] = df_elem['Ionic.radius.VI'][i,]
    }
  }

  
  # element attributes of A site, A1 site, B site and B1 site
  df_elem_attr <- data.frame(aw_A, ad_A, mp_A, fie_A, eln_A, ir_A,
                             aw_A1, ad_A1, mp_A1, fie_A1, eln_A1, ir_A1,
                             aw_B, ad_B, mp_B, fie_B, eln_B, ir_B,
                             aw_B1, ad_B1, mp_B1, fie_B1, eln_B1, ir_B1)
  
  colnames(df_elem_attr) <- c("aw_A", "ad_A", "mp_A", "fie_A", "eln_A", "ir_A",
                              "aw_A1", "ad_A1", "mp_A1", "fie_A1", "eln_A1", "ir_A1",
                              "aw_B", "ad_B", "mp_B", "fie_B", "eln_B", "ir_B",
                              "aw_B1", "ad_B1", "mp_B1", "fie_B1", "eln_B1", "ir_B1")
  
  
  # experimental conditions
  feature_names = c("temperature", "pH2O", "A_fraction",
                    "A1_fraction", "B_fraction", "B1_fraction")
  
  
  df_expe_selected = AA1BB1O3_sample[feature_names]
  colnames(df_expe_selected) = c("temperature", "pH2O", "fraction_A",
                                "fraction_A1", "fraction_B", "fraction_B1")
  
  #
  # combine the experimental conditions, the element attributes and calculated element descriptors
  # df_sample = cbind(df_expe_selected, df_elem_attr, df_descriptor)
  (df_sample = cbind(df_expe_selected, df_elem_attr))
  # 
  
  return(df_sample)  
}







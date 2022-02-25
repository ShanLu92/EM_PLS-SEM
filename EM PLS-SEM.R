rm(list=ls())

library(seminr)
library(DMwR2)
library(mice) 

# define the model structure as Seminr suggested
# here we are using the European Customer Satisfaction Index model
dd_mm <- constructs(
  composite("Image",        multi_items("IMAG", 1:5), weights = mode_A),
  composite("Expectation",  multi_items("CUEX", 1:3), weights = mode_A),
  composite("Quality",      multi_items("PERQ", 1:7), weights = mode_A),
  composite("Value",        multi_items("PERV", 1:2), weights = mode_A),
  composite("Satisfaction", multi_items("CUSA", 1:3), weights = mode_A),
  composite("Complaints",   single_item("CUSCO")),
  composite("Loyalty",      multi_items("CUSL", 1:3), weights = mode_A)
) 

dd_sm <- relationships(
  paths(from = "Image",        to = c("Expectation", "Satisfaction", "Loyalty")),
  paths(from = "Expectation",  to = c("Quality", "Value", "Satisfaction")),
  paths(from = "Quality",      to = c("Value", "Satisfaction")),
  paths(from = "Value",        to = c("Satisfaction")),
  paths(from = "Satisfaction", to = c("Complaints", "Loyalty")),
  paths(from = "Complaints",   to = "Loyalty")
)


n = 100
cs_times = 501 # iteration for EM PLS-SEM
threshold = 1e-04

rdm_times = 200

set.seed(121)

for (na_pcent in c(0.01, 0.02, .04, 0.08, 0.12, 0.16)){
  result = data.frame(matrix(ncol = rdm_times, nrow = 7))
  print (na_pcent)
  for (rdm in 1:rdm_times){
    
    ##### simulate datasets with missing values #####
    mobi = read.csv(paste0(toString(rdm), '_simul_esci_n=', toString(n), '.csv'))
    mobi = mobi[,-1]
    dd = scale(mobi)
    # head(dd)
    
    suppressMessages(invisible(capture.output(plsm <- estimate_pls(data = dd,
                         measurement_model = dd_mm,
                         structural_model = dd_sm))))
    
    mmMatrix = plsm$mmMatrix
    smMatrix = plsm$smMatrix
    
    scorenames = c("Image", "Expectation", "Quality", 
                   "Value", "Satisfaction", "Loyalty")
    
    na_n = ceiling(nrow(dd)*ncol(dd)*na_pcent)
    na_loc = sample(nrow(dd)*ncol(dd), na_n)
    
    na_row = na_loc%/%ncol(dd)+1
    na_col = na_loc%%ncol(dd)+1 # start from 1
    
    dd = data.frame(scale(mobi))
    
    ##############  EM PLS-SEM ########################
    # randomly assign values to empty cells from the available values
    start_time <- Sys.time()
      
    na_value_real = c()
    for (i in 1:length(na_n)){
      na_value_real = append(na_value_real, dd[na_row[i], na_col[i]])
      dd[na_row[i], na_col[i]] = sample(dd[, na_col[i]], 1)
    }
    
    # inner_mae_old = 100 # if use two consecutive inner mae diff<1e-06, then it alway yeils small mae in 1501 iteration
    for (cs in 1:cs_times){
      
      suppressMessages(invisible(capture.output(
        nan_plsm <- estimate_pls(data = dd,
                                 measurement_model = dd_mm,
                                 structural_model = dd_sm)
      )))
      
      for (i in 1:length(na_n)){
        
        na_colname_now = colnames(dd)[na_col[i]]
        na_row_now = na_row[i]
        na_scorename_now = mmMatrix[mmMatrix[,2]==na_colname_now, 'construct']
        mvars_name = mmMatrix[mmMatrix[,1]==na_scorename_now, 'measurement']
        na_score = nan_plsm$construct_scores[, na_scorename_now]
        
        if (length(mvars_name)>1){
          # fit jth col first, fixed others
          # variable names related to the scores, not include the missing value variable
          mvars_name = setdiff(mvars_name, na_colname_now) 
          # variable outer weights related to the scores, not include the missing value variable
          mvar_weight = nan_plsm$outer_weights[mvars_name, na_scorename_now]
          # variable outer weights of the missing value 
          na_weight = nan_plsm$outer_weights[na_colname_now, na_scorename_now]
          # in composite design, the w_ij * x_ij = scores -sum(other x_ij * outer_weights)
          na_fill = as.numeric(na_score[na_row_now])
          
          for (j in 1:length(mvars_name)){
            na_fill = na_fill - as.numeric(mvar_weight[j] * dd[na_row_now, mvars_name[j]])
          }
          na_fill = na_fill/na_weight
          dd[na_row_now, na_colname_now] = na_fill
        }else{#some LV has only one MV: nan_plsm$outer_weights[na_colname_now, na_scorename_now]=1
          na_fill = as.numeric(na_score[na_row_now])
        }
      }
      outer_mae_new = (1/nrow(mmMatrix))*sum(abs(nan_plsm$outer_weights - plsm$outer_weights))
      inner_mae_new = (1/nrow(smMatrix))*sum(abs(nan_plsm$path_coef - plsm$path_coef))
      
      # this round and last round inner weights difference < 1e06, then stop
      if (length(inner_mae_new)>1){
        if (abs(inner_mae_new) < 1e-06 ) break
      }
      # inner_mae_old = inner_mae_new
    }
      
    end_time <- Sys.time()
      
    na_value_impute = c()
    for (i in 1:length(na_n)){
      na_value_impute = append(na_value_impute, dd[na_row[i], na_col[i]])
    }
    
    # the mae of outer weights
    result[1,rdm] = outer_mae_new
    # the mae of inner weights
    result[2,rdm]= inner_mae_new
    result[3,rdm] = sum(abs(nan_plsm$rSquared[1,]-plsm$rSquared[1,]))
    result[4,rdm] = sum(abs(nan_plsm$rSquared[2,]-plsm$rSquared[2,]))
    result[5,rdm] = cs
    result[6,rdm] = (1/length(na_value_impute))*sum(abs(na_value_impute-na_value_real))
    result[7,rdm] = end_time - start_time
  }
  
  rownames(result) = c('outer_mae', 'inner_mae',
                       'rsqua', 'rsqua_adjust',
                       'cs', 'real-impute_mae',
                       'time')
  
  write.csv(result, paste0('EM_n=',toString(n),'_cs=',toString(cs), 
                           '_threshold=', toString(threshold),
                           '_percent=',toString(na_pcent),'_result.csv'))
  
}
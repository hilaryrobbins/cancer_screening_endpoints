#function

data_corrfigure <- function(datax, cancer_list, var_reduction){
  data_figure <- data.frame(cancer_type2=NA,Pearson=NA,Pearson_low=NA,Pearson_up=NA,
                            Linear_intercept=NA,Linear_x=NA,Linear_x_SE=NA, Linear_r2=NA, 
                            n_pub=NA)
  
  for (i in c(1:length(cancer_list)) ) {
    data_cancer <- datax %>% filter(cancer_type2==cancer_list[i])
    data_figure[i,"cancer_type2"]<-cancer_list[i]
    data_figure[i,"Pearson"]<-with(data_cancer, 
                                    cor(get(var_reduction), Reduction_CaMort_arm2vs1_perc, 
                                        method="pearson",use="complete.obs"))
    data_figure[i,"Pearson_low"]<-tryCatch(with(data_cancer, 
                                                 ci_cor(get(var_reduction),Reduction_CaMort_arm2vs1_perc,
                                                        probs=c(0.025,0.975), method="pearson", type="normal"))[["interval"]][1],
                                            error=function(e){return(-1)})
    data_figure[i,"Pearson_up"]<-tryCatch(with(data_cancer, 
                                                ci_cor( get(var_reduction),Reduction_CaMort_arm2vs1_perc,
                                                        probs=c(0.025,0.975), method="pearson", type="normal"))[["interval"]][2],
                                           error=function(e){return(1)})
    model_lm <- lm(Reduction_CaMort_arm2vs1_perc ~ get(var_reduction), data=data_cancer)
    data_figure[i,"Linear_intercept"]<-broom::tidy(model_lm)[broom::tidy(model_lm)$term=="(Intercept)","estimate", drop=T]
    data_figure[i,"Linear_x"]<-broom::tidy(model_lm)[grepl("var_reduction",broom::tidy(model_lm)$term),"estimate", drop=T]
    data_figure[i,"Linear_x_SE"]<-broom::tidy(model_lm)[grepl("var_reduction",broom::tidy(model_lm)$term),"std.error", drop=T]
    data_figure[i,"Linear_r2"]<-summary(model_lm)[["r.squared"]]
    data_figure[i,"n_pub"] <- with(datax %>% filter(!is.na(get(var_reduction))), sum(cancer_type2==cancer_list[i]))
  }
  data_figure <- data_figure %>%mutate(Pearson_all=paste0("Pearson \u03c1=",specify_decimal(Pearson,2)," (",specify_decimal(Pearson_low,2)," to ",specify_decimal(Pearson_up,2),")"),
                                        Linear_slope=paste0("Linear regression y=",specify_decimal(Linear_intercept,2),"+",specify_decimal(Linear_x,2),"x"),
                                        Linear_slope= str_replace(Linear_slope,"\\+-","\\-"),
                                        Linear_r2=paste0("R-squared=",specify_decimal(Linear_r2,2)),
                                        Linear_x_conf=paste0( specify_decimal(Linear_x,2)," (",specify_decimal(Linear_x-1.96*Linear_x_SE,2),"-",specify_decimal(Linear_x+1.96*Linear_x_SE,2),")"))
  
  return(data_figure)
}

#data_figure1<-data_corrfigure(mced.filtered, cancer_list, "Reduction_LateStage_arm2vs1_perc")


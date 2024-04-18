
#############
# Purpose: Cancer screening RCT endpoints systematic review
# Code authors: Xiaoshuang Feng, Hana Zahed, Hilary Robbins
# Citation for article: Feng et al, JAMA 2024
# !!! See README.txt for important details.
#############

rm(list=ls(all=TRUE)) 
setwd("") # insert path where data and the "function" script is saved
path_to_save<- "" # insert path to folder where results should be saved 

###Packages needed
packages <- c("confintr","tidyverse","readxl","ggpubr","ggrepel","peRspective","tagger","meta","aod","table1","peRspective")
lapply(packages, require, c = T)
devtools::install_github("eliocamp/tagger")

# Function to create "not in" operator
'%!in%' <- function(x,y)!('%in%'(x,y))

### Load data
mced.full <- read_excel("table_clean.xlsx")
dim(mced.full) # 63, with 60 papers 3 of which have 3 arms
source("functions_mcedtest.R")

## Data checks and harmonization 
str(mced.full) # some count variables considered as character -> change to numeric 

# 1. transform all count variables to numeric ####
mced.full <- mced.full %>% mutate(across(contains(c("Stage","deaths")),~as.numeric(.)))
# checks 
mced.full %>% select(contains(c("Stage", "deaths"))) %>% purrr::map(summary) #12 NA for stage 4 | 22 NA for all cause deaths : 12 papers did not have separate stage 4 info; 22 papers  did not have all cause mortality info 


#2. reclassify cancer types ####
table(mced.full$Cancer_type)
mced.full<-mced.full %>% mutate(cancer_type2=case_when(Cancer_type %in%c("Liver cancer","Nasopharyngeal carcinoma","Oral cancer")~"Other",
                                                       Cancer_type=="Breast cancer"~"Breast",Cancer_type=="Colorectal cancer"~"Colorectum",
                                                       Cancer_type=="Lung cancer"~"Lung",Cancer_type=="Ovarian cancer"~"Ovary",Cancer_type=="Prostate cancer"~"Prostate"))

# checks 
table(mced.full$Cancer_type,mced.full$cancer_type2)

#3. Define stage and cancer mortality reduction ####
mced.full <- mced.full %>%
  rowwise() %>%
  mutate(case_arm1      = sum(c(Stage_0_2_arm1, Stage_3_4_arm1), na.rm = T),
         case_arm2      = sum(c(Stage_0_2_arm2, Stage_3_4_arm2), na.rm = T),
       
         CumInc_Stage_0_2_arm1 = Stage_0_2_arm1/N_arm1, 
         CumInc_Stage_3_4_arm1 = Stage_3_4_arm1/N_arm1,
         CumInc_Stage_4_arm1 = Stage4_arm1/N_arm1,

         CumInc_Stage_0_2_arm2 = Stage_0_2_arm2/N_arm2, 
         CumInc_Stage_3_4_arm2 = Stage_3_4_arm2/N_arm2,
         CumInc_Stage_4_arm2 = Stage4_arm2/N_arm2,

         Cumstages_shift_arm1 = Stage_3_4_arm1/case_arm1,
         Cumstages_shift_arm2 = Stage_3_4_arm2/case_arm2,
         
         CumInc_CaMort_arm1 = Cancer_deaths_arm1/N_arm1,
         CumInc_CaMort_arm2 = Cancer_deaths_arm2/N_arm2,
         
         Reduction_CaMort_arm2vs1_perc = 100*((CumInc_CaMort_arm1 - CumInc_CaMort_arm2)/CumInc_CaMort_arm1),  # using control arm numbers as denominator
         Reduction_stageshift_arm2vs1_perc = 100*((Cumstages_shift_arm1 - Cumstages_shift_arm2)/Cumstages_shift_arm1),  # using control arm numbers as denominator
         Reduction_LateStage_arm2vs1_perc = 100*((CumInc_Stage_3_4_arm1 - CumInc_Stage_3_4_arm2)/CumInc_Stage_3_4_arm1),   # using control arm numbers as denominator
         Reduction_Stage4_arm2vs1_perc = 100*((CumInc_Stage_4_arm1 - CumInc_Stage_4_arm2)/CumInc_Stage_4_arm1),   # using control arm numbers as denominator
         
         total_N_arm2vs1 = N_arm1 + N_arm2) %>% 
  ungroup()

# checks 
summary(mced.full$Reduction_CaMort_arm2vs1_perc)
summary(mced.full$Reduction_LateStage_arm2vs1_perc)
summary(mced.full$Reduction_Stage4_arm2vs1_perc) # 15 NA 
mced.full %>% filter(is.na(Reduction_Stage4_arm2vs1_perc)) %>% 
  dplyr::select(`Citation (author, journal year)`,Stage4_arm1, Stage4_arm2, Reduction_Stage4_arm2vs1_perc )
# 12 of the 15 did not have separations of stage 4 vs stage 3-4 in the original dataset 
# 3 additional NA when calculating reduction in stage 4 incidence due to absence of stage 4 cancers (N of stage 4 cancer cases reported as 0 in the trial paper)
summary(mced.full$Reduction_stageshift_arm2vs1_perc)

## 3.a Add proportion tests comparing the 2 arms to the dataset (one for late-stage cancer, one for cancer mortality) ####
# compare cumulative incidence of late-stage cancers (stage 3+4) & using stage4 only as late stage cancer 
mced.full <- mced.full %>% rowwise() %>% 
  mutate(pval_stage3_4=round(prop.test(c(Stage_3_4_arm2, Stage_3_4_arm1),c(N_arm2,N_arm1))$p.value, 4), 
         pval_ca_mortality=round(prop.test(c(Cancer_deaths_arm2,Cancer_deaths_arm1),c(N_arm2,N_arm1))$p.value, 4),
         pval_stage3_4_p=round(prop.test(c(Stage_3_4_arm2, Stage_3_4_arm1),c(case_arm2, case_arm1))$p.value, 4),
         pval_stage4=tryCatch(round(prop.test(c(Stage4_arm2, Stage4_arm1),c(N_arm2, N_arm1))$p.value, 4),
                              error=function(e){return(NA)}),
         pval_stage4_p= tryCatch(round(prop.test(c(Stage4_arm2,Stage4_arm1),c(case_arm2, case_arm1))$p.value, 4),
                                 error=function(e){return(NA)})) %>% 
  ungroup()


# checks 
summary(mced.full$pval_ca_mortality)
summary(mced.full$pval_stage3_4)
summary(mced.full$pval_stage3_4_p)
summary(mced.full$pval_stage4)
summary(mced.full$pval_stage4)
summary(mced.full$pval_stage4_p)
summary(mced.full$pval_stage4_p)

## 3.b add labels for significance of effect ####
mced.full<-mced.full %>% mutate(sig_ca=ifelse(pval_ca_mortality<0.05&!is.na(pval_ca_mortality),1,0),
                                sig_stage3_4=ifelse(pval_stage3_4<0.05&!is.na(pval_stage3_4),1,0),
                                sig_stage4=case_when(pval_stage4<0.05&!is.na(pval_stage4)~1,
                                                     pval_stage4>=0.05&!is.na(pval_stage4)~0,
                                                     is.na(pval_stage4)~NA),
                                sig_stage3_4_p=ifelse(pval_stage3_4_p<0.05&!is.na(pval_stage3_4_p),1,0))

table(mced.full$sig_stage4,exclude = NULL)


#4. Create dummy variables to summarize trial size ####
mced.full<-mced.full %>% mutate(Trial_size=case_when(total_N_arm2vs1<=5000~"<=5k",
                                                     total_N_arm2vs1>5000&total_N_arm2vs1<=50000~"5k-50k",
                                                     total_N_arm2vs1>50000&total_N_arm2vs1<=100000~"50k-100k",
                                                     total_N_arm2vs1>100000~"100k+")) 
mced.full$Trial_size<-factor(mced.full$Trial_size,levels = c("<=5k","5k-50k","50k-100k","100k+"))
table(mced.full$Trial_size)

#5. Get data for main analysis ####
#!!!!!Data for main analysis (keep one paper per study)
mced.filtered<-mced.full %>% filter(primary_analysis=="yes")

#######################################################################################################
# Table 1 ####
#######################################################################################################

table1_unique<-table1(~ cancer_type2+Trial_size+factor(sig_stage3_4_p)+factor(sig_stage3_4)+factor(sig_ca)+factor(sig_stage4), 
                      data=mced.filtered,overall="Primary",na.is.category=F)

table1_multi<-table1(~ cancer_type2+Trial_size+factor(sig_stage3_4_p)+factor(sig_stage3_4)+factor(sig_ca)+factor(sig_stage4), 
                      data=mced.full,overall="Multiple",na.is.category=F)

quantile (mced.filtered$`Follow-up, years`,probs = seq(0, 1, 0.25), na.rm = FALSE)
quantile (mced.full$`Follow-up, years`,probs = seq(0, 1, 0.25), na.rm = FALSE)
#######################################################################################################
# Table2 ####
#######################################################################################################

function_sum_per <- function(x, n){
  sumN <- sum(x, na.rm = T)
  per <- round((sumN/n)*100, 1)
  all <- paste0(sumN, " (", per, ")" )
  return(all)
}
stage3_4red_all <- mced.filtered %>%
  dplyr::summarise(N=n(),
                   TP_05=function_sum_per(pval_stage3_4<0.05 & pval_ca_mortality<0.05,N),
                   TN_05=function_sum_per(pval_stage3_4>=0.05 & pval_ca_mortality>=0.05,N),
                   FP_05=function_sum_per(pval_stage3_4<0.05 & pval_ca_mortality>=0.05,N),
                   FN_05=function_sum_per(pval_stage3_4>=0.05 & pval_ca_mortality<0.05,N),
                   TP_1=function_sum_per(pval_stage3_4<0.1 & pval_ca_mortality<0.1,N),
                   TN_1=function_sum_per(pval_stage3_4>=0.1 & pval_ca_mortality>=0.1,N),
                   FP_1=function_sum_per(pval_stage3_4<0.1 & pval_ca_mortality>=0.1,N),
                   FN_1=function_sum_per(pval_stage3_4>=0.1 & pval_ca_mortality<0.1,N)) %>% 
  mutate(alt_outcome="reduction stage 3-4",cancer_type2="all")

stage3_4red_cancer <- mced.filtered %>% group_by(cancer_type2) %>% 
  dplyr::summarise(N=n(),
                   TP_05=as.character(sum(pval_stage3_4<0.05 & pval_ca_mortality<0.05)),
                   TN_05=as.character(sum(pval_stage3_4>=0.05 & pval_ca_mortality>=0.05)),
                   FP_05=as.character(sum(pval_stage3_4<0.05 & pval_ca_mortality>=0.05)),
                   FN_05=as.character(sum(pval_stage3_4>=0.05 & pval_ca_mortality<0.05)),
                   TP_1=as.character(sum(pval_stage3_4<0.1 & pval_ca_mortality<0.1)),
                   TN_1=as.character(sum(pval_stage3_4>=0.1 & pval_ca_mortality>=0.1)),
                   FP_1=as.character(sum(pval_stage3_4<0.1 & pval_ca_mortality>=0.1)),
                   FN_1=as.character(sum(pval_stage3_4>=0.1 & pval_ca_mortality<0.1)))%>% mutate(alt_outcome="reduction stage 3-4")


stage4red_cancer <- mced.filtered %>% 
  dplyr::summarise(N=sum(!is.na(pval_stage4)),
                   TP_05=function_sum_per(pval_stage4<0.05 & pval_ca_mortality<0.05, N),
                   TN_05=function_sum_per(pval_stage4>=0.05 & pval_ca_mortality>=0.05, N),
                   FP_05=function_sum_per(pval_stage4<0.05 & pval_ca_mortality>=0.05, N),
                   FN_05=function_sum_per(pval_stage4>=0.05 & pval_ca_mortality<0.05, N),
                   TP_1=function_sum_per(pval_stage4<0.1 & pval_ca_mortality<0.1, N),
                   TN_1=function_sum_per(pval_stage4>=0.1 & pval_ca_mortality>=0.1, N),
                   FP_1=function_sum_per(pval_stage4<0.1 & pval_ca_mortality>=0.1,N),
                   FN_1=function_sum_per(pval_stage4>=0.1 & pval_ca_mortality<0.1, N))%>% 
  mutate(alt_outcome="reduction 4",  cancer_type2="all")

stageshift_cancer <- mced.filtered %>%
  dplyr::summarise(N=n(),
                   TP_05=function_sum_per(pval_stage3_4_p<0.05 & pval_ca_mortality<0.05,N),
                   TN_05=function_sum_per(pval_stage3_4_p>=0.05 & pval_ca_mortality>=0.05,N),
                   FP_05=function_sum_per(pval_stage3_4_p<0.05 & pval_ca_mortality>=0.05,N),
                   FN_05=function_sum_per(pval_stage3_4_p>=0.05 & pval_ca_mortality<0.05,N),
                   TP_1=function_sum_per(pval_stage3_4_p<0.1 & pval_ca_mortality<0.1,N),
                   TN_1=function_sum_per(pval_stage3_4_p>=0.1 & pval_ca_mortality>=0.1,N),
                   FP_1=function_sum_per(pval_stage3_4_p<0.1 & pval_ca_mortality>=0.1,N),
                   FN_1=function_sum_per(pval_stage3_4_p>=0.1 & pval_ca_mortality<0.1,N))%>% 
  mutate(alt_outcome="stage shift",cancer_type2="all")



table2 <- bind_rows(stage3_4red_all, stage3_4red_cancer,stage4red_cancer,stageshift_cancer) %>% relocate(c("cancer_type2","alt_outcome"), .before = 1) 
write.csv(table2, paste0(path_to_save, "table2.csv"))

##################################################################################################
# Figure 1 ####
# main results, stage iii-iv and stage iv for correlation with cancer mortality reduction with no weighting
##################################################################################################
cancer_list<-unique(mced.filtered$cancer_type2)

data_figure1<-data_corrfigure(mced.filtered, cancer_list, "Reduction_LateStage_arm2vs1_perc")
  
#########here we have data for the slope confidence interval in the manuscript
data_figure1<-data_figure1 %>% mutate(across(c("Linear_slope","Linear_r2","Pearson_all"),
                                             ~ifelse(cancer_type2=="Other", "", .)))
                                      
mced.filtered$cancer_type2<-factor(mced.filtered$cancer_type2,levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other"))
data_figure1$cancer_type2<-factor(data_figure1$cancer_type2,levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other"))

data_segment<-data.frame(cancer_type2=c("Breast","Colorectum","Lung","Ovary","Prostate","Other"))
data_segment$cancer_type2<-factor(data_segment$cancer_type2,levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other"))

Figure1<-ggplot(data=mced.filtered, aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc)) + 
    geom_smooth(aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc),
                data=mced.filtered %>% filter(!cancer_type2 %in%c("Other") ),
                method="lm", se=T, fullrange=T, linewidth=0.3, linetype="dashed", color="black",fill = "gray88")+
    geom_point(aes(size=total_N_arm2vs1,colour=cancer_type2)) + 
    xlab("Reduction in stage III-IV cancer, percent") + ylab("Reduction in cancer mortality, percent") +
    facet_wrap(~cancer_type2)+
    scale_color_manual(values=c("#DF8F44FF","#00A1D5FF","#374E55FF","#B24745FF","#79AF97FF","#6A6599FF"))+
    geom_segment(aes(x = 0, y = -60, xend = 0, yend = 80),
                 arrow = arrow(length = unit(0.3, "cm")))+
    geom_segment(aes(x = -60, y = 0, xend = 80, yend = 0),
                 arrow = arrow(length = unit(0.3, "cm")))+
    theme_bw() +  
    coord_cartesian(xlim = c(-50, 75),ylim = c(-50, 75)) +
    guides(size="none") + 
    theme(axis.title = element_text(size=16), axis.text = element_text(size=16), 
          legend.position ="none",
          legend.background=element_rect(fill = alpha("white", 0.5)), 
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 15),strip.background = element_rect(colour="black",fill="white"))+
    geom_text(data = data_figure1,mapping = aes(x = 5,y = -30,label = Pearson_all, hjust=0))+
    geom_text(data = data_figure1,mapping = aes(x = 5,y = -40,label = Linear_slope, hjust=0 ))+
    geom_text(data = data_figure1,mapping = aes(x = 5,y = -50,label = Linear_r2, hjust=0 ))+
  geom_segment(data = data_segment %>% filter(cancer_type2=="Breast"), 
               aes(x = 3, y = -7, xend = 10, yend = -7),
               arrow = arrow(length = unit(0.1, "cm")))+
  geom_segment(data =  data_segment %>% filter(cancer_type2=="Breast"), 
               aes(x = 3, y = 3, xend = 3, yend =10 ),
               arrow = arrow(length = unit(0.1, "cm")))+
  geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
            mapping = aes(x = 13,y = -7,label = "Reduced stage III-IV cancer", hjust=0,fontface = "italic"),  size=3)+
  geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
          mapping = aes(x = 3,y = 42,label = "Reduced cancer mortality", angle=90,fontface = "italic"),  size=3)+
tag_facets(tag_levels  = c("A", "B","C","D","E","F"),tag_suffix = "")

Figure1
ggsave(Figure1,filename="Figure1_ABCDEF.pdf", path=path_to_save,device = cairo_pdf,
       height=8, width=13, units="in")

####Meta analysis 

result_meta<-metacor(cor=Pearson ,n=n_pub,studlab=cancer_type2,data=filter(data_figure1, cancer_type2!="Other"))

forest(result_meta,random = F)
  
pdf(paste0(path_to_save, "Figure1_meta",".pdf"), width = 9, height = 4)
meta::forest(result_meta,xlim=c(-1,1),random = F,
             leftcols =c("studlab", "n","effect", "ci", "w.common"),
             leftlabs =c("Cancer type","Studies, N","Correlation","95% CI","Weight"),
             xlab =c("Correlation (95% CI)"),
             layout = "JAMA",print.tau2 =T)
dev.off()
 

########Heterogeneity test
#Wald test is often used to determine if one or more predictor variables in a regression model are equal to zero.
#H0: Some set of predictor variables are all equal to zero.
#HA: Not all predictor variables in the set are equal to zero.
mced.filtered_noother<-mced.filtered %>% filter(cancer_type2!="Other")
table(mced.filtered_noother$cancer_type2)
model_cancer<-lm(Reduction_CaMort_arm2vs1_perc~Reduction_LateStage_arm2vs1_perc+Reduction_LateStage_arm2vs1_perc*cancer_type2,data =mced.filtered_noother )
summary(model_cancer)
wald.test(Sigma = vcov(model_cancer), b = coef(model_cancer), Terms = 7:10)
 #   Chi-squared test:
 # X2 = 15.6, df = 4, P(> X2) = 0.0036
 
 
##################################################################################################
# Figure 2 ####
# main results, stage iv for correlation with cancer mortality reduction 
##################################################################################################
with(filter(mced.filtered, !is.na(Reduction_Stage4_arm2vs1_perc)), table(cancer_type2))
# Breast cancer now only has 3 studies, the correlation confidence interval would not be precise

data_figure2<-data_corrfigure(mced.filtered, cancer_list, "Reduction_Stage4_arm2vs1_perc")


# given that we only have 2 publications for prostate , and we are not interested in the overall trend in other -> remove the results 
data_figure2 <- data_figure2 %>% mutate(across(c("Pearson_all","Linear_slope","Linear_r2"),
                                               ~ifelse(cancer_type2 %in% c("Prostate","Other"),"", .)))

data_figure2$cancer_type2<-factor(data_figure1$cancer_type2,levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other"))

Figure2<-ggplot(data=mced.filtered, aes(x=Reduction_Stage4_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc)) + 
   geom_smooth(aes(x=Reduction_Stage4_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc), data=mced.filtered %>% filter(!cancer_type2 %in%c("Prostate","Other") ),
               method="lm", se=T, fullrange=T, linewidth=0.3, linetype="dashed", color="black",fill = "gray88")+
   geom_point(aes(size=total_N_arm2vs1,colour=cancer_type2)) + 
   xlab("Reduction in stage IV cancer, percent") + ylab("Reduction in cancer mortality, percent") +
   facet_wrap(~cancer_type2)+
   scale_color_manual(values=c("#DF8F44FF","#00A1D5FF","#374E55FF","#B24745FF","#79AF97FF","#6A6599FF"))+
   geom_segment(aes(x = 0, y = -60, xend = 0, yend = 80),
                arrow = arrow(length = unit(0.3, "cm")))+
   geom_segment(aes(x = -60, y = 0, xend = 80, yend = 0),
                arrow = arrow(length = unit(0.3, "cm")))+
   theme_bw() +  
   coord_cartesian(xlim = c(-50, 75),ylim = c(-50, 75)) +
   guides(size="none") + 
   theme(axis.title = element_text(size=16), axis.text = element_text(size=16), 
         legend.position ="none",
         legend.background=element_rect(fill = alpha("white", 0.5)), 
         panel.grid.minor = element_blank(),
         strip.text = element_text(size = 15),strip.background = element_rect(colour="black",fill="white"))+
   geom_text(data = data_figure2,mapping = aes(x = 5,y = -20,label = Pearson_all, hjust=0))+
   geom_text(data = data_figure2,mapping = aes(x = 5,y = -30,label = Linear_slope, hjust=0 ))+
   geom_text(data = data_figure2,mapping = aes(x = 5,y = -40,label = Linear_r2, hjust=0 )) +
   geom_segment(data = data_segment %>% filter(cancer_type2=="Breast"), 
                aes(x = 3, y = -7, xend = 10, yend = -7),
                arrow = arrow(length = unit(0.1, "cm")))+
   geom_segment(data =  data_segment %>% filter(cancer_type2=="Breast"), 
                aes(x = 3, y = 3, xend = 3, yend =10 ),
                arrow = arrow(length = unit(0.1, "cm")))+
   geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
             mapping = aes(x = 13,y = -7,label = "Reduced stage IV cancer", hjust=0,fontface = "italic"),  size=3)+
   geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
             mapping = aes(x = 3,y = 42,label = "Reduced cancer mortality", angle=90,fontface = "italic"),  size=3) +
   tag_facets(tag_levels  = c("A", "B","C","D","E","F"),tag_suffix = "")
 
Figure2
ggsave(Figure2,filename="Figure2_ABCDEF.pdf", path=path_to_save,device = cairo_pdf,
        height=8, width=13, units="in")
 
#Meta analysis
result_meta<-metacor(cor=Pearson,n=n_pub,studlab=cancer_type2,data=filter(data_figure2, cancer_type2 %!in% c("Prostate","Other")))
forest(result_meta, random = F)
 
pdf(paste0(path_to_save, "Figure2_meta",".pdf"), width = 9, height = 4)
meta::forest(result_meta,xlim=c(-1,1),random = F,
             leftcols =c("studlab", "n","effect", "ci", "w.common"),
             leftlabs =c("Cancer type","Studies, N","Correlation","95% CI","Weight"),
             xlab =c("Correlation (95% CI)"),
             layout = "JAMA",print.tau2 =T)
dev.off()
 
 
##################################################################################################
##Sensitivity analysis
##################################################################################################

# sup figure 2 using the unique trials but use the last timepoint data.####
mced.full<-mced.full %>% mutate(year_pub= str_extract(`Citation (author, journal year)`, "\\d+"))
mced.full<-mced.full %>% mutate(year_pub=substr(year_pub,1,4))
 
mced.filtered.lastfup <- mced.full %>% 
  mutate(study_cancer_id_arm=paste0(`Study name`,"-","cancer_type2","-",Intervention), 
         cancer_type2=factor(cancer_type2, levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other"))) %>% 
   group_by(study_cancer_id_arm) %>%
   arrange(desc(year_pub)) %>%
   filter(row_number()==1) %>% 
  ungroup()

mced.filtered.lastfup <-  filter(mced.filtered.lastfup, !grepl( "ERSPC-Finland|ERSPC-Sweden|ERSPC-Rotterdam",`Study name` ))
#these will only be excluded from the "unique analysis" , but will participate in other sensitivity analysis and results

dim(mced.filtered.lastfup) #41
table(mced.filtered.lastfup$cancer_type2)

data_figure_sup2<-data_corrfigure(mced.filtered.lastfup, cancer_list, "Reduction_LateStage_arm2vs1_perc")

data_figure_sup2<-data_figure_sup2 %>%
  mutate(across(c("Linear_slope","Linear_r2","Pearson_all"),
                ~ifelse(cancer_type2=="Other","",.)),
         cancer_type2=factor(cancer_type2, levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other")))  
 

p_sup2<- ggplot(data=mced.filtered.lastfup, aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc)) + 
   geom_smooth(aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc), data=mced.filtered.lastfup %>% filter(cancer_type2!="Other"),
               method="lm", se=T, fullrange=T, linewidth=0.3, linetype="dashed", color="black",fill = "gray88")+
   geom_point(aes(size=total_N_arm2vs1,colour=cancer_type2)) + 
   xlab("Reduction in stage III-IV cancer, percent") + ylab("Reduction in cancer mortality, percent") +
   facet_wrap(~cancer_type2)+
   scale_color_manual(values=c("#DF8F44FF","#00A1D5FF","#374E55FF","#B24745FF","#79AF97FF","#6A6599FF"))+
   geom_segment(aes(x = 0, y = -60, xend = 0, yend = 80),
                arrow = arrow(length = unit(0.3, "cm")))+
   geom_segment(aes(x = -60, y = 0, xend = 80, yend = 0),
                arrow = arrow(length = unit(0.3, "cm")))+
   theme_bw() +  
   coord_cartesian(xlim = c(-50, 75),ylim = c(-50, 75)) +
   guides(size="none") + 
   theme(axis.title = element_text(size=16), axis.text = element_text(size=16), 
         legend.position ="none",
         legend.background=element_rect(fill = alpha("white", 0.5)), 
         panel.grid.minor = element_blank(),
         strip.text = element_text(size = 15),strip.background = element_rect(colour="black",fill="white"))+
   geom_text(data = data_figure_sup2,mapping = aes(x = 5,y = -20,label = Pearson_all, hjust=0))+
   geom_text(data = data_figure_sup2,mapping = aes(x = 5,y = -30,label = Linear_slope, hjust=0 ))+
   geom_text(data = data_figure_sup2,mapping = aes(x = 5,y = -40,label = Linear_r2, hjust=0 ))+
   geom_segment(data = data_segment %>% filter(cancer_type2=="Breast"), 
                aes(x = 3, y = -7, xend = 10, yend = -7),
                arrow = arrow(length = unit(0.1, "cm")))+
   geom_segment(data =  data_segment %>% filter(cancer_type2=="Breast"), 
                aes(x = 3, y = 3, xend = 3, yend =10 ),
                arrow = arrow(length = unit(0.1, "cm")))+
   geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
             mapping = aes(x = 13,y = -7,label = "Reduced stage III-IV cancer", hjust=0,fontface = "italic"),  size=3)+
   geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
             mapping = aes(x = 3,y = 42,label = "Reduced cancer mortality", angle=90,fontface = "italic"),  size=3)   

p_sup2
ggsave(p_sup2,filename="Figure_sup2.pdf", path=path_to_save,device = cairo_pdf,
       height=8, width=13, units="in")


#Meta result
result_meta<-metacor(cor=Pearson,n=n_pub,studlab=cancer_type2,data=filter(data_figure_sup2, cancer_type2!="Other"))
forest(result_meta, random = F)
 
pdf(paste0(path_to_save, "Figure_sup2_meta",".pdf"), width = 8, height = 4)
meta::forest(result_meta,xlim=c(-1,1),random = F,
             leftcols =c("studlab", "n","effect", "ci", "w.common"),
             leftlabs =c("Cancer type","Studies, N","Correlation","95% CI","Weight"),
             xlab =c("Correlation (95% CI)"),
             layout = "JAMA",print.tau2 =T)
dev.off()
 

############################################################
# sup figure 3 Using first data on stage and last data on mortality ####
# supplementary figure 2 using the unique trials but use the last timepoint data.
mced.filtered<-mced.filtered %>% mutate(year_pub= str_extract(`Citation (author, journal year)`, "\\d+"),
                                        year_pub=substr(year_pub,1,4))
 
mced.filtered <- mced.filtered %>% mutate(study_cancer_id_arm=paste0(`Study name`,"-","cancer_type2","-",Intervention),
                                          study_cancer_id_arm_year=paste0(`Study name`,"-","cancer_type2","-",Intervention,"-",year_pub)) 
 
mced.filtered.lastfup<- mced.filtered.lastfup %>% mutate(study_cancer_id_arm_year=paste0(`Study name`,"-","cancer_type2","-",Intervention,"-",year_pub)) 
unique(mced.filtered$study_cancer_id_arm_year)##41 fine!
unique(mced.filtered.lastfup$study_cancer_id_arm_year)##41 fine!
list1<-unique(mced.filtered.lastfup$study_cancer_id_arm)##41 fine!
list2<-unique(mced.filtered$study_cancer_id_arm)##41 fine!
list1[!list1 %in% list2]
list2[!list2 %in% list1]
# all same names, fine
 
data_same<-mced.filtered %>% filter(study_cancer_id_arm_year %in%mced.filtered.lastfup$study_cancer_id_arm_year)
data_same %>% dim()
#27 same, 14 need to update
mced.filtered.lastfup<-mced.filtered.lastfup %>% mutate(year_pub_later=year_pub)
mced.filtered.lastfup_add<-mced.filtered.lastfup %>% dplyr::select(study_cancer_id_arm,year_pub_later,
                                                                    Reduction_CaMort_arm2vs1_perc,pval_ca_mortality)
 
mced.filtered.first.last<-mced.filtered %>% dplyr::select(-c("Reduction_CaMort_arm2vs1_perc","pval_ca_mortality"))
 
 
mced.filtered.first.last<-merge(mced.filtered.first.last,mced.filtered.lastfup_add,by.x = "study_cancer_id_arm", 
                                 by.y = "study_cancer_id_arm")
dim(mced.filtered.first.last)
 
mced.filtered.first.last<-mced.filtered.first.last %>% ungroup()

data_figure_sup3<-data_corrfigure(mced.filtered.first.last, cancer_list, "Reduction_LateStage_arm2vs1_perc")


data_figure_sup3<-data_figure_sup3 %>% mutate(across(c("Linear_slope","Linear_r2","Pearson_all"), 
                                                     ~ifelse(cancer_type2=="Other", "",.)),
                                              cancer_type2=factor(cancer_type2,levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other")))
  
mced.filtered.first.last$cancer_type2<-factor(mced.filtered.first.last$cancer_type2,levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other"))
 
p_sup3<- ggplot(data=mced.filtered.first.last, aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc)) + 
   geom_smooth(aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc), data=mced.filtered.first.last %>% filter(cancer_type2!="Other"),
               method="lm", se=T, fullrange=T, linewidth=0.3, linetype="dashed", color="black",fill = "gray88")+
   geom_point(aes(size=total_N_arm2vs1,colour=cancer_type2)) + 
   xlab("Reduction in stage III-IV cancer, percent") + ylab("Reduction in cancer mortality, percent") +
   facet_wrap(~cancer_type2)+
   scale_color_manual(values=c("#DF8F44FF","#00A1D5FF","#374E55FF","#B24745FF","#79AF97FF","#6A6599FF"))+
   geom_segment(aes(x = 0, y = -60, xend = 0, yend = 80),
                arrow = arrow(length = unit(0.3, "cm")))+
   geom_segment(aes(x = -60, y = 0, xend = 80, yend = 0),
                arrow = arrow(length = unit(0.3, "cm")))+
   theme_bw() +  
   coord_cartesian(xlim = c(-50, 75),ylim = c(-50, 75)) +
   guides(size="none") + 
   theme(axis.title = element_text(size=16), axis.text = element_text(size=16), 
         legend.position ="none",
         legend.background=element_rect(fill = alpha("white", 0.5)), 
         panel.grid.minor = element_blank(),
         strip.text = element_text(size = 15),strip.background = element_rect(colour="black",fill="white"))+
   geom_text(data = data_figure_sup3,mapping = aes(x = 5,y = -20,label = Pearson_all, hjust=0))+
   geom_text(data = data_figure_sup3,mapping = aes(x = 5,y = -30,label = Linear_slope, hjust=0 ))+
   geom_text(data = data_figure_sup3,mapping = aes(x = 5,y = -40,label = Linear_r2, hjust=0 ))+
   geom_segment(data = data_segment %>% filter(cancer_type2=="Breast"), 
                aes(x = 3, y = -7, xend = 10, yend = -7),
                arrow = arrow(length = unit(0.1, "cm")))+
   geom_segment(data =  data_segment %>% filter(cancer_type2=="Breast"), 
                aes(x = 3, y = 3, xend = 3, yend =10 ),
                arrow = arrow(length = unit(0.1, "cm")))+
   geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
             mapping = aes(x = 13,y = -7,label = "Reduced stage III-IV cancer", hjust=0,fontface = "italic"),  size=3)+
   geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
             mapping = aes(x = 3,y = 42,label = "Reduced cancer mortality", angle=90,fontface = "italic"),  size=3)   

p_sup3
ggsave(p_sup3,filename="Figure_sup3.pdf", path=path_to_save,device = cairo_pdf,
       height=8, width=13, units="in")


###Meta 
result_meta<-metacor(cor=Pearson,n=n_pub,studlab=cancer_type2,data=filter(data_figure_sup3, cancer_type2!="Other"))
forest(result_meta, random = F)
 
pdf(paste0(path_to_save, "Figure_sup3_meta",".pdf"), width = 8, height = 4)
meta::forest(result_meta,xlim=c(-1,1),random = F,
              leftcols =c("studlab", "n","effect", "ci", "w.common"),
              leftlabs =c("Cancer type","Studies, N","Correlation","95% CI","Weight"),
              xlab =c("Correlation (95% CI)"),
              layout = "JAMA",print.tau2 =T)
dev.off()
 
 
# sup figure 4 using the multiple trials  ####
data_figure_sup4<-data_corrfigure(mced.full, cancer_list, "Reduction_LateStage_arm2vs1_perc")

data_figure_sup4<-data_figure_sup4 %>% mutate(across(c("Linear_slope","Linear_r2","Pearson_all"),
                                                     ~ifelse(cancer_type2=="Other","",.)),
                                              cancer_type2=factor(cancer_type2,levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other")))
  
mced.full$cancer_type2<-factor(mced.full$cancer_type2,levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other"))

  
p_sup4<-  ggplot(data=mced.full, aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc)) + 
    geom_smooth(aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc), data=mced.full %>% filter(cancer_type2!="Other"),
                method="lm", se=T, fullrange=T, linewidth=0.3, linetype="dashed", color="black",fill = "gray88")+
    geom_point(aes(size=total_N_arm2vs1,colour=cancer_type2)) + 
    xlab("Reduction in stage III-IV cancer, percent") + ylab("Reduction in cancer mortality, percent") +
    facet_wrap(~cancer_type2)+
    scale_color_manual(values=c("#DF8F44FF","#00A1D5FF","#374E55FF","#B24745FF","#79AF97FF","#6A6599FF"))+
    geom_segment(aes(x = 0, y = -60, xend = 0, yend = 80),
                 arrow = arrow(length = unit(0.3, "cm")))+
    geom_segment(aes(x = -60, y = 0, xend = 80, yend = 0),
                 arrow = arrow(length = unit(0.3, "cm")))+
    theme_bw() +  
    coord_cartesian(xlim = c(-50, 75),ylim = c(-50, 75)) +
    guides(size="none") + 
    theme(axis.title = element_text(size=16), axis.text = element_text(size=16), 
          legend.position ="none",
          legend.background=element_rect(fill = alpha("white", 0.5)), 
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 15),strip.background = element_rect(colour="black",fill="white"))+
    geom_text(data = data_figure_sup4,mapping = aes(x = 5,y = -20,label = Pearson_all, hjust=0))+
    geom_text(data = data_figure_sup4,mapping = aes(x = 5,y = -30,label = Linear_slope, hjust=0 ))+
    geom_text(data = data_figure_sup4,mapping = aes(x = 5,y = -40,label = Linear_r2, hjust=0 ))+
    geom_segment(data = data_segment %>% filter(cancer_type2=="Breast"), 
                 aes(x = 3, y = -7, xend = 10, yend = -7),
                 arrow = arrow(length = unit(0.1, "cm")))+
    geom_segment(data =  data_segment %>% filter(cancer_type2=="Breast"), 
                 aes(x = 3, y = 3, xend = 3, yend =10 ),
                 arrow = arrow(length = unit(0.1, "cm")))+
    geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
              mapping = aes(x = 13,y = -7,label = "Reduced stage III-IV cancer", hjust=0,fontface = "italic"),  size=3)+
    geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
              mapping = aes(x = 3,y = 42,label = "Reduced cancer mortality", angle=90,fontface = "italic"),  size=3)  
p_sup4
ggsave(p_sup4,filename="Figure_sup4.pdf", path=path_to_save,device = cairo_pdf,
       height=8, width=13, units="in")


#Meta result
result_meta<-metacor(cor=Pearson,n=n_pub,studlab=cancer_type2,data=filter(data_figure_sup4, cancer_type2!="Other"))
forest(result_meta,  random = F)
  
pdf(paste0(path_to_save, "Figure_sup4_meta",".pdf"), width = 8, height = 4)
meta::forest(result_meta,xlim=c(-1,1),random = F,
             leftcols =c("studlab", "n","effect", "ci", "w.common"),
             leftlabs =c("Cancer type","Studies, N","Correlation","95% CI","Weight"),
             xlab =c("Correlation (95% CI)"),
             layout = "JAMA",print.tau2 =T)
dev.off()
  

# sup Figure 5, force linear regression model intercept as 0 ####
data_nointercept<-data.frame(cancer_type2=NA, Linear_x=NA,Linear_r2=NA)

for (i in c(1:length(cancer_list)) ) {
  data_nointercept[i,1]<-cancer_list[i]
  model_lm <- lm(Reduction_CaMort_arm2vs1_perc ~ Reduction_LateStage_arm2vs1_perc+0, data=mced.filtered[mced.filtered$cancer_type2==cancer_list[i],])
  data_nointercept[i,2]<-summary(model_lm)[["coefficients"]][,"Estimate"]
  data_nointercept[i,3]<-summary(model_lm)[["r.squared"]]
}

data_nointercept
data_nointercept<-data_nointercept %>%
  mutate(Linear_slope=paste0("Linear regression y=",specify_decimal(Linear_x,2),"x"),
         Linear_r2=paste0("R-squared=",specify_decimal(Linear_r2,2)),
         across(c("Linear_slope","Linear_r2"),~ifelse(cancer_type2=="Other","",.)),
         cancer_type2=factor(cancer_type2,levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other")))

mced.filtered$cancer_type2<-factor(mced.filtered$cancer_type2,levels = c("Breast","Colorectum","Lung","Ovary","Prostate","Other"))

p_sup5<-ggplot(data=mced.filtered, aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc)) + 
  geom_point(aes(size=total_N_arm2vs1,colour=cancer_type2)) + 
  geom_smooth(formula = y~x+0,method="lm",data = mced.filtered %>% filter(!cancer_type2 %in%c("Other")),
              se=T, fullrange=T, linewidth=0.3, linetype="dashed", color="black",fill = "gray88")+
  xlab("Reduction in stage III-IV cancer, percent") + ylab("Reduction in cancer mortality, percent") +
  facet_wrap(~cancer_type2)+
  scale_color_manual(values=c("#DF8F44FF","#00A1D5FF","#374E55FF","#B24745FF","#79AF97FF","#6A6599FF"))+
  geom_segment(aes(x = 0, y = -60, xend = 0, yend = 80),
               arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = -60, y = 0, xend = 80, yend = 0),
               arrow = arrow(length = unit(0.3, "cm")))+
  theme_bw() +  
  coord_cartesian(xlim = c(-50, 75),ylim = c(-50, 75)) +
  guides(size="none") + 
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16), 
        legend.position ="none",
        legend.background=element_rect(fill = alpha("white", 0.5)), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),strip.background = element_rect(colour="black",fill="white"))+
  geom_text(data = data_nointercept,mapping = aes(x = 5,y = -30,label = Linear_slope, hjust=0 ))+
  geom_text(data = data_nointercept,mapping = aes(x = 5,y = -40,label = Linear_r2, hjust=0 )) +
  geom_segment(data = data_segment %>% filter(cancer_type2=="Breast"), 
               aes(x = 3, y = -7, xend = 10, yend = -7),
               arrow = arrow(length = unit(0.1, "cm")))+
  geom_segment(data =  data_segment %>% filter(cancer_type2=="Breast"), 
               aes(x = 3, y = 3, xend = 3, yend =10 ),
               arrow = arrow(length = unit(0.1, "cm")))+
  geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
            mapping = aes(x = 13,y = -7,label = "Reduced stage III-IV cancer", hjust=0,fontface = "italic"),  size=3)+
  geom_text(data =  data_segment %>% filter(cancer_type2=="Breast"), 
            mapping = aes(x = 3,y = 42,label = "Reduced cancer mortality", angle=90,fontface = "italic"),  size=3)

p_sup5
ggsave(p_sup5,filename="Figure_sup5.pdf", path=path_to_save,device = cairo_pdf,
       height=8, width=13, units="in")

# Supplementary Figure 6,CRC by screening tools
mced.filtered<-mced.filtered %>% 
  mutate(tools=ifelse(cancer_type2=="Colorectum",ifelse(Intervention%in%c("FOBT","1. Annual FOBT","2. Biennial FOBT"),"FOBT","Endoscopy" ),NA_character_))
table(mced.filtered$tools,exclude = NULL)
#Endoscopy      FOBT 
#4         7 

tool_list<-unique(mced.filtered$tools)[!is.na(unique(mced.filtered$tools))]

data_figure_crc <- data.frame(NULL)
for(i in c("Endoscopy", "FOBT")){
  data_crc_tool <- mced.filtered %>% filter(cancer_type2=="Colorectum", 
                                            tools==i)
  data_figure_crc_interim<-data_corrfigure(data_crc_tool, c("Colorectum"), "Reduction_LateStage_arm2vs1_perc") %>% 
    mutate(tools=i)
  data_figure_crc <- bind_rows(data_figure_crc,data_figure_crc_interim)
}
data_figure_crc

p_sup_crc<-ggplot(data=mced.filtered %>% filter(cancer_type2=="Colorectum"), aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc)) + 
  geom_smooth(aes(x=Reduction_LateStage_arm2vs1_perc, y=Reduction_CaMort_arm2vs1_perc), 
              method="lm", se=T, fullrange=T, linewidth=0.3, linetype="dashed", color="black",fill = "gray88")+
  geom_point(aes(size=total_N_arm2vs1,colour=tools)) + 
  xlab("Reduction in stage III-IV cancer, percent") + ylab("Reduction in cancer mortality, percent") +
  facet_wrap(~tools)+
  scale_color_manual(values=c("#4DBBD5FF","#3C5488FF"))+
  geom_segment(aes(x = 0, y = -60, xend = 0, yend = 80),
               arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = -60, y = 0, xend = 80, yend = 0),
               arrow = arrow(length = unit(0.3, "cm")))+
  theme_bw() +  
  coord_cartesian(xlim = c(-50, 75),ylim = c(-50, 75)) +
  guides(size="none") + 
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16), 
        legend.position ="none",
        legend.background=element_rect(fill = alpha("white", 0.5)), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),strip.background = element_rect(colour="black",fill="white"))+
  geom_text(data = data_figure_crc,mapping = aes(x = 5,y = -20,label = Pearson_all, hjust=0))+
  geom_text(data = data_figure_crc,mapping = aes(x = 5,y = -30,label = Linear_slope, hjust=0 ))+
  geom_text(data = data_figure_crc,mapping = aes(x = 5,y = -40,label = Linear_r2, hjust=0 )) +
  geom_segment(data = data.frame(tools="Endoscopy"), 
               aes(x = 3, y = -7, xend = 10, yend = -7),
               arrow = arrow(length = unit(0.1, "cm")))+
  geom_segment(data =  data.frame(tools="Endoscopy"), 
               aes(x = 3, y = 3, xend = 3, yend =10 ),
               arrow = arrow(length = unit(0.1, "cm")))+
  geom_text(data =  data.frame(tools="Endoscopy"), 
            mapping = aes(x = 13,y = -7,label = "Reduced stage III-IV cancer", hjust=0,fontface = "italic"),  size=3)+
  geom_text(data =  data.frame(tools="Endoscopy"), 
            mapping = aes(x = 3,y = 35,label = "Reduced cancer mortality", angle=90,fontface = "italic"),  size=3)   

p_sup_crc
ggsave(p_sup_crc,filename="Figure_sup6_crctools.pdf", path=path_to_save,device = cairo_pdf,
       height=5, width=9, units="in")

  
#####################################################################
#Side analysis not in the tables/plots
#Slope confidence interval in the manuscript, see data_figure1 results

# Manuscript
# total participants
mced.filtered %>% summarise(sum(total_N_arm2vs1))
#2796398  
# Remove two duplicated control arms numbers
mced.filtered %>% filter(`Citation (author, journal year)` %in% c("Mandel et al, N Engl J Med 199314**","Jacobs et al, Lancet 201643**")) %>% 
  dplyr::select("Citation (author, journal year)",year_pub,Intervention,N_arm2,N_arm1,total_N_arm2vs1)
#101299          
#15394           
#2796398  -101299 -15394  
#2679705
  
##########heterogeneity test for the slope by using the last time point

mced.filtered_noother_last<-mced.filtered.lastfup %>% filter(cancer_type2!="Other")
table(mced.filtered_noother_last$cancer_type2)
model_cancer_last<-lm(Reduction_CaMort_arm2vs1_perc~Reduction_LateStage_arm2vs1_perc+Reduction_LateStage_arm2vs1_perc*cancer_type2,
                      data =mced.filtered_noother_last )
summary(model_cancer_last)
wald.test(Sigma = vcov(model_cancer_last), b = coef(model_cancer_last), Terms = 7:10)
#Chi-squared test:
#  X2 = 13.9, df = 4, P(> X2) = 0.0075

#############################################################################
#############################################################################
###############################END HERE######################################
###############################END HERE######################################
###############################END HERE######################################
#############################################################################
#############################################################################





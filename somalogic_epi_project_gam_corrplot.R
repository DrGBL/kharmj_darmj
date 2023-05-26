#see here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7265853/

#also see here https://www.nature.com/articles/s41586-020-2700-3

#also look at pathfindR
#https://ccforum.biomedcentral.com/articles/10.1186/s13054-021-03503-x#Sec2
#https://cran.r-project.org/web/packages/pathfindR/index.html

#need to add IFNalpha2

library(tidyverse)
library(ggpubr)
library(grid)
library(ggthemes)
library(quantreg)
library(mgcv)
library(gratia)
#library(Cairo)
#library(visreg)
library(ggcorrplot)

setwd("/scratch/richards/guillaume.butler-laporte/somalogic_epi")

time_cutoff<-30

soma_raw<-read_csv("infe_417-soma_data=normalized-nat_log_tranf=FALSE-standardize=FALSE-remove_outliers=FALSE.csv") %>%
  filter(Days_symptom_update3<=time_cutoff) %>%
  mutate(A3=factor(A3,
                   labels=c("Control",
                            "Case"))) %>%
  mutate(sex=factor(sex)) %>%
  mutate(anonymized_patient_id=factor(anonymized_patient_id))
soma<-soma_raw
soma[,43:ncol(soma)]<-scale(log(soma_raw[,43:ncol(soma_raw)]))

soma_outlier<-scale(log(soma_raw[,43:ncol(soma_raw)]))
soma_outlier[which(soma_outlier>3 | soma_outlier<(-3))]<-NA
soma_outlier[which(soma_outlier<=3 & soma_outlier>=(-3))]<-1

soma[,43:ncol(soma)]<-soma[,43:ncol(soma)]*soma_outlier


######functions######
gam_fxn_full_only<-function(protein, name){
  gam_full<-soma %>%
  gam(as.formula(paste(protein,"~s(Days_symptom_update3, by=A3, bs=\"tp\")+
        s(Days_symptom_update3, by=sex, bs=\"tp\")+
        s(age_at_diagnosis, by=sex, bs=\"tp\")+
        s(age_at_diagnosis, by=A3, bs=\"tp\")+
        ti(Days_symptom_update3, age_at_diagnosis, by=sex, bs=\"tp\")+
        ti(Days_symptom_update3, age_at_diagnosis, by=A3, bs=\"tp\")+
        A3*sex")),
      data=.,
      method = "REML")
  return(gam_full)
}

gam_fxn_null_A3_only<-function(protein, name){
  gam_null<-soma %>%
    gam(as.formula(paste(protein,"~s(Days_symptom_update3, by=sex, bs=\"tp\")+
        s(age_at_diagnosis, by=sex, bs=\"tp\")+
        ti(Days_symptom_update3, age_at_diagnosis, by=sex, bs=\"tp\")+
        sex")),
        data=.,
        method = "REML")
  return(gam_null)
}

gam_fxn<-function(protein, name){
  loess_plot<-soma %>%
    filter(Days_symptom_update3<=time_cutoff) %>%
    mutate(A3=factor(A3,
                     labels=c("Control",
                              "Case"))) %>%
    ggplot(aes_string(x="Days_symptom_update3", 
                      y=protein,
                      group="A3",
                      colour="A3")) +
    geom_smooth() +
    xlab("Days since symptoms onset")+
    ylab(name)+
    scale_x_continuous(breaks=seq(0,30,5))+
    scale_color_manual(values=c("dodgerblue4", "red3"))+
    theme_bw()
  
  gam_full<-soma %>%
    filter(Days_symptom_update3<=time_cutoff) %>%
    mutate(A3=factor(A3,
                     labels=c("Control",
                              "Case"))) %>%
    mutate(sex=factor(sex)) %>%
    mutate(anonymized_patient_id=factor(anonymized_patient_id)) %>%
    gam(as.formula(paste(protein,"~s(Days_symptom_update3, by=A3, bs=\"tp\")+
        s(Days_symptom_update3, by=sex, bs=\"tp\")+
        s(age_at_diagnosis, by=sex, bs=\"tp\")+
        s(age_at_diagnosis, by=A3, bs=\"tp\")+
        ti(Days_symptom_update3, age_at_diagnosis, by=sex, bs=\"tp\")+
        ti(Days_symptom_update3, age_at_diagnosis, by=A3, bs=\"tp\")+
        A3*sex")),
        data=.,
        method = "REML")
  
  gam_null_outcome<-soma %>%
    filter(Days_symptom_update3<=time_cutoff) %>%
    mutate(A3=factor(A3,
                     labels=c("Control",
                              "Case"))) %>%
    mutate(sex=factor(sex)) %>%
    mutate(anonymized_patient_id=factor(anonymized_patient_id)) %>%
    gam(as.formula(paste(protein,"~s(Days_symptom_update3, by=sex, bs=\"tp\")+
        s(age_at_diagnosis, by=sex, bs=\"tp\")+
        ti(Days_symptom_update3, age_at_diagnosis, by=sex, bs=\"tp\")+
        sex")),
        data=.,
        method = "REML")
  
  gam_null_sex<-soma %>%
    filter(Days_symptom_update3<=time_cutoff) %>%
    mutate(A3=factor(A3,
                     labels=c("Control",
                              "Case"))) %>%
    mutate(sex=factor(sex)) %>%
    mutate(anonymized_patient_id=factor(anonymized_patient_id)) %>%
    gam(as.formula(paste(protein,"~s(Days_symptom_update3, by=A3, bs=\"tp\")+
        s(age_at_diagnosis, by=A3, bs=\"tp\")+
        ti(Days_symptom_update3, age_at_diagnosis, by=A3, bs=\"tp\")+
        A3")),
        data=.,
        method = "REML")
  
  #mean predict
  gam_plot_M<-predict_graph(data=test_data_M, model=gam_full, name=name)
  gam_plot_F<-predict_graph(data=test_data_F, model=gam_full, name=name)
  
  gam_appraisal<-appraise(gam_full)
  gam_anova_outcome<-anova(gam_null_outcome, gam_full, test="Chisq")
  gam_anova_sex<-anova(gam_null_sex, gam_full, test="Chisq")
  
  return(list(loess_plot, gam_full, gam_plot_M, gam_plot_F, gam_appraisal, gam_anova_outcome, gam_anova_sex))
}

predict_graph<-function(data, model, name){
  fits<-data %>%
    predict(model, newdata=., type="response", se=TRUE)
  
  predicts <- data.frame(data, fits) %>% 
    mutate(lower = fit - 1.96*se.fit,
           upper = fit + 1.96*se.fit)
  
  result_plot <- ggplot(aes(x=Days_symptom_update3,
                            y=fit,
                            group=A3,
                            colour=A3), 
                        data=predicts) +
    geom_ribbon(aes(ymin = lower, ymax=upper, group=A3, colour=A3, fill = A3),
                alpha=0.4) +
    geom_line(color="black")+    
    xlab("Days since symptoms onset")+
    ylab(name)+
    scale_x_continuous(breaks=seq(0,30,5))+
    scale_color_manual(values=c("dodgerblue4", "red3"))+
    scale_fill_manual(values=c("dodgerblue4", "red3"))+
    ylim(c(-3,3))+
    theme_bw()
  
  return(result_plot)
}

plot_output<-function(){
  #unadjusted
  unadj_plot<-ggarrange(il1a[[1]], 
                        il1b[[1]],
                        il2[[1]],
                        il4[[1]], 
                        il6[[1]],
                        il7[[1]],
                        il10[[1]],
                        il12[[1]],
                        il13[[1]],
                        il17a[[1]],
                        il17b[[1]],
                        il17c[[1]],
                        il17d[[1]],
                        il17f[[1]],
                        gcsf[[1]],
                        gmcsf[[1]],
                        mcsf[[1]],
                        tnfalpha[[1]],
                        ifngamma[[1]],
                        ip10[[1]],
                        mcp10[[1]],
                        labels = LETTERS[1:21],
                        ncol = 6, 
                        nrow = 4,
                        common.legend = TRUE, 
                        legend = "bottom") %>%
    annotate_figure(.,
                    fig.lab = "Unadjusted Cytokines Levels",
                    fig.lab.face = "bold",
                    fig.lab.pos="bottom.left"
    )
  
  #adjusted male
  i<-3
  adj_plot_male<-ggarrange(il1a[[i]], 
                           il1b[[i]],
                           il2[[i]],
                           il4[[i]], 
                           il6[[i]],
                           il7[[i]],
                           il10[[i]],
                           il12[[i]],
                           il13[[i]],
                           il17a[[i]],
                           il17b[[i]],
                           il17c[[i]],
                           il17d[[i]],
                           il17f[[i]],
                           gcsf[[i]],
                           gmcsf[[i]],
                           mcsf[[i]],
                           tnfalpha[[i]],
                           ifngamma[[i]],
                           ip10[[i]],
                           mcp10[[i]],
                           labels = LETTERS[1:21],
                           ncol = 6, 
                           nrow = 4,
                           common.legend = TRUE, 
                           legend = "bottom") %>%
    annotate_figure(.,
                    fig.lab = "Predicted Cytokines Levels, 65yo M",
                    fig.lab.face = "bold",
                    fig.lab.pos="bottom.left"
    )
  
  #adjusted female
  i<-4
  adj_plot_female<-ggarrange(il1a[[i]], 
                             il1b[[i]],
                             il2[[i]],
                             il4[[i]], 
                             il6[[i]],
                             il7[[i]],
                             il10[[i]],
                             il12[[i]],
                             il13[[i]],
                             il17a[[i]],
                             il17b[[i]],
                             il17c[[i]],
                             il17d[[i]],
                             il17f[[i]],
                             gcsf[[i]],
                             gmcsf[[i]],
                             mcsf[[i]],
                             tnfalpha[[i]],
                             ifngamma[[i]],
                             ip10[[i]],
                             mcp10[[i]],
                             labels = LETTERS[1:21],
                             ncol = 6, 
                             nrow = 4,
                             common.legend = TRUE, 
                             legend = "bottom") %>%
    annotate_figure(.,
                    fig.lab = "Predicted Cytokines Levels, 65yo F", 
                    fig.lab.face = "bold",
                    fig.lab.pos="bottom.left"
    )
  
  return(list(unadj_plot,adj_plot_male,adj_plot_female))
}

######run the full analyses #####

protein_list<-colnames(soma)[43:ncol(soma)]
protein_df<-data.frame(Protein=protein_list,
                       p_values=rep(NA_real_, length(protein_list)))
for(i in 1:length(protein_list)){
#for(i in 1:100){
  protein_df[i,2]<-anova(gam_fxn_null_A3_only(protein_df$Protein[i]), 
                         gam_fxn_full_only(protein_df$Protein[i]), 
                         test="Chisq")$`Pr(>Chi)`[2]
}

saveRDS(protein_df, file = "full_gam_p_values.rds")

#######correlation plots between cytokines######
corr <- soma %>%
  filter(Days_symptom_update3<=time_cutoff) %>%
  dplyr::select(c(IL1A.4851.25,
                  IL1B.3037.62,
                  IL2.3070.1,
                  IL4.2906.55,
                  IL6.4673.13,
                  IL7.4140.3,
                  IL10.2773.50,
                  IL12A.IL12B.10367.62,
                  IL13.3072.4,
                  IL17A.9170.24,
                  IL17B.3499.77,
                  IL17C.9255.5,
                  IL17D.4136.40,
                  IL17F.14026.24,
                  CSF3.4840.73,
                  CSF2.4697.59,
                  CSF1.3738.54,
                  TNF.5936.53,
                  IFNG.15346.31,
                  CXCL10.4141.79,
                  CCL2.2578.67)) %>%
  stats::cor(x=.,use="pairwise.complete.obs",
             method=c("spearman"))

colnames(corr) <- c("IL-1\u03b1",
                    "IL-1\u03b2",
                    "IL-2",
                    "IL-4",
                    "IL-6",
                    "IL-7",
                    "IL-10",
                    "IL-12",
                    "IL-13",
                    "IL-17A",
                    "IL-17B",
                    "IL-17C",
                    "IL-17D",
                    "IL-17F",
                    "G-CSF",
                    "GM-CSF",
                    "M-CSF",
                    "TNF-\u03b1",
                    "IFN-\u1d5e",
                    "CXCL10 / IP-10",
                    "CCL2 / MCP-1")
rownames(corr) <- colnames(corr)

ggcorrplot(corr, hc.order = TRUE, ggtheme = ggplot2::theme_bw(), colors= c("dodgerblue4","white","red3")) %>%
  ggsave(filename="correlation_cytokines.pdf", device=cairo_pdf, height = 8.5, width = 8.5,  units = "in", dpi = 600)



#######test data with mean age and all females######
sex_pred<-"M"
mean_age<-mean(soma$age_at_diagnosis)
quartiles_age<-quantile(soma$age_at_diagnosis, c(0.25, 0.75))
test_data_M<-soma %>%
  filter(Days_symptom_update3<=time_cutoff) %>%
  dplyr::select(A3, Days_symptom_update3) %>%
  mutate(A3=factor(A3, 
                   labels=c("Control",
                            "Case"))) %>%
  mutate(sex=factor(rep("M", nrow(.)))) %>%
  mutate(age_at_diagnosis=rep(mean_age, nrow(.)))

test_data_F<-soma %>%
  filter(Days_symptom_update3<=time_cutoff) %>%
  dplyr::select(A3, Days_symptom_update3) %>%
  mutate(A3=factor(A3, 
                   labels=c("Control",
                            "Case"))) %>%
  mutate(sex=factor(rep("F", nrow(.)))) %>%
  mutate(age_at_diagnosis=rep(mean_age, nrow(.)))


######run the analyses ######
#IL1a
il1a<-gam_fxn(protein="IL1A.4851.25", name="IL-1\u03b1") 

#IL1b
il1b<-gam_fxn(protein="IL1B.3037.62", name="IL-1\u03b2") 

#IL2
il2<-gam_fxn(protein="IL2.3070.1", name="IL-2") 

#IL4
il4<-gam_fxn(protein="IL4.2906.55", name="IL-4") 

#IL6
il6<-gam_fxn(protein="IL6.4673.13", name="IL-6")

#IL7
il7<-gam_fxn(protein="IL7.4140.3", name="IL-7")

#IL10
il10<-gam_fxn(protein="IL10.2773.50", name="IL-10")

#IL12
il12<-gam_fxn(protein="IL12A.IL12B.10367.62", name="IL-12")

#IL13
il13<-gam_fxn(protein="IL13.3072.4", name="IL-13")

#IL17a
il17a<-gam_fxn(protein="IL17A.9170.24", name="IL-17A")

#IL17B
il17b<-gam_fxn(protein="IL17B.3499.77", name="IL-17B")

#IL17C
il17c<-gam_fxn(protein="IL17C.9255.5", name="IL-17C")

#IL17D
il17d<-gam_fxn(protein="IL17D.4136.40", name="IL-17D")

#IL17F
il17f<-gam_fxn(protein="IL17F.14026.24", name="IL-17F")

#G-CSF
gcsf<-gam_fxn(protein="CSF3.4840.73", name="G-CSF")

#GM-CSF
gmcsf<-gam_fxn(protein="CSF2.4697.59", name="GM-CSF")

#M-CSF
mcsf<-gam_fxn(protein="CSF1.3738.54", name="M-CSF")

#TNF alpha
tnfalpha<-gam_fxn(protein="TNF.5936.53", name="TNF-\u03b1")

#Ifn gamma
ifngamma<-gam_fxn(protein="IFNG.15346.31", name="IFN-\u1d5e")

#IP-10
ip10<-gam_fxn(protein="CXCL10.4141.79", name="CXCL10 / IP-10")

#MCP-1
mcp10<-gam_fxn(protein="CCL2.2578.67", name="CCL2 / MCP-1")

###### plot things and p-values#####
fig_final<-plot_output()
ggsave(filename="unadjusted_cytokines.pdf", plot = fig_final[[1]], device=cairo_pdf, height = 8.5, width = 11,  units = "in", dpi = 600)
ggsave(filename="predicted_cytokines_65M.pdf", plot = fig_final[[2]], device=cairo_pdf, height = 8.5, width = 11,  units = "in", dpi = 600)
ggsave(filename="predicted_cytokines_65F.pdf", plot = fig_final[[3]], device=cairo_pdf, height = 8.5, width = 11,  units = "in", dpi = 600)


#p-values
p_df<-data.frame(Cytokine=c("il1a", 
                            "il1b",
                            "il2",
                            "il4", 
                            "il6",
                            "il7",
                            "il10",
                            "il12",
                            "il13",
                            "il17a",
                            "il17b",
                            "il17c",
                            "il17d",
                            "il17f",
                            "gcsf",
                            "gmcsf",
                            "mcsf",
                            "tnfalpha",
                            "ifngamma",
                            "ip10",
                            "mcp10"),
                 p_value_A3=c(il1a[[6]]$`Pr(>Chi)`[2], 
                           il1b[[6]]$`Pr(>Chi)`[2],
                           il2[[6]]$`Pr(>Chi)`[2],
                           il4[[6]]$`Pr(>Chi)`[2], 
                           il6[[6]]$`Pr(>Chi)`[2],
                           il7[[6]]$`Pr(>Chi)`[2],
                           il10[[6]]$`Pr(>Chi)`[2],
                           il12[[6]]$`Pr(>Chi)`[2],
                           il13[[6]]$`Pr(>Chi)`[2],
                           il17a[[6]]$`Pr(>Chi)`[2],
                           il17b[[6]]$`Pr(>Chi)`[2],
                           il17c[[6]]$`Pr(>Chi)`[2],
                           il17d[[6]]$`Pr(>Chi)`[2],
                           il17f[[6]]$`Pr(>Chi)`[2],
                           gcsf[[6]]$`Pr(>Chi)`[2],
                           gmcsf[[6]]$`Pr(>Chi)`[2],
                           mcsf[[6]]$`Pr(>Chi)`[2],
                           tnfalpha[[6]]$`Pr(>Chi)`[2],
                           ifngamma[[6]]$`Pr(>Chi)`[2],
                           ip10[[6]]$`Pr(>Chi)`[2],
                           mcp10[[6]]$`Pr(>Chi)`[2]),
                 p_value_sex=c(il1a[[7]]$`Pr(>Chi)`[2], 
                              il1b[[7]]$`Pr(>Chi)`[2],
                              il2[[7]]$`Pr(>Chi)`[2],
                              il4[[7]]$`Pr(>Chi)`[2], 
                              il6[[7]]$`Pr(>Chi)`[2],
                              il7[[7]]$`Pr(>Chi)`[2],
                              il10[[7]]$`Pr(>Chi)`[2],
                              il12[[7]]$`Pr(>Chi)`[2],
                              il13[[7]]$`Pr(>Chi)`[2],
                              il17a[[7]]$`Pr(>Chi)`[2],
                              il17b[[7]]$`Pr(>Chi)`[2],
                              il17c[[7]]$`Pr(>Chi)`[2],
                              il17d[[7]]$`Pr(>Chi)`[2],
                              il17f[[7]]$`Pr(>Chi)`[2],
                              gcsf[[7]]$`Pr(>Chi)`[2],
                              gmcsf[[7]]$`Pr(>Chi)`[2],
                              mcsf[[7]]$`Pr(>Chi)`[2],
                              tnfalpha[[7]]$`Pr(>Chi)`[2],
                              ifngamma[[7]]$`Pr(>Chi)`[2],
                              ip10[[7]]$`Pr(>Chi)`[2],
                              mcp10[[7]]$`Pr(>Chi)`[2])) %>%
  mutate(threshold_outcome=p_value_A3<0.05/21) %>%
  mutate(threshold_sex=p_value_sex<0.05/21)

######VEGF family######
 
# #VEGF-A
# vegfa<-soma %>%
#   filter(Days_symptom_update3<=time_cutoff) %>%
#   mutate(A3=factor(A3,
#                    labels=c("Control",
#                             "Case"))) %>%
#   ggplot(aes(x=Days_symptom_update3, 
#              y=	VEGFA.2597.8,
#              group=A3,
#              colour=A3)) +
#   geom_smooth() +
#   xlab("Days since symptoms onset")+
#   ylab("VEGF-A")+
#   scale_color_manual(values=c("dodgerblue4", "red3"))+
#   theme_bw()
# 
# 
# #VEGF-B
# vegfb<-soma %>%
#   filter(Days_symptom_update3<=time_cutoff) %>%
#   mutate(A3=factor(A3,
#                    labels=c("Control",
#                             "Case"))) %>%
#   ggplot(aes(x=Days_symptom_update3, 
#              y=	VEGFB.9453.12,
#              group=A3,
#              colour=A3)) +
#   geom_smooth() +
#   xlab("Days since symptoms onset")+
#   ylab("VEGF-B")+
#   scale_color_manual(values=c("dodgerblue4", "red3"))+
#   theme_bw()
# 
# 
# 
# #VEGF-C
# vegfc<-soma %>%
#   filter(Days_symptom_update3<=time_cutoff) %>%
#   mutate(A3=factor(A3,
#                    labels=c("Control",
#                             "Case"))) %>%
#   ggplot(aes(x=Days_symptom_update3, 
#              y=	VEGFC.3132.1,
#              group=A3,
#              colour=A3)) +
#   geom_smooth() +
#   xlab("Days since symptoms onset")+
#   ylab("VEGF-C")+
#   scale_color_manual(values=c("dodgerblue4", "red3"))+
#   theme_bw()
# 
# 
# #VEGF-D FIGF.13098.93
# vegfd1<-soma %>%
#   filter(Days_symptom_update3<=time_cutoff) %>%
#   mutate(A3=factor(A3,
#                    labels=c("Control",
#                             "Case"))) %>%
#   ggplot(aes(x=Days_symptom_update3, 
#              y=	FIGF.13098.93,
#              group=A3,
#              colour=A3)) +
#   geom_smooth() +
#   xlab("Days since symptoms onset")+
#   ylab("VEGF-D version 1")+
#   scale_color_manual(values=c("dodgerblue4", "red3"))+
#   theme_bw()
# 
# 
# #VEGF-D FIGF.14705.1
# vegfd2<-soma %>%
#   filter(Days_symptom_update3<=time_cutoff) %>%
#   mutate(A3=factor(A3,
#                    labels=c("Control",
#                             "Case"))) %>%
#   ggplot(aes(x=Days_symptom_update3, 
#              y=	FIGF.14705.1,
#              group=A3,
#              colour=A3)) +
#   geom_smooth() +
#   xlab("Days since symptoms onset")+
#   ylab("VEGF-D version 2")+
#   scale_color_manual(values=c("dodgerblue4", "red3"))+
#   theme_bw()
# 

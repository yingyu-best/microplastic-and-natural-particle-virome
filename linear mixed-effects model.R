###lmm
library(lme4)
library(car)
eco.sum<-read.csv("eco.function.csv")
lmm.input.diversity<-data.frame(scale(eco.sum[1:120,-c(1:4)]))

divs1<-sapply(2:ncol(lmm.input.diversity),function(j){
  message("Now j=",j," in ",ncol(lmm.input.diversity),". ",date())
  if (length(unique(lmm.input.diversity[,j]))<3){
    result<-rep(NA,18)
  } else {
    div<-data.frame(divtest=lmm.input.diversity[,j],temp=lmm.input.diversity[,1],eco.sum[1:120,1:4])
    div$type <- gsub("MP", "1", div$type)
    div$type <- gsub("NP", "0",div$type)
    fm1<-lmer(divtest~type*temp+(1|site),data=div) #river作为控制变量
    
    presult<-car::Anova(fm1,type=2)
    coefs<-coef(summary(fm1))[ , "Estimate"]  ##four coefs
    names(coefs)<-paste0(names(coefs),".mean")
    
    SEvalues<-coef(summary(fm1))[ , "Std. Error"] ##standard errors
    names(SEvalues)<-paste0(names(SEvalues),".se")
    
    tvalues<-coef(summary(fm1))[ , "t value"] ##t values
    names(tvalues)<-paste0(names(tvalues),".t")
    
    chisqP<-c(presult[,1],presult[,3])
    names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
    
    result<-c(coefs,tvalues,SEvalues,chisqP)}
  result
})
colnames(divs1)<-colnames(lmm.input.diversity[-1])

#visualization
library(dplyr)
library(stringr)
divs1.reshape=reshape2::melt(divs1)
divs1.reshape <- divs1.reshape %>%
  mutate(extracted = str_extract(divs1.reshape$Var1, "(?<=\\.).+"))

#keep mean  se
divs1.reshape <- divs1.reshape[divs1.reshape$extracted %in% c("mean", "se","P"), ]
divs1.reshape$Var1 <- gsub("\\.se|\\.mean|\\.P", "", divs1.reshape$Var1)
divs1.reshape$Var1 <- gsub("type1", "type", divs1.reshape$Var1)
library(dplyr)
library(tidyr)
# Transform the data frame
df_wide <-divs1.reshape %>%
  pivot_wider(names_from = extracted, values_from = value)
df_wide  <- df_wide [df_wide $Var1 != "(Intercept)", ]


ggplot(data.frame(df_wide),aes(y = Var1,
                               x = mean,
                               color =Var1))+
  theme_bw()+
  geom_errorbar(aes(xmin = mean - se, xmax = mean + se),position = position_dodge(0.7), width = 0.2) +#调整误差线长度
  #geom_bar(position = "stack",stat="identity",alpha = 1)+
  geom_point(size = 3) +
  facet_wrap(~Var2,scales = 'free_y',ncol = 1,strip.position = "left")+
  scale_color_manual(values=c("#EA7580FF" ,"#088BBEFF", "#F8CD9CFF"))+
  theme(axis.title = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        #axis.text.x = element_blank(), 
        axis.text = element_text(size = 10))+
  ylab("effect size")+xlab("")+
  #加显著性
  geom_text(aes(label = ifelse(P < 0.001, "***",ifelse(P < 0.01, "**",ifelse(P < 0.05, "*",ifelse(P < 0.1, "o", ""))))), 
            position = position_dodge(0.7), hjust = -0.01, color="black",size = 5)


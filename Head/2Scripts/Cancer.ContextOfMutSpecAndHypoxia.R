# our table: one line - one mutation
# (мутации мтДНК) ~ гипоксия (Buffa score AND/OR скорость деления клеток) + VAF (время появления) + ткань (тканеспецифичное что-то)

#################################################
##### 1: DERIVE TABLE
#################################################

rm(list=ls(all=TRUE))
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

Mut = read.table("../../Body/1Raw/mtDNA_snv_Oct2016.txt", head = TRUE, sep = '\t')  
Decode = read.table("../../Body/1Raw/PancancerInfoMetadata.txt", head = TRUE, sep = '\t')  
Decode = Decode[,c(2,4)]
Mut = merge(Mut,Decode, by.x = 'sample', by.y = 'submitter_donor_id', all.x = TRUE)
Hyp = read.table("../../Body/1Raw/hypoxicCancers.txt", head = TRUE, sep = '\t')  
HypMut = merge(Mut,Hyp, by = 'aliquot_id')

HypMut$Subst = paste(HypMut$ref,HypMut$var,sep = '>')
names(HypMut)
HypMut = select(HypMut,aliquot_id,sample,position,Subst,X,X.1,X.2,X.3,tissue,Tier2,tumor_var_freq,hypoxia_score_winter,hypoxia_score_ragnum,hypoxia_score_buffa)
HypMut$tumor_var_freq = as.numeric(gsub('%','',HypMut$tumor_var_freq))
nrow(HypMut) # 3110 mutations with known hypoxia
VecOfSamples = unique(HypMut$sample); length(VecOfSamples) # 828 samples

CancerTissue = c('Bladder','Bone/SoftTissue','Breast','Biliary','Cervix','Lymphoid','Myeloid','Colon/Rectum','Prostate','Esophagus','Stomach','CNS','Head/Neck','Kidney','Liver','Lung','Ovary','Pancreas','Skin','Thyroid','Uterus')  
TurnOverDays = c(200,5373,84.5,200,6,30,30,5,120,11,5.5,10000,16,1000,400,5143,11000,360,147,4138,4); length(TurnOverDays)
Turn = data.frame(CancerTissue,TurnOverDays)
Turn = Turn[order(Turn$TurnOverDays),]
HypMut = merge(HypMut,Turn,by.x = 'Tier2', by.y = 'CancerTissue')
unique(HypMut$Subst)
unique(HL$X.3)
unique(NL$X.3)
summary(HypMut$hypoxia_score_buffa)#median = 3
#################################################
####### ANALYSES:
#################################################
length(unique(HypMut$X.2))#148
length(unique(HypMut$X.3))#12
unique(HypMut$X.3)
HypMutContext = select(HypMut,Subst,Tier2,X.2,X.3,hypoxia_score_buffa)

### for light chain ###
JustLch = !stringr::str_detect(HypMutContext$X.3, 'n')
ContextLch = (JustLch = dplyr::filter(HypMutContext, JustLch))

#Which kind of contest mutspecter in normoxic cancer tissues?
NormL = (dplyr::filter(ContextLch, hypoxia_score_buffa < 3)) %>% group_by(Subst,Tier2,X.2,X.3) %>% summarise(N = n())
unique(NormL$Subst)
nomberOfmut = NormL[with(NormL, order(N)), ]
CTandTCnormLch = NormL %>% 
  filter(!stringr::str_detect(X.3,'pT>G|pC>G|pT>A|pC>A'))
unique(CTandTCnormLch$X.3)#"pC>T" "pT>C"

#for all tissues
CTandTCnormLch %>% 
  ggplot(mapping = aes(x = N, y = X.2)) +
  geom_col()

#for any tissue with context.mut N>5
CTandTCnormLch %>% 
  dplyr::filter(N >5) %>%
  ggplot(mapping = aes(x = X.2, y = N))+
  geom_col() +
  ggtitle('Контекст основных мутационных переходов в раковых тканях',
          subtitle = 'Нормоксия') +
  coord_flip() +
  facet_wrap(~ Tier2)

#Which kind of contecst mutspectr in hypoxic cancer tissues?
HypoxicL = (dplyr::filter(ContextLch, hypoxia_score_buffa > 3)) %>% group_by(Tier2,X.2,X.3) %>% summarise(N = n())
CTandTChypoxicLch = HypoxicL %>% 
  filter(!stringr::str_detect(X.3,'pT>G|pC>G|pT>A|pC>A'))
unique(CTandTChypoxicLch$X.3)

#for all tissues
CTandTChypoxicLch %>% 
  ggplot(mapping = aes(x = N, y = X.2)) +
  geom_col()
#for any tissue with c.mut N>5
CTandTChypoxicLch %>% 
  dplyr::filter(N >5) %>%
  ggplot(mapping = aes(x = X.2, y = N))+
  geom_col() +
  ggtitle('Контекст основных мутационных переходов в раковых тканях',
          subtitle = 'Гипоксия') +
  coord_flip() +
  facet_wrap(~ Tier2)

### for heavy chain ###
JustHch = !stringr::str_detect(HypMutContext$X.3, 'p')
ContextHch = (JustHch = dplyr::filter(HypMutContext, JustHch))

#Which kind of contecst mutspectr in normoxic cancer tissues Hch?
NormH = (dplyr::filter(ContextHch, hypoxia_score_buffa < 3)) %>% group_by(Subst,Tier2,X.2,X.3) %>% summarise(N = n())
unique(NormH$Subst)
unique(NormH$X.3)
CTandTCnormHch = NormH %>% 
  filter(!stringr::str_detect(X.3,'nT>G|nC>G|nT>A|nC>A'))
unique(CTandTCnormHch$X.3)
#for all tissues
CTandTCnormHch %>% 
  ggplot(mapping = aes(x = N, y = X.2)) +
  geom_col()
#for any tissue with c.mut N>5
CTandTCnormHch %>% 
  dplyr::filter(N >5) %>%
  ggplot(mapping = aes(x = X.2, y = N))+
  geom_col() +
  ggtitle('Контекст основных мутационных переходов в раковых тканях',
          subtitle = 'Нормоксия') +
  coord_flip() +
  facet_wrap(~ Tier2)

#in hypoxic cancer tissues?
HypoxicH = (dplyr::filter(ContextHch, hypoxia_score_buffa > 3)) %>% group_by(Tier2,X.2,X.3) %>% summarise(N = n())
CTandTChypoxicHch = HypoxicH %>% 
  filter(!stringr::str_detect(X.3,'nT>G|nC>G|nT>A|nC>A'))
unique(CTandTChypoxicHch$X.3)

#for all tissues
CTandTChypoxicHch %>% 
  ggplot(mapping = aes(x = N, y = X.2)) +
  geom_col()
#for any tissue with c.mut N>5
CTandTChypoxicHch %>% 
  dplyr::filter(N >5) %>%
  ggplot(mapping = aes(x = X.2, y = N))+
  geom_col() +
  ggtitle('Контекст основных мутационных переходов в раковых тканях',
          subtitle = 'Гипоксия') +
  coord_flip() +
  facet_wrap(~ Tier2)

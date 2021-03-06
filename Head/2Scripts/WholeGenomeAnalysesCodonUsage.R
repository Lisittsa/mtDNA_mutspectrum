rm(list=ls(all=TRUE))

############ read aminoacid and codon sequences
WholeGenomes = read.table("../../Body/2Derived/AllGenesCodons.csv", header = TRUE, sep = '\t')

Final = c()
for (step in 1:nrow(WholeGenomes))
{ # step = 3
  if (WholeGenomes$Quality[step] == 1)
  {
  species = as.character(WholeGenomes$Species[step])
  gene = as.character(WholeGenomes$Gene[step])
  ############ translate codons to aminoacids
  ### make triplets 
  Codons = as.character(WholeGenomes$Codons[step])
  AminoSeq = as.character(WholeGenomes$ModifiedAmino[step])
  AminoSeqVec = unlist(strsplit(AminoSeq,''))
  length(Codons)
  length(AminoSeq)
  CodonsVec = unlist(strsplit(Codons,''))
  StartNuc = 1
  CodonsVec1 = c()
  # if length(CodonsVec)/3 == integer
  for (i in 1:(length(CodonsVec)/3))
  {
    CodonsVec1 = c(CodonsVec1,paste(CodonsVec[StartNuc : (StartNuc+2)],collapse = ''))
    StartNuc = StartNuc+3
  }
  
  TranslationalTest = data.frame(CodonsVec1,AminoSeqVec)
  
  ### translate them
  Temp = as.character(TranslationalTest$CodonsVec1)
  Temp = gsub("TTT|TTC",'F',Temp)
  Temp = gsub("TTA|TTG|CTT|CTC|CTA|CTG",'L',Temp)
  Temp = gsub("ATT|ATC",'I',Temp)
  Temp = gsub("ATG|ATA",'M',Temp)
  Temp = gsub("GTT|GTC|GTA|GTG",'V',Temp)
  Temp = gsub("TCT|TCC|TCA|TCG|AGT|AGC",'S',Temp)
  Temp = gsub("CCT|CCC|CCA|CCG",'P',Temp)
  Temp = gsub("ACT|ACC|ACA|ACG",'T',Temp)
  Temp = gsub("GCT|GCC|GCA|GCG",'A',Temp)
  Temp = gsub("TAT|TAC",'Y',Temp)
  Temp = gsub("AGA|AGG|TAA|TAG",'*',Temp)
  Temp = gsub("CAT|CAC",'H',Temp)
  Temp = gsub("CAA|CAG",'Q',Temp)
  Temp = gsub("AAT|AAC",'N',Temp)
  Temp = gsub("AAA|AAG",'K',Temp)
  Temp = gsub("GAT|GAC",'D',Temp)
  Temp = gsub("GAA|GAG",'E',Temp)
  Temp = gsub("TGT|TGC",'C',Temp)
  Temp = gsub("TGG|TGA",'W',Temp)
  Temp = gsub("CGT|CGC|CGA|CGG",'R',Temp)
  Temp = gsub("GGT|GGC|GGA|GGG",'G',Temp)
  TranslationalTest$TranslatedAa = Temp
  
  nrow(TranslationalTest[TranslationalTest$TranslatedAa != TranslationalTest$AminoSeqVec,]) # should be zero
  
  VecOfAllCodons = c('TTT','TTC','TTA','TTG','CTT','CTC','CTA','CTG','ATT','ATC','ATG','ATA','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','AGT','AGC','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','TAT','TAC','AGA','AGG','TAA','TAG','CAT','CAC','CAA','CAG','AAT','AAC','AAA','AAG','GAT','GAC','GAA','GAG','TGT','TGC','TGG','TGA','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG')
  CodonUsage = c(VecOfAllCodons,as.character(TranslationalTest$CodonsVec1))
  CU = data.frame(table(CodonUsage)); CU$Freq = CU$Freq-1
  NamesOfCodons = as.character(CU$CodonUsage)
  
  Line = c(species, gene, CU$Freq)
  Final = rbind(Final,Line)
  }
}

Names = c('Species','Gene',NamesOfCodons)
Final = as.data.frame(Final)
names(Final) = Names

str(Final$AAA)
Final <- data.frame(lapply(Final, as.character), stringsAsFactors=FALSE)
str(Final)

## count the number of A T G C in neutral positions of each gene (8 synon fourlfold codons)
Final$NeutralA = as.numeric(as.character(Final$CTA[1])) + as.numeric(Final$GTA) + as.numeric(Final$TCA) + as.numeric(Final$CCA)  + as.numeric(Final$ACA)  + as.numeric(Final$GCA)  + as.numeric(Final$CGA)  + as.numeric(Final$GGA)
Final$NeutralT = as.numeric(Final$CTT) + as.numeric(Final$GTT) + as.numeric(Final$TCT) + as.numeric(Final$CCT)  + as.numeric(Final$ACT)  + as.numeric(Final$GCT)  + as.numeric(Final$CGT)  + as.numeric(Final$GGT)
Final$NeutralG = as.numeric(Final$CTG) + as.numeric(Final$GTG) + as.numeric(Final$TCG) + as.numeric(Final$CCG)  + as.numeric(Final$ACG)  + as.numeric(Final$GCG)  + as.numeric(Final$CGG)  + as.numeric(Final$GGG)
Final$NeutralC = as.numeric(Final$CTC) + as.numeric(Final$GTC) + as.numeric(Final$TCC) + as.numeric(Final$CCC)  + as.numeric(Final$ACC)  + as.numeric(Final$GCC)  + as.numeric(Final$CGC)  + as.numeric(Final$GGC)

## merge with taxonomy

Taxa = read.table("../../Body/2Derived/full_table.txt", header = TRUE, sep = '\t')
Taxa = Taxa[,c(1,2)]
names(Taxa)=c('Species','Taxonomy')
nrow(Final)
Final = merge(Final,Taxa)
nrow(Final)  # somebody is missing!!! check it!!

## derive class from taxonomy
Final$Class = 'SOLYANKA';
for (i in 1:nrow(Final))
{ 
  if (length(grep('Mammalia', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'Mammalia'; }
  if (length(grep('Amphibia', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'Amphibia'; }
  
  if (length(grep('Actinopterygii', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'Actinopterygii'; }
  if (length(grep('Aves', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'Aves'; }
  
  if (length(grep('Lepidosauria', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'Reptilia'; }
  if (length(grep('Testudines', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'Reptilia'; }
  if (length(grep('Crocodylia', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'Reptilia'; }
  
  if (length(grep('Dipnoi', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'AncientFish'; }
  if (length(grep('Coelacanthiformes', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'AncientFish'; }
  if (length(grep('Cyclostomata', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'AncientFish'; }
  if (length(grep('Chondrichthyes', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'AncientFish'; }
  if (length(grep('Cephalochordata', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'AncientFish'; }
  if (length(grep('Tunicata', Final$Taxonomy[i])) == 1) {Final$Class[i] = 'AncientFish'; }
}

table(Final$Class)

## print it out:

write.table(Final, "../../Body/3Results/AllGenesCodonUsage.txt", sep = '\t')  

#### how many species has all 13 genes with high quality?
test = data.frame(table(Final$Species))
nrow(test[test$Freq == 13,]) # 3834
nrow(test[test$Freq != 13,]) # 43
nrow(test[test$Freq < 13,]) # 22
nrow(test[test$Freq > 13,]) # 21

#### who are they?
less = test[test$Freq < 13,]$Var1 
# Aplidium_conicum          Asymmetron_inferum        Asymmetron_lucayanum      Botrylloides_leachii      Botrylloides_pizoni       Branchiostoma_belcheri    Branchiostoma_floridae    Branchiostoma_japonicum  
#[9] Branchiostoma_lanceolatum Chionodraco_myersi        Ciona_intestinalis        Datnioides_microlepis     Epigonichthys_cultellus   Epigonichthys_maldivensis Herdmania_momus           Myrichthys_maculosus     
#[17] Pelomedusa_subrufa        Podiceps_cristatus        Prioniturus_luconensis    Psittacus_erithacus       Sphenodon_punctatus       Styela_plicata           

table(Final[Final$Species %in% less,]$Class)
#Actinopterygii    AncientFish           Aves       Reptilia 
# 36             48             36             23 

more = test[test$Freq > 13,]$Var1 
#Aceros_waldeni             Ardeola_bacchus            Botaurus_stellaris         Butorides_striata          Champsocephalus_gunnari    Egretta_sacra              Gorsachius_magnificus     
#[8] Heteronotia_binoei         Hoplobatrachus_rugulosus   Ixobrychus_eurhythmus      Ixobrychus_sinensis        Nannopterum_brasilianus    Notiomystis_cincta         Penelopides_panini        
#[15] Phalacrocorax_carbo        Phoebastria_albatrus       Phoebastria_immutabilis    Phoebastria_nigripes       Thalassarche_melanophrys   Tropiocolotes_tripolitanus Turdus_philomelos 

table(Final[Final$Species %in% more,]$Class)  
#Actinopterygii       Amphibia           Aves       Reptilia 
#16             14            262             33 


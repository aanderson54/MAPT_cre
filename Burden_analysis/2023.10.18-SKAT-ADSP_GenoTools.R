args = commandArgs(trailingOnly=TRUE)

library(SKAT)

# Expects argument 1 to be a fam file, argument 2 to be a covariate file in plink format, and argument 3 to be a burden set
# The burden set as well as the phenotype designation should be binary and 0/1 coded
FAM_Cov<-Read_Plink_FAM_Cov(args[1],args[2],Is.binary=TRUE,flag1=1)
BurdenSet <- as.matrix(read.table(args[3], header=FALSE, sep ="\t"))

# Extract covariates and other information from FAM_Cov
PC1=FAM_Cov$PC1
PC2=FAM_Cov$PC2
PC3=FAM_Cov$PC3
PC4=FAM_Cov$PC4
PC5=FAM_Cov$PC5
PC6=FAM_Cov$PC6
PC7=FAM_Cov$PC7
PC8=FAM_Cov$PC8
PC9=FAM_Cov$PC9
PC10=FAM_Cov$PC10
SAS=FAM_Cov$SAS
MDE=FAM_Cov$MDE
EUR=FAM_Cov$EUR
EAS=FAM_Cov$EAS
CAS=FAM_Cov$CAS
AMR=FAM_Cov$AMR
AJ=FAM_Cov$AJ
AFR=FAM_Cov$AFR
AAC=FAM_Cov$AAC
Sex=FAM_Cov$Sex
Age_harmonized_wImpute=FAM_Cov$Age_harmonized_wImpute
IsBloodwImpute=FAM_Cov$IsBloodwImpute
Illumina_NovaSeq=FAM_Cov$Illumina_NovaSeq
Illumina_HiSeqX=FAM_Cov$Illumina_HiSeqX
Illumina_HiSeq_2000=FAM_Cov$Illumina_HiSeq_2000
USUHS=FAM_Cov$USUHS
WashU=FAM_Cov$WashU
USUHS_Miami=FAM_Cov$USUHS_Miami
Broad=FAM_Cov$Broad
Illumina=FAM_Cov$Illumina
Baylor=FAM_Cov$Baylor
GENENTECH=FAM_Cov$GENENTECH
MEDGENOME=FAM_Cov$MEDGENOME
NYGC=FAM_Cov$NYGC
MACROGEN=FAM_Cov$MACROGEN
H2=FAM_Cov$H2
Pheno<-FAM_Cov$Phenotype

# Create null model for SKAT
obj<-SKAT_Null_Model(Pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SAS + MDE + EUR + EAS + CAS + AMR + AJ + AFR + AAC + Sex + Age_harmonized_wImpute + IsBloodwImpute + Illumina_NovaSeq + Illumina_HiSeqX + Illumina_HiSeq_2000 + USUHS + Broad + Illumina + Baylor + GENENTECH + MEDGENOME + NYGC + MACROGEN + H2,out_type="D",Adjustment=FALSE)

# Perform SKAT with robust method
SKATout<-SKATBinary_Robust(BurdenSet,obj,method="Burden")
SKATout$p.value
SKATout$mac
SKATout$param$n.marker


# Create a matrix to store results
BigTable <- matrix(nrow = 1, ncol = 6)
BigTable[1,1:3] <- args[1:3]
BigTable[1,4:6] <- c(SKATout$p.value, SKATout$mac, SKATout$param$n.marker)

# SKAT part is done, now get effect size and run Fisher's Exact Test
TotalControls <- as.data.frame(table(Pheno))[1,2]
TotalCases <- as.data.frame(table(Pheno))[2,2]
PhenoControls <- (Pheno-1)*-1
CasesWithVar <- as.data.frame(table(Pheno*BurdenSet))[2,2]
CasesWithVar[is.na(CasesWithVar)] <- 0
CasesWithoutVar <- TotalCases-CasesWithVar
ControlsWithVar <- as.data.frame(table(PhenoControls*BurdenSet))[2,2]
ControlsWithVar[is.na(ControlsWithVar)] <- 0
ControlsWihtoutVar <- TotalControls-ControlsWithVar

# Function to run Fisher's Exact Test on a given row of data
row_fisher <- function(row, alt = 'two.sided', cnf = 0.95) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           p_val = f$p.value,
           or = f$estimate[[1]],
           or_ll = f$conf.int[1],
           or_ul = f$conf.int[2]))
}

# Prepare data for Fisher's Exact Test
TestData <- matrix(nrow = 1, ncol = 4)
TestData[1,1:4] = c(CasesWithVar, CasesWithoutVar, ControlsWithVar, ControlsWithoutVar)
colnames(TestData) <- c('CaW', 'CaWo', 'CoW', 'CoWo')

# Apply Fisher's test to data
FET <- t(apply(TestData, 1, row_fisher))

# Combine all results
OutputTable <- cbind(BigTable,FET)
OutputName <- paste(args[1],args[2],args[3],"_SKAT-output.txt",sep = "")

write.table(OutputTable,file=OutputName,quote=FALSE,sep="\t",row.names=FALSE,col.names = FALSE,na="")
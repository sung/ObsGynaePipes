#Used librarys
library("gplots") #function rich.colors
library("preprocessCore") #function normalize.quantiles
	
#Color used for graphes
my.col <- colorRampPalette(c("#FFFFFF", "black", "blue", "#FA8072","#00A2FF", "#00CC00", "#E0E0E0"))(7) #1:Backgroundcolor for all graphs, 2: Foregroundcolor for all graphs (E6E6E6), 3: Fill for histograms, 4: Red, for boxplots, 5: Blue, for boxplots, 6: Green, for boxplots, 7: Light gray

############
#Parameters#
############
###Set folders
	setwd("~/TauComparison")
	folder <- c("~/TauComparison/")

###Choose organism and data set for analysis
organism <- "Hum" #"Hum" or "Mus"
expDataSource <- "Fagerberg" #Brawand, ENCODE, Thorrez for Mouse; Brawand, Fagerberg, Ge for Human;
	add <- "" 

###Human
#####Fagerberg - RNA-seq, downloaded ArrayExpress, 27 tissues
#####Brawand - RNA-seq, Bgee processed, 8 tissues
#####Ge - Microarray, Bgee processed, 32 tissues

###Mouse
#####ENCODE - RNA-seq, self processed, 22 tissues
#####Brawand - RNA-seq, Bgee processed, 6 tissues
#####Thorrez - Microarray, Bgee processed, 19 tissues (placent missing)


###+++###
fTissueNames <- function(organism, dataSource)
{	
	if(organism == "Mus")
	{
		if (dataSource == "ENCODE") {
			tissuesNames <- c("cerebellum", "cortex", "heart", "kidney", "liver", "lung", "placenta", "smintestine", "spleen", "testis", "thymus", "adrenal", "bladder", "colon", "duodenum", "flobe", "gfat", "lgintestine", "mamgland", "ovary", "sfat", "stomach")
		} else if (dataSource == "Brawand") {
			tissuesNames <- c("brain","cerebellum","heart","kidney","liver", "testis")
		} else if (dataSource == "Thorrez") {
			tissuesNames <- c("diaphragm", "spleen", "muscle", "liver", "brain", "lung", "kidney", "adrenal", "marrow", "adipose", "pituitary", "sgland", "svesicle", "thymus", "testis", "heart", "smintestine", "eye", "fgonad")
		}
	}  else if (organism == "Hum") {
		if (dataSource == "Fagerberg") {
			tissuesNames <- c("colon","kidney", "liver", "pancreas", "lung", "prostate", "brain", "stomach", "spleen", "lymphnode", "appendix", "smint", "adrenal", "duodenum", "fat", "endometrium", "placenta", "testis", "gbladder", "ubladder", "thyroid", "esophagus", "heart", "skin", "ovary", "bonem", "sgland")	
		} else if (dataSource == "Brawand") {
			tissuesNames <- c("fcortex","pcortex","tlobe","cerebellum","heart","kidney","liver", "testis")	
		} else if (dataSource == "Ge") {
			tissuesNames <- c("heart", "thymus", "spleen", "fgonad", "kidney", "muscle", "pancreas", "prostate", "smintestine", "colon", "placenta", "ubladder", "mamgland", "uterus", "thyroid", "skin", "trachea", "cerebellum", "brain", "adrenal", "marrow", "amygdala", "nucleus", "callosum", "Ammons", "thalamus", "pituitary", "spinal", "testis", "liver", "stomach", "lung")
		} 	
	}			
	return(tissuesNames)
}	
###***###***###

###+++###
fTissuePrintNames <- function(organism, dataSource)
{	
	if(organism == "Mus")
	{
		if (dataSource == "ENCODE") {
			tissuesPrintNames <- c("Cerebellum", "Cortex", "Heart", "Kidney", "Liver", "Lung","Placenta","Small Intestine","Spleen","Testis", "Thymus", "Adrenal", "Bladder", "Colon", "Duodenum", "Frontal Lobe", "Genital Fat Pad", "Large Intestine", "Mammary Gland", "Ovary", "Subcutaneous Fat Pad", "Stomach") 
		} else if (dataSource == "Brawand") {
			tissuesPrintNames <- c("Brain","Cerebellum","Heart","Kidney","Liver", "Testis")
		} else if (dataSource == "Thorrez") {
			tissuesPrintNames <- c("Diaphragm", "Spleen", "Muscle", "Liver", "Brain", "Lung", "Kidney", "Adrenal", "Bone Marrow", "Adipose Tissue", "Pituitary Gland", "Salivary Gland", "Seminal Vesicle", "Thymus", "Testis", "Heart", "Small Intestine", "Eye", "Female Gonad")
		} 			
	}  else if (organism == "Hum") {
		if (dataSource == "Fagerberg") {
			tissuesPrintNames <- c("Colon","Kidney", "Liver", "Pancreas", "Lung", "Prostate", "Brain", "Stomach", "Spleen", "Lymph Node", "Appendix", "Small Intestine", "Adrenal", "Duodenum", "Fat", "Endometrium", "Placenta", "Testis", "Gall Bladder", "Urinary Bladder", "Thyroid", "Esophagus","Heart", "Skin", "Ovary", "Bone Marrow", "Salivary Gland")	
		} else if (dataSource == "Brawand") {
			tissuesPrintNames <- c("Frontal Cortex","Prefrontal Cortex","Temporal Lobe","Cerebellum","Heart","Kidney","Liver", "Testis")
		} else if (dataSource == "Ge") {
			tissuesPrintNames <- c("Heart", "Thymus", "Spleen", "Female Gonad", "Kidney", "Skeletal Muscle", "Pancreas", "Prostate Gland", "Small Intestine", "Colon", "Placenta", "Urinary Bladder", "Mammary Gland", "Uterus", "Thyroid Gland", "Skin", "Trachea", "Cerebellum", "Brain", "Adrenal Gland", "Bone Marrow", "Amygdala", "Caudate Nucleus", "Corpus Callosum", "Ammons Horn", "Thalamus", "Pituitary Gland", "Spinal Cord", "Testis", "Liver", "Stomach", "Lung"  )
		}	
	}	
		
	return(tissuesPrintNames)
}	
###***###***###	
	
	#Tissue names for different organisms and data sets
		tissuesEFNames <- c("brain", "heart", "kidney", "liver", "lung", "placenta", "smint", "spleen", "testis", "adrenal", "bladder", "colon", "duodenum", "ovary", "fat", "stomach") #common tissues between ENCODE and Fagerberg 16
		tissuesBrNames <- c("brain", "cerebellum", "heart", "kidney", "liver", "testis") #ENCODE, Brawand human and mouse
		tissuesMUNames <- c("brain", "heart", "kidney", "liver", "testis") #
		tissuesMUNamesCortex <- c("cortex", "heart", "kidney", "liver", "testis") #
		
		tissuesBrHNames <- c("fcortex", "cerebellum", "heart","kidney","liver", "testis")
		tissuesEFNamesMouse <- c("cortex", "heart", "kidney", "liver", 	"lung", "placenta", "smintestine", "spleen", "testis", "adrenal", 	"bladder", "colon", "duodenum", "ovary", "sfat", "stomach") #common tissues between ENCODE and Fagerberg 16 in mouse
		tissuesEFNamesHuman <- c("brain", "heart", "kidney", "liver", 	"lung", "placenta", "smint", "spleen", "testis", "adrenal", 	"ubladder", "colon", "duodenum", "ovary", "fat", "stomach") #common tissues between ENCODE and Fagerberg 16 in human	
		tissuesGTNames <- c("spleen", "muscle", "liver", "brain", "lung", "kidney", "adrenal", "marrow", "pituitary", "thymus", "testis", "heart", "smintestine", "fgonad")
		

############
#Input data#
############
fInputData <- function() #
{
		##Tissue expression
		if (organism == "Mus")
		{
			if (expDataSource == "ENCODE") {
				orgExpression <- read.table(paste(folder,"EncodeCshlAdult8wksEnsV68RNAseqGene.txt",sep=""), sep="\t", header=TRUE)
			} else if (expDataSource == "Brawand") {
				orgExpression <- read.table(paste(folder,"Mus_musculus_RNA-Seq_RPKM_GSE30352_Tissues.txt",sep=""), sep="\t", header=TRUE)	
			}  else if (expDataSource == "Thorrez") {
				orgExpression <- read.table(paste(folder,"Mus_musculus_probesets_GSE9954_A-AFFY-45_gcRMA_TissuesBgee.txt",sep=""), sep="\t", header=TRUE)
			}		
		}  else if (organism == "Hum") {
			if (expDataSource == "Fagerberg") {
				orgExpression <- read.table(paste(folder, "ArrayExpressHumanAdultEnsV69RNAseq.txt",sep=""), sep="\t", header=TRUE) 
				colnames(orgExpression) <-lapply(colnames(orgExpression),function(x){x <- unlist(strsplit(toString(x), split='_', fixed=TRUE))[1]})
			} else if (expDataSource == "Brawand") {
				orgExpression <- read.table(paste(folder,"Homo_sapiens_RNA-Seq_RPKM_GSE30352_Tissues.txt",sep=""), sep="\t", header=TRUE)
			}  else if (expDataSource == "Ge") {
				orgExpression <- read.table(paste(folder,"Homo_sapiens_probesets_GSE2361_A-AFFY-33_gcRMA_TissuesBgee.txt",sep=""), sep="\t", header=TRUE)
			}
		}	
return(orgExpression)
}
#######

###############
###Functions###				
###############

###+++###
#Function requires data frame with expression values
#Mean values between replicares are calculated
fReplicateMean <- function(x, source, organism, names)
{
	if(source == "Brawand")
	{
		if(organism == "Hum")
		{
			x$Averaged.RPKM.fcortex <- rowMeans(x[,regexpr("frontal.cortex",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.pcortex <- rowMeans(x[,regexpr("prefront.cortex",colnames(x))>0])
			x$Averaged.RPKM.tlobe <- x[,regexpr("temporal.lobe",colnames(x))>0]	
			x$Averaged.RPKM.cerebellum <- rowMeans(x[,regexpr("cerebellum",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.heart <- rowMeans(x[,regexpr("heart",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.kidney <- rowMeans(x[,regexpr("kidney",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.liver <- rowMeans(x[,regexpr("liver",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.testis <- rowMeans(x[,regexpr("testis",colnames(x))>0], na.rm=TRUE, dim=1)
			x <- x[,c("Ensembl.Gene.ID", names)]
		} else if(organism == "Mus") {
			x$Averaged.RPKM.brain <- rowMeans(x[,regexpr("brain",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.cerebellum <- rowMeans(x[,regexpr("cerebellum",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.heart <- rowMeans(x[,regexpr("heart",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.kidney <- rowMeans(x[,regexpr("kidney",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.liver <- rowMeans(x[,regexpr("liver",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.testis <- rowMeans(x[,regexpr("testis",colnames(x))>0], na.rm=TRUE, dim=1)
			x <- x[,c("Ensembl.Gene.ID", names)]
		}	
	} else if(source == "ENCODE"){
		if(organism == "Mus")
		{
			x$Averaged.RPKM.cerebellum <- rowMeans(x[,regexpr("Cbellum",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.cortex <- rowMeans(x[,regexpr("Cortex",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.heart <- rowMeans(x[,regexpr("Heart",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.kidney <- rowMeans(x[,regexpr("Kidney",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.liver <- rowMeans(x[,regexpr("Liver",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.lung <- rowMeans(x[,regexpr("Lung",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.placenta <- rowMeans(x[,regexpr("Plac",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.smintestine <- rowMeans(x[,regexpr("Smint",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.spleen <- rowMeans(x[,regexpr("Spleen",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.testis <- rowMeans(x[,regexpr("Testis",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.thymus <- rowMeans(x[,regexpr("Thymus",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.adrenal <- rowMeans(x[,regexpr("Adrenal",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.bladder <- rowMeans(x[,regexpr("Bladder",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.colon <- rowMeans(x[,regexpr("Colon",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.duodenum <- rowMeans(x[,regexpr("Duod",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.flobe <- rowMeans(x[,regexpr("Flobe",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.gfat <- rowMeans(x[,regexpr("Gfat",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.lgintestine <- rowMeans(x[,regexpr("Lgint",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.mamgland <- rowMeans(x[,regexpr("Mamg",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.ovary <- rowMeans(x[,regexpr("Ovary",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.sfat <- rowMeans(x[,regexpr("Sfat",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.stomach <- rowMeans(x[,regexpr("Stom",colnames(x))>0], na.rm=TRUE, dim=1)
			x <- x[,c("Ensembl.Gene.ID", names)]
		} 
	} else if(source == "Fagerberg"){
		x$Averaged.RPKM.colon <- rowMeans(x[,regexpr("colon",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.kidney <- rowMeans(x[,regexpr("kidney",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.liver <- rowMeans(x[,regexpr("liver",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.pancreas <- rowMeans(x[,regexpr("pancreas",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.lung <- rowMeans(x[,regexpr("lung",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.prostate <- rowMeans(x[,regexpr("prostate",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.brain <- rowMeans(x[,regexpr("brain",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.stomach <- rowMeans(x[,regexpr("stomach",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.spleen <- rowMeans(x[,regexpr("spleen",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.lymphnode <- rowMeans(x[,regexpr("lymphnode",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.appendix <- rowMeans(x[,regexpr("appendix",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.smint <- rowMeans(x[,regexpr("smallintestine",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.adrenal <- rowMeans(x[,regexpr("adrenal",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.duodenum <- rowMeans(x[,regexpr("duodenum",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.fat <- rowMeans(x[,regexpr("fat",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.endometrium <- rowMeans(x[,regexpr("endometrium",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.placenta <- rowMeans(x[,regexpr("placenta",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.testis <- rowMeans(x[,regexpr("testis",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.gbladder <- rowMeans(x[,regexpr("gallbladder",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.ubladder <- rowMeans(x[,regexpr("urinarybladde",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.thyroid <- rowMeans(x[,regexpr("thyroid",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.esophagus <- rowMeans(x[,regexpr("esophagus",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.heart <- rowMeans(x[,regexpr("heart",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.skin <- rowMeans(x[,regexpr("skin",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.ovary <- rowMeans(x[,regexpr("ovary",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.bonem <- rowMeans(x[,regexpr("bonem",colnames(x))>0], na.rm=TRUE, dim=1)
		x$Averaged.RPKM.sgland <- rowMeans(x[,regexpr("salivarygland",colnames(x))>0], na.rm=TRUE, dim=1)
		x <- x[,c("Ensembl.Gene.ID", names)]
	}  else if(source == "Thorrez"){
		if(organism == "Mus") {
		 	x$Averaged.RPKM.diaphragm <- rowMeans(x[,regexpr("diaphragm",colnames(x))>0], na.rm=TRUE, dim=1)
		 	x$Averaged.RPKM.spleen <- rowMeans(x[,regexpr("spleen",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.muscle <- rowMeans(x[,regexpr("muscle",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.liver <- rowMeans(x[,regexpr("liver",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.brain <- rowMeans(x[,regexpr("brain",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.lung <- rowMeans(x[,regexpr("lung",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.kidney <- rowMeans(x[,regexpr("kidney",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.adrenal <- rowMeans(x[,regexpr("adrenal",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.marrow <- rowMeans(x[,regexpr("marrow",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.adipose <- rowMeans(x[,regexpr("adipose",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.pituitary <- rowMeans(x[,regexpr("pituitary",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.sgland <- rowMeans(x[,regexpr("saliva",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.svesicle <- rowMeans(x[,regexpr("vesicle",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.thymus <- rowMeans(x[,regexpr("thymus",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.testis <- rowMeans(x[,regexpr("testis",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.heart <- rowMeans(x[,regexpr("heart",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.smintestine <- rowMeans(x[,regexpr("intestine",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.eye <- rowMeans(x[,regexpr("eye",colnames(x))>0], na.rm=TRUE, dim=1)
			x$Averaged.RPKM.fgonad <- rowMeans(x[,regexpr("gonad",colnames(x))>0], na.rm=TRUE, dim=1)
			x <- x[,c("Ensembl.Gene.ID", names)]
		 	}
	} else if(source == "Ge"){
		if(organism == "Hum") {
		 	x$Averaged.RPKM.heart <- x[,regexpr("heart",colnames(x))>0]
		 	x$Averaged.RPKM.thymus <-  x[,regexpr("thymus",colnames(x))>0]
			x$Averaged.RPKM.spleen <-  x[,regexpr("spleen",colnames(x))>0]
			x$Averaged.RPKM.fgonad <-  x[,regexpr("female",colnames(x))>0]
			x$Averaged.RPKM.kidney <-  x[,regexpr("kidney",colnames(x))>0]
			x$Averaged.RPKM.muscle <- x[,regexpr("muscle",colnames(x))>0]
		 	x$Averaged.RPKM.pancreas <-  x[,regexpr("pancreas",colnames(x))>0]
			x$Averaged.RPKM.prostate <-  x[,regexpr("prostate",colnames(x))>0]
			x$Averaged.RPKM.smintestine <-  x[,regexpr("intestine",colnames(x))>0]
			x$Averaged.RPKM.colon <-  x[,regexpr("colon",colnames(x))>0]
			x$Averaged.RPKM.placenta <-  x[,regexpr("placenta",colnames(x))>0]
			x$Averaged.RPKM.ubladder <- x[,regexpr("bladder",colnames(x))>0]
		 	x$Averaged.RPKM.mamgland <-  x[,regexpr("mammary",colnames(x))>0]
			x$Averaged.RPKM.uterus <-  x[,regexpr("uterus",colnames(x))>0]
			x$Averaged.RPKM.thyroid <-  x[,regexpr("thyroid",colnames(x))>0]
			x$Averaged.RPKM.skin <-  x[,regexpr("skin",colnames(x))>0]
			x$Averaged.RPKM.trachea <- x[,regexpr("trachea",colnames(x))>0]
		 	x$Averaged.RPKM.cerebellum <-  x[,regexpr("cerebellum",colnames(x))>0]
			x$Averaged.RPKM.brain <-  x[,regexpr("brain",colnames(x))>0]
			x$Averaged.RPKM.adrenal <-  x[,regexpr("adrenal",colnames(x))>0]
			x$Averaged.RPKM.marrow <-  x[,regexpr("marrow",colnames(x))>0]
			x$Averaged.RPKM.amygdala <- x[,regexpr("amygdala",colnames(x))>0]
		 	x$Averaged.RPKM.nucleus <-  x[,regexpr("nucleus",colnames(x))>0]
			x$Averaged.RPKM.callosum <-  x[,regexpr("callosum",colnames(x))>0]
			x$Averaged.RPKM.Ammons <-  x[,regexpr("horn",colnames(x))>0]
			x$Averaged.RPKM.thalamus <-  x[,regexpr("thalamus",colnames(x))>0]
			x$Averaged.RPKM.pituitary <- x[,regexpr("pituitary",colnames(x))>0]
		 	x$Averaged.RPKM.spinal <-  x[,regexpr("spinal",colnames(x))>0]
			x$Averaged.RPKM.testis <-  x[,regexpr("testis",colnames(x))>0]
			x$Averaged.RPKM.liver <-  x[,regexpr("liver",colnames(x))>0]
			x$Averaged.RPKM.stomach <-  x[,regexpr("stomach",colnames(x))>0]
			x$Averaged.RPKM.lung <-  x[,regexpr("lung",colnames(x))>0]
			x <- x[,c("Ensembl.Gene.ID", names)]
	 		}
	} 	
	return(x)
}
###***###***###

###+++###
fPlotExpression <- function(x, dataName, fileName, names) #names=tissuesPrintNames #dataName=string to be printed on the x-axis
{
	dev.new(height=9, width=12)
	par(cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
	palette(rev(rich.colors(ncol(x))))
	
	plot(density(x[,1],n=1000), main = "Expression values among different tissues", xlab=dataName,col=(1), lwd=3)
	for(i in c(2:length(names)))
	{	
		lines(density(x[,i],n = 1000), col=(i), lwd=3)
	}
	legend("topright", names, col=(1:length(names)), lty="solid", lwd=3)
	
	dev.copy2pdf(device=quartz, file=paste(folder, organism, "Expression", expDataSource, "", fileName, add, ".pdf", sep=""),onefile=TRUE)#,paper="A4r"
	#dev.off()

return()
}	
###***###***###

###+++###
#Function requires data frame to be normalized
#1. All 0 are set to NA, to exclude them from quatile normalization
#2. Data are quantile normalized
#3. 0 values (the one set to NA) are set back to 0
fQN <- function(x) #
{
	x[x==0] <- NA
	x_m <- as.matrix(x)
	x <- normalize.quantiles(x_m)
	x[is.na(x)] <- 0
	return(data.frame(x))
}	
###***###***###

###+++###	
#Function require a vector with expression of one gene in different tissues.
#Mean is calculated taking in account tissues with 0 expression. 2+0+4=2
fmean <- function(x)
	{
		if(!all(is.na(x)))
	 	{
	 		res <- mean(x, na.rm=TRUE)
	 	} else {
	 		res <- NA
	 	}
	 	return(res)
	}
###***###***###	

###+++###	
#Function require a vector with expression of one gene in different tissues.
#Max is calculated taking in account tissues with 0 expression. 2+0+4=2
fmax <- function(x)
	{
		if(!all(is.na(x)))
	 	{
	 		res <- max(x, na.rm=TRUE)
	 	} else {
	 		res <- NA
	 	}
	 	return(res)
	}
###***###***###	

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues
fTau <- function(x)
{
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(max(x)!=0)
 			{
 				x <- (1-(x/max(x)))
 				res <- sum(x, na.rm=TRUE)
 				res <- res/(length(x)-1)
 			} else {
 				res <- 0
 			}
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		} 
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	} 
 	return(res)
}
###***###***###

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fGini <- function(x)
{
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(sum(x!=0))
 			{
 				res <- gini(x)*(length(x)/(length(x)-1))
  			} else {
 				res <- 0
 			}
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}	 		
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fTsi <- function(x)
{
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(sum(x!=0))
 			{
 				res <- max(x) / sum(x)
 			} else {
 				res <- 0
 			} 	
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}	
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###	

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Function requires setting of a treshold (rpkm)	
fCounts <- function(x, rpkm)
{
	if(all(!is.na(x)))
 	{
 		res <- length(which(x > rpkm))	
 		if (res > 0)
 		{
 			res <- (1 - res/length(x))*(length(x)/(length(x)-1))  #Modification: To bring to normalized scale
 		}
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###		

###+++###
#Function require a data frame with expression data, and give back a vector with EEi values for each gene
#If expression for one tissue is not known, gene specificity for this gene is NA
fEe <- function(x)
{
	if(!all(is.na(x)))
 	{
 		x <- as.matrix(x)
 		x[x<0] <- NA
 		x <- cbind(x, r=rowSums(x, na.rm=FALSE))
		x <- rbind(x, c=colSums(x, na.rm=TRUE))	
 		x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)])
 						
 		res <- apply(x[-nrow(x),-ncol(x)], c(1), FUN=max)
 		res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale
 	} else {
 		res <- NA
 		print("No data avalable.")
 	}
  	return(res)
}
###***###***###

###+++###
#Function require a data frame with expression data, and give back a vector with PEM scores
#If expression for one tissue is not known, gene specificity for this gene is NA
fPem <- function(x)
{
	if(!all(is.na(x)))
 	{
 		x <- as.matrix(x)
 		x[x<0] <- NA
 		x <- cbind(x, r=rowSums(x, na.rm=FALSE)) #Add column with expression of gene per tissue
		x <- rbind(x, c=colSums(x, na.rm=TRUE))	#Add row with expression of all genes in a given tissue
 		x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)]) #calculate the score
 		
 		x[x<1] <- 1
 		x <- log10(x)
 		 		
 		x<- abs(x)				
 		res <- apply(x[-nrow(x),-ncol(x)], c(1), FUN=max) #choose only the maximal score for each gene
 		res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale from 0 to 1
 	} else {
 		res <- NA
 		print("No data avalable.")
 	}
  	return(res)
}
###***###***###

###+++###
#Hg entropy
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fHg <- function(x)
{	
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(sum(x) !=0)
 			{
 				p <- x / sum(x)
 				res <- -sum(p*log2(p), na.rm=TRUE)
 				res <- 1 - (res/log2(length(p))) #Modification: To bring to normalized scale
 			} else {
 				res <- 0
 			} 		
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###	

###+++###
#Z-score
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fZ <- function(x)
{	
	if(all(!is.na(x)))
 	{
 	res <-  apply(scale(t(x), center=TRUE, scale=TRUE),2,max)/((length(x[1,])-1)/sqrt(length(x[1,])))
 	res[is.na(res)] <- 0
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###	

###+++###
#SPM score from TISGED
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fSpm <- function(x)
{	
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{	
 			if(sum(x) !=0)
 			{
 				spm <- x^2/(x%*%x)
 				res <- max(spm) #Modification:To bring to normalized scale. Choose max
 			} else {
 				res <- 0
 			}	 		
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###

###+++###
#Function to calculate and draw correlation between two data sets
##p=parameter names #organisms = used organisms #datasets = used data sets #textPos = text position in smoothScatter, orth= file with orthology, same order as in organisms
fScatDataSet <- function(p, organisms, datasets, add, textPos, orth) 
{
	x1 <- read.table(paste(folder, organisms[1], datasets[1], "TScomparisonTable_9_", add[1],".txt", sep=""), header = TRUE, sep=" ")
	x2 <- read.table(paste(folder, organisms[2], datasets[2], "TScomparisonTable_9_",  add[2], ".txt", sep=""), header = TRUE, sep=" ")

	x1 <- x1[,c("Ensembl.Gene.ID", p)]	
	colnames(x1) <- c("Ensembl.Gene.ID", paste(add[1], p, sep="."))
	x2 <- x2[,c("Ensembl.Gene.ID", p)]	
	colnames(x2) <- c("Ensembl.Gene.ID", paste(add[2], p, sep="."))
	if (!is.na(orth)){
		orthologs <- read.table(paste(folder, orth, ".txt", sep=""), header = TRUE, sep=",")
		orthologs <- orthologs[regexpr("one2one", orthologs$Homology.Type)>0, ]
		orthologs <- orthologs[,c(1,2)]
		x <- merge(orthologs, x1, by=c("Ensembl.Gene.ID"))
		x <- x[,-1]
		colnames(x) <- c("Ensembl.Gene.ID", paste(add[1], p, sep="."))
		x <- merge(x, x2, by=c("Ensembl.Gene.ID"))

	} else {
		x <- merge(x1, x2, by=c("Ensembl.Gene.ID"))	
	}
	
	dev.new(height=12, width=12)
	par(mfrow=c(3,3), cex.main=2.5, cex.axis=1.2, cex.lab=1.4, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
	
	for(i in 1:length(p)) {
		smoothScatter(x[,paste(add[1],".", p[i], sep="")], x[,paste(add[2],".", p[i], sep="")],  xlab=paste(p[i], " in ", add[1], " tissues",sep=""), 			ylab=paste( p[i], " in ", add[2], " tissues",sep=""), nrpoint=Inf, cex=1, nbin=100, xlim=c(0,1), ylim=c(0,1))
		c <- cor.test(x[,paste(add[1], ".", p[i], sep="")], x[,paste(add[2], ".", p[i], sep="")], method="spearman")
		print(c)
		c <- round(c$estimate, digits=2)
		text(x=textPos[1]+0.02, y=textPos[2], pos=2, cex=1.2,labels=paste("rho = ", c,sep=""), col="red", font=2)
		c <- cor.test(x[,paste(add[1], ".", p[i], sep="")], x[,paste(add[2], ".", p[i], sep="")], method="pearson")
		print(c)
		c <- round(c$estimate, digits=2)
		text(x=textPos[1]+0.02, y=(textPos[2]+0.05), pos=2, cex=1.2,labels=paste("   R = ", c,sep=""), col="red", font=2)
	}
	#title(paste("Correlation between specificity parameters in ", organisms[1], " ", add[1], " and ", organisms[2], " ", add[2], " tissues", sep=""), outer=TRUE, line=-2)
	title(paste(" ", sep=""), outer=TRUE, line=-2)
	dev.copy2pdf(device=quartz, file=paste(folder, organisms[1], datasets[1], organisms[2], datasets[2], "TScomp_ScatPlot_9_", add[1], "_",add[2],".pdf", sep=""),onefile=TRUE)#,paper="A4r"
		#dev.off()
	
	xh <- as.matrix(x[,paste(add[1], p, sep=".")])
	xm <- as.matrix(x[,paste(add[2], p, sep=".")])
	xs <- cor(xh, xm, method="spearman")
	xp <- cor(xh, xm, method="pearson")
	capture.output(c("Spearman correlation"),file=paste(folder, organisms[1], datasets[1], organisms[2], datasets[2], "CorrelationTS_", add[1], "_", add[2], ".txt", sep="")) 
	capture.output(xs, append=TRUE, file=paste(folder, organisms[1], datasets[1], organisms[2], datasets[2], "CorrelationTS_", add[1], "_", add[2], ".txt", sep="")) 
	capture.output(c("Pearson correlation"), append=TRUE, file=paste(folder, organisms[1], datasets[1], organisms[2], datasets[2], "CorrelationTS_", add[1], "_", add[2], ".txt", sep="")) 
	capture.output(xp, append=TRUE, file=paste(folder, organisms[1], datasets[1], organisms[2], datasets[2], "CorrelationTS_", add[1], "_", add[2], ".txt", sep="")) 
	capture.output(c("Orthologous genes"), append=TRUE, file=paste(folder, organisms[1], datasets[1], organisms[2], datasets[2], "CorrelationTS_", add[1], "_", add[2], ".txt", sep="")) 
	capture.output(paste(length(x[,1])), append=TRUE, file=paste(folder, organisms[1], datasets[1], organisms[2], datasets[2], "CorrelationTS_", add[1], "_", add[2], ".txt", sep="")) 
	
	return()
}	
###***###***###

###+++###
#Function to calculate and draw correlation for one parameter (x axis) and 8 other parameters
##x = data set, p0 = paramerter to compare, p=parameters names #tPos = text position in smoothScatter(x,y,pos,offset), n = number of graphs in (a,b) format
fScatPlot <- function(x, p0, p, add, tPos) 
{
	dev.new(height=12, width=16)
	par(mfrow=c(3,3), cex.main=0.95, cex.axis=1.2, cex.lab=1.4, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])

	for(i in 1:length(p)) {
		smoothScatter(x[,p0], x[,p[i]], xlab=p0, ylab=p[i], nrpoint=Inf, cex=1, nbin=100, xlim=c(0,1), ylim=c(0,1))
		c <- cor(x[,p0], x[,p[i]], method="spearman")
		c <- round(c, digits=2)
		text(x=tPos[1], y=tPos[2], pos=tPos[3], offset=tPos[4], cex=1.2,labels=paste("rho = ", c,sep=""), col="red", font=2)
		c <- cor(x[,p0], x[,p[i]], method="pearson")
		c <- round(c, digits=2)
		text(x=tPos[1], y=tPos[2]+0.05, pos=tPos[3], offset=tPos[4], cex=1.2,labels=paste("   R = ", c,sep=""), col="red", font=2)
	}
	dev.copy2pdf(device=quartz, file=paste(folder, organism, expDataSource, "TScomp_ScatPlot_9_", p0,  "_", add,".pdf", sep=""),onefile=TRUE)#,paper="A4r"
		#dev.off()
	
	return()
}	
###***###***###

###+++###
#Function to calculate and draw correlation for one parameter (y axis) and 9 other parameters
##x = data set, p0 = paramerter to compare, p=parameters names #tPos = text position in smoothScatter(x,y,pos,offset), n = number of graphs in (a,b) format
fScatPlot2 <- function(x, p0, p, add, tPos, limX) 
{
	dev.new(height=12, width=12)
	par(mfrow=c(3,3), cex.main=0.95, cex.axis=1.2, cex.lab=1.4, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])

	for(i in 1:length(p)) {
		smoothScatter(x[,p0], x[,p[i]],  xlab=p0, ylab=p[i], nrpoint=Inf, cex=1, nbin=100, xlim=c(0,limX), ylim=c(0,1))
		c <- cor.test(x[,p0], x[,p[i]], method="spearman")
		print(c)
		c <- round(c$estimate, digits=2)
		text(x=tPos[1], y=tPos[2], pos=tPos[3], offset=tPos[4], cex=1.2,labels=paste("rho = ", c,sep=""), col="red", font=2)
	}
	dev.copy2pdf(device=quartz, file=paste(folder, organism, expDataSource, "TScomp_ScatPlot_9_", p0,  "_", add,".pdf", sep=""),onefile=TRUE)#,paper="A4r"
		#dev.off()
	
	return()
}	
###***###***###

###+++###
#Calculate and save tissue specificity parameters
#orgExpression = data set, rpkm = cutt off, add = number of tissues, tNames = tissues to use, tNamesNew = tissues to name, RNAseq = how to normalise (log_QN, QN, log, NA)
#Only genes with Ensembl IDs are used, or for Drosophila
#Normalization is done on all tissues, not dependent which tissues are selected later
#1. Data are normalized 
#2. All expression under rpkm is set to 0
#3. Replicates mean is calculated (fReplicateMean)
#4. Genes that not expressed in any tissue are removed
#5. Tissue specificity parameters are calculated
fTS <- function(orgExpression, rpkm, add, tNames, tNamesNew, RNAseq)
{	
	orgExpression <- orgExpression[regexpr("ENS", orgExpression$Ensembl.Gene.ID)>0 | regexpr("FBgn", orgExpression$Ensembl.Gene.ID)>0 | regexpr("PPAG", orgExpression$Ensembl.Gene.ID)>0, ]
	orgExpression <- na.omit(orgExpression)
	print(summary(orgExpression))
	if(RNAseq == "log_QN"){
		x <- orgExpression[,c(-1)]
		x[x < rpkm] <- 1
		x <- log2(x)
		rpkm <- log2(rpkm)
		orgExpression[,c(-1)] <- fQN(x)
	} else if (RNAseq == "QN")  {
		x <- orgExpression[,c(-1)]
		x[x < rpkm] <- 0
		orgExpression[,c(-1)] <- fQN(x)
	} else if (RNAseq == "log")  {
		x <- orgExpression[,c(-1)]
		x[x < rpkm] <- 1
		orgExpression[,c(-1)] <- log2(x)
		rpkm <- log2(rpkm)
	} else {
		x <- orgExpression[,c(-1)]
		x[x < rpkm] <- 0
		orgExpression[,c(-1)] <- x
	}
	
	orgExpression <- fReplicateMean(orgExpression, expDataSource, organism, paste("Averaged.RPKM.",tissuesNames, sep=""))
	orgExpression$Max <- apply(orgExpression[,c(-1)], c(1), fmax)
	orgExpression <- orgExpression[orgExpression$Max > rpkm,]
	orgExpression <- orgExpression[,c(-length(colnames(orgExpression)))]
	print(summary(orgExpression))
	fPlotExpression(orgExpression[,-1], paste("Normalized expression (cutoff", 2^rpkm, "RPKM)", sep=" "), paste("NormalizedQN_", 2^rpkm,"RPKM", sep=""), tissuesPrintNames)
	
	orgExpression <- orgExpression[,c("Ensembl.Gene.ID", paste("Averaged.RPKM.", tNames,sep="")) ]
	colnames(orgExpression) <- c("Ensembl.Gene.ID", paste("Averaged.RPKM.", tNamesNew, sep=""))
	nTissues <- length(tNamesNew)
	tissuesNames <- tNamesNew
	print(paste("Analysis done on", nTissues, "tissue:", sep=" "))
	print(tissuesNames)

	orgExpression$Tau <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fTau)
	orgExpression$Gini <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fGini)
	orgExpression$Tsi <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fTsi)
	orgExpression$Counts <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, function(x){x <- fCounts(x, rpkm)})
	orgExpression$Hg <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fHg)
	orgExpression$Zscore <- fZ(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))])
	orgExpression$Spm <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fSpm)
	orgExpression$Ee <- fEe(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))])
	orgExpression$Pem <- fPem(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))])
	
	orgExpression$Mean <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fmean)
	orgExpression$Max <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fmax)

	print(summary(orgExpression))
	
	p <- c("Tau", "Gini", "Tsi", "Counts", "Ee", "Hg", "Zscore", "Spm", "Pem")
	x <- as.matrix(orgExpression[,p])
	xs <- cor(x, method="spearman")
	xp <- cor(x, method="pearson")
	capture.output(c("Spearman correlation"),file=paste(folder, organism, expDataSource,"CorrelationTS_", add, ".txt", sep="")) 
	capture.output(xs, append=TRUE, file=paste(folder, organism, expDataSource, "CorrelationTS_", add, ".txt", sep="")) 
	capture.output(c("Pearson correlation"), append=TRUE, file=paste(folder, organism, expDataSource,"CorrelationTS_", add, ".txt", sep="")) 
	capture.output(xp, append=TRUE, file=paste(folder, organism, expDataSource,"CorrelationTS_", add, ".txt", sep="")) 

	dev.new(height=9, width=12)
		par(cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
		palette(rev(rich.colors(10)))
		#palette(rev(blues9))
	
		plot(density(orgExpression[,"Tau"],n=1000), main = " ", xlab="Tissue specificity",col=(1), lwd=4, lty=1
		,ylim=c(0,8), xlim=c(-0.1,1.1)
		)
		lines(density(orgExpression[,"Gini"],n = 1000), col=(2), lwd=4, lty=2)
		lines(density(orgExpression[,"Tsi"],n = 1000), col=(3), lwd=4, lty=1)
		lines(density(orgExpression[,"Counts"],n = 1000), col=(4), lwd=4, lty=2)
		lines(density(orgExpression[,"Ee"],n = 1000), col=(5), lwd=4, lty=1)
		lines(density(orgExpression[,"Hg"],n = 1000), col=(6), lwd=4, lty=2)
		lines(density(orgExpression[,"Zscore"],n = 1000), col=(7), lwd=4, lty=1)
		lines(density(orgExpression[,"Spm"],n = 1000), col=(8), lwd=4, lty=2)
		lines(density(orgExpression[,"Pem"],n = 1000), col=(9), lwd=4, lty=1)
				
		legend("topright",c("Tau", "Gini", "TSI", "Counts", "EE", "Hg", "Zscore", "SPM", "PEM"),col=(1:11), lwd=4, lty=c(1,2), bty="n", seg.len=4)
		
		dev.copy2pdf(device=quartz, file=paste(folder, organism, expDataSource, "TScomparison_9_",  add,".pdf", sep=""),onefile=TRUE)#,paper="A4r"
		#dev.off()

	write.table(orgExpression, file=paste(folder, organism, expDataSource,"TScomparisonTable_9_",  add,".txt",sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)
		
	fScatPlot(orgExpression, "Tau", c("Gini", "Tsi", "Counts", "Hg", "Zscore","Spm", "Ee", "Pem"), add, c(0, 0.94, 4, 0.5))
	fScatPlot2(orgExpression, "Mean", c("Tau", "Gini", "Tsi", "Counts", "Hg", "Zscore", "Spm", "Ee", "Pem"), add, c(ceiling(max(orgExpression$Mean)), 0.95, 2, 0.5), ceiling(max(orgExpression$Mean)))
	fScatPlot2(orgExpression, "Max", c("Tau", "Gini", "Tsi", "Counts", "Hg", "Zscore", "Spm", "Ee", "Pem"), add, c(ceiling(max(orgExpression$Max)), 0.95, 2, 0.5), ceiling(max(orgExpression$Max)))

	return()
}
###***###***###


###################################################
#Calculate gene specificity with different methods#
###################################################

###########
###Run#####
orgExpression <- fInputData()
rpkm <- 1

#For Mouse ENCODE
organism <- "Mus"
expDataSource <- "ENCODE"
tissuesNames <- fTissueNames (organism, expDataSource)
tissuesPrintNames <- fTissuePrintNames (organism, expDataSource)
nTissues <- length(tissuesNames)
orgExpression <- fInputData()
fTS(orgExpression, 1, "22", tissuesNames, tissuesNames, "log") #All ENCODE tissues 
fTS(orgExpression, 1, "16F", tissuesEFNamesMouse, tissuesEFNames, "log") #ENCODE tissues that are common with human Fagerberg  

fTS(orgExpression, 1, "22notLog", tissuesNames, tissuesNames, "") #All ENCODE tissues without log-transformation 
fTS(orgExpression, 1, "16FnotLog", tissuesEFNamesMouse, tissuesEFNames, "") #ENCODE tissues that are common with human Fagerberg without log-transformation 

fTS(orgExpression, 1, "22QN", tissuesNames, tissuesNames, "log_QN") #All ENCODE tissues with quantile normalisation 

#For Mouse Brawand
organism <- "Mus"
expDataSource <- "Brawand"
tissuesNames <- fTissueNames (organism, expDataSource)
tissuesPrintNames <- fTissuePrintNames (organism, expDataSource)
nTissues <- length(tissuesNames)
orgExpression <- fInputData()
fTS(orgExpression, 1, "6m", tissuesNames, tissuesNames, "log")	#All Brawand tissues for Mouse 

#For Mouse Microarray Thorrez
organism <- "Mus"
expDataSource <- "Thorrez"
tissuesNames <- fTissueNames (organism, expDataSource)
tissuesPrintNames <- fTissuePrintNames (organism, expDataSource)
nTissues <- length(tissuesNames)
orgExpression <- fInputData()
fTS(orgExpression, 2, "19mT", tissuesNames, tissuesNames, "") #All Microarray Thorrez tissues
fTS(orgExpression, 2, "14mT", tissuesGTNames, tissuesGTNames, "") #Microarray Thorrez tissues common with human Ge data set 


#For Human Fageberg
organism <- "Hum"
expDataSource <- "Fagerberg"
tissuesNames <- fTissueNames (organism, expDataSource)
tissuesPrintNames <- fTissuePrintNames (organism, expDataSource)
nTissues <- length(tissuesNames)
orgExpression <- fInputData()
fTS(orgExpression, 1, "27", tissuesNames, tissuesNames, "log")	#All Fagerberg tissues 
fTS(orgExpression, 1, "16E", tissuesEFNamesHuman, tissuesEFNames, "log") #Fagerberg tissues common with ENCODE mouse dataset

fTS(orgExpression, 1, "27notLog", tissuesNames, tissuesNames, "")	#All Fagerberg tissues without log-transformation
fTS(orgExpression, 1, "16EnotLog", tissuesEFNamesHuman, tissuesEFNames, "") #Fagerberg tissues common with ENCODE mouse dataset without log-transformation 

fTS(orgExpression, 1, "27QN", tissuesNames, tissuesNames, "log_QN")	#All Fagerberg tissues with quantile normalisation 

fTS(orgExpression, 0.1, "27notLog01RPKM", tissuesNames, tissuesNames, "")	#All Fagerberg tissues without log-transformation with 0.1 RPKM cutoff 


#For Human Brawand
organism <- "Hum"
expDataSource <- "Brawand"
tissuesNames <- fTissueNames (organism, expDataSource)
tissuesPrintNames <- fTissuePrintNames (organism, expDataSource)
nTissues <- length(tissuesNames)
orgExpression <- fInputData()
fTS(orgExpression, 1, "6h", tissuesBrHNames, tissuesBrNames, "log") #Brawand tissues common to other organisms Brawand 


#For Human Microarray Ge
organism <- "Hum"
expDataSource <- "Ge"
tissuesNames <- fTissueNames (organism, expDataSource)
tissuesPrintNames <- fTissuePrintNames (organism, expDataSource)
nTissues <- length(tissuesNames)
orgExpression <- fInputData()
fTS(orgExpression, 2, "32hG", tissuesNames, tissuesNames, "") #All Microarray Ge tissues 
fTS(orgExpression, 2, "14hG", tissuesGTNames, tissuesGTNames, "") #Microarray Ge tissues common to mouse Thorrez tissues 	


#Run for TScomparison
tsParameters <- c("Tau", "Gini", "Tsi", "Counts", "Hg", "Zscore", "Spm", "Ee", "Pem")

fScatDataSet(tsParameters, c("Hum", "Hum"), c("Fagerberg", "Fagerberg"), c("27", "16E"), c(0.95,0.05), NA) # F S8
fScatDataSet(tsParameters, c("Mus", "Mus"), c("ENCODE", "ENCODE"), c("22", "16F"), c(0.95,0.05), NA) #F S9
fScatDataSet(tsParameters, c("Hum", "Mus"), c("Fagerberg", "ENCODE"), c("16E", "16F"), c(0.95,0.05), "HumOrthologsMusEnsV75") #F 3
fScatDataSet(tsParameters, c("Hum", "Mus"), c("Brawand", "Brawand"), c("6h", "6m"), c(0.95,0.05), "HumOrthologsMusEnsV75") #F S14

fScatDataSet(tsParameters, c("Hum", "Hum"), c("Fagerberg", "Ge"), c("27", "32hG"), c(0.95,0.05), NA) #F 6	
fScatDataSet(tsParameters, c("Mus", "Mus"), c("ENCODE", "Thorrez"), c("22", "19mT"), c(0.95,0.05), NA) #F S55
fScatDataSet(tsParameters, c("Hum", "Hum"), c("Ge", "Ge"), c("32hG", "14hG"), c(0.95,0.05), NA) #F S50
fScatDataSet(tsParameters, c("Mus", "Mus"), c("Thorrez", "Thorrez"), c("19mT", "14mT"), c(0.95,0.05), NA) #F S51
fScatDataSet(tsParameters, c("Hum", "Mus"), c("Ge", "Thorrez"), c("14hG", "14mT"), c(0.95,0.05), "HumOrthologsMusEnsV75") #F S54
fScatDataSet(tsParameters, c("Mus", "Mus"), c("ENCODE", "Brawand"), c("22", "6m"), c(0.95,0.05), "MusOrthologsMusEnsV75") #F S57
fScatDataSet(tsParameters, c("Hum", "Hum"), c("Fagerberg", "Brawand"), c("27", "6h"), c(0.95,0.05), "HumOrthologsHumEnsV75") #F S56


fScatDataSet(tsParameters, c("Mus", "Mus"), c("ENCODE", "ENCODE"), c("22notLog", "16FnotLog"), c(0.95,0.05), NA) #F S70
fScatDataSet(tsParameters, c("Hum", "Hum"), c("Fagerberg", "Fagerberg"), c("27notLog", "16EnotLog"), c(0.95,0.05), NA) #F S69    
fScatDataSet(tsParameters, c("Hum", "Hum"), c("Fagerberg", "Fagerberg"), c("27", "27notLog"), c(0.95,0.05), NA) #F S72
fScatDataSet(tsParameters, c("Mus", "Mus"), c("ENCODE", "ENCODE"), c("22", "22notLog"), c(0.95,0.05), NA) #F S73


fScatDataSet(tsParameters, c("Hum", "Hum"), c("Fagerberg", "Fagerberg"), c("27", "27notLogQN"), c(0.95,0.05), NA) 
fScatDataSet(tsParameters, c("Mus", "Mus"), c("ENCODE", "ENCODE"), c("22", "22notLogQN"), c(0.95,0.05), NA) 

fScatDataSet(tsParameters, c("Hum", "Hum"), c("Fagerberg", "Fagerberg"), c("27", "27QN"), c(0.95,0.05), NA) #F S74
fScatDataSet(tsParameters, c("Mus", "Mus"), c("ENCODE", "ENCODE"), c("22", "22QN"), c(0.95,0.05), NA) #F S75

fScatDataSet(tsParameters, c("Hum", "Mus"), c("Fagerberg", "ENCODE"), c("16EnotLog", "16FnotLog"), c(0.95,0.05), "HumOrthologsMusEnsV75") #F S71


#Comparison of different count tresholds #TC
organism <- "Hum"
dataset <- "Fagerberg"
add <- "27notLog01RPKM"
tissuesNames <- fTissueNames (organism, expDataSource)
tissuesPrintNames <- fTissuePrintNames (organism, expDataSource)
nTissues <- length(tissuesNames)


x <- read.table(paste(folder, organism, dataset, "TScomparisonTable_9_", add,".txt", sep=""), header = TRUE, sep=" ")
x$Counts1 <- apply(x[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, function(x){x <- fCounts(x, 1)})
x$Counts10 <- apply(x[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, function(x){x <- fCounts(x, 10)})
x$Counts100 <- apply(x[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, function(x){x <- fCounts(x, 100)})

x <- x[,c("Counts", "Counts1", "Counts10", "Counts100")]
colnames(x) <- c("0.1", "1", "10", "100")
names <- c("0.1 RPKM", "1 RPKM", "10 RPKM", "100 RPKM")

dev.new(height=9, width=12)
	par(cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
	#palette(rev(rich.colors(ncol(x))))
	palette(rev(heat.colors(8)))
	
	plot(density(x[,1],n=1000), main = "", xlab="Counts",col=(2), lwd=3, ylim=c(0,15))
	for(i in c(2:length(names)))
	{	
		lines(density(x[,i],n = 1000), col=(i*2), lwd=3)
	}
	legend("topright", names, col=c(2,4,6,8), lty="solid", lwd=3)
	
	dev.copy2pdf(device=quartz, file=paste(folder, organism, "CountsDifferentThreshholds", expDataSource, "", add, ".pdf", sep=""),onefile=TRUE)#,paper="A4r"
	#dev.off()

###############
###############

#Comparison of RNA-seq and Microarray

organism <- "Hum"
dataset <- "Fagerberg"
add <- "27"

x <- read.table(paste(folder, organism, dataset, "TScomparisonTable_9_", add,".txt", sep=""), header = TRUE, sep=" ")
x2 <- read.table(paste(folder, "Hum", "Ge", "TScomparisonTable_9_", "32hG",".txt", sep=""), header = TRUE, sep=" ")

x3 <- read.table(paste(folder, "Mus", "ENCODE", "TScomparisonTable_9_", "22",".txt", sep=""), header = TRUE, sep=" ")
x4 <- read.table(paste(folder, "Mus", "Thorrez", "TScomparisonTable_9_", "19mT",".txt", sep=""), header = TRUE, sep=" ")

xnM <- x[!(x$Ensembl.Gene.ID %in% x2$Ensembl.Gene.ID),] #in human RNA-seq not Microarray
xM <- x[x$Ensembl.Gene.ID %in% x2$Ensembl.Gene.ID,] #in human RNA-seq and Microarray

xnMm <- x3[!(x3$Ensembl.Gene.ID %in% x4$Ensembl.Gene.ID),] #in mouse RNA-seq not Microarray
xMm <- x3[x3$Ensembl.Gene.ID %in% x4$Ensembl.Gene.ID,] #in mouse RNA-seq and Microarray

names <- c("Tau", "Gini")
xnM <- xnM[,names]
xM <- xM[,names]
xnMm <- xnMm[,names]
xMm <- xMm[,names]

dev.new(height=9, width=12)
	par(cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
	palette(rev(rich.colors(10)))
	
	plot(density(xM[,1],n=1000), main = "", xlab="Tau",col=(1), lwd=3, ylim=c(0,7.5))
		
	lines(density(xnM[,1],n = 1000), col=(2), lwd=3)
	lines(density(xMm[,1],n = 1000), col=(8), lwd=3)
	lines(density(xnMm[,1],n = 1000), col=(7), lwd=3)
	
	legend("topleft", c(paste("Human Tau, same genes as microarray data set, ", length(xM[,1]), " genes", sep=""), paste("Human Tau, genes not detected by microarray, ", length(xnM[,1]), " genes", sep=""), paste("Mouse Tau, same genes as microarray data set, ", length(xMm[,1]), " genes", sep=""), paste("Mouse Tau, genes not detected by microarray, ", length(xnMm[,1]), " genes",sep="")), col=c(1,2,8,7), lty="solid", lwd=3)
	
	dev.copy2pdf(device=quartz, file=paste(folder, "TauRNAseqMicroarrayComparison", ".pdf", sep=""),onefile=TRUE)#,paper="A4r"
	#dev.off()

###############
###############

fGO("Hum","Fagerberg", "27")
fGO("Mus", "ENCODE", "22")
#############################
##### Comparison GO terms
fGO <- function(organism, dataset, add) #TC
{
	x <- read.table(paste(folder, organism, dataset, "TScomparisonTable_9_", add, ".txt", sep=""), header = TRUE, sep=" ")

	xS <- read.table(paste(folder, organism, "EnsV75GenelistSperm", ".txt", sep=""), header = TRUE, sep=",") #0007283
	xN <- read.table(paste(folder, organism, "EnsV75GenelistNeuro", ".txt", sep=""), header = TRUE, sep=",") #0050877
	xX <- read.table(paste(folder, organism, "EnsV75GenelistXenobiotic", ".txt", sep=""), header = TRUE, sep=",") #0006805
	xP <- read.table(paste(folder, organism, "EnsV75GenelistProtein", ".txt", sep=""), header = TRUE, sep=",") #GO:0006457
	xM <- read.table(paste(folder, organism, "EnsV75GenelistMembrane", ".txt", sep=""), header = TRUE, sep=",") #GO:0061024
	xR <- read.table(paste(folder, organism, "EnsV75GenelistRNAsplicing", ".txt", sep=""), header = TRUE, sep=",") #GO:0008380
			
xS <- merge(xS, x, by="Ensembl.Gene.ID")
xN <- merge(xN, x, by="Ensembl.Gene.ID")
xX <- merge(xX, x, by="Ensembl.Gene.ID")
xP <- merge(xP, x, by="Ensembl.Gene.ID")
xM <- merge(xM, x, by="Ensembl.Gene.ID")
xR <- merge(xR, x, by="Ensembl.Gene.ID")
summary(xS)
summary(xN)
summary(xX)
summary(xP)
summary(xM)
summary(xR)

tsParameters <- c("Tau", "Gini", "Tsi", "Counts", "Hg", "Zscore", "Spm", "Ee", "Pem")
dev.new(height=9, width=12)
		par(mfrow=c(3,3), cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
for(p in tsParameters){
	
		palette(rev(rich.colors(10)))
			
		plot(density(xS[,p], n=100), main = p, xlab="Tissue specificity",col=(1), lwd=2, lty=1
		,ylim=c(0,10), xlim=c(-0.1,1.1)
		)
		lines(density(xR[,p], n=100), col=(8), lwd=2, lty=1)
		
		lines(density(xM[,p], n=100), col=(8), lwd=2, lty=3)
		lines(density(xP[,p], n=100), col=(8), lwd=2, lty=4)
		
		lines(density(xN[,p], n=100), col=(1), lwd=2, lty=4)
		lines(density(xX[,p], n=100), col=(1), lwd=2, lty=3)
		
		lines(density(x[,p], n=100), col=(10), lwd=2, lty=1)
				
		#legend("topright",c("Tau", "Gini", "TSI", "Counts", "EE", "Hg", "Zscore", "ZUscore", "SPM", "PEM", "PEMU"),col=(1:11), lty="solid", lwd=3)
		legend("topright",c("Spermatogenesis", "Neurological System", "Xenobiotic metabolism", "RNA splicing", "Protein folding", "Membrane organisation", "All"),col=c(1, 1, 1, 8, 8, 8, 10), lwd=2, lty=c(1, 4, 3, 1, 4, 3, 1), bty="n", seg.len=4)
	}	
		dev.copy2pdf(device=quartz, file=paste(folder, organism, "GOTermsComparison",  ".pdf", sep=""),onefile=TRUE)#,paper="A4r"
		#dev.off()
}
############################
###########################

#############################
x <- fSnum("Hum","Fagerberg", "27")
x <- fSnum("Mus", "ENCODE", "22")

#####Number of specific genes per tissues
fSnum <- function(organism, dataset, add) #TC
{
	x <- read.table(paste(folder, organism, dataset, "TScomparisonTable_9_", add, ".txt", sep=""), header = TRUE, sep=" ")
		#x$Counts <- apply(x[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, function(x){x <- fCounts(x, log2(10))})
	x[,c(grep("Averaged.RPKM.", colnames(x)))] <- ifelse(x[,c(grep("Averaged.RPKM.", colnames(x)))]==x$Max, 1, 0 )
	
	x$Organ <- apply(x[,c(grep("Averaged.RPKM.", colnames(x)))], c(1), function(x){x <- which(x>0)})
	x$Organ <- sapply(x$Organ, function(x){x <- as.numeric(x[1])})
	names <- gsub("Averaged.RPKM.", "", colnames(x[,c(grep("Averaged.RPKM.", colnames(x)))]))
	x$Organ <- names[x$Organ]

	
tsParameters <- c("Tau", "Gini", "Tsi", "Counts", "Hg", "Zscore", "Spm", "Ee", "Pem")
dev.new(height=9, width=12)
		trellis.par.set(list(background=list(col=my.col[1]), add.text=list(col=my.col[2], cex=1), axis.line=list(col=my.col[2]), axis.text=list(col=my.col[2], cex=0.6, font=1), 	par.main.text=list(col=my.col[2], cex=0.6), par.xlab.text=list(col=my.col[2], cex=0.6, font=1), par.ylab.text=list(col=my.col[2], cex=0.4, font=1), plot.line=list(col=my.col[2]), dot.line=list(lwd=1, lty=2, col="#4B4B4B"), strip.background=list(col=my.col[1]))) #trellis.par.get()
for(p in tsParameters){
	
		temp <- x[x[,p]>0.8,c(p,"Organ")]
		temp2 <- data.frame(Organ=names, Number=0)
		for(n in names){
			temp2[temp2[,1]==n,2] <- length(temp[temp$Organ==n,1])
		}
		#temp2 <- temp2[order(temp2[,2]),]
		palette(rev(rich.colors(10)))
		myplot <- barchart(temp2[,1]~temp2[,2], horizontal=TRUE, main=p, xlab="", xlim=c(0,2000), col="blue", border=c("transparent")
		, par.settings=list(layout.widths=list(ylab.axis.padding=0.1, xlab.axis.padding=0.1), axis.components=list(top=list(pad1=0, pad2=0), bottom=list(pad1=0.5, pad2=0)))
		)
		assign(paste("plot",p,sep=""),myplot)
			
		}	
		
		grid.arrange(plotTau, plotGini, plotTsi, plotCounts, plotHg, plotZscore, plotSpm, plotEe, plotPem, ncol=3, nrow=3)
				
		dev.copy2pdf(device=quartz, file=paste(folder, organism, "SpecificGeneNumbers",  ".pdf", sep=""),onefile=TRUE)#,paper="A4r"
		#dev.off()
	return(x)	
	
}
############################
###########################
#####################

f5comb("Hum","Fagerberg", "27")
f5comb("Mus", "ENCODE", "22")

f5comb("Hum", "Ge", "32hG")
f5comb("Mus", "Thorrez", "19mT")
#############################
#####Boxplot of all 5 combination of tissues

f5comb <- function(organism, dataset, add) #TC
{
	x <- read.table(paste(folder, organism, dataset, "TScomparisonTable_9_", add, ".txt", sep=""), header = TRUE, sep=" ")
	xOriginal <- x
	x <- x[,c(grep("Averaged.RPKM.", colnames(x)))]
	tsParameters <- c("Tau", "Gini", "Tsi", "Counts", "Hg", "Zscore", "Spm", "Ee", "Pem")
	sampleNumber <- 1000
		
	combTissues <- matrix(ncol = sampleNumber, nrow = 5)
	for(i in 1:sampleNumber){
		combTissues[,i] <- sample(length(colnames(x)), 5)
	}
		
	combRes <- data.frame("tissues"=0, "Tau"=0, "Gini"=0, "Tsi"=0, "Counts"=0, "Hg"=0, "Zscore"=0, "Spm"=0, "Ee"=0, "Pem"=0)
			
	for(t in 1:length(combTissues[1,])){
	
		x2 <- x[,combTissues[,t]]
		
		x2$Tau <- apply(x2[,c(grep("Averaged.RPKM.", colnames(x2)))], 1, fTau)
		x2$Gini <- apply(x2[,c(grep("Averaged.RPKM.", colnames(x2)))], 1, fGini)
		x2$Tsi <- apply(x2[,c(grep("Averaged.RPKM.", colnames(x2)))], 1, fTsi)
		x2$Counts <- apply(x2[,c(grep("Averaged.RPKM.", colnames(x2)))], 1, function(x){x <- fCounts(x, 0)})
		x2$Hg <- apply(x2[,c(grep("Averaged.RPKM.", colnames(x2)))], 1, fHg)
		x2$Zscore <- fZ(x2[,c(grep("Averaged.RPKM.", colnames(x2)))])
		x2$Spm <- apply(x2[,c(grep("Averaged.RPKM.", colnames(x2)))], 1, fSpm)
		x2$Ee <- fEe(x2[,c(grep("Averaged.RPKM.", colnames(x2)))])
		x2$Pem <- fPem(x2[,c(grep("Averaged.RPKM.", colnames(x2)))])
			
	
		combRes[t,1] <- t
		for(p in tsParameters){
			r <- cor(xOriginal[,p], x2[,p], method="pearson")
			combRes[t,p] <- r
		} 
	}

	write.table(combRes, file=paste(folder, organism, add,"_5TissuesCombinationsPSample", sampleNumber, ".txt",sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)
	write.table(combTissues, file=paste(folder, organism, add,"_5TissuesCombinationsTissuesPSample", sampleNumber, ".txt",sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)
	
	dev.new(height=8, width=16)
	boxplot(combRes[,-1], notch=TRUE, outline=TRUE, col="blue", ylim=c(-0.6,0.8))
	abline(h=0, lty=3)

	dev.copy2pdf(device=quartz, file=paste(folder, organism, add,"_5TissueCombinationsPSample", sampleNumber,  ".pdf", sep=""),onefile=TRUE)#,paper="A4r"
		#dev.off()	
}
############################
###########################


fTauComb("Hum", "Fagerberg", "27")
fTauComb("Mus", "ENCODE", "22")

fTauComb("Hum", "Ge", "32hG")
fTauComb("Mus", "Thorrez", "19mT")
#############################
#####Boxplot of all combinations of tissues for Tau
fTauComb <- function(organism, dataset, add) #TC
{
	x <- read.table(paste(folder, organism, dataset,"TScomparisonTable_9_", add, ".txt", sep=""), header = TRUE, sep=" ")
	xOriginal <- x
	x <- x[,c(grep("Averaged.RPKM.", colnames(x)))]
	tissues <- length(c(grep("Averaged.RPKM.", colnames(x))))
	tsParameters <- c(5:(as.numeric(tissues)-1))
	sampleNumber <- 1000
	
	combRes <- matrix(ncol = length(colnames(x))-4, nrow = sampleNumber)
	r <- 1
	for(n in c(5:(as.numeric(tissues)-1))){
		combTissues <- matrix(ncol = sampleNumber, nrow = n)
		for(i in 1:sampleNumber){
			combTissues[,i] <- sample(length(colnames(x)), n)
		}
		for(t in 1:length(combTissues[1,])){
			x2 <- x[,combTissues[,t]]
			x2$Tau <- apply(x2[,c(grep("Averaged.RPKM.", colnames(x2)))], 1, fTau)
				
			combRes[t,1] <- t
			res <- cor(xOriginal[,"Tau"], x2[,"Tau"], method="pearson")
			combRes[t,r+1] <- res	 
		}
		r <- r+1
		write.table(combTissues, file=paste(folder, organism, add,"_TauTissuesCombinations", n, "TissuesPSample", sampleNumber, ".txt",sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)
	}
	colnames(combRes) <- c("tissues", 5:(length(colnames(x))-1))
	write.table(combRes, file=paste(folder, organism, add,"_TauTissuesCombinationsPSample", sampleNumber, ".txt",sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)
	
	dev.new(height=8, width=16)
	
	boxplot(combRes[,-1], notch=TRUE, outline=TRUE, col="blue", ylim=c(0,1))

	dev.copy2pdf(device=quartz, file=paste(folder, organism, add,"_TauTissueCombinationsPSample", sampleNumber,  ".pdf", sep=""),onefile=TRUE)#,paper="A4r"
		#dev.off()
}
############################
###########################


fGOhm(c("Hum", "Mus"),c("Fagerberg", "ENCODE"), c("16E","16F"))
#############################
##### Comparison 
fGOhm <- function(organisms, datasets, adds) #TC
{
	xHum <- read.table(paste(folder, organisms[1], datasets[1], "TScomparisonTable_9_", adds[1], ".txt", sep=""), header = TRUE, sep=" ")
	xMus <- read.table(paste(folder, organisms[2], datasets[2], "TScomparisonTable_9_", adds[2], ".txt", sep=""), header = TRUE, sep=" ")
	
	xHumS <- read.table(paste(folder, organisms[1], "EnsV75GenelistSperm", ".txt", sep=""), header = TRUE, sep=",") #0007283
	xHumN <- read.table(paste(folder, organisms[1], "EnsV75GenelistNeuro", ".txt", sep=""), header = TRUE, sep=",") #0050877
	xHumX <- read.table(paste(folder, organisms[1], "EnsV75GenelistXenobiotic", ".txt", sep=""), header = TRUE, sep=",") #0006805
	xMusS <- read.table(paste(folder, organisms[2], "EnsV75GenelistSperm", ".txt", sep=""), header = TRUE, sep=",") #0007283
	xMusN <- read.table(paste(folder, organisms[2], "EnsV75GenelistNeuro", ".txt", sep=""), header = TRUE, sep=",") #0050877
	xMusX <- read.table(paste(folder, organisms[2], "EnsV75GenelistXenobiotic", ".txt", sep=""), header = TRUE, sep=",") #0006805

	
	xHumP <- read.table(paste(folder, organisms[1], "EnsV75GenelistProtein", ".txt", sep=""), header = TRUE, sep=",") #GO:0006457
	xHumM <- read.table(paste(folder, organisms[1], "EnsV75GenelistMembrane", ".txt", sep=""), header = TRUE, sep=",") #GO:0061024
	xHumR <- read.table(paste(folder, organisms[1], "EnsV75GenelistRNAsplicing", ".txt", sep=""), header = TRUE, sep=",") #GO:0008380
	xMusP <- read.table(paste(folder, organisms[2], "EnsV75GenelistProtein", ".txt", sep=""), header = TRUE, sep=",") #GO:0006457
	xMusM <- read.table(paste(folder, organisms[2], "EnsV75GenelistMembrane", ".txt", sep=""), header = TRUE, sep=",") #GO:0061024
	xMusR <- read.table(paste(folder, organisms[2], "EnsV75GenelistRNAsplicing", ".txt", sep=""), header = TRUE, sep=",") #GO:0008380

	xHumSpec <- rbind(xHumS, xHumN, xHumX)
	xHumSpec <- unique(xHumSpec)
	xHumSpec <- merge(xHumSpec, xHum, by="Ensembl.Gene.ID")
		xHumHK <- rbind(xHumP, xHumM, xHumR)
		xHumHK <- unique(xHumHK)
		xHumHK <- merge(xHumHK, xHum, by="Ensembl.Gene.ID")
			
	xMusSpec <- rbind(xMusS, xMusN, xMusX)
	xMusSpec <- unique(xMusSpec)
	xMusSpec <- merge(xMusSpec, xMus, by="Ensembl.Gene.ID")
		xMusHK <- rbind(xMusP, xMusM, xMusR)
		xMusHK <- unique(xMusHK)
		xMusHK <- merge(xMusHK, xMus, by="Ensembl.Gene.ID")


	write.table(xHumSpec, file=paste(folder, organisms[1], datasets[1], "TScomparisonTable_9_", adds[1], "Spec", ".txt", sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)
	write.table(xMusSpec, file=paste(folder, organisms[2], datasets[2], "TScomparisonTable_9_", adds[2], "Spec", ".txt", sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)
	fScatDataSet(tsParameters, c("Hum", "Mus"), c("Fagerberg", "ENCODE"), c("16ESpec", "16FSpec"), c(0.95,0.05), "HumOrthologsMusEnsV75")
	
	write.table(xHumHK, file=paste(folder, organisms[1], datasets[1], "TScomparisonTable_9_", adds[1], "HK", ".txt", sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)
	write.table(xMusHK, file=paste(folder, organisms[2], datasets[2], "TScomparisonTable_9_", adds[2], "HK", ".txt", sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)
	fScatDataSet(tsParameters, c("Hum", "Mus"), c("Fagerberg", "ENCODE"), c("16EHK", "16FHK"), c(0.95,0.05), "HumOrthologsMusEnsV75")

}
############################
###########################


#############################
fTvenn("Hum","Fagerberg", "27")
fTvenn("Mus", "ENCODE", "22")
	
fTvenn <- function(organism, dataset, add) #TC
{
	threshold <- 0.4
	x <- read.table(paste(folder, organism, dataset, "TScomparisonTable_9_", add, ".txt", sep=""), header = TRUE, sep=" ")
		#x$Counts <- apply(x[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, function(x){x <- fCounts(x, log2(10))})
	x[,c(grep("Averaged.RPKM.", colnames(x)))] <- ifelse(x[,c(grep("Averaged.RPKM.", colnames(x)))]==x$Max, 1, 0 )
	
	x$Organ <- apply(x[,c(grep("Averaged.RPKM.", colnames(x)))], c(1), function(x){x <- which(x>0)})
	x$Organ <- sapply(x$Organ, function(x){x <- as.numeric(x[1])})
	names <- gsub("Averaged.RPKM.", "", colnames(x[,c(grep("Averaged.RPKM.", colnames(x)))]))
	x$Organ <- names[x$Organ]
	
	tsParameters <- list(c("Tau", "Gini", "Counts", "Hg", "Pem") , c("Tau", "Tsi", "Zscore", "Spm", "Ee"), c("Tau", "Gini", "Counts", "Zscore", "Spm"))

	palette(rev(rich.colors(10)))
	
	for (set in c(1,2,3)){
		for(o in names){
			xV <- vector("list", length(tsParameters[[set]]))
			i<-1
			for(p in tsParameters[[set]]){
				temp <- x[x[,p]>threshold & x[,"Organ"]==o, c("Ensembl.Gene.ID")]
				xV[[i]] <- temp
				i <- i+1
			}
			names(xV) <- tsParameters[[set]]
			if(set==1){
				vc <- c(1,2,4,6,9)
			} else if(set==2) {
				vc <- c(1,3,7,8,5)
			} else {
				vc <- c(1,2,4,7,8)
			}
		
			print(venn.diagram(xV,filename=paste(organism,"_Venn_", o, "_set", set, "_", threshold,".png", sep=""), width=2000, height=2000, imagetype="png", col=vc, fill=vc, lwd=3, alpha=0.3, cex=1, cat.cex=1, cat.pos=c(1,-20,-100,120,30), cat.dist=c(0.2, 0.2, 0.27, 0.2, 0.2), main=o, main.pos=c(0.2,1), main.cex=1.5))	
		}
	}
}
############################
###########################
#####################

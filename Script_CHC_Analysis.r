#############################################
##                                         ##
##  @Autor: César Osvaldo Martínez Cantú   ##
##                                         ##
#############################################

#############################
##                         ##
##  Librerías a utilizar   ##
##                         ##
#############################
library(dplyr)
library(Biobase)
library(oligoClasses)
library(ArrayExpress)
library(oligo)
library(arrayQualityMetrics)
library(pd.mogene.2.0.st)
library(mogene20sttranscriptcluster.db)
library(limma)
library(ggplot2)
library(calibrate)
library(plotly)
library(pheatmap)
library(RColorBrewer)
library(mvtnorm)
library(gplots)
library(DAAG)

#########################################
##                                     ##
##  Definiendo ruta de la información  ##
##                                     ##
#########################################
setwd("C:/../../../DataGenomic/analysisPC-PD")
getwd()

##################################################################
##                                                              ##
##  Se carga archivo SDRF con los datos fenotipicos del raton.  ##
##                                                              ##
##################################################################
sdrf_location <- "SDRF.file" 
SDRF <- read.table(sdrf_location, header=T)
print(SDRF)
###########################################################################
##                                                                       ##
##    Source.Name Array.Data.File Factor.Value.phenotype Time.treatment  ##
## 1    PC181.CEL       PC181.CEL                Control        Control  ##
## 2    PC182.CEL       PC182.CEL                Control        Control  ##
## 3    PC183.CEL       PC183.CEL                Control        Control  ##
## 4    PD061.CEL       PD061.CEL                Tratado      6_Semanas  ##
## 5    PD062.CEL       PD062.CEL                Tratado      6_Semanas  ##
## 6    PD063.CEL       PD063.CEL                Tratado      6_Semanas  ##
## 7    PD101.CEL       PD101.CEL                Tratado     10_Semanas  ##
## 8    PD102.CEL       PD102.CEL                Tratado     10_Semanas  ##
## 9    PD103.CEL       PD103.CEL                Tratado     10_Semanas  ##
## 10   PD141.CEL       PD141.CEL                Tratado     14_Semanas  ##
## 11   PD142.CEL       PD142.CEL                Tratado     14_Semanas  ##
## 12   PD143.CEL       PD143.CEL                Tratado     14_Semanas  ##
## 13   PD181.CEL       PD181.CEL                Tratado     18_Semanas  ##
## 14   PD182.CEL       PD182.CEL                Tratado     18_Semanas  ##
## 15   PD183.CEL       PD183.CEL                Tratado     18_Semanas  ##
##                                                                       ##
###########################################################################

#########################################################
##                                                     ##
##  Se les asigna nombre a los renglones del archivo,  ##
##  tomando los valores de la columna Array.Data.File  ##
##                                                     ##
#########################################################
rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)
print(SDRF)
####################################################################################
##                                                                                ##
## An object of class 'AnnotatedDataFrame'                                        ##
##   rowNames: PC181.CEL PC182.CEL ... PD183.CEL (15 total)                       ##
##   varLabels: Source.Name Array.Data.File Factor.Value.phenotype Time.treatment ##
##   varMetadata: labelDescription                                                ##
##                                                                                ##
####################################################################################

#####################################
##                                 ##
##  Extraccion de caracteristicas  ##
##                                 ##
#####################################
raw_data_dir <- c(getwd())
raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir, 
                                SDRF$Array.Data.File), 
                                verbose = FALSE, phenoData = SDRF)
#####################################################################
##                                                                 ##
##  Colocar los archivos .CEL en la misma ruta del SDRF            ##
##  raw_Data contiene los archivos .CEL, los datos de expresiones  ##
##  y los datos fenotipicos en un mismo archivo.                   ##
##                                                                 ##
#####################################################################

###############################################################################
##                                                                           ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PC181.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PC182.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PC183.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD061.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD062.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD063.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD101.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD102.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD103.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD141.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD142.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD143.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD181.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD182.CEL              ##
## Reading in : C:/../../../DataGenomic/analysisPC-PD/PD183.CEL              ##
##                                                                           ##
###############################################################################

##################################################
##                                              ##
##  Se detiene si no es un objeto valido de R.  ##
##                                              ##
##################################################
stopifnot(validObject(raw_data)) 

#############################################
##                                         ##
##  Head me da solo los primeros 5 datos.  ##
##                                         ##
#############################################
head(Biobase::pData(raw_data)) 
###################################################################################
##                                                                               ##
##            Source.Name Array.Data.File Factor.Value.phenotype Time.treatment  ##
##  PC181.CEL   PC181.CEL       PC181.CEL                Control        Control  ##
##  PC182.CEL   PC182.CEL       PC182.CEL                Control        Control  ##
##  PC183.CEL   PC183.CEL       PC183.CEL                Control        Control  ##
##  PD061.CEL   PD061.CEL       PD061.CEL                Tratado      6_Semanas  ##
##  PD062.CEL   PD062.CEL       PD062.CEL                Tratado      6_Semanas  ##
##  PD063.CEL   PD063.CEL       PD063.CEL                Tratado      6_Semanas  ##
##                                                                               ##
###################################################################################

###################################################################################################
##                                                                                               ##
##  Conocer el número de muestra, el cual debe de coincidir con el número de muestras del SDRF.  ##
##                                                                                               ##
###################################################################################################
column_sample <- ncol(exprs(raw_data)) 
row_probes <- nrow(exprs(raw_data))

print(column_sample)
print(c("Número de muestras:",column_sample))
print(c("Número de sondas:",row_probes))

#####################################################################################################################
##                                                                                                                 ##
##   PC181.CEL PC182.CEL PC183.CEL PD061.CEL PD062.CEL PD063.CEL PD101.CEL PD102.CEL PD103.CEL PD141.CEL PD142.CEL ##
## 1      3108      3698      1976      4092      3112      4392      2415      3499      4518      3718      2506 ##
## 2       107       137        47        98        90        99        81        90       120       123        61 ##
## 3      2893      3630      1980      3925      3025      4605      2372      3128      4661      3255      2353 ##
## 4       101        62        57        83        96       103        87        54       101       137        52 ##
## 5       103       110        72        94        97       147        73        97       131        98        76 ##
## 6        49        63        49       103        54        48        67        68        65        85        61 ##
##   PD143.CEL PD181.CEL PD182.CEL PD183.CEL                                                                       ##
## 1      4129      2980      3031      2790                                                                       ##
## 2        78        60        65        50                                                                       ##
## 3      4158      2838      2854      2969                                                                       ##
## 4        85        57        50        54                                                                       ##
## 5       111        69        82        77                                                                       ##
## 6        68        84        94        59                                                                       ##
##                                                                                                                 ##
#####################################################################################################################

##############################################
##                                          ##
##  dim(exprs(raw_data))                    ##
##   renglones - 2598544     Columnas - 15  ##
##                                          ##
##############################################

#######################################
##                                   ##
##  Inicio del control de calidad.   ##
##                                   ##
#######################################

#######################################################
##                                                   ##
##  Biobase::exprs(raw_data)[1:5, 1:5] - (Opcional)  ##
##                                                   ##
#######################################################

#########################################################
##                                                     ##
##   PC181.CEL PC182.CEL PC183.CEL PD061.CEL PD062.CEL ##
## 1      3108      3698      1976      4092      3112 ##
## 2       107       137        47        98        90 ##
## 3      2893      3630      1980      3925      3025 ##
## 4       101        62        57        83        96 ##
## 5       103       110        72        94        97 ##
##                                                     ##
#########################################################
exp_raw <- log2(Biobase::exprs(raw_data)) 
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE) 

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], 
                    Phenotype = pData(raw_data)$Factor.Value.phenotype)

################################################################
##                                                            ##
##  Gráfico de Análisis de componentes principales.           ##
##  Agregar más colores en caso de tener más de 2 fenotipos.  ##
##                                                            ##
################################################################
jpeg("PCA_raw.jpg")
ggplot(dataGG, aes(PC1, PC2)) + geom_point(aes(colour = Phenotype)) + 
       ggtitle("PCA plot of the log-transformed raw expression data") + 
       xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) + 
       ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) + 
       theme(plot.title = element_text(hjust = 0.5)) + 
       coord_fixed(ratio = sd_ratio) + 
       scale_shape_manual(values = c(4,15)) + 
       scale_color_manual(values = c("darkorange2", "dodgerblue4")) + 
       theme_bw()
dev.off()  

########################
##                    ##
##  Gráfico de cajas  ##
##                    ##
########################
svg("Boxplot_raw.svg")
oligo::boxplot(raw_data, target = "core", 
               main = "Gráfico de cajas del log2 de las intensidades")
dev.off()

########################################
##                                    ##
##  Generar gráfico en 3D (Opcional)  ##
##                                    ##
########################################
###########################################################################################
##                                                                                       ##
##  dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2], PC3 = PCA_raw$x[,3],  ##
##                    Phenotype = pData(raw_data)$Factor.Value.phenotype,                ##
##                    Time = pData(raw_data)$Time.treatment)                             ##
##                                                                                       ##
##  p <- plot_ly(dataGG, x = ~PC1, y = ~PC2, z = ~PC3,                                   ##
##           color = ~Time, colors = c('#BF382A', '#0C4B8E','#D0F9B1'))                  ##
##                                                                                       ##
##  plot_ly(dataGG, x = ~PC1, y = ~PC2, z = ~PC3,                                        ##
##       color = ~Time, colors = c('#BF382A', '#0C4B8E','green'))                        ##
##                                                                                       ##
###########################################################################################

arrayQualityMetrics(expressionset = raw_data, 
                    outdir = "Reporte_Sin_Normalizar", 
                    force = TRUE , do.logtransform = TRUE, 
                    intgroup = c("Factor.Value.phenotype", "Time.treatment"))

##################################
##                              ##
##  Normalización de los datos  ##
##                              ##
##################################
carcinogen_eset <- oligo::rma(raw_data, target="core")  
#############################
##                         ##
##  Background correcting  ##
##  Normalizing            ##
##  Calculating Expression ##
##                         ##
#############################

############################
##                        ##
##  dim(carcinogen_eset)  ##
##  Features  Samples     ##
##    41345       15      ##
##                        ##
############################

exp_carcinogen_eset <- exprs(carcinogen_eset)

##########################################################
##                                                      ##
##  Indica número de renglones y columnas,              ##
##  llevar control del universo de probes que tenemos.  ##
##                                                      ##
##########################################################
dim(exp_carcinogen_eset) 
PCA_eset <- prcomp(t(exprs(carcinogen_eset)), scale = FALSE)

dataGG <- data.frame(PC1 = PCA_eset$x[,1], PC2 = PCA_eset$x[,2], 
                    Phenotype = pData(carcinogen_eset)$Factor.Value.phenotype)

################################################################
##                                                            ##
##  Grafico de PCA Normalizado                                ##
##  Agregar más colores en caso de tener más de 2 fenotipos.  ##
##                                                            ##
################################################################
jpeg("PCA_raw_Normalized.jpg")
ggplot(dataGG, aes(PC1, PC2)) + geom_point(aes(colour = Phenotype)) + 
        ggtitle("PCA plot of the log-transformed raw expression data") + 
        xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) + 
        ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        coord_fixed(ratio = sd_ratio) + 
        scale_shape_manual(values = c(4,15)) + 
        scale_color_manual(values = c("green", "dodgerblue4")) + 
        theme_bw()
dev.off()

########################################
##                                    ##
##  Generar gráfico en 3D (Opcional)  ##
##                                    ##
########################################
##############################################################################################
##                                                                                          ##
##  dataGG <- data.frame(PC1 = PCA_eset$x[,1], PC2 = PCA_eset$x[,2], PC3 = PCA_eset$x[,3],  ##
##                   Phenotype = pData(carcinogen_eset)$Factor.Value.phenotype,             ##
##                   Time = pData(carcinogen_eset)$Time.treatment)                          ##
##                                                                                          ##
##  p <- plot_ly(dataGG, x = ~PC1, y = ~PC2, z = ~PC3,                                      ##
##            color = ~Time, colors = c('#BF382A', '#0C4B8E','#D0F9B1'))                    ##
##                                                                                          ##
##  plot_ly(dataGG, x = ~PC1, y = ~PC2, z = ~PC3,                                           ##
##       color = ~Time, colors = c('#BF382A', '#0C4B8E','green'))                           ##
##                                                                                          ##
##############################################################################################

############################################
##                                        ##
##       Generando Mapa de Calor          ##
##  de las relaciones entre las muestras  ##
##                                        ##
############################################

dists <- as.matrix(dist(t(exp_carcinogen_eset), method = "manhattan")) 
colnames(dists) <- NULL
diag(dists) <- NA
rownames(dists) <- pData(carcinogen_eset)$Time.treatment
hmcol <- colorRampPalette(c("red","black","green"))(255)

png("Heatmap.png")
pheatmap(dists, col = rev(hmcol), clustering_distance_rows = "manhattan", 
        clustering_distance_cols = "manhattan")
dev.off()

############################
##                        ##
##  Generando Histograma  ##
##                        ##
############################

carcinogen_medians <- rowMedians(exprs(carcinogen_eset)) ## Medianas
svg("hist_normal.svg")
hist_res <- hist(carcinogen_medians, 100, freq=FALSE)
dev.off()

#########################
##                     ##
##  Filtrado de Genes  ##
##                     ##
#########################

emp_mu <- hist_res$breaks[which.max(hist_res$density)]
emp_sd <- mad(carcinogen_medians)/2
prop_cental <- 0.50 ##  Estandar de la distribucion
cut_val <- 0.05 / prop_cental ## el 5%  de la distribucion
thresh_median <- qnorm(0.05 / prop_cental, emp_mu, emp_sd)
no_of_samples <- table(paste0(pData(carcinogen_eset)$Factor.Value.disease., "_", 
                       pData(carcinogen_eset)$Factor.Value.phenotype.))
samples_cutoff <- min(no_of_samples)
idx_thresh_median <- apply(exprs(carcinogen_eset), 1, function(x){
sum(x > thresh_median) >= samples_cutoff } )

table(idx_thresh_median)
#########################
##                     ##
##  idx_thresh_median  ##
##      TRUE           ##
##      41345          ##
##                     ##
#########################

carcinogen_filtered <- subset(carcinogen_eset, idx_thresh_median)

carcinogen_filtered_anno <- AnnotationDbi::select(mogene20sttranscriptcluster.db, 
    keys=(featureNames(carcinogen_filtered)), 
    columns = c("SYMBOL", "GENENAME"), 
    keytype="PROBEID")
###################################################################
##                                                               ##
##  'select()' returned 1:many mapping between keys and columns  ##
##                                                               ##
###################################################################

carcinogen_filtered_anno <- subset(carcinogen_filtered_anno, !is.na(SYMBOL))
carcinogen_filtered_anno_grouped <- group_by(carcinogen_filtered_anno, PROBEID)
head(carcinogen_filtered_anno_grouped )
##########################################################################
##                                                                      ##
##  A tibble: 6 x 3                                                     ##
##  Groups:   PROBEID [6]                                               ##
##    PROBEID  SYMBOL  GENENAME                                         ##
##    <chr>    <chr>   <chr>                                            ##
##  1 17210855 Lypla1  lysophospholipase 1                              ##
##  2 17210869 Tcea1   transcription elongation factor A (SII) 1        ##
##  3 17210883 Gm16041 predicted gene 16041                             ##
##  4 17210887 Atp6v1h ATPase, H+ transporting, lysosomal V1 subunit H  ##
##  5 17210904 Oprk1   opioid receptor, kappa 1                         ##
##  6 17210912 Rb1cc1  RB1-inducible coiled-coil 1                      ##
##                                                                      ##
##########################################################################
dim(carcinogen_filtered_anno_grouped)  
#######################
##                   ##
##  [1] 36392     3  ##
##                   ##
#######################

anno_summarized <- dplyr::summarize(carcinogen_filtered_anno_grouped, 
                                    no_of_matches = n_distinct(SYMBOL))
anno_filtered <- filter(anno_summarized, no_of_matches > 1)

ids_to_exlude <- ((featureNames(carcinogen_filtered) %in% anno_filtered$PROBEID) | 
                   featureNames(carcinogen_filtered) %in% subset(carcinogen_filtered_anno, 
                   is.na(SYMBOL))$PROBEID)

carcinogen_final <- subset(carcinogen_filtered, !ids_to_exlude)
############################
##                        ##
## dim(carcinogen_final)  ##
## Features  Samples      ##
##   39724       15       ##
##                        ##
############################
fData(carcinogen_final)$PROBEID <- rownames(fData(carcinogen_final))
fData(carcinogen_final) <- left_join(fData(carcinogen_final), carcinogen_filtered_anno)
#############################
##                         ##
##  Joining, by = PROBEID  ##
##                         ##
#############################

rownames(fData(carcinogen_final)) <-fData(carcinogen_final)$PROBEID

#########################################
##                                     ##
##  Análisis de expresión diferencial  ##
##                                     ##
#########################################
carcinogen_exprs_final = exprs(carcinogen_final)
head(carcinogen_exprs_final)
fac_int <- pData(carcinogen_final)$Factor.Value.phenotype
design <- model.matrix(~0 + fac_int) ## Solamente funciona para un fenotipo binario, si se tienen más fenotipos, cambiar los contrastes
print(design)
########################################
##                                    ##
##     fac_intControl fac_intTratado  ##
##  1               1              0  ##
##  2               1              0  ##
##  3               1              0  ##
##  4               0              1  ##
##  5               0              1  ##
##  6               0              1  ##
##  7               0              1  ##
##  8               0              1  ##
##  9               0              1  ##
##  10              0              1  ##
##  11              0              1  ##
##  12              0              1  ##
##  13              0              1  ##
##  14              0              1  ##
##  15              0              1  ##
##  attr(,"assign")                   ##
##  [1] 1 1                           ##
##  attr(,"contrasts")                ##
##  attr(,"contrasts")$fac_int        ##
##  [1] "contr.treatment"             ##
##                                    ##
########################################

colnames(design) <- unique(factor(fac_int))
#####################################################
##                                                 ##
##  Se definen los grupos de acuerdo al fenotipo.  ##
##                                                 ##
#####################################################
print(design) 
##################################
##                              ##
##    Control Tratado           ##
##  1        1       0          ##
##  2        1       0          ##
##  3        1       0          ##
##  4        0       1          ##
##  5        0       1          ##
##  6        0       1          ##
##  7        0       1          ##
##  8        0       1          ##
##  9        0       1          ##
##  10       0       1          ##
##  11       0       1          ##
##  12       0       1          ##
##  13       0       1          ##
##  14       0       1          ##
##  15       0       1          ##
##  attr(,"assign")             ##
##  [1] 1 1                     ##
##  attr(,"contrasts")          ##
##  attr(,"contrasts")$fac_int  ##
##  [1] "contr.treatment"       ##
##                              ##
##################################

lmfit <- lmFit(carcinogen_final, design) ##Establecemos modelos lineales, correlaciones lineales entre los contrastes
################################################################
##                                                            ##
##  Cambiar aqui los contrastes, si se tienen mas fenotipos.  ##
##                                                            ##
################################################################
cont_tratado <- makeContrasts(Control-Tratado, levels = design)
print(cont_tratado)
##################################
##                              ##
##          Contrasts           ##
## Levels    Control - Tratado  ##
##   Control                 1  ##
##   Tratado                -1  ##
##                              ##
##################################

lmfit.cont <- contrasts.fit(lmfit, cont_tratado)
lmfit.cont.ebayes <- eBayes(lmfit.cont)

print(topTable(lmfit.cont.ebayes)) 
################################################################################################################
##                                                                                                            ##
##            PROBEID    SYMBOL                                                           GENENAME     logFC  ##
##  17335467 17335467    Cdkn1a                         cyclin-dependent kinase inhibitor 1A (P21) -4.471969  ##
##  17543396 17543396     Eda2r                                          ectodysplasin A2 receptor -2.911571  ##
##  17529764 17529764      Pls1                                              plastin 1 (I-isoform) -3.161245  ##
##  17261865 17261865     Ccng1                                                          cyclin G1 -1.936880  ##
##  17484068 17484068      Lhpp phospholysine phosphohistidine inorganic pyrophosphate phosphatase  2.801670  ##
##  17534385 17534385     Gria3                    glutamate receptor, ionotropic, AMPA3 (alpha 3) -2.981105  ##
##  17289794 17289794      Plk2                                                 polo like kinase 2 -2.044308  ##
##  17462437 17462437     Usp18                                    ubiquitin specific peptidase 18 -1.835164  ##
##  17301697 17301697 Tnfrsf10b             tumor necrosis factor receptor superfamily, member 10b -3.165906  ##
##  17394538 17394538     Sulf2                                                        sulfatase 2 -1.933202  ##
##             AveExpr         t      P.Value    adj.P.Val        B                                           ##
##  17335467  9.736155 -30.75620 1.229233e-15 4.881530e-11 22.83271                                           ##
##  17543396  5.782974 -19.57152 1.394459e-12 2.768838e-08 18.07583                                           ##
##  17529764  6.825638 -18.53864 3.201208e-12 4.237545e-08 17.42081                                           ##
##  17261865 10.216770 -17.92845 5.336955e-12 5.298529e-08 17.00953                                           ##
##  17484068  6.565091  17.09573 1.100483e-11 7.467486e-08 16.41676                                           ##
##  17534385  6.806257 -17.06770 1.128246e-11 7.467486e-08 16.39614                                           ##
##  17289794  7.906152 -16.63134 1.670781e-11 8.986179e-08 16.06939                                           ##
##  17462437  7.291381 -16.45519 1.962887e-11 8.986179e-08 15.93433                                           ##
##  17301697  5.808938 -16.38732 2.089462e-11 8.986179e-08 15.88180                                           ##
##  17394538  7.426638 -16.30111 2.262837e-11 8.986179e-08 15.81467                                           ##
##                                                                                                            ##
################################################################################################################
#################################
## ID --- Simbolo --- Nombre ####
#################################

lmfit.cont.ebayes.table <- topTable(lmfit.cont.ebayes, number = Inf)
lmfit.cont.ebayes.table2 = subset(lmfit.cont.ebayes.table,
                            lmfit.cont.ebayes.table$GENENAME != "NA")

###################################################################
##                                                               ##
##  Genera el resumen de la expresion diferencial de los genes,  ##
##  ordenado del mas significativo al menos significativo        ##
##                                                               ##
###################################################################
Diff_expresed_genes = topTable(lmfit.cont.ebayes, num = Inf) 
print(head(Diff_expresed_genes))
print(str(Diff_expresed_genes))
print(str(topTable(lmfit.cont.ebayes, num = Inf)))

print(str(Diff_expresed_genes))
######################################################################
##                                                                  ##
##  Filtramos los mas significativos, se puede modificar el logFC,  ##
##  Ejemplo: 3, 4 y el adj.P.value (5% de significancia)            ##
##                                                                  ##
######################################################################
Log_FC_Diff_expresed_genes = subset(Diff_expresed_genes, abs(logFC) > 2 & adj.P.Val < 0.05)

######################################################################
##                                                                  ##
##  Hacemos un dataframe de los valores de expresion del carcinoma  ##
##                                                                  ##
######################################################################
carcinogen_exprs_final_df = as.data.frame(carcinogen_exprs_final)

print(str(carcinogen_exprs_final_df))
##########################################################################
##                                                                      ##
##  Filtramos los datos de expresión con los genes más significativos,  ##
##  basado en el código de la sonda.                                    ##
##                                                                      ##
##########################################################################
carcinogen_exprs_final_dif_expres = carcinogen_exprs_final_df[
        rownames(carcinogen_exprs_final_df) %in% Log_FC_Diff_expresed_genes$PROBEID,
        ]

##############################################################
##                                                          ##
##  Reducimos el objeto a solo código de nombres y sondas.  ##
##                                                          ##
##############################################################
Log_FC_Diff_expresed_genes_names = Log_FC_Diff_expresed_genes[, c(1:2)]

############################################################################
##                                                                        ##
##  Generando el mapa de calor de todos los genes, con el FC de interes.  ##
##                                                                        ##
############################################################################
carcinogen_exprs_final_dif_expres_ord = subset(carcinogen_exprs_final_dif_expres)
carcinogen_exprs_final_dif_expres_ord = carcinogen_exprs_final_dif_expres_ord[
        rownames(carcinogen_exprs_final_dif_expres_ord) %in% Log_FC_Diff_expresed_genes$PROBEID,
        ]
carcinogen_exprs_final_dif_expres_ord$PROBEID = rownames(carcinogen_exprs_final_dif_expres_ord)


Log_FC_Diff_expresed_genes_merge_ord = merge(carcinogen_exprs_final_dif_expres_ord, Log_FC_Diff_expresed_genes_names)
Log_FC_Diff_expresed_genes_merge_final_ord = subset( Log_FC_Diff_expresed_genes_merge_ord, 
                                                     Log_FC_Diff_expresed_genes_merge_ord$SYMBOL != "NA")
Log_FC_Diff_expresed_genes_merge_final2_ord = Log_FC_Diff_expresed_genes_merge_final_ord[
        !duplicated(Log_FC_Diff_expresed_genes_merge_final_ord$SYMBOL),]

rownames(Log_FC_Diff_expresed_genes_merge_final2_ord) = Log_FC_Diff_expresed_genes_merge_final2_ord$SYMBOL

Log_FC_Diff_expresed_genes_merge_final2_ord$PROBEID = NULL 
Log_FC_Diff_expresed_genes_merge_final2_ord$SYMBOL = NULL 

diff_matrix2_ord = as.matrix(Log_FC_Diff_expresed_genes_merge_final2_ord)

png("HeatmapFC_AllGenes.png")
heatmap.2(diff_matrix2_ord, trace= "none", col = colorRampPalette(c("green","black","red"))(255),
         margins = c(5,10), lwid = c(5,15), lhei = c(3,15), density.info = "none" )
dev.off()

#############################################################################
##                                                                         ##
##  Calculamos la distribución de los cuartiles en el primer control,      ## 
##  Nota: modificar el objeto al nombre de sus controles                   ##
##  Modificar los valores de corte dependiendo de los cuartiles de salida  ##
##                                                                         ##
#############################################################################
print(summary(carcinogen_exprs_final_dif_expres$PC181.CEL))
#############################################################################
##                                                                         ##
##  Generando el mapa de calor del primer cuartil                          ##
##  Modificar los valores de corte dependiendo de los cuartiles de salida  ##
##                                                                         ##
#############################################################################
carcinogen_exprs_final_dif_expres_ord_1st = subset(carcinogen_exprs_final_dif_expres,PC181.CEL <= 3.828)
carcinogen_exprs_final_dif_expres_ord_1st2 = carcinogen_exprs_final_dif_expres_ord_1st[
        rownames(carcinogen_exprs_final_dif_expres_ord_1st) %in% Log_FC_Diff_expresed_genes$PROBEID,
        ]
carcinogen_exprs_final_dif_expres_ord_1st2$PROBEID = rownames(carcinogen_exprs_final_dif_expres_ord_1st2)


Log_FC_Diff_expresed_genes_merge_ord_1st = merge(carcinogen_exprs_final_dif_expres_ord_1st2, Log_FC_Diff_expresed_genes_names)
Log_FC_Diff_expresed_genes_merge_final_ord_1st = subset( Log_FC_Diff_expresed_genes_merge_ord_1st, 
                                                         Log_FC_Diff_expresed_genes_merge_ord_1st$SYMBOL != "NA")
Log_FC_Diff_expresed_genes_merge_final2_ord_1st = Log_FC_Diff_expresed_genes_merge_final_ord_1st[
        !duplicated(Log_FC_Diff_expresed_genes_merge_final_ord_1st$SYMBOL),]

rownames(Log_FC_Diff_expresed_genes_merge_final2_ord_1st) = Log_FC_Diff_expresed_genes_merge_final2_ord_1st$SYMBOL

Log_FC_Diff_expresed_genes_merge_final2_ord_1st$PROBEID = NULL 
Log_FC_Diff_expresed_genes_merge_final2_ord_1st$SYMBOL = NULL 

diff_matrix2_ord_1st = as.matrix(Log_FC_Diff_expresed_genes_merge_final2_ord_1st)

png("Heatmap_1st_cuartil.png")
heatmap.2(diff_matrix2_ord_1st, trace= "none", col = colorRampPalette(c("green","black","red"))(255),
         margins = c(5,10), lwid = c(5,15), lhei = c(3,15), density.info = "none"
       , Colv=FALSE )
dev.off()
##############################################################################
##                                                                          ##
##  Generando el mapa de calor del segundo  cuartil                         ##
##  Modificar los valores de corte dependiendo de los cuartiles de salida.  ##
##                                                                          ##
##############################################################################

carcinogen_exprs_final_dif_expres_ord_2st = subset(carcinogen_exprs_final_dif_expres,PC181.CEL > 3.828 & PC181.CEL < 6.423)
carcinogen_exprs_final_dif_expres_ord_2st2 = carcinogen_exprs_final_dif_expres_ord_2st[
        rownames(carcinogen_exprs_final_dif_expres_ord_2st) %in% Log_FC_Diff_expresed_genes$PROBEID,
        ]
carcinogen_exprs_final_dif_expres_ord_2st2$PROBEID = rownames(carcinogen_exprs_final_dif_expres_ord_2st2)


Log_FC_Diff_expresed_genes_merge_ord_2st = merge(carcinogen_exprs_final_dif_expres_ord_2st2, Log_FC_Diff_expresed_genes_names)
Log_FC_Diff_expresed_genes_merge_final_ord_2st = subset( Log_FC_Diff_expresed_genes_merge_ord_2st, 
                                                         Log_FC_Diff_expresed_genes_merge_ord_2st$SYMBOL != "NA")
Log_FC_Diff_expresed_genes_merge_final2_ord_2st = Log_FC_Diff_expresed_genes_merge_final_ord_2st[
        !duplicated(Log_FC_Diff_expresed_genes_merge_final_ord_2st$SYMBOL),]

rownames(Log_FC_Diff_expresed_genes_merge_final2_ord_2st) = Log_FC_Diff_expresed_genes_merge_final2_ord_2st$SYMBOL

Log_FC_Diff_expresed_genes_merge_final2_ord_2st$PROBEID = NULL 
Log_FC_Diff_expresed_genes_merge_final2_ord_2st$SYMBOL = NULL 


diff_matrix2_ord_2st = as.matrix(Log_FC_Diff_expresed_genes_merge_final2_ord_2st)

png("Heatmap_2nd_cuartil.png")
heatmap.2(diff_matrix2_ord_2st, trace= "none", col = colorRampPalette(c("green","black","red"))(255),
         margins = c(5,10), lwid = c(5,15), lhei = c(3,15), density.info = "none"
       , Colv=FALSE )
dev.off()
##############################################################################
##                                                                          ##
##  Generando el mapa de calor del tercer cuartil                           ##
##  Modificar los valores de corte dependiendo de los cuartiles de salida.  ##
##                                                                          ##
##############################################################################

carcinogen_exprs_final_dif_expres_ord_3st = subset(carcinogen_exprs_final_dif_expres, PC181.CEL >= 6.423)
carcinogen_exprs_final_dif_expres_ord_3st2 = carcinogen_exprs_final_dif_expres_ord_3st[
        rownames(carcinogen_exprs_final_dif_expres_ord_3st) %in% Log_FC_Diff_expresed_genes$PROBEID,
        ]
carcinogen_exprs_final_dif_expres_ord_3st2$PROBEID = rownames(carcinogen_exprs_final_dif_expres_ord_3st2)


Log_FC_Diff_expresed_genes_merge_ord_3st = merge(carcinogen_exprs_final_dif_expres_ord_3st2, Log_FC_Diff_expresed_genes_names)
Log_FC_Diff_expresed_genes_merge_final_ord_3st = subset( Log_FC_Diff_expresed_genes_merge_ord_3st, 
                                                         Log_FC_Diff_expresed_genes_merge_ord_3st$SYMBOL != "NA")
Log_FC_Diff_expresed_genes_merge_final2_ord_3st = Log_FC_Diff_expresed_genes_merge_final_ord_3st[
        !duplicated(Log_FC_Diff_expresed_genes_merge_final_ord_3st$SYMBOL),]

rownames(Log_FC_Diff_expresed_genes_merge_final2_ord_3st) = Log_FC_Diff_expresed_genes_merge_final2_ord_3st$SYMBOL

Log_FC_Diff_expresed_genes_merge_final2_ord_3st$PROBEID = NULL 
Log_FC_Diff_expresed_genes_merge_final2_ord_3st$SYMBOL = NULL 

diff_matrix2_ord_3st = as.matrix(Log_FC_Diff_expresed_genes_merge_final2_ord_3st)


png("Heatmap_3rd_cuartil.png")
heatmap.2(diff_matrix2_ord_3st, trace= "none", col = colorRampPalette(c("green","black","red"))(255),
         margins = c(5,10), lwid = c(5,15), lhei = c(3,15), density.info = "none"
        , Colv=FALSE )
dev.off()
#########################
##                     ##
##  Gráfico de Volcan  ##
##                     ##
#########################
# with(subset(lmfit.cont.ebayes.table2, adj.P.Val<.05 & abs(logFC)>2), 
#      points(logFC, -log10(P.Value), pch=20, col="green"))

png("volcano_plot.png")
with(lmfit.cont.ebayes.table2, plot(logFC, -log10(P.Value), 
    pch=20, main="Volcano plot", xlim=c(-7.5,6), ylim=c(0,12)))
with(subset(lmfit.cont.ebayes.table2, adj.P.Val<.05 & abs(logFC)>2), 
     points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(lmfit.cont.ebayes.table2, adj.P.Val<.05 & (logFC)< -2 ), 
     points(logFC, -log10(P.Value), pch=20, col="green"), pCutoff = 10e-32)
abline(v = 2, col ="red")
abline(v = -2, col ="darkgreen")
abline(h = 2, col ="black")

with(subset(lmfit.cont.ebayes.table2, adj.P.Val<.05 & abs(logFC)>3), 
    textxy(logFC, -log10(P.Value), labs=SYMBOL, cex=0.3, offset=0.5))

dev.off()

############# Dynamic Volcano Plot with 3 and 4
with(lmfit.cont.ebayes.table2, plot(logFC, -log10(P.Value), 
    pch=20, main="Volcano plot", xlim=c(-7.5,6), ylim=c(0,12)))
with(subset(lmfit.cont.ebayes.table2, adj.P.Val<.05 & abs(logFC)>4), 
     points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(lmfit.cont.ebayes.table2, adj.P.Val<.05 & (logFC)< -4 ), 
     points(logFC, -log10(P.Value), pch=20, col="green"), pCutoff = 10e-32)
abline(v = 4, col ="red")
abline(v = -4, col ="darkgreen")
abline(h = 2, col ="black")


results = decideTests(lmfit.cont.ebayes)
####################################
##                                ##
##  TestResults matrix            ##
##            Contrasts           ##
##             Control - Tratado  ##
##    17200001                 0  ##
##    17200003                 0  ##
##    17200005                 0  ##
##    17200007                 0  ##
##    17200009                 0  ##
##  39707 more rows ...           ##
##                                ##
####################################


DE_genes_php_pc <- subset(lmfit.cont.ebayes.table, adj.P.Val < 0.05)$PROBEID

back_genes_idx <- genefilter::genefinder(carcinogen_final, 
                  as.character(DE_genes_php_pc), method = "manhattan", scale = "none")

back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)
back_genes <- featureNames(carcinogen_final)[back_genes_idx]
intersect(back_genes, DE_genes_php_pc)

gene_IDs <- rownames(lmfit.cont.ebayes.table)
in_universe <- gene_IDs %in% c(DE_genes_php_pc, back_genes)
in_selection <- gene_IDs %in% DE_genes_php_pc
all_genes <- in_selection[in_universe]
all_genes <- factor(as.integer(in_selection[in_universe]) == 1)
names(all_genes) <- gene_IDs[in_universe]
differentially_espresed_IDs = all_genes[all_genes == TRUE]
differentially_expresss_genes_expression = exprs(carcinogen_final)[
    rownames(exprs(carcinogen_final)) %in% c(names(differentially_espresed_IDs)),]

###########################################
##                                       ##
##  Se generan los diferentes archivos.  ##
##                                       ##
###########################################

write.table(file="DE_express_control_tratados.txt", 
            differentially_expresss_genes_expression, 
            quote=F, sep="\t")

write.table(file="TopTable_genes.txt", 
            lmfit.cont.ebayes.table2, 
            quote=F, sep="\t")
            
write.table(file="TopTable_genes_underexpressed.txt",
            subset(lmfit.cont.ebayes.table2$SYMBOL, 
            lmfit.cont.ebayes.table2$adj.P.Val < 0.05 & 
            lmfit.cont.ebayes.table2$logFC <= 2), 
            quote=F)
write.table(file="TopTable_genes_overexpressed.txt",
            subset(lmfit.cont.ebayes.table2$SYMBOL, 
            lmfit.cont.ebayes.table2$adj.P.Val < 0.05 & 
            lmfit.cont.ebayes.table2$logFC >= 2), 
            quote=F)
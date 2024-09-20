setwd("G:/Mi unidad/1- ESPECIALIZACION/6 - PROCESAMIENTO_DATOS/ProyectoPMD/")
library(readr)
frecuancias_total <- read_csv("G:/Mi unidad/1- ESPECIALIZACION/6 - PROCESAMIENTO_DATOS/MicrobialData (5)/asv_table.csv")
id_LA <- read_csv("G:/Mi unidad/1- ESPECIALIZACION/6 - PROCESAMIENTO_DATOS/ProyectoPMD/id_LA.txt")
id_Chinactrl <- read_csv("G:/Mi unidad/1- ESPECIALIZACION/6 - PROCESAMIENTO_DATOS/ProyectoPMD/id_Chinactrl.txt")

taxonomia_total <- read_csv("G:/Mi unidad/1- ESPECIALIZACION/6 - PROCESAMIENTO_DATOS/MicrobialData (5)/taxonomy_table.csv")

extracted <- frecuancias_total[, c(id_LA$id, id_Chinactrl$id)]
extracted <- cbind(frecuancias_total$Species, extracted)
colnames(extracted)[colnames(extracted) == "frecuancias_total$Species"] <- "Species"
write.csv(extracted, "frec_LAyCHINA.csv", row.names = F, quote = F)



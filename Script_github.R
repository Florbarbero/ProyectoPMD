# Paso 1: con la otu_table, taxomony_table y la metadata creo el archivo phyloseq

#https://joey711.github.io/phyloseq/import-data
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
library(phyloseq)

#Importar a R otumat y taxmat

# Convertir el tibble a un data frame
otumat <- as.data.frame(otumat)

# Asignar los valores de la primera columna como nombres de fila
rownames(otumat) <- otumat[, 1]

# Convertir el tibble a un data frame
taxmat <- as.data.frame(taxmat)

# Asignar los nombres de fila
rownames(taxmat) <- rownames(otumat)

#eliminar la primera fila de ambas tablas
otumat <- otumat[, -1]
taxmat <- taxmat[, -1]


# Convertir otumat y taxmat de data frames a matrices
otumat <- as.matrix(otumat)
taxmat <- as.matrix(taxmat)

class(otumat)
class(taxmat)

# Crea una tabla de OTU (abundancia de especies)
OTU <- otu_table(otumat, taxa_are_rows = TRUE)

# Crea una tabla de taxonomía
TAX <- tax_table(taxmat)

physeq = phyloseq(OTU, TAX)
physeq

plot_bar(physeq, fill = "Phylum")

# Crear el objeto de metadatos sampledata
sample_data <- data.frame(
  País = c(
    # China samples
    rep("China", 33),
    
    # EE.UU samples
    rep("EE.UU", 109)
  )
)

# Asignar los nombres de fila
row_names <- c(
  "SRR7142468", "SRR7142469", "SRR7142470", "SRR7142478", "SRR7142479",
  "SRR7142486", "SRR7142487", "SRR7142490", "SRR7142491", "SRR7142492",
  "SRR7142493", "SRR7142494", "SRR7142495", "SRR7142503", "SRR7142504",
  "SRR7142508", "SRR7142511", "SRR7142512", "SRR7142513", "SRR7142518",
  "SRR7142524", "SRR7142471", "SRR7142472", "SRR7142473", "SRR7142474",
  "SRR7142475", "SRR7142489", "SRR7142496", "SRR7142501", "SRR7142509",
  "SRR7142510", "SRR7142514", "SRR7142515", "SRR8100127", "SRR8100131",
  "SRR8100143", "SRR8100147", "SRR8100150", "SRR8100152", "SRR8100156",
  "SRR8100160", "SRR8100161", "SRR8100165", "SRR8100169", "SRR8100173",
  "SRR8100177", "SRR8100181", "SRR8100185", "SRR8100187", "SRR8100194",
  "SRR8100198", "SRR8100204", "SRR8100208", "SRR8100214", "SRR8100216",
  "SRR8100218", "SRR8100223", "SRR8100228", "SRR8100229", "SRR8100233",
  "SRR8100243", "SRR8100247", "SRR8100249", "SRR8100253", "SRR8100257",
  "SRR8100262", "SRR8100266", "SRR8100270", "SRR8100274", "SRR8100277",
  "SRR8100280", "SRR8100284", "SRR8100289", "SRR8100293", "SRR8100297",
  "SRR8100303", "SRR8100307", "SRR8100311", "SRR8100315", "SRR8100320",
  "SRR8100324", "SRR8100326", "SRR8100328", "SRR8100332", "SRR8100336",
  "SRR8100340", "SRR8100350", "SRR8100352", "SRR8100356", "SRR8100363",
  "SRR8100367", "SRR8100369", "SRR8100371", "SRR8100373", "SRR8100377",
  "SRR8100381", "SRR8100383", "SRR8100387", "SRR8100394", "SRR8100398",
  "SRR8100400", "SRR8100403", "SRR8100407", "SRR8100412", "SRR8100415",
  "SRR8100419", "SRR8100421", "SRR8100423", "SRR8100427", "SRR8100434",
  "SRR8100438", "SRR8100442", "SRR8100446", "SRR8100450", "SRR8100454",
  "SRR8100456", "SRR8100458", "SRR8100462", "SRR8100465", "SRR8100469",
  "SRR8100473", "SRR8100477", "SRR8100481", "SRR8100483", "SRR8100486",
  "SRR8100490", "SRR8100498", "SRR8100502", "SRR8100508", "SRR8100512",
  "SRR8100517", "SRR8100520", "SRR8100521", "SRR8100524", "SRR8100529",
  "SRR8100531", "SRR8100535", "SRR8100546", "SRR8100550", "SRR8100554",
  "SRR8100555", "SRR8100558"
)



# Set the row names for the dataframe
rownames(sample_data) <- row_names

SAMPLEDATA <- sample_data(sample_data)

physeq1 = merge_phyloseq(physeq, SAMPLEDATA)
physeq1

# Paso 2: Cálculo de la diversidad alfa (Chao1 y Shannon)

# Instala el paquete si no lo tienes
install.packages("dplyr")

# Carga el paquete
library(dplyr)

diversity_indices <- estimate_richness(physeq1, measures = c("Chao1", "Shannon"))

# Mostrar los resultados
print(diversity_indices)

# Extraer el sample_data de physeq1
sample_data_df <- as.data.frame(sample_data(physeq1))

# Añadir la información de 'Lugar' a la tabla de índices de diversidad
diversity_indices$País <- sample_data_df$País

# Agrupar por 'Lugar' y calcular las medias de los índices de diversidad
diversity_summary <- diversity_indices %>%
  group_by(País) %>%
  summarise(Chao1_mean = mean(Chao1, na.rm = TRUE),
            Shannon_mean = mean(Shannon, na.rm = TRUE))

# Mostrar los resultados
print(diversity_summary)

# Paso 3: Calculo de Phylum Únicos por País
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(microbiomeMarker)

# Realizar el análisis LEfSe
lefse_results <- run_lefse(
  physeq1,
  group = "País", 
  lda_cutoff = 2.0, # Puedes ajustar este valor según tus necesidades
  norm = "CPM" # Normalización a usar (CPM, TSS, RLE, etc.)
)

# Generar el gráfico de barras de los efectos de los indicadores
plot_ef_bar(lefse_results) +
  ggtitle("Indicadores LEfSe por Tipo de Suelo") +
  theme_minimal()

# Extraer el slot marker_table del objeto lefse_results
marker_table <- lefse_results@marker_table@.Data

# Convertir a data.frame
marker_df <- data.frame(
  feature = marker_table[[1]],
  enriched_group = marker_table[[2]],
  LDA_score = marker_table[[3]],
  pvalue = marker_table[[4]],
  qvalue = marker_table[[5]]
)

# Extraer el segundo nivel del taxón (Phylum)
marker_df <- marker_df %>%
  mutate(
    Phylum = sapply(strsplit(feature, "\\|"), function(x) x[2]) # Extrae el segundo nivel, que debería ser el Phylum
  )

# Resumir por Phylum y grupo enriquecido
phylum_summary <- marker_df %>%
  group_by(Phylum, enriched_group) %>%
  summarise(LDA_score = mean(LDA_score), .groups = 'drop')

ggplot(phylum_summary, aes(x = reorder(Phylum, LDA_score), y = LDA_score, fill = enriched_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Indicadores LEfSe por País",
    x = "Phylum",
    y = "Puntuación LDA"
  ) +
  theme(legend.position = "bottom")

library(dplyr)
library(tidyr)
library(ggplot2)

# Filtrar para taxones exclusivos por País
exclusive_taxa <- marker_df %>%
  group_by(Phylum) %>%
  filter(n_distinct(enriched_group) == 1) %>%
  ungroup()

# Resumir por Phylum y grupo enriquecido (esto ya es redundante pero se incluye por claridad)
phylum_summary_exclusive <- exclusive_taxa %>%
  group_by(Phylum, enriched_group) %>%
  summarise(LDA_score = mean(LDA_score), .groups = 'drop')

# Crear el gráfico con los taxones exclusivos
ggplot(phylum_summary_exclusive, aes(x = reorder(Phylum, LDA_score), y = LDA_score, fill = enriched_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Phylum Exclusivos por País",
    x = "Phylum",
    y = "Puntuación LDA"
  ) +
  theme(legend.position = "bottom")





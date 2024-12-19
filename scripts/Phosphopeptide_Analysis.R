
# ============================================================================== Debut ===================================================================================================================#
# Fonction pour vérifier et installer les packages manquants
install_if_missing <- function(packages) {
  # Identifier les packages manquants
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  # Installer les packages manquants
  if (length(missing_packages) > 0) {
    install.packages(missing_packages, dependencies = TRUE)
  }
  # Charger tous les packages
  lapply(packages, library, character.only = TRUE)
}

# Liste des packages nécessaires
required_packages <- c(
  "ggplot2", 
  "ggrepel", 
  "dplyr",
  "readxl", 
  "VennDiagram", 
  "tidyverse", 
  "plotly",
  "ComplexHeatmap", 
  "circlize", 
  "magrittr", 
  "ggseqlogo", 
  "rstudioapi",
  "ontologyIndex"
)

# Vérification et installation des packages
install_if_missing(required_packages)

# ============================================================================== Functions ===================================================================================================================#

source("/Users/lamine/INEM/Projets/Mario/RSK2_analysis/scripts/Util_Phospo.R")
# Les fonction presentes dans le fichier Util.R :

#===============================================================================  Select_Diff()
# La fonction Select_Diff identifie les protéines ou peptides différentiellement exprimés
# en fonction des données de log2 fold-change (log2FC) et de p-value. Elle retourne un tableau 
# contenant les informations fusionnées avec les métadonnées, et classe les entrées en trois 
# catégories : UP (surexprimées), DOWN (sous-exprimées) ou NO (non significatives) selon des 
# seuils définis par l'utilisateur (cutoff.fc et cutoff.pv).
# 
# Paramètres :
# - log2FC : Dataframe contenant les valeurs log2 fold-change.
# - pvalue : Dataframe contenant les p-values correspondantes.
# - metadata : Tableau contenant les métadonnées des protéines/peptides.
# - condition : Nom de la condition d'intérêt (par défaut : "Crush+shControl : Control").
# - cutoff.fc : Seuil pour log2 fold-change (par défaut : 1.8).
# - cutoff.pv : Seuil pour les p-values (par défaut : 0.05).
# 
# Retourne :
# - Un dataframe annoté avec les colonnes fusionnées et une classification des directions (UP, DOWN, NO).

#===============================================================================  Plot_Volcano()
# La fonction Plot_Volcano génère un graphique de type volcan interactif pour visualiser les 
# protéines ou peptides différentiellement exprimés en fonction des valeurs de log2 fold-change 
# (log2FC) et de p-value.
#
# Paramètres :
# - data : Dataframe contenant les colonnes log2FC, pvalue, et Direction (UP, DOWN, NO).
# - cutoff.fc : Seuil pour le log2 fold-change à visualiser.
# - cutoff.pv : Seuil pour les p-values à visualiser.
# - title : Titre du graphique (par défaut : "Title").
#
# Retourne :
# - Un graphique interactif de type volcan (format plotly) avec des points colorés selon la direction 
#   (UP, DOWN, NO), des lignes de seuil pour log2FC et p-value, et des étiquettes pour les protéines/peptides.

#===============================================================================  plot_abundance_protein_site()
# La fonction plot_abundance_protein_site permet de visualiser l'abondance moyenne par site 
# pour une protéine donnée à travers différentes conditions expérimentales.
#
# Paramètres :
# - protein_name : Nom de la protéine d'intérêt (chaîne de caractères).
# - abundance_data : Dataframe contenant les données d'abondance pour les conditions et réplicats.
# - metadata_phospo_cleaned : Dataframe contenant les métadonnées associées, incluant les colonnes 'Gene' et 'Site'.
#
# Sortie :
# - Un graphique de type barplot, où les barres représentent l'abondance moyenne par site 
#   sous différentes conditions, colorées par condition.

#===============================================================================  Make_Rpos_Consensus()
# La fonction Make_Rpos_Consensus permet d'analyser les peptides phosphorylés et de générer un consensus
# basé sur les résidus "R" en positions -3 et -5 par rapport au site de phosphorylation (*).
#
# Paramètres :
# - data : Dataframe contenant les informations sur les peptides et protéines.
# - Protein_Name : Nom de la colonne contenant les noms des protéines (par défaut "Protein_Name").
# - Peptide_Localized : Nom de la colonne contenant la séquence localisée des peptides (par défaut "Peptide_Localized").
#
# Étapes principales :
# 1. Comptage des occurrences de peptides par protéine.
# 2. Attribution de catégories ("direct", "two", "three", "plus") selon le nombre de peptides.
# 3. Identification des résidus "R" aux positions -3 et -5 par rapport au site de phosphorylation.
# 4. Tri des résultats en mettant les substrats "direct" en priorité.
#
# Sortie :
# - Un dataframe enrichi avec des annotations des positions "R" et des catégories de substrats.

#===============================================================================  Plot_logo()
# Plot_logo
# Cette fonction génère un logo de séquence basé sur une liste de peptides alignés, permettant de visualiser 
# la conservation et les caractéristiques des acides aminés à différentes positions.
#
# Paramètres :
# - pep.seq : Vecteur contenant les séquences de peptides à analyser.
# - method : Méthode utilisée pour calculer le logo. Options :
#     - "prob" (par défaut) : Basée sur les fréquences des acides aminés (descriptive).
#     - "bits" : Utilise l'information en bits pour des analyses plus informatives.
# - title : Titre du graphique (par défaut "Motif").
#
# Fonctionnalités :
# 1. Nettoie les séquences en supprimant les astérisques (*) et les valeurs manquantes (NA).
# 2. Vérifie que toutes les séquences ont la même longueur.
# 3. Génère un logo de séquence à l'aide de ggseqlogo avec des positions dynamiques centrées sur 0.
#
# Sortie :
# - Un graphique ggplot représentant le logo de séquence.

#===============================================================================  Couleur R
colors_64 <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", 
               "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
               "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
               "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
               "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
               "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
               "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
               "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
               "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
               "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
               "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
               "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
               "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
               "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
               "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
colorLookup <- c(selTable = "#4DBBD5FF", selPlot = "#00A087FF", selPFC = "#E64B35FF", selGenes = "#3C5488FF")

# ======================================================================  Charger et processer les donnees  =================================================================================================##
# Charger le tableau 
Table1_INSERM_00052349_IMAC_FINAL_7_7_Original <- read_excel("Data/Table1_INSERM_00052349_IMAC_FINAL_7_7_Original.xlsx", 
                                                            # range = cell_cols("D:I"), 
                                                             sheet = "Summary", skip = 5)

# Définir explicitement les noms des colonnes à partir de la première ligne
Table1_Summary_Original <- Table1_INSERM_00052349_IMAC_FINAL_7_7_Original
colnames(Table1_Summary_Original) <- Table1_Summary_Original[4, ]

# Supprimer la ligne utilisée pour les noms de colonnes
phospho_data <- Table1_Summary_Original[-c(1:4), ]

# Conversion en dataframe
phospho_data <- as.data.frame(phospho_data)

# ================================================= Selectionner les donnees pour la suite ================
# Fonction pour selectionner les donnees
Select_Data <- function(phospho_data, column = 4:9){
  data = phospho_data %>%
    # Sélectionner les colonnes spécifiées (par défaut les colonnes 4 à 9)
    dplyr::select(column)
  
  # Convertir toutes les colonnes sélectionnées en numériques
  data <- data %>% mutate(across(everything(), as.numeric))
  
  # Définir les noms de lignes de la nouvelle donnée à partir de la colonne 'Index in Detail'
  rownames(data) <- phospho_data$`Index in Detail`
  
  return(data)  
}

# FC
Normalized_FC <- Select_Data(phospho_data, column = 4:9)
  
#Log2FC
Normalized_log2FC <- Select_Data(phospho_data, column = 98:103)

# Pvalue
Pvalue <- Select_Data(phospho_data, column = 10:15)

# Abundance
Abundance_Normalized <- Select_Data(phospho_data, column = 50:65)

# Remplacer les valeurs non numériques et convertir en numérique
Abundance <- as.data.frame(
  lapply(Abundance_Normalized, function(col) {
    numeric_col <- as.numeric(gsub("[^0-9\\.]", "", col))  # Nettoyage et conversion
    replace(numeric_col, is.na(numeric_col), 20)           # Remplacer NA par 20
  })
)
rownames(Abundance) <- rownames(Abundance_Normalized)
names(Abundance) <- names(Abundance_Normalized)

# Metdata
metadata_phospo <- phospho_data %>%
  dplyr::select(18, 19, 20, 2, 3, 23, 29, 26, 16, 17, 24)                             
# Définit les noms des lignes à partir de la colonne "Index in Detail"
rownames(metadata_phospo) <- phospho_data$`Index in Detail`


# Suppression des contenus après le point-virgule dans les colonnes Gene et Protein Name
metadata_phospo_cleaned <- metadata_phospo %>%
  mutate(
    Gene = str_extract(Gene, "^[^;]+"),                # Conserve tout avant le premier point-virgule
    `Protein Name` = str_extract(`Protein Name`, "^[^;]+") # Conserve tout avant le premier point-virgule
  )
rownames(metadata_phospo_cleaned) <- rownames(metadata_phospo)

# ======================================================================  Parametres  =================================================================================================##

cutoff.pv <- 0.05

cutoff.fc <- 1.8

# Permet de voir et selectionner en entier le nom de le condition (entrer un a un)
condition_1 <- Conditions_Selection(condition = NULL, Normalized_log2FC, Pvalue)

print(condition_1)

condition_2 <- Conditions_Selection(condition = NULL, Normalized_log2FC, Pvalue) 

print(condition_2)

# ======================================================================  Analyse differentiel  =================================================================================================##
# Differtial 
Condition_1_Diff <- Select_Diff(
  log2FC = Normalized_log2FC, 
  pvalue = Pvalue, 
  metadata = metadata_phospo_cleaned, 
  condition = condition_1,
  cutoff.fc = log2(cutoff.fc),
  cutoff.pv = cutoff.pv
)

Condition_2_Diff <- Select_Diff(
  log2FC = Normalized_log2FC, 
  pvalue = Pvalue, 
  metadata = metadata_phospo_cleaned, 
  condition = condition_2,
  cutoff.fc = log2(cutoff.fc),
  cutoff.pv = cutoff.pv
)

# ==============================================================================  Volcano-plots 

#Options du script : Possibilité de changer les seuils du Fold Change et de la p-value (par défaut fixer le Fold Change à |1.8| et la p-value à 0.05)
Condition_1_plot <- Plot_Volcano(Condition_1_Diff, cutoff.fc = cutoff.fc, cutoff.pv = cutoff.pv, title = condition_1, save = FALSE)
#Afficher le volcano
print(Condition_1_plot)

Condition_2_plot <- Plot_Volcano(Condition_2_Diff, cutoff.fc = cutoff.fc, cutoff.pv = cutoff.pv, title = condition_2)

print(Condition_2_plot)

# ==============================================================================  Diagramm Venn   
#Objectif : Dresser la liste des phopho-peptides communs/ nombre de protéines total commun à différentes conditions et affiner la liste par application de différents filtres

Select_Direction <- function(data, Direction = c("UP", "DOWN")){
  if (Direction == "UP") {  # Si la direction est 'UP'
    dir = data %>% 
      filter(Direction == "UP") %>%  # Filtrer les protéines upregulées
      pull(Protein_Name)  # Extraire les noms des protéines
  } else if (Direction == "DOWN") {  # Si la direction est 'DOWN'
    dir = data %>% 
      filter(Direction == "DOWN") %>%  # Filtrer les protéines downregulées
      pull(Protein_Name)  # Extraire les noms des protéines
  } else {  # Si la direction n'est ni 'UP' ni 'DOWN'
    stop("La direction sélectionnée n'est pas valable")
  }
  
  return(dir)  # Retourner les noms des protéines correspondant à la direction choisie
}

# Selectionner la direction UP ou DOWN pour chaque condition et comparer dans le Venn Diagramm
list_1 <- Select_Direction(Condition_1_Diff, Direction = "UP")
list_2 <- Select_Direction(Condition_2_Diff, Direction = "DOWN")

# Mettre la direction selectionner pour l'annotation du plot
dir_1 = "Up-regulated"
dir_2 = "Down-regulated"

draw.pairwise.venn(
  area1 = length(list_1), 
  area2 = length(list_2),
  cross.area = length(intersect(list_1, list_2)),
  fill = c("red", "blue"),
  category = c(paste0(condition_1, "\n", dir_1), paste0(condition_2, "\n", dir_2)), # Ajout de \n pour retour à la ligne
  lty = "blank",
  scaled = FALSE,
  cat.cex = 1.2, # Taille des légendes
  cat.fontfamily = "sans", # Police des légendes
  cat.pos = c(-30, 30), # Position des légendes
  cex = 1.5 # Taille du texte des cercles
)

# ==============================================================================  Heatmap 
# Calcul des moyennes des réplicats et réorganisation des données
Abundance_means <- Abundance %>% 
  mutate(
    Control_1 = rowMeans(select(., starts_with("Control 1"))),
    Control_2 = rowMeans(select(., starts_with("Control 2"))),
    RSK2_1 = rowMeans(select(., starts_with("RSK2 1"))),
    RSK2_2 = rowMeans(select(., starts_with("RSK2 2"))),
    Crush_shControl_1 = rowMeans(select(., starts_with("Crush+shControl 1"))),
    Crush_shControl_2 = rowMeans(select(., starts_with("Crush+shControl 2"))),
    Crush_shRSK2_1 = rowMeans(select(., starts_with("Crush+shRSK2 1"))),
    Crush_shRSK2_2 = rowMeans(select(., starts_with("Crush+shRSK2 2")))
  ) %>% 
  select(Control_1, Control_2, RSK2_1, RSK2_2, Crush_shControl_1, Crush_shControl_2, Crush_shRSK2_1, Crush_shRSK2_2)

# Annotations des états et conditions
Confidence_annotation <- data.frame(Confidence = metadata_phospo_cleaned$`Confidently Localized`, stringsAsFactors = FALSE, 
                                    row.names = rownames(Abundance_means))
Condition_annotation <- data.frame(Condition = c("Control", "Control", "RSK2", "RSK2", "Crush_shControl", "Crush_shControl", 
                                                 "Crush_shRSK2", "Crush_shRSK2"), stringsAsFactors = FALSE)
# Créer un vecteur pour les couleurs des états cellulaires
Confidence_colors <- structure(c("lightblue", "thistle"#, "pink", "lightgreen", "mediumslateblue"
                            ), 
                          names = unique(Confidence_annotation$Confidence))

Condition_annotation_colors <- structure(c("green", "orange", "pink", "purple"),
                                names = unique(Condition_annotation$Condition))

# Créer un vecteur pour les couleurs des états cellulaires
condition_colors <- structure(c("lightgreen", "lightblue", "thistle", "pink"), 
                              names = unique(Condition_annotation$Condition))

# Scale the data (Z-score scaling for rows, which are genes)
scaled_data <- scale(t(as.matrix(Abundance_means)))

# Créer un heatmap avec les états et conditions
Heatmap(
  scaled_data,
  name = "Expression",
  col = colorRamp2(c(min(scaled_data), 0 , max(scaled_data)), c("darkgreen", "lightgreen", "yellow")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = factor(Confidence_annotation$Confidence, levels = c("YES", "NO")),  # Split par groupe de gènes
  row_split = factor(Condition_annotation$Condition, levels = c("Control", "RSK2", "Crush_shControl", 
                                                                   "Crush_shRSK2")),
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = HeatmapAnnotation(
    Confidence_Localized = factor(Confidence_annotation$Confidence, levels = c("YES", "NO")),

    col = list(Confidence_Localized = Confidence_colors),
    show_legend = FALSE
  ),
  # bottom_annotation = HeatmapAnnotation(
  #   Condition = factor(Condition_annotation$Condition, levels = c("Control", "RSK2", "Crush_shControl", 
  #                                                                 "Crush_shRSK2")),
  #   col = list(Condition = condition_colors),
  #   show_legend = TRUE  # Désactiver la légende supplémentaire
  # ),
  left_annotation = rowAnnotation(
    Condition = factor(Condition_annotation$Condition, levels = c("Control", "RSK2", "Crush_shControl", 
                                                                    "Crush_shRSK2")),
    col = list(Condition = Condition_annotation_colors),
    show_legend = FALSE
  )

)

# ==============================================================================  Others plot  

# Utilisation de la fonction
plot_abundance_protein_site("Atat1", Abundance_means, metadata_phospo_cleaned)

# ======================================================================  Type de regulation & Consensus  ============================================================================================##

#Objectif : Identifier si la régulation est directe par phosphorylation (1 protéine = 1 phosphopeptide) versus Régulation transcription/traduction/synthèse protéique (1 protéine = plusieurs phospho-peptides). 
#Consensus séquence de RSK2 (séquence théorique dans laquelle chaque acide-aminé représente le résidu le plus fréquemment rencontré parmi les divers substrats de RSK2 connues) : RXRXXS*/T*  
# -	Position 0 = site de phosphorylation (S*/T* ; pour Sérine ou Thréonine ; identifié par une astérisque)
# -	Positions -1 et -2 = n’importe quels acides-aminés peuvent être présents (XX)
# -	Position -3 : arginine (R-3)
# -	Position -4 : n’importe quel acide-aminé peut être présent (X)
# -	Position -5 : arginine (R-5)
# ==============================================================================  Sur les donnes diff

# Rerchercher les substrats direct et les positions de R-3 et -5 dans les donnees differentiel des 2 conditions la fonction "Make_Rpos_Consensus"
Consensus_Condition_1 = Make_Rpos_Consensus(Condition_1_Diff, Protein_Name = "Protein_Name", Peptide_Localized = "Peptide_Localized")

Consensus_Condition_2 = Make_Rpos_Consensus(Condition_2_Diff, Protein_Name = "Protein_Name", Peptide_Localized = "Peptide_Localized")

# ==============================================================================  Sur intersection de venn Diagram
# Rechercher les substrats direct et les positions de R-3 et -5 dans la liste commune du Venn Diagramm

signatures <- list(UP = list_1,
                   DOWN = list_2
)
Inters.sign <- Intersect(signatures)

# Selectionner les proteines a l'intersection

phosppho_intersec <- metadata_phospo_cleaned %>%
  dplyr::filter(`Protein Name` %in% Inters.sign)

Consensus_Reg <- Make_Rpos_Consensus(phosppho_intersec, Protein_Name = "Protein Name", Peptide_Localized = "Localized -7/+7 Peptide")

# Étape 2 : Filtrer les protéines pour chaque catégorie
substrat_direct <- Consensus_Reg %>%
  dplyr::filter(substrat=="direct")

substrat_2 <- Consensus_Reg %>%
  dplyr::filter(substrat=="two")

substrat_3 <- Consensus_Reg %>%
  dplyr::filter(substrat=="three")

substrat_plus <- Consensus_Reg %>%
  dplyr::filter(substrat=="plus")


# ======================================================================  LOGO =================================================================================================##

#-------------------------
# Exemple d'utilisation de la fonction. Vous pouvez selectionner un, plusieurs peptide et dessiner leurs logo. Vous pouvez prendre directment les peptides d'un sous types de donnees egalement. 
# Plot logo pour  tous les substrats direct  
pep.seq <- substrat_direct$`Localized -7/+7 Peptide`
method <- "bits"
logo_plot_direct <- Plot_logo(pep.seq = pep.seq, method = method, title = "Motif")
print(logo_plot_direct)

# Plot logo pour  tous les substrats inderect (2)  
pep.seq <- substrat_2$`Localized -7/+7 Peptide`
method <- "bits"
logo_plot_2 <- Plot_logo(pep.seq = pep.seq, method = method, title = "Motif")
print(logo_plot_2)

# Plot logo pour  un seul peptides
logo_plot <- Plot_logo(pep.seq = "LPDPPLES*EDEDEEG", method = "prob", title = "Motif")
print(logo_plot)

# Plot logo pour  plusieurs peptides
logo_plot_mult <- Plot_logo(pep.seq = c("LPDPPLES*EDEDEEG", "KIPLHVFS*PKYTKLR", "AAKSTPGT*PLAKAKA") , method = "bits", title = "Motif")
print(logo_plot_mult)


# ======================================================================  Comparaison =================================================================================================##
# URL pour télécharger l'ontologie GO en format OBO
go_url <- "http://purl.obolibrary.org/obo/go.obo"

# Télécharger et charger l'ontologie dans R
go_ontology <- get_ontology(go_url, extract_tags = "minimal")

# URL AmiGO Solr API pour "axon regeneration" (GO:0031103)
axon_url <- "https://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity_label,bioentity_name,annotation_class_label,annotation_class,assigned_by,taxon,evidence_type,evidence_with,panther_family,type,bioentity_isoform,reference,date&fq=document_category:%22annotation%22&fq=isa_partof_closure:%22GO:0031103%22&q=*:*"

# Télécharger les données
axon_data <- read.csv(axon_url, sep = ",", header = TRUE)

Select_Common_Genes <- function(df, ref){
  # Comparer les gènes expérimentaux avec ceux de la base axon_data
  common_genes <- intersect(df$Gene, ref$bioentity_label)
  message("La liste des genes en commun est la suivante :")
  print(common_genes)
  # Créer un tableau contenant tous les gènes avec les informations correspondantes
  common_genes_df <- df %>%
    filter(Gene %in% common_genes)
  
  return(common_genes_df)
}

common_genes_df <- Select_Common_Genes(df = Consensus_Reg, ref = axon_data)

# ============================================================================== Fin ===================================================================================================================#


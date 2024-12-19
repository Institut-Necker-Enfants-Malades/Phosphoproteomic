
# Les fonctions 
#=============================================================================== Conditions_Selection()
Conditions_Selection <- function(condition = NULL, log2FC, pvalue){
  # Si aucune condition n'est spécifiée, afficher les conditions disponibles et demander une sélection
  if (is.null(condition)) {
    message("Aucune condition spécifiée. Veuillez choisir une condition parmi les suivantes :")
    available_conditions <- intersect(colnames(log2FC), colnames(pvalue))
    for (i in seq_along(available_conditions)) {
      message(paste0("[", i, "] ", available_conditions[i]))
    }
    
    # Interface pour la sélection
    repeat {
      user_choice <- as.numeric(readline("Entrez le numéro de la condition souhaitée : "))
      if (!is.na(user_choice) && user_choice >= 1 && user_choice <= length(available_conditions)) {
        condition <- available_conditions[user_choice]
        break
      } else {
        message("Sélection invalide. Veuillez entrer un numéro valide.")
      }
    }
    message(paste("Condition sélectionnée : ", condition))
  }
  
  # Vérification de l'existence de la condition dans log2FC et pvalue
  if (!(condition %in% colnames(log2FC)) || !(condition %in% colnames(pvalue))) {
    stop(paste("Erreur : La condition spécifiée ('", condition, "') n'existe pas dans les colonnes de log2FC ou pvalue. Veuillez vérifier les noms des colonnes.", sep = ""))
  }
return(condition)
}
#=============================================================================== Select_Diff()
Select_Diff <- function(log2FC, pvalue, metadata, condition = "Crush+shControl : Control", cutoff.fc = 1.8, cutoff.pv = 0.05) {
  # Vérification des colonnes requises dans metadata_phospo_cleaned
  required_columns <- c("Protein Name", "Gene", "Index in Detail", "Localized -7/+7 Peptide", "Confidently Localized", "URL")
  if (!all(required_columns %in% colnames(metadata_phospo_cleaned))) {
    stop("Les colonnes nécessaires sont manquantes dans metadata_phospo_cleaned : ", 
         paste(setdiff(required_columns, colnames(metadata_phospo_cleaned)), collapse = ", "))
  }
  
  # Vérification de l'existence de la condition dans log2FC et pvalue
  if (!(condition %in% colnames(log2FC)) || !(condition %in% colnames(pvalue))) {
    stop(paste("Erreur : La condition spécifiée ('", condition, "') n'existe pas dans les colonnes de log2FC ou pvalue. Veuillez vérifier les noms des colonnes.", sep = ""))
  }
  # Fusionner les données à partir de la condition spécifiée
  data <- data.frame(
    Protein_Name = metadata$`Protein Name`,
    Gene = metadata_phospo_cleaned$Gene,
    Index_in_Detail = metadata_phospo_cleaned$`Index in Detail`,
    log2FC = as.numeric(log2FC[[condition]]),          # Sélectionner la colonne de log2FC correspondant à la condition
    pvalue = pvalue[[condition]],                     # Sélectionner la colonne de pvalue correspondant à la condition
    Peptide_Localized = metadata_phospo_cleaned$`Localized -7/+7 Peptide`,
    Confidently_Localized = metadata_phospo_cleaned$`Confidently Localized`,
    URL = metadata_phospo_cleaned$URL,
    row.names = rownames(log2FC)                      # Utiliser les noms de lignes de log2FC
  )
  
  # Gestion des valeurs manquantes et conversion robuste
  data_cleaned <- data %>%
    mutate(pvalue = ifelse(pvalue == "-", NA, suppressWarnings(as.numeric(pvalue))))
  
  if (any(is.na(data_cleaned$log2FC))) {
    warning("Certaines valeurs de log2FC n'ont pas pu être converties en numérique.")
  }
  if (any(is.na(data_cleaned$pvalue))) {
    warning("Certaines valeurs de p-value n'ont pas pu être converties en numérique.")
  }
  
  # Ajouter une colonne "Direction" initialisée à "NO"
  data_cleaned$Direction <- "NO"
  
  # Définir "UP" pour les gènes avec log2FC > cutoff.fc et p-value < cutoff.pv
  data_cleaned$Direction[data_cleaned$log2FC > cutoff.fc & data_cleaned$pvalue < cutoff.pv] <- "UP"
  
  # Définir "DOWN" pour les gènes avec log2FC < -cutoff.fc et p-value < cutoff.pv
  data_cleaned$Direction[data_cleaned$log2FC < -cutoff.fc & data_cleaned$pvalue < cutoff.pv] <- "DOWN"
  
  # Résumé des résultats UP/DOWN/NO
  summary_results <- data.frame(
    Direction = c("UP", "DOWN", "NO"),
    Count = table(data_cleaned$Direction)
  )
  print(summary_results)
  
  # Retourner le tableau de données modifié
  return(data_cleaned)
}
#=============================================================================== Plot_Volcano()
Plot_Volcano <- function(data, cutoff.fc, cutoff.pv, title = "Title", save = FALSE) {
  # Créer le graphique volcan avec ggplot
  volcano_gplot <- ggplot(data=data, aes(x=log2FC, y=-log10(pvalue), col=Direction, label=Protein_Name)) + 
    geom_point() +
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-log2(cutoff.fc), log2(cutoff.fc)), col="red") +
    geom_hline(yintercept=-log10(cutoff.pv), col="red") +
    labs(title = title)
  
  # Vérifier si l'utilisateur souhaite sauvegarder le graphique
  if (save) {
    message("Vous avez choisi de sauvegarder le graphique.")
    
    # Ouvrir un volet pour choisir le dossier de sauvegarde
    selected_dir <- rstudioapi::selectDirectory(caption = "Sélectionnez un dossier pour sauvegarder le fichier")
    if (is.null(selected_dir) || selected_dir == "") {
      message("Aucun dossier sélectionné. Annulation de la sauvegarde.")
      return(volcano_gplot)
    }
    
    # Demander un nom de fichier
    file_name <- readline(prompt = "Entrez le nom du fichier (sans extension) : ")
    if (!nzchar(file_name)) {
      message("Nom de fichier invalide. Annulation de la sauvegarde.")
      return(volcano_gplot)
    }
    
    # Ajouter l'extension .png au nom de fichier
    file_name <- paste0(file_name, ".png")
    save_path <- file.path(selected_dir, file_name)
    
    # Sauvegarder le graphique
    tryCatch({
      ggsave(
        filename = save_path, 
        plot = volcano_gplot, 
        device = "png", 
        width = 10, 
        height = 7, 
        bg = "white"
      )
      message(paste("Graphique sauvegardé avec succès à l'emplacement :", save_path))
    }, error = function(e) {
      message("Erreur lors de la sauvegarde : ", e$message)
    })
    
    # Retourner NULL après la sauvegarde
    return(NULL)
  }
  
  # Convertir le graphique ggplot en graphique interactif avec plotly
  volcano_gplotly <- volcano_gplot %>% 
    ggplotly(tooltip = 'all')
  
  # Retourner le graphique interactif
  return(volcano_gplotly)
}


#=============================================================================== Intersect()
Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

#=============================================================================== plot_abundance_protein_site()
# Fonction de visualisation de l'abondance pour une protéine donnée
plot_abundance_protein_site <- function(protein_name, abundance_data, metadata_phospo_cleaned) {
  # Calcul des moyennes des réplicats
  replicate_means <- abundance_data %>%
    mutate(
      Control = rowMeans(select(., starts_with("Control")), na.rm = TRUE),
      RSK2 = rowMeans(select(., starts_with("RSK2")), na.rm = TRUE),
      Crush_shControl = rowMeans(select(., starts_with("Crush_shControl")), na.rm = TRUE),
      Crush_shRSK2 = rowMeans(select(., starts_with("Crush_shRSK2")), na.rm = TRUE)
    ) %>%
    select(Control, RSK2, Crush_shControl, Crush_shRSK2)
  
  # Ajout des métadonnées
  rownames(replicate_means) <- rownames(metadata_phospo_cleaned)
  merged <- replicate_means
  merged$Gene <- metadata_phospo_cleaned$Gene
  merged$Site <- metadata_phospo_cleaned$Site
  
  # Filtrer les données pour la protéine demandée
  protein_data <- merged %>% filter(Gene == protein_name)
  
  # Transformation en format long
  protein_data_long <- protein_data %>%
    pivot_longer(cols = c(Control, RSK2, Crush_shControl, Crush_shRSK2),
                 names_to = "Condition",
                 values_to = "Abundance")
  
  # Création du graphique
  ggplot(protein_data_long, aes(x = Site, y = Abundance, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    labs(
      title = paste("Abondance pour la protéine", protein_name),
      x = "Sites",
      y = "Abondance moyenne",
      fill = "Condition"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#=============================================================================== Make_Rpos_Consensus()
Make_Rpos_Consensus <- function(data, Protein_Name = "Protein_Name", Peptide_Localized = "Peptide_Localized") {
  # Vérification des colonnes
  if (!Protein_Name %in% colnames(data)) {
    stop(paste("La colonne", Protein_Name, "n'existe pas dans les données."))
  }
  if (!Peptide_Localized %in% colnames(data)) {
    stop(paste("La colonne", Peptide_Localized, "n'existe pas dans les données."))
  }
  
  # Étape 1 : Compter les occurrences de chaque protéine
  protein_counts <- data %>%
    dplyr::count(.data[[Protein_Name]], name = "Peptides_Count") %>%
    dplyr::arrange(desc(Peptides_Count))  # Trier par ordre décroissant
  
  # Étape 2 : Ajouter une colonne 'substrat' basée sur les catégories
  data <- data %>%
    dplyr::left_join(protein_counts, by = Protein_Name) %>%  # Ajouter la colonne 'Count'
    dplyr::mutate(
      substrat = dplyr::case_when(
        Peptides_Count == 1 ~ "direct",
        Peptides_Count == 2 ~ "two",
        Peptides_Count == 3 ~ "three",
        Peptides_Count > 3 ~ "plus",
        TRUE ~ NA_character_ # Par défaut, si aucune condition n'est remplie
      )
    )
  
  # Étape 3 : Identifier les positions R
  Consensus <- data %>%
    dplyr::mutate(
      R_positions = purrr::map_chr(.data[[Peptide_Localized]], function(peptide) {
        # Retirer les espaces pour éviter des erreurs
        peptide <- gsub(" ", "", peptide)
        
        # Identifier la position du site de phosphorylation (*)
        site_pos <- regexpr("\\*", peptide)
        
        # Vérifier les résidus autour du site, si la position est valide
        if (site_pos > 5) { # Vérifie qu'il y a au moins 5 AA à gauche
          left_seq <- substr(peptide, site_pos - 6, site_pos - 1)
          r_neg_3 <- substr(left_seq, 3, 3) == "R" # Résidu à la position -3 (4e AA)
          r_neg_5 <- substr(left_seq, 1, 1) == "R" # Résidu à la position -5 (2e AA)
          
          # Créer une annotation
          if (r_neg_3 && r_neg_5) {
            return("R-3 and R-5")
          } else if (r_neg_3) {
            return("R-3")
          } else if (r_neg_5) {
            return("R-5")
          } else {
            return("No R at -3/-5")
          }
        } else {
          return("Invalid peptide length")
        }
      })
    )
  
  # Étape 4 : Trier les résultats
  Consensus_sorted <- Consensus %>%
    arrange(desc(substrat == "direct"))  # Trier de manière décroissante sur "direct"
  
  return(Consensus_sorted)
}

#=============================================================================== Plot_logo()
Plot_logo <- function(pep.seq, method = c("prob", "bits"), title = "Motif") {
  # Method :
  #"prob" : Pour une analyse descriptive et intuitive des fréquences.
  #"bits" : Pour une analyse plus informative, en particulier pour identifier des régions hautement conservées ou spécifiques.
  #-----------------
  # Vérification et nettoyage des séquences
  sequences <- na.omit(pep.seq)  # Supprimer les séquences NA
  sequences <- gsub("\\*", "", sequences)  # Supprimer les astérisques (*)
  
  # Vérification des longueurs des séquences
  seq_lengths <- nchar(sequences)
  if (length(unique(seq_lengths)) > 1) {
    stop("Toutes les séquences doivent avoir la même longueur.")
  }
  
  # Sélection de la méthode (par défaut "prob")
  method <- match.arg(method)
  
  # Création du plot avec ggseqlogo
  ggseqlogo(list(Peptide = sequences), method = method) +
    ggtitle(title) +
    theme_minimal(base_size = 14) +
    scale_x_continuous(
      breaks = seq(1, seq_lengths[1]),  # Positions dynamiques
      labels = seq(-floor(seq_lengths[1] / 2), ceiling(seq_lengths[1] / 2) - 1)  # Centre 0
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
}


#=============================================================================== End

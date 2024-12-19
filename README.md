# `Analyse de Phosphoprotéomique RSK2`

## Description
Ce script R permet d'analyser des données de phosphoprotéomique, notamment :
- Analyse différentielle des phosphopeptides
- Génération de volcano plots
- Création de diagrammes de Venn pour comparer les conditions
- Visualisation par heatmap des niveaux d'expression
- Analyse des motifs de séquence (logos)
- Comparaison avec des données d'axon regeneration

## Prérequis

### `Installation des packages`
Le script inclut une fonction d'installation automatique des packages nécessaires :
```R
required_packages <- c(
  "ggplot2", "ggrepel", "dplyr", "readxl", "VennDiagram", 
  "tidyverse", "plotly", "ComplexHeatmap", "circlize", 
  "magrittr", "ggseqlogo", "rstudioapi", "ontologyIndex"
)

```

## `Structure des fichiers`
```R
project/
├── scripts/
│   ├── main_script.R
│   └── Util_Phospo.R
└── Data/
    └── Table1_INSERM_00052349_IMAC_FINAL_7_7_Original.xlsx

```
## `Utilisation`

### Préparation des données

Placez votre fichier Excel dans le dossier Data/
Le fichier doit contenir une feuille "Summary" avec les données de phosphoprotéomique


### Configuration initiale

Modifiez les paramètres selon vos besoins :

```R
Copycutoff.pv <- 0.05    # Seuil de p-value
cutoff.fc <- 1.8     # Seuil de Fold Change
```
### Analyses disponibles

- Analyse différentielle (Select_Diff)
- Volcano plots (Plot_Volcano)
- Diagrammes de Venn pour comparer les conditions
- Heatmap des abondances
- Analyse des motifs de séquence
- Recherche de consensus RSK2 (RXRXXS*/T*)
  
### Visualisations

- Volcano plots des résultats différentiels
- Diagramme de Venn des protéines régulées
- Heatmap des abondances normalisées
- Logo plots des motifs de séquence

### Format des données
Le fichier Excel d'entrée doit contenir :

- Données d'expression normalisées
- Log2 Fold Change
- P-values
- Métadonnées (Gene, Protein Name, etc.)

### Fonctions principales

```R
### 1. Conditions_Selection()

#Description : Permet d'afficher les conditions disponibles et demander une sélection pour eviter les erreurs de frappes.

### 2. Select_Diff()
Description : Identifie les protéines/peptides différentiellement exprimés basé sur le log2 fold-change et la p-value.

### 2. Plot_Volcano()
Description : Génère un volcano plot pour visualiser les changements d'expression.

### 3. Make_Rpos_Consensus()
Description : Analyse les séquences peptidiques pour identifier le motif consensus RSK2 (RXRXXS*/T*).

### 4. Plot_logo()
Description : Crée une visualisation de logo pour les motifs de séquence.

### Notes d'utilisation
- Les fonctions sont conçues pour être utilisées séquentiellement dans l'analyse
- Les paramètres par défaut sont optimisés pour l'analyse RSK2 mais peuvent être ajustés
- La visualisation peut être personnalisée via les paramètres des fonctions de plotting
```
## Notes importantes

- Les seuils de significance peuvent être ajustés (cutoff.pv et cutoff.fc)
- L'analyse de consensus RSK2 recherche le motif RXRXXS*/T*
- Les données d'axon regeneration sont téléchargées depuis [AmiGO](https://amigo.geneontology.org/amigo/term/GO:0031103) Solr API


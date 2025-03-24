### PARTIE 1. DETERMINANTS DES SALAIRES (PAR LES MOINDRES CARRES ORDINAIRES) 

# Installation et chargement des bibliothèques nécessaires
install.packages("haven")
install.packages("dplyr")

library(haven)
library(dplyr)

# Chargement de la base de données ECMOSS
data <- read_dta("Ecmoss_2006.dta")
selected_data <- data %>% select(s_net, duree, nbheur, age_r, sexe_r, statut, strate, qs25_r, qs26_r)

# Convertion des variables qualitatives en facteurs
selected_data <- selected_data %>%
  mutate(
    sexe_r = as.factor(sexe_r),
    statut = as.factor(statut),
    strate = as.factor(strate),
    qs25_r = as.factor(qs25_r),
    qs26_r = as.factor(qs26_r)
  )

# Stats des pour les variables quantis
quantitative_summary <- summary(selected_data[, sapply(selected_data, is.numeric)])
print(quantitative_summary)

# Modalités uniques pour les variables qualis
qualitative_summary <- sapply(selected_data[, sapply(selected_data, is.factor)], function(x) length(unique(x)))
print(qualitative_summary)

# Distribution du salaire net
hist(selected_data$s_net, 
     main = "Distribution du Salaire Net", 
     xlab = "Salaire Net (€)", 
     col = "lightblue", 
     breaks = 30)

# Distribution de la durée travaillée
hist(selected_data$duree, 
     main = "Distribution de la Durée Travaillée", 
     xlab = "Nombre de Jours Travaillés", 
     col = "lightgreen", 
     breaks = 30)

# Distribution de l'âge
hist(selected_data$age_r, 
     main = "Distribution de l'Âge des Salariés", 
     xlab = "Âge", 
     col = "lightcoral", 
     breaks = 30)

# Répartition par sexe
barplot(table(selected_data$sexe_r), 
        main = "Répartition par Sexe", 
        col = c("skyblue", "pink"), 
        xlab = "Sexe", 
        ylab = "Nombre d'Observations")

# Répartition par statut
barplot(table(selected_data$statut), 
        main = "Répartition par Statut", 
        col = "orange", 
        xlab = "Statut", 
        ylab = "Nombre d'Observations")

# Répartition par niveau de qualification
barplot(table(selected_data$qs25_r), 
        main = "Répartition par Niveau de Qualification", 
        col = "purple", 
        xlab = "Niveau de Qualification", 
        ylab = "Nombre d'Observations", 
        las = 2)

# Répartition par diplôme le plus élevé
barplot(table(selected_data$qs26_r), 
        main = "Répartition par Diplôme le Plus Élevé", 
        col = "blue", 
        xlab = "Diplôme", 
        ylab = "Nombre d'Observations", 
        las = 2)

# Pourcentage de valeurs manquantes pour chaque variable
missing_values <- sapply(selected_data, function(x) sum(is.na(x)) / length(x) * 100)
print(missing_values)

# Remplacement par la médiane
selected_data <- selected_data %>%
  mutate(
    duree = ifelse(is.na(duree), median(duree, na.rm = TRUE), duree),
    nbheur = ifelse(is.na(nbheur), median(nbheur, na.rm = TRUE), nbheur),
    age_r = ifelse(is.na(age_r), median(age_r, na.rm = TRUE), age_r),
    s_net = ifelse(is.na(s_net), median(s_net, na.rm = TRUE), s_net),
  )

# Vérification
missing_values <- sapply(selected_data, function(x) sum(is.na(x)))
print(missing_values)

# Modèlisation par les MCO
model <- lm(s_net ~ duree + nbheur + age_r + sexe_r + statut + strate + qs25_r + qs26_r, data = selected_data)
summary(model)

############
############
############

### PARTIE 2. CONSTRUCTION D’UNE RELATION BINAIRE ET FERMETURE TRANSITIVE

# Installation et chargement des bibliothèques nécessaires (pour celles non faites précédemment)
install.packages("haven")
install.packages("dplyr")

library(dplyr)
library(ggplot2)

# Chargement de la base de donnees
data <- read.csv("sgl-arbres-urbains-wgs84.csv")

# Remplacement des valeurs manquantes par la moyenne
data$hauteur <- ifelse(is.na(data$hauteur), mean(data$hauteur, na.rm = TRUE), data$hauteur)
data$diametre <- ifelse(is.na(data$diametre), mean(data$diametre, na.rm = TRUE), data$diametre)

# Histogrammes observer "hauteur" et "diametre"
ggplot(data, aes(x = hauteur)) + 
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) + 
  labs(title = "Distribution de la hauteur", x = "Hauteur (m)", y = "Fréquence")

ggplot(data, aes(x = diametre)) + 
  geom_histogram(binwidth = 0.1, fill = "green", color = "black", alpha = 0.7) + 
  labs(title = "Distribution du diamètre", x = "Diamètre (m)", y = "Fréquence")

# Definition des seuils 
s_h <- 600  # Pour la hauteur 
s_d <- 30   # Pour le diametre

# Construction de la relation binaire Q
relation_binaire <- function(df) {
  n <- nrow(df)
  Q <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        hauteur_cond <- abs(df$hauteur[i] - df$hauteur[j]) <= s_h
        diametre_cond <- abs(df$diametre[i] - df$diametre[j]) <= s_d
        
        if (hauteur_cond && diametre_cond) {
          Q[i, j] <- 1  
        } else if (df$hauteur[i] > df$hauteur[j] + s_h && df$diametre[i] > df$diametre[j] + s_d) {
          Q[i, j] <- 2  
        }
      }
    }
  }
  return(Q)
}

# Application de la fonction pour obtenir la matrice Q
matrice_Q <- relation_binaire(data)

# Fonction pour calculer la fermeture transitive
fermeture_transitive <- function(Q) {
  n <- nrow(Q)
  fermeture <- Q
  
  for (k in 1:n) {
    for (i in 1:n) {
      for (j in 1:n) {
        if (fermeture[i, k] > 0 && fermeture[k, j] > 0) {
          fermeture[i, j] <- max(fermeture[i, j], 1)  # Marque la relation transitive
        }
      }
    }
  }
  return(fermeture)
}

# Calcul de la fermeture transitive
closure_matrix <- fermeture_transitive(matrice_Q)

# Sauvegarde des matrices
write.csv(matrice_Q, "relation_binaire_Q.csv", row.names = FALSE)
write.csv(closure_matrix, "fermeture_transitive_Q.csv", row.names = FALSE)

############
############
############

### PARTIE 3. OPTIMISATION SOUS R

## Étape 1 : Définition de la fonction à minimiser

# Cette fonction représente notre objectif à minimiser :
# f(x, y, z) = x^2 + y^2 + z^2 - 2x - 3y + z
f <- function(x, y, z) {
  return(x^2 + y^2 + z^2 - 2*x - 3*y + z)
}

## Étape 2 : Calcul des sommets en résolvant les contraintes

# Les contraintes définissent la région admissible dans l’espace (x, y, z). Elles sont les suivantes :
# 1. x + 2y + z ≤ 4
# 2. x + y ≥ 1
# 3. x, y, z ≥ 0

# Cas 1 : Résolution des contraintes en égalité
# Nous utilisons les contraintes x + y = 1 et x + 2y + z = 4 pour exprimer x et z en fonction de y.

calc_sommet_cas1 <- function(y) {
  # Calcul de x à partir de la contrainte x + y = 1
  x <- 1 - y
  
  # Calcul de z à partir de la contrainte x + 2y + z = 4
  z <- 3 - y
  
  # Vérification des conditions de non-négativité (x, y, z >= 0)
  if (x >= 0 && z >= 0) {
    return(c(x, y, z))  # Retourne le point (x, y, z) si valide
  } else {
    return(NULL)  # Retourne NULL si les contraintes ne sont pas respectées
  }
}

# Générons les valeurs possibles pour y.
# y doit être compris entre 0 et 1 (limité par x >= 0 et z >= 0).
y_values_cas1 <- seq(0, 1, by = 1)  # On teste deux valeurs : y = 0 et y = 1

# Calculons les sommets pour le Cas 1 à l'aide de calc_sommet_cas1
# lapply applique la fonction à chaque valeur de y_values_cas1, et do.call(rbind, ...) combine les résultats dans un tableau.
sommets_cas1 <- do.call(rbind, lapply(y_values_cas1, calc_sommet_cas1))

# Cas 2 : Intersections avec les axes
# Ces points représentent les cas où certaines variables sont nulles.

# Intersection avec l'axe z (x = 0 et y = 0)
# Si x = 0 et y = 0, la contrainte x + 2y + z ≤ 4 donne z = 4.
sommet_cas2_1 <- c(0, 0, 4)

# Intersection avec l'axe x (y = 0 et z = 0)
# Si y = 0 et z = 0, la contrainte x + y ≥ 1 donne x = 1.
sommet_cas2_2 <- c(1, 0, 0)

# Intersection avec l'axe y (x = 0 et z = 0)
# Si x = 0 et z = 0, la contrainte x + y ≥ 1 donne y = 1.
sommet_cas2_3 <- c(0, 1, 0)

## Étape 3 : Combinaison de tous les sommets

# Nous regroupons les sommets calculés pour les différents cas dans un tableau unique pour faciliter les calculs ultérieurs.
sommets <- rbind(
  sommets_cas1,           # Sommets calculés pour le Cas 1
  sommet_cas2_1,          # Intersection avec l'axe z
  sommet_cas2_2,          # Intersection avec l'axe x
  sommet_cas2_3           # Intersection avec l'axe y
)

# Donnons des noms explicites aux colonnes
colnames(sommets) <- c("x", "y", "z")

# Transformons les sommets en un tableau de données pour ajouter des colonnes et manipuler plus facilement les données.
sommets <- as.data.frame(sommets)

## Étape 4 : Calculer f(x, y, z) pour chaque sommet

# Nous utilisons la fonction 'f' définie plus haut pour évaluer la valeur de la fonction objectif pour chaque sommet.

# La fonction 'mapply' applique 'f' ligne par ligne en prenant les colonnes x, y, et z comme arguments.
sommets$f_val <- mapply(f, sommets$x, sommets$y, sommets$z)

## Étape 5 : Identifier le sommet où la fonction est minimale

# La fonction 'which.min' retourne l'indice du sommet où la valeur de la fonction objectif (f_val) est minimale.
min_index <- which.min(sommets$f_val)

# Récupérons les coordonnées du sommet optimal
sommet_optimal <- sommets[min_index, ]

## Résultats

# Affichons tous les sommets calculés avec leurs valeurs de la fonction objectif
print(sommets)

# Affichons le sommet optimal et sa valeur minimale
print(sommet_optimal)

############
############
############

### PARTIE 4. SOUS-GRAPHE RECOUVRANT OPTIMAL DU METRO PARISIEN

# Chargement les données
metro_data <- read.csv("metro_distance_matrix_updated.csv", row.names = 1, check.names = FALSE)

# Convertion de la matrice en un data frame d'arêtes pondérées
edges <- data.frame()
stations <- rownames(metro_data)

for (i in 1:(nrow(metro_data) - 1)) {
  for (j in (i + 1):ncol(metro_data)) {
    edges <- rbind(edges, c(stations[i], stations[j], metro_data[i, j]))
  }
}

colnames(edges) <- c("Station1", "Station2", "Weight")
edges$Weight <- as.numeric(edges$Weight)

# Trie des arêtes par poids
edges <- edges[order(edges$Weight), ]

# Fonction Union-Find pour Kruskal
find <- function(parent, station) {
  if (parent[station] != station) {
    parent[station] <- find(parent, parent[station])
  }
  return(parent[station])
}

union <- function(parent, rank, station1, station2) {
  root1 <- find(parent, station1)
  root2 <- find(parent, station2)
  
  if (root1 != root2) {
    if (rank[root1] > rank[root2]) {
      parent[root2] <- root1
    } else if (rank[root1] < rank[root2]) {
      parent[root1] <- root2
    } else {
      parent[root2] <- root1
      rank[root1] <- rank[root1] + 1
    }
  }
}

# Implémentation de Kruskal
kruskal <- function(edges, stations) {
  parent <- setNames(stations, stations)
  rank <- setNames(rep(0, length(stations)), stations)
  mst <- data.frame(Station1 = character(), Station2 = character(), Weight = numeric())
  total_weight <- 0
  
  for (i in 1:nrow(edges)) {
    station1 <- edges$Station1[i]
    station2 <- edges$Station2[i]
    weight <- edges$Weight[i]
    
    if (find(parent, station1) != find(parent, station2)) {
      union(parent, rank, station1, station2)
      mst <- rbind(mst, edges[i, ])
      total_weight <- total_weight + weight
    }
  }
  
  return(list(mst = mst, total_weight = total_weight))
}

# Application de Kruskal
result <- kruskal(edges, stations)
mst <- result$mst
total_weight <- result$total_weight

# Analyse des degrés des stations dans le MST
degrees <- table(c(mst$Station1, mst$Station2))
importance <- degrees[c("La Défense (Grande Arche)", "Jussieu")]

# Affichage les résultats
print("Arêtes du MST :")
print(head(mst))
cat("\nPoids total du MST :", total_weight, "\n")
cat("\nImportance des stations :\n")
print(importance)


#Ensemble de fonctions permettant de tester si on peut fusionner deux noeuds du graphe.



#' Détermine les zones de z2 reouvrant plusieurs zones de z1
#'
#' Cette fonction permet de déterminer les colonnes de \code{m_crois}
#' correspondant aux zones de z2 à la frontière entre les zones i et j
#' de z1.
#'
#' @param i,j Entiers indiquant les zones de z1.
#' @param m_crois Matrice de croisement.
#' @param type Character entre "all", "multi" et "unique", pour
#' sélectionner le type de zones de z2 à prendre en compte.
#'
#' @return En sortie on a un vecteur d'entiers indiquant les numéros
#' de colonnes de \code{m_crois} impliqués.
zones_front <- function(i, j, m_crois, type = "all"){
  sum_col_i_j <- colSums(m_crois[c(i,j), ,drop = FALSE]) #somme des colonnes pour les lignes i et j
  col_sel_all <- which(0 < m_crois[i,] & m_crois[i,] < sum_col_i_j) #colonnes correspond aux zonages 2 communs à i et à j
  
  if(type == "all" | length(col_sel_all) == 0) return(col_sel_all)
  else{
    col_sel_unique <- which(colSums(m_crois[, col_sel_all, drop = FALSE]) == colSums(m_crois[c(i, j), col_sel_all, drop = FALSE]))
    col_sel_multi <- col_sel_all[! col_sel_all %in% col_sel_unique]
    
    if(type == "unique") return(col_sel_unique)
    if(type == "multi") return(col_sel_multi)
  }
}

#' Valeurs des arêtes reliant i et j
#'
#' i et j étant deux zones de z1, cette fonction donne la valeur
#' de l'arête reliant i et j et de celle reliant j à i.
#'
#' @param i,j Entiers indiquant les zones de z1.
#' @param m_crois Matrice de croisement.
#' @param ... d'autres paramètres hérités de \code{zones_front}.
#'
#' @return En sortie on obtient un vecteur de deux entiers
#' donnant la valeur de l'arête de i vers j et de celle de j
#' vers i.
valeurs_lien <- function(i, j, m_crois, ...){
  col_sel <- zones_front(i, j, m_crois, ...)
  if(length(col_sel) >= 2){
    res <- rowSums(m_crois[c(i,j), col_sel])
  }
  else if (length(col_sel) == 1){
    res <- m_crois[c(i, j), col_sel]
  }
  else{
    res <- c(0,0)
    warning("Look out, there is no link between i and j")
  }
  
  return(res)
}


#' Tester, selon la méthode 1, si on peut fusionner deux sommets
#'
#' Fonction permettant de tester si les sommets i et j peuvent
#' être fusionnés. La méthode 1 consiste à regarder les liens
#' entre les deux sommets et voir si ces deux liens ont une valeur
#' plus grande que le threshold de confidentialité. Si c'est le cas, on
#' peut fusionner les deux sommets.
#'
#' La fonction renvoie TRUE si on peut fusionner les deux sommets (
#' c'est-à-dire qu'il n'y a aucun problème de différenciation en considérant
#' un sommet sans considérer l'autre), et FALSE sinon.
#'
#' @param i,j Entiers indiquant deux lignes différentes de \code{m_crois}.
#' @param m_crois Matrice de croisement.
#' @param threshold threshold de confidentialité.
#'
#' @return En sortie, on a une liste contenant deux élements, fus et col.
#' fus vaut TRUE ou FALSE et indique s'il faut fusionner les deux lignes.
#' col est un vecteur d'entier qui indique les numéros de colonne à supprimer.
test_fus_m1 <- function(i, j, m_crois,threshold){
  sum_col_i_j <- colSums(m_crois[c(i,j), ,drop = FALSE]) #somme des colonnes pour les lignes i et j
  col_sel <- which(0 < m_crois[i,] & m_crois[i,] < sum_col_i_j) #colonnes correspond aux zonages 2 communs à i et à j
  
  if(length(col_sel) >= 2) { #les deux zones i et j partagent plusieurs zones de z2 en commun
    res <- sum(rowSums(m_crois[c(i,j), col_sel]) >= threshold) == 2 #res vaut TRUE si les deux liens sont au-dessus du threshold et FALSE sinon
    col_to_suppress <- col_sel[sum_col_i_j[col_sel] == colSums(m_crois[, col_sel])] #colonnes à retirer en cas de fusion (étape suivante)
    #ce sont les zones de z2 qui ne sont communes qu'à i et j (et pas d'autre zone de z1)
  }
  else if(length(col_sel) == 1) { #les deux zones i et j partagent une seule zone de z2 en commun
    res <- sum(m_crois[c(i,j), col_sel] >= threshold) == 2 #vaut TRUE ou FALSE
    col_to_suppress <- col_sel[sum_col_i_j[col_sel] == sum(m_crois[, col_sel])] #colonnes à retirer en cas de fusion (étape suivante)
  }
  else {
    warning(paste0(" zones ",i," and ",j," are not contiguous"))
    res <- FALSE
    col_to_suppress <- NULL
  }
  return(list(fus = res, col = col_to_suppress))
}


#' Tester s'il existe un chemin, selon la méthode 2
#'
#' Fonction permettant de tester s'il existe un chemin
#' entre les noeuds i et j où chaque arête a une valeur
#' plus grande que le threshold de confidentialité moins
#' la valeur du lien entre i et j.
#'
#' @param i,j Entiers indiquant deux sommets différents du graphe.
#' @param m_graph Matrice carré d'adjacence. Si on note a_ij l'élément
#' de cette matrice correspondant à la ième ligne et jème colonne, alors
#' a_ij=0 si les sommets i et j ne sont pas connectés et sinon a_ij est un entier
#' qui indique la valeur de l'arête entre i et j.
#' @param threshold threshold de confidentialité.
#' @param v_arete valeur du lien entre i et j
#'
#' @return Un booléen qui vaut TRUE si les deux sommets i et j sont
#' connectés, selon la méthode 2, et FALSE sinon.
is_connected_m2 <- function(i, j, m_graph, v_arete, threshold){
  
  if(v_arete >= threshold) return(TRUE)
  
  #______________________________________________
  #Transformation de la matrices d'adjacence:
  #______________________________________________
  #On met à 0 l'arête entre i et j
  m_graph[i, j] <- 0
  
  #On met à 0 les arêtes de valeur inférieur strictement au threshold - v
  #m_graph[m_graph < threshold - v]=0
  m_graph <- m_graph * (m_graph >= threshold - v_arete) #autre possibilité
  
  #______________________________________________
  #Tester s'il existe un chemin, dans cette nouvelle matrice
  #entre i et j
  #______________________________________________
  g <- igraph::graph_from_adjacency_matrix(m_graph, mode = "directed", weighted = TRUE)
  connected <- j %in% igraph::bfs(g, root = i,neimode = "out", unreachable = FALSE)$order
  
  return(connected)
}


#' Tester, selon la méthode 2, si on peut fusionner deux noeuds ou non.
#'
#' Permet de tester si les deux sommets i et j du graphe peuvent
#' être fusionnés.
#'
#' La méthode 2 s'applique après la méthode 1. Si les sommets i et j
#' sont reliés par une arête en-dessous du threshold de confidentialité
#' la méthode 1 conclut qu'on ne peut pas les fusionner. Mais la méthode
#' 2 regarde s'il existe un autre chemin entre i et j, avec des arêtes
#' ayant des valeurs suffisamment élevées (plus grande que le threshold - la valeur
#' du lien entre i et j). Si c'est le cas, on peut en déduire qu'il n'y
#' aura pas de problème de différentiation en considérant i et j séparément.
#' On peut donc fusionner ces deux zones.
#'
#' @param i,j Entiers indiquant les noeuds du graphe à tester.
#' @param m_crois Matrice de croisement.
#' @param threshold threshold de confidentialité.
#'
#' @return En sortie, on a une liste contenant deux élements, fus et col.
#' fus vaut TRUE ou FALSE et indique s'il faut fusionner les deux lignes.
#' col est un vecteur d'entier qui indique les numéros de colonne à supprimer.
test_fus_m2 <- function(i, j, m_crois, threshold){
  
  # Enlever de m_crois les zones de z2 multi entre i et j:
  col_multi <- zones_front(i, j, m_crois, type = "multi")
  m_graph <- matrix_graphe(m_crois[, -col_multi, drop = FALSE], multi = TRUE)
  v_aretes <- valeurs_lien(i, j, m_crois, type = "all")
  
  if(sum(v_aretes == 0) >= 1){
    res <- FALSE
    col_to_suppress <- integer(0)
  }
  
  else{
    
    res = is_connected_m2(i, j, m_graph, v_aretes[1], threshold) & is_connected_m2(j, i, m_graph, v_aretes[2], threshold)
    
    #______________________________________________
    #sélection des colonnes de m_crois à supprimer
    #(i.e. les carreaux à l'intérieur de la fusion i-j)
    #______________________________________________
    sum_col_i_j <- colSums(m_crois[c(i,j),, drop = FALSE]) #somme des colonnes pour les lignes i et j
    col_sel = 0 < m_crois[i,] & m_crois[i,] < sum_col_i_j #colonnes correspond aux zonages 2 communs à i et à j
    
    if(sum(col_sel) >= 2) {
      col_to_suppress <- which(col_sel)[sum_col_i_j[col_sel] == colSums(m_crois[,col_sel])] #colonnes à retirer en cas de fusion (étape suivante)
    }
    else if(sum(col_sel) == 1) {
      col_to_suppress <- which(col_sel)[sum_col_i_j[col_sel] == sum(m_crois[,col_sel])] #colonnes à retirer en cas de fusion (étape suivante)
    }
    else {
      warning(paste0("Zones ",i," and ",j," are not contiguous"))
      col_to_suppress <- integer(0)
    }
  }
  
  return(list(fus = res, col = col_to_suppress))
  
}


#' Tester si on peut fusionner deux zones de z1
#'
#' On teste selon la méthode 1, 2 ou les deux à la fois, si
#' deux zones i et j de z1 peuvent être fusionnées ou non
#'
#' @param i,j Entiers indiquant les zones de z1.
#' @param m_crois Matrice de croisement.
#' @param threshold threshold de confidentialité.
#' @param methode Méthode à choisir parmi "m1", "m2" et "both".
#'
#' @return En sortie, on a une liste contenant deux élements, fus et col.
#' fus vaut TRUE ou FALSE et indique s'il faut fusionner les deux lignes.
#' col est un vecteur d'entier qui indique les numéros de colonne à supprimer.
test_fus=function(i, j, m_crois, threshold, methode = "both"){
  if(methode == "m1") test_f <- test_fus_m1(i, j, m_crois, threshold)
  if(methode == "m2") test_f <- test_fus_m2(i, j, m_crois, threshold)
  if(methode == "both") {
    test_m1 <- test_fus_m1(i, j, m_crois, threshold)
    if(test_m1$fus == TRUE) test_f <- test_m1 #si la méthode 1 conclut à une fusion, pas besoin d'utiliser la méthode 2
    else test_f <- test_fus_m2(i, j, m_crois, threshold)
  }
  
  return(test_f)
}

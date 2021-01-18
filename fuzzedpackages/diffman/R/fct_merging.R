#Ensemble de fonctions pour realiser l'agregation du graphe (ou la reduction du nombre
#de zones de z1 par fusion de zones)




#' Fonction qui fusionne les zones numero i et numero j
#'
#' Fusion des zones i et j du zonage z1.
#'
#' Les colonnes a supprimer \code{col_to_suppress} sont determinees
#' lors du test de fusion. Si le test conclue qu'il faut fusionner
#' deux zones, alors on retourne en plus les zones de z2 qui sont
#' a l'interieur de cette zone fusionnee.
#'
#' @param i,j Entiers indiquant les zonages de z1.
#' @param m_crois Matrice de croisement.
#' @param col_to_suppress Vecteur d'entiers indiquant les colonnes a supprimer.
#'
#' @return En sortie, on a la matrice de croisement avec une ligne en moins,
#' et potentiellement des colonnes en moins egalement.
fusion <- function(i, j, m_crois, col_to_suppress){
  #on modifie le nom de la ligne i:
  rownames(m_crois)[i] <- paste0(c(rownames(m_crois)[i], rownames(m_crois)[j]), collapse = ".")
  
  
  if(length(col_to_suppress) %in% c(ncol(m_crois), ncol(m_crois) - 1)) {
    print("Everything is merged")
  }
  
  else if(length(col_to_suppress) > 0) {
    m_crois <- m_crois[, -col_to_suppress, drop = FALSE] #on supprime les colonnes
  }
  
  m_crois[i, ] <- m_crois[i, ] + m_crois[j, ] #on additionne les 2 lignes
  m_crois <- m_crois[-j, , drop = FALSE] #on retire la ligne j
  
  #On supprimer les colonnes avec un seul chiffre
  m_crois <- simplify_mcrois(m_crois)
  
  return(m_crois)
}



#' Fonction pour agreger une seule fois 2 zones de z1
#'
#' On test toutes les paires de zones de z1 jusqu'a trouver
#' deux zones qu'on peut fusionner. Alors on fusionne et on s'arrête la.
#'
#' @param m_crois Matrice de croisement.
#' @param threshold Entier indiquant le threshold de confidentialite.
#' @param shuffle Booleen indiquant s'il faut melanger les lignes
#' de \code{m_crois} apres avoir fusionne ou non. Vaut TRUE par defaut.
#' @param verbose Boolean. If TRUE, displays messages.
#' @param ... other parameters
#'
#' @return En sortie on a une liste de deux elements : \code{matrice} donne
#' la matrice de croisement avec une ligne en moins et \code{continue} indique
#' si l'agregation peut être poursuivie.
agregate_one <- function(m_crois, threshold, shuffle = TRUE, verbose = TRUE, ...){
  m_liens <- m_crois %*% t(m_crois) > 0 ; #on calcule la matrice de contiguite
  
  #===============================
  #Recherche d'une paire (i,j) a fusionner
  #===============================
  i <- 1 #indice de la zone etudiee
  continue_search <- TRUE
  while(i <= nrow(m_liens) & continue_search){
    vois <- which(m_liens[i, ]) #voisins de i dans le graphe (i inclus)
    vois <- vois[ -which(vois == i)] #voisins de i, sans i
    
    k <- 1 #indice du voisin
    while(k <= length(vois) & continue_search){
      test_res <- test_fus(i, vois[k], m_crois, threshold, ...)
      continue_search <- !test_res$fus  #on arrête de chercher en cas de fusion entre i et vois[k]
      k <- k + 1
    }
    i <- i + 1
  }
  
  
  #===============================
  #Fusion de la paire (i,j)
  #===============================
  if(!continue_search){
    m_crois <- fusion(i - 1, vois[k - 1], m_crois, test_res$col) #on fusionne
    if(shuffle) perm <- sample(1:nrow(m_crois)) else perm <- 1:nrow(m_crois) #pour melanger les lignes
    return(list(matrice = m_crois[perm, , drop = FALSE], continue = TRUE))
  }
  
  #===============================
  #Plus aucune fusion possible
  #===============================
  else {#dans ce cas, on est arrive a la fin des lignes de m_crois sans avoir trouve de paires a fusionner,
    #c'est donc qu'on ne peut plus agreger le graphe d'avantage.
    if(verbose) message("End of merging")
    return(list(matrice = m_crois, continue = FALSE))
  }
}


#' Fonction pour agreger le graphe
#'
#' On agrege le graphe en utilisant la fonction agregate_one un
#' grand nombre de fois, jusqu'a ce que plus aucune agregation ne soit possible
#' où jusqu'a un threshold limite d'agregations.
#'
#' @param m_crois Matrice de croisement.
#' @param threshold threshold de confidentialite.
#' @param kmax Entier indiquant le nombre d'agregtion maximale a effectuer.
#' Si \code{kmax} vaut 0, on agrege au maximum.
#' @param pas_pb Entier indiquant a quelle frequence afficher la barre de progression
#' (pb = progress bar). Par defaut, affiche tout les 2\% environ.
#' @param verbose Boolean. If TRUE, progress bar is displayed.
#' @param ... other parameters
#'
#' @return En sortie, on a la matrice de croisement apres agregation.
agregate <- function(m_crois, threshold, kmax = 0, pas_pb = floor(nrow(m_crois) / 50), verbose = TRUE,...){
  if(pas_pb == 0) pas_pb <- 1
  if(kmax == 0) kmax <- nrow(m_crois) #on veut tout agreger
  k <- 1 #permet de compter le nombre d'agregation realises
  continue <- TRUE #indique si on continue ou non a vouloir agreger le graphe
  
  if(verbose) message("     Start of the merging process")
  
  #Barre de progression (package progress)
  if(verbose){
    pb <- progress::progress_bar$new(
    format = "  merging [:bar] :percent in :elapsed , remain :eta",
    total = kmax, clear = FALSE, width= 50)
  
    pb$tick(0)
  }
  
  while(continue & k < kmax){
    res_ag <- agregate_one(m_crois, threshold, verbose = verbose, ...) #on agrege deux zones
    m_crois <- res_ag$matrice
    continue <- res_ag$continue
    if(k %% pas_pb == 0 & continue & verbose) pb$tick(pas_pb) #print(paste0(round(k * 100 / kmax, 1), "%"))
    k <- k+1
  }
  
  if(! is.null(dim(m_crois))) m_crois <- simplify_mcrois(m_crois)
  
  return(m_crois)
}
















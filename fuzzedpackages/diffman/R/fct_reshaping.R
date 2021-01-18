#Ce script rassemble les fonctions permettant de mettre en forme les données:
#



#' Enlever les zones de z2 entièrement incluses dans une zone de z1
#'
#' Cette fonction permet de retirer de la table individuelle, les observations
#' contenues dans une zone de z2 entièrement incluse dans une seule zone de z1.
#' Une zone de z2 est entièrement incluse dans une autre zone si toutes les
#' observations qu'elle contient sont également contenue dans cette autre zone.
#'
#' @param t_ind La table individuelle (format data.table) Chaque ligne représente une observation
#' et elle doit au moins contenir deux colonnes nommées 'z1' et 'z2' et indiquant
#' à quels zonages appartient l'observation.
#'
#' @return En sortie on récupère une table individuelle avec des observations en
#' moins
simplify_z2_rem <- function(t_ind){
  #define variables to avoid NOTE of "no visible binding for global
  # variable" when running R CMD check
  z1 = z2 = nb_z1 = NULL
  
  t_z2 <- t_ind[,.(nb_z1=sum(!duplicated(z1))),by=.(z2)] #t_z2 est une table ayant autant de ligne
  #qu'il y a de zones z2. Pour chaque zone de z2 elle indique le nombre de zones de z1, nb_z1, que la zone
  #z2 recouvre (au sens des observations et non de la surface).
  
  t_ind <- t_ind[z2 %in% t_z2[nb_z1>=2]$z2] #on ne garde que les observations situées dans des zones
  #de z2 "à la frontière" c'està-dire recouvrant plusieurs zones de z1
  
  return(t_ind)
}


#' Creer la matrice de croisement entre deux zonages
#' 
#' Compte le nombre d'observations aux intersections des 
#' deux zonages.
#' 
#' @param t_ind La table individuelle (format data.table) Chaque ligne représente une observation
#' et elle doit au moins contenir deux colonnes nommées 'z1' et 'z2' et indiquant
#' à quels zonages appartient l'observation.
#' 
#' @return table de croisement
tab_crois <- function(t_ind){
  #define variables to avoid NOTE of "no visible binding for global
  # variable" when running R CMD check
  z1 = z2 = NULL
  
  # on compte le nombre d'observations par croisement entre les zonages z1 et z2
  
  t_crois <- t_ind[,.(nb_obs=.N),by=.(z1=factor(z1),z2=factor(z2))]
  return(t_crois)
}


#' Fusionne les zones de z2 recouvrant les mêmes zones de z1
#'
#' Cette fonction permet de fusionner les zones de z2 dont les observations
#' sont réparties sur les mêmes zones de z1. Cela permet de diminuer le
#' nombre de zonages de z2.
#'
#' @param t_crois Table de croisement (c'est la table de fréquences lorsqu'on
#' croise les deux zonages z1 et z2).
#'
#' @return En sortie on obtient une table de croisement avec moins de classes
#' pour le zonage z2.
simplify_z2_fus <- function(t_crois){
  #define variables to avoid NOTE of "no visible binding for global
  # variable" when running R CMD check
  z1 = z2 = z2_b = nb_obs = NULL
  
  t_crois <- t_crois[,z2_b:=paste0(sort(z1),collapse="-"),by=.(z2)][ #z2_b est une "étiquette" qui donne les zones de z1 recouvertes par chaque zone de z2
    ,.(nb=sum(nb_obs)),by=.(z1,z2_b)][
      ,z2_b:=factor(z2_b)]
  
  #Remarque : les noms des zones de z2 ont changé avec l'application de cette fonction.
  #Avant, les noms étaient ceux du zonage z2 (noms des carreaux par exemple). Après
  #les noms sont basés sur les noms de z1 : une zone de z2 est nommée selon les zones
  #de z1 qu'elle recouvre.
  
  #on change les noms de colonnes
  colnames(t_crois)[colnames(t_crois)=="z2_b"] <- "z2"
  colnames(t_crois)[colnames(t_crois)=="nb"] <- "nb_obs"
  
  return(t_crois)
}


#' Obtenir une matrice (sparse) qui croise les zones de z1 et celles de z2
#'
#' La matrice de croisement est une matrice dont les lignes indiquent les zones
#' du zoange z1 et les colonnes les zones du zonage z2. Chaque élément de la
#' matrice donne le nombre d'observations situées à l'intersection entre une
#' zone de z1 et une zone de z2.
#'
#' @param t_crois Table de croisement (1 colonne z1 et une colonne z2).
#'
#' @return En sortie on obtient une matrice sparse de croisement. Les noms
#' des lignes correspondent aux noms des zones de z1 et les noms des colonnes
#' aux noms des zones de z2.
matrix_crois <- function(t_crois){
  #define variables to avoid NOTE of "no visible binding for global
  # variable" when running R CMD check
  z1 = z2 = NULL
  
  # On recalcul les facteurs, car des modalités ont pu disparaitre
  t_crois[, ":="(z1 = factor(z1), z2 = factor(z2))]
  
  # Matrice z1Xz2
  m_crois <- Matrix::sparseMatrix(i=as.numeric(t_crois$z1),
                                  j=as.numeric(t_crois$z2),
                                  x=t_crois$nb_obs,
                                  dimnames=list(levels(t_crois$z1),levels(t_crois$z2)))
  return(m_crois)
}


#' Creation de la matrice de croisement
#' 
#' Permet de creer directement la matrice
#' de croisement a partir de la table individuelle
#' des observations. Wrapper pour les fonctions
#' simplify_z2_rem, tab_crois, simplify_z2_fus,
#' matrix_crois et simplify_mcrois
#' 
#' @param t_ind La table individuelle (format data.table) Chaque ligne représente une observation
#' et elle doit au moins contenir deux colonnes nommées 'z1' et 'z2' et indiquant
#' à quels zonages appartient l'observation.
#' 
#' @return matrice sparse de croisement
matrix_crois_from_tind <- function(t_ind){
  t_ind <- as.data.table(t_ind)
  t_ind <- simplify_z2_rem(t_ind)
  t_crois <- tab_crois(t_ind)
  t_crois <- simplify_z2_fus(t_crois)
  m_crois <- matrix_crois(t_crois)
  m_crois <- simplify_mcrois(m_crois)
  
  return(m_crois)
}

#' Simplifier la matrice de croisement
#'
#' Simplifier la matrice de croisement revient à supprimer deux types de colones
#' (et donc supprimer des zones du zonage z2).
#' 1 - on supprime les colonnes vides, c'est-à-dire ne contenant que des 0.
#' 2 - on supprime les colonnes ne contenant qu'un seul élément différent de 0.
#'
#' @param m_crois Une matrice de croisement.
#'
#' @return En matrice on a une matrice de croisement simplifiée, c'est-à-dire
#' avec moins de colonnes.
simplify_mcrois <- function(m_crois){
  # col_to_delete_1 <- which(colSums(m_crois>0) == 1)
  # col_to_delete_2 <- which(colSums(m_crois>0) == 0)
  # col_to_delete <- c(col_to_delete_1, col_to_delete_2)
  
  #nouvelle façon:
  col_to_delete <- which(colSums(m_crois>0) <= 1)
  
  if(length(col_to_delete)>0)
    m_crois <- m_crois[, -col_to_delete, drop = FALSE]
  
  return(m_crois)
}



#' Créer la matrice de liens (ou matrice de contiguïté)
#'
#' Cette fonction permet de créer une matrice carré de booléens de taille
#' égale au nombre de zones du zonage z1.
#'
#' Un élément (i,j) de cette matrice
#' vaut TRUE (ou 1) si les zones i et j du zonage z1 sont contigues, c'est-à-dire
#' s'il existe au moins une zone de z2 recouvrant à la fois i et j. Les élements
#' de la diagonales portent la valeur FALSE.
#'
#' @param m_crois Matrice de croisement.
#'
#' @return En sortie on a une matrice carré de booléens.
matrix_liens <- function(m_crois){
  
  m_liens <- m_crois %*% t(m_crois) > 0
  diag(m_liens) <- FALSE
  
  return(m_liens)
}



#' Constuire la matrice d'adjacence du graphe
#'
#' Fonction permettant à partir de la matrice de croisement
#' de déterminer la matrice d'adjacence du graphe. Cette matrice
#' de graphe est pondérée et non symmétrique (ce qui correspond
#' à un graphe orienté).
#'
#' L'option \code{multi} permet de choisir si on prend en compte
#' ou non les zones de z2 recouvrant 3 zones de z1 ou plus. En effet
#' si on les prend en compte, alors certaines observations sont comptées
#' plusieurs fois dans le graphe, ce qui peut conduire à de mauvaises
#' interprétations.
#'
#' @param m_crois Matrice de croisement.
#' @param multi Booléen indiquant s'il faut considérer les zones de z2
#' recouvrant trois zones ou plus de z1.
#'
#' @return En sortie on obtient une matrice carré d'adjacence.
matrix_graphe <- function(m_crois, multi = TRUE){
  if(multi == FALSE) { #on enlève les carreaux sur 3 communes ou plus
    car_sel <- colSums(m_crois > 0) >= 3
    m_crois <- m_crois[, car_sel, drop = FALSE]
  }
  m_graph <- m_crois %*% (t(m_crois) > 0)
  diag(m_graph) <- 0
  return(m_graph)
}



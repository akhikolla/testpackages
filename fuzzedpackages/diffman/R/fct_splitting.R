#Ensemble de fonctions pour determiner les composantes connexes a tester exhaustivement.


#' Detecter un noeud central
#'
#' Un noeud central est un noeud qui lorsqu'on le retire
#' rend le graphe non connexe. Cette fonction permet de
#' detecter un noeud parmi les noeud centraux.
#'
#' @param m_crois Matrice de croisement.
#'
#' @return En sortie on a un entier qui indique le numero
#' du noeud central (c'est-a-dire le numero de la ligne).
#'Si aucun noeud central n'a ete trouve on retourne 0.
detect_central_node <- function(m_crois){
  
  m_liens <- matrix_liens(m_crois)
  g <- igraph::graph_from_adjacency_matrix(m_liens) #graphe simple (non oriente, non pondere) des liens entre zones de z1
  comp <- igraph::components(g) #les composantes connexes du graphe
  
  if(comp$no > 1) stop("Attention, la matrice de croisement donnee n'est pas connexe")
  
  # Classer les noeuds du graphe selon leur nombre de voisins.
  nb_vois <- rowSums(m_liens)
  ind_nd <- rev(order(nb_vois)) #indices des noeuds par ordre decroissant de nombre de voisins
  
  # On cherche un noeud qui lorsque retire, laisse le graphe en plus de trois composantes connexes
  i <- 1
  nd_res <- 0
  continue <- nb_vois[ind_nd[i]] >= 3 #on continue a chercher si le noeud en question a au moins 3 voisins
  i_max <- sum(nb_vois >= 3) #nombre de noeuds ayant au moins 3 voisins
  while(continue & i <= i_max){
    nd <- ind_nd[i] #numero du noeud
    m_liens_div <- m_liens[ -nd, -nd, drop = FALSE]
    g_div <- igraph::graph_from_adjacency_matrix(m_liens_div)
    nb_cc <- igraph::components(g_div)$no
    
    if(nb_cc >=3) {continue <- FALSE ; nd_res <- nd}
    else i <- i+1
  }
  
  return(nd_res)
}



#' Desagreger en composantes connexes
#'
#' Permet de separer la matrice de croisement en plusieurs
#' sous-matrices de croisement connexes.
#'
#' @param m_crois Matrice de croisement.
#'
#' @return On retourne une liste de matrices de croisement etant
#' chacune connexe.
comp_connexe_list <- function(m_crois){
  
  g <- igraph::graph_from_adjacency_matrix(matrix_liens(m_crois)) #graphe simple (non oriente, non pondere) des liens entre zones de z1
  comp <- igraph::components(g) #les composantes connexes du graphe
  #gpe <- groups(comp)
  
  i <- 1
  l_mat <- vector("list", sum(comp$csize > 1)) #liste des sous-matrices dont la representation sous forme de graphe est connexe.
  #on va completer cette liste au fur et a mesure.
  for (g in 1:comp$no){ # g indique le numero de la composante connexe en cours
    if(comp$csize[g] > 1){ #si la composante connexe est composee de plus de deux elements
      m_gpe <- m_crois[comp$membership == g,, drop = FALSE]
      m_gpe <- simplify_mcrois(m_gpe) # m_gpe est la matrice de croisement correspondant a la composante connexe g
      l_mat[[i]] <- m_gpe
      i <- i+1
    }
    
    #si comp$csize[g]==1, c'est-a-dire qu'il n'y a qu'un seul sommet sur la composante connexe
    #du graphe, alors il n'y a pas de probleme de differenciation sur ce sommet, et on ne l'ajoute
    #pas a la liste l_mat
  }
  
  return(l_mat)
}


#' Casser un graphe en plusieurs sous-graphes
#'
#' Cette fonction permet de decomposer un graphe
#' connexe en plusieurs sous-graphes connnexe
#' apres avoir detecter un noeud central.
#'
#' @param m_crois_connexe Matrice de croisement connexe.
#'
#' @return On retourne une liste de liste. Chaque sous-liste
#' est composee d'une matrice \code{m} et d'un \code{etat} parmi
#' "cassable" ou "incassable".
matrix_break <- function(m_crois_connexe){
  nd_central <- detect_central_node(m_crois_connexe)
  if(nd_central == 0){
    res <- list(list(m = m_crois_connexe, etat = "incassable"))
  }
  else{
    m <- m_crois_connexe[ -nd_central, , drop = FALSE]
    m_nd <- m_crois_connexe[nd_central,, drop = FALSE]
    g <- igraph::graph_from_adjacency_matrix(matrix_liens(m))
    comp <- igraph::components(g)
    
    compteur <- 1
    res <- list()
    for(k in 1:comp$no){
      m_k <- m[comp$membership == k, , drop = FALSE]
      rownames(m_nd) <- paste0(setdiff(rownames(m_crois_connexe),rownames(m_k)),collapse = ".")
      m_k <- rbind(m_k, m_nd)
      m_k <- simplify_mcrois(m_k)
      
      res[[compteur]] <- list(m = m_k, etat = "cassable")
      compteur <- compteur +1
    }
  }
  
  return(res)
}



#' Completement decomposer un graphe.
#'
#' Permet de decomposer un graphe, de façon iterative, jusqu'a
#' ce que tous les sous-graphes soient "incassables".
#'
#' @param m_crois_connexe Matrice de croisement connexe.
#' @param etat Caracteres parmi "cassable" ou "incassable" indiquant
#' l'etat de la matrice de croisement connexe.
#' @param taille_max Entier indiquant a partir de quelle taille de graphe
#' on decide de le decomposer (le "casser").
#'
#' @return En sortie on a une liste de matrices de croisement connexes soient
#' "incassables" soient ayant moins de lignes que \code{taille_max}.
matrix_atomize <- function(m_crois_connexe, etat = "cassable", taille_max){
  if(etat == "incassable" | nrow(m_crois_connexe) <= taille_max) l_res <- list(m_crois_connexe)
  else{
    l_res <- list()
    l_mat <- matrix_break(m_crois_connexe)
    for(i in 1:length(l_mat)){
      l_atom <- matrix_atomize(l_mat[[i]]$m, l_mat[[i]]$etat, taille_max)
      l_res <- c(l_res, l_atom)
    }
  }
  
  return(l_res)
}



#' Decomposer une matrice de croisement en composantes connexes
#'
#' Permet de decomposer totalement une matrice de croisement, qu'elle
#' soit connexe ou non, en plusieurs sous-matrices de croisements
#' connexes "incassables" ou ayant moins de lignes que \code{taille_max}
#'
#' @param m_crois Matrice de croisement.
#' @param taille_max Entier indiquant la taille a partir de laquelle
#' on decide de "casser" la matrice ou non. Vaut 10 par defaut.
decompose_m_crois <- function(m_crois, taille_max = 10){
  l_mat <- comp_connexe_list(m_crois)
  
  l_res <- list()
  for(i in 1:length(l_mat)){
    l_res <- c(l_res, matrix_atomize(l_mat[[i]], taille_max = taille_max))
  }
  
  return(l_res)
}

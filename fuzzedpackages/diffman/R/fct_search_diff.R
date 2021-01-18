#Ensemble des fonctions pour rechercher les problemes de differenciation


#' Tester toutes les combinaisons possibles (jusqua une certaine taille)
#'
#' Permet de tester, pour une liste de matrice de croisement, les agregats
#' de zones de z1 conduisant a un probleme de differenciation.
#'
#' @param list_m_crois Liste de matrices de croisement.
#' @param threshold threshold de confidentialite.
#' @param max_agregate_size Entier indiquant la taille maximale des agregats
#' a tester
#'
#' @return On retourne une liste d'agregats.
search_diff_agregate=function(list_m_crois, threshold, max_agregate_size = 20){
  
  a <- 1
  list_agregat <- list()
  
  #Barres de progressions :
  pb1 <- progress::progress_bar$new(
    format = "  composantes connexes [:bar] :percent, reste :eta",
    total = length(list_m_crois), clear = FALSE, width= 60,
    complete = "*", incomplete = ".")
  
  pb1$tick(0)
  
  
  for(i in 1:length(list_m_crois)){
    m <- list_m_crois[[i]]
    base <- as.matrix(m)
    if(ncol(base) == 1) colnames(base) = "1carreau"
    
    #### Preparation des donnees ####
    commune0 <- data.frame(z1 = rownames(base), men = rowSums(base), stringsAsFactors = FALSE)
    carreau0 <- data.frame(z2 = colnames(base), men = colSums(base), stringsAsFactors = FALSE)
    rownames(base) <- NULL
    colnames(base) <- NULL
    mat_conti_commune_carreau <- (base > 0) * 1L
    mat_conti_commune <- (mat_conti_commune_carreau %*% t(mat_conti_commune_carreau) > 0) * 1L
    diag(mat_conti_commune) <- 0L
    
    taille_test <- min(nrow(m), max_agregate_size)
    
    #### On test tous les agregats un a un ####
    for(taille_agregat in 1:taille_test){
      
      pb2 <- progress::progress_bar$new(
        format = "  agregats [:bar] :current ",
        total = taille_test, clear = TRUE, width= 50)
      
      #### appel de la detection de differenciation ####
      dfResultat <- differencierRcpp(iTailleCible = taille_agregat, 
                                     iSeuil = threshold,
                                     vNbObsTerritoire = commune0$men,
                                     vNbObsCarreaux = carreau0$men,
                                     mContiguiteT = mat_conti_commune,
                                     mContiguiteTC = mat_conti_commune_carreau)
      
      if(length(dfResultat)>0){
        for(k in 1:length(dfResultat)){
          num <- dfResultat[[k]][1:taille_agregat]
          list_agregat[[a]] <- commune0$z1[num]
          a <- a+1
        }
      }
      
      pb2$tick()
      
    }
    
    pb1$tick()
  }
  return(list_agregat)
}


#' Desagreger la liste d'agregat
#'
#' Permet de creer une liste de vecteur, chaque vecteur donnant
#' la composition en termes de zones de z1 de l'agregat.
#'
#' @param list_agregat Liste de chaînes de caractere. Chaque element de la liste
#' est une chaîne de caracteres donnant les nomes des zones de z1 composant l'agregat.
#' Chaque nom est separe dans la chaîne de caractere par le symbole ".".
#'
#' @return On retourne une liste de vecteur de chaînes de caracteres.
desagregate_list=function(list_agregat){
  lapply(list_agregat,function(x){ unlist(strsplit(x, ".", fixed = TRUE)) })
}


#' Recuperer les identifiants des observations a risque
#'
#' Permet de recuperer les identifiants des observations a
#' risque de differenciation. Permet egalement de donner le nombre
#' d'observations qu'il y a sur les zones "internes" et "externes"
#' lorsqu'on effectue la differenciation.
#'
#' Il se peut que certains agregats conduisent au même probleme de differenciation,
#' c'est-a-dire qu'ils permettent de deduire de l'information sur le même
#' groupe d'observation. Dans ce cas, cette fonction ne garde que l'agregat de plus
#' petite taille (en terme de nombre de zones de z1 impliques dans l'agregat). Il
#' se peut donc qu'en sortie, le nombre d'agregats soit inferieur au nombre d'agregats
#' en entree.
#'
#' @param list_agregat Liste de vecteurs de chaînes de caracteres donnant les noms
#' des zones de z1 constituant les differents agregats. Ces agregats ont dû
#' prealablement être identifies comme conduisant a un probleme de differenciation.
#' @param t_ind La table individuelle donnant pour chaque observation la zone de z1
#' et la zonne de z2 auquelle elle appartient.
#' @param threshold Entier indiquant le threshold de confidentialite.
#' @param verbose Boolean. If TRUE, progress bar is displayed.
#'
#' @return En sortie on obtient un data.frame/data.table donnant la liste
#' des observations a risque, sans doublons, et indiquant pour chaque observation
#' a risque, l'agregat sur lequel la differenciation est faite, le type de differenciation
#' (interne ou externe), la taille de l'agregat (nombre de zones de z1) et le nombre
#' d'observations a risque dans la même differenciation
find_id_obs_risque <- function(list_agregat, t_ind, threshold, verbose = TRUE){
  
  #define variables to avoid NOTE of "no visible binding for global
  # variable" when running R CMD check
  z1 = z2 = J = min_nb_obs = nb_obs = id_obs = min_taille_ag = agregat_size = NULL
  
  t_ind$id <- as.character(t_ind$id)
  
  taille_ag <- unlist(lapply(list_agregat, length)) # pour ne pas compter le "big agregat"
  list_agregat <- list_agregat[which(taille_ag <= length(unique(t_ind$z1)) / 2)] # pour ne pas compter le "big agregat"
  
  # obs_risque <- vector("list", length(list_agregat))
  # nb_risque <- vector("list", length(list_agregat))
  
  t_ind_2 <- copy(t_ind)
  setkey(t_ind, z1)
  setkey(t_ind_2, z2)
  
  taille_ag <- sapply(list_agregat, length) #taille des agregats
  
  
  pb <- progress::progress_bar$new(
    format = " recherche obs a risque [:bar] :current, reste :eta",
    total = length(list_agregat), clear = TRUE, width= 80,
    complete = "*", incomplete = ".")
  
  if(verbose) pb$tick(0)
  
  res <- data.frame()
  
  for(i in 1:length(list_agregat)){
    
    #On cherche les carreaux aux frontieres
    t_ind_in_ag <- t_ind[J(list_agregat[[i]]), nomatch = 0L] #individus dans l'agregat
    z2_i <- unique(t_ind_in_ag$z2) # zones des z2 (carreaux) dans l'agregat ou a la frontiere de l'agregat
    t_ind_f <- t_ind_2[J(z2_i) , nomatch = 0L ]
    t_ind_f <- t_ind_f[! z1 %in% list_agregat[[i]] ]
    z2_f <- unique(t_ind_f$z2) # zones des z2 (carreaux) a la frontiere de l'agregat
    
    #On retrouve les observations concernees
    id_obs_risque_in <- t_ind_in_ag[z2 %in% z2_f]$id #identifiants des observations a risque a l'interieur de l'agregat
    nb_obs_in <- length(id_obs_risque_in) #nombre d'observations a la frontiere DANS l'agregat
    id_obs_risque_out <- t_ind_f$id #identifiants des observations a risque a l'exterieur de l'agregat
    nb_obs_out <- length(id_obs_risque_out) #nombre d'observations a la frontiere HORS agregat
    
    #On complete le data.frame
    ag_c <- paste0(list_agregat[[i]], collapse = "-") #liste des noms de communes de l'agregat concatenes
    df_in = df_out = data.frame()
    if(nb_obs_in < threshold)
      df_in <- data.frame(id_obs = id_obs_risque_in,
                          agregat_z1 = ag_c,
                          agregat_size = taille_ag[[i]],
                          nb_obs = nb_obs_in,
                          type_diff = "internal")
    
    if(nb_obs_out < threshold)
      df_out <- data.frame(id_obs = id_obs_risque_out,
                           agregat_z1 = ag_c,
                           agregat_size = taille_ag[[i]],
                           nb_obs = nb_obs_out,
                           type_diff = "external")
    
    df <- rbind(df_in, df_out)
    res <- rbind(res, df)
    
    if(verbose) pb$tick()
  }
  
  
  #On ne retient qu'une ligne par observation (suppression des doublons):
  res <- as.data.table(res) #on transforme en data.table
  #On retient par ordre de priorite
  # 1 - le nombre minimal d'observations du groupe
  res[, min_nb_obs := min(nb_obs), by=.( id_obs)]
  res <- res[nb_obs == min_nb_obs, ]
  res[, min_nb_obs := NULL]
  # 2 -l'agregat de taille minimale
  res[, min_taille_ag := min(agregat_size), by=.( id_obs)]
  res <- res[agregat_size == min_taille_ag, ]
  res[, min_taille_ag := NULL]
  # 3 - pour ce qui reste on choisit au hasard
  res <- res[ !duplicated(id_obs), ]
  
  #On retourne le resultat
  return(res)
}



#' Perform all the process to detect risky observations
#'
#' Allow from a table of observations for which there are
#' two different nomenclatures (z1 and z2) to determine the 
#' observations at risk when using the differentiation technique
#'
#' Risky observations because of differentiation are the ones for
#' which information can be deduced on agregates smaller than the 
#' confidentiality threshold. For example, considering the confidentiality
#' threshold is 10 and if by making the difference
#' between some categories of z1 and some categories of z2 one can 
#' deduce the value of a variable for 5 observations, then those 5 
#' observations are considered as "risky".
#'
#' @param t_ind The table of observations (data.frame or data.table). Each row
#' correspond to an observtion and for each observation we must know in which
#' category of the z1 nomenclature it belongs and in which category of the z2
#' nomenclature.
#' @param threshold Strictly positive integer indicating the confidentiality
#' threshold. Observations are considered at risk if one can deduce information
#' on a agregate of n observations where n < threshold.
#' @param max_agregate_size Integer indicating the maximal size of agregates
#' which are tested exhaustively. If that number is too large (greater than 30), the
#' computations may not end because of the combinations number that can become very large.
#' Also the RAM can be overloaded.
#' @param save_file Character indicating the suffix of the name of the saved results.
#' If is null, results are not writing on the hardware. The path root is taken from the
#' working directory (getwd()).
#' @param simplify Boolean. If TRUE then the graph simplification (merging + splitting)
#' occures. Otherwise the exhaustive search is directly applied on the original graph.
#' @param verbose Boolean. If TRUE (default), the different steps of the process are 
#' notified and progress bars provide an estimation of time left.
#'
#' @return As an output there is a data.table or data.frame with five columns : 
#' \enumerate{
#' \item $id_obs for the observation at risk
#' \item $agregat for the agregate of categories from z1 nomenclature on which the 
#' differentiation is performed
#' \item $agregat_size indicating the number of categories composing the agregate
#' \item $nb_obs the number of observations on which information is deduced when 
#' the differentiation is computed (nb_obs must be stricly inferior to $threshold)
#' \item $type_diff the type of differentiation between "internal" or "external".
#' }
#' 
#' @examples 
#' res_diff <- find_pbm_diff(t_ex,threshold = 5,max_agregate_size = 15)
#' 
#' @export
find_pbm_diff <- function(t_ind, threshold, max_agregate_size, save_file = NULL, simplify = TRUE, verbose = TRUE){
  
  #define variables to avoid NOTE of "no visible binding for global
  # variable" when running R CMD check
  z1 = z2 = NULL
  
  if(verbose) message("******** Start of the process *********")
  
  to_save <- !is.null(save_file) #si save_file est renseignee, cela permet de sauvegarder les matrices croisees apres agregations
  
  t_ind <- as.data.table(t_ind) #conversion to data.table format
  
  if(verbose) message("< --- Creation of the crossing matrix --- >")
  t_ind <- simplify_z2_rem(t_ind)
  t_crois <- tab_crois(t_ind)
  t_crois <- t_crois[z1 != "blanchi" & z2 != "blanchi"]
  t_crois <- simplify_z2_fus(t_crois)
  m_crois <- matrix_crois(t_crois)
  

  #if(to_save) saveRDS(m_crois, paste0("Resultats_diffman/m_crois_",save_file,".RDS"))
  
  if(simplify){ #one can choose to skip these steps of graph reduction if desired
    
    if(verbose) message("< --- Merging method 1 --- > ")
    m_crois <- agregate(m_crois, threshold, methode = "m1", verbose = verbose)
    
    #if(to_save) saveRDS(m_crois, paste0("Resultats_diffman/m_crois_ag_m1_",save_file,".RDS"))
    
    if(verbose) message("< --- Merging methods 1 and 2 --- >")
    m_crois <- agregate(m_crois, threshold, methode = "both", verbose = verbose)
    
    #if(to_save) saveRDS(m_crois, paste0("Resultats_diffman/m_crois_ag_m2_",save_file,".RDS"))
    
    if(sum(dim(m_crois)==0)>0) {
      message("No differentiation problems detected !")
      return(NULL)
    }
    
    if(verbose) message("< --- Splitting the graph --- >")
    l_decomp <- decompose_m_crois(m_crois, max_agregate_size)
  }
  else{
    l_decomp <- comp_connexe_list(m_crois)
  }
  
  if(verbose) message("< --- Exhaustive search of differentiation problems --- >")
  l_ag <- search_diff_agregate(l_decomp, threshold, max_agregate_size)
  l_ag <- desagregate_list(l_ag)
  
  if(verbose) message("< --- Identification of risky observations --- >")
  obs_risque <- find_id_obs_risque(l_ag, t_ind, threshold)
  if(to_save) {
    saveRDS(obs_risque, paste0(save_file,".RDS"))
  }
  
  if(verbose) message("******** End of the process :) *********")
  
  message("There are ", 
          nrow(obs_risque), 
          " risky observations because of differentiation ",
          "when confidentiality threshold is set to ", threshold)
  
  return(obs_risque)
}

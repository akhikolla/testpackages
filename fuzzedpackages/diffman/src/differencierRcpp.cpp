// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <algorithm>
#include <list>        // std::list
#include <iostream>    // std::iostream
#include <iterator>
using namespace Rcpp;
using namespace std;

/**
 * Fonction permettant de détecter une différenciation pour un agrégat donné
 * 
 *  @param iPopAgregat : entier contenant le nombre d'observations de l'agrégat courant
 *  @param vNbLiensAgregat : vecteur d'entiers contenant, pour chaque carreau, le nombre de territoires contigüs pour l'agrégat
 *  @param vNbLiensTparC : vecteur d'entiers contenant, pour chaque carreau, le nombre de territoires contigüs
 *  @param iSeuil : entier indiquant le seuil en dessous duquel il y a rupture du secret statistique
 *  @param vNbObsCarreaux : vecteur d'entiers contenant le nombre d'observations pour chaque carreau
 *  
 *  @return vecteur de deux entiers. 
 *            Le premier indique le nombre d'observations impliquées dans une différenciation intérieure. 
 *            Le deuxième indique le nombre d'observations impliquées dans une différenciation exterieure.
 *  
 *  @author Arlindo Dos Santos / PSAR Analyse urbaine
 *  @version 0.1 06/04/2018 
 */
IntegerVector verifRcpp(const int iPopAgregat
                          , const IntegerVector vNbLiensAgregat
                          , const IntegerVector vNbLiensTparC
                          , const int iSeuil
                          , const IntegerVector vNbObsCarreaux)
{
  int iNbCarreaux = vNbLiensAgregat.size();
  int iNbObsExterieur = 0;
  int iNbObsInterieur = 0;
  
  for(int iCarreau = 0; iCarreau < iNbCarreaux; ++iCarreau)
  {
    if((vNbLiensAgregat(iCarreau) > 0))
    {
      if(vNbLiensAgregat(iCarreau)  < vNbLiensTparC(iCarreau))
        iNbObsExterieur += vNbObsCarreaux(iCarreau);
      else if(vNbLiensAgregat(iCarreau) == vNbLiensTparC(iCarreau))
        iNbObsInterieur += vNbObsCarreaux(iCarreau);
    }
  }
  
  return(IntegerVector::create(
      (abs(iPopAgregat - iNbObsInterieur) < iSeuil) * abs(iPopAgregat - iNbObsInterieur),
                                            (abs(iNbObsExterieur + iNbObsInterieur - iPopAgregat) < iSeuil) * abs(iNbObsExterieur + iNbObsInterieur - iPopAgregat)
  ));
}

/**
 * Fonction récursive générant tous les agrégats possibles (sans répétition) et vérifiant s'il y a différenciation
 * 
 *  @param vAgregat : vecteur d'entiers (de taille fixe) contenant l'indice de l'agrégat courant. Sa taille courante est précisée dans iIndiceFinAgregat
 *  @param iIndiceFinAgregat : entier indiquant la position de fin du tableau vAgregat (nombre de territoires de l'agrégat courant)
 *  @param iTailleCible : entiers indiquant la taille des agrégats à tester
 *  @param vTerritoiresDejaExplores : vecteur d'entiers indiquant pour chaque territoire s'il a déjà été exploré
 *  @param iPopAgregat : entier contenant le nombre d'observations de l'agrégat courant
 *  @param vNbLiensAgregat : vecteur d'entiers contenant, pour chaque carreau, le nombre de territoires contigüs pour l'agrégat
 *  @param iSeuil : entier indiquant le seuil en dessous duquel il y a rupture du secret statistique
 *  @param vNbLiensTparC : vecteur d'entiers contenant, pour chaque carreau, le nombre de territoires contigüs
 *  @param vNbObsTerritoire : vecteur d'entiers contenant le nombre d'observations pour chaque territoire
 *  @param vNbObsCarreaux : vecteur d'entiers contenant le nombre d'observations pour chaque carreau
 *  @param mContiguiteT : matrice d'entiers indiquant la contiguité territoire/territoire - la diagonale est supposée à zéro
 *  @param mContiguiteTC : matrice d'entiers indiquant la contiguité territoire/carreau - les territoires sont en ligne et les carreaux en colonne
 *  @param lResultat :  Paramètre passé par référence afin de recevoir le résultat produit
 *                      Une liste où chaque élément est une différenciation. 
 *                      Chaque différenciation est un vecteur contenant les indices des territoires constituant l'agrégat 
 *                      puis le nombre d'observations impliquées dans une différenciation intérieure 
 *                      et enfin le nombre d'observations impliquées dans une différenciation extérieure
 *  
 *  @return void - le résultat est retourné par référence via le paramètre lResultat
 *  
 *  @author Arlindo Dos Santos / PSAR Analyse urbaine
 *  @version 0.2 27/07/2018 
 */
void explorerRcpp(IntegerVector vAgregat
                    , int iIndiceFinAgregat
                    , const int iTailleCible
                    , IntegerVector vTerritoiresDejaExplores
                    , const int iPopAgregat
                    , IntegerVector vNbLiensAgregat
                    , const int iSeuil
                    , const IntegerVector vNbLiensTparC
                    , const IntegerVector vNbObsTerritoire
                    , const IntegerVector vNbObsCarreaux
                    , const IntegerMatrix mContiguiteT
                    , const IntegerMatrix mContiguiteTC
                    , std::list<IntegerVector>& lResultat)
{
  if(iIndiceFinAgregat < iTailleCible)
  {
    // aggrégat à agrandir
    IntegerVector vIndicesVoisins;
    int iNbterritoires = vNbObsTerritoire.size();
    
    // recherche de tous les voisins possibles
    if(iIndiceFinAgregat == 0)
    {
      // aggrégat vide
      vIndicesVoisins = IntegerVector(iNbterritoires);
      for (int iIndice = 0; iIndice < iNbterritoires; ++iIndice)
        vIndicesVoisins(iIndice) = iIndice;
    }
    else
    {
      // aggrégat incomplet
      int iIndiceAgregat;
      int iIndiceCheminAgregat;
      for(int iTerritoire = 0; iTerritoire < iNbterritoires; ++iTerritoire)
      {
        for(iIndiceCheminAgregat = 0; iIndiceCheminAgregat < iIndiceFinAgregat; ++iIndiceCheminAgregat)
        {
          iIndiceAgregat = vAgregat(iIndiceCheminAgregat);
          
          // vérifier que les territoires sont contigus
          if(mContiguiteT(iTerritoire, iIndiceAgregat) == 1)
          {
            // ne pas ajouter un territoire déjà exploré (ou déjà dans l'agrégat)
            if(vTerritoiresDejaExplores[iTerritoire] == 0)
              vIndicesVoisins.push_back(iTerritoire);
            
            // passer au territoire suivant car le territoire courant est contigu à au moins 1 territoire de l'agrégat
            // pas la peine de vérifier qu'il l'est aussi pour un autre territoire de l'agrégat
            break;
          }
        }
      }
    }
    
    int iNbVoisins = vIndicesVoisins.size();
    if (iNbVoisins == 0)
    {
      // Rcout << "aggrégat ne pouvant pas être agrandi jusqu'à la taille requise" << std::endl;
      return;
    }
    else
    {
      int iNbCarreaux = vNbObsCarreaux.size();
      int k;
      for (int j = 0; j < iNbVoisins; ++j)
      {
        IntegerVector vNbLiensAgregatNew(iNbCarreaux);
        int iIndiceVoisinCourant = vIndicesVoisins(j);
        
        for(k = 0; k < iNbCarreaux; ++k)
          vNbLiensAgregatNew(k) = vNbLiensAgregat(k) + mContiguiteTC(iIndiceVoisinCourant, k);
        
        vAgregat(iIndiceFinAgregat) = iIndiceVoisinCourant;
        ++iIndiceFinAgregat;
        
        vTerritoiresDejaExplores(iIndiceVoisinCourant) = 1;
        
        explorerRcpp(vAgregat
                       , iIndiceFinAgregat
                       , iTailleCible
                       , clone(vTerritoiresDejaExplores)
                       , iPopAgregat + vNbObsTerritoire(iIndiceVoisinCourant)
                       , vNbLiensAgregatNew
                       , iSeuil
                       , vNbLiensTparC
                       , vNbObsTerritoire
                       , vNbObsCarreaux
                       , mContiguiteT
                       , mContiguiteTC
                       , lResultat);
        
        vTerritoiresDejaExplores[iIndiceVoisinCourant] = 1;
        --iIndiceFinAgregat;
      }
    }
  }
  else
  {
    // l'aggrégat est de la taille voulue
    IntegerVector verification_locale = verifRcpp(iPopAgregat, vNbLiensAgregat, vNbLiensTparC, iSeuil, vNbObsCarreaux);
    
    // Rcout << "\nagregat: " << vAgregat[0] + 1 << " " << vAgregat[1] + 1 << " " << vAgregat[2] + 1 << " " << vAgregat[3] + 1;
    
    if ((verification_locale(0) + verification_locale(1)) > 0)
    {
      // // une différenciation a été détectée
      // Rcout << std::endl << "agregat";
      // for(int iDebug = 0; iDebug < iIndiceFinAgregat; ++iDebug)
      //   Rcout << " " << vAgregat(iDebug) + 1;
      
      IntegerVector nouvelleDifferenciation = vAgregat + 1;
      nouvelleDifferenciation.push_back(verification_locale(0));
      nouvelleDifferenciation.push_back(verification_locale(1));
      lResultat.push_back(nouvelleDifferenciation);
    }
  }
  return;
}

//' Recherche exhaustive des problemes de differentiation
//' 
//' Fonction exposée permettant de détecter les différenciations territoire/carreau
//' Elle prépare les paramètres et appelle la fonction récursive explorerRcpp
//' 
//' Remarque : Il est inutile d'explorer des tailles d'agrégats supérieures à nbTerritoires/2
//' 
//' @param iTailleCible : entiers indiquant la taille des agrégats à tester
//' @param iSeuil : entier indiquant le seuil en dessous duquel il y a rupture du secret statistique
//' @param vNbObsTerritoire : vecteur d'entiers contenant le nombre d'observations pour chaque territoire
//' @param vNbObsCarreaux : vecteur d'entiers contenant le nombre d'observations pour chaque carreau
//' @param mContiguiteT : matrice d'entiers indiquant la contiguité territoire/territoire - la diagonale est supposée à zéro
//' @param mContiguiteTC : matrice d'entiers indiquant la contiguité territoire/carreau - les territoires sont en ligne et les carreaux en colonne
//'  
//' @return liste de vecteurs dont chaque élément est une différenciation. 
//'          Chaque élément est un vecteur constitué 
//'          - des indices des territoires de l'agrégat 
//'          - du nombre d'observations différenciées à l'intérieur 
//'          - et enfin du nombre d'observations différenciées à l'extérieur
//'  
//' @author Arlindo Dos Santos / PSAR Analyse urbaine
// [[Rcpp::export]]
std::list<IntegerVector> differencierRcpp(const int iTailleCible
                                            , const int iSeuil
                                            , const IntegerVector vNbObsTerritoire
                                            , const IntegerVector vNbObsCarreaux
                                            , const IntegerMatrix mContiguiteT
                                            , const IntegerMatrix mContiguiteTC)
{
  std::list<IntegerVector> lResultat;
  
  const int iNbCarreaux = vNbObsCarreaux.size();
  const int iNbTerritoires = vNbObsTerritoire.size();
  
  int iLigne;
  IntegerVector vNbLiensTparC(iNbCarreaux);
  for(int iColonne = 0; iColonne < iNbCarreaux; ++iColonne)
    for(iLigne = 0; iLigne < iNbTerritoires; ++iLigne)
      vNbLiensTparC(iColonne) += mContiguiteTC(iLigne, iColonne);
  
  explorerRcpp(
    IntegerVector(iTailleCible)
    , 0L
    , iTailleCible
    , IntegerVector(iNbTerritoires)
    , 0L
    , IntegerVector(vNbObsCarreaux.size())
    , iSeuil
    , vNbLiensTparC
    , vNbObsTerritoire
    , vNbObsCarreaux
    , mContiguiteT
    , mContiguiteTC
    , lResultat
  );
  return(lResultat);
}

/*** R
# cat("\f")
# debut <- proc.time()
# 
# iTailleAgregatCible <- 1L
# 
# lResultat = differencierRcpp(iTailleAgregatCible
#                        , 11L
#                       , commune0$men
#                       , carreau0$men
#                       , mat_conti_commune
#                       , mat_conti_commune_carreau
#                       )
# 
# dfResultat <- t(as.data.frame(lResultat, fix.empty.names = FALSE))
# 
# if(nrow(dfResultat) > 0)
# {
#   nomColonnes <- c()
#   for(iCompteur in 1:(iTailleAgregatCible))
#   {
#     dfResultat[, iCompteur] <- commune0[dfResultat[, iCompteur], ]$depcom
#     nomColonnes <- c(nomColonnes, paste0("depcom", iCompteur))
#   }
#   nomColonnes <- c(nomColonnes, "interieur", "exterieur")
#   colnames(dfResultat) <- nomColonnes
# }
# 
# dfResultat
# proc.time() - debut

*/

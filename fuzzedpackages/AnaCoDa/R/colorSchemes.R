# colors for mixture elements, make as long as possible
.mixtureColors <- c("black", "red", "blue", "green4", "magenta", "darkorange", "deeppink3", "gray", "cyan", "brown", "yellow")
# for up to six codons
.codonColors <- list(GCA="blue", GCC="darkorange", GCG="purple", GCT="green4", #Ala
                     TGC="darkorange", TGT="green4", #Cys
                     GAC="darkorange", GAT="green4", #Asp
                     GAA="blue", GAG="purple", #Glu
                     TTC="darkorange", TTT="green4", #Phe
                     GGA="blue", GGC="darkorange", GGG="purple", GGT="green4", #Gly
                     CAC="darkorange", CAT="green4", #His
                     ATA="blue", ATC="darkorange", ATT="green4", #Ile
                     AAA="blue", AAG="purple", #Lys
                     CTA="blue", CTC="darkorange", CTG="purple", CTT="green4", TTA="darkturquoise", TTG="deeppink3", #Leu
                     ATG="purple",
                     AAC="darkorange", AAT="green4", #Asn
                     CCA="blue", CCC="darkorange", CCG="purple", CCT="green4", #Pro
                     CAA="blue", CAG="purple", #Gln
                     CGA="blue", CGC="darkorange", CGG="purple", CGT="green4", AGA="darkturquoise", AGG="deeppink3", #Arg
                     TCA="blue", TCC="darkorange", TCG="purple", TCT="green4", #Ser4
                     ACA="blue", ACC="darkorange", ACG="purple", ACT="green4", #Thr
                     GTA="blue", GTC="darkorange", GTG="purple", GTT="green4", #Val
                     TAC="darkorange", TAT="green4", #Tyr
                     AGC="darkorange", AGT="green4", #Ser2
                     TGG="blue",
                     TAA="blue", TAG="purple", TGA="darkturquoise") #Stop

.ribModelConstants <- list(
                    #FONSE/ROC "constants"
                    deltaM = "Mutation", deltaEta = "Selection", deltaOmega = "Selection",
                    #PA "constants
                    alpha = "Alpha",lambdaPrime = "LambdaPrime")


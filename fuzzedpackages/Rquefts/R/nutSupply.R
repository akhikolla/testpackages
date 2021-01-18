
nutSupply1 <- function(pH, SOC, Kex, Polsen, Ptotal=NA) {
	# SOC in g/kg, Kex in mmol/kg, P in mg/kg
	# Janssen et al., 1990. Geoderma 46: 299-318, Table 2

	# for recycling
	d <- cbind(as.vector(pH), as.vector(SOC), as.vector(Kex), as.vector(Polsen), as.vector(Ptotal))
	colnames(d) <- c("pH", "SOC", "Kex", "Polsen", "Ptotal")
	
	fN <- 0.25 * (d[,"pH"] - 3)
	N_supply <- fN * 6.8 * d[,"SOC"]
	
	fP <- 1 - 0.5 * (d[,"pH"] - 6)^2

	P_supply <- fP * 0.35 * d[,"SOC"] + 0.5 * d[,"Polsen"]

	i <- which(!is.na(d[,"Ptotal"]))
	if (length(i) > 0) {
		P_supply[i] <- fP[i] * 0.014 * d[,"Ptotal"][i] + 0.5 * d[,"Polsen"][i]
	}
	
	fK <- 0.625 * (3.4 - 0.4 * d[,"pH"])
	#Kex <- d[,3]
	K_supply <- (fK * 400 * d[,"Kex"]) / (2 + 0.9 * d[,"SOC"])
	
	x <- cbind(N_base_supply=N_supply, P_base_supply=P_supply, K_base_supply=K_supply)

	x[x<0] <- 0
	x
}



..getfP <- function(pH) {
# this version is not used!
# Sattari et al
	k <- is.na(pH)
	pH <- pH[!k]
	fP <- rep(NA, length(pH))

	i <- j <- pH < 4.7
	fP[i] <- 0.02 

	i <- (!j) & pH < 6
	fP[i] <- 1 - 0.5 * (pH[i] - 6)^2
	j <- i | j

	i <- (!j) & pH < 6.7
	fP[i] <- 1
	j <- i | j

	i <- (!j) & pH < 8
	fP[i] <- 1 - 0.25 * (pH[i] - 6.7)^2
	j <- i | j

	fP[!j] <- 0.02
	
	fPna <- rep(NA, length(k))
	fPna[!k] <- fP
	fPna
}

.getfN <- function(pH) {
# this version is not used!
# Sattari et al
	k <- is.na(pH)
	pH <- pH[!k]

	fN <- rep(NA, length(pH))

	i <- j <- pH < 4.7
	fN[i] <- 0.4 

	i <- (!j) & pH < 7
	fN[i] <- 0.25 * (pH[i] - 3)
	j <- i | j

	fN[!j] <- 1

	fNna <- rep(NA, length(k))
	fNna[!k] <- fN
	fNna
}

.getfK <- function(pH) {
# this version is not used!
# Sattari et al
	k <- is.na(pH)
	pH <- pH[!k]
	i <- as.integer(cut(pH, c(-1, 4.5,  6.8, 15)))
	
	fK <- rep(1, length(pH))
	fK[i==2] <- 6.1 * pH[i==2]^(-1.2)
	fK[i==3] <- 0.6

	fKna <- rep(NA, length(k))
	fKna[!k] <- fK
	fKna
}


# these are used
# identical for N, slight approx for P and K
.get_fn <- approxfun(c(4.7,7), c(0.4, 1), rule=2)

.get_fp <- approxfun(c(4.7, 5, 5.5, 6, 6.7, 7, 7.5, 8), c(0.02, .5, 0.875, 1, 1, 0.9775, 0.84, 0.02), rule=2)

.get_fk <- approxfun(c(4.5, 5, 6, 6.8), c(1, 0.8842312, 0.71, 0.6), rule=2)
#.old_fk <- function(ph) 0.625 * (3.4 - 0.4 * ph)


nutSupply2 <- function(temp, pH, SOC, Kex, Polsen, Ptotal) {
# SOC in g/kg, Kex in mmol/kg, P in mg/kg
# NPK supply after Sattari et al., 2017
	
	d <- cbind(pH=as.vector(pH), SOC=as.vector(SOC), Kex=as.vector(Kex), Polsen=as.vector(Polsen), Ptotal=as.vector(Ptotal), temp=as.vector(temp))

# temperature adjustment
# TeffN replaces fixed value 6.8 
	TeffN = 2 * 2 ^((d[,"temp"]-9)/9)
	fN <- .get_fn(d[,"pH"])
	N_supply <- fN * TeffN * d[,"SOC"]

	fP <- .get_fp(d[, "pH"])
	P_supply <- fP * 0.014 * d[,"Ptotal"] + 0.5 * d[,"Polsen"]
	
	fK <- .get_fk(d[,"pH"])
	K_supply <- (fK * 400 * d[,"Kex"]) / (2 + 0.9 * d[,"SOC"])
	
	cbind(N_base_supply=N_supply, P_base_supply=P_supply, K_base_supply=K_supply)
}


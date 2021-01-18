GetNoise <- function(data, noisePercent = "median") {
    if (is.null(noisePercent)) noisePercent = "median"
    if (noisePercent == "median") {
        sds <- matrixStats::colSds(data) # Vectorize the calulation of column mean
        noise <- median(sds)
    } else {
        sds <- matrixStats::colSds(data)
        sds <- sort(sds)
        ind <- round(length(sds) * noisePercent/100)
        noise <- sds[ind]
    }
    noise
}

AddNoisePerturb <- function(data, noise) {
    rowNum <- nrow(data)
    colNum <- ncol(data)
    dat = rnorm(rowNum * colNum, mean = 0, sd = noise)
    epsilon<- matrixStats::allocMatrix(rowNum,colNum,value = 0.0)
    epsilon <- dat
    

 
    list(
        data = data + epsilon,
        ConnectivityMatrixHandler = function(connectivityMatrix, ...) {
            connectivityMatrix
        }
    )
}

SubSampling <- function(data, percent = 80) {
    N = nrow(data)
    randCount = ceiling(N * percent/100)
    randOrder = sample(N, N)
    
    randMatrix = data[randOrder[1:randCount], ]
    list(
        data = randMatrix,
        ConnectivityMatrixHandler = function(connectivityMatrix, ...) {
            S = matrix(0, N, N)
            S[1:randCount, 1:randCount] <- connectivityMatrix
            for(i in 1:N){
                S[N,N] = 1
            }
            S[randOrder, randOrder]
        }
    )
}
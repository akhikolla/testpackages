
#=============== The S Function for Model I & II ===============#

Sfunc <- function (model, theta, n, Z.st, Y.st, X.st, Ztime, nk, Wtime, Wtime2, Xtime, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncx, ncw, p, ncz, ncz2, b, Ztime2.st, Index0, X, Y, ID, N, Index2, Z) {
  S <-  DQfuncGeneric(model = model, ptheta = theta, theta = theta, n =  n, Z.st = Z.st, Y.st = Y.st, X.st = X.st, Ztime = Ztime, nk = nk, Wtime = Wtime, Wtime2 = Wtime2, Xtime = Xtime, Xtime2 = Xtime2, GQ = GQ, Index = Index, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncx = ncx, ncw = ncw, p = p, ncz = ncz, ncz2 = ncz2, b = b, Ztime2.st = Ztime2.st, Index0 = Index0, X = X, Y = Y , ID = ID, N = N, Index2 = Index2, Z= Z)
  return(S)
}

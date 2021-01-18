
#=============== The S Function for Model I & II of Multiplicative Joint Model ===============#

SfuncMult <- function (model, theta, B.st, n, Y.st, b, Btime, Btime2, Index, Index0, Ztime, Ztime2, nknot, nk, Index1, rho, d, wGQ, ncz, ncb, N, Y, B, ID, Index2) {
  S <-  DQfuncMultGeneric(ptheta = theta, theta = theta, B.st = B.st, n = n, Y.st = Y.st, b = b, model = model, Btime = Btime, Btime2 = Btime2, Index = Index, Index0 = Index0, Ztime = Ztime, Ztime2 = Ztime2, nknot = nknot, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncz = ncz, ncb = ncb, N = N, Y = Y, B = B, ID = ID, Index2 = Index2)
  return(S)
}

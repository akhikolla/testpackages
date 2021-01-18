iopt <- function() {
  df <- data.frame(key=rep(NA, 23), value=rep("", 23), stringsAsFactors=FALSE)
  df[1, ] <- c(1, "NAD_83(2011/CORS96/2007) (North American plate fixed)")
  df[2, ] <- c(2, "NAD_83(PA11/PACP00) (Pacific plate fixed)")
  df[3, ] <- c(3, "NAD_83(MA11/MARP00) (Mariana plate fixed)")
  df[5, ] <- c(5, "WGS_84(transit) (NAD_83(2011) used)")
  df[6, ] <- c(6, "WGS_84(G730) (ITRF91 used)")
  df[7, ] <- c(7, "WGS_84(G873) (ITRF94 used)")
  df[8, ] <- c(8, "WGS_84(G1150) (ITRF2000 used)")
  df[9, ] <- c(9, "WGS_84(G1674) (ITRF2008 used)")
  df[10, ] <- c(10, "WGS_84(G1762) (IGb08 used)")
  df[11, ] <- c(11, "SIO/MIT_92 (ITRF91 used)")
  df[12, ] <- c(12, "ITRF88")
  df[13, ] <- c(13, "ITRF89")
  df[14, ] <- c(14, "ITRF90 or (PNEOS90/NEOS90)")
  df[15, ] <- c(15, "ITRF91")
  df[16, ] <- c(16, "ITRF92")
  df[17, ] <- c(17, "ITRF93")
  df[18, ] <- c(18, "ITRF94")
  df[19, ] <- c(19, "ITRF96")
  df[20, ] <- c(20, "ITRF97 or IGS97")
  df[21, ] <- c(21, "ITRF2000 or IGS00/IGb00")
  df[22, ] <- c(22, "ITRF2005 or IGS05")
  df[23, ] <- c(23, "ITRF2008 or IGS08/IGb08")
  df <- df[-c(4),]
  df
}



## Example pedigree with missing (unknown) birth year for some individuals
ped0 <- data.frame(     id=c( 1, 2, 3,  4, 5, 6, 7,  8, 9, 10, 11, 12, 13, 14),
                        fid=c( 0, 0, 0,  1, 1, 1, 3,  3, 3,  5,  4,  0,  0, 12),
                        mid=c( 0, 0, 0,  2, 0, 2, 2,  2, 5,  0,  0,  0,  0, 13),
                        birth_dt=c(NA, 0, 1, NA, 3, 3, 3, 3, 4, 4, 5, NA, 6, 6) + 2000)

## First run - using information from children
ped1 <- pedFixBirthYear(x=ped0, interval=1)

## Second run - using information from parents
ped2 <- pedFixBirthYear(x=ped1, interval=1, down=TRUE)

## Third run - using information from children, but with no success
ped3 <- pedFixBirthYear(x=ped2, interval=1)

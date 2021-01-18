## Example pedigree
ped <- data.frame(      id=1:10,
                       fid=c(0, 0, 0, 1, 1, 1, 3, 3, 3, 5),
                       mid=c(0, 0, 0, 2, 0, 2, 2, 2, 5, 0),
                  birth_dt=c(0, 0, 1, 2, 3, 3, 3, 4, 4, 5) + 2000)

## Set base population as those individuals that were born after year 2002
pedSetBase(x=ped, keep=ped$birth_dt > 2002, unknown=0)


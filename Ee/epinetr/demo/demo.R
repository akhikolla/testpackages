# Construct a population
pop <- Population(popSize = 200, map = map100snp, QTL = 20,
                  alleleFrequencies = runif(100),
                  broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40)

# Attach additive effects to all SNPs from a normal distribution
pop <- addEffects(pop, distrib = rnorm)

# Attach a random pairwise epistatic network to QTL
pop <- attachEpiNet(pop)

# Plot the epistatic network
plot(getEpiNet(pop))

# Attach a random epistatic network to QTLs with 3-way and 4-way
# interactions
pop <- attachEpiNet(pop, k=3:4)

# Plot the epistatic network
plot(getEpiNet(pop))

# Attach a scale-free pairwise epistatic network to QTLs
pop <- attachEpiNet(pop, scaleFree=TRUE)

# Plot the epistatic network
plot(getEpiNet(pop))

# Attach a scale-free pairwise epistatic network to QTLs with minimum
# degree of 2
pop <- attachEpiNet(pop, scaleFree=TRUE, m=2)

# Plot the epistatic network
plot(getEpiNet(pop))

# Attach a scale-free epistatic network to QTLs with 3-way and 4-way
# interactions
pop <- attachEpiNet(pop, scaleFree=TRUE, k=3:4)

# Plot the epistatic network
plot(getEpiNet(pop))

# Attach a scale-free epistatic network to QTLs with 3-way and 4-way
# interactions and a minimum degree of 2
pop <- attachEpiNet(pop, scaleFree=TRUE, k=3:4, m=2)

# Plot the epistatic network
plot(getEpiNet(pop))

# Run the simulation
pop2 <- runSim(pop, generations=150)

# Plot the run
plot(pop2)

# Retrieve pedigree from simulation
ped <- getPedigree(pop2)

# Run the pedigree dropper
pop3 <- runSim(pop, pedigree=ped)

# Plot the pedigree dropper run
plot(pop3)

#!/home2/bt/Rscript
#
# Run an imputation permutation test for a particular phenotype in an R/qtl
# cross file (assumed to be named input.csv).  Results are written to an output
# file via dput() so that multiple runs may be combined later.  The run count
# should be a serial integer identifying this run.
#

library(qtl)

args <- commandArgs(TRUE)
if (length(args) != 3)
    stop("usage: run_perm_test PHENO_COLUMN RUN_COUNT  NUM_PERMS")

pheno.col <- as.integer(args[1])
n.draws <- 64
run.count <- as.integer(args[2])
n.perm <- as.integer(args[3])


cross <- read.cross("csv", file="input.csv", genotypes=c("A", "B"))

# Specify that the cross data was generated from selfed RILs
class(cross)[1] <- "riself"

# Some markers are at the same coordinates, so apply a jittermap
cross <- jittermap(cross)

# Calculate the probabilities of the true underlying genotypes given the
# observed multipoint marker data, with allowance for genotyping errors
cross <- calc.genoprob(cross, step=1, map.function="kosambi")

# Create a grid showing the recombination fractions for all pairs of markers,
# or of the LOD scores between pairs of markers
cross <- est.rf(cross)

cross.hk <- sim.geno(cross, step=1, map.function="kosambi", error.prob=0.001,
                      n.draws=n.draws)

results <- scantwo(cross.hk, pheno.col=pheno.col, method="hk", n.perm=n.perm,
                   verbose=FALSE)

output.filename <- paste("perm_part_", pheno.col, "_", run.count, ".Rdat", sep="")
dput(results, output.filename)

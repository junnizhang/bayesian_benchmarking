
kBaseDir <- "/Users/johnbryant/Documents/bm/smoking"
kCalcDir <- file.path(kBaseDir, "Paper - calculations")
kOutDir <- kCalcDir
kPaperDir <- file.path(kBaseDir, "Paper")

library(latticeExtra)

setwd(kBaseDir)
load("census.RData")  ## loads data.frame 'census'


###############################################################################
## Set values that remain constant across simulations #########################
###############################################################################

## True finite-population quantities

smokes.all <- xtabs(~ age + sex + income,
                    data = census,
                    subset = smokes == 1L)
smokes.all <- Counts(smokes.all)
popn.all <- xtabs(~ age + sex + income,
                  data = census)
popn.all <- Counts(popn.all)
prob.true <- smokes.all / popn.all


## Finite population with region dimension

popn.all.reg <- xtabs(~ age + sex + income + region,
                      data = census)
popn.all.reg <- Counts(popn.all.reg)


## Size of sample

n.sample <- 60000


## Strata counts

region.popn <- table(census$region)
region.weight <- sqrt(region.popn)
strata.counts <- round(n.sample * prop.table(region.weight))


## Finite population correction

strata.popn <- table(census$region)
census$fpc <- strata.popn[match(census$region, names(strata.popn))]


## Multiplier for 'medium' benchmark variant

kMedium <- 0.5


## Headers

n.cell <- length(popn.all)
n.bench <- length(collapseDimension(popn.all, margin = "income"))

colnames.M <- paste("M", seq_len(n.cell), sep = ".")
colnames.W <- paste("W", seq_len(n.cell), sep = ".")
colnames.C <- paste("C", seq_len(n.cell), sep = ".")
colnames.D <- paste("D", seq_len(n.bench), sep = ".")


makeHeader <- function(cols) {
    id.vars <- c("iteration", "model.variant")
    ans <- c(id.vars, cols)
    ans <- paste(ans, collapse = ",")
    ans <- paste(ans, "\n", sep = "")
    ans
}

header.diagnostics <- makeHeader(c("accept.prob", "accept.bench",
                                   "autocorr", "psrf", "sigma",
                                   "n.young.rich", "n.young.rich.smokes"))
header.perform.cell <- makeHeader(c(colnames.M, colnames.W, colnames.C))
header.perform.bench <- makeHeader(colnames.D)


## Filenames

results.file.diagnostics <- "smoking_results_diagnostics.csv"
results.file.perform.cell <- "smoking_results_perform_cell.csv"
results.file.perform.bench <- "smoking_results_perform_bench.csv"



#################################################################################
## Functions for extracting results #############################################
#################################################################################


## These functions are closures, in that they include values fixed
## in simulation in their definitions.  It is essential to
## run the code in the "Set values that remain constant in simulations"
## section before using these functions.

makePSRF <- function(res) {
    mcmc <- MCMC(res,
                 where = c("model", "likelihood", "prob"),
                 sample = seq_len(n.cell))
    ans <- gelman.diag(mcmc, autoburnin = FALSE, multivariate = FALSE)
    max(ans$psrf[ , "Point est."])
}

makeSigma <- function(res) {
    sigma <- fetch(res, where = c("model", "prior", "sd"))
    mean(sigma)
}

diagnostics <- function(res) {
    metropolis <- metropolis(res)
    accept.prob <- as.numeric(metropolis["model.likelihood.prob", "acceptance"])
    uses.bench <- "model.benchmarks.values" %in% rownames(metropolis)
    if (uses.bench)
        accept.bench <- as.numeric(metropolis["model.benchmarks.values", "acceptance"])
    else
        accept.bench <- NA
    autocorr <- metropolis["model.likelihood.prob", "autocorr"]
    psrf <- makePSRF(res)
    sigma <- makeSigma(res)
    pop <- fetch(res, where = "exposure")
    n.young.rich <- sum(subarray(pop, age == "15-24" & income == "100001+"))
    y <- fetch(res, where = "y")
    n.young.rich.smokes <- sum(subarray(y, age == "15-24" & income == "100001+"))
    c(accept.prob,
      accept.bench,
      autocorr,
      psrf,
      sigma,
      n.young.rich,
      n.young.rich.smokes)
}

## Results based on cells obtained after aggregating away region,
## since region is not of substantive interest, but included only
## to account for sample design.
performanceCells <- function(res) {
    kProb = c(0.025, 0.975)
    finite.y <- finiteY(res, population = popn.all.reg)
    finite.y <- collapseDimension(finite.y, dimension = "region")
    finite.prob <- finite.y / popn.all
    mean.finite.prob <- collapseIterations(finite.prob, FUN = mean)
    quantiles.finite.prob <- collapseIterations(finite.prob, prob = kProb)
    lower <- slice(quantiles.finite.prob, dimension = "quantile", elements = 1L)
    upper <- slice(quantiles.finite.prob, dimension = "quantile", elements = 2L)
    error.sq <- (mean.finite.prob - prob.true)^2
    widths <- upper - lower
    in.credible.interval <- (lower <= prob.true) & (prob.true <= upper)
    c(as.numeric(error.sq),
      as.numeric(widths),
      100 * as.numeric(in.credible.interval))
}

performanceBenchmarks <- function(res, meanDirect, sdDirect) {
    finite.y <- finiteY(res, population = popn.all.reg)
    finite.y <- collapseDimension(finite.y, 
                                  margin = c(names(meanDirect), "iteration"))
    popn.all <- collapseDimension(popn.all, margin = names(meanDirect))
    finite.prob <- finite.y / popn.all
    mean.finite.prob <- collapseIterations(finite.prob, FUN = mean)
    as.numeric(abs(mean.finite.prob - meanDirect) / sdDirect)
}

writeToFile <- function(variant, output, filename) {
  line <- c(variant, output)
  line <- paste(line, collapse = ",")
  line <- paste(line, "\n", sep = "")
  cat(x = line, file = filename, append = TRUE)
}



#############################################################################
## Control values ###########################################################
#############################################################################

nBurnin <- 50000
nSim <- 50000
nChain <- 6
nThin <- 250

n.iteration <- 100

jump.theta.nonbench <- 0.4
jump.theta.bench <- 0.25
jump.bench <- 0.015



##############################################################################
## Simulation ################################################################
##############################################################################

set.seed(10)

setwd(kOutDir)
cat(header.diagnostics, file = results.file.diagnostics)
cat(header.perform.cell, file = results.file.perform.cell)
cat(header.perform.bench, file = results.file.perform.bench)

cat("Starting simulation at", format(Sys.time()), "\n")

for (iteration in seq_len(n.iteration)) {
    
    cat("Starting iteration", iteration, "\n")
    
    ## Draw sample ##########################################################
    
    s <- stratsample(strata = census$region, counts = strata.counts)
    sample <- census[s, ]
    smokes.sample <- xtabs(~ age + sex + region + income, 
                           data = sample, 
                           subset = smokes == 1L)
    smokes.sample <- Counts(smokes.sample)
    popn.sample <- xtabs(~ age + sex + region + income, 
                         data = sample)
    popn.sample <- Counts(popn.sample)
        
    
    ## Make benchmarks ######################################################
    
    design <- svydesign(id = ~1,
                        strata = ~region,
                        fpc = ~fpc,
                        data = sample)
    direct <- svyby(formula = ~smokes,
                    by = ~income, 
                    design = design,
                    FUN = svymean)
    mean.direct <- xtabs(smokes ~ income, data = direct)
    sd.direct <- xtabs(se ~ income, data = direct)
    mean.direct <- Values(mean.direct)
    sd.direct <- Values(sd.direct)
    bench.hard <- Benchmarks(mean = mean.direct, 
                             weights = popn.all.reg)
    bench.soft <- Benchmarks(mean = mean.direct, 
                             sd = sd.direct,
                             jump = jump.bench,
                             weights = popn.all.reg)
    bench.medium <- Benchmarks(mean = mean.direct, 
                               sd = kMedium * sd.direct,
                               jump = jump.bench,
                               weights = popn.all.reg)
    
    
    ## Closure for collecting results

    writeResults <- function(res, variant) {
        diag <- diagnostics(res)
        perform.cell <- performanceCells(res)
        perform.bench <- performanceBenchmarks(res = res,
                                               meanDirect = mean.direct,
                                               sdDirect = sd.direct)
        writeToFile(variant = variant,
                    output = diag,
                    filename = results.file.diagnostics)
        writeToFile(variant = variant,
                    output = perform.cell,
                    filename = results.file.perform.cell)
        writeToFile(variant = variant,
                    output = perform.bench,
                    filename = results.file.perform.bench)
    }
    
    
    ## Estimate models ###############################################################
    
    ## Not benchmarked

    variant <- c(iteration, "None")
    setwd(kOutDir)
    res <- estimateModel(Model(y ~ Binomial(prob ~ age * income + sex + region,
                                            jump = jump.theta.nonbench),
                               age ~ Exch(),
                               age:income ~ Exch()),
                         y = smokes.sample,
                         exposure = popn.sample,
                         filename = tempfile(),
                         nBurnin = nBurnin,
                         nSim = nSim,
                         nChain = nChain,
                         nThin = nThin)
    psrf <- makePSRF(res)
    cat(paste(variant, collapse = " "), "psrf = ", psrf, "\n")
    writeResults(res = res, variant = variant)


    ## Hard benchmarks
    
    variant <- c(iteration, "Hard")
    setwd(kOutDir)
    res <- estimateModel(Model(y ~ Binomial(prob ~ age * income + sex + region,
                                            benchmarks = bench.hard,
                                            jump = jump.theta.bench),
                               age ~ Exch(),
                               age:income ~ Exch()),
                         y = smokes.sample,
                         exposure = popn.sample,
                         filename = tempfile(),
                         nBurnin = nBurnin,
                         nSim = nSim,
                         nChain = nChain,
                         nThin = nThin)
    psrf <- makePSRF(res)
    cat(paste(variant, collapse = " "), "psrf = ", psrf, "\n")
    writeResults(res = res, variant = variant)


    ## Soft benchmarks
    
    variant <- c(iteration, "Soft")
    setwd(kOutDir)
    res <- estimateModel(Model(y ~ Binomial(prob ~ age * income + sex + region,
                                            benchmarks = bench.soft,
                                            jump = jump.theta.bench),
                               age ~ Exch(),
                               age:income ~ Exch()),
                         y = smokes.sample,
                         exposure = popn.sample,
                         filename = tempfile(),
                         nBurnin = nBurnin,
                         nSim = nSim,
                         nChain = nChain,
                         nThin = nThin)
    psrf <- makePSRF(res)
    cat(paste(variant, collapse = " "), "psrf = ", psrf, "\n")
    writeResults(res = res, variant = variant)


    ## Medium benchmarks
    
    variant <- c(iteration, "Medium")
    setwd(kOutDir)
    res <- estimateModel(Model(y ~ Binomial(prob ~ age * income + sex + region,
                                            benchmarks = bench.medium,
                                            jump = jump.theta.bench),
                               age ~ Exch(),
                               age:income ~ Exch()),
                         y = smokes.sample,
                         exposure = popn.sample,
                         filename = tempfile(),
                         nBurnin = nBurnin,
                         nSim = nSim,
                         nChain = nChain,
                         nThin = nThin)
    psrf <- makePSRF(res)
    cat(paste(variant, collapse = " "), "psrf = ", psrf, "\n")
    writeResults(res = res, variant = variant)

}

cat("Finishing simulation at", format(Sys.time()), "\n")



###################################################################################
### Examine the results ###########################################################
###################################################################################


## Performance measures ############################################################

colnames.M <- paste("M", 1:84, sep = ".")
colnames.W <- paste("W", 1:84, sep = ".")
colnames.C <- paste("C", 1:84, sep = ".")
colnames.D <- paste("D", 1:7, sep = ".")
n.cell <- 84
n.bench <- 7


setwd(kBaseDir)

pc <- read.csv("smoking_results_perform_cell.csv")
pb <- read.csv("smoking_results_perform_bench.csv")
pc <- aggregate(pc[c(colnames.M, colnames.W, colnames.C)],
                pc["model.variant"],
                mean)
pb <- aggregate(pb[colnames.D],
                pb["model.variant"],
                mean)
pc <- reshape(pc,
              varying = list(colnames.M, colnames.W, colnames.C),
              v.names = c("M", "W", "C"),
              idvar = "model.variant",
              timevar = "id",
              times = seq_len(n.cell),
              direction = "long")
pb <- reshape(pb,
              varying = list(colnames.D),
              v.names = "D",
              idvar = "model.variant",
              timevar = "id",
              times = seq_len(n.bench),
              direction = "long")
pc$M <- 1000 * sqrt(pc$M)
pc$W <- pc$W * 1000
pc <- reshape(pc,
              varying = list(c("M", "W", "C")),
              v.names = "value",
              idvar = c("model.variant", "id"),
              timevar = "measure",
              times = c("M", "W", "C"),
              direction = "long")
pb$measure <- "D"
names(pb)[match("D", names(pb))] <- "value"
pcb <- rbind(pc, pb)
rownames(pcb) <- NULL
pcb$measure <- factor(pcb$measure, levels = c("D", "M", "W", "C"))
pcb$model.variant <- factor(pcb$model.variant,
                            levels = c("Hard", "Medium", "Soft", "None"))
panel.special <- function(x, y, ...) {
    panel.bwplot(x, y, ...)
    medians <- tapply(x, y, median)
    format <- if (all(medians < 2)) "%3.2f" else "%3.1f"
    labels <- sprintf(format, medians)
    x <- as.numeric(medians)
    y <- match(names(medians), levels(y))
    panel.text(labels = labels, x = x, y = y + 0.4, cex = 0.8)
}
strip.special <- strip.custom(factor.levels =
                                  c(expression(Discrepancy~~italic(D[l])),
                                    expression(sqrt(MSE[italic(asl)])%*%1000),
                                    expression(Width~~italic(W[asl]^0.95)%*%1000),
                                    expression(Coverage~~italic(C[asl]^0.95))))
p <- bwplot(model.variant ~ value | measure,
            data = pcb,
            box.ratio = 0.9,
            par.settings = list(box.rectangle = list(col = "black"),
                box.umbrella = list(col = "black"),
                plot.symbol = list(col = "black", pch = "|"),
                fontsize = list(text = 8, points = 6),
                strip.background = list(col = "grey90"),
                layout.heights=list(strip=1.3)),
            scales = list(tck = 0.4,
                          x = list(relation = "free"),
                          y = list(labels = c(expression(Exact),
                                              expression(Inexact: lambda==1),
                                              expression(Inexact: lambda==0.5),
                                              expression(None)))),
            strip = strip.special,
            pch = "|",
            xlab = "",
            layout = c(4, 1),
            panel = panel.special)
graphics.off()
setwd(kBaseDir)
pdf(file = "smoking_performance.pdf", width = 6, height = 2.1)
plot(p)
dev.off()



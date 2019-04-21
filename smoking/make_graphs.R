
n.cell <- 84
n.bench <- 7
colnames.M <- paste("M", seq_len(n.cell), sep = ".")
colnames.W <- paste("W", seq_len(n.cell), sep = ".")
colnames.C <- paste("C", seq_len(n.cell), sep = ".")
colnames.D <- paste("D", seq_len(n.bench), sep = ".")

pc <- read.csv("smoking/smoking_results_perform_cell.csv")
pb <- read.csv("smoking/smoking_results_perform_bench.csv")
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
                            levels = c("Soft", "Medium", "Hard", "None"))
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
                                  c(expression(Consistency~~italic(D[j])),
                                    expression(sqrt(MSE[italic(i)])%*%1000),
                                    expression(Width~~italic(W[i])%*%1000),
                                    expression(Coverage~~italic(C[i]))))
p <- bwplot(model.variant ~ value | measure,
            data = pcb,
            box.ratio = 0.9,
            par.settings = list(box.rectangle = list(col = "black"),
                box.umbrella = list(col = "black"),
                plot.symbol = list(col = "black", pch = "|"),
                fontsize = list(text = 8, points = 6),
                strip.background = list(col = "grey90")),
            scales = list(tck = 0.4, x = list(relation = "free")),
            strip = strip.special,
            pch = "|",
            xlab = "",
            layout = c(4, 1),
            panel = panel.special)
graphics.off()
setwd(kCalcDir)
pdf(file = "smoking_performance.pdf", width = 6.6, height = 2.1)
plot(p)
dev.off()
file.copy(from = file.path(kCalcDir, "smoking_performance.pdf"),
          to = file.path(kPaperDir, "smoking_performance.pdf"),
          overwrite = TRUE)



p <- dplot(jitter(value) ~ income | age * sex, 
      data = finiteY(res, population = popn.all.reg)/popn.all,
      weights = popn.all.reg,
      par.settings = list(fontsize = list(text = 9)),
      scales = list(x = list(rot = 90), y = list(rot = 0)),
      overlay = list(values = smokes.all / popn.all, col = "black", lwd = 2))

graphics.off()
pdf(file = "Modelled vs actual - smoking simulation.pdf", paper = "a4r", h = 0, w = 0)
plot(p)
dev.off()



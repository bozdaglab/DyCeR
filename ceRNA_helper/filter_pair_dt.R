#### --- variables -------##########
pair.dt = pair.dt
cor.threshold = cor.threshold
adj.pvalue.threshold = adj.pvalue.threshold

#### --- analysis -------##########
filtered.pair.dt = pair.dt[pair.dt$cor > cor.threshold & pair.dt$partial.corr > 0 & pair.dt$sensitivity.cor > 0]
cat("   number of pair after correlation filtering:", nrow(filtered.pair.dt), "\n")

filtered.pair.dt$adj.cor.pvalue = p.adjust(p=filtered.pair.dt$cor.pvalue,method="BH")
filtered.pair.dt$adj.partial.cor.pvalue = p.adjust(p=filtered.pair.dt$partial.cor.pvalue,method="BH")
filtered.pair.dt$adj.hypergeometric.pvalue = p.adjust(p=filtered.pair.dt$hypergeometric.pvalue,method="BH")

filtered.pair.dt = filtered.pair.dt[filtered.pair.dt$adj.cor.pvalue < adj.pvalue.threshold 
                                    & filtered.pair.dt$adj.partial.cor.pvalue < adj.pvalue.threshold
                                    & filtered.pair.dt$adj.hypergeometric.pvalue < adj.pvalue.threshold]
cat("   number of pair after pvalue filtering:", nrow(filtered.pair.dt), "\n")

#### --- clean up -------##########

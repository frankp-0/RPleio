library(data.table)
library(ggplot2)
library(magrittr)

rpleio <- fread("../results/pgc.txt.gz")
pleio <- fread("../results/pgc_pleio.txt.gz")
vc <- data.table(rpstat = rpleio$pleio_stat,
                 pstat = pleio$pleio_stat)

png("../results/vc_compare.png")
vc %>% ggplot(aes(x = pstat, y = rpstat)) + geom_abline(slope = 1, intercept = 0) + geom_point() + ggtitle(expression(paste("Comparison of Test Statistics for ", tau^2, " > 0"))) + xlab(expression("S"["PLEIO"])) + ylab(expression("S"["RPleio"]))
dev.off()

ris <- fread("../results/pgc.is")
pis <- fread("../results/pgc_pleio.isf")
setorder(ris, theta)
setorder(pis, V1)
pis <- pis[order(ris$p),V2]
ris <- ris[order(p),p]
is <- data.table(pis = log(pis),
                 ris = log(ris))

png("../results/is_compare.png")
is %>% ggplot(aes(x = pis, y = ris)) + geom_abline(slope = 1, intercept = 0) + geom_point() + ggtitle("Comparison of Importance Sampling Log P-Values") + xlab(expression("Log P"["PLEIO"])) + ylab(expression("Log P"["RPleio"]))
dev.off()

p <- pleio$pleio_p
rp <- rpleio$pleio_p
p <- data.table(p = as.numeric(p), rp = rp)

png("../results/interp.png")
p %>% ggplot(aes(x = log(p), y = log(rp))) + geom_point() + geom_abline(slope = 1, intercept = 0) + ggtitle("Comparison of Interpolated Log P-values") + xlab("Log P (PLEIO)") + ylab("Log P (R_Pleio)")
dev.off()

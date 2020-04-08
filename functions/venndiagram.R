BiocManager::install("VennDiagram", version = "3.8")
library(VennDiagram)

grid.newpage()
draw.double(area1 = 1738, area2 = 3997, n12 = 1673, n13 = 4 category = c("African", "Other","test"), lty = "blank", 
                 fill = c("skyblue", "pink1","skyblue"))


?VennDiagram
?draw.triple.venn



BiocManager::install("VennDiagram", version = "3.8")
library(VennDiagram)

grid.newpage()
draw.triple.venn(area1 = 2109+1673+1738, area2 = 4040+1673+3997, area3 = 2053+3997+1738, n12 = 1673, n23 = 3997, n13 = 1738, 
                 n123 = 357, category = c("T Cells", "M Cells", "B Cells"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))3


draw.pairwise.venn(area1, area2, cross.area, category = rep("", 2),
                   euler.d = TRUE, scaled = TRUE, inverted = FALSE,
                   ext.text = TRUE, ext.percent = rep(0.05, 3), lwd =
                     rep(2, 2), lty = rep("solid", 2), col = rep("black",
                                                                 2), fill = NULL, alpha = rep(0.5, 2), label.col =
                     rep("black", 3), cex = rep(1, 3), fontface =
                     rep("plain", 3), fontfamily = rep("serif", 3), cat.pos
                   = c(-50, 50), cat.dist = rep(0.025, 2), cat.cex =
                     rep(1, 2), cat.col = rep("black", 2), cat.fontface =
                     rep("plain", 2), cat.fontfamily = rep("serif", 2),
                   cat.just = rep(list(c(0.5, 0.5)), 2), cat.default.pos
                   = "outer", cat.prompts = FALSE, ext.pos = rep(0, 2),
                   ext.dist = rep(0, 2), ext.line.lty = "solid",
                   ext.length = rep(0.95, 2), ext.line.lwd = 1,
                   rotation.degree = 0, rotation.centre = c(0.5, 0.5),
                   ind = TRUE, sep.dist = 0.05, offset = 0, cex.prop =
                     NULL, print.mode = "raw", sigdigs = 3, ...)

grid.newpage()
draw.pairwise.venn(area1 = 285, area2 = 145, cross.area = 134, category = c("African", "Other"), fill = c("red", "blue"))

draw.pairwise.venn(area1 = 285, area2 = 145, cross.area = 134)






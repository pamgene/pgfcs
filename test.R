library(dplyr)
library(reshape2)

X = matrix(nrow = 8, ncol =50, data = rnorm(8*50))
grp = factor(c(1,1,1,1,2,2,2,2))
pair = factor(c(1,2,3,4,1,2,3,4))
d = data.frame(X, grp, pair)
dm = melt(d, id.vars = c("grp", "pair"))
Mp = acast(dm, variable ~ pair ~ grp)
F = Mp[,,2]-Mp[,,1]


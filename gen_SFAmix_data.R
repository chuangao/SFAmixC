gen_SFAmix_data <- function (std = 1) {
    nf.s <- 10
    nf.d <- 5
    ng <- 10000
    ns <- 500
    nf <- nf.s + nf.d
    lams <- matrix(0, nrow = ng, ncol = nf.s)
    lamd <- matrix(rnorm(ng * nf.d, 0, std), nrow = ng, ncol = nf.d)
    ex <- matrix(rnorm(nf * ns, 0, 1), nrow = nf)
    block <- ng/nf.s
    for (i in 1:nf.s) {
        ne <- sample(50:100, 1)
        lams[(i - 1) * block + sample(1:block, ne, replace = F), 
            i] = rnorm(ne, 0, std)
    }
    err <- matrix(rnorm(ng * ns), nrow = ng, ncol = ns)
    lam <- as.matrix(cbind(lams, lamd))
    lam <- lam[, sample(1:nf, nf, replace = F)]
    y <- lam %*% ex + err
    return(list(y = y, lams = lams, lamd = lamd, ex = ex))
}

data <- gen_SFAmix_data(std=2)
system("mkdir -p ./data/")
write.table(data$y,"./data/gexp.txt",row.names=F,col.names=F,sep="\t",quote=F)


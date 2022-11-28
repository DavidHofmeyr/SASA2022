if(!"fastICA"%in%installed.packages()) install.packages("fastICA")
library(fastICA)

X0 <- cbind(images[,sample(1:8, 3)], rnorm(nrow(images)), runif(nrow(images)))
X <- X0%*%matrix(rnorm(25), 5, 5)

ica1 <- fastICA(X, 5)
ica2 <- fk_ICA(X, 5, nbin = 10000)
ica3 <- myica(X, 5)
pca <- list(S = X%*%eigen(cov(X))$vectors)

for(i in 1:5) ica1$S[,i] <- ica1$S[,i]*sign(cor(ica1$S[,i], X0)[which.max(abs(cor(ica1$S[,i], X0)))])
for(i in 1:5) ica2$S[,i] <- ica2$S[,i]*sign(cor(ica2$S[,i], X0)[which.max(abs(cor(ica2$S[,i], X0)))])
for(i in 1:5) ica3$S[,i] <- ica3$S[,i]*sign(cor(ica3$S[,i], X0)[which.max(abs(cor(ica3$S[,i], X0)))])
for(i in 1:5) pca$S[,i] <- pca$S[,i]*sign(cor(pca$S[,i], X0)[which.max(abs(cor(pca$S[,i], X0)))])


par(mfrow = c(6, 5))
par(mar = c(0, 0, 0, 0))
for(j in 1:5) plotimage(X0[,j])
for(j in 1:5) plotimage(X[,j])
for(j in 1:5) plotimage(ica1$S[,j])
for(j in 1:5) plotimage(ica2$S[,j])
for(j in 1:5) plotimage(ica3$S[,j])
for(j in 1:5) plotimage(pca$S[,j])


sum(apply(ica1$S, 2, function(x) max(abs(cor(x, X0)))))
sum(apply(ica2$S, 2, function(x) max(abs(cor(x, X0)))))
sum(apply(ica3$S, 2, function(x) max(abs(cor(x, X0)))))
sum(apply(pca$S, 2, function(x) max(abs(cor(x, X0)))))

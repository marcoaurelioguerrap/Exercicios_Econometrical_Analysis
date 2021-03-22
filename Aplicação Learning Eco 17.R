library(Matrix)

setwd("C:/Marco/Mestrado/Matérias/Econometria I/exercicios")

dados = read.csv( "Tabela 17.csv", sep = ";", dec=',')

lnq1 = log(dados[,5])
lnq2 = log(dados[,6])
lnq3 = log(dados[,7])

lnp1 = log(dados[,1])
lnp2 = log(dados[,2])
lnp3 = log(dados[,3])

lny = log(dados[,4])

X1 = cbind(1,lnp1 , lny)
X2 = cbind(1,lnp2 , lny)
X3 = cbind(1,lnp3 , lny)
x = cbind(X1,X2,X3)
colnames(X1) <- c('const','lnp1','lny')
colnames(X2) <- c('const','lnp2','lny')
colnames(X3) <- c('const','lnp3','lny')

y = c(lnq1,lnq2,lnq3)
X = matrix(data = 0, nrow=90 , ncol=9)

for (j in c(0,3,6)){
  for (i in 1:30){
    X[i+(j*10),1:3+j] = x[i,1:3+j]
  }
}

b_sem_sur = solve(crossprod(X)) %*% t(X) %*% y
Y_est_sem_sur = X %*% solve(crossprod(X)) %*% (t(X) %*% y)
e= y - Y_est_sem_sur
SSE_sem = sum(e*e)
SSE_sur = matrix(data = NA, nrow=3 , ncol=3)

for (j in c(0,3,6)){
  for (i in c(0,3,6)){
    SSE_sur[(i/3)+1,(j/3)+1] = sum(e[1:30+i*10]*e[1:30+j*10])
  }
}

SSE_sem_teste = e %*% t(e)
S_2_sem = SSE_sem /(nrow(X)-ncol(X))
W = S_2_sem * solve(crossprod(X))
b_sem_sur = solve(t(X)%*%solve(W)%*%X) %*% t(X) %*%solve(W)%*% y

cov.wt(X)

b_ols = matrix(data = NA, nrow=3 , ncol=3)
colnames(b_ols) <- c('const','lnp','lny')
Y_est = matrix(data = NA, nrow=30 , ncol=3)
SSE = matrix(data = NA, nrow=3 , ncol=3)
S_2 = matrix(data = NA, nrow=3 , ncol=3)
VarCov = array(NA, c(3,3,9))
colnames(VarCov) <- c('const','lnp','lny')
rownames(VarCov) <- c('const','lnp','lny')
Resultado = array(NA, c(3,4,3),dimnames=list(c( 'const','lnp','lny'),
                                             c( 'Coef','DesPad','b/DP','P-valor' ),
                                             c( '1','2','3')))
a=0
for (i in 1:3){
  b_ols[i,] = solve(t(X[[i]]) %*% X[[i]]) %*% (t(X[[i]]) %*% Y[,i])
  Y_est[,i] = X[[i]] %*% solve(crossprod(X[[i]])) %*% (t(X[[i]]) %*% Y[,i])
  e = Y - Y_est #erro ta aqui
}
for(i in 1:3){
  for(j in 1:3){
    SSE[i,j] = sum(e[,i]*e[,j])
    S_2[i,j] = SSE[i,j]/ (nrow(X[[i]])-ncol(X[[i]]))
    a=a+1
    print(a)
    VarCov[,,a] = S_2[i,j] * solve(t(X[[i]])%*%X[[j]])
  }
  Resultado[,,i] = cbind(b_ols[i,],diag(VarCov[,,i]^(1/2)),b_ols[i,]/diag(VarCov[,,i]^(1/2)),2*pt(-abs(b_ols[i,]/diag(VarCov[,,i]^(1/2))),df = nrow(dados)))
}
print(Resultado)
print(VarCov)

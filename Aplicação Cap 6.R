library(Matrix)

setwd("C:/Marco/Mestrado/Matérias/Econometria I/exercicios")

dados = read.csv( "TableF4-4.csv", sep = ",")

#define a matriz de variaveis
Y = matrix(data = NA, nrow=nrow(dados) , ncol= 1 )
X = matrix(data = NA, nrow=nrow(dados) , ncol=15)
X_rest = matrix(data = NA, nrow=nrow(dados) , ncol=10)
# ln C = B[0] + B[1]*ln Q + G[0]*ln Pk + G[1]*ln Pl + G[2]*ln Pf +
Y = log(dados$COST)
X[,1] = 1
X[,2] = log(dados$Q)
X[,3] = log(dados$PK)
X[,4] = log(dados$PL)
X[,5] = log(dados$PF)
#         + ??kk*(1/2)*(ln Pk)^2 + ??ll*(1/2)*(ln Pl)^2 + ??ff*(1/2)*(ln Pf)^2
# pode ta faltando fazer o log em tudo
X[,6] = (1/2)*(log(dados$PK)^2)
X[,7] = (1/2)*(log(dados$PL)^2)
X[,8] = (1/2)*(log(dados$PF)^2)
#         + ??kl*(ln Pk)*(ln Pl ) + ??kf*(ln Pk)*(ln Pf ) + ??lf*(ln Pl )*( ln Pf )
X[,9] = (log(dados$PK))*(log(dados$PL))
X[,10] = (log(dados$PK))*(log(dados$PF))
X[,11] = (log(dados$PL))*(log(dados$PF))
#         + ?? (1/2)*(ln Q)^2
X[,12] = (1/2)*(log(dados$Q)^2)
#         + ??Qk*(ln Q)*(ln Pk) + ??Ql*(ln Q)*(ln Pl ) + ??Qf (ln Q)*(ln Pf ) + ??.
X[,13] = (log(dados$Q))*(log(dados$PK))
X[,14] = (log(dados$Q))*(log(dados$PL))
X[,15] = (log(dados$Q))*(log(dados$PF))
print('A) Write out the R matrix and q vector in (5-8) that are needed to impose the restriction of linear homogeneity in prices.')

#Restrição
q= matrix(data=0 , ncol=1,nrow=5)
R = matrix(data=0 , ncol=15,nrow=5)
#??k + ??l + ?? f = 1, 
R[1,3]=1
R[1,4]=1
R[1,5]=1
q[1]=1
#??kk + ??kl + ??kf = 0
R[2,6]=1
R[2,9]=1
R[2,10]=1
q[2]=0
#??kl + ??ll + ??lf = 0
R[3,7]=1
R[3,9]=1
R[3,11]=1
q[3]=0
#??kf + ??lf + ??ff = 0
R[4,10]=1
R[4,11]=1
R[4,8]=1
q[4]=0
#??QK + ??Ql + ??Qf = 0.
R[5,13]=1
R[5,14]=1
R[5,15]=1
q[5]=0
print('R')
print(R)
print('q')
print(q)

print('B)')

# teste F ( 5-16) , nao deu certo 
B = solve(crossprod(X)) %*% t(X) %*% Y
M_1 = diag(dim((X %*% solve(crossprod(X))) %*% t(X))[1]) - ((X %*% solve(crossprod(X))) %*% t(X))
e =  M_1 %*% Y
AsyCov = as.numeric((t(e) %*% (e))/(nrow(X)-ncol(X))) * solve(crossprod(X))
F = (t(R %*% B - q) %*% solve( R %*% AsyCov %*% t(R) ) %*% (R %*% B - q))/(nrow(X)-ncol(X))


# modelo restrito
Y_rest = log(dados$COST / dados$PF)
X_rest[,1] = 1
X_rest[,2] = log(dados$Q)
X_rest[,3] = log(dados$PK/dados$PF)
X_rest[,4] = log(dados$PL/dados$PF)
X_rest[,5] = (1/2)*(log(dados$PK/dados$PF )^2)
X_rest[,6] = (1/2)*(log(dados$PL/dados$PF )^2)
X_rest[,7] = log(dados$PK/dados$PF)*log(dados$PL/dados$PF)
X_rest[,8] = (1/2)*(log(dados$Q )^2)
X_rest[,9] = log(dados$Q)*log(dados$PK / dados$PF)
X_rest[,10] = log(dados$Q)*log(dados$PL/ dados$PF)
B_rest = solve(crossprod(X_rest)) %*% t(X_rest) %*% Y_rest
M_1_rest = diag(dim((X_rest %*% solve(crossprod(X_rest))) %*% t(X_rest))[1]) - ((X_rest %*% solve(crossprod(X_rest))) %*% t(X_rest))
e_rest =  M_1_rest %*% Y_rest
#colnames(X) = c ('cons','Q','0.5*ln(Q)^2','PK','PL','PF')

#To conform to the underlying theory of production, it is necessary
#to impose the restriction that the cost function be homogeneous of degree one in
#the three prices. This is done with the restriction Bk + Bl + Bf = 1, or Bf = 1 ??? Bk ??? Bl .

Y_rest = log(dados$co / dados$PF)
X_rest[,1] = 1
X_rest[,2] = log(dados$Q)
X_rest[,3] = (1/2)*(log(dados$Q)^2)
X_rest[,4] = log(dados$PK / dados$PF)
X_rest[,5] = log(dados$PL / dados$PF)

colnames(X_rest) = c ('cons','Q','0.5*ln(Q)^2','ln(PK/PF)','ln(PL/PF)')


print( 'A) Using all 158 observations, compute the estimates of the parameters
       in the cost function and the estimate of the asymptotic covariance matrix.' )
B_rest = solve(crossprod(X_rest)) %*% (t(X_rest) %*% Y_rest)
rownames(B_rest) = colnames(X_rest)
colnames(B_rest) = ( 'Coef' )
print(B_rest)
#asymptotic covariance matrix
M_1 = diag(dim((X_rest %*% solve(crossprod(X_rest))) %*% t(X_rest))[1]) - ((X_rest %*% solve(crossprod(X_rest))) %*% t(X_rest))
e_rest =  M_1 %*% Y_rest
print ('AsyCov = (tee/n-k) * (tXX)^-1')

AsyCov = as.numeric((t(e_rest) %*% (e_rest))/(nrow(X_rest)-ncol(X_rest))) * solve(crossprod(X_rest))
colnames(AsyCov) <- NULL
rownames(AsyCov) <- NULL
print (AsyCov)


print( 'B) Note that the cost function does not provide a direct estimate of b_f . Compute
       this estimate from your regression results, and estimate the asymptotic standard error.' )
b_f = 1- B_rest[4]- B_rest[5]
EP = (AsyCov[4,4]+AsyCov[5,5]+2*AsyCov[4,5])^(1/2)
est = b_f / EP
print('b_f:')
print(b_f)
print('Erro Padrão:')
print(EP)
print('Est:')
print(est)

print( 'C) Compute an estimate of Q??? using your regression results and then form a confidence
       interval for the estimated efficient scale.' )

Q_2 = exp((1-B_rest[2])/B_rest[3])
#Deterina a Jacobiana da regressão f(b[2],b[3]) = exp((1-B_rest[2])/B_rest[3])
Jac_1 = - (exp((1-B_rest[2])/B_rest[3])/B_rest[3])
Jac_2 = -((1-B_rest[2])*exp((1-B_rest[2])/B_rest[3]))/(B_rest[3]^2)
Jac = cbind(Jac_1,Jac_2)

V=AsyCov[2:3,2:3]
B = cbind(B_rest[2],B_rest[3])
#Calcula a variancia aproximada de f(b[2],b[3])
V0= Jac %*% V %*% t(Jac)

Result = cbind(Q_2,V0^(1/2) , Q_2 + 1.96*V0^(1/2),Q_2 - 1.96*V0^(1/2))
colnames(Result)= c('Média','Des.Pad','Int.Inf','Int.Sup')
print(Result)
print('D) Examine the raw data and determine where in the sample the efficient scale lies.
That is, determine how many firms in the sample have reached this scale, and
whether, in your opinion, this scale is large in relation to the sizes of firms in
the sample.')

c=0

for (i in 1:158)
{
  if((dados$Q[i]>=Q_2))
  {c = c + 1 }else{next}
}

print(paste('A escala eficiente de produção é ',round(Q_2, digits=2),'. São no total',c,'firmas que produzem acima da escala eficiente. Considerando que são 158 firmas, são poucas as empresas que alcançaram esse patamar de eficiência.'))

library(Matrix)

setwd("C:/Marco/Mestrado/Matérias/Econometria I/exercicios")

dados = read.csv( "TableF4-4.csv", sep = ",")

#define a matriz de variaveis
Y = matrix(data = NA, nrow=nrow(dados) , ncol= 1 )
X = matrix(data = NA, nrow=nrow(dados) , ncol=6)
X_rest = matrix(data = NA, nrow=nrow(dados) , ncol=5)

Y = dados$COST
X[,1] = 1
X[,2] = dados$Q
X[,3] = (1/2)*(log(dados$Q)^2)
X[,4] = dados$PK
X[,5] = dados$PL
X[,6] = dados$PF
colnames(X) = c ('cons','Q','0.5*ln(Q)^2','PK','PL','PF')

#To conform to the underlying theory of production, it is necessary
#to impose the restriction that the cost function be homogeneous of degree one in
#the three prices. This is done with the restriction Bk + Bl + Bf = 1, or Bf = 1 ??? Bk ??? Bl .

Y_rest = log(Y / dados$PF)
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


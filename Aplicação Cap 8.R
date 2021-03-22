library(Matrix)

setwd("C:/Marco/Mestrado/Matérias/Econometria I/exercicios")

dados = read.csv( "TableF8-1.csv", sep = ",")

attach(dados)


#define a matriz de variaveis
Y = matrix(data = NA, nrow=nrow(dados) , ncol= 1 )
X = matrix(data = NA, nrow=nrow(dados) , ncol=10)
Z = matrix(data = NA, nrow=nrow(dados) , ncol=11)


Y = LWAGE
colnames(X) <- c('const','Ed','Exp' , 'Exp^2','Occ','Ind', 'South','SMSA','Blk','WKS')
colnames(Z) <- c('const','Ed','Exp' , 'Exp^2','Occ','Ind', 'South','SMSA','Blk','UNION','Fem')
X[,1] = 1
X[,2] = ED
X[,3] = EXP
X[,4] = EXP * EXP
X[,5] = OCC
X[,6] = IND
X[,7] = SOUTH
X[,8] = SMSA
X[,9] = BLK
X[,10] = WKS

for (i in 1:9){
  Z[,i] = X[,i]
}
Z[,10] = UNION
Z[,11] = FEM

#OLS
b_ols = solve(t(X) %*% X) %*% (t(X) %*% Y)
#2SLS
PX = Z %*% solve(crossprod(Z)) %*% (t(Z) %*% X)
b_2sls = solve(t(PX) %*% X) %*% (t(PX) %*% Y)
v = X - PX
#A) Qual é a função de demanda?
print('A) Qual é a função de demanda?')
print('LWage = B[1] + B[2]*Ed + B[3]*Exp + B[4]*Exp^2 + B[5]*Occ + B[6]*Ind + B[7]*South + B[8]*SMSA + B[9]*Blk + B[10]*WKS + e')

#B) Compare o resultado da estimação por OLS e 2SLS

Y_est = X %*% solve(crossprod(X)) %*% (t(X) %*% Y)
e = Y - Y_est
SSE = sum(e*e)
S_2 = SSE/ (nrow(X)-ncol(X))
DesPad = S_2^(1/2)
VarCov = S_2 * solve(crossprod(X))
DesPadM = S_2^(1/2) * solve(crossprod(X))
Resultado = cbind(b_ols,diag(VarCov^(1/2)),b_ols/diag(VarCov^(1/2)),2*pt(-abs(b_ols/diag(VarCov^(1/2))),df = nrow(dados)))
colnames(Resultado) =c( 'Coef','DesPad','b/DP','P-valor' )
print(round(Resultado,digits = 6))

Y_est_2sls = PX %*% solve(crossprod(PX)) %*% (t(PX) %*% Y)
e_2sls = Y - Y_est_2sls
SSE_2sls = sum(e_2sls*e_2sls)
S_2_2sls = SSE_2sls/ (nrow(PX)-ncol(PX))
DesPad_2sls = S_2_2sls^(1/2)
VarCov_2sls = S_2_2sls * solve(crossprod(PX))
DesPadM_2sls = S_2_2sls^(1/2) * solve(crossprod(PX))
Resultado_2sls = cbind(b_2sls,diag(VarCov_2sls^(1/2)),b_2sls/diag(VarCov_2sls^(1/2)),2*pt(-abs(b_2sls/diag(VarCov_2sls^(1/2))),df = nrow(dados)))
colnames(Resultado_2sls) =c( 'Coef','DesPad','b/DP','P-valor' )
print(round(Resultado_2sls,digits = 6))

#C) os instrumentos utilizados na equação são relevantes?
print('C) os instrumentos utilizados na equação são relevantes?')

print('A hipótese nula é de que os instrumentos são não significativos. Isso é feito na primeira regressão dos modelo de dois estágios. o teste usado é o F de significância conjunta do modelo. Pelo teste F é possível notar que os instrumentos são válidos.')

b_wks = solve(t(Z) %*% Z) %*% (t(Z) %*% X[,10])
X_est_wks = Z %*% solve(crossprod(Z)) %*% (t(Z) %*% X[,10])
e_wks = X[,10] - X_est_wks
SSE_wks = sum(e_wks*e_wks)
S_2_wks = SSE_wks/ (nrow(Z)-ncol(Z))
DesPad_wks = S_2_wks^(1/2)
VarCov_wks = S_2_wks * solve(crossprod(Z))
DesPadM_wks = S_2_wks^(1/2) * solve(crossprod(Z))
Resultado_wks = cbind(b_wks,diag(VarCov_wks^(1/2)),b_wks/diag(VarCov_wks^(1/2)),2*pt(-abs(b_wks/diag(VarCov_wks^(1/2))),df = nrow(dados)))
colnames(Resultado_wks) =c( 'Coef','DesPad','b/DP','P-valor' )
print(round(Resultado_wks,digits = 6))


R = matrix(0, nrow=10,ncol=11)
for (i in 1:10){
R[i,1+i] = 1
}


q = matrix(0, nrow=10,ncol=1)
F = (t(R %*% b_wks - q) %*% solve(R %*% (S_2_wks * solve(crossprod(Z))) %*% t(R)) %*% (R %*% b_wks - q))/(ncol(Z)-1)
teste_F = cbind( F,df(F, df1=10,df2=(nrow(dados)-ncol(Z))))
colnames(teste_F) =cbind( 'Teste F','P-valor' )
print(round(teste_F,digits=4))


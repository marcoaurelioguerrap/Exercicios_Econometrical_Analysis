library(Matrix)

setwd("C:/Marco/Mestrado/Matérias/Econometria I/exercicios")

dados = read.csv( "Koop-Tobias.csv", sep = ",")

#dados = dados_raw[!duplicated(dados_raw$PERSONID) & !dados_raw$PERSONID > 15,]

#define a matriz de variaveis
Y = matrix(data = NA, nrow=nrow(dados) , ncol= 1 )
X = matrix(data = NA, nrow=nrow(dados) , ncol=ncol(dados))
X_1 = matrix(data = NA, nrow=nrow(dados) , ncol=4)
X_2 = matrix(data = NA, nrow=nrow(dados) , ncol=4)

#Ed0 = matrix(data = 0, nrow = 31 , ncol = 1 )
#Ed = matrix(data = 0, nrow = 31 , ncol = 1 )

Y = dados$LOGWAGE
#X1 = [constant, education, experience, ability]
X_1[,1] = 1
X_1[,2] = dados$EDUC
X_1[,3] = dados$POTEXPER
X_1[,4] = dados$ABILITY
#X2 = [mother's education,father's education, broken home, number of siblings]
X_2[,1] = dados$MOTHERED
X_2[,2] = dados$FATHERED
X_2[,3] = dados$BRKNHOME
X_2[,4] = dados$SIBLINGS
X = cbind( X_1 , X_2 )

print( 'A) Compute the full regression of log wage on X1 and X2 and report all results.' )

B = solve(crossprod(X)) %*% (t(X) %*% Y)
Y_est = t(B) %*% t(X)
e = Y - Y_est
SSE = sum(e*e)
S_2 = SSE/ (nrow(dados)-ncol(dados)) 
DesPad = S_2^(1/2)
VarCov = S_2 * solve(crossprod(X))
Resultado = cbind(B,diag(VarCov^(1/2)),B/diag(VarCov^(1/2)),round(2*pt(-abs(B/diag(VarCov^(1/2))),df = nrow(dados)),digits=4))

rownames(Resultado) = c ('Cons','Educ','Exp','Habi','Educ_Mae','Educ_Pai','Broken Home','Irmaos')
colnames(Resultado) =c( 'Coef','DesPad','B/DP','P-valor' )
print(round(Resultado,digits = 4))

print( 'B). Use an F test to test the hypothesis that all coefficients except the constant term are zero.' )

R = diag(8)
R[1,1] = 0
q = matrix(0, nrow=8,ncol=1)
F = (t(R %*% B - q)[2:8] %*% solve(R[2:8,2:8] %*% (S_2 * solve(crossprod(X)))[2:8,2:8] %*% t(R)[2:8,2:8]) %*% (R %*% B - q)[2:8])/(ncol(X)-1)
teste_F = cbind( F,df(F, df1=7,df2=(nrow(dados)-ncol(dados))))
colnames(teste_F) =cbind( 'Teste F','P-valor' )
print(round(teste_F,digits=4))
print( 'C) Use an F statistic to test the joint hypothesis that the coefficients on the four household variables in X2 are zero.' )

F_X2 = (t(R %*% B - q)[5:8] %*% solve(R[5:8,5:8] %*% (S_2 * solve(crossprod(X)))[5:8,5:8] %*% t(R)[5:8,5:8]) %*% (R %*% B - q)[5:8])/(4)
teste_F_X2 = cbind( F_X2,df(F_X2, df1=7,df2=(nrow(dados)-ncol(dados))))
colnames(teste_F_X2) =cbind( 'Teste F','P-valor' )
print(round(teste_F_X2,digits=4))

print("D) Use a Wald test to carry out the test in part c.")

W = t(B[5:8])%*% solve(R[5:8,5:8] %*% (S_2 * solve(crossprod(X)))[5:8,5:8] %*% t(R)[5:8,5:8]) %*% B[5:8]
colnames(W) = ('Wald Test')
print(round(W,digits=4))



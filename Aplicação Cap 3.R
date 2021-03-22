library(Matrix)

setwd("C:/Marco/Mestrado/Matérias/Econometria I/exercicios")

dados_raw = read.csv( "Koop-Tobias.csv", sep = ",")

dados = dados_raw[!duplicated(dados_raw$PERSONID) & !dados_raw$PERSONID > 15,]

#define a matriz de variaveis
Y = matrix(data = NA, nrow=nrow(dados) , ncol= 1 )
X = matrix(data = NA, nrow=nrow(dados) , ncol=ncol(dados))
X_1 = matrix(data = NA, nrow=nrow(dados) , ncol=4)
X_2 = matrix(data = NA, nrow=nrow(dados) , ncol=3)

#Ed0 = matrix(data = 0, nrow = 31 , ncol = 1 )
#Ed = matrix(data = 0, nrow = 31 , ncol = 1 )

Y = dados$LOGWAGE
X_1[,1] = 1
X_1[,2] = dados$EDUC
X_1[,3] = dados$POTEXPER
X_1[,4] = dados$ABILITY
X_2[,1] = dados$MOTHERED
X_2[,2] = dados$FATHERED
X_2[,3] = dados$SIBLINGS
X = cbind( X_1 , X_2 )
print( 'A) Coeficientes de MQO da regressão de Y em X_1' )
B_1 = solve(crossprod(X_1)) %*% (t(X_1) %*% Y)
rownames(B_1) = c ('Cons','Educ','Exp','Habi')
colnames(B_1) =( 'Coef' )
print(B_1)
print( 'B) Coeficientes de MQO da regressão de Y em X_1 e X_2 (X = X_1 e X_2)' )
B = solve(crossprod(X)) %*% (t(X) %*% Y)
rownames(B) = c ('Cons','Educ','Exp','Habi','Educ_Mae','Educ_Pai','Irmaos')
colnames(B) =( 'Coef' )
print(B)
print( 'C) Coeficientes de MQO da regressão de X_2[k] em X_1 e , k = { Educ_Mae, Educ_Pa, Irmaos }' )

M_1 = diag(dim((X_1 %*% solve(crossprod(X_1))) %*% t(X_1))[1]) - ((X_1 %*% solve(crossprod(X_1))) %*% t(X_1))

Medias_X2_em_X1 = matrix(data = NA , nrow = 1, ncol = 3)
rownames(Medias_X2_em_X1) = c ('Medias')
colnames(Medias_X2_em_X1) = c ('Educ_Mae','Educ_Pai','Irmaos')
for (k in 1:3){
  Medias_X2_em_X1[,k] = mean( M_1 %*% X_2[,k])
}
print(Medias_X2_em_X1)
print( 'D) Calcule R2 para regressão de y em X_1 e X_2. Repita o calculo para o caso de se omitir o termo constante de X_1, O que acontece com  R2?')

i = matrix(data = 1, nrow=nrow(dados) , ncol= 1)

M_0 = diag(dim(i)[1]) - (i %*% solve(crossprod(i)) %*% t(i))

R2 = (t(B) %*% t(X) %*% M_0 %*% X %*% B) / (t(Y) %*% M_0 %*% Y)
X_sem_cons = X[,-1]
B_sem_cons = solve(crossprod(X_sem_cons)) %*% (t(X_sem_cons) %*% Y)
R2_sem_cons = (t(B_sem_cons) %*% t(X_sem_cons) %*% M_0 %*% X_sem_cons %*% B_sem_cons) / (t(Y) %*% M_0 %*% Y)

R_quadrados = cbind(R2,R2_sem_cons)
rownames(R_quadrados) = c ('R2')
colnames(R_quadrados) = c ( 'Com const', 'Sem const')
print( R_quadrados)

print('E) Calcule o R2 ajustado.')

R2_ajus = 1 - (( 1 - R2 ) * ( nrow(X) - 1 ) / ( nrow(X)- ncol(X)))
print(R2_ajus)

print('')

print('F) y em X_1 e M_1 %*% X_2 ')

X_f = cbind( X_1 , M_1 %*% X_2 )
B_f = solve(crossprod(X_f)) %*% (t(X_f) %*% Y)
rownames(B_f) = c ('Cons','Educ','Exp','Habi','Educ_Mae','Educ_Pai','Irmaos')
colnames(B_f) =( 'Coef' )
print(B_f)

library(Matrix)

options(warn = -1)

setwd("C:/Marco/Mestrado/Matérias/Econometria I/exercicios")

dados = read.csv( "Koop-Tobias.csv", sep = ",")

#define a matriz de variaveis
Y = matrix(data = NA, nrow=nrow(dados) , ncol= 1 )
X = matrix(data = NA, nrow=nrow(dados) , ncol=ncol(dados))
X_1 = matrix(data = NA, nrow=nrow(dados) , ncol=4)
X_2 = matrix(data = NA, nrow=nrow(dados) , ncol=4)
X_3 = matrix(data = NA, nrow=nrow(dados) , ncol=2)
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
#X_3 = [grad, col ]
#dados = dados_raw[!duplicated(dados_raw$PERSONID) & !dados_raw$PERSONID > 15,]

X1 = cbind( X_1 , X_2[,-3] )

print( 'A) Compute o regressão de minimos quadrados. Qual é o valor mariginal, em $/hrs , de uma adição de um ano de educação para alguem com 12 anos quando todas as outras variaveis estão na média e BRKHOME = 0')

B1 = solve(crossprod(X1)) %*% (t(X1) %*% Y)


print( paste('Salario_Marginal[educação| educ = 12 & Brkhome = 0] =',round(B1[2],digits=7),'$/hrs'))

print(exp(t(colMeans(X1)) %*% B1)*B1[2])

print( 'B)Substitua Educ por (Col, Grad), Qual é o valor marginal de uma graduação no colegio?' )

hist(dados$EDUC,breaks = seq.int(8,20,1))

X_3 = matrix(data = NA, nrow=nrow(dados) , ncol=2)
X_3[,1] =  as.numeric(dados$EDUC > 12 & dados$EDUC <=16)
X_3[,2] =  as.numeric(dados$EDUC > 16) 

X_D = cbind( X_1[,1], X_3, X_1[,-1:-2] , X_2  )

B_D = solve(crossprod(X_D)) %*% (t(X_D) %*% Y)
rownames(B_D) <- c( 'const','Col','Grad','PExp','Habil','Mae','Pai','BRKHOME','Irmaos')
colnames(B_D) <- c('Coef')
print(B_D)
print( paste('variação percentual do Salario[educação| diploma colegio] =',round(B_D[2],digits=7),'$/hrs'))

print( 'C)Educação ao quadrado' )

X2 = cbind( X_1 , X_2 ,X_1[,2]^2 )

B_D2 = solve(crossprod(X2)) %*% (t(X2) %*% Y)
Y_est = t(B_D2) %*% t(X2)
e = Y - Y_est
SSE = sum(e*e)
S_2 = SSE/ (nrow(X2)-ncol(X2)) 
DesPad = S_2^(1/2)
VarCov = S_2 * solve(crossprod(X2))
Resultado = cbind(B_D2,diag(VarCov^(1/2)),B_D2/diag(VarCov^(1/2)),round(2*pt(-abs(B_D2/diag(VarCov^(1/2))),df = nrow(dados)),digits=4))

rownames(Resultado) = c ('Cons','Educ','Exp','Habi','Educ_Mae','Educ_Pai','Broken Home','Irmaos','Educ^2')
colnames(Resultado) =c( 'Coef','DesPad','B/DP','P-valor' )
print(round(Resultado,digits = 4))

print(paste('H0: B[Educ^2] = 0 é rejeitado. B[Educ^2] é diferente de zero com um nivel de significancia de  1%'))


a = t(colMeans(X2)) %*% B_D2
curve(a + B_D2[2]*x + B_D2[9]*x^2,0,20 , ylab = 'LogWage',xlab = 'Schooling',main="LogWage X Schooling")


print('D) Levando em conta a interação de Educ x Habilidade, qual é o valor adcional de um ano de estudos? Compute um intervalo de confiança para o impacto marinal no ln do salario dado um ano a mais na educação')

X_4 = matrix(data = NA, nrow=nrow(dados) , ncol=1)
X_4[,1] =  dados$EDUC * dados$ABILITY

X3 = cbind( X_1 , X_2 , X_4)

B3 = solve(crossprod(X3)) %*% (t(X3) %*% Y)

Y_est = t(B3) %*% t(X3)
e = Y - Y_est
SSE = sum(e*e)
S_2 = SSE/ (nrow(X3)-ncol(X3)) 
DesPad = S_2^(1/2)
VarCov = S_2 * solve(crossprod(X3))
Resultado = cbind(B3,diag(VarCov^(1/2)),B3/diag(VarCov^(1/2)),round(2*pt(-abs(B3/diag(VarCov^(1/2))),df = nrow(dados)),digits=4))

rownames(Resultado) = c ('Cons','Educ','Exp','Habi','Educ_Mae','Educ_Pai','Broken Home','Irmaos','Educ*Habil')
colnames(Resultado) =c( 'Coef','DesPad','B/DP','P-valor' )
print(round(Resultado,digits = 4))

media2 = B3[2] + B3[9]*mean(dados$ABILITY)

EP = (VarCov[2,2]+(mean(dados$ABILITY)^2)*VarCov[9,9] + 2*(mean(dados$ABILITY))*VarCov[2,9])^(1/2)

lim_Inf = media2 - 1.96 * EP

lim_Sup = media2 + 1.96 * EP

resultado = cbind(media2, lim_Inf,lim_Sup)
rownames(resultado)= c('intervalo de conf.')
print(resultado)

print("E) regressão c) e d)")

X4 = cbind( X_1 , X_2 , X_4 , X_1[,2]^2 )

B4 = solve(crossprod(X4)) %*% (t(X4) %*% Y)

Y_est = t(B4) %*% t(X4)
e = Y - Y_est
SSE = sum(e*e)
S_2 = SSE/ (nrow(X4)-ncol(X4)) 
DesPad = S_2^(1/2)
VarCov = S_2 * solve(crossprod(X4))
Resultado = cbind(B4,diag(VarCov^(1/2)),B4/diag(VarCov^(1/2)),round(2*pt(-abs(B4/diag(VarCov^(1/2))),df = nrow(dados)),digits=4))

rownames(Resultado) = c ('Cons','Educ','Exp','Habi','Educ_Mae','Educ_Pai','Broken Home','Irmaos','Educ*Habil','Educ^2')
colnames(Resultado) =c( 'Coef','DesPad','B/DP','P-valor' )
print(round(Resultado,digits = 4))

media_baixo = mean(dados$ABILITY[dados$ABILITY < mean(dados$ABILITY)])
media_alto = mean(dados$ABILITY[dados$ABILITY > mean(dados$ABILITY)])
lim_Inf_baixo = media_baixo - 1.96 * EP
lim_Sup_baixo = media_baixo + 1.96 * EP
lim_Inf_alto = media_alto - 1.96 * EP
lim_Sup_alto = media_alto + 1.96 * EP

a2 = t(colMeans(X4[,c(-2,-4,-9,-10)])) %*% B4[c(-2,-4,-9,-10)]
al = a2 + media_baixo*B4[2]
ah = a2 + media_alto*B4[2]
curve(ah + B4[2]*x + B4[10]*x^2 + B4[9]*media_alto*x,7,22,ylim=c(1.2,3) 
      , ylab = 'LogWage',xlab = 'Schooling', col="black")
curve(al + B4[2]*x + B4[10]*x^2 + B4[9]*media_baixo*x,7,22
      ,ylab = 'LogWage',xlab = 'Schooling',add=TRUE,lty=2)
legend(8,3, legend=c("media alta", "media baixa"),
       lty=1:2, cex=0.5)

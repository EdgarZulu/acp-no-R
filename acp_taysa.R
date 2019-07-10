rm(list=ls())
ls()
setwd("E:/Documentos/_TRABALHO/IFMT/MESTRADO/Tayza")
dados2=read.table("dados_acp2.txt",h=T,dec=",")
dados=dados2[,c(-1)]
cor(dados)
which((is.na(dados)))

# dados=scale(dados)
dados
str(dados)
attach(dados)
attach(dados2)
########################################
# install.packages("Factoshiny")
require(Factoshiny)
PCAshiny(dados)




##########################################################
pairs(dados,lower.panel = NULL,
      col=as.numeric(dados2$var),pch=7,font.labels=2)

dados.acp=princomp(dados2[,-1])
biplot(dados.acp)

require(vegan)
m1.rda=rda(dados2[,-1])
biplot(m1.rda,
       type=c("text","points"),
       display=c("sites","species"))

spp.names=levels(dados2$var)

ordihull(m1.rda,
         groups=dados2$var)
legend("bottomright",
       col=c(3,4,5),
       lty=1,legend=spp.names)



## matriz de correla??o
matcor=cor(dados)
print(matcor,digits=3)

acpcor=prcomp(dados,scale=T,retx=T)
summary(acpcor)
sum(acpcor$sdev^2)



plot(1:ncol(dados), acpcor$sdev^2, type = "b", xlab = "Componente",
     ylab = "Vari?ncia", pch = 20, cex.axis = 1.3, cex.lab = 1.3)

acpcor$rotation[, 1:2]

print(acpcor$sdev[1:2] * t(acpcor$rotation[, 1:2]), digits = 3)

# par(mfrow=c(1,2))
biplot(acpcor, xlab = "CP1", ylab = "CP2",cex.lab = 1.5, cex.axis = 1.5)

#########################################################################
# install.packages("vegan")

require(vegan)
resultado<-prcomp(dados)
resultado
summary(resultado)

# install.packages("xtable")
require(xtable)
xtable(resultado)

resultado$x
biplot(resultado)

### cruster hier?rquico
library(vegan)

tratamento.tipo<-vegdist(dados,"euclidean")
### outros metodos ## ward, single, complete, averange, mcquitty, median, centroid
hclust(tratamento.tipo,method="single")
cluster<-hclust(tratamento.tipo,method="complete")
plot(cluster)


### cruster aglomerativo
tratamentos<-vegdist(dados,"euclidean")
clusterUPGMA<-hclust(tratamentos,method="average") 
###  outros metodos ## centrid, mcquitty, median, averange
plot(clusterUPGMA)

########################################
# install.packages(c("factoextra","FactoMineR","ggplot2","ggcorrplot","psych"))
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(ggcorrplot)
library(psych)

boxplot(dados)
boxplot(scale(dados)) ### padronizados
matcor=round(cor(dados),4)
matcor


### matriz de correla??o colorida
ggcorrplot(matcor, hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           method="circle", 
           colors = c("tomato2", "white", "springgreen3"), 
           title="Correlograma", 
           ggtheme=theme_bw)

res.pca.cor <- PCA(dados, scale.unit = T, graph = T)  ### mapa do ACP
round(res.pca.cor$eig,3) ## autovalores
round(res.pca.cor$svd$V,3)  ## autovetores

# Quanto mais pr?xima uma vari?vel for do c?rculo de correla??es, melhor sua representa??o no mapa fatorial (e mais importante ? a vari?vel para a interpreta??o desses componentes)
# As vari?veis pr?ximas ao centro do gr?fico s?o menos importantes para os primeiros componentes.
# No gr?fico abaixo os componentes s?o coloridas de acordo com os valores do coseno quadrado:

fviz_pca_var(res.pca.cor, col.var="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.6) + theme_minimal()


# Coordenadas de vari?veis
round(res.pca.cor$var$coord,2)

# Cos2: ? uma medida que indica a qualidade da representa??o para vari?veis no mapa fatorial
round(res.pca.cor$var$cos2,2)

# A contribui??o das vari?veis pode ser extra?da da seguinte forma:
round(res.pca.cor$var$contrib,2)

# Quanto maior o valor da contribui??o, mais a vari?vel contribui para o componente.
# As vari?veis mais importantes associadas a um determinado PC podem ser visualizadas, usando a fun??o fviz_contrib () [factoextra package], da seguinte forma:
# Contribui??es de vari?veis no PC1
fviz_contrib(res.pca.cor, choice = "var", axes = 1)+ theme_minimal()

# Controle as cores das vari?veis usando suas contribui??es
# a cor representa a contribui??o conjunta dim1-dim2
fviz_pca_var(res.pca.cor, col.var="contrib")+ theme_minimal()

res.desc <- dimdesc(res.pca.cor, axes = c(1,2))
# Descri??o da dimens?o 1
res.desc$Dim.1

# Gr?fico de escores (indiv?duos ou pontos objetos)
# As coordenadas dos escores nos componentes principais s?o:
round(res.pca.cor$ind$coord,2)

fviz_pca_ind(res.pca.cor)+ theme_minimal()

fviz_pca_ind(res.pca.cor, col.ind="cos2") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=0.50) + theme_minimal()

# Contribui??es de escores para PC1
fviz_contrib(res.pca.cor, choice = "ind", axes = 1)+ theme_minimal()

# Contribui??es de escores para PC2
fviz_contrib(res.pca.cor, choice = "ind", axes = 2)+ theme_minimal()

# Contribui??o total em PC1 e PC2
fviz_contrib(res.pca.cor, choice = "ind", axes = 1:2)+ theme_minimal()

fviz_pca_biplot(res.pca.cor) + theme_minimal()

######### outra forma
pc=princomp(dados,cor=TRUE,score=TRUE)
summary(pc)
plot(pc)
plot(pc,type='l')
biplot(pc)
attributes(pc)
pc$loadings

#####################################################################
#####################################################################
#####################################################################
### outra forma - a que eu mais gostei
require(FactoMineR)
require(factoextra)

# viasualizar os dados
head(dados[1:3,])
## caso tenha NA
dados=na.omit(dados)
## escala padroniza??o
dados3=scale(dados)
dados3
### gerando a ACP
res.acp=PCA(dados3,graph=F)
### autovalores e autovetores
res.val=get_eigenvalue(res.acp)
res.val
## grafico
fviz_eig(res.acp,label=T,ylab=c(0,110))
### extraindo os autovalores e autovetores
autovar=get_pca_var(res.acp)
autovet=get_pca_ind(res.acp)
# novo gr?fico
fviz_pca(res.acp,col.var="blue")
## criando um grupo
grupo=as.factor(dados2[,1])
## novo grafico
fviz_pca_biplot(res.acp,habillage=grupo,
                title="GR?FICO DE ACP DOS TRATAMENTOS",col.var="red",ellipse.level=0.95)


# install.packages("corrplot")
require(corrplot)
## checar a qualidade
autovar$cos2
## grafico de correla??o
corrplot(autovar$cos2)


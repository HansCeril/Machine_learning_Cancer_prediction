#===================================================Packages==========================================================
######################################################################################################################

library(fields)
library(spam)
library(grid)
library(maps)
library(dendextend)
library(gplots)
library(glmnet)
library(pROC)
library(FactoMineR)
library(class)
library(pROC)
library(rpart)
library(randomForest)
library(vegan)

#***********************************************Chemin repertory*******************************************************
setwd("/Users/ceril/Desktop/Machine_learning_projet_finale")
#**********************************************************************************************************************

#=================================================ReadFile==========================================================
#*******************************************************************************************************************
X = read.table("xtrain.txt",row.names=1)
Y = read.table("ytrain.txt")
X = t(X)
print(X)
Y = Y$V1

#==================================================Data_Plot=========================================================
image.plot(t(X))
print (ncol(X))

#====================================================================================================================
#+++++++++++++++++++++++++++++++++++++++++++Coefficient de corrélation+++++++++++++++++++++++++++++++++++++++++++++++
# On effectue une corrélation de l'expréssion des gène et des labels.

CXY = NULL
for (k in 1:ncol(X)){
  CXY[k] = abs(cor(X[,k],as.numeric(Y)))
}
hist(CXY) # On effectue un histogramme de distribution.
#===================================================================================================================

#===================================================Kmeans=========================================================
#===================================================================================================================
# On commence par le k-means car il est facile à mettre en œuvre et à comprendre. 
#Je mets k à une valeur élevée et j'exécute k-means de 2 clusters à k clusters et ensuite je trace le graphique 
# du coude bien connu pour avoir une idée du nombre de clusters dont j'ai besoin.

# Dans ce cas je choisis k = 50.

# -------------------Fonction MULTIKMEANS----------------

multiKmeans = function(data,max.clusters,iter){
  # cluster sum of squares
  css = list()
  km = list()
  for(i in seq(2:max.clusters)){
    cl = NULL
    cl = kmeans(data,centers=i,iter.max=iter)
    css[i]=cl$tot.withinss   
    km[[i]]=cl
  } 
  df = data.frame(ss=matrix(unlist(css), nrow=max.clusters-1, byrow=T))
  cluster = seq(2,max.clusters)
  df = cbind(df,as.data.frame(cluster))
  return(list(css=df,cluster=km))  
}
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

#---------------------------ELBOWGRAPH-------------------------------------------------
elbowGraph = function(multiKmeansObj){
  require(ggplot2)
  p = NULL
  p = qplot(x=as.numeric(cluster),y=ss,data=multiKmeansObj)
  p = p + geom_point(colour = "red", size = 2) + geom_text(aes(label=as.numeric(cluster),size=1,hjust=-0.5, vjust=0))
  p = p + labs(title = "Elbow Graph") + xlab("# of Clusters") + ylab("Sum of Squares") + theme_bw()
  p = p +  theme(legend.position = "none") 
  return(p)
}
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

#-------------------------PlotMUltiKmean---------------
plot.kmeans2 =
  function (x, data = NULL, class = NULL, legend.position = c("right",
                                                              "bottom", "left", "top", "none"), title = "K-Means Results",
            xlab = "Principal Component 1", ylab = "Principal Component 2",
            ...)
  {
    require(useful)
    toPlot <- fortify(model = x, data = data)
    legend.position <- match.arg(legend.position)
    ggplot(toPlot, aes(x = .x, y = .y, colour = .Cluster, label= .Cluster)) +
      geom_point(aes_string(shape = class)) + scale_color_discrete("Cluster") +
      theme(legend.position = legend.position) + labs(title = title, x = xlab, y = ylab) + geom_text(size=3,hjust=-0.5, vjust=0)
  }
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°


cl = multiKmeans(X, 50, 200)
elbow = elbowGraph(cl$css)
elbow

# A partir du graphe obtenue on pourrait dire que K est compris entre 2 et 6 clusters.

# Le k-means initial donne 2 <= k <= 6 grappes et je pense que 4 ou 5 grappes semblent correctes mais 
# je vais utiliser le groupement hiérarchique pour mieux comprendre le nombre de grappes. 
# Veuillez noter que j'utilise d'autres algorithmes de clustering pour avoir un aperçu du nombre de grappes 
#que je choisirai pour mes k-means précédents.

w = which(CXY>.29)

ILDs = scale(X)
dis = vegdist(X[,w],"euclidean")
hc = hclust(dis, "ward.D2")
plot(hc)

#Dans le graphique ci-dessus, ce que nous cherchons est pour les grandes tiges et il semble que 3,6 ou 7 grappes
#sont de grands candidats.

#Voyons comment ces différentes options de clusters ressemblent.

plot(hc, labels = FALSE) 
rect.hclust(hc, k = 3, border = "red") # 3 cluster

plot(hc, labels = FALSE)
rect.hclust(hc, k = 6, border = "blue") # 6 cluster

plot(hc, labels = FALSE)
rect.hclust(hc, k = 12, border = "yellow") # 12 cluster

# 3 clusters seem like a potential solution. Let’s see how looks like 3 clusters
clus3 = cl$cluster[[3]]
print (clus3)
print(X)
plot(clus3, X)
print (clus3$centers)

clus12 = cl$cluster[[12]]
plot(clus12, X)

clus8 = cl$cluster[[6]]
plot(clus8, X)
a = (clus8$centers)

# il semble que 3 et 8 sont de bonnes options pour le nombre de clusters
#cependant avec 8 clusters, il y a beaucoup de chevauchement entre les clusters.

#============================================Lasso===========================================================
library(glmnet)
n.train = 170
set.seed(Seed)
test <- pred <-NULL
B=100

for (b in 1:B){
  i.TRAIN = sample(1:n,n.train) ; 
  i.TEST = setdiff(1:n,i.TRAIN) ; 
  test = c(test,i.TEST)
  rho = NULL
  for (k in 1:ncol(X)){
    rho[k] = abs(cor(X[i.TRAIN,k],as.numeric(Y[i.TRAIN])))
  }
  w = which(abs(rho)>.22)
  TRAIN = data.frame(X=scale(data[i.TRAIN,w]),Y=data$Y[i.TRAIN])
  Xtest = scale(data[i.TEST,w],center=apply(data[i.TRAIN,w],2,mean),scale=apply(data[i.TRAIN,w],2,sd))
  TEST = data.frame(X=Xtest,Y=data$Y[i.TEST])
  cv = cv.glmnet(as.matrix(TRAIN[,-dim(TRAIN)[2]]),TRAIN$Y,family="binomial",type.measure="auc",nfolds=3)
  mod.lr = glmnet(as.matrix(TRAIN[,-dim(TRAIN)[2]]),TRAIN$Y,lambda=cv$lambda.1se,family="binomial")
  pred = c(pred,predict(mod.lr,as.matrix(TEST[,-dim(TRAIN)[2]]),type="response"))
}

colnames(TRAIN)[mod.lr$beta@i]

boxplot(TRAIN[,mod.lr$beta@i[1]]~data$Y[i.TRAIN])
boxplot(TEST[,mod.lr$beta@i[1]]~data$Y[i.TEST])

boxplot(TRAIN[,mod.lr$beta@i[3]]~data$Y[i.TRAIN])
boxplot(TEST[,mod.lr$beta@i[3]]~data$Y[i.TEST])

roc.lr = roc(data[test,dim(data)[2]],pred) #Area under the curve: 0.904
plot(roc.lr)
grid()
print(roc.lr)

#================================================KNN===========================================================
library(class)
library(pROC)

Seed = "1223"
set.seed(Seed)

n = nrow(X) ; n.train = round(n*2/3)
data = data.frame(cbind(X,Y))
data$Y = as.factor(data$Y)

set.seed(Seed)
test <- prob1 <- prob3 <- prob5 <- prob7<- prob9 <-NULL
B = 100

for (b in 1:B){
  i.TRAIN = sample(1:n,n.train) ; 
  i.TEST = setdiff(1:n,i.TRAIN) ; 
  test = c(test,i.TEST)
  rho = NULL
  for (k in 1:ncol(X)){
    rho[k] = abs(cor(X[i.TRAIN,k],as.numeric(Y[i.TRAIN])))
  }
  w = which(abs(rho)>.25)
  TRAIN = data[i.TRAIN,c(w,dim(data)[2])]
  TEST = data[i.TEST,c(w,dim(data)[2])]
  knn1 = knn(TRAIN[,-(length(w)+1)],TEST[,-(length(w)+1)],TRAIN$Y,k=1,prob=TRUE)
  prob1 = c(prob1,attributes(knn1)$prob)
  knn3 = knn(TRAIN[,-(length(w)+1)],TEST[,-(length(w)+1)],TRAIN$Y,k=3,prob=TRUE)
  prob3 = c(prob3,attributes(knn3)$prob)
  knn5 = knn(TRAIN[,-(length(w)+1)],TEST[,-(length(w)+1)],TRAIN$Y,k=5,prob=TRUE)
  prob5 = c(prob5,attributes(knn5)$prob)
  knn7 = knn(TRAIN[,-(length(w)+1)],TEST[,-(length(w)+1)],TRAIN$Y,k=7,prob=TRUE)
  prob7 = c(prob7,attributes(knn7)$prob)
  knn9 = knn(TRAIN[,-(length(w)+1)],TEST[,-(length(w)+1)],TRAIN$Y,k=9,prob=TRUE)
  prob9 = c(prob9,attributes(knn7)$prob)
}

Ytest = data[test,dim(data)[2]]

boxplot(prob7~Ytest,ylab="Predicted probabilities", xlab = "label") # Prob7 semble le mieux adapté
boxplot(prob5~Ytest,ylab="Predicted probabilities")
boxplot(prob3~Ytest,ylab="Predicted probabilities")
boxplot(prob1~Ytest,ylab="Predicted probabilities")

#***********************ROC_Curve************************

roc5 = roc(data[test,dim(data)[2]],prob5)
plot(roc5)
roc9 = roc(data[test,dim(data)[2]],prob9)
plot(roc9,add=TRUE,col="gray")
roc7 = roc(data[test,dim(data)[2]],prob7)
plot(roc7,add=TRUE,col="red",lty=2)
legend("bottomright",legend=c("k=5","k=7","k=9"),lty=c(1,2,1),col=c("black","blue","red"))

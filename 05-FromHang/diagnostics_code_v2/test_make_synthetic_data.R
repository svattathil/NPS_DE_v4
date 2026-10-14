set.seed(7); n<-73; G<-300; K<-7
ids<-sprintf("S%03d",1:n)
z<-rnorm(n)
Y<-sapply(1:K,function(k) rbinom(n,1,plogis(-1.3+0.25*k+0.8*z)))
colnames(Y)<-c("agit","anx","apa","depd","irr","psych","nite")
ph<-data.frame(protsample=ids,Y,age=1)
write.table(ph,"2_Pipeline/02-Prepare_analysis_data/phenos_cleaned_SVs_TEST.txt",sep="\t",quote=FALSE,row.names=FALSE)
X<-matrix(rnorm(G*n),G,n); X[sample(length(X),0.03*length(X))]<-NA
X[1:5,]<- X[1:5,] + 0.9*matrix(rep(z,5),5,byrow=TRUE)  # 5 true signals
colnames(X)<-ids
res<-data.frame(protein=sprintf("P%04d",1:G),X,check.names=FALSE)
write.table(res,"2_Pipeline/02-Prepare_analysis_data/resid_regress_ignoreSVs_TEST.txt",sep="\t",quote=FALSE,row.names=FALSE)

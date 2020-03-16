######## CV and fitting for GDSC data, fully penalized model
# code for creating CV loops taken from
# https://github.com/zhizuio/IPFStructPenalty

library(glmnet)
library(rrpack)

load("GDSC_my1.RData")
source("admm.R")

n <- dim(GDSC$y)[1]
q <- dim(GDSC$y)[2]
p <-  GDSC$p
num.nonpen <- GDSC$num.nonpen
MultSeed <- c(2067, 6253, 8009, 10893, 32940, 43024, 59327, 72, 82702, 148392)

# Consider keeping for loop commented out, and run in parallel
#for(s in 1:10){
  s = 1
  cat("Seed ", s, ": ", MultSeed[s], "\n")
  set.seed(MultSeed[s])
  # select train dataset of 80% cell lines of each tissue
  trainID <- NULL
  for(i.tissue in 1:num.nonpen) trainID <- c(trainID, sample(which(GDSC$x[,i.tissue]==1),round(sum(GDSC$x[,i.tissue]==1)*0.8)))
  # make sure at least one mutated cell lines of each cancer tissue for every mutation features
  repeat{
    if(min(colSums(GDSC$x[trainID, num.nonpen+(1+p[1]+p[2]):sum(p)]))>=1) break
    trainID <- NULL
    for(i.tissue in 1:num.nonpen) trainID <- c(trainID, sample(which(GDSC$x[,i.tissue]==1),round(sum(GDSC$x[,i.tissue]==1)*0.8)))
  }
  x_test <- GDSC$x[-trainID,]
  x <- GDSC$x[trainID,]
  y_test <- GDSC$y[-trainID,]
  y <- GDSC$y[trainID,]
  
  x_test[,num.nonpen+1:p[1]] <- log(GDSC$x[-trainID, num.nonpen+1:p[1]])
  x[,num.nonpen+1:p[1]] <- log(GDSC$x[trainID, num.nonpen+1:p[1]])
  
  foldid <- sample(rep(seq(5),length=dim(x)[1]))
  nfolds = max(foldid)
  
  ################################################################################################
  # iRRR
  ################################################################################################
  
  lambda.seq = 10^seq(-2,0,by=0.2)
  mspe.seq = matrix(NA,ncol=length(lambda.seq),nrow=nfolds)
  for (k in 1:nfolds){ 
    indx = which(foldid!=k)
    ind = which(foldid==k)
    # Setting up the data for training
    X0 = x[indx,1:num.nonpen]
    X1 = scale(x[indx,(num.nonpen+1):(num.nonpen+p[1])])
    X2 = scale(x[indx,(num.nonpen+p[1]+1):(num.nonpen+p[1]+p[2])])
    X3 = x[indx,(num.nonpen+p[1]+p[2]+1):(num.nonpen+sum(p))]
    Xlist = list(X0,X1,X2,X3)
    ww = rep(NA,length(Xlist))
    for (i in 1:length(Xlist)){
      ww[i] = svd(Xlist[[i]])$d[1]*(sqrt(q) + sqrt(as.numeric(rankMatrix(Xlist[[i]]))))/dim(Xlist[[i]])[1]
    }
    # Setting up the data for testing
    X0 = x[ind,1:num.nonpen]
    X1 = scale(x[ind,(num.nonpen+1):(num.nonpen+p[1])])
    X2 = scale(x[ind,(num.nonpen+p[1]+1):(num.nonpen+p[1]+p[2])])
    X3 = x[ind,(num.nonpen+p[1]+p[2]+1):(num.nonpen+sum(p))]
    X_test = cbind(X0,X1,X2,X3)
    for (l in 1:length(lambda.seq)){
      fit <- iRRR(y[indx,],X=Xlist,lam1 = lambda.seq[l],paramstruct = list(maxrho = 500,varyrho=1,noreg = c(1),
                                                                           rho = 0.0001,
                                                                           Niter = 2000,
                                                                           weight = ww,
                                                                           Tol = 1e-6))
      mu = matrix(rep(fit$mu,dim(x[ind,])[1]),nrow=dim(x[ind,])[1],byrow=T)

      yhat = X_test%*%fit$C+mu
      corr = cor(as.vector(yhat),as.vector(y[ind,]))
      MSPE = mean((y[ind,]-yhat)^2)
      mspe.seq[k,l] = MSPE
      print("local low rank")
      print(paste("lam1:",lambda.seq[l],"fold:",k))
      print(paste("MSPE:",round(MSPE,4),"corr:",round(corr,4),"rank X0:",rankMatrix(fit$A[[1]]),
      "rank X1:",rankMatrix(fit$A[[2]]),"rank X2:",rankMatrix(fit$A[[3]]),
      "rank X3:",rankMatrix(fit$A[[4]])))
      
    }
  }
  print("local low rank")
  print(mspe.seq)
  print(lambda.seq)
  meanMSPE = apply(mspe.seq,2,mean)
  lam1.cvOpt = lambda.seq[which.min(meanMSPE)]
  
  
  
  # Finally fit iRRR with optimal lambda value to the full training set, predict on test set
  # Setting up data
  X0 = x[,1:num.nonpen]
  X1 = scale(x[,(num.nonpen+1):(num.nonpen+p[1])])
  X2 = scale(x[,(num.nonpen+p[1]+1):(num.nonpen+p[1]+p[2])])
  X3 = x[,(num.nonpen+p[1]+p[2]+1):(num.nonpen+sum(p))]
  Xlist = list(X0,X1,X2,X3)
  ww = rep(NA,length(Xlist))
  for (i in 1:length(Xlist)){
    ww[i] = svd(Xlist[[i]])$d[1]*(sqrt(q) + sqrt(as.numeric(rankMatrix(Xlist[[i]]))))/dim(Xlist[[i]])[1]
  }
  X0_test = x_test[,1:num.nonpen]
  X1_test = scale(x_test[,(num.nonpen+1):(num.nonpen+p[1])])
  X2_test = scale(x_test[,(num.nonpen+p[1]+1):(num.nonpen+p[1]+p[2])])
  X3_test = x_test[,(num.nonpen+p[1]+p[2]+1):(num.nonpen+sum(p))]
  X_test = cbind(X0_test,X1_test,X2_test,X3_test)






  # Local Low Rank
  fit <- iRRR(y,X=Xlist,lam1 = lam1.cvOpt,paramstruct = list(maxrho = 500,varyrho=1,noreg = c(1),
                                                                            rho = 0.0001,
                                                                            Niter = 2000,
                                                                            weight = ww,
                                                                            Tol = 1e-6))
  mu = matrix(rep(fit$mu,dim(x_test)[1]),nrow=dim(x_test)[1],byrow=T)


  yhat = X_test%*%fit$C+mu
  corr = cor(as.vector(yhat),as.vector(y_test))
  MSPE = mean((y_test-yhat)^2)
  print("local low rank")
  print(paste("lam1:",lam1.cvOpt,"fold: test sample"))
  print(paste("MSPE:",round(MSPE,4),"corr:",round(corr,4),"rank X0:",rankMatrix(fit$A[[1]]),
              "rank X1:",rankMatrix(fit$A[[2]]),
              "rank X2:",rankMatrix(fit$A[[3]]),
              "rank X3:",rankMatrix(fit$A[[4]])))
  
  
  iRRRfinalfit = list(fit=fit,yhat=yhat,y=y_test,lam1=lam1.cvOpt,corr=corr,MSPE=MSPE,
                      ranks = c(rankMatrix(fit$A[[1]]),rankMatrix(fit$A[[2]]),
                                rankMatrix(fit$A[[3]]),rankMatrix(fit$A[[4]])))
  
  
  ################################################################################################
  # GLMNET
  ################################################################################################
  alpha.seq = seq(0,1,length.out=11)
  params = matrix(NA,ncol=3,nrow=length(alpha.seq))
  colnames(params) = c("cvm","lambda","alpha")
  cvglmnetfits = list()
  for (i in 1:length(alpha.seq)){
    cvGLMNET = cv.glmnet(x,y,alpha=alpha.seq[i],nfolds=5,foldid=foldid,family="mgaussian",penalty.factor=c(rep(0,num.nonpen),rep(1,sum(p))))
    cvglmnetfits[[i]] = list(fit=cvGLMNET,yhat=yhat,y=y_test,corr=corr,MSPE=MSPE)
    params[i,] = c(cvGLMNET$cvm[which.min(cvGLMNET$cvm)],cvGLMNET$lambda[which.min(cvGLMNET$cvm)],alpha.seq[i])
  }
  # Fit using best values of lambda and alpha
  optInd = which.min(params[i,1])
  glmnetoptfit = cvglmnetfits[[optInd]]$fit
  
  # Predict on new dataset
  yhat = predict(glmnetoptfit,newx=x_test)
  corr = cor(as.vector(yhat),as.vector(y_test))
  MSPE = mean((as.vector(y_test)-as.vector(yhat))^2)
  glmnetfinalfit = list(fit=glmnetoptfit,yhat=yhat,y=y_test,lambda=glmnetoptfit$lambda.min,
                        alpha=alpha.seq[optInd],corr=corr,MSPE=MSPE)

  ################################################################################################
  # RRR
  ################################################################################################
  
  cvfit <- cv.rrr(Y=y,X=x,nfold=5,norder=foldid)
  yhat = x_test%*%coef(cvfit)
  corr = cor(as.vector(yhat),as.vector(y_test))
  MSPE = mean((as.vector(y_test)-as.vector(yhat))^2)
  rrrfinalfit = list(fit=cvfit,yhat=yhat,y=y_test,rank=cvfit$rank,corr=corr,MSPE=MSPE)
  
#}
  



  save.image(file=paste0("results_nonpen_s=",MultSeed[s],".Rdata"))  






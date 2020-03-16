require(Matrix)
require(rrpack)
require(glmnet)

library(doParallel)
library(foreach)

source('admm.R')


### Wrap all of this in a function
numericalSimulation = function(setting){
  
  # setting is a 3-vector containing rhoY (rho_epsilon in paper), rhoX and s (the setting)
  
  # Parameters of data simulation
  p0 = 0     # dimension of clinical data
  p1 = p2 = 100
  p = p1 + p2 + p0
  r1 = 4
  r2 = 6
  q = 24      # n. of drugs 
  n = 90     # samples
  
  s = setting[3]
  
  ### correlation among the drugs, i.e. the noises are correlation w.r.t q
  rhoY = setting[1]
  Omega = matrix(rhoY, ncol = q, nrow = q); diag(Omega) = 1
  Sigma = qr.solve(Omega)        # generate covariance matrix
  out = eigen(Sigma, symmetric = TRUE)
  S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
  # generate noise
  noise = matrix(rnorm(q*n), nrow = n, ncol = q) %*% S.sqrt
  
  ### correlation with each X_k
  rhoX = setting[2]
  Omega1 = matrix(rhoX, ncol = p1, nrow = p1);Omega2 = matrix(rhoX, ncol = p2, nrow = p2)
  diag(Omega1) = 1; diag(Omega2) = 1
  # generate covariance matrix
  Sigma1 = qr.solve(Omega1); Sigma2 = qr.solve(Omega2)
  out1 = eigen(Sigma1, symmetric = TRUE) ; out2 = eigen(Sigma2, symmetric = TRUE)
  S.sqrt1 = out1$vectors %*% diag(out1$values^0.5) %*% t(out1$vectors)
  S.sqrt2 = out2$vectors %*% diag(out2$values^0.5) %*% t(out2$vectors)
  
  ### simulate ss times
  msee =  mspe = mse1 = mse2 = mse3 = rank1 = rank2 = rank0 =
    mspe.rrr = msee.rrr = enet.mse = enet.msp  = c()
  n.iters = 50
  for (ss in 1:n.iters) {
    print(ss)
    if (rhoX==0 & rhoY==0){
      Xlist = list( # need to be centered
        X1 = scale(matrix(rnorm(n*p1), nr = n, nc = p1) , scale=F) ,
        X2 = scale(matrix(rnorm(n*p2), nr = n, nc = p2) , scale=F) )
      XlistTest = list(
        X1 = scale(matrix(rnorm(n*p1), nr = n, nc = p1) , scale=F) ,
	X2 = scale(matrix(rnorm(n*p2), nr = n, nc = p2) , scale=F) )

    } else{
      # X is correlated
      Xlist = list( # need to be centered
        # X0 = matrix( rbinom(n*p0,1,0.5), nr = n, nc = p0),
        X1 = scale(matrix(rnorm(n*p1), nr = n, nc = p1) %*% S.sqrt1 , scale=F) ,
        X2 = scale(matrix(rnorm(n*p2), nr = n, nc = p2) %*% S.sqrt2 , scale=F) )
      XlistTest = list( # need to be centered
        # X0 = matrix( rbinom(n*p0,1,0.5), nr = n, nc = p0),
        X1 = scale(matrix(rnorm(n*p1), nr = n, nc = p1) %*% S.sqrt1 , scale=F) ,
        X2 = scale(matrix(rnorm(n*p2), nr = n, nc = p2) %*% S.sqrt2 , scale=F) )

    }
    #
    # Setting S1: local low-rank
    if (s==1){
      b = list(
        beta0 = matrix(rnorm(q*p0), nr = q, nc = p0),
        beta1 = matrix(rnorm(q*r1), nr = q, nc = r1) %*% matrix(rnorm(p1*r1), nr = r1,nc = p1),
        beta2 = matrix(rnorm(q*r2), nr = q, nc = r2) %*% matrix(rnorm(p2*r2), nr = r2,nc = p2)
      )
      beta0 = do.call(cbind, b)
    }
    
    # 
    # Setting S2: low-rank + sparsity
    if (s==2){
      b = list(
        beta0 = matrix(rnorm(q*p0), nr = q, nc = p0),
        beta1 = matrix(rnorm(q*r1), nr = q, nc = r1) %*% matrix(rnorm(p1*r1), nr = r1,nc = p1),
        beta2 = matrix(rbinom(q*p2, 1, 0.5) , nr = q ,nc = p2)*rnorm(q*p2,0)
      )
      beta0 = do.call(cbind, b)
    }
    
    
    # Setting S3: GLOBAL LOW-RANK
    if (s==3){
      r0 = 6
      beta0 = t(matrix(rnorm(p*r0), nr = p, nc = r0) %*% matrix(rnorm(q*r0), nr = r0,nc = q))
    }
    
    # Setting S4: sparsity
    if (s==4){
      beta0 = matrix(rbinom(q*p, 1, 0.2) , nr = q, nc = p)* matrix( rnorm(p*q) ,nr = q,nc = p)
    }
    
    # Generating data
    X = do.call(cbind, Xlist)
    X_star = do.call(cbind,XlistTest)
    noise = matrix(rnorm(q*n), nrow = n, ncol = q) %*% S.sqrt
    noise_star = matrix(rnorm(q*n), nrow = n, ncol = q) %*% S.sqrt
    Y = X%*%t(beta0) +  noise
    Y_star = X_star%*%t(beta0) + noise_star
    
    # We need to do CV
    nfolds = 5
    foldid = sample(rep(seq(5),length=n))
    # lambda.seq = seq(0.01,.5,by=.05)
    lambda.seq = 10^seq(-3,1,length.out=11)
    mspe_irrr_cv = matrix(NA,ncol=nfolds,nrow=length(lambda.seq))
    
    for (fold in 1:nfolds){
      indx = which(foldid!=fold)
      ind = which(foldid==fold)
      X_train = list(scale(Xlist$X1[indx,],scale=F),scale(Xlist$X2[indx,],scale=F))
      Y_train = Y[indx,]
      X_test = list(scale(Xlist$X1[ind,],scale=F),scale(Xlist$X2[ind,],scale=F))
      Y_test = Y[ind,]
      # Calculate weights
      ww = rep(NA,length(X_train))
      for (i in 1:length(X_train)){
        ww[i] = svd(X_train[[i]])$d[1]*(sqrt(q) + sqrt(as.numeric(rankMatrix(X_train[[i]]))) )/length(indx)
      }
      
      for (l in 1:length(lambda.seq)){
        # Reshape X as list
        fit = iRRR(Y_train, X_train, lam1= lambda.seq[l],
                   paramstruct = list(maxrho = 50,noreg = c(0),varyrho=1,
                                      rho = .0001,
                                      Niter = 2000,
                                      weight = ww,
                                      Tol = 1e-6))
        # Predict on X_test
        mu = matrix(rep(fit$mu,dim(Y_test)[1]),nrow=dim(Y_test)[1],byrow=T)
        ha = do.call(rbind,fit$A)
        A = fit$A
        
        mspe_irrr_cv[l,fold] = mean((Y_test - (do.call(cbind,X_test)%*%ha))^2)
      }
    }
    lam1opt = lambda.seq[which.min(apply(mspe_irrr_cv,1,mean))]
    plot(log10(lambda.seq),apply(mspe_irrr_cv,1,mean),type="l")
    points(log10(lam1opt),min(apply(mspe_irrr_cv,1,mean)),col="red")
    
    ww = rep(NA,length(Xlist))
    for (i in 1:length(Xlist)){
      ww[i] = svd(Xlist[[i]])$d[1]*(sqrt(q) + sqrt(as.numeric(rankMatrix(Xlist[[i]]))))/n
    }
    
    fit = iRRR(Y, Xlist, lam1= lam1opt,
               paramstruct = list(maxrho = 50,noreg = c(0),varyrho=1,
                                  rho = .0001,
                                  Niter = 2000,
                                  weight = ww,
                                  Tol = 1e-6))
    # Predict
    mu = matrix(rep(fit$mu,dim(Y)[1]),nrow=dim(Y)[1],byrow=T)
    ha = do.call(rbind,fit$A)
    A = fit$A
    
    msee[ss] = mean((beta0 - t(ha))^2)
    mspe[ss] = mean((Y_star - (X_star%*%ha+mu))^2) 
    
    # global low-rank
    fit.rr = cv.rrr(Y, X, norder=foldid,nfold=5)
    be.rrr = fit.rr$coef
    mspe.rrr[ss] = mean((Y_star - X_star%*%be.rrr)^2)
    msee.rrr[ss] = mean((t(beta0) - be.rrr)^2)
    
    # glmnet
    alpha.seq = seq(0,1,length.out=6)
    mspealpha = rep(NA,length(alpha.seq))
    cv_enet <- list()
    cv_enet_cvm = c()
    for (zz in 1:length(alpha.seq)){
      cv_enet[[zz]] <- cv.glmnet(X, Y, 
                                 intercept = F, 
                                 alpha = alpha.seq[zz], 
                                 standardize.response = F,
                                 foldid = foldid,
                                 family="mgaussian")
      cv_enet_cvm[zz] <- min(cv_enet[[zz]]$cvm)
    }
    fit.enet <- cv_enet[[which.min(cv_enet_cvm)]]
    enet.pred =  predict(fit.enet, newx = X_star,s="lambda.min")
    enet.msp[ss] = mean((Y_star-enet.pred[,,1])^2)
    enet.b = t(do.call(cbind,coef(fit.enet,s="lambda.min"))[-1,])
    enet.mse[ss] = mean((beta0 - enet.b)^2)
  }
  
  results = list(s=s,rhoY=setting[1],rhoX=setting[2],mspe=mspe,mspe.rrr=mspe.rrr,
                 enet.msp=enet.msp,msee=msee,msee.rrr=msee.rrr,enet.mse=enet.mse)
}

# For parallel processing
#cl<-makeCluster(24)
#registerDoParallel(cl)

# Set up list of simulation settings
sim.settings = cbind(c(0,0.3,0.6,0,0,0),c(0,0,0,0.3,0.6,0.9))
sim.settings = cbind(rbind(sim.settings,sim.settings,sim.settings,sim.settings),
                     c(rep(1,nrow(sim.settings)),
                       rep(2,nrow(sim.settings)),
                       rep(3,nrow(sim.settings)),
                       rep(4,nrow(sim.settings))))
colnames(sim.settings) = c("rhoY","rhoX","s")


# For parallel processing
#res = foreach(i=1:nrow(sim.settings),
#              .combine=list,
#              .multicombine=T,
#              .packages = c("Matrix","rrpack","glmnet","MASS")) %dopar% numericalSimulation(sim.se#ttings[i,])
#
#save.image(file="numSim2.Rdata")
#stopImplicitCluster()

# For running a single setting
xxx = 1
res = numericalSimulation(sim.settings[xxx,])

save(res,file=paste0("results_",xxx,".Rdata"))



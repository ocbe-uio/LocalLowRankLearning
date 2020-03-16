##############################################################################
############ Integrative Multi-View Regression:  #############################
## A function to implement the ADMM algorithm for the iRRR method.
##############################################################################
## Code translated from the MATLAB code provided by 
## https://github.com/reagan0323/iRRR



# Libraries
library(MASS)
library(parallel)
library(Matrix)

iRRR = function(Y,X,lam1,paramstruct=NULL){
  # Inputs:
  # Y: (n x q) continous response matrix
  # X: list of length K, where each element in the list
  #    contains a (n*p_i) predictor data matrix
  # lam1: >0, tuning for nuclear norm
  # paramstruct: a list with  lam0: tuning for ridge penalty
  #                           weights: (Kx1) weight vector
  #                           noreg: a vector of view numbers that are not to be regularized
  
  # Initial setup
  K = length(X)
  weight = rep(1,K)
  Tol = 1e-3
  Niter = 500
  varyrho = 1
  rho = 0.1
  lam0 = 0
  maxrho = 5
  randomstart = 0
  fig = 1
  if (!is.null(paramstruct)){
    if ("lam0"%in%names(paramstruct)){
      lam0 = paramstruct$lam0
    }
    if("weight"%in%names(paramstruct)){
      weight = list(paramstruct)[[1]]$weight
      #print(weight)
    }
    if("Tol"%in%names(paramstruct)){
      Tol = paramstruct$Tol
    }
    if ("Niter"%in%names(paramstruct)){
      Niter = paramstruct$Niter
    }
    if("randomstart"%in%names(paramstruct)){
      randomstart=paramstruct$randomstart
    }
    if("varyrho"%in%names(paramstruct)){
      varyrho = paramstruct$varyrho
    }
    if ("rho"%in%names(paramstruct)){
      rho = paramstruct$rho
    }
    if("fig"%in%names(paramstruct)){
      fig=paramstruct$fig
    }
    if("noreg"%in%names(paramstruct)){
      noreg=paramstruct$noreg
    }
    if("maxrho"%in%names(paramstruct)){
      maxrho = paramstruct$maxrho
    }
  }
  
  # Initialization
  n = nrow(Y)
  q = ncol(Y)
  ## Centering Y
  #meanY = apply(Y,2,mean)
  #Y = scale(Y,scale=F)
  p = rep(0,K)
  cX = c()
  meanX = c()
  for (i in 1:K){
    ni = dim(X[[i]])[1]
    p[i] = dim(X[[i]])[2]
    if (ni!=n){
      print("Error, samples do not match")
    }
    # Column center Xs
    meanX = c(meanX,apply(X[[i]],2,mean,na.rm=T))
    X[[i]] = scale(X[[i]],scale=F)
    # Normalize (wrt to weights)
    X[[i]] = X[[i]]/weight[i]
    cX = cbind(cX,X[[i]])
  }
  
  # Initial parameter estimates
  mu = apply(Y,2,mean,na.rm=T)
  # Majorize Y to get a working Y
  wY = Y
  temp = rep(1,n)%*%t(mu)
  wY[which(is.na(wY),arr.ind=T)] = temp[which(is.na(wY),arr.ind=T)] #wY is a complete matrix
  mu = apply(wY,2,mean)
  wY1 = scale(wY,scale=F)
  
  B = list() # Lagrange params for B
  Theta = list() # vertically concatenated B
  cB = c()
  for (i in 1:K){
    if (randomstart){
      B[[i]] = matrix(rnorm(p[i]*q),ncol=q)
    } else{
      B[[i]] = ginv(t(X[[i]])%*%X[[i]]) %*% t(X[[i]])%*%wY1 # OLS with pseudoinverse
    }
    Theta[[i]] = matrix(0,nrow=p[i],ncol=q)
    cB = rbind(cB,B[[i]])
  }
  A = B #low-rank alias
  cA = cB
  cTheta = matrix(0,nrow=sum(p),ncol=q)
  
  tmp = svd(1/sqrt(n)*cX)
  D_cX = tmp$d
  V_cX = tmp$v
  
  if (!varyrho){ # i.e. fixed rho
    DeltaMat = V_cX%*% diag(1/(D_cX^2 + lam0 + rho)) %*%t(V_cX)
          +  (diag(sum(p))-V_cX%*%t(V_cX))/(lam0+rho) # inv(1/n*X^TX+(lam0+rho)I)
  }
  
  # Check objective value
  obj = ObjValue1(Y,X,mu,A,lam0,lam1)
  obj_ls = ObjValue1(Y,X,mu,A,0,0)
  
  ####################################################################################
  ###### ADMM   ######################################################################
  ####################################################################################
  niter = 0
  diff = Inf
  rec_obj = c(obj,obj_ls)
  rec_Theta = c()
  rec_primal = c()
  rec_dual = c()
  while (niter<Niter && abs(diff)>Tol){
    niter = niter+1
    cB_old = cB
    
    ###############     Majorization      #########################################
    Eta = rep(1,n)%*%t(mu) + cX%*%cB
    wY = Y
    wY[which(is.na(wY),arr.ind=T)] = Eta[which(is.na(wY),arr.ind=T)] #working response
    mu=(apply(wY,2,mean))
    wY1 = scale(wY,scale=F) # column centered
    
    # estimate concatenated B
    if (varyrho){
      DeltaMat = V_cX%*% diag(1/(D_cX^2 + lam0 + rho)) %*%t(V_cX)+ ## NB MAYBE TRANSPOSED?
        (diag(sum(p))-V_cX%*%t(V_cX))/(lam0+rho) # inv(1/n*X^TX+(lam0+rho)I)
    }
    cB = DeltaMat%*%((1/n)*t(cX)%*%wY1+rho*cA+cTheta)
    # 
    begin = 1
    for (i in 1:K){
      end = begin + p[i]-1
      #print(c(begin,end))
      B[[i]] = cB[begin:end,]
      begin = end + 1
    }
    
    # # estimate A_i in parallel
    # # update Theta_i in parallel right after estimating A_i
    # # NB! This means a cluster will have to be set up, let's call it cl
    # # This solution is not the most elegant one, but it works
    # clusterExport(cl=cl, varlist=c("A","B","Theta","rho","SoftThres","lam1"),envir = environment())
    # #print(environment())
    # listreturn = parLapply(cl,1:K,parUpdate)
    # begin = 1
    # for (i in 1:K){
    #   end = begin+p[i]-1
    #   A[[i]] = listreturn[[i]]$Anew
    #   cA[begin:end,] = A[[i]] # Reshape
    #   Theta[[i]] = listreturn[[i]]$Thetanew
    #   cTheta[begin:end,] = Theta[[i]] # Reshape
    #   begin = end+1
    # }
    
    # Estimate A_i without parallelization
    for (i in 1:K){
      temp = B[[i]]-Theta[[i]]/rho
      tmpSVD = svd(temp)
      # Check if view number is in noreg, then there is no soft-tresholding
      # and hence no low-rankness
      if (i%in%noreg){
        A[[i]] = temp
      } else{
        A[[i]] = tmpSVD$u%*%SoftThres(tmpSVD$d,(lam1/rho))%*%t(tmpSVD$v)
      }
      Theta[[i]] = Theta[[i]]+rho*(A[[i]]-B[[i]])
    }
    # Update cA and cB
    begin = 1
    for (i in 1:K){
      end = begin+p[i]-1
      cA[begin:end,] = A[[i]] # Reshape
      cTheta[begin:end,] = Theta[[i]] # Reshape
      begin = end+1
    }
    
    # Update rho
    if (varyrho){
      rho = min(maxrho,1.1*rho) # increasing rho
    }
    
    # Stopping rule ###########
    
    # Primal and dual residuals
    primal = norm((cA-cB),type="F")
    rec_primal = c(rec_primal,primal)
    dual = norm(cB-cB_old,type="F")
    rec_dual = c(rec_dual,dual)
    
    # Objective function value
    obj = ObjValue1(Y,X,mu,A,lam0,lam1)
    obj_ls = ObjValue1(Y,X,mu,A,0,0)
    rec_obj = rbind(rec_obj,c(obj,obj_ls))
    
    # Stopping rule
    #diff = mean((cA-cB)^2)
    diff = primal
    
    #print(paste("rho:",round(rho,4),"diff:",round(diff,4)))
  }
  if (niter==Niter){
    #print(paste("iRRR does NOT converge after",Niter,"iterations!"))
  } else{
    #print(paste("iRRR converges after",niter,"iterations!"))
  }
  # Outputs
  # rescale parameter estimate, and add back mean
  C = c()
  for (i in 1:K){
    A[[i]] = A[[i]]/weight[i] # NB division or multiplication?
    B[[i]] = B[[i]]/weight[i]
    C = rbind(C,A[[i]])
  }
  mu = t(t(mu)-meanX%*%C)
  return(list(C=C,mu=mu,A=A,B=B,Theta=Theta,primal=rec_primal,dual=rec_dual,rho=rho))
}

ObjValue1 = function(Y,X,mu,B,lam0,lam1){
  # Calc 1/(2n) ||Y-sum(X_i*B_i)||^2+lam0/2*sum(|B_i|^2_F)+lam1*sum(|B_i|_*)
  # with column centered Xi's and (potentially non-centered and missing) Y
  n = dim(Y)[1]
  q = dim(Y)[2]
  K = length(X)
  obj = 0
  pred = rep(1,n)%*%t(mu)
  for (i in 1:K){
    pred = pred + X[[i]]%*%B[[i]]
    obj = obj + 0.5*lam0*norm(B[[i]],type="F")^2 + lam1*sum(svd(B[[i]])$d)
  }
  obj = obj+1/(2*n)*sum(sum((Y-pred)^2,na.rm=T),na.rm=T)
  return(obj)
}

parUpdate = function(i){
  Anew = list()
  Thetanew = list()
  temp = B[[i]]-Theta[[i]]/rho
  tmpSVD = svd(temp)
  #Anew = tmpSVD$u%*%SoftThres(tmpSVD$d,(lam1/rho))%*%t(tmpSVD$v)
  Anew = tmpSVD$u%*%SoftThres(tmpSVD$d,(lam1/rho))%*%tmpSVD$v
  Thetanew = Theta[[i]]+rho*(Anew-B[[i]])
  list("Anew"=Anew,"Thetanew"=Thetanew,i=i)
}

SoftThres = function(Din,lam){
  d = diag(Din)
  d[which(d>0,arr.ind=T)] = pmax(d[d>0]-lam,0)
  d[which(d<0,arr.ind=T)] = pmin(d[d<0]+lam,0)
  return(d)
}


iRRR = compiler::cmpfun(iRRR)





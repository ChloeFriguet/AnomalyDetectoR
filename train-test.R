library(kernlab)


######################
### training phase ###
######################

# @input
# - nominalData : nominal dataset
# - sig : parameter for the RBF kernel, not as in <kernlab:rbfdot> (inverse kernel width parameter) => sigma = 1/(2*sig^2) )
# Param for K^-1 computation:
# - reg : method for regularisation / default :"none" - other : "tikhonov"
# - param : parameter for Tikhonov regularisation method or the % of the total eigenvalues kept in truncated SVD / default : 0
# - inv : method for matrix inverse computation / default : "solve" - other: "svd"; "svdTrunc"; "moorePenrose"; "cholesky"
train = function(nominalData, sig, param = 0, inv = "solve", reg="none"){
  # settings
  n=dim(nominalData)[1] # sample size
  rbfkernel = rbfdot(sigma =  1/(2*sig^2) )
  
  # compute K^-1
  K = kernelMatrix(kernel = rbfkernel, x = nominalData)
  
  # Tikhonov regularisation if needed
  if (reg =="tikhonov") K = K + param*diag(1, n, n) 

  # tool function for matrix inversion methods-----------------
  mat.inverse = function(mat, type){
    switch(type,
           cholesky = {mat.chol = chol(x = mat); chol2inv(x = mat.chol)},
           svd = {s=svd(mat); D=diag(1/s$d); U=s$u; U %*% D %*% t(U)},
           svdTrunc ={ s=svd(mat); ss <- cumsum(s$d); len <- length(ss); ss <- ss/ss[len];  rk <- min(which(pmin((ss - param), (rep(0, len))) == 0));D <- as.matrix(diag(1/s$d)[1:rk, 1:rk]); U <- s$u[, 1:rk]; U %*% D %*% t(U)},
           moorePenrose = ginv(mat),
           solve = solve(mat)
    ) 
  } # end mat.inverse-------------------------------------------
  
  # matrix inversion
  K.inv = mat.inverse(mat = K,type = inv)
  # verif : OK
  # image(K.inv %*% K)
  # diag(t(K) %*% K.inv %*% K) # theoretical values : all = 1
  
  tau.values = rep(0,n) # gives the tau value associated to set {S\xi}
  #remove one element by one element to update K^{-1}_{S\xi}

  for (i in 1:n){
    #construct correct matrices:  Ki.inv=K^{-1}_{S\xi} = E - F H^-1 G
    F = K.inv[-i,i]; 
    G = K.inv[i,-i];
    H = K.inv[i,i];
    E = K.inv[-i,-i];
    Ki.inv = E - ((F) %*% t(G))/H
    Ki = kernelMatrix(kernel = rbfkernel, x = nominalData[-i,], y = t(nominalData[i,]))
    tau.values[i] = t(Ki) %*% Ki.inv %*% Ki
  }
  #plot(tau.vec, pch=20, ylim=c(0,1))

  return(list(tau.values = tau.values, K.inv=K.inv, sig=sig, nominalData=nominalData))
  
} # end train

##################
### test phase ###
##################

# @input
# - res.train : results from the train() function on nominal data (nominal dataset, K^-1 matrix,  tau vector, sig value)
# - eta : (possibly a vector of) query point(s) to be tested
test = function(res.train, eta){
  # settings
  q=dim(eta)[1] # nb of query points
  rbfkernel = rbfdot(sigma =  1/(2*res.train$sig^2) )
  
  # projection
  Keta = kernelMatrix(kernel = rbfkernel, x = res.train$nominalData, y = eta)
  eta.values = diag(t(Keta) %*% res.train$K.inv %*% Keta)
  
  # p-value
  eta.pval = unlist(lapply(1:q, function(i){
    1-mean(eta.values[i]<res.train$tau.values)
  }))
  return(list(eta.values=eta.values, eta.pval = eta.pval))
}# end test


######################
### general method ###
######################

# @input
# - nominalData : nominal dataset
# - sig : parameter for the RBF kernel, not as in <kernlab:rbfdot> (inverse kernel width parameter) => sigma = 1/(2*sig^2) )
# Param for K^-1 computation:
# - reg : method for regularisation / default :"none" - other : "tikhonov"
# - param : parameter for Tikhonov regularisation method or the % of the total eigenvalues kept in truncated SVD / default : 0
# - inv : method for matrix inverse computation / default : "solve" - other: "svd"; "svdTrunc"; "moorePenrose"; "cholesky"
# - eta : (possibly a vector of) query point(s) to be tested / if NULL : only the tau values are returned

KPCA.AnomalyDetection = function(nominalData, eta=NULL, sig, param = 0, inv = "solve", reg="none") {
  
  res.train = train(nominalData = nominalData, sig = sig,param = param, reg = reg, inv = inv) 
  res.test=NULL
  if (is.null(eta)==FALSE) {
    res.test = test(res.train = res.train,eta = eta)  
  }
  return(list(nominal.values = res.train$tau.values, eta.values=res.test$eta.values, eta.pval = res.test$eta.pval))
  
}



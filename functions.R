######################
### tool functions ###
######################

# for inversion method  
mat.inverse = function(mat, type){
  switch(type,
         cholesky = {mat.chol = chol(x = mat); chol2inv(x = mat.chol)},
         svd = {s=svd(mat); D=diag(1/s$d); U=s$u; U %*% D %*% t(U)},
         svdTrunc ={ s=svd(mat); ss <- cumsum(s$d); len <- length(ss); ss <- ss/ss[len];  rk <- min(which(pmin((ss - param), (rep(0, len))) == 0));D <- as.matrix(diag(1/s$d)[1:rk, 1:rk]); U <- s$u[, 1:rk]; U %*% D %*% t(U)},
         moorePenrose = ginv(mat),
         solve = solve(mat)
  ) 
} # end mat.inverse  
# verif : OK
# mat.inverse(mat=K,type="cholesky")
# mat.inverse(mat=K,type="moorePenrose") 
# mat.inverse(mat=K,type="svdTrunc")   
# mat.inverse(mat=K,type="svd") 
# mat.inverse(mat=K,type="solve") 
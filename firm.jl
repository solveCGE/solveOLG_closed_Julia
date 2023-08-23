# production function and marginal products
function fY(K,L,TFP) TFP.*K.^alpha.*L.^(1-alpha); end
function MPK(K,L,TFP) TFP.*alpha.*(K./L).^(alpha-1); end
function MPKinv(MPKin,L,TFP) L.*(MPKin./(TFP.*alpha)).^(1/(alpha-1)); end
function MPL(K,L,TFP) TFP.*(1-alpha).*(K./L).^alpha; end

# invert the firm problem: finds r for given V (and LS, TFP and taxes)
function rdemand(assetsupply, maxiter = 20, tol = 1e-6, verbose = false)

  K2      = copy(K)
  uck2    = copy(uck)
  r2      = copy(r)
  qTob2   = copy(qTob)
  
  error   = Inf
  iter    = 0
  
  while true
    
    iter           += 1
    error_old      = error
    qTob_old       = qTob2
    
    K2             = assetsupply./qTob2
    #K2[1]          = K[1] #predetermined
    uck2           = MPK.(K2,LD,TFP)
    qTob2          = (1.0.-tauprof).*uck2 .+ tauprof*delta .+ (1-delta)
    
    error = sum(abs.(qTob2.-qTob_old))
    
    if verbose cat(paste("Iteration:\t",iter,"\t\tError:\t",error,"\n")); end
    
    if iter > maxiter
      if verbose println("No convergence!!") end
  	  break
    end
    if error < tol
      if verbose println("Convergence!!") end
      break
    end
    if error > error_old
      if verbose println("Increasing error: stop at previous step!") end
    qTob2 = qTob_old
    break
    end
    
  end
  
  r2[1:(tend-1)] = qTob2[2:tend].-1.0;
  r2[tend] = r2[tend-1];
  
  return r2
end


# computes demographic transition (and updates intervivo transfers accordingly)
function compdemo()
  
  global Nv, Nz, N, Nc, ivv, ivz
  
  # compute demography transition
  for tt in 2:tend
    Nv[1,tt]  = NB[tt]
    for i in 2:nag
      Nv[i,tt] = Nv[i-1,tt-1]*gamv[i-1,tt-1]
    end
  end
  
  Nz[:] = per2coh(Nv) # write INTO variable, do not OVERWRITE it (then no repacking required)
  N[:]  = aggcoh2per(Nz)
  Nc[:] = sum(Nv[1:(fag-1),:],dims=1)
  
  # Compute neutral intervivo-transfers by rescaling received transfers
  for tt in 1:tend
    ivgiven          = -sum(Nv[:,tt].*ivv[:,tt].*(ivv[:,tt].<0))
    ivreceived       = sum(Nv[:,tt].*ivv[:,tt].*(ivv[:,tt].>0))
    
    ivv[ivv[:,tt].>0,tt] = ivv[ivv[:,tt].>0,tt]*(ivgiven/ivreceived)
    
    if abs(sum(ivv[:,tt].*Nv[:,tt]))>1e-10
	  println("ERROR IN RECOMPDEMO: Unbalanced intervivo transfers!")
    end
  end
  
  ivz[:,:] = per2coh(ivv)

end

# computes the further life expectancy
function lifeexpect(gamv)
  
  local nag = length(gamv)
  
  lifeexpectageN = zeroscol(nag)
  for a in (nag-1):-1:1
    lifeexpectageN[a] = (lifeexpectageN[a+1]+1)*gamv[a]+gamv[a+1]*(1-gamv[a])
  end
  
  return lifeexpectageN
end

# main routine that solves the transition path of the full model
function solveOLG(starttime = 1, maxiter = 200, tol = 1e-4, damping_budget = 1.0, damping_assets = 1.0, damping_ab = 1.0, damping_r = 0.5, damping_new_assets = 0.7)
  
  global uck, K, Inv, qTob, Y, w, wz, V, TaxF, Div
  global Cons, LS, A, ab, iv, Nw, Nr, P, tauW, TaxP, Taxl, Rev, CG, Exp, PB
  global edy, edg, edl, eda, ediv, edab, edw
  global HH_nonconvz
  global tauWv, tauWz, tauF, tauC, tauCv, tauCz, taul, taulv, taulz, tauprof, cGv, CG
  global r, rz, rv, abv, abz, LD
  global Av, Consv, lambdav, Savv, dis_totv, ellv, ev, wv, pcv, yv 
  
  println("\nRunning Tatonnement Algorithm for Transition:\n");
  
  tic_loop = DateTime(now())
  
  scaleA            = 1.0;           # initialize
  scaleab           = onesrow(tend); # initialize
  
  #===== demography ======#
  compdemo(); # recomputes demographic transition

  for iter in 1:maxiter
    tic_iter = DateTime(now())

    #===== solve the firm problem for given labor demand ======#
    uck[(starttime+1):tend]     = (r[starttime:(tend-1)]+delta*(1.0.-tauprof[(starttime+1):tend]))./(1.0.-tauprof[(starttime+1):tend])
    K[(starttime+1):tend]       = MPKinv(uck[(starttime+1):tend],LD[(starttime+1):tend],TFP[(starttime+1):tend])
    Inv[starttime:(tend-1)]     = K[(starttime+1):tend] - (1-delta)*K[starttime:(tend-1)]

    Inv[tend]       = delta*K[tend]
    qTob[:]         = (1.0.-tauprof).*MPK.(K,LD,TFP) .+ tauprof*delta .+ (1-delta)

    Y[:]            = fY.(K,LD,TFP)
    w[:]            = MPL.(K,LD,TFP)./(1.0.+tauF)
    wv[:,:]         = kron(w,onescol(nag))
    wz[:,:]         = per2coh(wv)
    V[:]            = qTob.*K
    TaxF[:]         = tauprof.*(Y.-(1.0.+tauF).*w.*LD.-delta*K)
    Div[:]          = Y.-(1.0.+tauF).*w.*LD.-Inv.-TaxF

    #===== solve the households' problem for given prices and tax rates ======#
    HHall(starttime, (iter == 1), scaleA)

    #===== aggregation ======#
    Cons[:]      = aggcoh2per(Consz.*Nz)
    LS[:]        = aggcoh2per(notretz.*ellz.*thetaz.*Nz)
    A[:]         = aggcoh2per(Az.*Nz)
    ab[:]        = aggcoh2per(abz.*Nz)
    iv[:]        = aggcoh2per(ivz.*Nz) # should be 0 by construction
    Nw[:]        = aggcoh2per(notretz.*Nz)
    Nr[:]        = aggcoh2per((1.0.-notretz).*Nz)

    # government budget
    P[:]           = aggcoh2per((1.0.-notretz).*pz.*Nz)
    tauW[:]        = aggcoh2per(tauWz.*notretz.*ellz.*thetaz.*Nz)./LS
    TaxP           = aggcoh2per((1.0.-notretz).*tauWz.*pz.*Nz)
    Taxl           = aggcoh2per(taulz.*Nz)
    Rev[:]         = TaxF+(tauF.*LD+tauW.*LS).*w.+Taxl.+tauC.*Cons.+TaxP
    CG[:]          = sum(cGv.*Nv,dims=1)
    Exp[:]         = CG.+P

    # follow given debt-path
    PB[starttime:(tend-1)]  = DG[starttime:(tend-1)].-DG[(starttime+1):tend]./(1.0.+r[starttime:(tend-1)])
    PB[tend]                = r[tend]*DG[tend]/(1+r[tend])

    #===== excess demands ======# 
    edy[:]       = Inv.+Cons.+CG.-Y
    edg[:]       = Rev.-Exp.-PB
    edl[:]       = LD.-LS
    eda[:]       = DG.+V.-A
    ediv[:]      = -iv
    edab[:]      = aggcoh2per((1.0.-gamz).*Savz.*Nz).-ab
    edw[:]       = 1.0.*edy .+ w.*edl .+ ediv .+ edab .+ edg .+ eda - [eda[2:tend];eda[tend]]'./(1.0.+r) # Walras' Law
  
    # check Walras' Law: this always has to hold (even out of equilibrium)! If not there is something wrong with accounting in the model
    if maximum(abs.(edw[starttime:(tend-1)]))> 1e-10 error("Error: Walras Law does not hold!"); end

    toc_iter = DateTime(now())
    dur_iter = toc_iter - tic_iter
    
    #===== checking error and breaking loop ======# 	
    err             = sum(abs.(edy[starttime:tend]))+sum(abs.(edg[starttime:tend]))+sum(abs.(edl[starttime:tend]))+sum(abs.(eda[starttime:tend]))+sum(abs.(ediv[starttime:tend]))+sum(abs.(edab[starttime:tend]));
    err2            = log(err/tol);
    
    println("Iteration: ", @sprintf("%3d",iter) ,"  scaleA: ", @sprintf("%.6f",scaleA), "  scaleab: ", @sprintf("%.6f",mean(scaleab)), "  non-conv.HH: ", @sprintf("%2d",sum(HH_nonconvz)), "  Time: ",  @sprintf("%5d",dur_iter.value), " ms  log of err/tol: ", @sprintf("%2.8f",err2))
    
    if (err2 < 0.0) 
      println(repeat(" ", 102), "Convergence!\n\n")
      break
    end
    if (iter == maxiter)
      println(repeat(" ", 102), "No Convergence!\n\n")
      break
    end
    
    HH_nonconvz[:,:] .= 0 # reset convergence counter
    
    #======= updating for next iteration =======#
    # budget rules
    budget_surplus  = edg*damping_budget
    
    if (budget_bal == 1)
      tauWv[:,:]     = tauWv .- kron(budget_surplus./(w.*LS),onescol(nag))
      tauWz[:,:]     = per2coh(tauWv)
    end
    if (budget_bal == 2)
      tauF[:]        = tauF .- budget_surplus./(w.*LD)
    end
    if (budget_bal == 3)
      tauC[:]        = tauC .- budget_surplus./Cons
      tauCv[:,:]     = kron(tauC,onescol(nag))
      tauCz[:,:]     = per2coh(tauCv)
    end
    if (budget_bal == 4)
      taul[:]            = taul .- budget_surplus./(N.-Nc)
      taulv[fag:nag,:]   = kron(taul,onescol(nag-fag+1))
      taulz[:,:]         = per2coh(taulv)
    end
    if (budget_bal == 5)
      tauprof[:]    = tauprof - budget_surplus./(Y.-(1.0.+tauF).*w.*LD.-delta*K);
    end
    if (budget_bal == 6)
      cGv[:,:]      = cGv .+ kron(budget_surplus./N,onescol(nag))
      CG[:]         = sum(cGv.*Nv,dims=1)
    end
    
    # price updating
    newassets       = damping_new_assets*(A.-DG) .+ (1-damping_new_assets).*V
    r_new           = rdemand(newassets)
    r[:]            = damping_r*r_new .+ (1-damping_r)*r
    rv[:,:]         = kron(r,onescol(nag))
    rz[:,:]         = per2coh(rv)

    scaleab         = 1.0.+(aggcoh2per((1.0.-gamz).*Savz.*Nz)./ab.-1.0)*damping_ab;
    abv[:,:]        = abv.*kron(scaleab,onescol(nag))
    abz[:,:]        = per2coh(abv)
    LD[:]           = LS
    scaleA          = 1+((DG[starttime]+V[starttime])/A[starttime]-1)*damping_assets;

  end
  
  # convert cohort-view variables back to period-view variables
  # (those where only cohort-view variables were altered in solveOLG)
  Av[:,:]       = coh2per(Az)
  Consv[:,:]    = coh2per(Consz)
  lambdav[:,:]  = coh2per(lambdaz)
  Savv[:,:]     = coh2per(Savz)
  dis_totv[:,:] = coh2per(dis_totz)
  ellv[:,:]     = coh2per(ellz)
  rv[:,:]       = coh2per(rz)
  wv[:,:]       = coh2per(wz)
  pcv[:,:]      = coh2per(pcz)
  yv[:,:]       = coh2per(yz)

  toc_loop = DateTime(now())
  dur_loop = toc_loop-tic_loop
  println("Computation time:\t", @sprintf("%.4f",dur_loop.value/1000), " sec")
  println("CHECK SOLUTION:\t\t", @sprintf("%.16f",sum(abs.(edy).+abs.(edl).+abs.(edg).+abs.(eda).+abs.(ediv).+abs.(edab))))

end
println("Start calibration:\n")

# parameters
delta                    = 0.05                 # depreciation rate
r0                       = 0.04                 # real interest rate
sigma                    = 0.9                  # elasticity of inter-temporal substitution
sigL                     = 0.3                  # labor supply elasticity

# note: ages are off-set by 1 year, e.g. age group 1 contains 0-year olds
fag                      = 14                   # first economically active age-group (age 15)
rag0                     = 61.3                 # retirement age group (retirement age 62.3), non-whole numbers allowed
iag0                     = 51                   # first age group giving inter-vivo transfers

ivpc                     = 0.2                  # intervivo transfer received per capita

# some normalizations
N0                       = 100.0                # population
Y0                       = 100.0                # GDP
L0                       = 30.0                 # total labor supply in efficiency units
w0                       = 2.0                  # wage rate
Cons0                    = 0.55*Y0              # consumption share (calibrated using taul0)

### DEMOGRAPHY ###
for i in 1:nag
  gamv0[i] = 1-0.89^(nag-i+1)                   # some simple profile
end

# survival of last age group is 0
gamv0[nag] = 0
  
# compute demography
Nv0[1]  = 1
for i in 2:nag
  Nv0[i] = Nv0[i-1]*gamv0[i-1]
end

# rescale population
NB0             = 1/sum(Nv0)*N0
Nv0             = Nv0/sum(Nv0)*N0

avage0 = sum(Nv0.*collect(0:(nag-1)))/N0
report("REPORT: Average age:",avage0)
lifeexpect0 = lifeexpect(gamv0)
report("REPORT: Life-expectancy at age 0:", lifeexpect0[1])
report("REPORT: Life-expectancy at age 65:", lifeexpect0[66])

### AGE PROFILES ###

# indicator for not-retired
notretv0[1:floor(Int,rag0)]   .= 1                       # not retired
notretv0[floor(Int,rag0)+1]   =  rag0-floor(rag0)        # partly retired

# intervivo-transfers
ivv0[iag0:nag]            = -seq(ivpc,ivpc*2,nag-iag0+1) # some increasing profile (from ivpc to 2*ivpc)
ivv0[fag:(iag0-1)]        = -sum(ivv0[iag0:nag].*Nv0[iag0:nag])/sum(Nv0[fag:(iag0-1)])*onescol(iag0-fag)

iv0                       = sum(ivv0.*Nv0)
if abs(iv0)>1e-10 
  error("ERROR: UNBALANCED INTERVIVO TRANSFERS!")
end

thetav0                     = zeroscol(nag)                               # labor productivity parameters
theta_peak                  = floor(Int,rag0)-10                          # assumption: productivity peaks 10 years before retirement
thetav0[fag:theta_peak]     = seq(0.7,1.0,theta_peak-fag+1)
thetav0[(theta_peak+1):nag] = seq(1.0,0.1,nag-theta_peak)

ellv0                       = L0/sum(Nv0.*thetav0.*notretv0)*onescol(nag)   # labor supply

# partition of population
Nc0     = sum(Nv0[1:(fag-1)])        # number of children
Nw0 	= sum(notretv0.*Nv0)-Nc0 		  # number of workers
Nr0 	= sum((1.0.-notretv0).*Nv0) 	    # number of retirees
report("REPORT: Old-age dependency ratio:",sum(Nv0[66:nag])/sum(Nv0[16:65]))
report("REPORT: Economic dependency ratio:",(Nc0+Nr0)/Nw0)
report("CHECK: Newborns - deaths:", sum((1.0.-gamv0).*Nv0)-NB0)
report("CHECK: Children + workers + retriees - pop.:", Nc0+Nw0+Nr0-N0)

### POLICY PARAMETERS ###
tauWv0                   = 0.15*onescol(nag)            # wage tax rate worker & retiree
tauF0                    = 0.2                          # payroll tax rate
tauC0                    = 0.2                          # consumption tax rate
tauprof0                 = 0.1                          # profit tax rate
pv0                      = 0.65*sum(w0.*ellv0.*thetav0.*Nv0)/N0*onescol(nag)  # old-age pension (65% of average wage earnings)
DG0                      = 0.6*Y0                       # government debt level (60% of GDP)

# cGv0 is used to balance budget in calibration
cGv0_profile             = 0.2*onescol(nag)
cGv0_profile[1:25]       = seq(0.4,0.2,25)
cGv0_profile[55:nag]     = seq(0.2,1.0,nag-55+1) # some U-shaped profile

# price of consumption and age specific prices and tax rates (but the same for all age groups)
pc0     = 1+tauC0
tauCv0  = tauC0*onescol(nag)
pcv0    = pc0*onescol(nag)
wv0     = w0*onescol(nag)
rv0     = r0*onescol(nag)

LS0     = sum(notretv0.*ellv0.*thetav0.*Nv0)   # aggregate labor supply
LD0     = LS0
uck0    = (r0+delta*(1-tauprof0))/(1-tauprof0) # user-cost of capital
K0      = (Y0-(1+tauF0)*w0*LD0)/uck0
Inv0    = delta*K0
alpha   = K0*uck0/(K0*uck0+LS0*((1+tauF0)*w0))
qTob0   = (1-tauprof0)*alpha*Y0/K0 + tauprof0*delta + (1-delta) # = 1+r0
TFP0    = Y0/((K0^alpha)*(LS0^(1-alpha)))
#LD0     = ((1-alpha)*TFP0/((1+tauF0)*w0))^(1/alpha)*K0 # also true
TaxF0   = tauprof0*(Y0-(1+tauF0)*w0*LD0-(delta*K0))
Div0    = Y0-(1+tauF0)*w0*LD0-Inv0-TaxF0
V0      = (1+r0)*Div0/r0

# MATCH CALIBRATION TARGETS;
xcalib0 = [0.01, 0.3719, 0.40, 13, 1] # starting guesses for nlsolve()

xout = nlsolve(calibfind,xcalib0,ftol=1e-12)
if maximum(abs.(calibfind(xout.zero))) > 1e-6
	error("NEWTON METHOD DID NOT CONVERGE!\n")
end

### CALIBRATION OF LABOR SUPPLY MARGINS ###

# set parl0 in order to reproduce ell0, FOC ell0
parlv0   = (wv0.*(1.0.-tauWv0).*thetav0./pcv0).*(ellv0.^(-1/sigL)).*(Consv0.^(-1/sigma)); parlv0[1:(fag-1)] .= 0
# set parl1 in order to normalize disutility of labor to 0
parlv1   = (sigL/(1+sigL)).*parlv0.*(ellv0.^((1+sigL)/sigL))
dis_totv0= (sigL/(1+sigL)).*parlv0.*(ellv0.^((1+sigL)/sigL)).-parlv1

report("REPORT: Asset-to-output ratio:", A0/Y0)

checkA0         = sum(Av0.*Nv0)-A0
checkAv0        = Av0[nag]+yv0[nag]+ivv0[nag]+abv0[nag]-pc0*Consv0[nag] # end of period assets of last age group are zero
checkN0         = sum(Nv0)-N0

chkcalib = [edy0,edl0,edg0,ediv0,eda0,edab0,checkA0,checkAv0,checkN0]

report("CHECK: Calibration:",sum(chkcalib));

# fill time-dependent variables with calibration values
Cons                          = Cons0*onesrow(tend)
DG                            = DG0*onesrow(tend)
Inv                           = Inv0*onesrow(tend)
LD                            = LD0*onesrow(tend)
LS                            = LS0*onesrow(tend)
K                             = K0*onesrow(tend)
N                             = N0*onesrow(tend)
NB                            = NB0*onesrow(tend)
PB                            = PB0*onesrow(tend)
TFP                           = TFP0*onesrow(tend)
ab                            = ab0*onesrow(tend)
pc                            = pc0*onesrow(tend)
r                             = r0*onesrow(tend)
rag                           = rag0*onesrow(tend)
tauC                          = tauC0*onesrow(tend)
tauF                          = tauF0*onesrow(tend)
tauW                          = tauW0*onesrow(tend)
taul                          = taul0*onesrow(tend)
tauprof                       = tauprof0*onesrow(tend)
uck                           = uck0*onesrow(tend)

# fill time-dependent and age-dependent variables with calibration values
Av                            = kron(Av0, onesrow(tend))
Az                            = kron(Av0, onesrow(ncoh))
Consv                         = kron(Consv0, onesrow(tend))
Consz                         = kron(Consv0, onesrow(ncoh))
Nv                            = kron(Nv0, onesrow(tend))
Nz                            = kron(Nv0, onesrow(ncoh))
Savv                          = kron(Savv0, onesrow(tend))
Savz                          = kron(Savv0, onesrow(ncoh))
abv                           = kron(abv0, onesrow(tend))
abz                           = kron(abv0, onesrow(ncoh))
cGv                           = kron(cGv0, onesrow(tend))
cGz                           = kron(cGv0, onesrow(ncoh))
ellv                          = kron(ellv0, onesrow(tend))
ellz                          = kron(ellv0, onesrow(ncoh))
gamv                          = kron(gamv0, onesrow(tend))
gamz                          = kron(gamv0, onesrow(ncoh))
ivv                           = kron(ivv0, onesrow(tend))
ivz                           = kron(ivv0, onesrow(ncoh))
lambdav                       = kron(lambdav0, onesrow(tend))
lambdaz                       = kron(lambdav0, onesrow(ncoh))
notretv                       = kron(notretv0, onesrow(tend))
notretz                       = kron(notretv0, onesrow(ncoh))
pv                            = kron(pv0, onesrow(tend))
pz                            = kron(pv0, onesrow(ncoh))
tauCv                         = kron(tauCv0, onesrow(tend))
tauCz                         = kron(tauCv0, onesrow(ncoh))
tauWv                         = kron(tauWv0, onesrow(tend))
tauWz                         = kron(tauWv0, onesrow(ncoh))
taulv                         = kron(taulv0, onesrow(tend))
taulz                         = kron(taulv0, onesrow(ncoh))
thetav                        = kron(thetav0, onesrow(tend))
thetaz                        = kron(thetav0, onesrow(ncoh))
rv                            = kron(rv0, onesrow(tend))
rz                            = kron(rv0, onesrow(ncoh))
wv                            = kron(wv0, onesrow(tend))
wz                            = kron(wv0, onesrow(ncoh))

## some optional plots of the calibration
if genplots
  plot(seq(0,nag-1), Av0, label = "", xlabel="age", ylabel="assets")

  plot(seq(0,nag-1), Consv0, label = "consumption", xlabel = "age")
  plot!(seq(0,nag-1), notretv0.*ellv0.*thetav0.*wv0.*(1.0.-tauWv0), label = "net labor income")
  plot!(seq(0,nag-1), (1.0.-notretv0).*pv0.*(1.0.-tauWv0), label = "net pension income")
  plot!(seq(0,nag-1), cGv0, label = "public consumption")
end

;

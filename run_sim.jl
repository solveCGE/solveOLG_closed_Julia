#------------------------------------------------------------------------------#
# "solveOLG_closed_Julia"
# 
# Solves an AK-OLG-Model, closed economy, with income effects in Julia
# Philip Schuster, August, 2023
#
# run_sim.jl: main script, shock definition, then computes transition path
#------------------------------------------------------------------------------#

println("\nSIMPLE AUERBACH-KOTLIKOFF CLOSED ECONOMY MODEL IN JULIA\n")

# control center
const tend            = 300   # number of periods
const nag             = 100   # number of age groups (nag = x => max age = x-1)
const budget_bal      = 3     # budget closing instrument (1.. tauW, 2.. tauF, 3.. tauC, 4.. taul, 5.. tauprof, 6.. cG) to fulfill given debt path
const genplots        = false # generate plots

# initializes all variables globally
include("initdata.jl")

# run calibration routine
include("calib.jl")   # <- change parameters in here

##########################
# POLICY SHOCK SECTION
# ========================
# Note: it is typically easiest to introduce shocks to period-view variables 
#       and then convert them to cohort-view variables using per2coh()

## tax shocks
#tauprof = tauprof*0+0.15                    # profit tax rate is set to 15%
#tauprof[10:tend] = tauprof[10:tend]*0+0.15  # profit tax rate is set to 15%, starting in 10 years (announced today)
#tauprof[1:10] = tauprof[1:10]*0+0.15        # profit tax rate is set to 15% temporarily for next 10 years
#tauWv = tauWv*1.02; tauWz = per2coh(tauWv)  # wage tax is increased by 2%

## delayed retirement (uncomment whole block)
#rag[1:10] = seq(rag0,rag0+2,10); rag[11:tend].=rag0+2 
#notretv[:,:] .= 0
#for tt in 1:tend
#  notretv[1:floor(Int,rag[tt]),tt] .= 1
#  notretv[floor(Int,rag[tt])+1,tt] = rag[tt]-floor(rag[tt])
#end
#notretz = per2coh(notretv)                  # effective retirement age is increased linearly by 2 years over the next 10 years
 
## pension cut
#pv = pv*0.95; pz = per2coh(pv) # pensions are cut by 5%

## productivity shocks
#thetav = thetav*1.02; thetaz = per2coh(thetav)                       # individual productivity increases permanently by 2%
#thetav[:,1:30] = thetav[:,1:30]*1.02; thetaz = per2coh(thetav)       # individual productivity increases by 2% in the next 30 years
#TFP = TFP*1.02                                                       # total factor productivity increases permanently by 2%

## fertility shocks
#NB = NB*1.02             # 2% more newborns every year
NB[1:30] = NB[1:30]*1.02 # 2% more newborns every year over next 50 years

## mortality shocks
gamv[60:nag,:] = 1.0.-(1.0.-gamv[60:nag,:])*0.9; gamv[nag,:].=0; gamz = per2coh(gamv)  # reduction of old-age mortality by 10%

## shock to the initial capital stock
#K[1] = K0*0.99                                                       # 1% of capital stock is lost in first period

## change in targeted debt path
#DG[5:20] = seq(DG0,DG0*0.9,20-5+1); DG[21:tend] = DG0*0.9 # reduce public debt by 10% over 15 years starting in 5 years

##########################

# Solve transition path to new steady state
solveOLG(1, 200, 1e-2)

# some transition plots
if genplots
  plot(seq(0,tend), [r0;r], label = "", xlabel = "time", ylabel = "real interest rate")

  plot(seq(0,tend), [w0;w], label = "", xlabel = "time", ylabel = "wage rate")

  y_min = minimum([N/N0 Y/Y0 Inv/Inv0 Cons/Cons0 CG/CG0 A/A0 P/P0])
  y_max = maximum([N/N0 Y/Y0 Inv/Inv0 Cons/Cons0 CG/CG0 A/A0 P/P0])

  plot(seq(0,tend), [1 N/N0]'*100, label = "population", xlabel = "time", ylabel = "period 0 = 100")
  ylims!(y_min*100,y_max*100)
  plot!(seq(0,tend), [1 Y/Y0]'*100, label = "GDP")
  plot!(seq(0,tend), [1 Inv/Inv0]'*100, label = "investment")
  plot!(seq(0,tend), [1 Cons/Cons0]'*100, label = "consumption")
  plot!(seq(0,tend), [1 CG/CG0]'*100, label = "public consumption")
  plot!(seq(0,tend), [1 A/A0]'*100, label = "aggregate assets")
  plot!(seq(0,tend), [1 P/P0]'*100, label = "pension expend.")

end

;

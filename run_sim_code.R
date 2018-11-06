### Run Simulation Code ###
# To run simulation program:

# 1) Run Model4_FluSim_TrueVE_Program using input file
  trueve<-model4.trueve("inputfile.csv",vaccov=.495)

# 2) Run Model4_FluSim_Program using input file and calling values of true VE
#    from trueve
  model4.flusim("inputfile.csv",ve.tsi=trueve[3,],ve.tmai=trueve[4,])

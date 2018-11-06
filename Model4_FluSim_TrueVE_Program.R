##################################################
### Model 4 Flu Simulation Program for True VE ###
### Assumes random vaccination    
##################################################
model4.trueve<-function(params.file="params.csv",vaccov){
  
  ### Call neccessary packages #############################################
  library(reshape2)
  library(data.table)
  
  ### Read in .csv file with parameters ####################################
  params<-read.csv(params.file, header=T)
  
  ### Input parameters #####################################################
  title<-names(params)[2]
  seed<-params[2,2]                     # seed
  N<-params[15,2]                       # population size
  weeks<-params[13,2]                   # number of weeks in the study
  # joint probability of health status and health awareness
  pi.11<-params[16,2]                   # X=1, U=1 
  pi.10<-params[17,2]                   # X=1, U=0 
  pi.01<-params[18,2]                   # X=0, U=1 
  pi.00<-params[19,2]                   # X=0, U=0 
  
  # probability of vaccination
  alpha.11<-params[20,2]                # vaccine coverage for X=1, U=1 
  alpha.10<-params[21,2]                # vaccine coverage for X=1, U=0 
  alpha.01<-params[22,2]                # vaccine coverage for X=0, U=1 
  alpha.00<-params[23,2]                # vaccine coverage for X=0, U=0 
  
  # probability of influenza infection/FARI in each week
  gamma.01<-params[(25+1):(25+weeks),2] # standard person (V=0,X=1)
  
  # multipliers
  phi.gamma<-params[(25+weeks)+1,2]     # multiplier for gamma when X=0             
  theta.gamma<-params[(25+weeks)+2,2]   # multiplier for gamma when V=1
  
  # probability of NFARI in each week
  z1<-(25+weeks)+3
  beta.01<-params[(z1+1):(z1+weeks),2]  # standard person (V=0,X=1,U=1)     
  
  # multipliers  
  phi.beta<-params[(z1+weeks)+1,2]      # multiplier for beta when X=0             
  theta.beta<-params[(z1+weeks)+2,2]    # multiplier for beta when V=1
  
  # probability of medical visit for FARI
  z2<-(z1+weeks)+3
  delta2.01<-params[z2+1,2]             # standard person (V=0,U=1)
  # multipliers 
  theta.delta2<-params[z2+2,2]          # multiplier for delta2 when V=1
  mu.delta2<-params[z2+3,2]             # multiplier for delta2 when U=0
  
  # probability of medical visit for FARI
  delta1.01<-params[z2+5,2]             # standard person (V=0,U=1)
  
  # multipliers  
  theta.delta1<-params[z2+6,2]          # multiplier for delta2 when V=1
  mu.delta1<-params[z2+7,2]             # multiplier for delta2 when U=0
  
  # probability of testing positive for flu
  tau1<-params[z2+10,2]                # given NFARI
  tau2<-params[z2+9,2]                 # given FARI
  #########################################################################
  
  ### Initialize Simulation ###############################################
  sim<-params[1,2] 
  set.seed(seed)
  
  ve.tn<-c(rep(NA,sim))
  ve.tcc<-c(rep(NA,sim))
  
  ve.tsi<-c(rep(NA,sim))
  ve.tmai<-c(rep(NA,sim))
  
  # Create progress bar to monitor simulation progress
  pb <- txtProgressBar(min = 0, max = sim, style = 3)
  
  ### Begin Simulation ####################################################    
  for(s in 1:sim){
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, s)
    #########################################################################
    
    ### Calculate probabilities based on input parameters ###################
    # probability of influenza infection/FARI in each week
    gamma.00<-gamma.01*phi.gamma              # (V=0,X=0)
    gamma.10<-gamma.01*phi.gamma*theta.gamma  # (V=1,X=0)
    gamma.11<-gamma.01*theta.gamma            # (V=1,X=1)
    
    # cummulative probability of influenza infection/FARI up to week j
    c.gamma.01<-gamma.01
    c.gamma.00<-gamma.00
    c.gamma.10<-gamma.10
    c.gamma.11<-gamma.11
    
    for (j in 2:weeks){
      c.gamma.01[j]<-sum(gamma.01[1:j])
      c.gamma.00[j]<-sum(gamma.00[1:j])
      c.gamma.10[j]<-sum(gamma.10[1:j])
      c.gamma.11[j]<-sum(gamma.11[1:j])
    }
    
    # conditional probability of influenza infection/FARI in week j given no FARI prior
    cond.gamma.01<-gamma.01
    cond.gamma.00<-gamma.00
    cond.gamma.10<-gamma.10
    cond.gamma.11<-gamma.11
    
    for (j in 2:weeks){
      cond.gamma.01[j]<-gamma.01[j]/(1-sum(gamma.01[1:j-1]))
      cond.gamma.00[j]<-gamma.00[j]/(1-sum(gamma.00[1:j-1]))
      cond.gamma.10[j]<-gamma.10[j]/(1-sum(gamma.10[1:j-1]))
      cond.gamma.11[j]<-gamma.11[j]/(1-sum(gamma.11[1:j-1]))
    }
    
    # probability of medical visit for FARI
    delta2.11<-delta2.01*theta.delta2                      # (V=1,U=1)
    delta2.00<-delta2.01*mu.delta2                         # (V=0,U=0)
    delta2.10<-delta2.01*theta.delta2*mu.delta2            # (V=1,U=0)
    
    # probability of NFARI in each week
    beta.00<-beta.01*phi.beta              # (V=0,X=0)
    beta.10<-beta.01*phi.beta*theta.beta   # (V=1,X=0)
    beta.11<-beta.01*theta.beta            # (V=1,X=1)
    
    # probability of medical visit for NFARI
    delta1.11<-delta1.01*theta.delta1                      # (V=1,U=1)
    delta1.00<-delta1.01*mu.delta1                         # (V=0,U=0)
    delta1.10<-delta1.01*theta.delta1*mu.delta1            # (V=1,U=0)
    #########################################################################
    
    ### Expected True VE ###################################################
    # VE against SI
    e.ve.tsi<-1-(sum(gamma.10*(pi.00+pi.01) + gamma.11*(pi.10+pi.11))/sum(gamma.00*(pi.00+pi.01) + gamma.01*(pi.10+pi.11)))
    # VE against MAI    
    e.ve.tmai<-1-(sum((gamma.10*(delta2.10*pi.00 + delta2.11*pi.01)+gamma.11*(delta2.10*pi.10 + delta2.11*pi.11)))/sum((gamma.00*(delta2.00*pi.00 + delta2.01*pi.01)+gamma.01*(delta2.00*pi.10 + delta2.01*pi.11))))
    ##########################################################################
    
    ### Generate Population ##################################################
    # generate population and begin with everyone standard (V=0,X=1,U=1) and susceptible (Y=0)
    id=c(seq(1,N))                        # vector of participant ID numbers
    Y=matrix(c(rep(0,weeks*N)),nrow=N)    # matrix of participant ARI status (Y=0-susceptible, 1-NFARI, 2-FARI)
    # columns are weeks, rows are individuals
    V<-c(rep(0,N))                        # vector of participant vaccination status (all start unvaccinated)
    X<-c(rep(1,N))                        # vector of participant health status (X=0-frail, 1-healthy)
    U<-c(rep(1,N))                        # vector of participant health awareness (U=0-low, 1-high)
    M<-matrix(c(rep(0,weeks*N)),nrow=N)   # matrix of whether participant sought medical care for their ARI (M=0-no, 1-yes)
    # columns are weeks, rows are individuals
    t<-matrix(c(rep(999,weeks*N)),nrow=N) # matrix of whether a participant tested positive for their ARI (M=0-no, 1-yes)
    # columns are weeks, rows are individuals
    
    # intitialize week parameter
    w<-0
    # Assign individuals to have (X,U) - note: everyone starts out with X=1, U=1
    nx1u0<-ceiling(N*pi.10)          # number of people with X=1, U=0
    nx0u1<-ceiling(N*pi.01)          # number of people with X=0, U=1
    nx0u0<-ceiling(N*pi.00)          # number of people with X=0, U=0
    nx1u1<-N-sum(nx1u0,nx0u1,nx0u0)  # number of people with X=1, U=1
    
    px1u0<-sample(id,nx1u0)          # IDs of people with X=1, U=0
    px0u1<-sample(id,nx0u1)          # IDs of people with X=0, U=1
    px0u0<-sample(id,nx0u0)          # IDs of people with X=0, U=0
    
    X[c(px0u1,px0u0)]<-0             # change X from 1 to 0 for people in px0u1 and px0u0
    U[c(px1u0,px0u0)]<-0             # change U from 1 to 0 for people in px1u0 and px0u0 
    
    # Determine vaccinated individuals (depends on X and U) - random vaccination
    # Vaccinations prior to the season
    #nv.00<-ceiling(nx0u0*alpha.00)   # (X=0,U=0)
    #nv.01<-ceiling(nx0u1*alpha.01)   # (X=0,U=1)
    #nv.10<-ceiling(nx1u0*alpha.10)   # (X=1,U=0)
    #nv.11<-ceiling(nx1u1*alpha.11)   # (X=1,U=1)
    nv<-ceiling(vaccov*N)
    
    # select IDs to be vaccinated
    #pv.00<-sample(px0u0,nv.00)       # (X=0,U=0)
    #pv.01<-sample(px0u1,nv.01)       # (X=0,U=1) 
    #pv.10<-sample(px1u0,nv.10)       # (X=1,U=0) 
    #pv.11<-sample(id[X==1 & U==1],nv.11) # (X=1,U=1)
    pv<-sample(id,nv)
    
    # assign IDs to be vaccinated
    #V[c(pv.00,pv.01,pv.10,pv.11)]<-1
    V[pv]<-1
    
    # set week of vaccination to 0
    wV<-c(rep(999,N))                    # vector of week of vaccination (999-not vaccinated)
    #wV[c(pv.00,pv.01,pv.10,pv.11)]<-0
    wV[pv]<-0
    # Determine who will get FARI
    flu.ind<-c(rep(999,N))               # empty vector
    for (i in 1:N){ 
      #cat("Person ",i,"\n")
      #rand.num<-runif(1)                 # select random uniform number
      for (j in 1:weeks){
        rand.num<-runif(1)
        # if rand.num<=gamma_{1vx} then person i will have an FARI in week j
        if (V[i]==0 & X[i]==0 & rand.num<=cond.gamma.00[j]){flu.ind[i]<-j; break} # skip to next person if current person gets flu in week j
        if (V[i]==0 & X[i]==1 & rand.num<=cond.gamma.01[j]){flu.ind[i]<-j; break}
        if (V[i]==1 & X[i]==0 & rand.num<=cond.gamma.10[j]){flu.ind[i]<-j; break}
        if (V[i]==1 & X[i]==1 & rand.num<=cond.gamma.11[j]){flu.ind[i]<-j; break}
      }  
      
      # if rand.num<=gamma_{jvx} & rand.num>gamma_{(j-1)vx} then person i will have an FARI in week j
      #for (j in 2:weeks){
      #  if (V[i]==0 & X[i]==0 & rand.num<=c.gamma.00[j] & rand.num>c.gamma.00[j-1]){flu.ind[i]<-j}
      #  if (V[i]==0 & X[i]==1 & rand.num<=c.gamma.01[j] & rand.num>c.gamma.01[j-1]){flu.ind[i]<-j}
      #  if (V[i]==1 & X[i]==0 & rand.num<=c.gamma.10[j] & rand.num>c.gamma.10[j-1]){flu.ind[i]<-j}
      #  if (V[i]==1 & X[i]==1 & rand.num<=c.gamma.11[j] & rand.num>c.gamma.11[j-1]){flu.ind[i]<-j}
      #  if(flu.ind[i]==j){break}       # skip to next person if current person gets flu in week j
      #}  
      #}                                  # end loop over weeks
    }                                    # end FARI loop over people
    ##########################################################################
    # Empty vectors to indicate cases and controls  
    cases.si<-c(rep(0,N))                 # cases of SI
    cases.mai<-c(rep(0,N))                # cases of MAI
    # Empty "switch" vectors
    switch1<-c(rep(0,N))                  # Switch 1 - first NFARI (switches to 1 once a person has first NFARI)          
    switch2<-c(rep(0,N))                  # Switch 2 - first (only) FARI (switches to 1 once a person has FARI)
    switch3<-c(rep(0,N))                  # Switch 3 - first time seeking care (switches to 1 once a person first seeks medical care)
    ### Begin main simulation ################################################
    # Simulate over weeks
    for (j in 1:weeks){
      # Simulate over people
      for (i in 1:N){
        # FARI
        if(flu.ind[i]==j){
          Y[i,j]<-2                   # Set Y to 2
          cases.si[i]<-1
          switch2[i]<-1               # turn Switch 2 on
          
          # determine if person seeks care for FARI
          if (Y[i,j]==2){
            rand.num2<-runif(1)     # generate random number
            if (V[i]==0 & U[i]==0 & rand.num2<=delta2.00){M[i,j]<-1; cases.mai[i]<-1}
            if (V[i]==0 & U[i]==1 & rand.num2<=delta2.01){M[i,j]<-1; cases.mai[i]<-1}
            
            if (V[i]==1 & U[i]==0 & rand.num2<=delta2.10){M[i,j]<-1; cases.mai[i]<-1}
            if (V[i]==1 & U[i]==1 & rand.num2<=delta2.11){M[i,j]<-1; cases.mai[i]<-1}
          }
        } # end if flu.ind==j
      } # end loop over people
    } # end loop over weeks
    
    # Risk ratio
    a.si<-sum(cases.si[cases.si==1 & V==1])/(vaccov*N)     # True case of SI (V=1)
    b.si<-sum(cases.si[cases.si==1 & V==0])/(N-vaccov*N)   # True case of SI (V=0)
    
    a.mai<-sum(cases.mai[cases.mai==1 & V==1])/(vaccov*N)  # True case of MAI (V=1)
    b.mai<-sum(cases.mai[cases.mai==1 & V==0])/(N-vaccov*N)# True case of MAI (V=0)
    
    ve.tsi[s]<-1-(a.si/b.si)      # True VE against SI
    ve.tmai[s]<-1-(a.mai/b.mai)   # True VE against MAI
    
    ### End of simulation ####################################################
  }
  
  # Close progress bar
  #  close(pb)
  ##########################################################################
  
  ### Output ###############################################################
  # Print True VE and Estimates
  cat("\n")
  ve.true<-data.frame(E_VE_TSI=e.ve.tsi,E_VE_TMAI=e.ve.tmai,
                      VE_TSI=ve.tsi,VE_TMAI=ve.tmai)
  ve.true.mean<-apply(ve.true,2,mean) # Mean VE over all simulations
  ve.true.mean<-data.frame(ve.true.mean)
  return(ve.true.mean)
  ##########################################################################
  
  ### End function #########################################################
}
##########################################################################  

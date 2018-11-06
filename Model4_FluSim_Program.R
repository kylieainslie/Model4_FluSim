######################################
### Model 4 Flu Simulation Program ###
######################################
model4.flusim<-function(params.file="params.csv",ve.tsi,ve.tmai){
  
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
  
  # Empty vectors to store values in each simulation
  # frequencies
  a.asc<-c(rep(NA,sim))
  b.asc<-c(rep(NA,sim))
  a.psc<-c(rep(NA,sim))
  b.psc<-c(rep(NA,sim))
  a<-c(rep(NA,sim))
  b<-c(rep(NA,sim))
  c.tn<-c(rep(NA,sim))
  d.tn<-c(rep(NA,sim))
  c.tcc<-c(rep(NA,sim))
  d.tcc<-c(rep(NA,sim))
  
  # odds ratios and VE estimates
  ve.asc<-c(rep(NA,sim))
  ve.psc<-c(rep(NA,sim))
  or.tn<-c(rep(NA,sim))
  or.tcc<-c(rep(NA,sim))
  ve.tn<-c(rep(NA,sim))
  ve.tcc<-c(rep(NA,sim))
  # observed proportion of vaccination
  obs.vac<-c(rep(NA,sim))
  
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
    
    # Determine vaccinated individuals (depends on X and U)
    # Vaccinations prior to the season
    nv.00<-ceiling(nx0u0*alpha.00)   # (X=0,U=0)
    nv.01<-ceiling(nx0u1*alpha.01)   # (X=0,U=1)
    nv.10<-ceiling(nx1u0*alpha.10)   # (X=1,U=0)
    nv.11<-ceiling(nx1u1*alpha.11)   # (X=1,U=1)
    
    # select IDs to be vaccinated
    pv.00<-sample(px0u0,nv.00)       # (X=0,U=0)
    pv.01<-sample(px0u1,nv.01)       # (X=0,U=1) 
    pv.10<-sample(px1u0,nv.10)       # (X=1,U=0) 
    pv.11<-sample(id[X==1 & U==1],nv.11) # (X=1,U=1)
    
    # assign IDs to be vaccinated
    V[c(pv.00,pv.01,pv.10,pv.11)]<-1
    
    # Observed vaccination coverage
    obs.vac[s]<-sum(nv.00,nv.01,nv.10,nv.11)
    # set week of vaccination to 0
    wV<-c(rep(999,N))                 # vector of week of vaccination (999-not vaccinated)
    wV[c(pv.00,pv.01,pv.10,pv.11)]<-0
    
    # Determine who will get FARI
    flu.ind<-c(rep(999,N))            # empty vector
    for (i in 1:N){ 
      #cat("Person ",i,"\n")
      #rand.num<-runif(1)             # select random uniform number
      
      #cat("Week ",j,"\n")
      for (j in 1:weeks){
        rand.num<-runif(1)
        # if rand.num<=gamma_{1vx} then person i will have an FARI in week 1
        if (V[i]==0 & X[i]==0 & rand.num<=cond.gamma.00[j]){flu.ind[i]<-j; break} # skip to next person if current person gets flu in week j
        if (V[i]==0 & X[i]==1 & rand.num<=cond.gamma.01[j]){flu.ind[i]<-j; break}
        if (V[i]==1 & X[i]==0 & rand.num<=cond.gamma.10[j]){flu.ind[i]<-j; break}
        if (V[i]==1 & X[i]==1 & rand.num<=cond.gamma.11[j]){flu.ind[i]<-j; break}
      }
      #for (j in 2:weeks){
      #  if rand.num<=gamma_{jvx} & rand.num>gamma_{(j-1)vx} then person i will have an FARI in week j
      #   if (V[i]==0 & X[i]==0 & rand.num<=c.gamma.00[j] & rand.num>c.gamma.00[j-1]){flu.ind[i]<-j; break}
      #   if (V[i]==0 & X[i]==1 & rand.num<=c.gamma.01[j] & rand.num>c.gamma.01[j-1]){flu.ind[i]<-j; break}
      #   if (V[i]==1 & X[i]==0 & rand.num<=c.gamma.10[j] & rand.num>c.gamma.10[j-1]){flu.ind[i]<-j; break}
      #   if (V[i]==1 & X[i]==1 & rand.num<=c.gamma.11[j] & rand.num>c.gamma.11[j-1]){flu.ind[i]<-j; break}
      # }                                  # end loop over weeks
    }                                    # end FARI loop over people
    ##########################################################################
    # take random subsets of the population to be the two cohorts
    # ASC cohort
    asc.cohort<-sample(id,5000)
    asc.vac<-sum(V[asc.cohort])
    # PSC cohort
    psc.cohort<-sample(id[-asc.cohort],10000) # sample from people not in ASC cohort
    psc.vac<-sum(V[psc.cohort])
    # create indicator for participation in a cohort study (0=no, 1=yes)
    cohort<-ifelse(id%in%c(asc.cohort,psc.cohort),1,0)
    
    ##########################################################################
    # Empty vectors to indicate cases and controls in case-control studies
    cases<-c(rep(0,N))                    # cases (1=ASC, 2=PSC, 3=case-control)
    tncontrols<-c(rep(0,N))               # TN controls
    tcccontrols<-c(rep(0,N))              # TCC controls
    true.cases.si<-c(rep(0,N))            # True cases of SI
    true.cases.mai<-c(rep(0,N))           # True cases of MAI
    # Empty "switch" vectors
    switch1<-c(rep(0,N))                  # Switch 1 - first NFARI (switches to 1 once a person has first NFARI)          
    switch2<-c(rep(0,N))                  # Switch 2 - first (only) FARI (switches to 1 once a person has FARI)
    switch3<-c(rep(0,N))                  # Switch 3 - first time seeking care (switches to 1 once a person first seeks medical care)
    # Indicator of being included in a study as either a case or control
    included<-c(rep(0,N))
    ### Begin main simulation ################################################
    # Simulate over weeks
    for (j in 1:weeks){
      # Simulate over people
      for (i in 1:N){
        if(switch2[i]==1){next} # if a person has already tested positive for FARI skip
        if(switch1[i] & switch3[i]==1){next} # if a person has already sought medical care and tested negative skip
        if(included[i]>0){next} # if a person has already been included as a case or control skip
        ### FARI ###
        if(flu.ind[i]==j){
          Y[i,j]<-2                   # Set Y to 2
          true.cases.si[i]<-1         # True cases
          # test individuals in ASC cohort
          if (Y[i,j]==2 & i %in% asc.cohort){
            rand.num1<-runif(1)      # generate random number
            if (rand.num1<tau2){t[i,j]<-1
            switch2[i]<-1               # turn Switch 2 on
            }
            if (rand.num1>tau2){t[i,j]<-0
            switch1[i]<-1
            }              
          }
          # determine if person seeks care for FARI
          if (Y[i,j]==2){
            rand.num2<-runif(1)     # generate random number
            if (V[i]==0 & U[i]==0 & rand.num2<=delta2.00){M[i,j]<-1}
            if (V[i]==0 & U[i]==1 & rand.num2<=delta2.01){M[i,j]<-1}
            
            if (V[i]==1 & U[i]==0 & rand.num2<=delta2.10){M[i,j]<-1}
            if (V[i]==1 & U[i]==1 & rand.num2<=delta2.11){M[i,j]<-1}
          }
          # determine if person tests positive for FARI
          if (M[i,j]==1 & i %in% id[-asc.cohort]){
            switch3[i]<-1             # turn Switch 3 on
            true.cases.mai[i]<-1      # True case of MAI
            rand.num2a<-runif(1)      # generate random number
            if (rand.num2a<tau2){t[i,j]<-1
            switch2[i]<-1}
            if (rand.num2a>tau2){t[i,j]<-0
            switch1[i]<-1}
          }
        } # end flu.ind==j
        
        ### NFARI ###
        if(flu.ind[i]!=j){
          rand.num3<-runif(1)          # generate random number
          # if rand.num<=gamma_{1vx} then person i will have an FARI in week 1
          if (V[i]==0 & X[i]==0 & rand.num3<=beta.00[j]){Y[i,j]<-1}
          if (V[i]==0 & X[i]==1 & rand.num3<=beta.01[j]){Y[i,j]<-1}
          if (V[i]==1 & X[i]==0 & rand.num3<=beta.10[j]){Y[i,j]<-1}
          if (V[i]==1 & X[i]==1 & rand.num3<=beta.11[j]){Y[i,j]<-1}
          # test ASC cohort members for NFARI - false positives will be considered cases
          if (Y[i,j]==1 & i%in%asc.cohort){
            rand.num1a<-runif(1)
            if (rand.num1a<tau1){t[i,j]<-1
            switch2[i]<-1}
            if (rand.num1a>tau1){t[i,j]<-0
            switch1[i]<-1}              
          }
          # determine if person seeks care for NFARI
          if(Y[i,j]==1){
            switch1[i]<-1             # turn Switch 1 on
            rand.num4<-runif(1) # generate random number
            if (V[i]==0 & U[i]==0 & rand.num4<=delta1.00){M[i,j]<-1}
            if (V[i]==0 & U[i]==1 & rand.num4<=delta1.01){M[i,j]<-1}
            
            if (V[i]==1 & U[i]==0 & rand.num4<=delta1.10){M[i,j]<-1}
            if (V[i]==1 & U[i]==1 & rand.num4<=delta1.11){M[i,j]<-1}
          }
          # determine if person with NFARI tests positive for FARI (false positive)
          if (M[i,j]==1 & i %in% id[-asc.cohort]){
            switch3[i]<-1
            rand.num4a<-runif(1)
            if (rand.num4a<tau1){t[i,j]<-1
            switch2[i]<-1}
            if (rand.num4a>tau1){t[i,j]<-0
            switch1[i]<-1}
          }  
        } # end of if flu.ind!=j
        
        ### Determine cases and controls from test-results 
        # ASC cases
        if (i %in% asc.cohort & t[i,j]==1) {cases[i]<-1; included[i]=1}
        # PSC cases
        if (i %in% psc.cohort & M[i,j]==1 & t[i,j]==1){cases[i]<-2; included[i]=2}
        # case-control study cases
        if (cohort[i]==0){
          if (M[i,j]==1 & t[i,j]==1 & tncontrols[i]!=1){cases[i]<-3; included[i]=3}
          # TN controls
          if (M[i,j]==1 & t[i,j]==0 & tncontrols[i]!=1 & cases[i]<1){tncontrols[i]<-1; included[i]=4}
        }
      } # end loop over people
    } # end loop over weeks
    
    # TCC controls (no ARI throughout study & matched by X to cases)
    noari.x0<-sample(which(cohort==0 & switch1==0 & switch2==0 & included==0 & X==0),length(cases[cases==3 & X==0]))
    noari.x1<-sample(which(cohort==0 & switch1==0 & switch2==0 & included==0 & X==1),length(cases[cases==3 & X==1]))
    tcccontrols[c(noari.x0,noari.x1)]<-1
    
    # don't match by X
    #noari<-sample(which(switch1==0 & switch2==0),length(cases[cases==1]))
    #tcccontrols[noari]<-1
    
    # Final 2x2 table (summed over all weeks)
    a.asc[s]<-sum(cases[cases==1 & V==1])/asc.vac                           # ASC Case (V=1)
    b.asc[s]<-sum(cases[cases==1 & V==0])/(length(asc.cohort)-asc.vac)      # ASC Case (V=0)
    a.psc[s]<-length(cases[cases==2 & V==1])/psc.vac                        # PSC Case (V=1)
    b.psc[s]<-length(cases[cases==2 & V==0])/(length(psc.cohort)-psc.vac)   # PSC Case (V=0)
    
    a[s]<-length(cases[cases==3 & V==1])                     # Case (V=1)
    b[s]<-length(cases[cases==3 & V==0])                     # Case (V=0)
    c.tn[s]<-sum(tncontrols[tncontrols==1 & V==1])           # TN Control (V=1)
    d.tn[s]<-sum(tncontrols[tncontrols==1 & V==0])           # TN Control (V=0)
    c.tcc[s]<-sum(tcccontrols[tcccontrols==1 & V==1])        # TCC Control (V=1)
    d.tcc[s]<-sum(tcccontrols[tcccontrols==1 & V==0])        # TCC Control (V=0)
    
    # Odds Ratio
    or.tn[s]<-(a[s]*d.tn[s])/(b[s]*c.tn[s])                  # OR for TN study
    or.tcc[s]<-(a[s]*d.tcc[s])/(b[s]*c.tcc[s])               # OR for TCC study
    
    # VE Estimates
    ve.asc[s]<-1-(a.asc[s]/b.asc[s])                         # VE for ASC study for sim s
    ve.psc[s]<-1-(a.psc[s]/b.psc[s])                         # VE for PSC study for sim s
    ve.tn[s]<-1-or.tn[s]                                     # VE for TN study for sim s
    ve.tcc[s]<-1-or.tcc[s]                                   # VE for TCC study for sim s
    
    ### Study dataset ########################################################
    # wide dataframe
    mydata<-data.frame(c(rep(s,N)),id,V,X,U,Y,M,t,cases,tncontrols,tcccontrols)
    Ynames<-paste("Y",seq(1:weeks),sep="")
    Mnames<-paste("M",seq(1:weeks),sep="")
    tnames<-paste("T",seq(1:weeks),sep="")
    names(mydata)<-c("Sim","ID","V","X","U",Ynames,Mnames,tnames,"Cases","TN Control","TCC Controls")
    
    # Convert dataframe from wide to long
    long<- melt(setDT(mydata), measure.vars=list(c(6:(6+weeks-1)), 
                                                 c((6+weeks):(6+(2*weeks)-1)),
                                                 c((6+(2*weeks)):(6+(3*weeks)-1))), 
                variable.name='Week', value.name=c('Y', 'M', 'T'))
    long1<-long[order(long$Sim,long$ID),]
    long2<-data.frame(long1[,1:5],long1[,9:12],long1[,6:8])
    
    # Output outcomes files - write to .csv files
    outcomes.ind<-params[9,2]
    timestamp.ind<-params[10,2]
    if (outcomes.ind==1 & s==1){ # include column names
      if(timestamp.ind==0){
        outcomes_wide<-write.table(mydata,file=paste("outcomes_wide_",params[1,2],"_",".csv",sep=""),append=T,sep=",",row.names=F)
        outcomes_long<-write.table(long2,file=paste("outcomes_",params[1,2],"_",".csv",sep=""),append=T,sep=",",row.names=F)
      }
      if(timestamp.ind==1){ # add timestamp to outcomes file name
        outcomes_wide<-write.table(mydata,file=paste("outcomes_wide_",params[1,2],"_",Sys.Date(),".csv",sep=""),append=T,sep=",",row.names=F)
        outcomes_long<-write.table(long2,file=paste("outcomes_",params[1,2],"_",Sys.Date(),".csv",sep=""),append=T,sep=",",row.names=F)
      }  
    }
    if (outcomes.ind==1 & s>1){ # exclude column names
      if(timestamp.ind==0){
        outcomes_wide<-write.table(mydata,file=paste("outcomes_wide_",params[1,2],"_",".csv",sep=""),append=T,sep=",",row.names=F,col.names=F)
        outcomes_long<-write.table(long2,file=paste("outcomes_",params[1,2],"_",".csv",sep=""),append=T,sep=",",row.names=F,col.names=F)
      }
      if(timestamp.ind==1){ # add timestamp to outcomes file name
        outcomes_wide<-write.table(mydata,file=paste("outcomes_wide_",params[1,2],"_",Sys.Date(),".csv",sep=""),append=T,sep=",",row.names=F,col.names=F)
        outcomes_long<-write.table(long2,file=paste("outcomes_",params[1,2],"_",Sys.Date(),".csv",sep=""),append=T,sep=",",row.names=F,col.names=F)
      }  
    }
    
    ##########################################################################
    
    ### End of simulation ####################################################
  }
  
  # Close progress bar
  close(pb)
  ##########################################################################
  ### Output ###############################################################
  # Print mean vaccination coverage (for calculation of true VE)
  cat("\n","Observed Vaccination Coverage","\n",mean(obs.vac)/N,"\n")  
  # Print True VE and Estimates
  cat("\n","True VE against SI","\n",ve.tsi,"\n","True VE against MAI","\n",ve.tmai,"\n")  
  
  ve.all<-data.frame(VE_ASC=ve.asc,VE_PSC=ve.psc,VE_TN=ve.tn,VE_TCC=ve.tcc)
  ve.mean<-apply(ve.all,2,mean) # Mean VE over all simulations 
  
  cohort.table<-data.frame(A_ASC=a.asc,B_ASC=b.asc,A_PSC=a.psc,B_PSC=b.psc)
  table2x2<-data.frame(A=a,B=b,C_TN=c.tn,D_TN=d.tn,C_TCC=c.tcc,D_TCC=d.tcc)
  table2x2.mean<-apply(table2x2,2,mean)
  table2x2.prob<-data.frame(A=table2x2.mean[1]/N,B=table2x2.mean[2]/N,
                            C_TN=table2x2.mean[3]/N,D_TN=table2x2.mean[4]/N,
                            C_TCC=table2x2.mean[5]/N,D_TCC=table2x2.mean[6]/N)
  
  out<-data.frame(VE_ASC=ve.mean[1],VE_PSC=ve.mean[2],VE_TN=ve.mean[3],VE_TCC=ve.mean[4],
                  Bias_ASC_SI=ve.mean[1]-ve.tsi,Bias_PSC_SI=ve.mean[2]-ve.tsi,Bias_TN_SI=ve.mean[3]-ve.tsi,Bias_TCC_SI=ve.mean[4]-ve.tsi,
                  Bias_ASC_MAI=ve.mean[1]-ve.tmai,Bias_PSC_MAI=ve.mean[2]-ve.tmai,Bias_TN_MAI=ve.mean[3]-ve.tmai,Bias_TCC_MAI=ve.mean[4]-ve.tmai)
  return(list(cohort.table,table2x2.prob,out))
  ##########################################################################
  
  ### End function #########################################################
}
##########################################################################   

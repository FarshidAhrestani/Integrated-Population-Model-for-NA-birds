# JAGS-code-for-IPM-model-for-NA-birds
JAGS code used to model population dynamics of N.American birds by integrating Breeding Bird Survey (BBS) and Monitoring Avian Productivity and Survivorship (MAPS) data
    model{
    
    #######################################################################
    # Integrated Population Model for BBS & MAPS data                     #
    #                                                                     #
    # this version has spatially-constant demographic rates               #
    #                                                                     #
    #######################################################################
    
    #######################################################################
    # 1. Define priors
    #######################################################################
    
    #----------------------------------------------------------------------
    # Priors for BBS model
    #----------------------------------------------------------------------
    
    # observer effects
    for(i in 1 : nobservers) {
    	obs[i] ~ dnorm(0.0,tauobs)
    } #i
    
    # start-up effects
    eta ~ dnorm(0.0,1.0E-6)
    tauobs ~ dgamma(0.001,0.001)
    sdobs <- 1 / pow(tauobs, 0.5)

    # initial population size
    for (s in 1:nstrata){
    	n1[s] ~ dunif(0,50)
    	ntot[1,s] ~ dpois(n1[s])
    } #s
    
    # overdispersion
    for (k in 1:ncounts){
    	noise[k] ~ dnorm(0.0, taunoise)
    }

    taunoise ~ dgamma(0.001,0.001)
    sdnoise <- 1 / pow(taunoise, 0.5)

    #-----------------------------------------------------------------------
    # Priors for MAPS CJS model 
    #-----------------------------------------------------------------------

    #### year-specific intercepts for survival (phi), residency (pi),
    #### prob.of predetermining a resident (rho), and recapture prob. (p)
    for(t in 1:nyears-1){
    	p0[t] ~ dunif(0,1) # Priors for mean recapture probability
    	lp0[t] <- log(p0[t]/(1-p0[t])) # Logit transformation of recapture probability 
    	phi0[t] ~ dunif(0,1) # Priors for mean survival probability
    	lphi0[t] <- log(phi0[t]/(1-phi0[t])) # Logit transformation of survival probability
    	rho0[t] ~ dunif(0,1) # Priors for residency probability                       
    	lrho0[t] <- log(rho0[t]/(1-rho0[t]))
    	pi0[t] ~ dunif(0,1) 
    	lpi0[t] <- log(pi0[t]/(1-pi0[t]))
    }
    
    #### random station effects for p model ####
    sigma.p ~ dunif(0, 10)
    tau.p  <- pow(sigma.p, -2)
    for(j in 1:nsta){
    	alpha[j] ~ dnorm(0,tau.p)
    }
    
    #------------------------------------------------------------
    # Recruitment (not informed by data here) 
    #------------------------------------------------------------
    
    for (t in 1:nyears-1) {
    	f[t] ~ dunif(0,10)
    }
    
    #############################################################
    # 2. Likelihood for the BBS data 
    #############################################################
    
    for (k in 1:ncounts){
    	log(lambda[k]) <- log(ntot[year[k],strat[k]]) + obs[obser[k]] + eta*firstyr[k] + noise[k]
    	count[k] ~ dpois(lambda[k])
    }
    
    ############################################################
    # 3. Likelihood for the CJS-MAPS data
    ############################################################
    
    for(i in 1:nind){
    	for(t in 1:first[i]){
    		## Define latent state at first capture
    		z[i,t] ~ dbern(1)
    	}
    	logit(rho[i]) <- lrho0[first[i]]
    	R[i] ~ dbern(pi[i,first[i]])
    	mu[i] <- R[i]*rho[i]
    	r[i] ~ dbern(mu[i])
    	for(t in first[i]:nyears-1){
    		logit(p[i,t]) <- lp0[t]  + alpha[sta[i]] 
    		logit(pi[i,t]) <- lpi0[t] 
    		logit(phi[i,t]) <- lphi0[t] 
    	}
    	## State process
    	for (t in first[i]+1:nyears){
    		mu2[i,t] <- z[i,t-1]*phi[i,t-1]*R[i]
    		z[i,t] ~ dbern(mu2[i,t])
    		# Observation process
    		mu1[i,t] <- p[i,t-1]*z[i,t]
    		y[i,t] ~ dbern(mu1[i,t])
    	}
    }

    ############################################################
    # 3. Population process model
    ############################################################

    for (s in 1:nstrata){
    	for (t in 2:nyears){
    		#ntot[t,s] ~ dpois(mn.tot[t,s]) 
    		#mn.tot[t,s] <- nsurv[t,s] + nrecr[t,s]
    		ntot[t,s] <- nsurv[t,s] + nrecr[t,s]
    		nsurv[t,s] ~ dbin(phi0[t-1], ntot[t-1,s])
    		nrecr[t,s] ~ dpois(mn.recr[t,s])
    		mn.recr[t,s] <- f[t-1]*ntot[t-1,s]
    	}
    }

    ############################################################
    # 5. Derived parameters
    ############################################################

    totareaweight <- sum(areaweight[1:nstrata]) 
    for (i in 1:nstrata){
    	for( t in 1 : nyears ) {
    		n[i,t] <- nonzeroweight[i]*exp(log(ntot[t,i]))
    		# + 0.5*sdnoise*sdnoise + 0.5*sdobs*sdobs
    		N[i,t] <- areaweight[i]*n[i,t]/totareaweight
    	}
    }
    for( i in 1 : nstrata ) {
    	B[i] <- pow(n[i,nyears]/n[i,1],1/(nyears-1))
    }
    for( t in 1 : nyears ) {
    	CompIndex[t] <- sum(N[1:nstrata,t])
    }
    
    Bbar <- pow(CompIndex[nyears]/CompIndex[1],1/(nyears-1))

    }

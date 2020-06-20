##MEND with two P pools
rm(list = ls ())

library(deSolve)
##solving MEND model with B, D, P, Q, M, EP, EM and IC (CO2)

Decom <- function(times, states, parameters, flux_function) {
  with(as.list(c(states, parameters)), {
    F1 <- (1/Ec) * (Vd + mR) * B *D /( kD + D)  ## microbial uptake 
    F2p1 <- Vp1* EP1 * P1 / (kP1 + P1) # POM decomposition
    F2p2 <- Vp2* EP2 * P2 / (kP2 + P2)
    F3 <- Vm * EM * M / (kM + M) #MOAM decomposition
    F4 <- (1/Ec -1) * Vd * B * D /( kD + D) #microbial growth respiration 
    F5 <- (1/Ec -1) * mR * B * D /( kD + D) ## microbial maintenance respiration
    F6 <- Kads * D *(1- Q/ Qmax) ##adsorption
    F7 <- Kdes * Q/Qmax  ##desorption
    F8 <- (1- pEP-pEP- pEM) * mR * B   #microbial biomass decay
    F9ep1 <- P1/(P1+P2)*pEP * mR * B  #enzyme production
    F9ep2 <- P2/(P1+P2)*pEP * mR * B  
    F9em <- pEM * mR *B
    F10ep1 <- rEP * EP1  #enzyme decay
    F10ep2 <- rEP * EP2
    F10em <- rEM *EM
    
    dP1 <- Ip1 + (1 - gD) * F8- F2p1
    dP2 <- Ip2 - F2p2
    dM <- (1 - fD) * (F2p1+F2p2) - F3
    dB <- F1- (F4 + F5) - F8 - (F9em + F9ep1+F9ep2)
    dD <- Id + fD * (F2p1+F2p2) + gD * F8 + F3 + (F10em + F10ep1+F10ep2)- F1 - (F6 - F7)
    dQ <- F6 - F7
    dEP1 <- F9ep1 - F10ep1
    dEP2 <- F9ep2 - F10ep2
    dEM <- F9em - F10em
    dIC <- F4 + F5   #CO2 fllux
    dTot <- Ip1 + Ip2 + Id -(F4 + F5) 
    return(list(c(dP1, dP2, dM, dB, dD, dQ, dEP1,dEP2, dEM, dIC, dTot)))
  })
}



times <- seq(0, 8760, 24) ##per hour 24hour by 365 days =8760

B<-0.2
D<-0.5
P1<-4
P2<-2
Q<-0.3
M<-10
EP1<-0.00001
EP2<-0.00001
EM<- 0.00001
IC<- 0
Tot<- B+D+P1+P2+M+Q+EP1+EP2+EM

states <- c(P1 = P1, P2=P2,  M = M, B=B, D=D, Q = Q,EP1 = EP1, EP2=EP2,EM = EM, IC = IC, Tot= Tot) 
parameters <- c(Ec = 0.47, Vd = 5e-4, mR = 2.8e-4, kD = 0.26, Vp1 = 2.5,  kP1 = 50, Vp2=2.5, kP2=5, 
                Vm = 0.1, kM = 250,Kads = 0.006, Kdes = 0.001, Qmax = 1.7,
                pEP=0.01, pEM=0.01, rEP=1e-3, rEM=1e-3, 
                Ip1 = 4e-5, Ip2= 8e-5,Id = 8e-5, fD=0.5, gD=0.5)



DecomOut = ode(y = states, parms = parameters, times = times, func = Decom, method = "radau")
print("done")
plot(DecomOut)
write.csv(DecomOut, "MENDout.csv")
import math 

def SABR(alpha,beta,rho,nu,F,K,time): # all variables are scalars
    
    if K <= 0:   # negative rates' problem, need to shift the smile
        VOL = 0
        
    elif F == K: # ATM formula
        V = (F*K)**((1-beta)/2.)
        logFK = math.log(F/K)
        A = 1 + ( ((1-beta)**2*alpha**2)/(24.*(V**2)) + (alpha*beta*nu*rho)/(4.*V) + ((nu**2)*(2-3*(rho**2))/24.) ) * time
        B = 1 + (1/24.)*(((1-beta)*logFK)**2) + (1/1920.)*(((1-beta)*logFK)**4)
        VOL = (alpha/V)*A
        
    elif F != K: # not-ATM formula
        V = (F*K)**((1-beta)/2.)
        logFK = math.log(F/K)
        z = (nu/alpha)*V*logFK
        x = math.log( ( math.sqrt(1-2*rho*z+z**2) + z - rho ) / (1-rho) )
        A = 1 + ( ((1-beta)**2*alpha**2)/(24.*(V**2)) + (alpha*beta*nu*rho)/(4.*V) + ((nu**2)*(2-3*(rho**2))/24.) ) * time
        B = 1 + (1/24.)*(((1-beta)*logFK)**2) + (1/1920.)*(((1-beta)*logFK)**4)
        VOL = (nu*logFK*A)/(x*B)
    

    return VOL 

def shift(F,K):
    shift = 0.001 - K[0]
    for j in range(len(K)):
        K[j] = K[j] + shift
        F = F + shift   
   


def smile(alpha,beta,rho,nu,F,K,time): # F, time and the parameters are scalars, K and MKT are vectors, i is the index for tenor/expiry label

    VOLS = list()
    K = [F + 0.0001*k for k in K] 
    for j in range(len(K)):
        
        if K[0] <= 0:
            shift(F,K)
        VOLS.append(SABR(alpha,beta,rho,nu,F,K[j],time))

    return VOLS



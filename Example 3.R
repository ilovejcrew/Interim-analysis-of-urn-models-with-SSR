#mua and sigmaa is the parameter values for treatment1
#mub and sigmab is the parameter values for treatment2
#n1 is the number of enrollment at the first stage
#n2 is the number of enrollment at the second stage
#n3_0 is the original planned number of enrollment at the third stage
#m is the replication number
#c1,c2,c3 are O'Brien Fleming boundaries for the interim analysis, respectively
#ssrfl is the indicator for whether we will combine the SSR in the simulation

simulation=function(mua,mub,sigmaa,sigmab,n1,n2,n3_0,m,c1,c2,c3,ssrfl)
{

# Basic setting

ntotal=n1+n2+n3_0

nmax = 500

#nmax = 300


pcut=0.9

t1=0.2  #information time
t2=0.5

failure=NULL          #total failure number
failureratio=NULL     # failure rate
bmaxind=0              #number of the limit of increase of sample size is reached.

########## Adaptive Randomization with DBCD ##########

# values kept after m replications
number=0		      # reject number
numofssr=0              # of ssr was implemented
rho1<-NULL			#allocation proportion of treatment1
rho2<-NULL
urn1<-NULL
urn2<-NULL
reject1=0               # number of reject at first look
reject2=0               # number of reject at second look
reject3=0               # number of reject at third look
SS=NULL                 # final sample size 
SSplus=NULL             #if SSR will be done, the increase of sample size
for  ( i in 1:m ){
        ball1=5         # number of type 1 ball in the urn
        ball2=5         # number of type 2 ball in the urn
        
        xx1<-NULL       # response of patients in treatment 1
        xx2<-NULL       # response of patients in treatment 2

        N1=0            # number of patients in treatment 1
        N2=0            # number of patients in treatment 2

        xx1<-c(xx1,rnorm(10,mua,sigmaa))
        xx2<-c(xx2,rnorm(10,mub,sigmab))
        N1=10
        N2=10

        for (j in 1:(n1-20)){
          	Rho1=ball1/(ball1+ball2) # probability that drawing type 1
	          x<-runif(1,0,1)	
		        if (x>=0 & x<Rho1) {
			                  N1<-N1+1
                        new=rnorm(1,mua,sigmaa) #generate the new response outcome based on the distribution of the response for treatment 1
			                  xx1<-c(xx1,new)	#update the response vector
			                  }
	         	if (x>=Rho1 & x<=1) {
			                  N2<-N2+1
                        new=rnorm(1,mub,sigmab)
			                  xx2<-c(xx2,new)
			                  }
            ball1=ball1+sd(xx1)/(sd(xx1)+sd(xx2)) 	#sd(xx1)/(sd(xx1)+sd(xx2)) balls of type1 are added to the urn
         	  ball2=ball2+sd(xx2)/(sd(xx1)+sd(xx2))      #sd(xx2)/(sd(xx1)+sd(xx2)) balls of type2 are added to the urn
		        }
        
        stat=abs((mean(xx1)-mean(xx2))/sqrt(var(xx1)/N1+var(xx2)/N2)) #calculate the sequential test statistics

        if (stat>c1) {
         number=number+1    #if the test statistics is larger than the critical value, then the reject number increases by 1
         reject1=reject1+1  #if the test statistics is larger than the critical value, then number of reject at first look increases by 1
        } else {           #else the trial goes to the second stage

        for (j in 1:n2){
		            Rho1=ball1/(ball1+ball2) 
	              x<-runif(1,0,1)	
		            if (x>=0 & x<Rho1) {
		                   	N1<-N1+1
                        new=rnorm(1,mua,sigmaa)
			                  xx1<-c(xx1,new)
			                  }
		            if (x>=Rho1 & x<=1) {
			                  N2<-N2+1
                        new=rnorm(1,mub,sigmab)
			                  xx2<-c(xx2,new)
			                  }
                ball1=ball1+sd(xx1)/(sd(xx1)+sd(xx2))
         	      ball2=ball2+sd(xx2)/(sd(xx1)+sd(xx2))
             }
        
        stat2=abs((mean(xx1)-mean(xx2))/sqrt(var(xx1)/N1+var(xx2)/N2)) #calculate the sequential test statistics at the second interim analysis
                                
        if (stat2>c2) {
         number=number+1 #if the test statistics is larger than the critical value, then the reject number increases by 1
         reject2=reject2+1 #if the test statistics is larger than the critical value, then number of reject at second look increases by 1
      
        } else {

          cpfl = FALSE  # indicator of whether SSR will be done based on observed treatment effects
          mu_1 = (mean(xx2)-mean(xx1))/sqrt(var(xx1)/(N1/(N1+N2))+var(xx2)/(N2/(N1+N2)))  #calculate the observed treatment effect
         
          cp_1 =   1-pnorm((c3 - stat2*sqrt(t2)- mu_1*sqrt(ntotal)*sqrt(1-t2))/sqrt(1-t2) ) #calculate the conditional power

          if(0.01<cp_1 & cp_1<pcut) cpfl = TRUE
          
          fx = function(ntotal, pcut0 = pcut ){
                1-pnorm( (c3 - stat2*sqrt(t2)- mu_1*sqrt(ntotal)*sqrt(1-t2))/sqrt(1-t2) )- pcut0           
         }
          
          if(ssrfl & cpfl){
    		    ncp = floor(uniroot(fx,c(n3_0,1000000))$root) - n1 - n2 #if we decide to do the SSR, ncp is the sample size of the third stage to satisfy the power
            n3 = min(nmax, max(n3_0, ncp)) #n3 is the new updated sample size for stage 3
            SSplus=c(SSplus,n3-n3_0)       #if SSR will be done, the increase of sample size
                  } else {
            n3=n3_0
          }  
          if (n3==n3_0) numofssr=numofssr+1   #number of SSR
          if (n3==nmax) bmaxind=bmaxind+1     #number of the limit of increase of sample size is reached.
   
         for (j in 1:n3){
		            Rho1=ball1/(ball1+ball2) 
	              x<-runif(1,0,1)	
	             	if (x>=0 & x<Rho1) {
		                  	N1<-N1+1
                        new=rnorm(1,mua,sigmaa)
			                  xx1<-c(xx1,new)
			                  }
		            if (x>=Rho1 & x<=1) {
			                  N2<-N2+1
                        new=rnorm(1,mub,sigmab)
			                  xx2<-c(xx2,new)
			                  }
                ball1=ball1+sd(xx1)/(sd(xx1)+sd(xx2))
         	      ball2=ball2+sd(xx2)/(sd(xx1)+sd(xx2))
		 
		        }
        
        
        
        b = n3/n3_0
#stat30 is the original test statistics for stage 3, and stat3 is the our new proposed test statistics for stage 3
      	stat30=(mean(xx2)-mean(xx1))/sqrt(var(xx1)/N1+var(xx2)/N2) 
        stat3 = sqrt((n1+n2+n3)/(n1+n2+n3_0))*(mean(xx2)-mean(xx1))/sqrt(var(xx1)/N1+var(xx2)/N2)*sqrt(1-(n1+n2)/(n1+n2+n3_0))/sqrt(b*(n3_0/(n1+n2+n3_0)))-
           stat2*sqrt((n1+n2)/(n1+n2+n3_0))*sqrt(n3_0/(n1+n2+n3_0))/sqrt(b*(n3_0/(n1+n2+n3_0))) +
           stat2*sqrt((n1+n2)/(n1+n2+n3_0))
        
        
        if(ssrfl & cpfl){	
           stat = stat3  #if SSR is done, we will use stat3
        }else{
          stat = stat30 #if SSR is not done, we will use stat30
      }
          


        if (stat>c3) {
        number=number+1
        reject3=reject3+1
       
        } 
       }
       }
       rho1=c(rho1,N1/(N1+N2)) #allocation proportion of treatment1
       rho2=c(rho2,N2/(N1+N2)) #allocation proportion of treatment2
       urn1=c(urn1,ball1/(ball1+ball2)) #propotion of balls of type1 in the urn
       urn2=c(urn2,ball2/(ball1+ball2)) #propotion of balls of type2 in the urn
       SS=c(SS,N1+N2) #total sample size
       }
result=c(number/m,mean(rho1),sd(rho1),mean(urn1),sd(urn1),mean(SS),sd(SS),reject1+reject2,bmaxind)

return(result)

rm(list = ls(all = TRUE))
}


simulation=function(pa,pb,n1,n2,n3_0,m,c1,c2,c3,ssrfl)
{

# Basic setting

ntotal=n1+n2+n3_0

nmax = 500

#nmax = 300

pcut=0.9

t1=0.2
t2=0.5

failure=NULL
failureratio=NULL
bmaxind=0

########## Adaptive Randomization with DBCD ##########

# values kept after m replications
number=0		# reject number
numofssr=0              # of ssr was implemented
rho1<-NULL
rho2<-NULL
urn1<-NULL
urn2<-NULL
reject1=0 # number of reject at first look
reject2=0
reject3=0
SS=NULL                # final sample size 
SSplus=NULL            #if SSR will be done, the increase of sample size
for  ( i in 1:m ){
        ball1=5  # number of type 1 ball in the urn
        ball2=5
        
        xx1<-NULL  # response of patients go to treatment 1
        xx2<-NULL

        N1=0  # number of patients go to treatment1
        N2=0

        for (j in 1:n1){
          	Rho1=ball1/(ball1+ball2) # probability that drawing type 1
	
                x<-runif(1,0,1)	
		if (x>=0 & x<Rho1) {
			N1<-N1+1
                        new=rbinom(1,1,pa)
			xx1<-c(xx1,new)
                                                
		}
		if (x>=Rho1 & x<=1) {
			N2<-N2+1
                        new=rbinom(1,1,pb)
			xx2<-c(xx2,new)
                       
		}
                p1hat=(sum(xx1)+1)/(N1+1)
                p2hat=(sum(xx2)+1)/(N2+1)
                ball1=ball1+sqrt(p1hat)
         	ball2=ball2+sqrt(p2hat)
		        }
        
        p1hat=(sum(xx1)+1)/(N1+1)
        p2hat=(sum(xx2)+1)/(N2+1)
        
        stat=abs((p1hat-p2hat)/sqrt(p1hat*(1-p1hat)/N1+p2hat*(1-p2hat)/N2))	

        if (stat>c1) {
        number=number+1
        reject1=reject1+1
        } else {

        for (j in 1:n2){
		Rho1=ball1/(ball1+ball2) # probability that drawing type 1
	
                x<-runif(1,0,1)	
		if (x>=0 & x<Rho1) {
			N1<-N1+1
                        new=rbinom(1,1,pa)
			xx1<-c(xx1,new)
                        
		}
		if (x>=Rho1 & x<=1) {
			N2<-N2+1
                        new=rbinom(1,1,pb)
			xx2<-c(xx2,new)
                       
		}
                p1hat=(sum(xx1)+1)/(N1+1)
                p2hat=(sum(xx2)+1)/(N2+1)
		ball1=ball1+sqrt(p1hat)
         	ball2=ball2+sqrt(p2hat)
		        }
        
        p1hat=(sum(xx1)+1)/(N1+1)
        p2hat=(sum(xx2)+1)/(N2+1)
       
        stat2=abs((p1hat-p2hat)/sqrt(p1hat*(1-p1hat)/N1+p2hat*(1-p2hat)/N2))	
                               
        if (stat2>c2) {
        number=number+1
         reject2=reject2+1
       
        } else {

          cpfl = FALSE
          mu_1 = (mean(xx2)-mean(xx1))/sqrt(mean(xx1)*(1-mean(xx1))/(N1/(N1+N2))+mean(xx2)*(1-mean(xx2))/(N2/(N1+N2)))
         
         cp_1 =   1-pnorm((c3 - stat2*sqrt(t2)- mu_1*sqrt(ntotal)*sqrt(1-t2))/sqrt(1-t2) )

          if(0.01<cp_1 & cp_1<pcut) cpfl = TRUE
          
          fx = function(ntotal, pcut0 = pcut ){
                1-pnorm( (c3 - stat2*sqrt(t2)- mu_1*sqrt(ntotal)*sqrt(1-t2))/sqrt(1-t2) )- pcut0           
         }
          
          if(ssrfl & cpfl){
            ncp = floor(uniroot(fx,c(n3_0,1000000))$root) - n1 - n2
            n3 = min(nmax, max(n3_0, ncp))
            SSplus=c(SSplus,n3-n3_0)                   
          } else {
            n3=n3_0
          }  
          if (n3==n3_0) numofssr=numofssr+1
          if (n3==nmax) bmaxind=bmaxind+1
  
         for (j in 1:n3){
		Rho1=ball1/(ball1+ball2) # probability that drawing type 1
	
                x<-runif(1,0,1)	
		if (x>=0 & x<Rho1) {
			N1<-N1+1
                        new=rbinom(1,1,pa)
			xx1<-c(xx1,new)
                       
		}
		if (x>=Rho1 & x<=1) {
			N2<-N2+1
                        new=rbinom(1,1,pb)
			xx2<-c(xx2,new)
                        
		}
                p1hat=(sum(xx1)+1)/(N1+1)
                p2hat=(sum(xx2)+1)/(N2+1)
                ball1=ball1+sqrt(p1hat)
         	ball2=ball2+sqrt(p2hat)
		
		        }
        
        
        
 b = n3/n3_0
      	stat30=(mean(xx2)-mean(xx1))/sqrt(mean(xx1)*(1-mean(xx1))/N1+mean(xx2)*(1-mean(xx2))/N2) 
        stat3 = sqrt((n1+n2+n3)/(n1+n2+n3_0))*(mean(xx2)-mean(xx1))/sqrt(mean(xx1)*(1-mean(xx1))/N1+mean(xx2)*(1-mean(xx2))/N2)*sqrt(1-(n1+n2)/(n1+n2+n3_0))/sqrt(b*(n3_0/(n1+n2+n3_0)))-
           stat2*sqrt((n1+n2)/(n1+n2+n3_0))*sqrt(n3_0/(n1+n2+n3_0))/sqrt(b*(n3_0/(n1+n2+n3_0))) +
           stat2*sqrt((n1+n2)/(n1+n2+n3_0))
        
        
        if(ssrfl & cpfl){	
          stat = stat3
        }else{
          stat = stat30
        }
          


        if (stat>c3) {
        number=number+1
        reject3=reject3+1
       
        } 
       }
       }
       rho1=c(rho1,N1/(N1+N2))
       rho2=c(rho2,N2/(N1+N2))
       urn1=c(urn1,ball1/(ball1+ball2))
       urn2=c(urn2,ball2/(ball1+ball2))
       SS=c(SS,N1+N2)
       failure=c(failure,length(xx1)+length(xx2)-sum(xx1)-sum(xx2))
       failureratio=c(failureratio,(length(xx1)+length(xx2)-sum(xx1)-sum(xx2))/(length(xx1)+length(xx2)))
}
result=c(number/m,mean(rho1),sd(rho1),mean(urn1),sd(urn1),mean(SS),sd(SS),mean(failure),sd(failure),mean(failureratio),sd(failureratio),reject1,reject2,reject3,bmaxind)

return(result)

rm(list = ls(all = TRUE))
}
rm(list = ls())

library("bife")
library("reshape")
library("car")
library("grid")
library("gridExtra")

N <- 1000
Time <- c(4,8,10,12,16,20,50,100,200)
N_star <- c(100,1000)


#true value
beta  <- 1 
delta <- 1

set.seed(666)

for(i in 1:length(N_star)){  
  
  N0 <- N_star[i]
  UCL.param.mean <- BCL.param.mean <- CL.param.mean <- matrix(NA,length(Time),2)
  UCL.param.sd <- BCL.param.sd <- CL.param.sd <- matrix(NA,length(Time),2)
  UCL.param.se <- BCL.param.se <- CL.param.se <- matrix(NA,length(Time),2)
  table_2A <- table_2B <- c()
  
  
  pdf("Figure 2.pdf", width = 12, height = 12, onefile = TRUE)
  par(oma=c(2,2,2,2),mfrow=c(3,3))

  for(j in 1:length(Time)){
 
    T0 <- Time[j]
    BCL.param <- UCL.param <- CL.param <- matrix(NA,N,2)
    
    for(k in 1:N){
      
      #DGP
      x <- matrix(rnorm(N0*T0),N0,T0) 
      h <- matrix(rnorm(N0*T0),N0,T0) 
      u <- matrix(rnorm(N0*T0),N0,T0)
      a <- matrix(rnorm(N0),N0,1)
      y <- matrix(0,N0,T0)
      v <- matrix(0,N0,T0)
      d <- matrix(0,N0,T0)
      w <- matrix(0,N0,T0)
      alpha <- matrix(0,N0,1)
      
      for(i in 1:N0){
        for(t in 1:T0){
          
          v[i,t] <- log(u[i,t]/(1-u[i,t]))
          
          if(x[i,t] + h[i,t] > 0){
            d[i,t] <- 1
          }
        }
        alpha[i] <- sqrt(T0)*rowMeans(x)[i] + a[i]
      }
      
      for(i in 1:N0){
        for(t in 1:T0){
          w[i,t] <- alpha[i] + beta*x[i,t] + delta*d[i,t]
        }
      }  
      
      for(i in 1:N0){
        for(t in 1:T0){
          
          if(w[i,t] + v[i,t] > 0){
            y[i,t] <- 1
          }
          
        }
      }
      
      x1 <- as.data.frame(cbind(c(1:nrow(x)),x)) 
      x1.long <- reshape(x1, varying = list(names(x1)[2:ncol(x1)]),
                         idvar = "V1", timevar = "time", times = c(1:T0), direction = "long")
      d1 <- as.data.frame(cbind(c(1:nrow(d)),d)) 
      d1.long <- reshape(d1, varying = list(names(d1)[2:ncol(d1)]),
                         idvar = "V1", timevar = "time", times = c(1:T0), direction = "long")
      y1 <- as.data.frame(cbind(c(1:nrow(y)),y)) 
      y1.long <- reshape(y1, varying = list(names(y1)[2:ncol(y1)]),
                         idvar = "V1", timevar = "time", times = c(1:T0), direction = "long")
      
      alpha_panel <- c()
      for(t in 1:T0){
        alpha_panel <- append(alpha_panel,alpha)
      }
      alpha_panel <- as.data.frame(alpha_panel)
      
      data0 <- cbind(x1.long,y1.long[,3],d1.long[,3],alpha_panel)
      colnames(data0) <- c("id","time","X","Y","D","alpha")
      row.names(data0) <-c()
      remove("x1","x1.long","y1","y1.long","d1","d1.long","alpha_panel")
    
      #Estimation
      CL <- glm(Y ~ X + D + alpha, data = data0, family = binomial())
      CL.param[k,] <- CL$coefficients[-c(1,4)]
      
      BCL <- bife(Y ~ X + D |id, data = data0, bias_corr = "ana")
      BCL.param[k,] <- coef(BCL)
      
      UCL <- bife(Y ~ X + D |id, data = data0, bias_corr = "no")
      UCL.param[k,] <- coef(UCL)
    }
    
    #mean
    BCL.param.mean[j,] <- colMeans(BCL.param)
    UCL.param.mean[j,] <- colMeans(UCL.param)
    CL.param.mean[j,] <- colMeans(CL.param)
    
    #sd
    BCL.param.sd[j,] <- apply(BCL.param,2,sd)
    UCL.param.sd[j,] <- apply(UCL.param,2,sd)
    CL.param.sd[j,] <- apply(CL.param,2,sd)
    
    #se
    se <- function(x){
      sd(x)/sqrt(length(x))
    }
    
    BCL.param.se[j,] <- apply(BCL.param,2,se)
    UCL.param.se[j,] <- apply(UCL.param,2,se)
    CL.param.se[j,] <- apply(CL.param,2,se)
    
    #Figure 2
    heading = paste("T=",T0) 
    d1 <- density(UCL.param[,1])
    d2 <- density(BCL.param[,1])
    d3 <- density(CL.param[,1])
    plot(d3,col = "blue", xlab="Beta_hat",main=heading)
    lines(d2,col = "orange")
    lines(d1,col = "red")
    
    #table 2
    UVB_beta <- ks.test(UCL.param[,1],BCL.param[,1])
    UVC_beta <- ks.test(UCL.param[,1],CL.param[,1])
    BVC_beta <- ks.test(BCL.param[,1],CL.param[,1])
    KS.Test_beta <- c(UVB_beta$p.value, UVC_beta$p.value, BVC_beta$p.value)
    
    UVB_delta <- ks.test(UCL.param[,2],BCL.param[,2])
    UVC_delta <- ks.test(UCL.param[,2],CL.param[,2])
    BVC_delta <- ks.test(BCL.param[,2],CL.param[,2])
    KS.Test_delta <- c(UVB_delta$p.value, UVC_delta$p.value, BVC_delta$p.value)
    
    temp1 <- cbind(T0,t(KS.Test_beta))
    temp2 <- cbind(T0,t(KS.Test_delta))
    table_2A <- rbind(table_2A, temp1)
    table_2B <- rbind(table_2B, temp2)
    colnames(table_2A) <- c("(beta)","UCL vs. BCL","UCL vs. CL","BCL vs. CL")
    colnames(table_2B) <- c("(delta)","UCL vs. BCL","UCL vs. CL","BCL vs. CL")
    remove(temp1,temp2)
  
    cat("Currently completed T=: ", T0,"\n")
  }
  dev.off()  
  
  table_1 <- cbind(Time,UCL.param.mean,BCL.param.mean,CL.param.mean)[,c(1,2,4,6,3,5,7)]
  colnames(table_1) <- c("T","UCL.beta","BCL.beta","CL.beta","UCL.delta","BCL.delta","CL.delta")
  
  #Figure 1
  pdf(sprintf("Figure 1_(N_star_%s).pdf", N0), width = 12, height = 12, onefile = TRUE)
  attach(as.data.frame(table_1))
  UCL.plot <- plot(T,UCL.beta, xlab="T",ylab="Beta_hat", type = "l", col = "red", ylim=c(0.6,2))
  lines(T, BCL.beta, col="orange")
  lines(T, CL.beta, col="blue")
  detach(as.data.frame(table_1))
  dev.off()
  
  pdf(sprintf("Table 1_(N_star_%s).pdf", N0), width = 12, height = 12, onefile = TRUE)
  grid.table(format(round(table_1,digits = 4), nsmall = 4))
  dev.off() 
  
  pdf(sprintf("Table 2_(N_star_%s).pdf", N0), width = 12, height = 12, onefile = TRUE)
  par(oma=c(2,2,2,2),mfrow=c(2,1))
  grid.table(format(round(table_2A,digits = 4), nsmall = 4))
  grid.newpage()
  grid.table(format(round(table_2B,digits = 4), nsmall = 4))
  dev.off() 
  
  cat("Completed N*=: ", N0,"\n")
  
}  

 


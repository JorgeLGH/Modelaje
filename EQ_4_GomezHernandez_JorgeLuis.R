library(deSolve)
##
##
SIR <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), { 
    N<- S+I+R
    dS<- nu*N-beta*S*I/(S+I+R)-mu*S
    dI<- beta*S*I/(S+I+R) -gama*I - mu*I
    dR<- gama*I - mu*R
    list(c(dS, dI, dR))
  })
}
#beta<gama
pars <- c(beta= 2, gama = 4, mu=0.1, nu=0.1 )  
condiciones_iniciales <- c(S =999, I=1,R=0)
tiempo <- seq(0, 20, by = 0.001)
out <- ode(condiciones_iniciales, tiempo, SIR, pars) 
head(out)
matplot(out[ , 1], out[ , 2:4], type = "l", xlab = "tiempo", ylab = "Población",
        main = "SIS", lwd = 2)
legend("topright", c("Susceptible", "Infectado", "Recuperado"), col = 1:3,lty=1:3,cex=0.5)
#beta>gama
pars <- c(beta= 4, gama = 2, mu=0.1, nu=0.1 )
condiciones_iniciales <- c(S =999, I=1,R=0)
tiempo <- seq(0, 20, by = 0.001)
out <- ode(condiciones_iniciales, tiempo, SIR, pars) 
head(out)
matplot(out[ , 1], out[ , 2:4], type = "l", xlab = "tiempo", ylab = "Población",
        main = "SIS", lwd = 2)
legend("topright", c("Susceptible", "Infectado", "Recuperado"), col = 1:3,lty=1:3,cex=0.5)
#beta=gama
pars <- c(beta= 4, gama = 4, mu=0.1, nu=0.1 )
condiciones_iniciales <- c(S =999, I=1,R=0)
tiempo <- seq(0, 20, by = 0.001)
out <- ode(condiciones_iniciales, tiempo, SIR, pars) 
head(out)
matplot(out[ , 1], out[ , 2:4], type = "l", xlab = "tiempo", ylab = "Población",
        main = "SIS", lwd = 2)
legend("topright", c("Susceptible", "Infectado", "Recuperado"), col = 1:3,lty=1:3,cex=0.5)



####
####
SEIR <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), { 
    N<- S+E+I+R
    dS<- nu*N-beta*S*I/(S+I+R)-mu*S
    dE<- beta*S*I/(S+I+R) -mu*E -delta*E
    dI<- delta*E -gama*I - mu*I
    dR<- gama*I - mu*R
    list(c(dS, dI, dR,dE))
  })
}
#beta>gama
pars <- c(beta= 4, gama = 2, mu=0.1, nu=0.1, delta=2 )
condiciones_iniciales <- c(S =999, I=1,R=0,E=0)
tiempo <- seq(0, 20, by = 0.001)
out <- ode(condiciones_iniciales, tiempo, SEIR, pars) 
head(out)
matplot(out[ , 1], out[ , 2:5], type = "l", xlab = "tiempo", ylab = "Población",
        main = "SEIR", lwd = 2)
legend("topright", c("Susceptible", "Infectado", "Expuesto", "Recuperado"), col = 1:5,lty=1:3,cex=0.5)
#beta<gama
pars <- c(beta= 2, gama = 4, mu=0.1, nu=0.1, delta=2 )
condiciones_iniciales <- c(S =999, I=1,R=0,E=0)
tiempo <- seq(0, 20, by = 0.001)
out <- ode(condiciones_iniciales, tiempo, SEIR, pars) 
head(out)
matplot(out[ , 1], out[ , 2:5], type = "l", xlab = "tiempo", ylab = "Población",
        main = "SEIR", lwd = 2)
legend("topright", c("Susceptible", "Infectado", "Expuesto", "Recuperado"), col = 1:5,lty=1:3,cex=0.5)
#beta=gama
pars <- c(beta= 4, gama = 4, mu=0.1, nu=0.1, delta=2 )
condiciones_iniciales <- c(S =999, I=1,R=0,E=0)
tiempo <- seq(0, 20, by = 0.001)
out <- ode(condiciones_iniciales, tiempo, SEIR, pars) 
head(out)
matplot(out[ , 1], out[ , 2:5], type = "l", xlab = "tiempo", ylab = "Población",
        main = "SEIR", lwd = 2)
legend("topright", c("Susceptible", "Infectado", "Expuesto", "Recuperado"), col = 1:5,lty=1:3,cex=0.5)



#Hice la evaluación de todo con distintos valores de gama y beta porque afectan de manera directa el valor del R_0, de modo
#que ya está implícito que se cambia el R_0 dependiendo de la intensidad de estos cambios, pero por fines prácticos solo evalué
#cuando R_0 es menor, mayor o igual a 1
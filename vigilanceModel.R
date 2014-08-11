library(deSolve); library(ggplot2); library(reshape); library(gridExtra); 
library(reshape2)

#---------------------------------------#
#         Model Specification           #
#---------------------------------------#

# Model parameters and values
parameters <- c(r = 0.2,                     #intrinsic rate of growth
                v = 0.105,                   #vigilance
                dN = 0.1,                    #prey death rate
                K = 100,                     #carrying capacity
                c = 0.1,                     #conversion rate
                dP = 0.1,                    #predator death rate
                m = 0.1,                     #encounter rate
                k = 1,                       #inverse of predator lethality
                b = 5,                       #value of vigilance
                z = expression(m/(k+(b*v))), #predation risk
                r0 = expression((1-v)*r))    #prey birth rate 


# Initial conditions
max <- with(as.list(parameters), return(K*(1-(dN/r)))) #maximum starting N & isocline range
state <- c(N = max, P = 5)                             #starting population sizes
times <- expression(seq(0, 1000, by = 1))              #number of iterations to run

# Isocline range
isoRange <- expression(seq(0, max, 1))

# Run the models!
fixedVigilance(state, times, parameters, "fixedDynamics", isoRange)
flexibleVigilance(state, times, parameters, "flexibleDynamics", isoRange)

# Find the equilibrium abundance of P in the fixed vigilance model (w/ vigilance)
predEquilibrium(parameters, vigilance = TRUE)

# plot Sweeps
sweep <- sweepParameters(state, times, parameters)
plotSweepNP(sweep)
plotSweepV(sweep)

# G-functions
Gvalues <- getGP(state, times, parameters, "fixed")
Gvalues <- getGP(state, times, parameters, "flexible")
contourPlot(Gvalues, "k", "b") #choose m+b, k+b, or k+m

#---------------------------------------#
#              Functions                #
#---------------------------------------#

#---DYNAMICS---#

# with fixed vigilance, calculates change in N and P over a time step
fixedDynamics <- function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
        z <- eval(z)                           #evaluate z
        r0 <- eval(r0)                         #evaluate r0
        dN.dt <- r0*N*((K-N)/K)-(dN*N)-(z*N*P) #determine change in N
        dP.dt <- ((z*c*N)-dP)*P                #determine change in P
        list(c(dN.dt, dP.dt))      
    })
}

# with flexible vigilance, calculates change in N and P over a time step
flexibleDynamics <- function(t, state, parameters) {
    with(as.list(c(state, parameters)),{   
        v <- (((m*K*P)/(r*b*(K-N)))^0.5)-(k/b) #determine vigilance, overwrite value in parameters
        if (v < 0 | v > 1 | is.na(v)) v <- 0   #if v is less than 0, or greater than 1, set it to 0
        z <- eval(z)                           #evaluate z
        r0 <- eval(r0)                         #evaluate r0
        dN.dt <- r0*N*((K-N)/K)-(dN*N)-(z*N*P) #determine change in N
        dP.dt <- ((z*c*N)-dP)*P                #determine change in P
        list(c(dN.dt, dP.dt))                  
    })
} 

# returns data frame of N and P for a given time sequence (then tacks on v)
getDynamics <- function(state, times, parameters, modelType) {
    with(as.list(parameters),{
        dynamics <- ode(y = state, times = eval(times), func = match.fun(modelType), parms = parameters)
        if (modelType == "flexibleDynamics") {
            v <- (((m*K*dynamics[,"P"])/(r*b*(K-dynamics[,"N"])))^0.5)-(k/b)
            v[v < 0 | is.na(v)] <- 0
        }
        return(data.frame(dynamics, v)) #this just tacks on v to the output, can't find a better way to do this
    })
}

# plots dynamics of N and P over time
plotDynamics <- function(dynamics) {
    return(ggplot(melt(as.data.frame(dynamics[,c(1:3)]), id="time"), 
                  aes(time, value, group = variable, color = variable)) + 
               geom_line(size = 1.5) + labs(y = "population") + theme_bw() + 
               theme(legend.direction = "horizontal", legend.position=c(0.5, 0.9)))
}

# Find the equilibrium predator abundance in the fixed model with/without vigilance
predEquilibrium <- function(parameters, vigilance) {
    with(as.list(parameters),{
        ifelse(vigilance, {
            z <- eval(z)
            r0 <- eval(r0)
            P <- (c*K*z*(r0-dN)-(r0*dP))/(c*K*z^2)
        }, P <- ((k*(r-dN))/m)-((r*(k^2)*dP)/c*K*m^2))
        return(P)
    })
}

#---ISOCLINES---#

# return data frame of N and P isocline values (fixed vigilance)
fixedIsoclines <- function(isoRange, parameters) {
    with(as.list(parameters),{
        z <- eval(z)                                     #evaluate z
        r0 <- eval(r0)                                   #evaluate r0
        N <- eval(isoRange)                              #find the range of N values
        P <- (r0-dN)/z - (r0*N)/(z*K)                    #find the prey isocline value (P) for each N
        Niso <- data.frame(N = N, P = P, isocline = "N") #write N and P to a dataframe    
        P <- 0.25*N                                      #set the range of P values to 25% N
        N <- dP/(c*z)                                    #find the predator isocline value (N) for each P
        Piso <- data.frame(N = N, P = P, isocline = "P") #write N and P to a dataframe    
        return(rbind(Niso, Piso))
    })
}

# return data frame of N and P isocline values (flexible vigilance)
flexibleIsoclines <- function(isoRange, parameters) {
    with(as.list(parameters),{
        N <- eval(isoRange)                                       #evaluate the range of N to test
        P <- b*K*((r*(1+k/b)*((K-N)/K)-dN)^2)/(4*r*m*(K-N))       #evaluate prey isocline (P) for N range (Abdel's version)
        #P <- ((K*b*dN-K*b*r-K*k*r+N*b*r+N*k*r)^2)/(4*b*r*K*m*(K-N)) #Paul's version (same result as Abdel's)
        v <- (((m*K*P)/(r*b*(K-N)))^0.5)-(k/b)                    #evaluate v for those N and P values
        Niso <- data.frame(N = N, P = P, v = v, isocline = "N")   #write N, P, v to a data frame
        Niso[Niso$v < 0, "P"] <- (k*(r-dN))/m - ((r*k*Niso[Niso$v < 0, "N"])/(m*K)) #if v negative, recalculate isocline (P)
        P <- (((c*N)/dP)^2)*(r*m*(K-N)/(b*K))                     #evaluate predator isocline (P) for N range 
        v <- (((m*K*P)/(r*b*(K-N)))^0.5)-(k/b)                    #evaluate v for the new N and P values
        Piso <- data.frame(N = N, P = P, v = v, isocline = "P")   #write N, P, v to a data frame
        Piso[Piso$v < 0, "N"] <- (k*dP)/(c*m)                     #if v negative, recalculate isocline (N)
        P <- (r*k^2)/(b*m) - (r*k^2*N)/(b*m*K)                    #plot the isoleg where v is 0
        isoLeg <- data.frame(N = N, P = P, v = 0, isocline = "leg")
        return(rbind(Niso, Piso, isoLeg))
    })
}

# plot predator and prey isoclines
plotIsoclines <- function(isoRange, parameters, isoType) {
    return(ggplot(match.fun(isoType)(isoRange, parameters), 
                  aes(N, P, group = isocline, color = isocline)) +
               geom_line(size = 1.5)  + theme_bw() + 
               theme(legend.direction = "horizontal", legend.position=c(0.5, 0.9)))
}

#---CONSOLIDATION---#

# calculate and print output for the fixed vigilance model
fixedVigilance <- function(state, times, parameters, modelType, isoRange) {
    dynamics <- getDynamics(state, times, parameters, modelType)
    cat("Model: Fixed vigilance\nParameters:\n")
    print(unlist(as.list(parameters[1:9])))
    cat("Ending values:\n")
    print(round(dynamics[nrow(dynamics),], 3))
    return(grid.arrange(plotDynamics(dynamics), 
                        plotIsoclines(isoRange, parameters, fixedIsoclines),
                        ncol=2))
}

# calculate and print output for the flexible vigilance model
flexibleVigilance <- function(state, times, parameters, modelType, isoRange) {
    dynamics <- getDynamics(state, times, parameters, modelType)
    cat("Model: Flexible vigilance\nParameters:\n")
    print(unlist(as.list(parameters[1:9])))
    cat("Ending values:\n")
    print(round(dynamics[nrow(dynamics),], 3))
    return(grid.arrange(plotDynamics(dynamics), 
                        plotIsoclines(isoRange, parameters, flexibleIsoclines),
                        ggplot(dynamics, aes(N, v)) + geom_point(size = 1.5) + theme_bw(),
                        ggplot(dynamics, aes(P, v)) + geom_point(size = 1.5) + theme_bw()))
}

#---PARAMETER SWEEPS---#

# return the appropriate testing range for a given parameter
getParameterRange <- function(parameter) {
    ifelse(parameter == "b" | parameter == "k", range <- seq(0, 10, 1), {
        ifelse(parameter == "m", range <- seq(0, 0.45, 0.05), range <- seq(0, 0.5, 0.05))
    })
    return(range)
}

# return the equilibrium N, P, and v values for a given parameter range
getParameterValues <- function(parameter, times, parameters, last) {
    max <- with(as.list(parameters), return(K*(1-(dN/r))))
    endingValues <- data.frame()
    model <- c(rep("fixed", 2), rep("flexible", 3))
    group <- c(rep(c("N", "P"), 2), "v")
    range <- getParameterRange(parameter)
    for (value in range) {
        parameters[parameter] <- value
        if (parameter == "dN" | parameter == "r") {
            max <- with(as.list(parameters), return(K*(1-(dN/r))))
            if (is.na(max) | max < 0 | is.infinite(max)) next
        }
        state <- c(N = max, P = 5) 
        fixed <- getDynamics(state, times, parameters, fixedDynamics)[last,]
        flexible <- getDynamics(state, times, parameters, flexibleDynamics)[last,]
        values <- c(fixed[, 2], fixed[, 3], flexible[, 2], flexible[, 3], flexible[, 4])    
        result <- data.frame(parameter = parameter, parValue = value,
                             model = model, group = group, value = values)
        endingValues <- rbind(endingValues, result)
    }
    return(endingValues)
}

# return a data frame of equilibrium N, P, and v values for all parameter ranges
sweepParameters <- function(state, times, parameters) {
    last <- length(eval(times))
    parameters["v"] <- getDynamics(state, times, parameters, flexibleDynamics)[last, 4] #overwrite parameter value with v
    evalParameters <- names(unlist(as.list(parameters[c(1:3, 5:9)]))) #create a list of eight parameters to sweep
    sweep <- lapply(evalParameters, getParameterValues, times = times, parameters = parameters, last = last)
    return(do.call("rbind", sweep))
}

# plot the N and P dynamics for the parameter sweep
plotSweepNP <- function(sweep) {
    return(ggplot(sweep[sweep$group != "v",], aes(parValue, value, color = group)) + 
               geom_line(aes(linetype = model), size = 1)  +
               facet_wrap(~parameter+group, scale = "free", ncol = 4) + 
               labs(x = "parameter value", y = "population") + theme_bw() + 
               theme(strip.text = element_text(size = 14, face = "bold")))
}

# plot the v dynamics for the parameter sweep
plotSweepV <- function(sweep) {
    return(ggplot(sweep[sweep$group == "v",], aes(parValue, value)) + 
               geom_line(size = 1)  +
               facet_wrap(~parameter, scale = "free_x", ncol = 4) + 
               labs(x = "parameter value", y = "v*") + theme_bw() + 
               theme(strip.text = element_text(size = 14, face = "bold")))
}

#---G-FUNCTIONS---#

# return data frame of G-function values for range in m, k, and b
getGP <- function(state, times, parameters, modelType){
    last <- length(eval(times))   
    mRange <- seq(0, 0.50, 0.05)
    kRange <- bRange <- range <- seq(1, 10, 1)    
    Gvalues <- data.frame()   
    ifelse(modelType == "fixed", 
           vals <- getDynamics(state, times, parameters, fixedDynamics)[last,],
           vals <- getDynamics(state, times, parameters, flexibleDynamics)[last,])
    N <- vals[,2]
    P <- vals[,3]
    parameters["v"] <- vals[,4]
    with(as.list(parameters),{
        for (m in mRange) {
            for (k in kRange) {
                for (b in bRange) {
                    z <- m/(k+b*v)
                    ifelse(modelType == "fixed", 
                           G <- z*c*N - dP,
                           G <- c*N*(((r*m*(K-N))/(b*K*P))^0.5) - dP)
                    result <- data.frame(m=m, k=k, b=b, G=G)
                    Gvalues <- rbind(Gvalues, result)
                }
            }
        }
        return(Gvalues)
    }) 
}

# create a filled contour plot of G-function values for pairing of m, k, and b
contourPlot <- function(Gvalues, xpar, ypar) {
    xlim <- ylim <- 1:10
    par <- c(xpar, ypar)
    ifelse(xpar == "m", {xlim <- unique(Gvalues$m); xlab <- "m"}, 
           ifelse(xpar == "k", xlab <- "k", xlab <- "b"))
    ifelse(ypar == "m", {ylim <- unique(Gvalues$m); ylab <- "m"}, 
           ifelse(ypar == "k", ylab <- "k", ylab <- "b"))
    ifelse(!("m" %in% par), exclude <- "m", 
           ifelse(!("k" %in% par), exclude <- "k", exclude <- "b"))
    ifelse(exclude == "k",
           Y <- acast(Gvalues[Gvalues$k == 1, -2], m ~ b),
           ifelse(exclude == "m",
                  Y <- acast(Gvalues[Gvalues$m == 0.45, -1], k ~ b),
                  Y <- acast(Gvalues[Gvalues$b == 1, -3], k ~ m)))
    return(filled.contour(xlim, ylim, Y, color.palette = heat.colors, nlevels = 3,
                          xlab = xlab, ylab = ylab, 
                          main = paste("predator G-function with change in", xpar, "and", ypar),  
                          plot.axes={ axis(1, xlim); axis(2, ylim); 
                                      contour(x = xlim, y = ylim, z = Y, 
                                              drawlabels = FALSE, add = T,
                                              nlevels = 3) }))
}
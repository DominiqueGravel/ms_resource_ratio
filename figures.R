##################################################################
##################################################################
# Preliminary functions
##################################################################
##################################################################
library(rootSolve)

# The model
model = function(t,y,pars) {
	with(as.list(c(y,pars)), {

		Gx = min(a1*N1x/q1,a2*N2x/q2)
		Gy = min(a1*N1y/q1,a2*N2y/q2)	

		dN1x 	= e*S1x - e*N1x - q1*Gx*Px + q1*phi*d*Dx - DN1*(N1x-N1y)
		dN2x 	= e*S2x - e*N2x - q2*Gx*Px + q2*phi*d*Dx - DN2*(N2x-N2y)
		dPx 	= Gx*Px - m*Px - DP*(Px-Py)
		dDx 	= m*Px - d*Dx - DD*(Dx-Dy)
		
		dN1y 	= e*S1y - e*N1y - q1*Gy*Py + q1*phi*d*Dy - DN1*(N1y-N1x)
		dN2y 	= e*S2y - e*N2y - q2*Gy*Py + q2*phi*d*Dy - DN2*(N2y-N2x)
		dPy 	= Gy*Py - m*Py - DP*(Py-Px)
		dDy 	= m*Py - d*Dy - DD*(Dy-Dx)
			
		list(c(dN1x,dN2x,dN1y,dN2y,dPx,dPy,dDx,dDy))
		})
	}

# Function calculating equilibrium densities
eq_fn = function(S1x,S2x,S1y,S2y,DN1,DN2,DP,DD,q1) {
	q2=1-q1
	m 	= 0.1
	a1 	= m*q1/N1s
	a2 	= m*q2/N2s

	pars = c(e=e,a1=a1,a2=a2,q1=q1,q2=1-	q1,m=m,S1x=S1x,S2x=S2x,S1y=S1y,S2y=S2y,phi=phi,d=d,DN1=DN1,DN2=DN2,DP=DP,DD=DD)	
	y = c(N1x=S1x,N2x=S2x,N1y=S1y,N2y=S2y,Px=0.01,Py = 0.01,Dx=0,Dy=0)
	eq = as.list(runsteady(y=y, func=model, parms=pars)[[1]])

	Gx = min(a1*eq$N1x/q1,a2*eq$N2x/q2)
	Gy = min(a1*eq$N1y/q1,a2*eq$N2y/q2)	

	S1x_net = q1*Gx*eq$Px/e + eq$N1x
	S2x_net = q2*Gx*eq$Px/e + eq$N2x
	S1y_net = q1*Gy*eq$Py/e + eq$N1y
	S2y_net = q2*Gy*eq$Py/e + eq$N2y

	return(list(eq,c(S1x_net,S2x_net,S1y_net,S2y_net)))
	}

##################################################################
##################################################################
# Figure 1: Graphical representation of the resource-ratio theory
##################################################################
##################################################################

quartz(height = 10, width = 5)
par(mar=c(6,6,2,1), mfrow = c(2,1))

# Parameters
N1s = 5
N2s = 5
e 	= 0.1
phi = 0.3
d 	= 0.1
S1x = 50
S2x = 30
S1y = 0
S2y = 0
DN1 = 0
DN2 = 0
DP = 0
DD = 0
q1 = 0.5

# Equibrium densities
eq = eq_fn(S1x,S2x,S1y,S2y,DN1,DN2,DP,DD,q1)
eqp = eq_fn(S1x,S2x,S1y,S2y,DN1=10,DN2=10,DP,DD,q1)

eqN1 = eq[[1]]$N1x
S1xnet = eq[[2]][1]
c1 = S1xnet - eqN1
v1 = eqN1 - c1*.15

eqN2 = eq[[1]]$N2x
S2xnet = eq[[2]][2]
c2 = S2xnet - eqN2
v2 = eqN2 - c2*.15

eqN1p = eqp[[1]]$N1x
S1p = eqp[[2]][1]
v1p = eqN1p - c1*.15

eqN2p = eqp[[1]]$N2x
S2p = eqp[[2]][2]
v2p = eqN2p - c2*.15

##################################################################
# Panel a) Extension to resource-ratio
##################################################################
# Supply points
plot(S1x,S2x,xlim=c(0,50),ylim=c(0,50), xlab = "Nutrient 1", ylab = "Nutrient 2",pch=19,cex = 1.5, cex.axis = 1.5, cex.lab = 1.75)
points(S1p,S2p, pch = 21, bg = "white", cex = 1.5)

# ZNGI
arrows(x0 = -100, y0 = N1s, x1 = 100,  y1 = N2s, length = 0, lwd = 0.5, lty = 3)
arrows(x0 = N1s, y0 = -100, x1 = N1s , y1 = 100, length = 0, lwd = 0.5, lty = 3)

arrows(x0 = N1s, y0 = N2s, x1 = 100, y1 = N2s, length = 0, lwd=2)
arrows(x0 = N1s, y0 = N2s, x1 = N1s , y1 = 100, length = 0, lwd = 2)


# Consumption vectors
arrows(x0 = eqN1, y0 = eqN2, x1 = v1, y1 = v2, length = 0.1, lwd=1)
arrows(x0 = eqN1, y0 = eqN2, x1 = S1x, y1 = S2x, length = 0, lwd= 0.5, lty = 3)

arrows(x0 = eqN1p, y0 = eqN2p, x1 = v1p, y1 = v2p, length = 0.1, lwd=1)
arrows(x0 = eqN1p, y0 = eqN2p, x1 = S1p, y1 = S2p, length = 0, lwd= 0.5, lty = 3)
points(S1p,S2p, pch = 21, bg = "white", cex = 1.5)

# Equilibrium points 
points(c(eqN1,eqN1p),c(eqN2,eqN2p),pch = 21, bg = c("black","white"))

# Calculation of the slope
q2 = 1-q1
alpha = q2/q1
lines(c(35,35+10),c(15,15),lty = 3, lwd = 0.5)
lines(c(35+10,35+10),c(15,15+alpha*10),lty = 3, lwd = 0.5)

text(x = S1x, y = S2x*1.1, labels = "S", pos = 3)
text(x = S1p, y = S2p*1.1, labels = "S'", pos = 3)

text(x = 48.5, y = 15, labels = expression(alpha==frac(q[2],q[1])),cex = 0.75)

arrows(x0 = S1x, y0 = S2x*1.1, x1 = S1p , y1 = S2p*1.1, length = 0.1, col= "grey")

mtext(text = "A)", side = 3, line = 0.5, adj = 0, cex = 1.75)

##################################################################
# Panel b) Conditions for coexistence
##################################################################

# ZNGIs
plot(x = NULL, y = NULL, xlim=c(0,50),ylim=c(0,50), xlab = "Nutrient 1", ylab = "Nutrient 2",pch=19,cex = 1.5, cex.axis = 1.5, cex.lab = 1.75)

# ZNGI
N1as = 5
N2as = 15

N1bs = 15
N2bs = 5

arrows(x0 = N1as, y0 = N2as, x1 = 100,   y1 = N2as, length = 0, lwd=2)
arrows(x0 = N1as, y0 = N2as, x1 = N1as , y1 = 100, length = 0, lwd = 2)

arrows(x0 = N1bs, y0 = N2bs, x1 = 100,   y1 = N2bs, length = 0, lwd=2)
arrows(x0 = N1bs, y0 = N2bs, x1 = N1bs , y1 = 100, length = 0, lwd = 2)

# Consumption vectors
arrows(x0 = N2as, y0 = N1bs, x1 = 100,   y1 = 100*1.5, length = 0, lwd=1, lty = 3)
arrows(x0 = N2as, y0 = N1bs, x1 = 100,   y1 = 100*0.5, length = 0, lwd=1, lty = 3)

# Example of supply points
points(x=45,y=20, pch = 21, bg = "black", cex = 1.5)
points(x=35,y=30, pch = 21, bg = "white", cex = 1.5)

text(x = 35, y = 30*1.1, labels = "S'", pos = 3)
text(x = 45, y = 20*1.1, labels = "S", pos = 3)

arrows(x0 = 45, y0 = 20*1.1, x1 = 35 , y1 = 30*1.1, length = 0.1, col= "grey")



mtext(text = "B)", side = 3, line = 0.5, adj = 0, cex = 1.75)

setwd("/Users/DGravel/Documents/Manuscripts/Inprep/ms_RR_space")
dev.copy2eps(file = "R-Ratio_theory.eps")

##################################################################
##################################################################
# Figure 2: Effect of nutrient diffusion
##################################################################
##################################################################

quartz(height = 5, width = 5)
par(mar=c(6,6,2,1))

# The ZNGI, supply points and consumption vector
N1s = 5
N2s = 5
e 	= 0.1
phi = 0.3
d 	= 0.3	

plot(x=NULL,y=NULL,xlim=c(0,50),ylim=c(0,50),xlab="Nutrient 1",ylab= "Nutrient 2",pch=10,col=c("black","darkgrey"),cex = 2,cex.axis = 1.5, cex.lab = 1.75,cex.main = 2)
arrows(x0 = N1s, y0 = N2s, x1 = 100, y1 = N2s, length = 0, lwd=2)
arrows(x0 = N1s, y0 = N2s, x1 = N1s , y1 = 100, length = 0, lwd = 2)
arrows(x0 = N1s, y0 = N2s, x1 = N1s + 100, y1 = N2s + 100*(1-q1)/q1,lty=3,lwd = 0.75)

# Equilibrium points in absence of dispersal + projection vectors
S1x=10
S2x=35
S1y=30
S2y=20

eq = eq_fn(S1x,S2x,S1y,S2y,DN1=0,DN2=0,DP=0,DD=0,q1=0.5)
arrows(x0=eq[[1]]$N1x,y0=eq[[1]]$N2x,x1=eq[[2]][1],y1=eq[[2]][2],length=0,lwd=0.5,lty=3,col="black")
arrows(x0=eq[[1]]$N1y,y0=eq[[1]]$N2y,x1=eq[[2]][3],y1=eq[[2]][4],length=0,lwd=0.5,lty=3,col="darkgrey")
points(c(eq[[1]]$N1x,eq[[1]]$N1y),c(eq[[1]]$N2x,eq[[1]]$N2y),pch=21,cex=2,bg=c("black","darkgrey"))
points(c(eq[[2]][1],eq[[2]][3]),c(eq[[2]][2],eq[[2]][4]),pch=21,cex=2,bg=c("black","darkgrey"))

# Equilibrium points in presence of dispersal
eq = eq_fn(S1x,S2x,S1y,S2y,DN1=0.1,DN2=0.1,DP=0,DD=0,q1=0.5)
arrows(x0=eq[[1]]$N1x,y0=eq[[1]]$N2x,x1=eq[[2]][1],y1=eq[[2]][2],length=0,lwd=0.5,lty=3,col="black")
arrows(x0=eq[[1]]$N1y,y0=eq[[1]]$N2y,x1=eq[[2]][3],y1=eq[[2]][4],length=0,lwd=0.5,lty=3,col="darkgrey")
points(c(eq[[1]]$N1x,eq[[1]]$N1y),c(eq[[1]]$N2x,eq[[1]]$N2y),pch=21,cex=2,col=c("black","darkgrey"),bg="white")
points(c(eq[[2]][1],eq[[2]][3]),c(eq[[2]][2],eq[[2]][4]),pch=21,cex=2,col=c("black","darkgrey"),bg="white")

setwd("/Users/DGravel/Documents/Manuscripts/Inprep/ms_RR_space")
dev.copy2pdf(file = "NutrientDiffusion.pdf")

##################################################################
##################################################################
# Figure 3: Effect of detritus diffusion
##################################################################
##################################################################

quartz(height = 5, width = 5)
par(mar=c(6,6,2,1))

# The ZNGI, supply points and consumption vector
N1s = 5
N2s = 5
e 	= 0.1
phi = 0.5
d 	= 0.3	

plot(x=NULL,y=NULL,xlim=c(0,50),ylim=c(0,50),xlab="Nutrient 1",ylab= "Nutrient 2",pch=10,col=c("black","darkgrey"),cex = 2,cex.axis = 1.5, cex.lab = 1.75,cex.main = 2)
arrows(x0 = N1s, y0 = N2s, x1 = 100, y1 = N2s, length = 0, lwd=2)
arrows(x0 = N1s, y0 = N2s, x1 = N1s , y1 = 100, length = 0, lwd = 2)
arrows(x0 = N1s, y0 = N2s, x1 = N1s + 100, y1 = N2s + 100*(1-q1)/q1,lty=3,lwd = 0.75)

# Equilibrium points in absence of dispersal + projection vectors
S1x=10
S2x=35
S1y=30
S2y=20

eq = eq_fn(S1x,S2x,S1y,S2y,DN1=0,DN2=0,DP=0,DD=0,q1=0.5)
arrows(x0=eq[[1]]$N1x,y0=eq[[1]]$N2x,x1=S1x,y1=S2x,length=0,lwd=0.5,lty=3,col="black")
arrows(x0=eq[[1]]$N1y,y0=eq[[1]]$N2y,x1=S1y,y1=S2y,length=0,lwd=0.5,lty=3,col="darkgrey")
points(c(eq[[1]]$N1x,eq[[1]]$N1y),c(eq[[1]]$N2x,eq[[1]]$N2y),pch=21,cex=2,bg=c("black","darkgrey"))
points(c(eq[[2]][1],eq[[2]][3]),c(eq[[2]][2],eq[[2]][4]),pch=21,cex=2,bg=c("black","darkgrey"))

# Equilibrium points in presence of dispersal
eq = eq_fn(S1x,S2x,S1y,S2y,DN1=0,DN2=0,DP=0,DD=1,q1=0.5)
arrows(x0=eq[[1]]$N1x,y0=eq[[1]]$N2x,x1=eq[[2]][1],y1=eq[[2]][2],length=0,lwd=0.5,lty=3,col="black")
arrows(x0=eq[[1]]$N1y,y0=eq[[1]]$N2y,x1=eq[[2]][3],y1=eq[[2]][4],length=0,lwd=0.5,lty=3,col="darkgrey")
points(c(eq[[1]]$N1x,eq[[1]]$N1y),c(eq[[1]]$N2x,eq[[1]]$N2y),pch=21,cex=2,col=c("black","darkgrey"),bg="white")
points(c(eq[[2]][1],eq[[2]][3]),c(eq[[2]][2],eq[[2]][4]),pch=21,cex=2,col=c("black","darkgrey"),bg="white")

setwd("/Users/DGravel/Documents/Manuscripts/Inprep/ms_RR_space")
dev.copy2pdf(file = "DetritusDiffusion.pdf")

##################################################################
##################################################################
# Figure 4: Effect of producer diffusion
##################################################################
##################################################################

quartz(height = 5, width = 5)
par(mar=c(6,6,2,1))

# The ZNGI, supply points and consumption vector
N1s = 5
N2s = 5
e 	= 0.1
phi = 0.3
d 	= 0.3	

plot(x=NULL,y=NULL,xlim=c(0,50),ylim=c(0,50),xlab="Nutrient 1",ylab= "Nutrient 2",pch=10,col=c("black","darkgrey"),cex = 2,cex.axis = 1.5, cex.lab = 1.75,cex.main = 2)
arrows(x0 = N1s, y0 = N2s, x1 = 100, y1 = N2s, length = 0, lwd=2)
arrows(x0 = N1s, y0 = N2s, x1 = N1s , y1 = 100, length = 0, lwd = 2)
arrows(x0 = N1s, y0 = N2s, x1 = N1s + 100, y1 = N2s + 100*(1-q1)/q1,lty=3,lwd = 0.75)

# Equilibrium points in absence of dispersal + projection vectors
S1x=10
S2x=35
S1y=30
S2y=20

eq = eq_fn(S1x,S2x,S1y,S2y,DN1=0,DN2=0,DP=0,DD=0,q1=0.5)
arrows(x0=eq[[1]]$N1x,y0=eq[[1]]$N2x,x1=S1x,y1=S2x,length=0,lwd=0.5,lty=3,col="black")
arrows(x0=eq[[1]]$N1y,y0=eq[[1]]$N2y,x1=S1y,y1=S2y,length=0,lwd=0.5,lty=3,col="darkgrey")
points(c(eq[[1]]$N1x,eq[[1]]$N1y),c(eq[[1]]$N2x,eq[[1]]$N2y),pch=21,cex=2,bg=c("black","darkgrey"))
points(c(eq[[2]][1],eq[[2]][3]),c(eq[[2]][2],eq[[2]][4]),pch=21,cex=2,bg=c("black","darkgrey"))

# Equilibrium points in presence of dispersal
eq = eq_fn(S1x,S2x,S1y,S2y,DN1=0,DN2=0,DP=10,DD=0,q1=0.5)
arrows(x0=eq[[1]]$N1x,y0=eq[[1]]$N2x,x1=eq[[2]][1],y1=eq[[2]][2],length=0,lwd=0.5,lty=3,col="black")
arrows(x0=eq[[1]]$N1y,y0=eq[[1]]$N2y,x1=eq[[2]][3],y1=eq[[2]][4],length=0,lwd=0.5,lty=3,col="darkgrey")
points(c(eq[[1]]$N1x,eq[[1]]$N1y),c(eq[[1]]$N2x,eq[[1]]$N2y),pch=21,cex=2,col=c("black","darkgrey"),bg="white")
points(c(eq[[2]][1],eq[[2]][3]),c(eq[[2]][2],eq[[2]][4]),pch=21,cex=2,col=c("black","darkgrey"),bg="white")

setwd("/Users/DGravel/Documents/Manuscripts/Inprep/ms_RR_space")
dev.copy2pdf(file = "ProducerDiffusion.pdf")

##################################################################
##################################################################
# Figure 5: Alternate stable states
##################################################################
##################################################################

# Step 1: the model
two_sp_model = function(t,y,pars) {
	with(as.list(c(y,pars)), {

		Gix = min(a1i*N1x/q1i,a2i*N2x/q2i)
		Giy = min(a1i*N1y/q1i,a2i*N2y/q2i)	

		Gjx = min(a1j*N1x/q1j,a2j*N2x/q2j)
		Gjy = min(a1j*N1y/q1j,a2j*N2y/q2j)	

		dN1x 	= e*S1x - e*N1x - q1i*Gix*Pix - q1j*Gjx*Pjx + q1i*phi*d*Dix + q1j*phi*d*Djx - DN1*(N1x-N1y)
		dN2x 	= e*S2x - e*N2x - q2i*Gix*Pix - q2j*Gjx*Pjx + q2i*phi*d*Dix + q2j*phi*d*Djx - DN2*(N2x-N2y)

		dPix 	= Gix*Pix - m*Pix - DP*(Pix-Piy)
		dPjx 	= Gjx*Pjx - m*Pjx - DP*(Pjx-Pjy)

		dDix 	= m*Pix - d*Dix - DD*(Dix-Diy)
		dDjx 	= m*Pjx - d*Djx - DD*(Djx-Djy)
				
		dN1y 	= e*S1y - e*N1y - q1i*Giy*Piy - q1j*Gjy*Pjy + q1i*phi*d*Diy + q1j*phi*d*Djy - DN1*(N1y-N1x)
		dN2y 	= e*S2y - e*N2y - q2i*Giy*Piy - q2j*Gjy*Pjy + q2i*phi*d*Diy + q2j*phi*d*Djy - DN2*(N2y-N2x)

		dPiy 	= Giy*Piy - m*Piy - DP*(Piy-Pix)
		dPjy 	= Gjy*Pjy - m*Pjy - DP*(Pjy-Pjx)

		dDiy 	= m*Piy - d*Diy - DD*(Diy-Dix)
		dDjy 	= m*Pjy - d*Djy - DD*(Djy-Djx)
			
		list(c(dN1x,dN2x,dN1y,dN2y,dPix,dPjx,dPiy,dPjy,dDix,dDjx,dDiy,dDjy))
		})
	}

e=0.1
m = 0.1
N1is = 5
N2is = 10
slope_i = 2
q1i = 1/(slope_i+1)
q2i = 1-q1i
a1i = m*q1i/N1is
a2i = m*q2i/N2is

N1js = 10
N2js = 5
slope_j = 0.5
q1j = 1/(slope_j+1)
q2j = 1-q1j
a1j = m*q1j/N1js
a2j = m*q2j/N2js

phi = 0.9
d = 0.1
DN1 = 0
DN2 = 0
DP = 0
DD = 1

S1x = 17
S2x = 32
S1y = 32
S2y = 17

DD = 1

e = 0.1
phi = 0.4

pars = c(e=e,a1i=a1i,a2i=a2i,q1i=q1i,q2i=1-q1i,a1j=a1j,a2j=a2j,q1j=q1j,q2j=1-q1j,m=m,S1x=S1x,S2x=S2x,S1y=S1y,S2y=S2y,phi=phi,d=d,DN1=DN1,DN2=DN2,DP=DP,DD=DD)

y = c(N1x=S1x,N2x=S2x,N1y=S1y,N2y=S2y,Pix=0.01,Pjx=0.01,Piy=0.01,Pjy=0.01,Dix=0,Djx=0,Diy=0,Djy=0)

eq = as.list(runsteady(y=y, func=two_sp_model, parms=pars)[[1]])

# Ça marche, la dispersion favorise la coexistence. Par contre, elle semble pour l'instant très stable

Gix = min(a1i*eq$N1x/q1i,a2i*eq$N2x/q2i)
Giy = min(a1i*eq$N1y/q1i,a2i*eq$N2y/q2i)	

Gjx = min(a1j*eq$N1x/q1j,a2j*eq$N2x/q2j)
Gjy = min(a1j*eq$N1y/q1j,a2j*eq$N2y/q2j)	

S1x_net = (q1i*Gix*eq$Pix + q1j*Gjx*eq$Pjx)/e + eq$N1x
S2x_net = (q2i*Gix*eq$Pix + q2j*Gjx*eq$Pjx)/e + eq$N2x
S1y_net = (q1i*Giy*eq$Piy + q1j*Gjy*eq$Pjy)/e + eq$N1y
S2y_net = (q2i*Giy*eq$Piy + q2j*Gjy*eq$Pjy)/e + eq$N2y
	
par(mar=c(6,6,2,1))	
plot(c(S1x,S1y),c(S2x,S2y),xlim=c(0,50),ylim=c(0,50),xlab = "Nutrient 1", ylab = "Nutrient 2",pch=10,cex = 2,cex.axis = 1.5, cex.lab = 1.75, cex.main = 2)

# ZNGI species i
arrows(x0 = N1is, y0 = N2is, x1 = 100, y1 = N2is, length = 0, lwd = 2,col="darkblue")
arrows(x0 = N1is, y0 = N2is, x1 = N1is , y1 = 100, length = 0, lwd = 2,col="darkblue")
	
# ZNGI species j
arrows(x0 = N1js, y0 = N2js , x1 = 100, y1 = N2js, length = 0, lwd=2,col="darkred")
arrows(x0 = N1js, y0 = N2js, x1 = N1js , y1 = 100, length = 0, lwd = 2,col="darkred")
	
# Projection vectors
arrows(x0 = N1js, y0 = N2is, x1 = N1js + 100, y1 = N2is + 100*(1-q1i)/q1i,lty=3,lwd = 0.75, col = "darkblue")
arrows(x0 = N1js, y0 = N2is, x1 = N1js + 100, y1 = N2is + 100*(1-q1j)/q1j,lty=3,lwd = 0.75,col = "darkred")

# Net supply points
arrows(x0 = S1y, y0 = S2y, x1 = S1y_net, y1 = S2y_net, lwd = 0.5, length = 0.1)
arrows(x0 = S1x, y0 = S2x, x1 = S1x_net, y1 = S2x_net, lwd = 0.5, length = 0.1)
arrows(x0 = S1y, y0 = S2y, x1 = eq$N1y, y1 = eq$N2y, lwd = 0.5, length = 0.1,lty=3)

points(c(S1x_net,S1y_net),c(S2x_net,S2y_net),pch=21,bg="white",cex=2,lty=3)

# Equilibrium
points(N1js,N2is,pch=19,cex=2)

# In monocultures
# Species i
y = c(N1x=S1x,N2x=S2x,N1y=S1y,N2y=S2y,Pix=0.01,Pjx=0.0,Piy=0.01,Pjy=0.0,Dix=0,Djx=0,Diy=0,Djy=0)
eq = as.list(runsteady(y=y, func=two_sp_model, parms=pars)[[1]])
Gix = min(a1i*eq$N1x/q1i,a2i*eq$N2x/q2i)
Giy = min(a1i*eq$N1y/q1i,a2i*eq$N2y/q2i)	
Gjx = min(a1j*eq$N1x/q1j,a2j*eq$N2x/q2j)
Gjy = min(a1j*eq$N1y/q1j,a2j*eq$N2y/q2j)	
S1x_net = (q1i*Gix*eq$Pix + q1j*Gjx*eq$Pjx)/e + eq$N1x
S2x_net = (q2i*Gix*eq$Pix + q2j*Gjx*eq$Pjx)/e + eq$N2x
S1y_net = (q1i*Giy*eq$Piy + q1j*Gjy*eq$Pjy)/e + eq$N1y
S2y_net = (q2i*Giy*eq$Piy + q2j*Gjy*eq$Pjy)/e + eq$N2y
arrows(x0 = S1y, y0 = S2y, x1 = S1y_net, y1 = S2y_net, lwd = 0.5, length = 0.1,col="darkblue")
arrows(x0 = S1x, y0 = S2x, x1 = S1x_net, y1 = S2x_net, lwd = 0.5, length = 0.1,col="darkblue")
arrows(x0 = S1y, y0 = S2y, x1 = eq$N1y, y1 = eq$N2y, lwd = 0.5, length = 0.1,lty=3,col="darkblue")
points(c(S1x_net,S1y_net),c(S2x_net,S2y_net),pch=21,cex=2,col="darkblue",bg="white")
points(c(eq$N1x,eq$N1y),c(eq$N2x,eq$N2y),pch=21,cex=2,bg="darkblue")	

# Species j
y = c(N1x=S1x,N2x=S2x,N1y=S1y,N2y=S2y,Pix=0.0,Pjx=0.01,Piy=0.0,Pjy=0.01,Dix=0,Djx=0,Diy=0,Djy=0)
eq = as.list(runsteady(y=y, func=two_sp_model, parms=pars)[[1]])
Gix = min(a1i*eq$N1x/q1i,a2i*eq$N2x/q2i)
Giy = min(a1i*eq$N1y/q1i,a2i*eq$N2y/q2i)	
Gjx = min(a1j*eq$N1x/q1j,a2j*eq$N2x/q2j)
Gjy = min(a1j*eq$N1y/q1j,a2j*eq$N2y/q2j)	
S1x_net = (q1i*Gix*eq$Pix + q1j*Gjx*eq$Pjx)/e + eq$N1x
S2x_net = (q2i*Gix*eq$Pix + q2j*Gjx*eq$Pjx)/e + eq$N2x
S1y_net = (q1i*Giy*eq$Piy + q1j*Gjy*eq$Pjy)/e + eq$N1y
S2y_net = (q2i*Giy*eq$Piy + q2j*Gjy*eq$Pjy)/e + eq$N2y
arrows(x0 = S1y, y0 = S2y, x1 = S1y_net, y1 = S2y_net, lwd = 0.5, length = 0.1,col="darkred")
arrows(x0 = S1x, y0 = S2x, x1 = S1x_net, y1 = S2x_net, lwd = 0.5, length = 0.1,col="darkred")
arrows(x0 = S1y, y0 = S2y, x1 = eq$N1y, y1 = eq$N2y, lwd = 0.5, length = 0.1,lty=3,col="darkred")
points(c(S1x_net,S1y_net),c(S2x_net,S2y_net),pch=21,cex=2,col="darkred",bg="white")
points(c(eq$N1x,eq$N1y),c(eq$N2x,eq$N2y),pch=21,cex=2,bg="darkred")
		
setwd("/Users/DGravel/Documents/Manuscripts/Inprep/ms_RR_space")
dev.copy2pdf(file = "AlternateStates.pdf")		

		
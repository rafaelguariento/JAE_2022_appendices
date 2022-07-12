# Code for generating the results of "Prey defense phenotype mediates multiple-predator effects in tri-trophic food-webs"
# After running one of the provided scripts that generate the results, run this script to produce the individual panels represented in Figure 2.

library(plot3D)
library(RColorBrewer)

corA<-brewer.pal(5, "Greens")
corB<-brewer.pal(5, "Reds")
corC<-brewer.pal(5, "Blues")
corD<-brewer.pal(5, "Purples")


hist3D (z = result_V12_r2,  border = "black",zlim = c(0, 51),clim = c(0,51), col=corD,add = F,colkey = F,alpha=0.05,d = 10,theta = 120,phi = 20,bty = "b", scale = T,ticktype = "detailed", opaque.top = T)
hist3D (z = result_A_r2, border = "black", zlim = c(0, 51),clim=c(0,51),col=corA,add = T,colkey = F,alpha=0.05,d = 10,theta = 120,phi = 20,bty = "b", scale = T,ticktype = "detailed", opaque.top = T)
hist3D (z = result_V1_r2,  border = "black",zlim = c(0, 51),clim = c(0,51), col=corB,add = T,colkey = F,alpha=0.05,d = 10,theta = 120,phi = 20,bty = "b", scale = T,ticktype = "detailed", opaque.top = T)
hist3D (z = result_V2_r2,  border = "black", zlim = c(0, 51), clim = c(0,51),col=corC,add = T,colkey = F,alpha=0.05,opaque.top = T)


hist3D (z = result_V1_r3,  border = "black",zlim = c(0, 51),clim = c(0,51), col=corB,add = F,colkey = F,alpha=0.05,d = 10,theta = 120,phi = 20,bty = "b", scale = T,ticktype = "detailed", opaque.top = T)
hist3D (z = result_A_r3, border = "black", zlim = c(0, 51),clim=c(0,51),col=corA,add = T,colkey = F,alpha=0.05,d = 10,theta = 120,phi = 20,bty = "b", scale = T,ticktype = "detailed", opaque.top = T)
hist3D (z = result_V2_r3,  border = "black", zlim = c(0, 51), clim = c(0,51),col=corC,add = T,colkey = F,alpha=0.05,opaque.top = T)
hist3D (z = result_V12_r3,  border = "black",zlim = c(0, 51),clim = c(0,51), col=corD,add = T,colkey = F,alpha=0.05,d = 10,theta = 120,phi = 20,bty = "b", scale = T,ticktype = "detailed", opaque.top = T)


hist3D (z = result_V1_r5,  border = "black",zlim = c(0, 51),clim = c(0,51), col=corB,add = F,colkey = F,alpha=0.05,d = 10,theta = 120,phi = 20,bty = "b", scale = T,ticktype = "detailed", opaque.top = T)
hist3D (z = result_A_r5, border = "black", zlim = c(0, 51),clim=c(0,51),col=corA,add = T,colkey = F,alpha=0.05,d = 10,theta = 120,phi = 20,bty = "b", scale = T,ticktype = "detailed", opaque.top = T)
hist3D (z = result_V2_r5,  border = "black", zlim = c(0, 51), clim = c(0,51),col=corC,add = T,colkey = F,alpha=0.05,opaque.top = T)
hist3D (z = result_V12_r5,  border = "black",zlim = c(0, 51),clim = c(0,51), col=corD,add = T,colkey = F,alpha=0.05,d = 10,theta = 120,phi = 20,bty = "b", scale = T,ticktype = "detailed", opaque.top = T)


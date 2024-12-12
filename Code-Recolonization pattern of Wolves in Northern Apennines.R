
# Code for the analyses described in:

# Recolonization pattern of Wolves in Northern Apennines, Central Italy: 
# A Bayesian Analysis Using Opportunistic and Systematic Data

# OPPORTUNISTIC DATA ANALYSIS

# Load DB.opp.csv
rawdata.opp <- file.choose() 
rawdata.opp <- read.csv(rawdata.opp,h=T,sep=",",stringsAsFactors = F)

# Load Covs.opp.csv
covs.opp <- file.choose() 
covs.opp <- read.csv(covs.opp,h=T,sep=",",stringsAsFactors = F)

#create 3d arrays for modelling

# Data
data.opp <- as.matrix(rawdata.opp[,grepl("V",colnames(rawdata.opp))])

nsites.opp <- length(data.opp[,1])
nsurveys.opp <- length(data.opp[1,])/8
nyears.opp <- 8

db.opp <- array(NA,dim=c(nsites.opp,nsurveys.opp,nyears.opp))

db.opp[,,1] <- data.opp[,1:11]
db.opp[,,2] <- data.opp[,12:22]
db.opp[,,3] <- data.opp[,23:33]
db.opp[,,4] <- data.opp[,34:44]
db.opp[,,5] <- data.opp[,45:55]
db.opp[,,6] <- data.opp[,56:66]
db.opp[,,7] <- data.opp[,67:77]
db.opp[,,8] <- data.opp[,78:88]

# Effort
effort.opp <- as.matrix(covs.opp[,grepl("V",colnames(covs.opp))])

eff.opp <- array(NA,dim=c(nsites.opp,nsurveys.opp,nyears.opp))

eff.opp[,,1] <- effort.opp[,1:11]
eff.opp[,,2] <- effort.opp[,12:22]
eff.opp[,,3] <- effort.opp[,23:33]
eff.opp[,,4] <- effort.opp[,34:44]
eff.opp[,,5] <- effort.opp[,45:55]
eff.opp[,,6] <- effort.opp[,56:66]
eff.opp[,,7] <- effort.opp[,67:77]
eff.opp[,,8] <- effort.opp[,78:88]

#Create visitation data

Visits <-db.opp
#convert all 0s in the original data to 1s and all the NAs to 0s:
Visits[Visits == 0] <- 1
Visits[Visits > 0] <- 1
Visits[is.na(Visits)] <- 0


#Check if data and effort correspond
length(which(db.opp == 0 | db.opp== 1))
length(which(effort.opp > 0))

library(dplyr)
if_else(which(db.opp== 0 | db.opp== 1) == which(eff.opp > 0),1,0)
#All 1s fine!

#Check if the three arrays match: pick random portions of the three
# eff.opp and Visits have 0s where db.opp has NAs

head(db.opp)[,,1]
head(eff.opp)[,,1]
head(Visits)[,,1]

head(db.opp)[,,4]
head(eff.opp)[,,4]
head(Visits)[,,4]

(db.opp)[10:20,,4]
(eff.opp)[10:20,,4]
(Visits)[10:20,,4]

head(db.opp)[,,8]
head(eff.opp)[,,8]
head(Visits)[,,8]

#Seems fine

################################################################################

#create a list with all the data necessary for modelling procedure

str(bdata.opp<-list(y=db.opp,
                    nsites=nsites.opp,nsurveys=nsurveys.opp,nyears=nyears.opp,
                    effort=eff.opp,Visit=Visits,
                    latitude=covs.opp$Latitude,
                    longitude=covs.opp$Longtiude,
                    cluster=covs.opp$cluster,nCat.site=length(unique(covs.opp$cluster)),
                    year=c(-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5)))
                    
#set the initial starting parameters for the modelling procedure to all 1s
zst.opp <- array(1, dim=c(nsites.opp,nyears.opp))
inits.opp <- function() {list(z=zst.opp)}


#Choose the parameters to monitor during MCMC process
params.opp <- c("a.psi","b1","b2","b3","psi","avg.psi",
                "a.detp","b4","p","avg.p",
                "alpha_cluster","b5","ps1","ps2",
                "rdm.site","mu","sigma1","rdm.res",
                "rdm.year","diff.rdm.year","sigma2",
                "Chi2","FT","Chi2rep","FTrep",
                "Chi2.chat","FT.chat","Chi2.bpv","FT.bpv")

library(jagsUI)

#run the model (Rhat < 1.1 indicates convergence)
(out.opp <-jagsUI(bdata.opp, inits.opp, params.opp, "Model Opportunistic Data.txt", n.chains=3, 
                    n.burnin=25000, n.thin=10, n.iter=75000, n.cores=6, parallel=T, DIC=T))

################################################################################
################################################################################

# SYSTEMATIC DATA ANALYSIS

# Load DB.sys.csv
rawdata.sys <- file.choose() 
rawdata.sys <- read.csv(rawdata.sys,h=T,sep=",",stringsAsFactors = F)

# Load Covs.opp.csv
covs.sys <- file.choose() 
covs.sys <- read.csv(covs.sys,h=T,sep=",",stringsAsFactors = F)

#create 3d arrays for modelling

# Data
data.sys <- as.matrix(rawdata.sys[,grepl("V",colnames(rawdata.sys))])

nsites.sys <- length(data.sys[,1])
nsurveys.sys <- length(data.sys[1,])/2
nyears.sys <- 2

db.sys <- array(NA,dim=c(nsites.sys,nsurveys.sys,nyears.sys))

db.sys[,,1] <- data.sys[,1:5]
db.sys[,,2] <- data.sys[,6:10]

# Effort
effort.sys <- as.matrix(covs.sys[,grepl("V",colnames(covs.sys))])

eff.sys <- array(NA,dim=c(nsites.sys,nsurveys.sys,nyears.sys))

eff.sys[,,1] <- effort.sys[,1:5]
eff.sys[,,2] <- effort.sys[,6:10]

################################################################################

#create a list with all the data necessary for modelling procedure
str(bdata.sys <- list(y=db.sys,
                      nsites=nsites.sys,nsurveys=nsurveys.sys,nyears=nyears.sys,
                      effort=eff.sys,
                      latitude=covs.sys$Latitude,
                      longitude=covs.sys$Longtiude))

#set the initial starting parameters for the modelling procedure to all 1s
zst.sys <- array(1,dim=c(nsites.sys,nyears.sys))    # Avoid data/model/inits conflict
inits.sys <- function(){list(z = zst.sys)}

#Choose the parameters to monitor during MCMC process
params.sys <- c("psi.int","b1","b2","b3","avg.psi","psi",
                "a.detp.int","b4","avg.p","p",
                "rdm.site","alpha.mu","alpha.sd",
                "rdm.year","sd.rdm.year",
                "Chi2","FT","Chi2rep","FTrep",
                "Chi2.chat","FT.chat","Chi2.bpv","FT.bpv")

out.sys <- jagsUI(bdata.sys, inits.sys, params.sys, "Model Systematic Data.txt", n.chains=3, 
                    n.burnin=25000, n.thin=10, n.iter=75000, parallel=T,n.cores=6, DIC=T)

out.sys

################################################################################

# Figure 2

# Create a vector of reproductive units scaled to fit in the graph
# Pup response to wh: 2011= NA, 2012= NA, 2013= 0, 2014= 1, 2015= 1, 2016= 1, 2017= 2, 2018= 4, 2019= 2, 2020= 2

pup.dets <- c(NA,NA,0,0.2,0.2,0.2,0.4,0.8,0.4,0.4)

#OCCUPANCY PREDICTED
#Plot average occupancy for each year and BCI
windowsFonts(Times=windowsFont("Times New Roman"))

# Set up plotting parameters, specify Times New Roman, and adjust label positions
par(family="Times", mfrow=c(1,1), mar=c(6, 4, 4, 4) + 0.5, mgp=c(1.75, 0.5, 0))

avg.occu <- cbind(out.opp$sims.list$avg.psi,out.sys$sims.list$avg.psi)
avg.p <- cbind(out.opp$sims.list$avg.p,out.sys$sims.list$avg.p)

#OCCUPANCY MODEL WITH COVS
boxplot(avg.occu, ylab="Probability",xlab="Year",ylim=c(0,1),xlim=c(1,11),axes=F,outline=F,boxwex=0.15,range=1.5,
        col=c(rep("#66bb6a",8),rep("#ff7043",2)),at=c(0.85,1.85,2.85,3.85,4.85,5.85,6.85,7.85,8.85,9.85),cex.lab=1.5)
title(main="Average occupancy and detection probabilities",cex.main=2.5)
axis(1,at=c(1:10),labels=c("2011","2012","2013","2014","2015","2016","2017","2018","2019","2020"),cex.axis=1.25)
axis(2,at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),cex.axis=1.25)
boxplot(avg.p, ylab="Probability",xlab="Year",ylim=c(0,1),xlim=c(1,11),axes=F,outline=F,boxwex=0.15,range=1.5,
        col=c(rep("#bdbdbd",8),rep("white",2)),add=T,at=c(1.15,2.15,3.15,4.15,5.15,6.15,7.15,8.15,9.15,10.15))
axis(4,at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,1,2,3,4,5), line=-7.75,,cex.axis=1.25)
points(pup.dets,pch=18,type="b",col="black",cex=2, lty=9)
mtext("Reproductive units", side=4, line=-6, cex=1.5)
legend("bottom", inset=-0.2, xpd=TRUE, horiz=TRUE, bty="n",
       legend=c("Occupancy (Opportunistic)","Occupancy (Systematic)", 
                "Detection (Opportunistic)","Detection (Systematic)", "Reproductive units"),
       fill=c("#66bb6a", "#ff7043", "#bdbdbd", "white", NA), 
       border=c("black", "black", "black", "black", NA), 
       pch=c(NA, NA, NA, NA, 18), 
       col=c(NA, NA, NA, NA, "black"),
       pt.cex=1,x.intersp=0.25,text.width=c(1.5,1.3,1.4,1.2,1),cex=1.2)

################################################################################

# Figure 3

# Load WH pup responses.csv
pups.db <- file.choose()
pups.db <- read.csv(pups.db,h=T,sep=",",stringsAsFactors = F)

# Create datasets with occupancy estimates at each site and site coordinates

Occu.opp <- data.frame(out.opp$mean$psi,Latitude=covs.opp$Y,Longitude=covs.opp$X)
Occu.sys <- data.frame(out.sys$mean$psi,Latitude=covs.sys$Y,Longitude=covs.sys$X)

library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

# Check available fonts
windowsFonts()

# Create each plot with its specific aesthetics
p1 <- ggplot(Occu.opp, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X1, color = X1)) + 
  theme_bw() + theme_minimal() + xlim(587000,615000) + ylim(4865000,4891000) +
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  theme(legend.position= "none",axis.title.x= element_text(size=16,family="serif"),axis.text.x = element_blank(),
        axis.title.y= element_blank(),axis.text.y = element_text(size=10,family="serif")) +
  xlab("2011")

p2 <- ggplot(Occu.opp, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X2, color = X2)) + 
  theme_bw() + theme_minimal() + xlim(587000,615000) + ylim(4865000,4891000) +
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  theme(legend.position= "none",axis.title.x= element_text(size=16,family="serif"),axis.text.x = element_blank(),
        axis.title.y= element_blank(),axis.text.y = element_blank()) +
  xlab("2012")

p3 <- ggplot(Occu.opp, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X3, color = X3)) + 
  theme_bw() + theme_minimal() + xlim(587000,615000) + ylim(4865000,4891000) + 
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  theme(legend.position= "none",axis.title.x= element_text(size=16,family="serif"),axis.text.x = element_blank(),
        axis.title.y= element_blank(),axis.text.y = element_blank()) +
  xlab("2013")

p4 <- ggplot(Occu.opp, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X4, color = X4)) + 
  theme_bw() + theme_minimal() + xlim(587000,615000) + ylim(4865000,4891000) + 
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  theme(legend.position= "none",axis.title.x= element_text(size=16,family="serif"),axis.text.x = element_blank(),
        axis.title.y= element_blank(),axis.text.y = element_blank()) +
  xlab("2014") +
  geom_point(data=subset(pups.db,year == "2014"),aes(x=x,y=y), pch=18,cex=3.75,col="#03a9f4") # add pup response

p5 <- ggplot(Occu.opp, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X5, color = X5)) + 
  theme_bw() + theme_minimal() + xlim(587000,615000) + ylim(4865000,4891000) + 
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  theme(legend.position= "none",axis.title.x= element_text(size=16,family="serif"),axis.text.x = element_blank(),
        axis.title.y= element_blank(),axis.text.y = element_blank()) +
  xlab("2015") +
  geom_point(data=subset(pups.db,year == "2015"),aes(x=x,y=y), pch=18,cex=3.75,col="#03a9f4") # add pup response

p6 <- ggplot(Occu.opp, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X6, color = X6)) + 
  theme_bw() + theme_minimal() + xlim(587000,615000) + ylim(4865000,4891000) + 
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  theme(legend.position= "none",axis.title.x= element_text(size=16,family="serif"),axis.text.x= element_text(size=10,family="serif"),
        axis.title.y= element_blank(),axis.text.y= element_text(size=10,family="serif")) +
  xlab("2016") +
  geom_point(data=subset(pups.db,year == "2016"),aes(x=x,y=y), pch=18,cex=3.75,col="#03a9f4") # add pup response

p7 <- ggplot(Occu.opp, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X7, color = X7)) + 
  theme_bw() + theme_minimal() + xlim(587000,615000) + ylim(4865000,4891000) + 
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  theme(legend.position= "none",axis.title.x= element_text(size=16,family="serif"),axis.text.x= element_text(size=10,family="serif"),
        axis.title.y= element_blank(),axis.text.y = element_blank()) +
  xlab("2017") +
  geom_point(data=subset(pups.db,year == "2017"),aes(x=x,y=y), pch=18,cex=3.75,col="#03a9f4") # add pup response

p8 <- ggplot(Occu.opp, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X8, color = X8)) + 
  theme_bw() + theme_minimal() + xlim(587000,615000) + ylim(4865000,4891000) + 
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  theme(legend.position= "none",axis.title.x= element_text(size=16,family="serif"),axis.text.x= element_text(size=10,family="serif"),
        axis.title.y= element_blank(),axis.text.y = element_blank()) +
  xlab("2018") +
  geom_point(data=subset(pups.db,year == "2018"),aes(x=x,y=y), pch=18,cex=3.75,col="#03a9f4") # add pup response

p9 <- ggplot(Occu.sys, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X1, color = X1)) + 
  theme_bw() + theme_minimal() + xlim(587000,615000) + ylim(4865000,4891000) + 
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  theme(legend.position= "none",axis.title.x= element_text(size=16,family="serif"),axis.text.x= element_text(size=10,family="serif"),
        axis.title.y= element_blank(),axis.text.y = element_blank()) +
  xlab("2019") +
  geom_point(data=subset(pups.db,year == "2019"),aes(x=x,y=y), pch=18,cex=3.75,col="#03a9f4") # add pup response

p10 <- ggplot(Occu.sys, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X2, color = X2)) + 
  theme_bw() + theme_minimal() + xlim(587000,615000) + ylim(4865000,4891000) + 
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  theme(legend.position= "none",axis.title.x= element_text(size=16,family="serif"),axis.text.x= element_text(size=10,family="serif"),
        axis.title.y= element_blank(),axis.text.y = element_blank()) +
  xlab("2020") +
  geom_point(data=subset(pups.db,year == "2020"),aes(x=x,y=y), pch=18,cex=3.75,col="#03a9f4") # add pup response

# Define a "ghost" plot to extract the legend
plot.legend <- ggplot(Occu.sys, aes(x = Longitude, y = Latitude)) + 
  geom_point(aes(size = X1, color = X1)) +
  scale_size(range = c(1, 5), limits = c(0, 1)) +
  scale_color_viridis_c(limits = c(0, 1), option = "rocket") +
  guides(
    size = guide_legend(title = "Occupancy probability", keywidth = 3, keyheight = 1.5), 
    color = guide_legend(title = "Occupancy probability", keywidth = 3, keyheight = 1.5)
  ) +
  theme_void() +  
  theme(
    legend.position = "bottom",  
    legend.direction = "horizontal",  
    plot.margin = margin(t = 10, b = 10, l = 10, r = 10),
    legend.text = element_text(size = 18, family="serif"), # Adjust the text size of legend
    legend.title = element_text(size = 18, family="serif"), # text size of legend title
    legend.spacing.x = unit(0.5, "lines"), 
    legend.key.width = unit(1, "lines") 
  )

g_legend <- ggplot_gtable(ggplot_build(plot.legend))$grobs[[which(sapply(ggplot_gtable(ggplot_build(plot.legend))$grobs, function(x) x$name) == "guide-box")]]

latitude_label <- textGrob("Latitude", gp = gpar(fontsize = 22, fontfamily = "serif"),  rot = 90) # axis name rotated vertically
longitude_label <- textGrob("Longitude", gp = gpar(fontsize = 22, fontfamily = "serif"))
plot_title <- textGrob("Wolf Occupancy", gp = gpar(fontsize = 28, fontfamily = "serif", fontface = "bold"))

# Arrange the plots and add the custom labels, including the title
grid.arrange(
  arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol = 5),
  bottom = longitude_label,
  left = latitude_label,
  g_legend,
  top = plot_title,  
  nrow = 2,
  heights = c(11, 0.6)  
)

################################################################################

# Appendix 2

windowsFonts(Times=windowsFont("Times New Roman"))

par(family="Times",mfrow = c(2, 2))

# First plot
plot(x = out.opp$sims.list$Chi2, y = out.opp$sims.list$Chi2rep, xlab = " ", ylab = "Simulated", cex.lab = 1.5,cex.axis=1.25)
title(main = "Opportunistic data model", line = 2, cex.main = 2)  # Move title up
text(x = 80, y = 170, labels = "bpv= 0.25, c.hat=1.15", cex = 1.25)
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("Chi-squared Discrepancy", side = 3, line = 0.5, cex = 1.25)

# Second plot
plot(x = out.sys$sims.list$Chi2, y = out.sys$sims.list$Chi2rep, xlab = " ", ylab = " ", cex.main = 2,cex.axis=1.25)
title(main = "Systematic data model", line = 2, cex.main = 2)  # Move title up
text(x = 30, y = 102, labels = "bpv= 0.05, c.hat=1.48", cex = 1.25)
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("Chi-squared Discrepancy", side = 3, line = 0.5, cex = 1.25)

# Third plot
plot(x = out.opp$sims.list$FT, y = out.opp$sims.list$FTrep, xlab = "Observed", ylab = "Simulated", cex.lab = 1.5,cex.axis=1.25)
title(main = " ", line = 2, cex.main = 1.5)  # Empty title if you want to skip
text(x = 22, y = 62, labels = "bpv= 0.23, c.hat=1.14", cex = 1.25)
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("Freeman-Tukey Discrepancy", side = 3, line = 0.5, cex = 1.25)

# Fourth plot
plot(x = out.sys$sims.list$FT, y = out.sys$sims.list$FTrep, xlab = "Observed", ylab = " ", cex.lab = 1.5,cex.axis=1.25)
title(main = " ", line = 2, cex.main = 1.5)  # Empty title if you want to skip
text(x = 9, y = 24.5, labels = "bpv= 0.23, c.hat=1.44", cex = 1.25)
abline(a = 0, b = 1, col = "red", lwd = 2)
mtext("Freeman-Tukey Discrepancy", side = 3, line = 0.5, cex = 1.25)

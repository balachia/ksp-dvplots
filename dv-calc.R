library(data.table)
library(ggplot2)
library(gtools)
library(lattice)
library(latticeExtra)
library(RColorBrewer)

rm(list=ls())

# kerbin gravity
kg <- 9.81

# lines of constant cargo/stage weight ratio, as a function of cargo weight
iso.csr <- function(mass, csr, eng.mass=0, full.mass=1) {
        # returns the number of tanks that will produce a given cargo/stage ratio
        (mass/csr - eng.mass) / full.mass
    }

# returns number of tanks that will produce a given thrust/weight ratio, for a given cargo mass
iso.twr <- function(mass, g=kg, eng.thr=1, eng.mass=0, full.mass=1, twr=1) {
        # returns number of tanks that will produce thrust/weight ratio of 1 in gravity
        ((eng.thr / (g * twr)) - eng.mass - mass) / full.mass
    }

# the rocket equation, with isp in seconds
dv.calc <- function(mass, isp, g=kg, eng.mass=0, ntanks=1, full.mass=1, empty.mass=1) {
        isp * kg * log((eng.mass + mass + ntanks * full.mass) / (eng.mass + mass + ntanks * empty.mass))
    }

# load engines file
edt <- fread('engines.csv')
setkey(edt,name)

# set some defaults
# set default gravity to kerbin (9.81 m/s)
edt[is.na(g), g := kg]

engs <- edt[,name]
#engs <- engs[1]
# engs <- c('xenon')
#eng.name <- engs[1]

# desired iso-csr curves, with highlighting of the 1:2 ratio
csr.levels <- c(1/10, 1/5, 1/4, 1/3, 1/2, 1, 1.5, 2)
csr.levels <- c(csr.levels, rep(1/2,4))
twr.levels <- seq(1,2.4,0.2)
twr.levels <- c(1,1.2,1.5,2,2.5,3,4,5)

# start the reactor
# free mars
res <- lapply(engs, function (eng.name) {
        print(eng.name)

        # load engine characteristics and set graph properties
        eng <- as.list(edt[eng.name])
        topm.max <- 14
        topm.max <- eng$topm
        tank.seq <- seq(eng$tank.min,eng$tank.max,eng$tank.gap)

        # fix tick marks along x axis
#         xbreak <- max(1, signif(topm.max/10,1))
        xbreak <- signif(topm.max/10,1)
        topm.seq <- seq(0,topm.max,xbreak/100)

        dva.twr <- function(x, twr) { 
                ntanks.twr <- iso.twr(x, g=eng$g, eng.thr=eng$thr, eng.mass=eng$eng.mass,
                    full.mass=eng$mf, twr=twr)
                dv.calc(x, isp=eng$ispa, g=eng$g, eng.mass=eng$eng.mass, ntanks=ntanks.twr,
                    full.mass=eng$mf, empty.mass=eng$me)
            }
        dvv.twr <- function(x, twr) { 
                ntanks.twr <- iso.twr(x, g=eng$g, eng.thr=eng$thr, eng.mass=eng$eng.mass,
                    full.mass=eng$mf, twr=twr)
                dv.calc(x, isp=eng$ispv, g=eng$g, eng.mass=eng$eng.mass, ntanks=ntanks.twr,
                    full.mass=eng$mf, empty.mass=eng$me)
            }
        dva.csr <- function(x, csr) {
                ntanks.csr <- iso.csr(x, csr, eng.mass=eng$eng.mass, full.mass=eng$mf)
                dv.calc(x, isp=eng$ispa, g=eng$g, eng.mass=eng$eng.mass, ntanks=ntanks.csr,
                    full.mass=eng$mf, empty.mass=eng$me)
            }
        dvv.csr <- function(x, csr) {
                ntanks.csr <- iso.csr(x, csr, eng.mass=eng$eng.mass, full.mass=eng$mf)
                dv.calc(x, isp=eng$ispv, g=eng$g, eng.mass=eng$eng.mass, ntanks=ntanks.csr,
                    full.mass=eng$mf, empty.mass=eng$me)
            }
        dva.twr1 <- function(x) dva.twr(x,1)
        dvv.twr1 <- function(x) dvv.twr(x,1)

        # topweight graphs
        twdt <- as.data.table(expand.grid(ntanks=tank.seq,
                                          env=c('ispa','ispv'), 
                                          topm=seq(0,topm.max,xbreak/100)))
        twdt[, isp := sapply(as.character(env), function(x) eng[[x]])]
        twdt[, mass.f := eng$eng.mass + ntanks*eng$mf]
        twdt[, mass.e := eng$eng.mass + ntanks*eng$me]
        twdt[, dv := isp * kg * log((mass.f + topm) / (mass.e + topm))]
        twdt[, twr := eng$thr / (mass.f + topm) / eng$g]
        twdt[, cargo.stage.ratio := topm/mass.f]

#         print(twdt)
        ymax <- twdt[, max(dv)]
        ymax.digits <- floor(log10(ymax)) - 1
        ymax <- trunc(ymax / (10^ymax.digits)) + 1 
        ymax <- ymax * 10^ymax.digits

        # alternatively set y-max to theoretical dv limit
        dv.max <- eng$ispv * kg * log(eng$mf / eng$me) 
        ymax <- dv.max

        ybreak <- max(1, signif(ymax/10,1))

        twdt[, linelab := ifelse(topm %% xbreak == 0, as.character(round(twr,2)), '')]
        twdt[, linelab := ifelse(topm %% xbreak == 0, as.character(round(cargo.stage.ratio,2)), '')]

        # make isoquants - fun math
        # one version for atmosphere, one for vacuum
        isodt <- CJ(topm=seq(0,topm.max,xbreak/25),
                    dv=seq(0,dv.max,ybreak/25))
        isodt[, Ka := exp(dv / (kg * eng$ispa))]
        isodt[, Kv := exp(dv / (kg * eng$ispv))]
        isodt[, na := ((topm + eng$eng.mass)*(Ka-1))/(eng$mf - Ka*eng$me)]
        isodt[, nv := ((topm + eng$eng.mass)*(Kv-1))/(eng$mf - Kv*eng$me)]

        # csr isoquants
        isodt[, csra := topm / (eng$eng.mass + na*eng$mf)]
        isodt[, csrv := topm / (eng$eng.mass + nv*eng$mf)]

        # twr isoquants
        # TODO: fix to account for TWR of relevant bodies
        # TODO: something is still wrong here, we appear to be off by a factor of something
        #isodt[, twra := (eng$thr * (eng$ispa/eng$ispv)) / (kg*(topm + eng$eng.mass + na*eng$mf))]
        #isodt[, twrv := eng$thr / (kg*(topm + eng$eng.mass + nv*eng$mf))]
        isodt[, twra := (eng$thr * (eng$ispa/eng$ispv)) / (eng$g*(topm + eng$eng.mass + na*eng$mf))]
        isodt[, twrv := eng$thr / (eng$g*(topm + eng$eng.mass + nv*eng$mf))]

        #print(eng$g)
        #print(eng$mf)
        #print(isodt[topm==0 & (dv %% ybreak == 0)])


        ############################################################
        # lattice plots
        ############################################################

        # set palettes and legends
        ncolor <- (eng$tank.max / eng$tank.gap) - (eng$tank.min / eng$tank.gap) + 1
        palbase <- colorRampPalette(brewer.pal(9,'Blues'))(2 * ncolor)
        palfunc <- function(x) palbase[ncolor + x/eng$tank.gap]
        twrpal <- brewer.pal(length(twr.levels),'YlGn')
        tankkey <- list(text=list(as.character(seq(eng$tank.min,eng$tank.max,eng$tank.gap))),
                      lines=list(col=palbase[1:ncolor + ncolor],lwd=4),
                      space='right',title=paste0('Tanks (',eng$mf,')'),
                      cex.title=1.4)
        twrkey <- list(col=twrpal,at=twr.levels,
                       space='top')
        plotleg <- list(right=c(fun=draw.key(tankkey)),
                        top=c(fun=draw.colorkey(twrkey)))
        plotleg <- list(right=c(fun=draw.key(tankkey)))
        plotleg <- list(right=list(fun=draw.key,args=list(key=tankkey)),
                        top=list(fun=draw.colorkey,args=list(key=twrkey)))

        # start by plotting dv-lines per stage mass
        plotatm <- xyplot(dv ~ topm, data=twdt[env=='ispa'],
                          groups=ntanks, col=palfunc(twdt[env=='ispa',ntanks]),
                          type='l', lwd=4)
        plotvac <- xyplot(dv ~ topm, data=twdt[env=='ispv'],
                          groups=ntanks, col=palfunc(twdt[env=='ispv',ntanks]),
                          type='l', lwd=4)
        # add csr labels to lines
        plotatm <- plotatm + layer(panel.text(topm,dv,labels=linelab,pos=3),
                                   data=twdt[env=='ispa'])
        plotvac <- plotvac + layer(panel.text(topm,dv,labels=linelab,pos=3),
                                   data=twdt[env=='ispv'])
        # add isoquants
        csratm <- contourplot(csra~topm*dv,data=isodt,
                    at=csr.levels,labels=FALSE,lwd=2,
                    col=brewer.pal(3,'Accent')[2],
                    colorkey=FALSE)
        csrvac <- contourplot(csrv~topm*dv,data=isodt,
                    at=csr.levels,labels=FALSE,lwd=2,
                    col=brewer.pal(3,'Accent')[2],
                    colorkey=FALSE)
        twratm <- levelplot(twra~topm*dv,data=isodt,
                    at=twr.levels,col.regions=twrpal,alpha.regions=0.5,
                    contour=TRUE,
                    colorkey=FALSE)
        twrvac <- levelplot(twrv~topm*dv,data=isodt,
                    at=twr.levels,col.regions=twrpal,alpha.regions=0.5,
                    contour=TRUE,
                    colorkey=FALSE)
        plotatm <- twratm + as.layer(csratm) + as.layer(plotatm)
        plotvac <- twrvac + as.layer(csrvac) + as.layer(plotvac)
        #plotatm <- plotatm + as.layer(csratm)

        # add key
        #plotatm <- update(plotatm,key=tankkey)
        #plotatm <- update(plotatm,key=twrkey)
        plotatm <- update(plotatm,legend=plotleg,
                          main=paste0(eng$name,
                              ' atm (thr:mass ', round(eng$thr * eng$ispa/eng$ispv), ' : ', eng$eng.mass,
                              ', isp ', eng$ispa, ', ',  round(eng$g/kg,2), ' g)'),
                          scales=list(x=list(at=seq(0,topm.max,xbreak)),
                                      y=list(at=seq(0,ymax,ybreak))))
        plotvac <- update(plotvac,legend=plotleg,
                          main=paste0(eng$name,
                              ' vac (thr:mass ', eng$thr, ' : ', eng$eng.mass,
                              ', isp ', eng$ispv, ', ',  round(eng$g/kg,2), ' g)'),
                          scales=list(x=list(at=seq(0,topm.max,xbreak)),
                                      y=list(at=seq(0,ymax,ybreak))))
#         plotatm <- contourplot(csra~topm*dv,data=isodt) + as.layer(plotatm)


        ############################################################
        # ggplots
        ############################################################

        ggp.wrap <- defmacro(ggp.base, csr.levels, expr={
                ggp.base +
                scale_x_continuous(breaks=seq(0,topm.max,xbreak)) + 
                scale_y_continuous(breaks=seq(0,ymax,ybreak), limits=c(0,ymax)) +
                geom_text(hjust=0, vjust=0) +
                geom_point(alpha=0.1) + 
                scale_color_gradient(name=paste0('Tank Masses (',eng$mf,')')) + 
                geom_hline(yintercept=dv.max)
            })

        ggpa <- ggp.wrap({
                ggplot(twdt[env=='ispa'], aes(topm, dv, color=ntanks, label=linelab)) + 
                labs(title=paste0(eng$name,
                    ' (thr:mass ', eng$thr, ' : ', eng$eng.mass,
                    ', atm ', eng$ispa, ', ',  round(eng$g/kg,2), ' kg)')) +
                stat_function(fun=dva.twr1, size=0, geom='area', alpha=0.3, color='black') + 
                stat_function(fun=dva.csr, args=list(csr=csr.levels),
                    geom='point', size=0.5, alpha=0.5, n=1000*length(csr.levels))
            }, csr.levels)

        ggpv <- ggp.wrap({
                ggplot(twdt[env=='ispv'], aes(topm, dv, color=ntanks, label=linelab)) + 
                labs(title=paste0(eng$name,
                    ' (thr:mass ', eng$thr, ' : ', eng$eng.mass,
                    ', vac ', eng$ispv, ', ',  round(eng$g/kg,2), ' kg)')) +
                stat_function(fun=dvv.twr1, size=0, geom='area', alpha=0.3, color='black') + 
                stat_function(fun=dvv.csr, args=list(csr=csr.levels),
                    geom='point', size=0.5, alpha=0.5, n=1000*length(csr.levels)) 
            }, csr.levels)
            

        ggpa <- ggplot(twdt[env=='ispa'], aes(topm, dv, color=ntanks, label=linelab)) + 
            scale_x_continuous(breaks=seq(0,topm.max,xbreak)) + scale_y_continuous(breaks=seq(0,ymax,ybreak), limits=c(0,ymax)) +
            geom_text(hjust=0, vjust=0) +
            stat_function(fun=dva.twr1, size=0, geom='area', alpha=0.3, color='black') + 
            geom_point(alpha=0.1) + 
            scale_color_gradient(name=paste0('Tank Masses (',eng$mf,')')) + 
            geom_hline(yintercept=dv.max) + 
            labs(title=paste0(eng$name, ' (t:m ', eng$thr, ' : ', eng$eng.mass, ', atm ', eng$ispa, ', ',  round(eng$g/kg,2), ' kg)'))

        ggpv <- ggplot(twdt[env=='ispv'], aes(topm, dv, color=ntanks, label=linelab)) + 
            scale_x_continuous(breaks=seq(0,topm.max,xbreak)) + scale_y_continuous(breaks=seq(0,ymax,ybreak), limits=c(0,ymax)) +
            geom_text(hjust=0, vjust=0) +
            stat_function(fun=dvv.twr1, size=0, geom='area', alpha=0.3, color='black') + 
            geom_point(alpha=0.1) + 
            scale_color_gradient(name=paste0('Tank Masses (',eng$mf,')')) + 
            geom_hline(yintercept=dv.max) + 
            labs(title=paste0(eng$name, ' (t:m ', eng$thr, ' : ', eng$eng.mass, ', vac ', eng$ispv, ', ',  round(eng$g/kg,2), ' kg)'))

        ggpa
        ggpv

        png(paste0('dvplots/', eng$name, '-atm.png'), width=800, height=600)
        print(ggpa)
        dev.off()

        png(paste0('dvplots/', eng$name, '-vac.png'), width=800, height=600)
        print(ggpv)
        dev.off()

        png(paste0('dvplots2/', eng$name, '-atm.png'), width=800, height=600)
        print(plotatm)
        dev.off()

        png(paste0('dvplots2/', eng$name, '-vac.png'), width=800, height=600)
        print(plotvac)
        dev.off()



        res <- list()
        res$name <- eng.name
#         res$atm <- ggpa
#         res$vac <- ggpv
        res
    })

# dva
# ggplot(twdt, aes(topm, dva)) + geom_point()

# dvv
# ggplot(twdt, aes(topm, dvv)) + geom_point()

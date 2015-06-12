library(data.table)
library(ggplot2)
library(gtools)
library(lattice)
library(latticeExtra)

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
# engs <- engs[1]
# engs <- c('xenon')
# eng.name <- engs[1]

# desired iso-csr curves, with highlighting of the 1:2 ratio
csr.levels <- c(1/10, 1/5, 1/4, 1/3, 1/2, 1, 1.5, 2)
csr.levels <- c(csr.levels, rep(1/2,4))

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

        # cargo/stage weight ratio contour
#         csrdt <- as.data.table(expand.grid(topm=seq(0,topm.max, xbreak / 10), dv=seq(0,dv.max, ybreak/10)))
#         csrdt[, Ka := exp(dv / (kg*eng$ispa))]
#         csrdt[, Kv := exp(dv / (kg*eng$ispv))]
#         csrdt[, na := ((eng$eng.mass + topm)*(Ka-1))/(eng$mf - Ka*eng$me)]
#         csrdt[, csra := topm/(eng$eng.mass + na * eng$mf)]
#         csrdt[, nv := ((eng$eng.mass + topm)*(Kv-1))/(eng$mf - Ka*eng$me)]
#         csrdt[, csrv := topm/(eng$eng.mass + nv * eng$mf)]
#         csrdt[csra < 0, csra := 0]
#         csrdt[csrv < 0, csrv := 0]

        #plotatm <- xyplot(dv ~ topm, data=twdt, groups=ntanks, type='l')

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
        #print(plotatm)
        print(ggpa)
        dev.off()

        png(paste0('dvplots/', eng$name, '-vac.png'), width=800, height=600)
        print(ggpv)
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

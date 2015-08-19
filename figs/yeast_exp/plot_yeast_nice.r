library("lattice")
library("colorRamps")
library("RColorBrewer")
library("grid")
source("/home/bram/R/src/bramlib.r")

heights <- c(0.3, 1, 0.4, 1, 0.4, 1, 0.3)
widths <- c(0.2,1,0.2)

line.lwd <- 0.25

lwd <- 0.5
tick.cex <- 0.55
label.cex <- 0.7
legend.cex <- 0.55
in.plot.label.cex <- 0.55
plot.tck <- -0.4
xyplot.lwd <- 0.5
points.lwd <- 0.3

if (!exists("dat"))
{
    dat <-  read.table("summary_yeast.csv",sep=";",header=T) 
}

lo <- grid.layout(
                    ncol=length(widths),
                    nrow=length(heights),
                    heights=heights,
                    widths=widths
                    )

spores <- c("spore1","spore2","spore3","spore4")
colors <- brewer.pal(4,"YlOrBr")

panel.xh <- function(x,y,the.data,...)
{
    nmito <- sort(unique(the.data$nmito_min))
    pars <- trellis.par.get()

    pars$box.rectangle$col <- "black"
    pars$box.rectangle$lwd <- line.lwd
    pars$box.umbrella$col <- "black"
    pars$box.umbrella$lwd <- line.lwd
    pars$box.umbrella$lty <- 1
    pars$plot.symbol$pch <- 20
    pars$plot.symbol$col <- "gray"
    pars$plot.symbol$fill <- "gray"
    pars$plot.symbol$lwd <- 0
    pars$plot.symbol$alpha <- 0.5
    pars$plot.symbol$cex <- 0.3

    trellis.par.set(pars)


    for (nmito_i in nmito)
    {
        for(j in 1:length(spores))
        {
            sub.dat <- the.data[the.data$nmito_min == nmito_i,spores[[j]]]
                        
            panel.bwplot(x=
                as.character(
                        rep(nmito_i + j-2.5,times=length(sub.dat))
                ),
                y=sub.dat,
                horizontal=F,
                pch="|",
                fill=colors[[j]]
            )
                            
        }
    }
}

block <- function(
        row, 
        col, 
        dataset,
        label="A",
        label.main="",
        xlab="",
        ylab="")
{
    xp <- xyplot(c(0,1) ~ c(7,44),
                    panel=panel.xh,
                    the.data=dataset)

    pushViewport(viewport(layout.pos.row=row,
                            layout.pos.col=col,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        do.call("panel.xh",trellis.panelArgs(xp,1))

        grid.text(x=0.5,y=1.1,label=label.main,gp=gpar(cex=label.cex))

        #grid.rect(gp=gpar(lwd=lwd,fill="transparent"))
        grid.lines(x=unit(c(0,1),"npc"),y=unit(c(0,0),"npc"),gp=gpar(lwd=line.lwd))
        grid.lines(x=unit(c(0,0),"npc"),y=unit(c(0,1),"npc"),gp=gpar(lwd=line.lwd))

        grid.text(x=0.02,y=1.1,label=label,gp=gpar(cex=label.cex*1.3))
    upViewport()
    
    pushViewport(viewport(layout.pos.row=row+1,
                            layout.pos.col=col,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$x.limits,
                        side="top",
                        labels=T,
                        distance=0.5,
                        cex=tick.cex,
                        labelcex=label.cex,
                        tck=plot.tck,
                        lwd=line.lwd,
                        y.text.off=0.4,
                        nsub=0,
                        text=xlab)
    upViewport()
    
    pushViewport(viewport(layout.pos.row=row,
                            layout.pos.col=col-1,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            expr <- ylab
        single.axis(range=xp$y.limits,
                        side="right",
                        tck=plot.tck,
                        distance=0.5,
                        labels=T,
                        cex=tick.cex,
                        lwd=line.lwd,
                        labelcex=label.cex,
                        x.text.off=0.3,
                        nsub=5,
                        text=expr)
    upViewport()
}


# loop through the different values of ch
for (type.i in sort(unique(dat$type)))
{
    type.desc <- c("selection against heteroplasmy",
                    "selection favoring parent-2 mitochondria",
                    "selection favoring heteroplasmy")

    type.desc.brief <- c("sel_homoplasmy",
                        "sel_parent_2",
                        "sel_heteroplasmy")


    for (ch.val in sort(unique(dat$ch)))
    {
        for (error.val in sort(unique(dat$nmito_error)))
        {
            init.plot(filename=
                    paste("fig_boxplots_eM_",error.val,"_ch",ch.val,"_",
                            type.desc.brief[[type.i+1]],sep=""),
                    font="myriad",width=500,height=500
                    )

            pushViewport(viewport(layout=lo))

           
                type.label <- type.desc[[type.i+1]]
                if (ch.val == 0)
                {
                    type.label <- "no selection"
                }
                grid.text(x=0.5,y=0.98,label=
                        eval(substitute(
                                        expression(
                                                paste(type_label," ",
                                                        epsilon[italic(M)]," = ",error_val,"; ",
                                                        italic(s)," = ",ch_val
                                                )
                                        ),list(type_label=type.label,
                                                error_val=error.val,
                                                ch_val=ch.val)
                                )
                        )
                )

                block(row=2,col=2,
                        dataset=dat[dat$ascus_slope == 0 
                                    & dat$nmito_error == error.val 
                                    & dat$ch == ch.val 
                                    & dat$type == type.i,],
                                    #            xlab=expression(paste("number of mitochondria per individual, ",italic(M))),
                        ylab=expression(paste("frequency parent-1 mitochondria, ",italic(p))),
                        label.main="no mixing gradient (random)"
                        )

                
                pushViewport(viewport(layout.pos.row=2,
                                        layout.pos.col=3
                                        ))

                    draw.my.bw.box(x.left.pos=0.1,
                                    box.width=0.05,
                                    box.height=0.1,
                                    lwd=line.lwd,
                                    box.col=colors[[1]],
                                    whiskers=F)
                    grid.text(x=0.19,y=.95,just="left",label="spore 1",gp=gpar(cex=legend.cex))

                    draw.my.bw.box(x.left.pos=0.1,
                                    y.top.pos=0.85,
                                    box.width=0.05,
                                    box.height=0.1,
                                    lwd=line.lwd,
                                    box.col=colors[[2]],
                                    whiskers=F)
                    grid.text(x=0.19,y=0.8,just="left",label="spore 2",gp=gpar(cex=legend.cex))
                    
                    draw.my.bw.box(x.left.pos=0.1,
                                    y.top.pos=.7,
                                    box.width=0.05,
                                    box.height=0.1,
                                    lwd=line.lwd,
                                    box.col=colors[[3]],
                                    whiskers=F)
                    grid.text(x=0.19,y=0.65,just="left",label="spore 3",gp=gpar(cex=legend.cex))
                    
                    draw.my.bw.box(x.left.pos=0.1,
                                    y.top.pos=0.55,
                                    box.width=0.05,
                                    box.height=0.1,
                                    lwd=line.lwd,
                                    box.col=colors[[4]],
                                    whiskers=F)
                    grid.text(x=0.19,y=0.5,just="left",label="spore 4",gp=gpar(cex=legend.cex))
                upViewport()
                
                block(row=4,col=2,
                        label="B",
                        dataset=dat[dat$ascus_slope == 0.5 
                                    & dat$nmito_error == error.val 
                                    & dat$ch == ch.val 
                                    & dat$type == type.i,],
                                    #            xlab=expression(paste("number of mitochondria per individual, ",italic(M))),
                        ylab=expression(paste("frequency parent-1 mitochondria, ",italic(p))),
                        label.main="slight mixing gradient"
                        )
                
                block(row=6,col=2,
                        label="C",
                        dataset=dat[dat$ascus_slope == 0.9 
                                    & dat$nmito_error == error.val 
                                    & dat$ch == ch.val 
                                    & dat$type == type.i,],
                        xlab=expression(paste("average number of mitochondria per cell, ",italic(bar(M)))),
                        ylab=expression(paste("frequency parent-1 mitochondria, ",italic(p))),
                        label.main="strong mixing gradient"
                        )


            upViewport()
            exit.plot()
        }
    }
}

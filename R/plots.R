bySpeciesComposition <- function(dataset,cnt){
    cols=palette.colors(length(unique(dataset$new_txgroups)),"Pastel 2")
    ngroups=length(unique(na.omit(dataset$new_area)))
    par(mfrow=c(1,ngroups),oma=c(3,1,0,0),mar=c(1,1,1,1))
    for(i in 1:ngroups){
        region=dataset[dataset$new_area == i,]
        region.freq <- t(tapply(region[[cnt]],list(region$new_periods,region$new_txgroups),function(i)sum(i,na.rm=T)))
       region.freq <- apply(region.freq,2,function(i)i/sum(i,na.rm=T))
        barplot(region.freq[,ncol(region.freq):1],legend=(i==ngroups),args.legend=list(bg="white",cex=.75),space=0,main=paste("Area",i),col=cols,lwd=.1,border=.6)
    }
}

countTotal<- function(dataset,cnt){
    ngroups=length(unique(na.omit(dataset$new_area)))
    areacol <- palette.colors(length(unique(dataset$new_area)),"Pastel 1",recycle=T)
    counts <- tapply(dataset[[cnt]], list(dataset$new_periods,dataset$new_area ),sum,na.rm=T)
    plot(counts[,1],type="n",cex=3,lwd=3,ylim=range(counts,na.rm=T),col=areacol[1],ylab=cnt,xaxt="n",xlab="")
    for(i in 1:ncol(counts))lines(rev(counts[,i]),type="o",cex=3,lwd=3,ylim=range(counts),col=areacol[i])
    axis(1,    at=seq_along(levels(dataset$new_periods)),label=rev(levels(dataset$new_periods)))
    legend("topright",lwd=3,pch=1,col=areacol,legend=paste0("area",1:ngroups))
}


plot2dim <- function(ca.res,ngroup){
           cares <- ca.res$row$coord
           if((is.null(dim(cares))))
               par(mfrow=c(1,1))
           else
               par(mfrow=c(2,1))
           par(oma=c(4,4,2,1),mar=c(0,.5,.5,0),xpd=NA)

           if((is.null(dim(cares))))
               an <- names(cares)
           else
               an <- rownames(cares)


           perarea <- do.call("rbind",strsplit(an,"-"))
           colnames(perarea) <- c("new_periods","new_area")
           cares <- cbind.data.frame(perarea,cares)

           areacol <- palette.colors(ngroup,"Pastel 1",recycle=T)
           if((is.null(dim(ca.res$row$coord)))){
               colnames(cares)[ncol(cares)]="Dim 1"
               dim=1
               plot(1,1,type="n",main="",xlab="",ylab=paste("Dim",dim),xlim=c(1,length(unique(na.omit(cares$new_periods)))),ylim=range(cares[,paste("Dim",dim)]),xaxt="n")
               lapply(1:ngroup,function(area) lines(sort(cares$new_periods[cares$new_area == area]),cares[cares$new_area == area,paste("Dim",dim)][order(cares$new_periods[cares$new_area == area])],col=areacol[area],lwd=3,type="o"))
           }
           else{
               for(dim in 1:2){
                   plot(1,1,type="n",main="",xlab="",ylab=paste("Dim",dim),xlim=c(1,length(unique(na.omit(cares$new_periods)))),ylim=range(cares[,paste("Dim",dim)]),xaxt="n")
                   lapply(1:ngroup,function(area) lines(seq_along(cares$new_periods[cares$new_area == area]),rev(cares[cares$new_area == area,paste("Dim",dim)]),col=areacol[area],lwd=3,type="o"))
               }
           }
           axis(1, at=1+seq_along(unique(cares$new_periods)),label=rev(unique(cares$new_periods)))
}

plotPhase <- function(data){
    gms <- unique(st_drop_geometry(data[,c("GMM","GMS")]))
    gms <- gms[order(gms[,1],decreasing=T),]
    abs_periods <- seq(8500,2500,-500)
    plot(1,1,type="n",xlim=range(-9000,-2000),ylim=range(nrow(gms):1))
    sapply(seq(1,length(abs_periods)-1,2),function(i)rect(-abs_periods[i], -1000, -abs_periods[i + 1], 5500, col = adjustcolor('lightblue',.2), border = NA))
    segments(x0=-gms[,1]-gms[,2],y0=nrow(gms):1,x1=-gms[,1]+gms[,2],y1=nrow(gms):1,col=adjustcolor(1,.3))
    points(-gms[,1],nrow(gms):1,cex=.2)
    abline(v=-abs_periods)
}

plotCAarrows <- function(cares){
    an <- rownames(cares$row$coord)
    gp <- as.numeric(gsub("(\\d+)-(\\d)","\\2",an))
    # Function to calculate a shorter arrow length
    shorten_arrows <- function(gp.coords, shorten_length = 0.03) {
        max_index <- nrow(gp.coords)
        xy0=gp.coords[max_index:2, 1:2]
        xy1=gp.coords[(max_index-1):1, 1:2]
        x0=xy0[,1]
        y0=xy0[,2]
        x1=xy1[,1]
        y1=xy1[,2]

        # Calculate direction vectors
        dx <- x1 - x0
        dy <- y1 - y0
        distance <- sqrt(dx^2 + dy^2)

        # Normalize and scale the direction
        x1_new <- x1 - dx/distance * shorten_length
        y1_new <- y1 - dy/distance * shorten_length
        x0_new <- x0 + dx/distance * shorten_length
        y0_new <- y0 + dy/distance * shorten_length

        return(cbind(x0_new, y0_new,x1_new, y1_new))
    }

    # Prepare your arrows with the shortened lengths

    plot(cares$row$coord[,1:2],col=areapal[gp+1],pch=20,type="n")
    abline(v=0,lty=2,lwd=1.2)#(cts.ca$row$coord[,1:2],col=areapal[gp+1],pch=20)
    abline(h=0,lty=2,lwd=1.2)#(cts.ca$row$coord[,1:2],col=areapal[gp+1],pch=20)
    points(cares$row$coord[,1:2],col=areapal[gp+1],pch=20)
    text(cares$row$coord[,1:2],an,col=areapal[gp+1],pos=3)

    for(g in 1:2){
        gp.coor=cares$row$coord[gp==g,1:2]
        gp.coor <- shorten_arrows(gp.coor,.03)
        arrows( x0=gp.coor[,1],x1=gp.coor[,3],y0=gp.coor[,2],y1=gp.coor[,4],col=areapal[g+1],lw=1,length=.1)
    }

    text(cares$col$coord[,1:2],rownames(cares$col$coord),col="dark green",font=3)
}

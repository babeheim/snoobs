
save_temp <- FALSE

library(RColorBrewer)
library(bbmle)
library(rethinking)
library(tictoc)

set.seed(1912)

scaffold <- FALSE
save_temp <- FALSE

dir_init <- function(path, verbose = FALSE, overwrite = TRUE) {
  if (substr(path, 1, 2) != "./") stop("path argument must be formatted
    with './' at beginning")
  contents <- dir(path, recursive = TRUE)
  if (dir.exists(path)) {
    if (overwrite) {
      if (verbose) {
        if (length(contents) == 0) print(paste("folder ", path, " created.", sep = ""))
        if (length(contents) > 0) print(paste("folder ", path,
          " wiped of ", length(contents), " files/folders.", sep = ""))
      }
      if (dir.exists(path)) unlink(path, recursive = TRUE)
      dir.create(path)
    }
  } else {
    if (verbose) {
      print(paste("folder ", path, " created.", sep = ""))
    }
    dir.create(path)
  }
}

write_log <- function(title, path, start_time) {
  tic.log(format = TRUE)
  msg_log <- unlist(tic.log())
  msg_log <- paste0("- ", msg_log)
  if (!exists("start_time")) start_time <- NA
  header <- c(
    title,
    paste("start_time:", start_time),
    paste("finish_time:", Sys.time()),
    "events:")
  msg_log <- c(header, msg_log)
  writeLines(msg_log, path)
  print("tictoc log written to file")
}

age.pyramid <- function(ages, males, by.unit, sitename=""){

	if(length(ages)==0 | length(males)==0) return()
	
	ages[ages==0] <- 1

	plot.data <- table(ceiling(ages/by.unit)*by.unit, males)
	
	nonzero.groups <- as.numeric(rownames(plot.data))
	all.groups <- seq(by.unit, max(nonzero.groups), by=by.unit)
	
	zeros.data <- matrix(0, nrow=length(all.groups), ncol=2)
	colnames(zeros.data) <- c("0", "1")
	rownames(zeros.data) <- all.groups
	
	zeros.data[all.groups %in% nonzero.groups,] <- plot.data
	
	plot.data <- zeros.data
	
	if(dim(plot.data)[2]==1) plot.data <- cbind(plot.data, rep(0, dim(plot.data)[1]))
	
	# break up the men and women's jobs separately...

	yt <- as.numeric(rownames(plot.data))
	yb <- c((min(yt)-by.unit),yt[-length(yt)])
	xl <- -plot.data[,1]
	xr <- plot.data[,2]

	yt <- seq(by.unit, max(yt), by=by.unit)
	yb <- seq(0, max(yt)-by.unit, by=by.unit)
	
	xspan <- max(abs(min(xl)), max(xr))
	
	# ymax <- max(yt)
	ymax <- 90

	plot(c(-xspan, xspan), c(0, ymax), type="n", yaxt="n", xaxt="n", main=paste(sitename, " age pyramid (n=", sum(!is.na(ages)), ")", sep=""), xlab="Individuals", ylab="", axes=F)

	age.lab <- 1:(ymax/10)*10
	abline(h=age.lab, col="lightblue")
	text( rep(.9*max(xr), length(age.lab)), age.lab, labels=paste("Age", age.lab), pos=3, cex=.7, col="lightblue", offset=.2)

	rect(-xr, yb, rep(0, length(yt)), yt, col="NA", border="blue", lty=2)
	rect(rep(0, length(yt)), yb, -xl, yt, col="NA", border="pink", lty=2)
	rect(xl, yb, rep(0, length(yt)), yt, col="pink", border=NA)
	rect(rep(0, length(yt)), yb, xr, yt, col="blue", border=NA)
	rect(xl, yb, xr, yt, col=NA, border="white")

	tks <- min(xl):max(xr)
	tks[tks%%10==0]
	axis(1, at=tks[tks%%10==0], labels=abs(tks[tks%%10==0]))

}





# need to make sure there's no census 0...

decomp.grapher <- function(data, delta=T, line=FALSE, name="Spectral", cell.outline=1, my.order="float"){

	og.data <- data
	
	
	data <- og.data
	
	phibar <- data[,1]
	data <- data[,-1] # drop first column
	
	n.terms <- ncol(data) 
	n.cen <- nrow(data)
	my.colors <- sample(brewer.pal(n.terms, name))
	my.term.names <- colnames(data)
	

	## Ordering Schemes
	
	
	# Manual
	if(my.order == "manual"){
		my.I.order <- c(2,4,3,1)
		my.D.order <- c(1,3,4,2)
		
		positive.data <- data[,my.I.order]
		positive.data[positive.data < 0] <- 0
		negative.data <- data[,my.D.order]
		negative.data[negative.data > 0] <- 0
		
		color.matrix.I <- data
		color.matrix.D <- data
		for(i in 1:n.cen){
			color.matrix.I[i,] <- my.colors[my.I.order]
			color.matrix.D[i,] <- my.colors[my.D.order]
		}
	}
	
	
	
	# floating order scheme - order bars largest to smallest magnitude going out
	if(my.order == "float"){
		orders <- t(apply(data, 1, function(z) order(z, decreasing=T)))
		
		# modal.order <- apply(orders, 2, modal)
		# orders <- matrix(rep(modal.order, nrow(orders)), ncol=4, nrow=nrow(orders), byrow=T)
		
		data.ordered <- matrix(NA, nrow=n.cen, ncol=n.terms)
		for(i in 1:n.cen){
			current.values <- data[i,]
			data.ordered[i,] <- current.values[orders[i,]]
		}  # the data to plot is now organized in biggest to smallest, going left to right...but the columns no longer correspond to particular terms...only the 

		color.matrix.I <- orders
		for(i in 1:n.cen){
			color.matrix.I[i,] <- my.colors[orders[i,]]
		}
		
		positive.data <- data.ordered
		positive.data[positive.data < 0] <- 0
		
		negative.data <- data.ordered[,n.terms:1]
		negative.data[negative.data > 0] <- 0
		
		color.matrix.D <- color.matrix.I[,n.terms:1]
	}
	
	
	
	
	
	
	
	
	# according to the ordering scheme, populate the y.starts matrix
	
	censuses <- as.numeric(rownames(data))
	x.starts <- censuses - 0.5
	x.stops <- censuses + 0.5
	
	y.stops.I <- t(apply(positive.data, 1, cumsum))
	y.stops.I <- cbind(rep(0, n.cen), y.stops.I)
	
	y.stops.D <- t(apply(negative.data, 1, cumsum))
	y.stops.D <- cbind(rep(0, n.cen), y.stops.D)
	

	# right now, for each census, it decides of the four terms which has the largest mag, then the next largest, then the next, and calls their values from data.ordered in that order.  It does this for both positive and negative.  What I WANT it to do is always use the same order for positive and negative mags...

	
	delta.phibar <- diff(phibar)
	
	dev.new()
	par(mfrow=c(2,2))

	# First Graph: Frequency over Time
	
	scalar <- 7
	
	plot(censuses, phibar, type="l", las=1, ylim=c(0, 1), yaxt="n", ylab="", xlab="year", col="black", main="frequency of trait among individuals", frame.plot=FALSE, xaxt="n")
	axis(2, at=c(seq(0,1, by=0.125)), tcl=-.25, labels=FALSE)
	axis(2, at=c(seq(0,1, by=0.25)), las=1)
	axis(1, at=c(seq(0, n.cen, by=5)), label=c(seq(0, n.cen, by=5)*5))
	axis(2, at=c(0.5+0.01, 0.5+0.01*scalar), las=1)
	
	
	cell.outline=1
	
	y.stops.I2 <- scalar*y.stops.I
	y.stops.D2 <- scalar*y.stops.D
	
	for(i in 1:n.terms){
		rect(xleft=x.starts, xright=x.stops, ybottom=phibar+c(y.stops.I2[-1,i], 0), ytop=phibar+c(y.stops.I2[-1,i+1], 0), lty=cell.outline, col=color.matrix.I[,i], border="gray4")
		rect(xleft=x.starts, xright=x.stops, ybottom=phibar+c(y.stops.D2[-1,i], 0), ytop=phibar+c(y.stops.D2[-1,i+1], 0), lty=cell.outline, col=color.matrix.D[,i], border="gray4")
	}	

	points(censuses, phibar, type="p", pch=20, cex=1)
	
	text(locator(1), label="Bar heights magnified by a factor of 7")
	
	
	o <- order(colMeans(data, na.rm=T), decreasing=T)

	ordered.names <- my.term.names[o]
	ordered.colors <- my.colors[o]
	
	textbox.col <- "black"
	if(cell.outline==0) textbox.col <- "white"
	
	for(i in 1:n.terms){
		loc <- locator(1)
		if(i ==1) record <- loc
		my.pos <- xy.coords(loc$x, record$y)
		points(my.pos, pch=22, bg=ordered.colors[i], col=textbox.col, cex=1.5)
		text(my.pos, labels=ordered.names[i], pos=4)
	}
	

	
	
	
	
	
	
	# Second Graph: Delta Decomposition
	
	cell.outline=1
	
	plot(as.numeric(names(delta.phibar)), delta.phibar, type="l", las=1, ylab="change in trait frq.", xlab="year", col="white", main="frequency delta decomposition", frame.plot=FALSE, ylim=c(min(y.stops.D, na.rm=T), max(y.stops.I, na.rm=T)), xaxt="n")
	axis(1, at=c(seq(0, n.cen, by=5)), label=c(seq(0, n.cen, by=5)*5))

	for(i in 1:n.terms){
		rect(xleft=x.starts, xright=x.stops, ybottom=y.stops.I[,i], ytop=y.stops.I[,i+1], lty=cell.outline, col=color.matrix.I[,i], border="gray4")
		rect(xleft=x.starts, xright=x.stops, ybottom=y.stops.D[,i], ytop=y.stops.D[,i+1], lty=cell.outline, col=color.matrix.D[,i], border="gray4")
	}	


	o <- order(colMeans(data, na.rm=T), decreasing=T)

	ordered.names <- my.term.names[o]
	ordered.colors <- my.colors[o]
	
	textbox.col <- "black"
	if(cell.outline==0) textbox.col <- "white"
	
	for(i in 1:n.terms){
		loc <- locator(1)
		if(i ==1) record <- loc
		my.pos <- xy.coords(loc$x, record$y)
		points(my.pos, pch=22, bg=ordered.colors[i], col=textbox.col, cex=1.5)
		text(my.pos, labels=ordered.names[i], pos=4)
	}
	

	
	# Third Plot
	
	o <- order(colMeans(data, na.rm=T), decreasing=T)

	data.plot <- data[,o]
	
	# boxplot.border.colors <- ifelse(cell.outline==1, rep("black", length(my.colors)), my.colors[o])
	boxplot.border.colors <- my.colors[o]
	if(cell.outline==1) boxplot.border.colors <- "black"
	my.box <- boxplot(data.plot, col=my.colors[o], border=boxplot.border.colors, pch=20, main="average effect on trait frequency", xaxt="n", xlab="", las=1, frame.plot=FALSE, cex=1, ylim=c(lb1, ub2), lty=1)
# 	boxplot(data.plot, col=my.colors[o], pch=20, main="Average Effect on Trait Frequency", xaxt="n", xlab="", border=my.colors[o], las=1)
	
	medians <- my.box$stats[3,]
	
	abline(h=0, lty=2, col=gray(.35))
	boxplot(data.plot, col=my.colors[o], pch=20, main="average effect on trait frequency", add=T,  xlab="", border=boxplot.border.colors, las=1, frame.plot=FALSE, xaxt="n", cex=1, ylim=c(lb1, ub2), lty=1)
	axis(1, at=c(1:n.terms), labels=structure.names[o], tick=FALSE)
	
	width <- 0.4
	x.box.starts <- 1:n.terms - width
	x.box.stops <- 1:n.terms + width
	
	ordered.colors <- my.colors[o]
	
	for(i in 1:n.terms){
	  points(c(x.box.starts[i], x.box.stops[i]), c(medians[i], medians[i]), type="l", col=ordered.colors[i], lwd=3)
	}
	
	
	means <- colMeans(data.plot, na.rm=T)
	
	for(i in 1:n.terms){
		points(c(x.box.starts[i], x.box.stops[i]), c(means[i], means[i]), type="l", col="black", lwd=2)
	}
	
	
	
	
	
	# Fourth Plot
	
	
	
	
	force.percent.data <- abs(data)/rowSums(abs(data), na.rm=T)
	
	data <- force.percent.data
	
	orders <- t(apply(data, 1, function(z) order(z, decreasing=T)))
	
	o <- order(colMeans(force.percent.data, na.rm=T), decreasing=T)
	
	
	for(i in 1:dim(orders)[1]){
		orders[i,] <- o
	}
	
	data.plot <- matrix(NA, nrow=dim(data)[1], ncol=dim(data)[2])
	for(i in 1:dim(orders)[1]){
		current.values <- data[i,]
		data.plot[i,] <- current.values[orders[i,]]
	}
	
	color.matrix <- orders
	for(i in 1:dim(orders)[1]){
		color.matrix[i,] <- my.colors[orders[i,]]
	}
	

	index <- 2:((n.terms)*2)
	for(i in index){
		pos <- which(data.plot[,my.colnames[i]] >= 0)
		if(i %% 2 == 0){
			if(length(pos) > 0){
				y.starts[pos,i:max(index)] <- y.starts[pos,i-1]+data.plot[pos,my.colnames[i]]
			}
		}
	}
	
 
	plot(censuses, graph.data, type="l", las=1, ylim=c(0,1), ylab="%", xlab="census", col="white", main="relative size of force categories", yaxt="n", frame.plot=FALSE)
	ats <- seq(0, 1, by=0.2)
	axis(2, at=ats, labels=100*ats, las=1) 
	# axis(1, at=c(censuses))
	index <- 1:(2*n.terms)
	index <- index[index%%2==1]
	for(i in index){
		rect(xleft=x.starts, xright=x.stops, ybottom=y.starts[,i], ytop=y.starts[,i+1], col=color.matrix[,my.colnames[i]], lty=cell.outline)
# 		rect(xleft=x.starts, xright=x.stops, ybottom=float+y.starts[,(i+2*n.terms)], ytop=float+y.starts[,(i+2*n.terms+1)], col=color.matrix[,my.colnames[max.index+1-(i+2*n.terms)]], lty=0)
	}	



}





































decomp.grapher <- function(data, delta=T, line=FALSE, name="Spectral", cell.outline=1){

	censuses <- as.numeric(rownames(data))
	p <- data[,1]
	
	x.starts <- censuses - 0.5
	x.stops <- censuses + 0.5

	n.parts <- dim(data)[2] - 1
	
	n.I <- n.parts
	n.D <- n.parts

	y.starts <- matrix(0, nrow=length(censuses), ncol=(n.parts*4))
	rownames(y.starts) <- censuses

	my.colnames <- rep(sort(rep(1:n.parts, 2)), 2)
	colnames(y.starts) <- my.colnames
	
	data.red <- data[,-1] # sort it how you want
		
	structure.names <- colnames(data.red)
	colnames(data.red) <- structure.names
	
	orders <- t(apply(data.red, 1, function(z) order(z, decreasing=T)))
	
	# modal.order <- apply(orders, 2, modal)
	# orders <- matrix(rep(modal.order, nrow(orders)), ncol=4, nrow=nrow(orders), byrow=T)
	
	data.plot <- matrix(NA, nrow=dim(data.red)[1], ncol=dim(data.red)[2])
	for(i in 1:dim(orders)[1]){
		current.values <- data.red[i,]
		data.plot[i,] <- current.values[orders[i,]]
	}
	# the data to plot is now organized in biggest to smallest, going left to right...but the columns no longer correspond to particular terms
	
	# color.vector <- sample(c("blue", "darkblue", "darkred", "red", "cadetblue", "chartreuse2", "deeppink4", "goldenrod1"))
	# my.colors <- color.vector[1:n.parts]
	
	my.colors <- sample(brewer.pal(n.parts, name))

	
	color.matrix <- orders
	for(i in 1:dim(orders)[1]){
		color.matrix[i,] <- my.colors[orders[i,]]
	
	}
	

	index <- 2:((n.parts)*2)
	for(i in index){
		pos <- which(data.plot[,my.colnames[i]] > 0)
		if(i %% 2 == 0){
			if(length(pos) > 0){
				y.starts[pos,i:max(index)] <- y.starts[pos,i-1]+data.plot[pos,my.colnames[i]]
			}
		}
	}
	
	index <- 2*n.parts+2:((n.parts)*2)
	max.index <- max(index)
	for(i in index){
		neg <- which(data.plot[,my.colnames[max.index+1-i]] < 0)
		if(i %% 2 == 0){
			if(length(neg) > 0){
				y.starts[neg,i:max(index)] <- y.starts[neg,i-1]+data.plot[neg,my.colnames[max.index+1-i]]
			}
		}	
	}


	# right now, for each census, it decides of the four terms which has the largest mag, then the next largest, then the next, and calls their values from data.plot in that order.  It does this for both positive and negative.  What I WANT it to do is always use the same order for positive and negative mags...


	
	
	
	
	upper <- max(y.starts[,(1:n.parts)*2])
	lower <- min(y.starts[,(1:n.parts)*2 + 2*n.parts])

	lag <- diff(p)
	graph.data <- c(lag, 0)
	float <- 0
	if(delta==FALSE) float <- p	
	

	censuses.plot <- as.numeric(names(graph.data))

	
	# par(bg="#FCFCFC")  # aka gray(.9) 	 	 
	dev.new()
	par(mfrow=c(2,2))

	lb1 <- min(float)+min(y.starts)
	ub2 <- (max(float)+max(y.starts))
	
	
	# First Graph: Frequency Decomposition
	
	plot(censuses.plot, p, type="l", las=1, ylim=c(0, 1), yaxt="n", ylab="", xlab="census", col="black", main="frequency of trait among individuals", frame.plot=FALSE, xlim=c(min(na.omit(censuses.plot)), max(na.omit(censuses.plot))))
	#abline(h=1, col="gray")
	#abline(h=0, col="gray")
	
	axis(2, at=c(seq(0,1, by=0.125)), tcl=-.25, labels=FALSE)
	axis(2, at=c(seq(0,1, by=0.25)), las=1)
	# axis(1, at=c(censuses.plot))
	
	# Second Graph: Delta Decomposition
	
	plot(censuses.plot, graph.data, type="l", las=1, ylim=c(lb1, ub2*1.3), ylab="change in trait frq.", xlab="census", col="white", main="frequency delta decomposition", frame.plot=FALSE)
# 	axis(1, at=c(censuses.plot))
	index <- 1:(2*n.parts)
	index <- index[index%%2==1]
	for(i in index){
		rect(xleft=x.starts, xright=x.stops, ybottom=float+y.starts[,i], ytop=float+y.starts[,i+1], col=color.matrix[,my.colnames[i]], lty=cell.outline)
		rect(xleft=x.starts, xright=x.stops, ybottom=float+y.starts[,(i+2*n.parts)], ytop=float+y.starts[,(i+2*n.parts+1)], col=color.matrix[,my.colnames[max.index+1-(i+2*n.parts)]], lty=cell.outline)
	}	

	abline(h=0, col="white", lwd=1)

	if(delta==T & line==T) points(rowSums(data.plot, na.rm=T), pch=20, col="white", type="o")
	if(delta==FALSE & line==T) points(p, pch=20, col="white", type="o")

	
	
	
	o <- order(colMeans(data.red, na.rm=T), decreasing=T)

	ordered.names <- structure.names[o]
	ordered.colors <- my.colors[o]
	
	textbox.col <- "black"
	if(cell.outline==0) textbox.col <- "white"
	
	for(i in 1:n.parts){
		loc <- locator(1)
		if(i ==1) record <- loc
		my.pos <- xy.coords(loc$x, record$y)
		points(my.pos, pch=22, bg=ordered.colors[i], col=textbox.col, cex=1.5)
		text(my.pos, labels=ordered.names[i], pos=4)
	}
	

#	points(censuses.plot, graph.data, type="l", col="white", lwd=1)

	# ub <- max(data[,2:(n.parts+1)], na.rm=T)
	# lb <- min(data[,2:(n.parts+1)], na.rm=T)

	# plot(censuses.plot, data[,2], type="l", col=my.colors[1], ylim=c(lb, ub), las=1, ylab="Change in Trait Frq.", xlab="", xlim=c(min(na.omit(censuses.plot), na.rm=T), 2+max(censuses.plot, na.rm=T)))
	# abline(h=0)
	# points(censuses.plot, data[,3], type="l", col=my.colors[2])
	# points(censuses.plot, data[,4], type="l", col=my.colors[3])
	# points(censuses.plot, data[,5], type="l", col=my.colors[4])
	# points(censuses.plot, data[,6], type="l", col=my.colors[5])
	# points(censuses.plot, data[,7], type="l", col=my.colors[6])
	# points(censuses.plot, data[,8], type="l", col=my.colors[7])

	# means <- colMeans(data, na.rm=T)
	# means <- means[-1]
	# axis(4, at=means, labels=NA, las=1, tck=0.03)
	
	# points(rep(2+max(censuses.plot, na.rm=T), length(means)), means, pch=22, bg=my.colors)
	

	
	
	data.red <- data[,-1]
	colnames(data.red) <- structure.names
	
	o <- order(colMeans(data.red, na.rm=T), decreasing=T)

	data.red.plot <- data.red[,o]
	
	# boxplot.border.colors <- ifelse(cell.outline==1, rep("black", length(my.colors)), my.colors[o])
	boxplot.border.colors <- my.colors[o]
	if(cell.outline==1) boxplot.border.colors <- "black"
	my.box <- boxplot(data.red.plot, col=my.colors[o], border=boxplot.border.colors, pch=20, main="average effect on trait frequency", xaxt="n", xlab="", las=1, frame.plot=FALSE, cex=1, ylim=c(lb1, ub2), lty=1)
# 	boxplot(data.red.plot, col=my.colors[o], pch=20, main="Average Effect on Trait Frequency", xaxt="n", xlab="", border=my.colors[o], las=1)
	
	medians <- my.box$stats[3,]
	
	abline(h=0, lty=2, col=gray(.35))
	boxplot(data.red.plot, col=my.colors[o], pch=20, main="average effect on trait frequency", add=T,  xlab="", border=boxplot.border.colors, las=1, frame.plot=FALSE, xaxt="n", cex=1, ylim=c(lb1, ub2), lty=1)
	axis(1, at=c(1:n.parts), labels=structure.names[o], tick=FALSE)
	
	width <- 0.4
	x.box.starts <- 1:n.parts - width
	x.box.stops <- 1:n.parts + width
	
	ordered.colors <- my.colors[o]
	
	for(i in 1:n.parts){
	  points(c(x.box.starts[i], x.box.stops[i]), c(medians[i], medians[i]), type="l", col=ordered.colors[i], lwd=3)
	}
	
	
	means <- colMeans(data.red.plot, na.rm=T)
	
	for(i in 1:n.parts){
		points(c(x.box.starts[i], x.box.stops[i]), c(means[i], means[i]), type="l", col="black", lwd=2)
	}
	
	
	
	
	
	
	
	
	
	
	
	data.red <- data[,-1]
	
	force.percent.data <- abs(data.red)/rowSums(abs(data.red), na.rm=T)
	
	data.red <- force.percent.data
	
	orders <- t(apply(data.red, 1, function(z) order(z, decreasing=T)))
	
	o <- order(colMeans(force.percent.data, na.rm=T), decreasing=T)
	
	
	for(i in 1:dim(orders)[1]){
		orders[i,] <- o
	}
	
	data.plot <- matrix(NA, nrow=dim(data.red)[1], ncol=dim(data.red)[2])
	for(i in 1:dim(orders)[1]){
		current.values <- data.red[i,]
		data.plot[i,] <- current.values[orders[i,]]
	}
	
	color.matrix <- orders
	for(i in 1:dim(orders)[1]){
		color.matrix[i,] <- my.colors[orders[i,]]
	}
	

	index <- 2:((n.parts)*2)
	for(i in index){
		pos <- which(data.plot[,my.colnames[i]] >= 0)
		if(i %% 2 == 0){
			if(length(pos) > 0){
				y.starts[pos,i:max(index)] <- y.starts[pos,i-1]+data.plot[pos,my.colnames[i]]
			}
		}
	}
	
 
	plot(censuses.plot, graph.data, type="l", las=1, ylim=c(0,1), ylab="%", xlab="census", col="white", main="relative size of force categories", yaxt="n", frame.plot=FALSE)
	ats <- seq(0, 1, by=0.2)
	axis(2, at=ats, labels=100*ats, las=1) 
	# axis(1, at=c(censuses.plot))
	index <- 1:(2*n.parts)
	index <- index[index%%2==1]
	for(i in index){
		rect(xleft=x.starts, xright=x.stops, ybottom=y.starts[,i], ytop=y.starts[,i+1], col=color.matrix[,my.colnames[i]], lty=cell.outline)
# 		rect(xleft=x.starts, xright=x.stops, ybottom=float+y.starts[,(i+2*n.parts)], ytop=float+y.starts[,(i+2*n.parts+1)], col=color.matrix[,my.colnames[max.index+1-(i+2*n.parts)]], lty=0)
	}	



}















# need to include an age cutoff for the new one




decomposer <- function(data, check=F, drop.empty=T){

	# if grouping is not NA, then decomposer needs to modify the states to reflect transfers!  

	grouping <- NA
	
	if("group.id" %in% colnames(data)) grouping <- "group.id"
	
	full.census.data <- as.numeric(data[,"census"])
	full.phenotype.data <- as.numeric(data[,"phi"])
	full.census.list <- sort(unique(full.census.data))
	census.phenotype.frqs <- mean(full.phenotype.data[full.census.data==full.census.list[1]])
	n.censuses <- length(full.census.list)
	full.state.data <- as.numeric(data[,"state"])

	output.table <- matrix(NA, ncol=6+3, nrow=n.censuses)
	colnames(output.table) <- c("p", "rho", "emig", "immig", "RS", "delta", "mortality", "predicted", "actual")

	
	
	if(!is.na(grouping)){
	
		full.group.data <- as.numeric(data[,grouping])
		full.group.list <- sort(unique(full.group.data))
		n.groups <- length(full.group.list)
		
		group.table <- matrix(NA, ncol=n.groups, nrow=17)
		rownames(group.table) <- c("lp.n", "tp.n", "lp.phi.bar", "tp.phi.bar", "Delta.phi.bar", "lp.c", "tp.c", "Delta.c", "G", "i", "e", "b", "d", "a", "rho.bar.a", "r", "rho.bar")
		colnames(group.table) <- full.group.list
		
		output.table <- matrix(NA, ncol=11+3, nrow=n.censuses)
		colnames(output.table) <- c("p", "r.s", "r.t", "I", "i", "E", "e", "RS", "rs", "D", "d", "de", "predicted", "actual")

	}
		
		
	for(t in 2:n.censuses){

		lp.rows <- which(full.census.data==full.census.list[t-1])
		lp.phi.data <- as.numeric(data[lp.rows,"phi"])
		lp.state <- as.numeric(data[lp.rows, "state"])
		
		lp.active <- which(lp.state!=2 & lp.state!=3)
		lp.n <- length(na.omit(lp.phi.data[lp.active]))

		tp.rows <- which(full.census.data==full.census.list[t])
		tp.phi.data <- as.numeric(data[tp.rows,"phi"])
		tp.state <- as.numeric(data[tp.rows, "state"])
		
		tp.active <- which(tp.state!=2 & tp.state!=3)
		tp.n <- length(na.omit(tp.phi.data[tp.active]))
		
		lp.phi.bar <- mean(lp.phi.data[lp.active])
		tp.phi.bar <- mean(tp.phi.data[tp.active])
		census.phenotype.frqs <- c(census.phenotype.frqs, tp.phi.bar)

		Delta.phi.bar <- tp.phi.bar - lp.phi.bar
		G <- tp.n/lp.n

		r <- sum(tp.state==0 | tp.state==1)/lp.n
		tp.rhos <- as.numeric(data[tp.rows,"i.change"])
		tp.rho.bar <- mean(tp.rhos[tp.active], na.rm=T)
		
		e <- sum(tp.state==2)/lp.n
		phi.bar.e <- mean(tp.phi.data[tp.state==2])
		d <- sum(tp.state==3)/lp.n
		phi.bar.d <- mean(tp.phi.data[tp.state==3])
		i <- sum(tp.state==4)/lp.n
		phi.bar.i <- mean(tp.phi.data[tp.state==4])
		b <- sum(tp.state==5)/lp.n
		phi.bar.b <- mean(tp.phi.data[tp.state==5])
		lp.b.data <- as.numeric(data[lp.rows,"n.kids"])/2
		phi.bar.bu <- sum(lp.b.data*lp.phi.data)/sum(tp.state==5)
		
		tp.m.deltas <- as.numeric(data[tp.rows, "m.delta"])
		tp.f.deltas <- as.numeric(data[tp.rows, "f.delta"])
		tp.delta.bar <- sum(c(tp.m.deltas, tp.f.deltas), na.rm=T)/sum(tp.state==5)
			
			
		# single-level BDICE
		immigration <- i*(phi.bar.i - lp.phi.bar)
		emigration <- e*(phi.bar.e - lp.phi.bar)*(-1)
		reproductive.success <- b*(phi.bar.bu - lp.phi.bar)
		biased.transmission <- b*tp.delta.bar/2
		mortality <- d*(phi.bar.d - lp.phi.bar)*(-1)
		individual.change <- r*tp.rho.bar
		
		predicted <- sum(immigration, emigration, reproductive.success, biased.transmission, mortality, individual.change, na.rm=T)*(1/G)
		
		if(is.na(grouping)){
			output.table[t,] <- c(0, individual.change, emigration, immigration, reproductive.success, biased.transmission, mortality, predicted, Delta.phi.bar)
		}
		# c("i.change", "emig", "immig", "RS", "trans.bias", "mortality", "predicted", "actual")
		
		if(!is.na(grouping)){
				
			for(j in 1:n.groups){
						
				glp.rows <- which(full.census.data==full.census.list[t-1] & full.group.data==full.group.list[j])
				gtp.rows <- which(full.census.data==full.census.list[t] & full.group.data==full.group.list[j])
				glp.state <- as.numeric(data[glp.rows, "state"])
				gtp.state <- as.numeric(data[gtp.rows, "state"])
				
				glp.active <- which(glp.state!=2 & glp.state!=3)
				gtp.active <- which(gtp.state!=2 & gtp.state!=3)
				
				glp.phi.data <- as.numeric(data[glp.rows,"phi"])
				gtp.phi.data <- as.numeric(data[gtp.rows,"phi"])
				
				glp.n <- length(na.omit(glp.phi.data[glp.active]))
				gtp.n <- length(na.omit(gtp.phi.data[gtp.active]))
				g.G <- gtp.n/glp.n
				
				glp.phi.bar <- mean(glp.phi.data[glp.active], na.rm=T)		
				gtp.phi.bar <- mean(gtp.phi.data[gtp.active], na.rm=T)
				g.Delta.phi.bar <- gtp.phi.bar - glp.phi.bar
						
				glp.c <- glp.n/lp.n
				gtp.c <- gtp.n/tp.n
				Delta.c <- gtp.c - glp.c
				
	# 0 - remaining, no group transfer
	# 1 - remaining, transfered groups
	# 2 - emigrated 
	# 3 - died 
	# 4 - immigrated 
	# 5 - born 
		
				g.r <- sum(gtp.state==0)/glp.n
				g.a <- sum(gtp.state==1)/glp.n
		
				gtp.rhos <- as.numeric(data[gtp.rows,"i.change"])
				gtp.rho.bar.s <- mean(gtp.rhos[gtp.state==0], na.rm=T)
				gtp.rho.bar.a <- mean(gtp.rhos[gtp.state==1], na.rm=T)
				
				g.e <- sum(gtp.state==2)/glp.n
				g.d <- sum(gtp.state==3)/glp.n
				g.i <- sum(gtp.state==4)/glp.n
				g.b <- sum(gtp.state==5)/glp.n
					
				group.table[,j] <- c(glp.n, gtp.n, glp.phi.bar, gtp.phi.bar, g.Delta.phi.bar, glp.c, gtp.c, Delta.c, g.G, g.i, g.e, g.b, g.d, g.a, gtp.rho.bar.a, g.r, gtp.rho.bar.s)
				# c("lp.n", "tp.n", "lp.phi.bar", "tp.phi.bar", "Delta.phi.bar", "lp.c", "tp.c", "Delta.c", "G", "i", "e", "b", "d", "a", "rho.bar.a", "r", "rho.bar")
			}
						
			# The Baldini Equation
			# sum(group.table["Delta.c",]*group.table["lp.phi.bar",]) + sum(group.table["tp.c",]*group.table["Delta.phi.bar",])
					
			# multilevel, my version
			sessile.i.change <- sum(group.table["r",]*group.table["rho.bar",]*group.table["lp.c",], na.rm=T)
			transfer.i.change <- sum(group.table["a",]*group.table["rho.bar.a",]*group.table["lp.c",], na.rm=T)
			
			cow.immig <- sum(group.table["i",]*group.table["lp.phi.bar",]*group.table["lp.c",], na.rm=T) - i*lp.phi.bar
			ew.immig <- i*phi.bar.i - sum(group.table["i",]*group.table["lp.phi.bar",]*group.table["lp.c",], na.rm=T)
			
			cow.emig <- sum(group.table["e",]*group.table["lp.phi.bar",]*group.table["lp.c",], na.rm=T) - e*lp.phi.bar
			ew.emig <- e*phi.bar.e - sum(group.table["e",]*group.table["lp.phi.bar",]*group.table["lp.c",], na.rm=T)
			
			cow.birth <- sum(group.table["b",]*group.table["lp.phi.bar",]*group.table["lp.c",], na.rm=T) - b*lp.phi.bar
			ew.birth <- b*phi.bar.bu - sum(group.table["b",]*group.table["lp.phi.bar",]*group.table["lp.c",], na.rm=T)
			
			cow.death <- sum(group.table["d",]*group.table["lp.phi.bar",]*group.table["lp.c",], na.rm=T) - d*lp.phi.bar
			ew.death <- d*phi.bar.d - sum(group.table["d",]*group.table["lp.phi.bar",]*group.table["lp.c",], na.rm=T)
				
			predicted <- sum(sessile.i.change, transfer.i.change, cow.immig, ew.immig, -1*cow.emig, -1*ew.emig, cow.birth, ew.birth, -1*cow.death, -1*ew.death, biased.transmission, na.rm=T)*(1/G)	
			
			output.table[t,] <- c(0, sessile.i.change, transfer.i.change, cow.immig, ew.immig, -1*cow.emig, -1*ew.emig, cow.birth, ew.birth, -1*cow.death, -1*ew.death, biased.transmission, predicted, Delta.phi.bar)	
			# c("p", "s.change", "t.change", "cow.i", "ew.i", "cow.e", "ew.e", "cow.b", "ew.b", "cow.d", "ew.d", "trans.bias", "predicted", "actual")
		}
		
		print(full.census.list[t])
		
	}
			
	output.table[,"p"] <- census.phenotype.frqs
	
	if(check==F){
		cols <- dim(output.table)[2]
		output.table <- output.table[,-c(cols-1, cols)]
	}
	
	rownames(output.table) <- full.census.list
	if(any(!apply(output.table, 2, function(z) any(!is.na(z)))) & drop.empty){
		output.table <- output.table[,-which(!apply(output.table, 2, function(z) any(!is.na(z))))]
	}
	if(any(!apply(output.table, 2, function(z) any(z!=0 & !is.na(z)))) & drop.empty){
		output.table <- output.table[,-which(!apply(output.table, 2, function(z) any(z!=0 & !is.na(z))))]
	}
	
	output.table
	
}












# ah, a problem.  I need to classify everyone who was "born" but whose parents were not in the population at the last census as an immigrant!
# another problem: individuals who in the population as deaths or emigrants but were not present in the previous time census should be dropped

BDICE.sim.decomp.prepper <- function(data, phenotype="snoob", grouping=NA){
			
	# 00. age cutoffs have to go in here...
	# need to ID a group variable here, and then DON'T specify one in the actual decomposer (the prepper will be the only function to call a group)

		# needs to turn the simulation data into the prepped data needed for decomposer

		# BDICE sim: 
		# census id state mom dad age last.im last.em died male mate counter zygotes h.gene height snoob language t3 
		
		# desired variables: 
		# census id state mid fid m.delta f.delta i.change n.kids phenotype   and a grouping variable
		
	# 0. specify a variable to be called "phenotype", 
		
	# 1. state has to be translated: 
		
	# State Codes:
	# Dead - 0
	# Active - 1
	# Pregnant - 2
	# Emigrant - 11

	# 0 - remaining, no group transfer
	# 1 - remaining, transfered groups
	# 2 - emigrated 
	# 3 - died 
	# 4 - immigrated 
	# 5 - born 
	
	age.cutoff <- 0
	drop.rows <- which(data[,"age"] < age.cutoff*365)
	if(length(drop.rows)>0) data <- data[-drop.rows,]
	
	id.data <- data[,"id"]
	id.list <- sort(unique(id.data))
	census.data <- data[,"census"]
	census.list <- sort(unique(census.data))
	state.data <- data[,"state"]
		
	ghost.rows.to.drop <- integer(0)
	for(j in 2:length(census.list)){
		this.census.names <- id.data[census.data == census.list[j]]
		this.census.states <- state.data[census.data == census.list[j]]
		this.census.outfluxers <- this.census.names[this.census.states==11 | this.census.states==0]
		last.census.names <- id.data[census.data == census.list[j-1]]
		ghosts <- this.census.outfluxers[!this.census.outfluxers %in% last.census.names]
		if(length(ghosts)>0){
			ghost.rows.to.drop <- c(ghost.rows.to.drop, which(id.data %in% ghosts & census.data == census.list[j]))
		}
	}
	
	if(length(ghost.rows.to.drop) > 0) data <- data[-ghost.rows.to.drop,]
	
	id.data <- data[,"id"]
	id.list <- sort(unique(id.data))
	census.data <- data[,"census"]
	state.data <- data[,"state"]
	
	
	
	
	
	new.state.data <- state.data
	new.state.data[state.data==1] <- 0
	
	mentioned.matrix <- table(id.data, census.data)
	census.list <- as.numeric(colnames(mentioned.matrix))
	first.mention.list <- apply(mentioned.matrix, 1, function(z) min(as.numeric(names(z[z==1]))))

	
	if(!is.na(grouping)){
		
		grouping.data <- data[,grouping]
			
		for(j in 2:length(census.list)){
		
			this.census.names <- id.data[census.data == census.list[j]]
			last.census.names <- id.data[census.data == census.list[j-1]]
			this.census.groups <- grouping.data[census.data == census.list[j]]
			last.census.groups.data <- grouping.data[census.data == census.list[j-1]]
			
			last.census.groups <- last.census.groups.data[match(this.census.names, last.census.names)]

			insert.state.vec <- rep(0, length(this.census.names))
			
			insert.state.vec[which(last.census.groups != this.census.groups)] <- 1
			
			new.state.data[census.data==census.list[j]] <- insert.state.vec

#			if(j %% 10 == 0) print(j)
		}
		
	}
	
	
	
	for(j in 1:length(census.list)){
		
		new.arrivals <- id.list[first.mention.list == census.list[j]]
		
		data.subset <- data[census.data==census.list[j],]
		
		immigrated.vec <- data.subset[match(new.arrivals, data.subset[,"id"]),"mom"] == 0
			
		insert.vec <- immigrated.vec
		insert.vec[immigrated.vec==1] <- 4
		insert.vec[immigrated.vec==0] <- 5
		
		men.and.women.last.census <- data[census.data==census.list[j-1], "id"]
		
		recorded.moms.of.newborns <- data.subset[match(new.arrivals, data.subset[,"id"]), "mom"]
		recorded.dads.of.newborns <- data.subset[match(new.arrivals, data.subset[,"id"]), "dad"]
	
		recorded.moms.of.newborns[which(!(recorded.moms.of.newborns %in% men.and.women.last.census))]
		new.arrivals[which(!(recorded.moms.of.newborns %in% men.and.women.last.census))]
		recorded.dads.of.newborns[which(!(recorded.dads.of.newborns %in% men.and.women.last.census))]
		new.arrivals[which(!(recorded.dads.of.newborns %in% men.and.women.last.census))]
	
		if(any(!(recorded.moms.of.newborns %in% men.and.women.last.census))){
			insert.vec[which(!(recorded.moms.of.newborns %in% men.and.women.last.census))] <- 4
		}

		if(any(!(recorded.dads.of.newborns %in% men.and.women.last.census))){
			insert.vec[which(!(recorded.dads.of.newborns %in% men.and.women.last.census))] <- 4
		}
				
		new.state.subvec <- new.state.data[census.data==census.list[j]]
		
		new.state.subvec[match(new.arrivals, id.data[census.data==census.list[j]])] <- insert.vec
		
		new.state.data[census.data==census.list[j]] <- new.state.subvec 
				
#		if(j %% 10 == 0) print(j)
	}

	# in case their first census was also their last census, we need to have this done after the birth/immigrant coding
	new.state.data[state.data==0] <- 3
	new.state.data[state.data==11] <- 2 

	# New state codes:
	# 0 - remaining, no group transfer
	# 1 - remaining, transfered...only works if grouping is not NA
	# 2 - emigrated 
	# 3 - died 
	# 4 - immigrated - easy enough...first time they show up find out if they were born or not
	# 5 - born 

	m.delta.data <- rep(0, length(id.data))
	f.delta.data <- m.delta.data
	phenotype.data <- data[,as.character(phenotype)]
	n.kids.data <- rep(0, length(census.data))

	# i should probably re-code the "new state data" so that i'm not relying on the state codes from the simulation...this could prove problematic in the long run...

	# 3. i.change also to be calculated for second census and so forth

	i.change.data <- rep(NA, length(id.data))

	for(j in 2:length(census.list)){
		
		census.remaining.rows <- census.data==census.list[j] & (new.state.data==0 | new.state.data == 1)
		census.remainers <- id.data[census.remaining.rows]
		
		remainers.current.phenotype <- phenotype.data[ id.data %in% census.remainers & census.data == census.list[j] ] 
		remainers.old.phenotype <- phenotype.data[ id.data %in% census.remainers & census.data == census.list[j-1] ] 
				
		individual.change <- remainers.current.phenotype - remainers.old.phenotype
		
		i.change.data[id.data %in% census.remainers & census.data == census.list[j]] <- individual.change
		
#		if(j %% 10 == 0) print(j)
		
	}

	# he he 
	# barplot(tapply(i.change.data, c(census.data), function(z) sum(z, na.rm=T)))

	# 2. m.delta and f.delta have to be calculated for each child's census of birth, but NA all other times.  
	# 4. n.kids calculated as # of kids born to that person THAT time census, not cumulative

	for(j in 1:(length(census.list)-1)){
			
		next.census.birth.rows <- census.data==census.list[j+1] & new.state.data==5
		next.census.births <- id.data[next.census.birth.rows]
		next.census.moms.data <- data[next.census.birth.rows,"mom"]
		next.census.dads.data <- data[next.census.birth.rows,"dad"]
		next.census.moms <- sort(unique(data[next.census.birth.rows,"mom"]))
		next.census.dads <- sort(unique(data[next.census.birth.rows,"dad"]))
		
		mom.n.kids <- sapply(next.census.moms, function(z) sum(next.census.moms.data==z))
		dad.n.kids <- sapply(next.census.dads, function(z) sum(next.census.dads.data==z))
		
		n.kids.data[census.data==census.list[j] & id.data %in% next.census.moms] <- mom.n.kids
		n.kids.data[census.data==census.list[j] & id.data %in% next.census.dads] <- dad.n.kids	
			
			
		kid.phenotypes <- phenotype.data[next.census.birth.rows]
			
		mom.phenotypes <- phenotype.data[census.data==census.list[j] & id.data %in% next.census.moms]
		names(mom.phenotypes) <- next.census.moms
		
		dad.phenotypes <- phenotype.data[census.data==census.list[j] & id.data %in% next.census.dads]
		names(dad.phenotypes) <- next.census.dads
		
		m.deltas <- kid.phenotypes - mom.phenotypes[as.character(next.census.moms.data)]
		f.deltas <- kid.phenotypes - dad.phenotypes[as.character(next.census.dads.data)]
		
		m.delta.data[next.census.birth.rows] <- m.deltas
		f.delta.data[next.census.birth.rows] <- f.deltas
		
#		if(j %% 10 == 0) print(j)
		
	}


	if(is.na(grouping)){
		prepped.data <- cbind(census.data, id.data, new.state.data, data[,"mom"], m.delta.data, data[,"dad"], f.delta.data, i.change.data, n.kids.data, phenotype.data)
		colnames(prepped.data) <- c("census", "id", "state", "mid", "m.delta", "fid", "f.delta", "i.change", "n.kids", "phi")
	} else {
		prepped.data <- cbind(census.data, id.data, new.state.data, data[,"mom"], m.delta.data, data[,"dad"], f.delta.data, i.change.data, n.kids.data, phenotype.data, grouping.data)
		colnames(prepped.data) <- c("census", "id", "state", "mid", "m.delta", "fid", "f.delta", "i.change", "n.kids", "phi", "group.id")	
	}
	
	
	# FIgured out the bug I need to fix: for each time census, identify those who are now dead or emigrated, and make sure the phenotype, group, etc. for them is actually the 
	# data from the previous time census; strictly speaking we do not know the values recorded at the time of death, and they cannot appear in the calculations!

	
	final.id.data <- unlist(prepped.data[,"id"])
	final.state.data <- unlist(prepped.data[,"state"])
	final.phi.data <- unlist(prepped.data[,"phi"])
	if(!is.na(grouping)) final.group.data <- unlist(prepped.data[,"group.id"])
	outflux.ids <- final.id.data[final.state.data==3 | final.state.data==2]
	
	
	for(i in 1:length(outflux.ids)){
		this.person.rows <- which(final.id.data==outflux.ids[i])
		if(final.phi.data[this.person.rows[length(this.person.rows)]] != final.phi.data[this.person.rows[length(this.person.rows)-1]]) print(paste(outflux.ids[i], " was fixed.", sep=""))
		prepped.data[this.person.rows[length(this.person.rows)], "phi"] <- final.phi.data[this.person.rows[length(this.person.rows)-1]]
		if(!is.na(grouping)){
			prepped.data[this.person.rows[length(this.person.rows)], "group.id"] <- final.group.data[this.person.rows[length(this.person.rows)-1]]
		}
	}
	
	prepped.data


	

}












# workhorse functions

daily.immigration <- function(crude.immigration.rate.f, data){
	daily.immigration.rate.per.1000 <- (crude.immigration.rate.f/1000)/365
	pop.size <- length(active.rows)
	mean.daily.n.immigrants <- daily.immigration.rate.per.1000*pop.size 
	todays.n.immigrants <- rpois(1, mean.daily.n.immigrants)		
	unoccupied.rows <- sum(is.na(data[,1]))
	return(todays.n.immigrants)
}

daily.emigration <- function(crude.emigration.rate.f, data){
	daily.emigration.rate.per.1000 <- (crude.emigration.rate.f/1000)/365
	pop.size <- length(active.rows)
	mean.daily.n.emigrants <- daily.emigration.rate.per.1000*pop.size
	todays.n.emigrants <- rpois(1, mean.daily.n.emigrants)
	return(todays.n.emigrants)		
}

emigrant.picker <- function(todays.n.emigrants.f, data){
	baseline.prob.emig <- 1/length(active.rows)
	t1s <- pop.reg[active.rows, "snoob"]
	alpha <- log(baseline.prob.emig/(1- baseline.prob.emig))
	emig.prob.vec <- exp(alpha + emigration.bias.snoob*t1s)/(1+exp(alpha + emigration.bias.snoob*t1s))
	emig.prob.vec <- emig.prob.vec/sum(emig.prob.vec)
	emigrants <- sample(active.rows, todays.n.emigrants.f, prob=emig.prob.vec)
	return(emigrants)
}

conception <- function(asf, data){
	
	age.cats <- data[fecund.women.rows, "age.cat"]
	snoob <- data[fecund.women.rows,"snoob"]
	
	baseline <- (asf[as.character(age.cats)]/1000)/365
	alpha = log(baseline/(1- baseline))
	
	daily.pr.conception <- (1-(1/(1+exp(alpha + fertility.bias.snoob*snoob))))
	got.pregnant <- rbinom(length(fecund.women.rows), 1, daily.pr.conception)
	return(got.pregnant)

}

grim.reaper <- function(asm, data){

	male <- data[active.rows,"male"]
	snoob <- data[active.rows,"snoob"]
	age.cats <- data[active.rows, "age.cat"]
	
	baseline <- (asm[as.character(age.cats)]/1000)/365
	alpha = log(baseline/(1- baseline))

	daily.pr.death <- (1-(1/(1+exp(alpha + mortality.bias.male*male + mortality.bias.snoob*snoob))))
	died <- rbinom(length(active.rows), 1, daily.pr.death)
	return(died)
	
}

age.binner <- function(age.vec){

	zeros <- which(age.vec < 365)
	hundreds <- which(age.vec >= 36500)

	age.years <- floor(age.vec/365)
	age.cats <- age.years - age.years %% 5
	age.cats[age.cats==0] <- 1
	age.cats[zeros] <- 0
	age.cats[hundreds] <- 100

	return(age.cats)

}



immigrant.maker <- function(data){ 
	state <- 1
	mom <- 0
	dad <- 0
	age <- sample(1:(365*50),1)
	age.cat <- age.binner(age)
	last.im <- age
	last.em <- NA
	died <- NA
	male <- rbinom(1,1,prob=0.5)
	mate <- NA
	counter <- NA
	
	h.gene <- sample( 80:120 , 1)
	height <- height.fun(age=age, h.gene=h.gene) 
	
	baseline.immigration <- mean(data[active.rows, "snoob"])
	alpha <- log(baseline.immigration/(1-baseline.immigration))
	
	snoob <- rbinom(1, 1,prob=(exp(alpha + immigration.bias.snoob)/(1+exp(alpha + immigration.bias.snoob))))
	
	language <- rbinom(1, 1, mean(data[active.rows, "language"]))
	
	aiy <- age/365
	a <- 0.18
	if(snoob==1) a <- a*snoob.less.risky.factor
	b <- 0.07
	risk <- a*aiy/exp(b*aiy)
	
	normal.traits <- c(snoob, language, risk)
	
	trait.vec <- c(h.gene, height, normal.traits)
		
	immigrant <- c(state, mom, dad, age, age.cat, last.im, last.em, died, male, mate, counter, trait.vec)
	immigrant
}


risk.setter <- function(data){
	aiy <- pop.reg[active.rows,"age"]/365
	snoob <- pop.reg[active.rows,"snoob"]
	a <- 0.18
	b <- 0.07
	a.vec <- rep(a, length(active.rows))
	a.vec[snoob==1] <- snoob.less.risky.factor*a
	risk <- a.vec*aiy/exp(b*aiy)
	return(risk)
}

# I should set the risk scoring to be (a) stochastic and (b) peaking at the same time as RS, which is around 25 yrs old.  

# a <- 0.18
# b <- 0.07
# curve(100*a*x/exp(b*x), from=0, to=90, add=T)
# abline(v=1/b, lty=2, col="gray")

baby.maker <- function(mom.id, dad.id, data){
	state <- 1
	age <- 0
	age.cat <- 0
	last.im <- NA
	last.em <- NA
	died <- NA
	male <- rbinom(1,1,prob=0.5)
	mate <- NA
	counter <- NA
	
	# snoob inheritance
	mom.snoob <- data[mom.id, "snoob"]
	dad.snoob <- data[dad.id, "snoob"]
	midparent.snoob <- mean(c(mom.snoob, dad.snoob))  # GOD DAMN that was a hard bug to solve...
	if(midparent.snoob == 1) midparent.snoob <- 0.9999
	if(midparent.snoob == 0) midparent.snoob <- 0.0001
	baseline.trans.bias <- midparent.snoob
	alpha <- log(baseline.trans.bias/(1-baseline.trans.bias))
	prob.kid.is.snoob <- exp(alpha + transmission.bias.snoob)/(1+ exp(alpha + transmission.bias.snoob))
	snoob <- rbinom(1,1,prob=prob.kid.is.snoob)
	
	# english inheritance
	mom.english <- data[mom.id, "language"]
	dad.english <- data[dad.id, "language"]
	parent.frq.english <- mean(c(mom.english, dad.english))
	parent.alpha <- 0.99
	pop.frq.english <- mean(pop.reg[active.rows,"language"])
	pr.acquire.english <- parent.alpha*parent.frq.english + (1-parent.alpha)*pop.frq.english
	language <- rbinom(1, 1, pr.acquire.english)
	
	risk <- 0
	
	normal.traits <- c(snoob, language, risk)
	
	h.gene <- sample( c(data[mom.id, "h.gene"], data[mom.id, "h.gene"]) , 1)
	height <- height.fun(age=age, h.gene=h.gene) 
	
	trait.vec <- c(h.gene, height, normal.traits)
	
	
	newborn <- c(state, mom.id, dad.id, age, age.cat, last.im, last.em, died, male, mate, counter, trait.vec)
	newborn
}


mate.finder <- function(mom.id, data){  # this seems the most inefficient of the simulator's code...it has to run mom by mom, which is quite a lot of them...
	
	current.mate <- NA
	
	moms.snoob <- data[mom.id, "snoob"]
	
	last.mate <- data[mom.id, "mate"]
	if(!is.na(last.mate)){
		if(data[last.mate, "state"] == 1 & data[last.mate, "age"] > 15*365 & data[last.mate, "age"] < 60*365){
			current.mate <- last.mate
		}
	}
	# if the last mate exists, is alive and within the age range, he's your man
	
	if(is.na(current.mate)){	
		# "available" means within age range, alive, male
		# male.mating.bias <- 1
		
		available.men <- active.rows[which( data[active.rows, "male"]==1 & data[active.rows, "age"] >= reproductive.min.age*365 & data[active.rows, "age"] < 60*365 )]
					
		if(length(available.men)>0){
			current.mate <- available.men[1]
			if(length(available.men) > 1){
				available.men.traits <- data[available.men,"snoob"]
				n.men <- length(available.men)
			
				prob.choose.snoob <- mate.similarity.bias.snoob*moms.snoob + (1-moms.snoob)*(1-mate.similarity.bias.snoob)
				prob.choose.nonsnoob <- 1 - prob.choose.snoob
	
				prob.is.chosen <- rep(0.5, n.men)
				prob.is.chosen[available.men.traits==1] <- prob.choose.snoob
				prob.is.chosen[available.men.traits==0] <- prob.choose.nonsnoob
								
				current.mate <- sample(available.men, 1, prob=prob.is.chosen)
			}
		} else {
			current.mate <- NA
		}
	}
	current.mate	

}



# height function: 

height.fun <- function(age, h.gene){
	if(any(age>100)) age <- age/365
	heights <- h.gene*exp(-3 + .2*age)/(1 + exp(-3+.2*age))
	heights
}



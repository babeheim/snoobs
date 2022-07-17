

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

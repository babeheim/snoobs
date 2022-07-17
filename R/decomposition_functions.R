
decompose_prepped_census <- function(data, check=F, drop.empty=T){

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

prep_census_data <- function(data, phenotype="snoob", grouping=NA){
      
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

#      if(j %% 10 == 0) print(j)
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
        
#    if(j %% 10 == 0) print(j)
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
    
#    if(j %% 10 == 0) print(j)
    
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
    
#    if(j %% 10 == 0) print(j)
    
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
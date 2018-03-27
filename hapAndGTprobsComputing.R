#' Returns a matrix of haplotype (rows) pobabilities in different single cells (columns).
#' 
#' @param hapStatus A \code{vector} of decoded halotype status (strings).
#' @param counts A one row \code{data.frame} containing the read counts of the segment.
#' @param chrCellTypes Cell types in the chromosome where the segment lies.
#' @inheritParams getPossibleCNs 
#' @param haplotypeMode TODO ...
#' @author Maryam Ghareghani
#' @export
#' 

newgetCellStatProbabilities = function(hapStatus, counts, chrCellTypes, p, chrCellsDispPars, binLength, alpha, haplotypeMode = FALSE)
{
  numCells = length(chrCellTypes)
  pr = matrix(, nrow = length(hapStatus), ncol = numCells)
  rownames(pr) = hapStatus
  colnames(pr) = 1:numCells
  segLen = as.integer(counts[,3]) - as.integer(counts[,2])
  
  for (i in 1:length(hapStatus))
  {
    for (j in 1:numCells)
    {
      segType = getSegType(chrCellTypes[j], hapStatus[i])
      disp = dispersionPar(segType, chrCellsDispPars[j], segLen, binLength, alpha)
      pr[i,j] = dnbinom(as.integer(counts[,2*j+2]), size = disp[1], prob = p)*dnbinom(as.integer(counts[,2*j+3]), size = disp[2], prob = p)
    }
  }
  
  if (! haplotypeMode)
  {
    for (i in 1:length(hapStatus))
    {
      sisterHapSt = sisterHaplotype(hapStatus[i])
      i2 = match(sisterHapSt, hapStatus)
      if (i < i2)
      {
        pr[i,] = (pr[i,]+pr[i2,])/2
        pr[i2,] = pr[i,]
      }
    }
  }
  
  # for (j in 1:numCells)
  # {
  #   if (sum(pr[,j]) > 0)
  #     pr[,j] = pr[,j]/sum(pr[,j])
  # }
  
  pr
}


#TODO remove or at least rename it later
#' Compute a noramlInvCN probability table given the haplotype probability table.
#' 
#' @param haplotypeProbTable A haplotype probability table with the set of decoded (non-binary) haplotypes for rows and cells for columns.
#' @author Maryam Ghareghani
#' @export
#' 

getGenotypeProbTable = function(haplotypeProbTable)
{
  haplotypes = rownames(haplotypeProbTable)
  genotypes = NULL
  hapIdx = list()
  
  for (i in 1:length(haplotypes))
  {
    #genotypeStatus = getGenotypeStatus(decodeStatus(haplotypes[i]))
    genotypeStatus = getGenotypeStatus(haplotypes[i])
    idx = which(genotypes == genotypeStatus)
    if (length(idx) == 0)
    {
      genotypes = c(genotypes, genotypeStatus)
      hapIdx[[length(hapIdx)+1]] = i
    }
    else
    { # length(idx) should be 1
      hapIdx[[idx]] = c(hapIdx[[idx]], i)
    }
  }
  
  genotypeProbTable = matrix(, nrow = length(genotypes), ncol = ncol(haplotypeProbTable))
  rownames(genotypeProbTable) = genotypes
  
  for (i in 1:length(genotypes))
  {
    if (length(hapIdx[[i]]) == 1)
    {
      genotypeProbTable[i,] = haplotypeProbTable[hapIdx[[i]],]
    } else {
      genotypeProbTable[i,] = colSums(haplotypeProbTable[hapIdx[[i]],])#/length(hapIdx[[i]]) # I think we shouldn't devide it by the number of haplotypes corresponding to the genotype
    }
  }
  
  genotypeProbTable
}


#' Compute a genotype probability table given the haplotype probability table.
#' 
#' @param haplotypeProbTable A haplotype probability table with the set of decoded (non-binary) haplotypes for rows and cells for columns.
#' @author Maryam Ghareghani
#' @export
#' 

newgetGenotypeProbTable = function(haplotypeProbTable)
{
  haplotypes = rownames(haplotypeProbTable)
  genotypes = NULL
  hapIdx = list()
  
  #for (i in 1:length(haplotypes))
  #{
  #  haplotypes[i] = substr(haplotypes[i], 2, nchar(haplotypes[i]))
  #}
  
  for (i in 1:length(haplotypes))
  {
    sisterHap = sisterHaplotype(haplotypes[i])
    sisPos = match(sisterHap, haplotypes)
    if (sisPos > i-1)
    {
      hapIdx[[length(hapIdx)+1]] = unique(c(sisPos,i))
      genotypes = c(genotypes, haplotypes[i])
    }
  }
  
  genotypeProbTable = matrix(, nrow = length(genotypes), ncol = ncol(haplotypeProbTable))
  rownames(genotypeProbTable) = genotypes
  
  for (i in 1:length(genotypes))
  {
    if (length(hapIdx[[i]]) == 1)
    {
      genotypeProbTable[i,] = haplotypeProbTable[hapIdx[[i]],]
    } else {
      genotypeProbTable[i,] = colSums(haplotypeProbTable[hapIdx[[i]],])
    }
  }
  
  genotypeProbTable
}


#' Takes as input a haplotype probTable and outout the genotype probTable.
#' 
#' @param hapProbDF A \code{data.frame} containing haplotype probability table.
#' @author Maryam Ghareghani
#' @export

getGenotypeProbDataFrame = function(hapProbDF)
{
  GTs = list()
  # blacklist the names of non probability columns
  haplotypeNames = colnames(hapProbDF[8:ncol(hapProbDF)])
  maximumCN = length(which(grepl("CN*", haplotypeNames)))-1
  # blacklist copy number names
  haplotypeNames = haplotypeNames[(maximumCN+2):length(haplotypeNames)]
  GTnames = NULL
  
  if (startsWith(haplotypeNames[1],"X"))
  {
    haplotypeNames = sapply(haplotypeNames, function(x) substr(x, 2, nchar(x)))
    names(haplotypeNames) = NULL
  }
  
  for (i in 1:length(haplotypeNames))
  {
    sisterHap = sisterHaplotype(haplotypeNames[i])
    sisPos = match(sisterHap, haplotypeNames)
    if (sisPos > i-1)
    {
      GTs[[length(GTs)+1]] = unique(c(sisPos,i))
      GTnames = c(GTnames, haplotypeNames[i])
    }
  }
  
  GTprobDF = hapProbDF[,1:(maximumCN+8)]
  for (i in 1:length(GTs))
  {
    if (length(GTs[[i]]) == 1)
    {
      GTprobDF = cbind(GTprobDF, hapProbDF[,maximumCN+8+GTs[[i]]])
    } else
    {
      GTprobDF = cbind(GTprobDF, hapProbDF[,maximumCN+8+GTs[[i]][1]]+hapProbDF[,maximumCN+8+GTs[[i]][2]])
    }
  }
  
  colnames(GTprobDF)[(maximumCN+9):ncol(GTprobDF)] = GTnames
  
  GTprobDF
}


#' Normalize a probability table to a table in which each column sums up to 1.
#' 
#' @param probTable A probability table.
#' @author Maryam Ghareghani
#' @export
#' 

normalizeProbTable = function(probTable)
{
  colsum = colSums(probTable)
  newProbTable = probTable
  for (j in 1:ncol(probTable))
  {
    if (colsum[j] != 0)
    {
      newProbTable[,j] = newProbTable[,j]/colsum[j]
    }
  }
  
  newProbTable
}


#' We assume that the posteriori probability for different status is uniform with probability regFactor and a probability computed based on the NB model with probability (1-regFactor).
#' 
#' Regularize the probability table according to the mentioned assumption.
#'
#' @param probTable A non-regularized probability table in which every row corresponds to a status and each column corresponds to a single cell.
#' @param regFactor The regularization factor.
#' @author Maryam Ghareghani
#' @export
#' 

regularizeProbTable = function(probTable, regFactor = 1e-10)
{
  newProbTable = probTable
  for (j in 1:ncol(probTable))
  {
    newProbTable[,j] = regFactor/nrow(probTable) + (1-regFactor)*newProbTable[,j]
    newProbTable[,j] = newProbTable[,j]/sum(newProbTable[,j])
  }
  newProbTable
}


#' For each segment, computes the probabilities of having the same haplotype for all single cells
#' 
#' @param probTable.l The list of segments probability table
#' @author Maryam Ghareghani
#' @export
#' 



jumpProbs <- function(probTable.l)
{
  chroms <- sapply(probTable.l, function(x) x$chr[1])
  probTable.l.chroms <- split(probTable.l, factor(chroms, levels = unique(chroms)))
  jump.probs <- list()
  for (k in 1:length(probTable.l.chroms)) {
    jump.probs[[k]] <- list()
    for (i in 1:(length(probTable.l.chroms[[k]])-1)) {
      prod.probs <- probTable.l.chroms[[k]][[i]][,(maximumCN+9):n] * probTable.l.chroms[[k]][[i+1]][,(maximumCN+9):n]
      jump.probs[[k]][[i]] <- rowSums(prod.probs)
    }
  }
  
  names(jump.probs) <- names(probTable.l.chroms)
  return(jump.probs)
}

#' returns the genotyoe class (CN_loss, CN_gain, normal_CN2, or inv_CN2)
#' 
#' @param genotype.name The genotype or haplotype coding
#' @author Maryam Ghareghani
#' @export
#' 

genotype_class <- function(genotype.name)
{
  geno.class <- ""
  if (startsWith(genotype.name, "X"))
  {
    genotype.name <- substr(genotype.name, 2, nchar(genotype.name))
  }
  genotype.int.vec <- sapply(1:nchar(genotype.name), function(x) as.numeric(substr(genotype.name,x,x)))
  CN <- sum(genotype.int.vec)
  
  if (CN < 2)
  {
    geno.class <- "CN_loss"
  } else if (CN > 2 )
  {
    geno.class <- "CN_ gain"
  } else # CN = 2
  {
    if (genotype.int.vec[2] + genotype.int.vec[4] > 0)
    {
      geno.class <- "inv_CN2"
    } else
    {
      geno.class <- "normal_CN2"
    }
  }
  
  return(geno.class)
}

#' returns the genotyoe class (CN_loss, CN_gain, normal_CN2, or inv_CN2) probability dataframe
#' 
#' @param prob.tab The genotype or haplotype probability dataframe
#' @author Maryam Ghareghani
#' @export
#' 

get_GT_class_prob_table <- function(prob.tab)
{
  # getting the index of the last column before CN probs
  last.CN.col.idx <- max(which(startsWith(names(prob.tab), "CN")))
  
  # getting genotype names
  gt.names <- names(prob.tab)[(last.CN.col.idx+1):ncol(prob.tab)]
  
  # getting genotype classes
  gt.classes <- sapply(gt.names, genotype_class)
  
  prob.gt.class <- prob.tab[,1:7]
  for (g in unique(gt.classes))
  {
    g.idx <- which(gt.classes == g)+last.CN.col.idx
    prob.gt.class <- cbind(prob.gt.class, data.frame(rowSums(prob.tab[g.idx])))
  }
  names(prob.gt.class)[8:ncol(prob.gt.class)] <-unique(gt.classes)
  
  return(prob.gt.class)
}

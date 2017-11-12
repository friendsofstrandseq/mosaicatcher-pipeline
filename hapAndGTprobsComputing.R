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

newgetCellStatProbabilities = function(hapStatus, counts, chrCellTypes, p, chrCellsDispPars, binLen, alpha, haplotypeMode = FALSE)
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
      disp = dispersionPar(segType, chrCellsDispPars[j], segLen, binLen, alpha)
      pr[i,j] = dnbinom(as.integer(counts[,2*j+2]), size = disp[1], prob = p)*dnbinom(as.integer(counts[,2*j+3]), size = disp[2], prob = p)
    }
  }
  
  if (! haplotypeMode)
  {
    for (i in 1:length(hapStates))
    {
      sisterHapSt = sisterHaplotype(hapStates[i])
      i2 = match(sisterHapSt, hapStates)
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
  haplotypeNames = colnames(hapProbDF[14:ncol(hapProbDF)])
  GTnames = NULL
  
  for (i in 1:length(haplotypeNames))
  {
    haplotypeNames[i] = substr(haplotypeNames[i], 2, nchar(haplotypeNames[i]))
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
  
  GTprobDF = hapProbDF[,1:13]
  for (i in 1:length(GTs))
  {
    if (length(GTs[[i]]) == 1)
    {
      GTprobDF = cbind(GTprobDF, hapProbDF[,13+GTs[[i]]])
    } else
    {
      GTprobDF = cbind(GTprobDF, hapProbDF[,13+GTs[[i]][1]]+hapProbDF[,13+GTs[[i]][2]])
    }
  }
  
  colnames(GTprobDF)[14:ncol(GTprobDF)] = GTnames
  
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

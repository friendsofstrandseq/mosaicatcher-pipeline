#' returns the hapotyoe name
#' 
#' @param hap.code The haplotype coding
#' @author Maryam Ghareghani
#' @export
#' 

get_hap_name <- function(hap.code)
{
  hap.codes <- c("1010", "0010", "1000", "0000", "0110", "1001", "0101", "2010", "1020", "2020")
  hap.names <- c("ref_hom", "del_h1", "del_h2", "del_hom", # ref and del
                 "inv_h1", "inv_h2", "inv_hom", # inv
                 "dup_h1", "dup_h2", "dup_hom") # dup
  
  hap.idx <- match(hap.code, hap.codes)
  
  if (!is.na(hap.idx))
  {
    return(hap.names[hap.idx])
  }
  
  return("complex")
}

#' converts the haplotype to the genotype name
#' 
#' @param hap.name The haplotype name
#' @author Maryam Ghareghani
#' @export
#' 

haplo_to_geno_name <- function(hap.name)
{
  geno.name <- gsub("h1","het",hap.name)
  geno.name <- gsub("h2","het",geno.name)
  
  return(geno.name)
}
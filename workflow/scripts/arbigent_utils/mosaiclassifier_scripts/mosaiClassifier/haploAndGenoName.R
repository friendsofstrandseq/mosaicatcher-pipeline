#' returns the hapotyoe name
#' 
#' @param hap.code The haplotype coding
#' @author Maryam Ghareghani
#' @export
#' 

get_hap_name <- function(hap.code)
{
  hap.codes <- c("1010", "0010", "1000", "0000", "0110", "1001", "0101", "2010", "1020", "2020", "1110", "1011")
  hap.names <- c("ref_hom", "del_h1", "del_h2", "del_hom", # ref and del
                 "inv_h1", "inv_h2", "inv_hom", # inv
                 "dup_h1", "dup_h2", "dup_hom", # dup
                 "idup_h1", "idup_h2") # inv-dup
  
  hap.idx <- match(hap.code, hap.codes)
  
  if (!is.na(hap.idx))
  {
    return(hap.names[hap.idx])
  }
  
  return("complex")
}


#' returns the genotype code given the genotyoe name
#' 
#' @param geno.name genotype name
#' @author Maryam Ghareghani
#' @export
#' 

geno_name_to_code <- function(geno.name)
{
  geno.codes <- c("1010", "1010", 
                 "0010", "0000", # del
                 "0110", "0101", # inv
                 "2010", "2020", # dup
                 "1110") # inv-dup
  
  geno.names <- c("hom_ref", "false_del",
                 "het_del", "hom_del", # del
                 "het_inv", "hom_inv", # inv
                 "het_dup", "hom_dup", # dup
                 "inv_dup") # inv-dup
  
  geno.idx <- match(geno.name, geno.names)
  
  if (!is.na(geno.idx))
  {
    return(geno.codes[geno.idx])
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

#' translates the haplotype code to the corresponding classe of genotypes (normal, inv, CN loss, CN gain)
#' 
#' @param hap.code The haplotype coding
#' @author Maryam Ghareghani
#' @export
#' 

haplo_code_to_geno_class <- function(hap.code)
{
  if (length(hap.code) < 1) return (character())
  
  dd = as.data.table(str_split(hap.code,"", simplify = T, n = 4))
  dd = dd[, lapply(.SD, as.integer)]
  dd[, state := ifelse(V1+V2+V3+V4 != 2, 
                       ifelse(V1+V2+V3+V4<2, "loss", "gain"), 
                       ifelse(V2+V4>0, "inv", "ref") )]
  return(dd$state)
}

#' get CN from haplotype
#' 
#' @param hap.code The haplotype coding
#' @author Maryam Ghareghani
#' @export
#' 

haplo_code_to_geno_class <- function(hap.code)
{
  if (length(hap.code) < 1) return (character())
  
  dd = as.data.table(str_split(hap.code,"", simplify = T, n = 4))
  dd = dd[, lapply(.SD, as.integer)]
  dd[, CN:=V1+V2+V3+V4, by=1:nrow(dd)]
  return(dd$CN)
}

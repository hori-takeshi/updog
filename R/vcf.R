## Wrapper functions for loading in vcf file and fitting mupdog.

#' Wrapper for \code{\link[VariantAnnotation]{readVcf}},
#' \code{\link[VariantAnnotation]{writeVcf}}, and
#' \code{\link{flexdog}}.
#'
#' Reads in a VCF file with \code{\link[VariantAnnotation]{readVcf}},
#' fits \code{\link{flexdog}} to each SNP in the VCF file, then
#' writes the output as a VCF file using
#' \code{\link[VariantAnnotation]{writeVcf}}.
#'
#' @param input The input file. This file should be in VCF format.
#' @param output The ouptut file. This file will be in VCF format.
#' @param ploidy The ploidy of the species.
#' @param nc The number of cores to be use if parallelization
#'     is desired. We use
#'     the \code{\link[parallel]{parallel}} package
#'     for parallel computing.
#' @param alt_field The name of the genotype field in \code{input}
#'     that contains the alternative counts.
#' @param size_field The name of the genotype field in \code{input}
#'     that contains the read-depth.
#' @param ... Further arguements to pass to \code{\link{flexdog}}.
#'
#' @return The file \code{output} is created that
#'     contains the new INFO fields:
#' \describe{
#' \item{BIAS}{The allele (ascertainment) bias. The probability of ascertaining
#'     a reference read divided by the probability of ascertaining an alternative
#'     read.}
#' \item{SEQ}{The sequencing error rate.}
#' \item{OD}{The overdispersion parameter.}
#' \item{PM}{The posterior proportion of individuals misgenotyped.}
#' \item{PRIOR}{The estimated prior distribution. Comma separated.}
#' }
#' The \code{output} file also contains the new GENO fields:
#' \describe{
#' \item{GT}{The (unphased) gentotype. In the format 1/1/0/0 where 1 represents the
#'     alternative allele and 0 represents the reference allele.}
#' \item{DS}{The posterior mean dosage. Defined between 0 and \code{ploidy}.}
#' \item{GQ}{The posterior probability of misgenotyping.}
#' \item{GP}{The posterior probability of each genotype.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
#' @import VariantAnnotation
#' @importFrom S4Vectors DataFrame
#'
#' @references Gerard, David, Luis Felipe Ventorim Ferrao,
#'     Antonio Augusto Franco Garcia, and Matthew Stephens. 2018.
#'     "Harnessing Empirical Bayes and Mendelian Segregation
#'     for Genotyping Autopolyploids from Messy Sequencing Data."
#'     \emph{bioRxiv}. Cold Spring Harbor Laboratory.
#'     doi:10.1101/281550.
#'
#'     Valerie Obenchain, Michael Lawrence,
#'     Vincent Carey, Stephanie Gogarten, Paul Shannon,
#'     Martin Morgan; VariantAnnotation : a Bioconductor
#'     package for exploration and annotation of
#'     genetic variants , \emph{Bioinformatics}, Volume 30,
#'     Issue 14, 15 July 2014, Pages 2076-2078,
#'     \url{https://doi.org/10.1093/bioinformatics/btu168}.
#'
#' @seealso \code{\link{flexdog}} for the underlying fitting function.
#'
#'
#'
vcfdog <- function(input,
                   output,
                   ploidy,
                   nc         = 1,
                   alt_field  = "AA",
                   size_field = "DP",
                   ...) {

  ## Check input ---------------------------
  assertthat::are_equal(length(input), length(output),
                        length(nc), length(alt_field),
                        length(size_field), length(ploidy),
                        1)
  assertthat::assert_that(nc >= 1)
  assertthat::assert_that(ploidy >= 1)
  assertthat::are_equal(ploidy %% 1, 0)
  assertthat::are_equal(nc %% 1, 0)
  assertthat::assert_that(is.character(input))
  assertthat::assert_that(is.character(output))
  assertthat::assert_that(is.character(alt_field))
  assertthat::assert_that(is.character(size_field))

  ## Read in vcf file ----------------------
  vout <- readVcf(file = input, genome = "updog")

  altmat  <- geno(vout)[[alt_field]]
  sizemat <- geno(vout)[[size_field]]

  nsnps <- nrow(sizemat)
  nind  <- ncol(sizemat)

  ## Fit flexdog ----------------------------
  if (nc == 1) {
    foreach::registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(nc)
    doParallel::registerDoParallel(cl = cl)
    stopifnot(foreach::getDoParWorkers() == nc)
  }
  index <- 1 ## stupid hack to get around CRAN note.
  retlist <- foreach::foreach(index = seq_len(nsnps),
                   .export = c("flexdog",
                               "dosage_to_gt",
                               "prob_to_phred")) %dopar% {
    uout <- flexdog(refvec  = altmat[index, ],
                    sizevec = sizemat[index, ],
                    ploidy  = ploidy,
                    verbose = FALSE,
                    ...)
    retlist <- list()
    retlist$GT <- vapply(X         = uout$geno,
                         FUN       = dosage_to_gt,
                         FUN.VALUE = "character",
                         ploidy    = ploidy)
    retlist$DS    <- uout$postmean
    retlist$GQ    <- prob_to_phred(p = 1 - uout$maxpostprob)
    retlist$GP    <- prob_to_phred(uout$postmat)
    retlist$BIAS  <- uout$bias
    retlist$SEQ   <- uout$seq
    retlist$OD    <- uout$od
    retlist$PM    <- uout$prop_mis
    retlist$PRIOR <- paste(format(prob_to_phred(uout$gene_dist), digits = 2, trim = TRUE, scientific = FALSE), collapse = ",")
    retlist
  }
  if (nc > 1) {
    parallel::stopCluster(cl)
  }

  ## Collapse data ---------------------------------------
  GTmat    <- t(vapply(X = retlist, FUN = function(x){x$GT}, FUN.VALUE = character(nind)))
  DSmat    <- t(vapply(X = retlist, FUN = function(x){x$DS}, FUN.VALUE = numeric(nind)))
  GParray  <- aperm(vapply(X = retlist, FUN = function(x){x$GP}, FUN.VALUE = matrix(numeric(1), nrow = nind, ncol = ploidy + 1)), c(3, 1, 2))
  GQmat    <- t(vapply(X = retlist, FUN = function(x){x$GQ}, FUN.VALUE = numeric(nind)))
  BIASvec  <- vapply(X = retlist, FUN = function(x){x$BIAS}, FUN.VALUE = numeric(1))
  SEQvec   <- vapply(X = retlist, FUN = function(x){x$SEQ}, FUN.VALUE = numeric(1))
  ODvec    <- vapply(X = retlist, FUN = function(x){x$OD}, FUN.VALUE = numeric(1))
  PMvec    <- vapply(X = retlist, FUN = function(x){x$PM}, FUN.VALUE = numeric(1))
  PRIORvec <- vapply(X = retlist, FUN = function(x){x$PRIOR}, FUN.VALUE = character(1))

  ## Write and save a vcf object -------------------------
  newdf <- DataFrame(Number = c(1, 1, 1, ploidy + 1),
                      Type   = c("String", "Float", "Float", "Float"),
                      Description = c("Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype",
                                      "Posterior mean genotype",
                                      "Phred-scaled posterior probability of misgenotyping",
                                      "Phred-scaled posterior probability of each dosage"),
                     row.names = c("GT", "DS", "GQ", "GP"))

  suppressWarnings({
    geno(header(vout)) <- newdf
    for (current_name in names(geno(vout))) {
      geno(vout)[[current_name]] <- NULL
    }
  })

  VariantAnnotation::geno(vout)$GT <- GTmat
  VariantAnnotation::geno(vout)$DS <- DSmat
  VariantAnnotation::geno(vout)$GQ <- GQmat
  VariantAnnotation::geno(vout)$GP <- GParray

  newdf <- DataFrame(Number = c(1, 1, 1, 1, ploidy + 1),
                      Type   = c("Float", "Float", "Float", "Float", "Float"),
                      Description = c("Allele/ascertainment bias. Probability of ascertaining reference divided by probability of ascertaining alternative.",
                                      "Sequencing error rate.",
                                      "Overdispersion parameter.",
                                      "Posterior proportion of individuals mis-genotyped",
                                      "Phred-scaled genotype distribution."),
                     row.names = c("BIAS", "SEQ", "OD", "PM", "PRIOR"))

  suppressWarnings({
    info(header(vout)) <- newdf
    for (current_name in names(info(vout))) {
      info(vout)[[current_name]] <- NULL
    }
  })

  info(vout)$BIAS  <- BIASvec
  info(vout)$SEQ   <- SEQvec
  info(vout)$OD    <- ODvec
  info(vout)$PM    <- PMvec
  info(vout)$PRIOR <- PRIORvec

  function_call <- utils::capture.output(print(match.call()))
  newdf <- DataFrame(Value = c(as.character(Sys.time()),
                               as.character(utils::packageVersion("updog")),
                               function_call),
                     row.names = c("fileDate",
                                   "source",
                                   "commandline"))

  meta(header(vout))$META <- newdf

  writeVcf(obj = vout, filename = output)
}

#' Convert the dosage to the common (unphased) GT format in VCF files.
#'
#' @param x The allele dosage (proportion of alternative alleles).
#' @param ploidy The ploidy of the species.
#'
#' @return A character in the form 1/1/0/0 where 1 represents the alternative
#'     allele and 0 represents the reference allele.
#'
#' @author David Gerard
dosage_to_gt <- function(x, ploidy) {
  if (is.na(x)) {
    retval <- paste(rep(x = ".", times = ploidy),
                    collapse = "/")
  } else {
    retval <- paste(c(rep(x = "1", times = x),
                      rep(x = "0", times = ploidy - x)),
                    collapse = "/")
  }
  return(retval)
}

#' Converts a proportion to a phred-scaled probability.
#'
#' @param p A proportion. Between 0 and 1.
#'
#' @return \code{-10 * log(p, base = 10)}.
#'
#' @author David Gerard
#'
prob_to_phred <- function(p) {
  -10 * log(x = p, base = 10)
}

#' A function to turn non-corrected into corrected labelled peptide/protein fractions (Corr LPFs)
#'
#' This function allows users to input their non-corrected LPF matrix from "AnnotateProteins.R" and correct the fractional enrichment in individual peptides using the enrichment percentages in amino acid soluble pools from paired treatments. The function returns the corrected matrix, which can be directly used to calculate fractional synthesis rates.
#' @param InputNonCorrMat
#' @param files2correct Defaults to NULL. List of file directories, where each file contains the enrichments in soluble amino acids per treatment evaluated.
#' @param AA4correction Defaults to NULL. Character vector with the names of the amino acids to be used for the calculations. The amino acids must be present in the files from the previous parameter.
#' @param AAinterprtFileDir File directory to the interpretation file that contains in one column amino acid names and in the second colum the single letter code for those amino acids.
#' @param ProtPTMs Character vector with the PTM codes that are part of peptides in the provided dataset.
#' @param CorrectMeans Defaults to FALSE. Allows users to return corrected treatment means instead of individual replicates.
#' @param EnrBoundary Defaults to one. Allows to define the boundarie of fractional enrichment accepted to be returned. 1 equals to 100%.
#' @param GroupPeptides Defaults to FALSE. Allows to group all peptides from a single protein into mean protein enrichment.
#' @keywords LCMS kProteomics Tracer Dynamics
#' @export
#' @examples
#'
#' ...

LPFcorrection <- function(InputNonCorrMat,
                          files2correct = NULL,
                          AA4correction = NULL,
                          AAinterprtFileDir,
                          ProtPTMs,
                          CorrectMeans = FALSE,
                          EnrBoundary = 1,
                          GroupPeptides = T) {

  ### needed functions

  LPFcorr <- function(Col2correct,
                      CorrectionMat,
                      PTMs,
                      CorrectMeans = CorrectMeans){

    ##### write a small function to remove PTMs from peptide sequences before counting AAs (use what was already written in one of the other functions)

    if(CorrectMeans == T) {

      corrected_LPF <- c()

    } else {

      corrected_LPF <- matrix(NA, nrow = nrow(Col2correct), ncol = ncol(Col2correct))

    }

    #count = 0

    for (j in 1:length(rownames(Col2correct))) {

      #count = count +1

      #cat(count, " ")

      peptide2color <- rownames(Col2correct)[j]

      if (length(strsplit(peptide2color, "-")[[1]]) > 1){

        PTMs <- PTMs

        split_PTM_peptide <- strsplit(peptide2color, "-")[[1]]

        out_vec <- c()

        for (i in 1:length(split_PTM_peptide)) {

          test <- split_PTM_peptide[i]

          test1 <- strsplit(split_PTM_peptide[i], "")[[1]]

          if(length(test1) > 0 & length(unique(test != PTMs)) == 1) {

            out_vec <- c(out_vec, split_PTM_peptide[i])

          }
        }

        peptide2color <- paste0(out_vec, collapse = "")

        pep_sequence <- strsplit(strsplit(peptide2color, "_")[[1]][2], "")[[1]]

      } else {

        pep_sequence <- strsplit(strsplit(peptide2color, "_")[[1]][2], "")[[1]]

      }

      ##### pep_sequence is already stripped off PTMs - ready to match enriched AAs

      AAs_MeanEnrichment <- c()

      for (i in 1:nrow(file_runner_corr_adjusted)) {

        AA_i <- rownames(file_runner_corr_adjusted)[i]

        test_AA_in_sequence <- grep(pattern = AA_i, x = pep_sequence)

        if(length(test_AA_in_sequence) > 0){

          AA_i_number <- length(grep(pattern = AA_i, x = pep_sequence))

          AA_i_MeanEnr <- mean(apply(file_runner_corr_adjusted[AA_i,], 2, as.numeric))

          AAs_MeanEnrichment <- c(AAs_MeanEnrichment, AA_i_MeanEnr)

        }
      }

      ##### correct means with means and not rep with reps

      if(CorrectMeans == T) {

        ColMeans2correct <- rowMeans(apply(Col2correct, MARGIN = 2, FUN = as.numeric))

        corrected_LPF[j] <- ColMeans2correct[j]/mean(AAs_MeanEnrichment)

      } else {

        corrected_LPF[j,] <- apply(Col2correct, MARGIN = 2, FUN = as.numeric)[j,]/mean(AAs_MeanEnrichment)

      }
    }

    if(CorrectMeans == T) {

      names(corrected_LPF) <- rownames(Col2correct)

    } else {

      colnames(corrected_LPF) <- colnames(Col2correct)
      rownames(corrected_LPF) <- rownames(Col2correct)

    }

    return(corrected_LPF)

  }


  ### main

  #### Import files

  if(!is.null(files2correct)) {

    test_annotateProt_NL <- InputNonCorrMat

    out_LPF_corrected_list <- list()

    out_LPF_corrected_names_list <- list()

    for (i in 1:length(files2correct)) {

      file_runner_corr <- read.csv(files2correct[i])

      if(!is.null(AA4correction)) {

        test_names_AA_peaks <- grep(pattern = paste0(AA4correction, collapse = "|"),
                                    file_runner_corr$X)

        if(length(test_names_AA_peaks)>0) {

          file_runner_corr <- file_runner_corr[test_names_AA_peaks,]

        } else {

          stop("ERROR: Amino acid peak names in the correction files do not match the names in AA4correction param")

        }
      }

      cat(paste0(paste(colnames(test_annotateProt_NL),
                       1:length(colnames(test_annotateProt_NL)),
                       sep = " "),
                 collapse = "\n"))

      columns2correct <- readline(paste0("which columns equal to these three samples (x): ",
                                         paste0(colnames(file_runner_corr), collapse = ",") , " ? (e.g.: 1,2,3)"))

      colIndexes2correct <- as.numeric(strsplit(columns2correct, ",")[[1]])

      subset_mat2correct <- test_annotateProt_NL[,colIndexes2correct]
      rownames(subset_mat2correct) <- test_annotateProt_NL[,"PeptideID"]

      #### Correct individual replicates from k-LC/MS with individual paired replicates from k-GC/MS
      #### When impossible dispose peptides as low quality data

      AAinterprtFile <- read.csv(AAinterprtFileDir)

      test_exist_AA <- grep(pattern = paste0(AA4correction, collapse = "|"), AAinterprtFile$Amino.acid)

      if (length(test_exist_AA) > 0) {

        ##### loop necessary to have the AA in the right order

        runner_name_single_letter <- c()

        for (j in 1:length(AA4correction)) {

          runner_name_single_letter[j] <- AAinterprtFile$Single.letter.abbreviation[grep(pattern = AA4correction[j],
                                                                                         AAinterprtFile$Amino.acid)]
        }

      } else {

        stop("ERROR: given amino acid names in AA4correction param do not coincide with any amino acid in AAinterprtFileDir")

      }

      file_runner_corr_adjusted <- file_runner_corr[,2:ncol(file_runner_corr)]
      rownames(file_runner_corr_adjusted) <- runner_name_single_letter

      out_LPF_corrected_mat <- LPFcorr(Col2correct = subset_mat2correct,
                                       CorrectionMat = file_runner_corr_adjusted,
                                       PTMs = ProtPTMs,
                                       CorrectMeans = CorrectMeans)

      #### deleting peptides with impossible enrichment percentages that are product of low quality data

      if(CorrectMeans == T) {

        test_zero_peps <- which(out_LPF_corrected_mat > EnrBoundary)

        if(length(test_zero_peps) == 0) {

          out_LPF_corrected_list[[i]] <- out_LPF_corrected_mat

          out_LPF_corrected_names_list[[i]] <- names(out_LPF_corrected_mat)

        } else {

          out_LPF_corrected_list[[i]] <- out_LPF_corrected_mat[-c(which(out_LPF_corrected_mat > EnrBoundary))]

          out_LPF_corrected_names_list[[i]] <- names(out_LPF_corrected_mat[-c(which(out_LPF_corrected_mat > EnrBoundary))])

        }

      } else {

        out_LPF_corrected_mat <- out_LPF_corrected_mat[which(rowMeans(out_LPF_corrected_mat) < EnrBoundary),]

        out_LPF_corrected_list[[i]] <- out_LPF_corrected_mat

        out_LPF_corrected_names_list[[i]] <- rownames(out_LPF_corrected_mat)

      }

    }

    good_quality_peptides <- Reduce(f = intersect, out_LPF_corrected_names_list)


    if(CorrectMeans == T) {

      out_mat <- matrix(NA, nrow = length(good_quality_peptides), ncol = length(files2correct))

      for (i in 1:length(out_LPF_corrected_list)) {

        out_mat[,i] <- out_LPF_corrected_list[[i]][good_quality_peptides]

      }

      colnames(out_mat) <- paste0("Corr_LPF", list2DF(strsplit(files2correct, "\\.csv")))
      rownames(out_mat) <- good_quality_peptides

    } else {

      count_ncol <- 0

      for (i in 1:length(out_LPF_corrected_list)) {

        count_ncol <- count_ncol + dim(out_LPF_corrected_list[[i]])[2]

      }

      out_mat <- suppressWarnings(matrix(NA, nrow = length(good_quality_peptides), ncol = 0))

      colnames_bound <- c()

      for (i in 1:length(out_LPF_corrected_list)) {

        out_mat <- cbind(out_mat, out_LPF_corrected_list[[i]][good_quality_peptides,])

        colnames_bound <- c(colnames_bound, rep(paste0("Corr_LPF", list2DF(strsplit(files2correct[i], "\\.csv"))), dim(out_LPF_corrected_list[[i]])[2]))

      }

      colnames(out_mat) <- colnames_bound

      rownames(out_mat) <- good_quality_peptides

    }
  }

  #### calling the fully annotated matrix to supplement it with the new data

  rownames(test_annotateProt_NL) <- test_annotateProt_NL[,"PeptideID"]

  LPF_plus_mat <- cbind(out_mat, test_annotateProt_NL[good_quality_peptides,])

  #### group proteins into a single entry to unify information, conserve multiple peptides in the appropriate slots

  if(GroupPeptides == T) {

    ProteinID_vector <- unique(LPF_plus_mat[,"HORVU"])

  } else {

    ProteinID_vector <- unique(LPF_plus_mat[,"PeptideID"])

  }

  consolidated_mat <- matrix(NA, nrow = length(ProteinID_vector), ncol = ncol(LPF_plus_mat))

  for (i in 1:length(ProteinID_vector)) {

    if(GroupPeptides == T) {

      runner_consolid_mat <- as.matrix(LPF_plus_mat[grep(pattern = ProteinID_vector[i], x = LPF_plus_mat[,"HORVU"]),])

    } else {

      runner_consolid_mat <- as.matrix(LPF_plus_mat[grep(pattern = ProteinID_vector[i], x = LPF_plus_mat[,"PeptideID"]),])

    }


    if(ncol(runner_consolid_mat) == 1){

      vector_consolid <- as.character(runner_consolid_mat)

      names(vector_consolid) <- colnames(LPF_plus_mat)

    } else {

      if(CorrectMeans == T) {

        vector_consolid <- c(colMeans(apply(runner_consolid_mat[,1:2], 2, as.numeric)),
                             unique(runner_consolid_mat[,3:4]),
                             paste0(runner_consolid_mat[,5], collapse = ";"),
                             unique(runner_consolid_mat[,6]),
                             paste0(runner_consolid_mat[,7], collapse = ";"),
                             colMeans(apply(runner_consolid_mat[,8:ncol(runner_consolid_mat)], 2, as.numeric)))

        names(vector_consolid) <- colnames(LPF_plus_mat)

      } else {

        vector_consolid <- c(colMeans(apply(runner_consolid_mat[,1:6], 2, as.numeric)),
                             unique(runner_consolid_mat[,7:8]),
                             paste0(runner_consolid_mat[,9], collapse = ";"),
                             unique(runner_consolid_mat[,10]),
                             paste0(runner_consolid_mat[,11], collapse = ";"),
                             colMeans(apply(runner_consolid_mat[,12:ncol(runner_consolid_mat)], 2, as.numeric)))

        names(vector_consolid) <- colnames(LPF_plus_mat)

      }

    }

    consolidated_mat[i,] <- vector_consolid

  }

  colnames(consolidated_mat) <- colnames(LPF_plus_mat)

  return(consolidated_mat)

}

#' A function to identify and compile protein enrichment information from a fed peptide subset
#'
#' This function allows users to identify the parent proteins from peptide subsets, highlight the coverage of the peptides in the protein sequence, retrieve and visualize the noise-corrected enrichments (non-corrected LPFs) and return an output table that is necessary for the next function in the workflow, which corrects LPFs.
#' @param PeptideVector Character vector with peptide identifiers as delivered by the EnrichmentSet.R function.
#' @param Treatment Factor vector with the containing sample information in the same order as samples are outlined in the enrichment file.
#' @param FileName Output file name. The output file is a PDF that contains boxplots featuring individual peptide non-corrected LPFs across samples and highlighted peptide sequences in the overall protein sequence provided in the FASTA file.
#' @param Path2FASTA FASTA file used for the analysis of the dataset in MaxQuant.
#' @param Path2MQev MaxQuant "evidence.txt" file.
#' @param LabelFactor Defaults to one treatment and one control, each with a labelled counterpart and triplicated. Needs to be defined as a factor and in the same order as the Treatment vector.
#' @param EnrichmentFileDir Parent directory were the enrichment file is contained.
#' @param ProtPTMs Character vector with the PTM codes that are part of peptides in the provided dataset.
#' @param SeqCharLength Defaults to 70. Defines the length of each line in the output protein sequences.
#' @param GroupPeptides Defaults to FALSE. Allows to group all peptides from a single protein into mean protein enrichment.
#' @param verbose Defaults to FALSE. When TRUE allows to visualize the progression through the function.
#' @param CorrectLab Defaults to TRUE. Corrects labelling percentages using the mean residual "labelling" in non-labelled samples. This ensures that only noise "labelling" percentages with low standard deviation are successfully moved to zero.
#' @keywords LCMS kProteomics Tracer Dynamics
#' @export
#' @examples
#' 
#' ...


AnnotateProteins <- function(PeptideVector,
                             Treatment = as.factor(c(rep("Control_NL",3),
                                                     rep("Control_L",3),
                                                     rep("Treatment_NL",3),
                                                     rep("Treatment_L",3))),
                             FileName = "test",
                             Path2FASTA = paste0(system.file("extdata", package = "ProtSynthesis"), "/", "160517_Hv_IBSC_PGSB_r1_proteins_HighConf_REPR_annotation.fasta"),
                             Path2MQev = paste0(system.file("extdata", package = "ProtSynthesis"), "/", "evidence.txt"),
                             LabelFactor = as.factor(rep(c(rep("Control", 3),
                                                           rep("Labelled", 3)),2)),
                             EnrichmentFileDir,
                             cexSeq = 1,
                             ProtPTMs = c("OX", "AC"),
                             SeqCharLength = 70,
                             GroupPeptides = FALSE,
                             verbose = FALSE,
                             CorrectLab = TRUE) {

    ## functions needed
  
  ProteinSequencePlot <- function(protein_sequence,
                                  lineCharLen = 50,
                                  peptide2color,
                                  cex_prot_seq = 2,
                                  PTMs = c("OX", "AC")){
    
    split_protein_sequence <- strsplit(protein_sequence, "")[[1]]
    
    protein_length <- length(split_protein_sequence)
    
    if (length(strsplit(peptide2color, "-")[[1]]) > 1){
      
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
      
    }
    
    peptide_location <- stringr::str_locate(pattern = tolower(peptide2color), string = protein_sequence)
    
    sequence_list <- list(sequence_before_pep = split_protein_sequence[1:(peptide_location[,"start"]-1)],
                          peptide_itself = split_protein_sequence[peptide_location[,"start"]:peptide_location[,"end"]],
                          sequence_after_pep = split_protein_sequence[(peptide_location[,"end"]+1):protein_length])
    
    Fragment_out_list <- list()
    
    for (j in 1:length(sequence_list)) {
      
      Fragment_sequence <- sequence_list[[j]]
      
      Fragment_length <- length(Fragment_sequence)
      
      lines_Nr <- ceiling(Fragment_length/lineCharLen)
      
      NewProtSeq <- c()
      
      for (i in 1:lines_Nr) {
        
        if (i < lines_Nr) {
          
          ## all except last line (which contains less characters)
          
          Char_Nr <- ((i-1)*lineCharLen):(i*lineCharLen)
          
          NewProtSeq <- c(NewProtSeq, c(Fragment_sequence[Char_Nr], "\n"))
          
        } else {
          
          ## last line
          
          Char_Nr <- ((i-1)*lineCharLen):Fragment_length
          
          NewProtSeq <- c(NewProtSeq, c(Fragment_sequence[Char_Nr]))
          
        }
      }
      
      Fragment_out_list[[j]] <- NewProtSeq
      
    }
    
    ## sequence plot
    
    plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n", xaxt = "n", yaxt = "n");
    unikn::mark(x = 0, y = 8.5, labels = c(paste0(Fragment_out_list[[1]], collapse = ""),
                                           paste0(Fragment_out_list[[2]], collapse = ""), 
                                           paste0(Fragment_out_list[[3]], collapse = "")),
                col_bg = "white",
                col = c("black", "red", "black"),
                y_layout = "flush", cex = cex_prot_seq)
    
  }
  
  ## main
  
  ### Data needed
  
  #### EnrichmentFile
  
  EnrichmentFile <- read.csv(file = EnrichmentFileDir, header = T, row.names = 1)
  
  EnrichmentFile <- as.matrix(EnrichmentFile)
  
  if(CorrectLab == T) {
    
    ### correction of labelled samples with mean residual enrichment in non-labelled counterparts
    
    Treatment_order <- as.character(aggregate(as.numeric(EnrichmentFile[1,]) ~ Treatment, FUN = max)[,1])
    
    Pep_mean <- t(apply(X = EnrichmentFile,
                        MARGIN = 1,
                        function(x){as.numeric(aggregate(as.numeric(x) ~ Treatment, FUN = mean)[,2])}))
    
    colnames(Pep_mean) <- Treatment_order
    
    correction_factor_index_NL <- grep(paste0(unique(Treatment[grep("Control", LabelFactor)]), collapse = "|"), colnames(Pep_mean), value = T)
    
    correction_factor_index <- apply(list2DF(strsplit(as.character(Treatment), "_"))[1,], 2, as.character)
    
    for (i in 1:length(unique(correction_factor_index))) {
      
      index_samples_i <- which(correction_factor_index == unique(correction_factor_index)[i]) #split(grep("Labelled", LabelFactor), f = as.character(Treatment[grep("Labelled", LabelFactor)]))[i]
      
      Matrix2substractFrom <- EnrichmentFile[,index_samples_i]
      
      SubtractionMatrix <- t(tcrossprod(rep(1,length(index_samples_i)), Pep_mean[,correction_factor_index_NL[i]]))
      
      EnrichmentFile[,index_samples_i] <- Matrix2substractFrom-SubtractionMatrix
      
    }
    
    EnrichmentFile[which(EnrichmentFile < 0)] <- 0
    
  }
  
  #### FASTA
  
  Hv_FASTAs <- seqinr::read.fasta(Path2FASTA)
  
  #### Evidence file col "Leading.proteins"
  
  evidence_file <- read.table(Path2MQev, header = T, sep = "\t")
  
  ### modifying peptide entries from the evidence file to make the search of modified peptides compatible
  
  evidence_file[,"Modified.sequence"] <- as.character(evidence_file$Modified.sequence)
  
  subset_PTMs <- evidence_file[grep(evidence_file$Modified.sequence, pattern = "\\)"),]
  
  peptide_i_tr <- c()
  
  for (i in 1:dim(subset_PTMs)[1]) {
    
    peptide_runner <- as.character(subset_PTMs$Modified.sequence[i])
    
    peptide_i_tr[i] <- paste0(replace(strsplit(peptide_runner, "")[[1]],
                                      c(which(strsplit(peptide_runner, "")[[1]] == "("),
                                        which(strsplit(peptide_runner, "")[[1]] == ")")),
                                      "-"),
                              collapse = "")
    
  }
  
  rep_df <- cbind(peptide_i_tr, as.character(subset_PTMs$Modified.sequence))
  
  ### replacing the transformed entries back into the evidence file
  
  if(verbose == T){
    
    cat("\n")
    cat("\n")
    cat(paste0("\n", "Working on the evidence file to match peptides, match: ", "...\n"))
    
  }
  
  
  count = 0
  
  for (i in 1:dim(rep_df)[1]) {
    
    count = count + 1
    
    if(verbose == T){
      
      cat(paste0("\n", count, " of ", dim(rep_df)[1], "...\n"))
      
    }
    
    test_entry <- grep(rep_df[i,2], evidence_file$Modified.sequence, fixed = T)
    
    if (length(test_entry) > 0) {
      
      evidence_file[grep(rep_df[i,2],
                         evidence_file$Modified.sequence,
                         fixed = T), "Modified.sequence"] <- rep_df[i,1]
      
    }
  }
  
  ### matching peptides and protein IDs, this procedure works for a single named peptide vector fed into the function....
  
  if (verbose == T) {
    
    cat("\n")
    cat("\n")
    cat("\n")
    cat("Matching peptides with protein FASTA IDs...")
    cat("...\n")
    
  }
  
  names <- c()
  
  list_out <- list()
  
  test <- EnrichmentFile[PeptideVector,]
  
  for (i in 1:length(rownames(test))) {
    
    if (verbose == T) {
      
      cat("\n")
      cat("running on peptide #", i, ": ", rownames(test)[i])
      cat("...\n")
      
    }
    
    if (length(unique(rownames(test))) == length(rownames(test))) {
      
      peptide_i <- paste0("_", strsplit(rownames(test)[i], "_")[[1]][3],"_")
      
      runner_singlePep <- evidence_file[grep(peptide_i,
                                             evidence_file$Modified.sequence,
                                             ignore.case = T),]
      
      #### protein alone will cause agglomerations in the list
      
      identified_proteins <- strsplit(as.character(unique(runner_singlePep$Leading.razor.protein)),
                                      ";")[[1]] ##### we take the leading razor protein as the one identifying the peptide...
      
      if (length(identified_proteins) > 1) {
        
        ##### multiple prots identified
        
        for (j in 1:length(identified_proteins)) {
          
          ## look for FASTA ID
          
          ### resulting object is class "SeqFastadna"
          
          matches <- list(Hv_FASTAs[grep(identified_proteins[j],
                                         names(Hv_FASTAs))][[1]])
          
          list_out <- c(list_out, matches)
          
        }
        
        names <- c(names,
                   rep(peptide_i, length(identified_proteins)))
        
      } else {
        
        ##### one prot identified
        
        ##### look for FASTA ID
        
        ###### resulting object is class "SeqFastadna"
        
        matches <- list(Hv_FASTAs[grep(identified_proteins, names(Hv_FASTAs))][[1]])
        
        list_out <- c(list_out, matches)
        
        names <- c(names, peptide_i)
        
      }
      
    } else {
      
      stop("ERROR: Please remove repeated peptides in your file")
      
    }
  }
  
  names(list_out) <- names
  
  ### take objects apart to build a data frame of annotations in order to recognize which proteins are there
  
  
  if (verbose == T) {
    
    cat("\n")
    cat("\n")
    cat("\n")
    cat("Building annotation matrix and graphical outputs...")
    cat("...\n")
    
  }
  
  out_mat <- matrix(data = NA, nrow = length(list_out), ncol = 4)
  
  for (i in 1:length(list_out)) {
    
    Sequence <- paste0(list_out[[i]][1:length(list_out[[i]])], collapse = "")
    
    FASTA <- attr(list_out[[i]], "Annot")
    
    PeptideID <- names(list_out)[i]
    
    HORVU <- strsplit(FASTA, "\\.")[[1]][1]
    
    out_mat[i,] <- c(Sequence, FASTA, PeptideID, HORVU)
    
  }
  
  colnames(out_mat) <- c("Sequence", "FASTA", "PeptideID", "HORVU")
  
  #### removing fully duplicated entries
  
  out_mat <- out_mat[-c(which(duplicated(x = as.data.frame(out_mat)))),]
  
  ### getting enrichment (%), sd and sequence coverage for proteins
  
  #### loop through out_mat getting peptides from common HORVU codes into a single entry in order to calculate enrichment, sd values and highlight sequence coverage per protein coding gene.
  
  #### Two outputs; 1- graphical output & 2- table with stats
  
  pdf(paste0(FileName, ".pdf"))
  
  if(GroupPeptides == T) {
   
    Vector2Loop <- unique(out_mat[,"HORVU"]) 
    
  } else {
    
    Vector2Loop <- unique(out_mat[,"PeptideID"]) 
    
  }
  
  ### the errors are in the information compilation of this loop now..........
  
  out_mat_pep <- matrix(NA, nrow = 0, ncol = (dim(out_mat)[2]+1+length(Treatment)+(length(levels(Treatment))*2))) ##### FASTA info + pep ID + Enrichments + means + St.Devs
  
  for (i in 1:length(Vector2Loop)) {
    
    if(GroupPeptides == T) {
      
      if (verbose == T) {
        
        cat("\n")
        cat("running on protein #", i, ": ", Vector2Loop[i])
        cat("...\n")
        
      }
      
      runner_mat <- as.matrix(out_mat[which(Vector2Loop[i] == out_mat[,"HORVU"]),])
      
    } else {
      
      if (verbose == T) {
        
        cat("\n")
        cat("running on peptide #", i, ": ", Vector2Loop[i])
        cat("...\n")
        
      }
      
      runner_mat <- as.matrix(out_mat[which(Vector2Loop[i] == out_mat[,"PeptideID"]),])
    }
    
    if (dim(runner_mat)[2] == 1) {
      
      runner_mat = t(runner_mat)
      
      ##### only one peptide per protein (one plot)
      
      ##### main pep runner table
      
      Pep_ernichments <- as.matrix(test[grep(pattern = as.character(runner_mat[,"PeptideID"]), x = rownames(test)),])
      
      if (dim(Pep_ernichments)[2]==1) {
        
        Pep_ernichments <- t(Pep_ernichments)
        
      }
      
      ##### meta information from main table derived for out_table
      
      pep_IDs <- paste0(rownames(test)[grep(pattern = as.character(runner_mat[,"PeptideID"]), x = rownames(test))], collapse = ";")
      names(pep_IDs) <- "Peptide_MetaInfo"
      
      full_enrichments <- colMeans(Pep_ernichments)
      
      means_pep <- aggregate(colMeans(Pep_ernichments) ~ Treatment, FUN = mean)[,2]
      names(means_pep) <- paste0("Mean_", aggregate(colMeans(Pep_ernichments) ~ Treatment, FUN = mean)[,1])
      
      Sds_pep <- aggregate(colMeans(Pep_ernichments) ~ Treatment, FUN = sd)[,2]
      names(Sds_pep) <- paste0("Sd_", aggregate(colMeans(Pep_ernichments) ~ Treatment, FUN = sd)[,1])
      
      out_vec <- c(pep_IDs, full_enrichments, means_pep, Sds_pep)
      
      ##### boxplots (Tukey HSD boxplots)
      
      par(mar = c(10,10,10,10))
      
      boxplot((colMeans(Pep_ernichments)*100) ~ Treatment,
              xlab = NULL, ylab = "non-corrected LPF (%)",
              main = c(strsplit(runner_mat[,"FASTA"], "\\|")[[1]][4],
                       as.character(runner_mat[,"PeptideID"])), las = 2);
      stripchart((colMeans(Pep_ernichments)*100) ~ Treatment,
                 method = "jitter", pch = 19, col = 1,
                 vertical = TRUE, add = TRUE) 
      
      ##### protein sequence peptide highlights plot
      
      protein_sequence <- as.character(unique(runner_mat[, "Sequence"]))
      
      Peptide2Highlight <- strsplit(runner_mat[,"PeptideID"], "_")[[1]][2]
      
      par(mar = c(0,0,0,0))
      
      ProteinSequencePlot(protein_sequence = protein_sequence,
                          lineCharLen = SeqCharLength, 
                          peptide2color = Peptide2Highlight, 
                          cex_prot_seq = cexSeq,
                          PTMs = ProtPTMs)
      
      if (dim(Pep_ernichments)[2] > 1) {
        ##### multiple peptides from the same species (different charges); information for out_table
        Peptide_Nr <- dim(Pep_ernichments)[2]
        
      } else {
        ##### single charged peptide species; information for out_table
        Peptide_Nr <- 1
        
      }
      ##### Paste new information in out_table (means, sds, multiple or single charges peptide Nr)
      
      out_runner_mat_vec <- c(runner_mat[1,], out_vec) ## this vector is to be pasted in out_mat...
      
      out_mat_pep <- rbind(out_mat_pep, out_runner_mat_vec)
      
    } else {
      
      ##### multiple peptides per protein (multiple plots)
      
      iter_pep <- runner_mat[,"PeptideID"]
      
      for (j in 1:length(iter_pep)) {
        
        ##### procedures for single peptides
        
        ##### main pep runner table
        
        Pep_ernichments <- as.matrix(test[grep(pattern = iter_pep[j], x = rownames(test)),])
        
        pep_IDs <- paste0(rownames(Pep_ernichments), collapse = ";")
        names(pep_IDs) <- "Peptide_MetaInfo"
        
        full_enrichments <- colMeans(Pep_ernichments)
        
        means_pep <- aggregate(colMeans(Pep_ernichments) ~ Treatment, FUN = mean)[,2]
        names(means_pep) <- paste0("Mean_", aggregate(colMeans(Pep_ernichments) ~ Treatment, FUN = mean)[,1])
        
        Sds_pep <- aggregate(colMeans(Pep_ernichments) ~ Treatment, FUN = sd)[,2]
        names(Sds_pep) <- paste0("Sd_", aggregate(colMeans(Pep_ernichments) ~ Treatment, FUN = sd)[,1])
        
        out_vec <- c(pep_IDs, full_enrichments, means_pep, Sds_pep)
        
        ##### boxplots (Tukey HSD boxplots)
        
        par(mar = c(10,10,10,10))
        
        boxplot((colMeans(Pep_ernichments)*100) ~ Treatment,
                xlab = NULL, ylab = "non-corrected LPF (%)",
                main = c(strsplit(runner_mat[,"FASTA"], "\\|")[[1]][4],
                         iter_pep[j]), las = 2);
        stripchart((colMeans(Pep_ernichments)*100) ~ Treatment,
                   method = "jitter", pch = 19, col = 1,
                   vertical = TRUE, add = TRUE) 
        
        ##### protein sequence peptide highlights plot
        
        protein_sequence <- as.character(unique(runner_mat[, "Sequence"]))
        
        Peptide2Highlight <- strsplit(iter_pep[j], "_")[[1]][2]
        
        par(mar = c(0,0,0,0))
        
        ProteinSequencePlot(protein_sequence = protein_sequence,
                            lineCharLen = SeqCharLength, 
                            peptide2color = Peptide2Highlight, 
                            cex_prot_seq = cexSeq,
                            PTMs = ProtPTMs)
        
        if (dim(Pep_ernichments)[1] > 1) {
          ##### multiple peptides from the same species (different charges); information for out_table
          Peptide_Nr <- dim(Pep_ernichments)[1]
          
        } else {
          ##### single charged peptide species; information for out_table
          Peptide_Nr <- 1
          
        }
        ##### Paste new information in out_table (means, sds, multiple or single charges peptide Nr)
        
        out_runner_mat_vec <- c(runner_mat[j,], out_vec) ## this vector is to be pasted in out_mat...
        
        out_mat_pep <- rbind(out_mat_pep, out_runner_mat_vec)
        
      }
    }
  }
  
  rownames(out_mat_pep) <- NULL
  
  dev.off()
  
  return(out_mat_pep)
  
}

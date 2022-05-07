#' A function to optimize peptide isotopolog intensities and transform them into IsoCorrectoR input files
#'
#' This function allows users to grab the output "intensities.tsv" table from isotopeEnrichment.py and transform it into the correct input format for IsoCorrectoR while selecting and optimal number of isotopolog peaks per peptide entry. Optimization of the isotopolog number relies on estimating first how many isotopologs can be expected from the atomic composition of each peptide and their natural isotopic abundance and secondly from the labelling percentage in soluble amino acid pools and the number of labelled amino acid residues in each peptide sequence.
#' @param PyResultsDir Parent directory were the XX file is contained.
#' @param returnCSV Deafults to TRUE. returns the needed IsoCorrectoR MeasurementFile.csv and MoleculeFile.csv to the working directory if set to TRUE.
#' @param verbose Defaults to FALSE. Allows to monitor the progression through peptides in order to spot mistakes in the input files if any. 
#' @param rmIsotopologs Defaults to 0, as an example 3 will keep 3 isotopolog peaks per peptide and discard all others. Allows users to keep a specific number of isotopolog peaks from all input peptides if beyond this number all other peaks should be detrimental to data analysis, e.g., noisy mass spectra.
#' @param OptimizeIsotopologNr Defaults to FALSE. Parameter to replace the rmIsotopologs parameter by a customization that calculates the optimal number of isotopolog peaks per peptide taking into account its molecular formula and the enrichment percentage in soluble amino acid precursors. To use it in TRUE mode all the information need to be supplied in the subsequent parameters.
#' @param files2correct Defaults to NULL. Needs to be used if OptimizeIsotopologNr is set to TRUE. List of file directories, where each file contains the enrichments in soluble amino acids per treatment evaluated.
#' @param AA4correction Defaults to NULL. Needs to be used if OptimizeIsotopologNr is set to TRUE. Character vector with the names of the amino acids to be used for the calculations. The amino acids must be present in the files from the previous parameter.
#' @param AAinterprtFileDir File directory to the interpretation file that contains in one column amino acid names and in the second colum the single letter code for those amino acids.
#' @param ElementalFileDir Directory to the elemental file to be used in IsoCorrectoR, which contains the natural isotopic abundance (NIA) of chemical elements.
#' @param ProtPTMs Character vector with the PTM codes that are part of peptides in the provided dataset.
#' @param LabelledSamplesNr Integer reflecting the total number of labelled samples in the dataset.
#' @keywords LCMS kProteomics Tracer Dynamics
#' @export
#' @examples
#' 
#' ...

isotopeEnrichment <- function(PyResultsDir = paste0(system.file("extdata", package = "ProtSynthesis"), "/", "intensities.tsv"),
                              returnCSV = TRUE,
                              verbose = FALSE,
                              rmIsotopologs = 0,
                              OptimizeIsotopologNr = FALSE,
                              files2correct = list.files(path = system.file("extdata", package = "ProtSynthesis"), pattern = "Mean_Norm_Factor"),
                              AA4correction = c("Serine", "Glycine"),
                              AAinterprtFileDir = paste0(system.file("extdata", package = "ProtSynthesis"), "/", "AminoAcidNames2SingleLetters.csv"),
                              ElementalFileDir = paste0(system.file("extdata", package = "ProtSynthesis"), "/", "ElementFile.csv"),
                              ProtPTMs = c("OX", "AC"),
                              LabelledSamplesNr = 6) {
  
  ## Functions needed
  ### list to data.frame
  
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  }
  
  ### Isotopolog number optimization function
  
  IsotopologNumberCalculation <- function(files2correct,
                                          AA4correction,
                                          AAinterprtFileDir,
                                          singlesingle_pep,
                                          ElementalFileDir,
                                          ProtPTMs){
    
    #### ceiling(((max(%EnrSolS) * nS) + (max(%EnrSolG) * nG)) + nC[atoms] * 0.01[13C-NIA%]) ## do the generalization to include all elements of the molecular formula
    
    ##### getting max(%EnrSolS) & max(%EnrSolG)
    
    out_correctedAA_list <- list()
    
    list_names_AAs <- list()
    
    names_AAs <- c()
    
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
      
      list_names_AAs[[i]] <- rownames(file_runner_corr_adjusted)
      
      names_AAs <- c(names_AAs, rownames(file_runner_corr_adjusted))
      
      out_correctedAA_list[[i]] <- file_runner_corr_adjusted
    }
    
    ### have to input the final amino acid enrichments instead of the linear modelled data because we have proof of that being an underestimation.
    
    if(unique(Reduce(f = union, x = list_names_AAs) == unique(names_AAs)) == T) {
      
      names_AAs <- unique(names_AAs)
      
      MAX_AA_j_i <- matrix(NA, nrow = length(names_AAs), ncol = length(out_correctedAA_list))
      
      names_out_treat <- c()
      
      for (i in 1:length(out_correctedAA_list)) {
        
        runner_mat <- out_correctedAA_list[[i]]
        
        names_out_treat[i] <- paste0(colnames(out_correctedAA_list[[i]]), collapse = " ")
        
        MAX_AA_j <- c()
        
        for (j in 1:length(names_AAs)) {
          
          MAX_AA_j[j] <- max(max(runner_mat[names_AAs[j],]))
          
        }
        
        names(MAX_AA_j) <- names_AAs
        
        MAX_AA_j_i[,i] <- MAX_AA_j
        
      }
      
      rownames(MAX_AA_j_i) <- names_AAs
      colnames(MAX_AA_j_i) <- names_out_treat
      
      
      AA_max_4_formula <- matrixStats::rowMaxs(MAX_AA_j_i)
      
      names(AA_max_4_formula) <- rownames(MAX_AA_j_i)
      
    } else {
      
      stop("ERROR: The function only works using the same amino acid set for all treatments... Please modify input files accordingly")
      
    }
    
    #### getting nAA_i atoms
    
    pep_sequence_PTM <- unique(singlesingle_pep$Sequence)
    
    ##### removed PTMs
    
    peptide2color <- pep_sequence_PTM
    
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
      
      pep_sequence <- strsplit(peptide2color, "")[[1]]
      
    } else {
      
      pep_sequence <- strsplit(peptide2color, "")[[1]]
      
    }
    
    ##### pep_sequence is already stripped off PTMs - ready to match enriched AAs
    
    Sum_AAs_lab_percentages <- c()
    
    for (i in 1:length(AA_max_4_formula)) {
      
      Sum_AAs_lab_percentages[i] <- length(grep(pattern = names(AA_max_4_formula)[i], pep_sequence)) * AA_max_4_formula[i]
      
    }
    
    nAA_isotopolog_Nr <- ceiling(sum(Sum_AAs_lab_percentages))
    
    #### getting ceiling(nX[atoms] * %[X-NIA])) according to molecular formula
    
    ##### molecular formula
    
    input_formula <- strsplit(unique(singlesingle_pep$Formula), " ")[[1]]
    
    ##### getting the number of atoms in each element
    
    Nr_of_Atoms <- as.numeric(gsub("[[:alpha:]]","",input_formula))
    
    ##### getting the NIA abundances from the atoms present in the molecular formula
    
    ElementalFile <- read.csv(ElementalFileDir)
    
    Letters <- gsub("[[:digit:]]","",input_formula)
    
    NIA_atoms_form <- c()
    
    for (i in 1:length(Letters)) {
      
      NIA_atoms_form[i] <- 1-as.numeric(strsplit(grep(x = strsplit(ElementalFile$Isotope.abundance_Mass.shift[ElementalFile$Element == Letters[i]],
                                                                   "/")[[1]],
                                                      pattern = "_0", value = T), "_0")[[1]])
      
    }
    
    ##### getting the number of isotopolog
    
    iso_number_NIAbased <- ceiling(sum(NIA_atoms_form))
    
    #### final sum that takes into account labelling in sol AAs and NIA to determine the optimal isotopolog number for the correction:
    
    isotopolog_number <- nAA_isotopolog_Nr + iso_number_NIAbased
    
    return(isotopolog_number)
    
  }
  
  
  ### get isotopologues function
  
  get_isotopologues <- function(singlesingle_pep,
                                Treatments,
                                files2correct,
                                AA4correction,
                                AAinterprtFileDir,
                                ElementalFileDir, 
                                ProtPTMs,
                                OptimizeIsotopologNr,
                                LabelledSamplesNr) {
    
    TreatmentNo = length(Treatments)
    
    dim_test <- na.omit(t(singlesingle_pep[1,12:dim(singlesingle_pep)[2]]))
    
    out_mat <- matrix(NA, nrow = dim(dim_test)[1], ncol = TreatmentNo)
    
    colnames(out_mat) <- Treatments
    
    IsoIDs <- c()
    
    count = 0
    
    for (i in 1:dim(singlesingle_pep)[1]) {
      
      count = count + 1
      #print(count)
      
      df_single_row <- na.omit(t(singlesingle_pep[i,12:dim(singlesingle_pep)[2]]))
      
      alocation <- which(singlesingle_pep$File[i] == colnames(out_mat))
      
      out_mat[,alocation] <- as.numeric(df_single_row[1:dim(dim_test)[1],1])
      
    }
    
    IsoIDs <- paste0(rep(as.character(unique(singlesingle_pep$Peptide.Key)), dim(dim_test)[1]),
                     "_",
                     apply(list2df(strsplit(rownames(df_single_row), "\\.")), 2, as.character)[2,])
    
    return_mat <- cbind("Measurements/Samples" = IsoIDs, out_mat)
    
    #### introducing and using the isotopolog number calculation
    
    if(OptimizeIsotopologNr == T){
      
      Isotopologs2return <- IsotopologNumberCalculation(files2correct = files2correct,
                                                        AA4correction = AA4correction,
                                                        AAinterprtFileDir = AAinterprtFileDir,
                                                        singlesingle_pep = singlesingle_pep,
                                                        ElementalFileDir = ElementalFileDir, 
                                                        ProtPTMs = ProtPTMs)
      
      return_mat = return_mat[1:(Isotopologs2return+1),]
      
      ### getting rid of rows with more zeros than non-labelled samples, assumes a balanced design
      
      return_mat <- return_mat[which(apply(return_mat, 1, function(x){length(which(as.numeric(x[2:length(x)]) == 0))}) <= LabelledSamplesNr),]
      
    }
    
    return(return_mat)
    
  }
  
  ## main
  ### get data
  
  data <- read.table(file = PyResultsDir, header = T, sep = "\t", fill = T)
  
  ### removing peptides without Special.Residues 
  #### (these peptides can be tested for non-specific enrichment from other amino acids)
  
  rm_pep <- which(data$Special.Residue.Count == 0)
  
  if(length(rm_pep) > 0) {
    
    data = data[-c(rm_pep),]  
    
  }
  
  rownames(data) <- paste0(rep("X", length(rownames(data))), rownames(data))
  
  #### Treatments
  
  Treatments <- as.character(unique(data$File))
  
  test_Treatments <- c()
  
  for (i in 1:length(Treatments)) {
    
    test_Treatments[i] <- length(strsplit(Treatments, split = "")[[i]]) > 0
  }
  
  Treatments = Treatments[grep(T, test_Treatments)]
  
  TreatmentNo = length(Treatments)
  
  #### Measurement file
  
  ##### First step must be subset the matrix to group of identical peptides in bins
  
  peptides <- unique(data$Peptide.Key)
  
  cat("...", "\n")
  cat(paste0(length(peptides), " Peptides found"))
  cat("...", "\n")
  
  MeasurementFile <- matrix(NA, nrow = 0, ncol = (TreatmentNo + 1))
  
  Molecule <- c()
  
  for (i in 1:length(peptides)) {
    
    test_sequence <- length(strsplit(as.character(strsplit(as.character(peptides[i]),
                                                           "_")[[1]][1]),
                                     "")[[1]])
    
    if (test_sequence > 0){
      
      ###### grabing each peptide
      
      single_pep_matches <- which(peptides[i] == data$Peptide.Key)
      
      single_pep <- data[single_pep_matches,]
      
      ###### filtering out repeated peptide peaks
      
      single_pep = single_pep[!duplicated(single_pep$File),]
      
      ###### grabbing identical peptides but differently charged separately
      
      test_charges <- unique(single_pep$Charge)
      
      if (verbose == T) {
        
        cat("...", "\n")
        cat(paste0("Peptide # ", i, " : ", peptides[i]))
        cat("...", "\n")
        
      }
      
      if (length(test_charges) > 1){
        
        MeasurMat <- matrix(NA, nrow = 0, ncol = (TreatmentNo + 1))
        
        for (j in 1:length(test_charges)) {
          
          singlesingle_pep <- single_pep[grep(test_charges[j], single_pep$Charge),]
          
          ###### grabbing molecular formula and potentially labelled residue number
          
          runner <- paste0(gsub(" ", x = unique(singlesingle_pep$Formula),
                                replacement = "", fixed = T),
                           "LabN",
                           unique(singlesingle_pep$Special.Residue.Count),
                           collapse = "")
          
          Molecule <- c(Molecule, runner)
          
          ###### grabbing the isotopologues
          
          MeasurMat <- rbind(MeasurMat,
                             get_isotopologues(singlesingle_pep,
                                               Treatments = as.character(unique(data$File)[1:TreatmentNo]),
                                               files2correct = files2correct,
                                               AA4correction = AA4correction,
                                               AAinterprtFileDir = AAinterprtFileDir,
                                               ElementalFileDir = ElementalFileDir, 
                                               ProtPTMs = ProtPTMs,
                                               OptimizeIsotopologNr = OptimizeIsotopologNr,
                                               LabelledSamplesNr = LabelledSamplesNr))  
        }
        
      } else {
        
        ###### peptides with only one charge
        
        singlesingle_pep <- single_pep
        
        ###### grabbing molecular formula and potentially labelled residue number
        
        runner <- paste0(gsub(" ", x = unique(singlesingle_pep$Formula),
                              replacement = "", fixed = T),
                         "LabN",
                         unique(singlesingle_pep$Special.Residue.Count),
                         collapse = "")
        
        Molecule <- c(Molecule, runner)
        
        ###### grabbing the isotopologues
        
        MeasurMat <- get_isotopologues(singlesingle_pep = singlesingle_pep,
                                       Treatments = as.character(unique(data$File)[1:TreatmentNo]),
                                       files2correct = files2correct,
                                       AA4correction = AA4correction,
                                       AAinterprtFileDir = AAinterprtFileDir,
                                       ElementalFileDir = ElementalFileDir, 
                                       ProtPTMs = ProtPTMs,
                                       OptimizeIsotopologNr = OptimizeIsotopologNr,
                                       LabelledSamplesNr = LabelledSamplesNr)
        
      } 
    }
    
    MeasurementFile <- rbind(MeasurementFile, MeasurMat)
  }
  
  ###### identifying and dealing with duplicated peptides 
  ###### (Enrichment does not allow duplicated names)
  
  dpl <- which(duplicated(x = MeasurementFile, MARGIN = 1) == T)
  
  if (length(dpl) > 0) {
    
    MeasurementFile = MeasurementFile[-c(dpl),]
    
  } else {
    
    cat("\n", "No Duplicated peptides found")
    
  }
  
  MeasurementFile[,1] = stringr::str_replace_all(MeasurementFile[,1], pattern = "\\(", replacement = "-")
  MeasurementFile[,1] = stringr::str_replace_all(MeasurementFile[,1], pattern = "\\)", replacement = "-")
  
  #### Molecule file
  
  Molecule_names <- list2df(strsplit(MeasurementFile[,1], "_"))
  
  Sequence_part <- apply(Molecule_names[3,], 2, as.character)
  
  charges_part <- apply(Molecule_names[5,], 2, as.character)
  
  specialResidue_part <- apply(Molecule_names[1,], 2, as.character)
  
  MoleculeFile <- as.data.frame(cbind(Molecule = unique(paste0(specialResidue_part,
                                                               "__",
                                                               Sequence_part,
                                                               "__", 
                                                               charges_part)),
                                      "MS ion or MS/MS product ion" = Molecule,
                                      "MS/MS neutral loss" = NA))
  
  ###### Removing any Isotopolog peaks ??
  
  if(rmIsotopologs > 0) {
    
    rm_big_isotopologs <- MeasurementFile
    
    list_of_all <- unique(apply(list2df(strsplit(rm_big_isotopologs[,1], "_"))[6,], MARGIN = 2, FUN = as.character))
    
    iter_vector2rm <- as.numeric(list_of_all[which(as.numeric(list_of_all) > rmIsotopologs)])
    
    for (i in 1:length(iter_vector2rm)) {
      
      cat("Removing Isotopolog ", iter_vector2rm[i], " of ", iter_vector2rm[length(iter_vector2rm)], "\n...")
      
      rm_big_isotopologs = rm_big_isotopologs[-c(which(apply(list2df(strsplit(rm_big_isotopologs[,1], "_"))[6,],
                                                             MARGIN = 2, FUN = as.character) == iter_vector2rm[i])),]
      
      
    }
    
    MeasurementFile = rm_big_isotopologs
    
  }
  
  MeasurementFile <- as.matrix(MeasurementFile)
  
  MeasurementFile[which(is.na(MeasurementFile))] <- 0
  
  MeasurementFile <- as.data.frame(MeasurementFile)
  
  
  if (returnCSV == T) {
    
    write.csv(x = MoleculeFile, file = "MoleculeFile.csv", quote = F, row.names = F)
    write.csv(x = MeasurementFile, file = "MeasurementFile.csv", quote = F, row.names = F)
    
  }
  
  return(list(MeasurementFile,
              MoleculeFile))
}

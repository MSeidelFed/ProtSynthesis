#' A function to apply thresholds and test statistically peptide enrichment in order to obtain good quality labelled peptide subsets 
#'
#' The statistical filters applied in this section are meant to remove falsely interpreted isotopolog abundances derived from "labelled" controls, this phenomenon can result from peptide coelution and contamination if the heavy isotopolog peaks with different peptides. This in turn results in an increased relative isotope abundance that does not come from the labelling experiment. This special scenario exemplifies the utility of having non-labelled controls in your samples. This function allows users to apply thresholds on residual labelling, noise and multiple parameters that leverage on the quality of the information that is delivered. The function returns a list of peptide subsets that fit the selected criteria and may be fed as input to the subsequent function in order to annotate their parent protein identities in the experimental dataset.
#' @param EnrichmentFileDir Parent directory were the enrichment file is contained.
#' @param Treatment Factor vector with the containing sample information in the same order as samples are outlined in the enrichment file.
#' @param LabelFactor Defaults to one treatment and one control, each with a labelled counterpart and triplicated. Needs to be defined as a factor and in the same order as the Treatment vector.
#' @param NLMAX Defaults to 5% (0.05). Allowed boundaries for noise, "false enrichment" in non-labelled / control samples.
#' @param LEnrLack Defaults to 2% (0.02). Below the selected boundary is TRUE lack of enrichment, above is FALSE.
#' @param SigLab Defaults to 5% (0.05). Below the selected boundary the FDR corrected P values comparing labeling between control and treatment are significantly different, i.e., TRUE.
#' @param SigLabControl Defaults to 5% (0.05). Below the selected boundary the FDR corrected P values comparing non-labelled and labelled control are significantly different, i.e., TRUE.
#' @param SigLabTreatment Defaults to 5% (0.05). Below the selected boundary the FDR corrected P values comparing non-labelled and labelled treatment are significantly different, i.e., TRUE.
#' @param CorrectLab Defaults to TRUE. Corrects labelling percentages using the mean residual "labelling" in non-labelled samples. This ensures that only noise "labelling" percentages with low standard deviation are successfully moved to zero.
#' @keywords LCMS kProteomics Tracer Dynamics
#' @export
#' @examples
#' 
#' ...

EnrichmentSet <- function(EnrichmentFileDir,
                          Treatment = as.factor(c(rep("Control_NL",3),
                                                  rep("Control_L",3),
                                                  rep("Treatment_NL",3),
                                                  rep("Treatment_L",3))),
                          LabelFactor = as.factor(rep(c(rep("Control", 3),
                                                        rep("Labelled", 3)),2)),
                          NLMAX = 0.05,
                          LEnrLack = 0.02,
                          SigLab = 0.05,
                          SigLabControl = 0.05,
                          SigLabTreatment = 0.05,
                          CorrectLab = TRUE){
  
  ## Needed functions
  ### autoscale function
  
  autoscale <- function(mat = distribution_test_mat(),
                        center_fun = colMeans,
                        scale_fun = matrixStats::colSds) {
    
    scaled_mat <- scale(mat, center = center_fun(mat), scale = scale_fun(mat))
    
    return(scaled_mat)
    
  }
  
  ### produce heatmap function
  
  ProduceHeatmap <- function(EnrichmentFileMat,
                             ListOfPeptides,
                             Autoscale = T,
                             Km,
                             Treatment,
                             main){
    
    EnrichmentFile <- EnrichmentFileMat
    
    HM_mat <- as.matrix(EnrichmentFile[ListOfPeptides,])
    
    if (Autoscale == T) {
      
      HM_mat[which(HM_mat == 0)] <- abs(rnorm(n = length(which(HM_mat == 0)), mean = 0.00000000001, sd = 0.000000001))
      
      HM_mat = t(autoscale(t(HM_mat)))
      
    }
    
    test_means <- t(apply(HM_mat, 1, function(x){aggregate(x ~ Treatment, FUN = mean)[,2]}))
    
    colnames(test_means) <- levels(Treatment)
    
    type = gsub("s\\d+_", "", lapply(levels(Treatment), FUN = as.character))
    
    ha = ComplexHeatmap::HeatmapAnnotation(df = data.frame(type = type))
    
    Heatmap_type1 <- ComplexHeatmap::Heatmap(matrix = test_means,
                                             clustering_distance_columns = "pearson",
                                             clustering_method_columns = "average"  ,
                                             name = "Intensities",
                                             km = Km,
                                             col = circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
                                                                        c("darkgoldenrod1",
                                                                          "yellow",
                                                                          "white",
                                                                          "violet",
                                                                          "purple")),
                                             top_annotation = ha,
                                             show_column_names = F,
                                             show_row_names = F,
                                             cluster_columns = F,
                                             cluster_rows = T,
                                             clustering_method_rows = "average",
                                             clustering_distance_rows = "pearson",
                                             column_title = main)
    
    return(Heatmap_type1)
    
  }

  ## main
  ### generate stats (append the stats as new columns in EnrichmentFile)
  
  EnrichmentFile <- read.csv(file = EnrichmentFileDir, header = T, row.names = 1)
  
  EnrichmentFile <- as.matrix(EnrichmentFile)
  
  Treatment_order <- as.character(aggregate(as.numeric(EnrichmentFile[1,]) ~ Treatment, FUN = max)[,1])
  
  Pep_max <- t(apply(X = EnrichmentFile,
                     MARGIN = 1,
                     function(x){as.numeric(aggregate(as.numeric(x) ~ Treatment, FUN = max)[,2])}))
  
  colnames(Pep_max) <- Treatment_order
  
  Pep_mean <- t(apply(X = EnrichmentFile,
                      MARGIN = 1,
                      function(x){as.numeric(aggregate(as.numeric(x) ~ Treatment, FUN = mean)[,2])}))
  
  colnames(Pep_mean) <- Treatment_order
  
  Pep_sd <- t(apply(X = EnrichmentFile,
                    MARGIN = 1,
                    function(x){as.numeric(aggregate(as.numeric(x) ~ Treatment, FUN = sd)[,2])}))
  
  colnames(Pep_sd) <- Treatment_order
  

  ### add subtraction of non-labelled residual labeling in each treatment from NL-samples to L-samples
  
  if(CorrectLab == T) {
    
    ### correction of labelled samples with mean residual enrichment in non-labelled counterparts
    
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
  
  FALSE_TRUE_Enr_Max <- ifelse(Pep_max > NLMAX, yes = T, no = F)
  
  FALSE_TRUE_Enr_Mean <- ifelse(Pep_mean > LEnrLack, yes = T, no = F)

  
  ## RandoDiStats omics testing
  
  EnrichmentFile[which(EnrichmentFile == 0)] <- abs(rnorm(n = length(which(EnrichmentFile == 0)), mean = 0.0000000001, sd = 0.0000001))
  
  omics_tests <- OmicsUnivariateStats(class_comparison_mat = t(EnrichmentFile),
                                      Factor1 = Treatment,
                                      TukeyReturns = "MeanComparisons",
                                      ReturnTukeyPlots = T,
                                      TukeyPDFName = "Tukey_tests",
                                      returnObject = "OmicsTests")
  
  print(colnames(omics_tests))
  
  cat("..\n")
  
  LvsL <- readline(prompt = "which column compares labelling across the two desired conditions? (1, 2..., n): ")
  
  cat("..\n")
  
  L20vsCon20 <- readline(prompt = "which column compares labelled vs non-labelled states in the control? (1, 2..., n): ")
  
  cat("..\n")
  
  L4vsCon4 <- readline(prompt = "which column compares labelled vs non-labelled states in the treatment? (1, 2..., n): ")
  
  cat("..\n")
  
  
  ### extract peptide subsets once the filters are built
  
  ### Filters for all
  
  cat("..\n")
  
  print(colnames(FALSE_TRUE_Enr_Max))
  
  cat("..\n")
  
  Test_Con_NL <- readline(prompt = "which column represents the non-labelled control? (1, 2..., n): ")
  
  cat("..\n")
  
  Test_Treat_NL <- readline(prompt = "which column represents the non-labelled treatment? (1, 2..., n): ")
  
  cat("..\n")
  
  ##### NL_20°C_MAX < 0.05 in autoscaled vectors, i.e., REAL NO NOISE IN NL CONTROL SAMPLES
  
  out_test_Con_NL <- names(which(FALSE_TRUE_Enr_Max[,colnames(FALSE_TRUE_Enr_Max)[as.numeric(Test_Con_NL)]] == F))
  
  ##### NL_4°C_MAX < 0.05 in autoscaled vectors, i.e., REAL NO NOISE IN NL TREATMENT SAMPLES
  
  out_test_Treat_NL <- names(which(FALSE_TRUE_Enr_Max[,colnames(FALSE_TRUE_Enr_Max)[as.numeric(Test_Treat_NL)]] == F))
  
  #### Cold vs Control sig. labelled peptides & Control vs Cold sig. labelled peptides
  
  ##### T.test cold control (< 0.05): significant = TRUE, i.e., REAL SIGNIFICANT DIFF BETWEEN CONTROL AND TREATED SAMPLES
  
  test_lab_20vs4 <- names(ifelse(test = omics_tests[,colnames(omics_tests)[as.numeric(LvsL)]] < SigLab,
                           yes = T, no = F))
  
  cat("..\n")
  
  print(colnames(FALSE_TRUE_Enr_Mean))
  
  cat("..\n")
  
  Test_Con_L_mean <- readline(prompt = "which column represents the labelled control? (1, 2..., n): ")
  
  cat("..\n")
  
  Test_Treat_L_mean <- readline(prompt = "which column represents the labelled treatment? (1, 2..., n): ")
  
  cat("..\n")
  
  ##### L_20°C rm Lack of Enrichment: L_20°C_AV > 0.02, i.e., REAL LABELLING IN L CONTROL SAMPLES
  
  out_test_Con_L <- names(which(FALSE_TRUE_Enr_Mean[,colnames(FALSE_TRUE_Enr_Mean)[as.numeric(Test_Con_L_mean)]] == T))
  
  ##### L_4°C rm Lack of Enrichment: L_4°C_AV > 0.02, i.e., REAL LABELLING IN L treated SAMPLES
  
  out_test_Treat_L <- names(which(FALSE_TRUE_Enr_Mean[,colnames(FALSE_TRUE_Enr_Mean)[as.numeric(Test_Treat_L_mean)]] == T))
  
  
  #### Cold.Sig.Lab.Peptides
  
  ##### T.test cold L vs NL (< 0.05): significant = TRUE, i.e., REAL SIGNIFICANT DIFF BETWEEN TREATMENT L AND NL SAMPLES
  
  test_lab_4Lvs4Con <- names(which(ifelse(test = omics_tests[,colnames(omics_tests)[as.numeric(L4vsCon4)]] < SigLabTreatment,
                                          yes = T, no = F) == T, useNames = T))
  
  
  #### Con.Sig.Lab.Peptides
  
  ##### T.test control L vs NL (< 0.05): significant = TRUE, REAL SIGNIFICANT DIFF BETWEEN CONTROL L AND NL SAMPLES
  
  test_lab_20Lvs20Con <- names(which(ifelse(test = omics_tests[,colnames(omics_tests)[as.numeric(L20vsCon20)]] < SigLabControl,
                                            yes = T, no = F) == T, useNames = T))
  
  
  #### NL.peptides
  
  
  test_pep_subsets <- list(out_test_Con_NL,
                           out_test_Treat_NL,
                           out_test_Con_L,
                           out_test_Treat_L,
                           test_lab_20vs4,
                           test_lab_4Lvs4Con,
                           test_lab_20Lvs20Con)
  
  names(test_pep_subsets) <- c("Con_NL",
                               "Treat_NL",
                               "Con_L",
                               "Treat_L",
                               "SIG_Tr_L_Con_L",
                               "SIG_Tr_L_Tr_NL",
                               "SIG_Con_L_Con_NL")
  
  ### Lists of peptides and intersections
  
  #### Venn Diagrams
  
  ##### Cold vs Control sig. labelled peptides UNION Control vs Cold sig. labelled peptides
  
  VennDiagram::venn.diagram(x = test_pep_subsets[1:5],
                            filename = "Sig_Lab_Between_Conditions.png",
                            col = c("yellow", "purple","lightblue", "grey", "pink"),
                            cex = 0.4,
                            cat.cex = 0.4,
                            cat.dist = 0.3,
                            cat.col =	c("yellow", "purple","lightblue", "grey", "pink"),
                            margin = 1)
  
  Union_TrCon_SigLab <- Reduce(f = union, x = test_pep_subsets[c(1:4)])
  
  last_index <- length(test_pep_subsets)+1
  
  test_pep_subsets[[last_index]] <- Union_TrCon_SigLab
  
  names(test_pep_subsets)[last_index] <- "Union_TrCon_SigLab"
  
  ##### Cold.Sig.Lab.Peptides
  
  VennDiagram::venn.diagram(x = test_pep_subsets[c(1,2,4,5)],
                            filename = "Sig_Lab_Cold.png",
                            col = c("yellow", "purple","lightblue", "grey"),
                            cex = 0.4,
                            cat.cex = 0.4,
                            cat.dist = 0.3,
                            cat.col =	c("yellow", "purple","lightblue", "grey"),
                            margin = 1)
  
  Tr_SigLab <- Reduce(f = intersect, x = test_pep_subsets[c(1,2,4)])
  
  last_index <- length(test_pep_subsets)+1
  
  test_pep_subsets[[last_index]] <- Tr_SigLab
  
  names(test_pep_subsets)[last_index] <- "Tr_TRUE_Lab"
  
  # Reduce(f = intersect, x = list(Tr_SigLab, Union_TrCon_SigLab))
  
  ##### Con.Sig.Lab.Peptides
  
  VennDiagram::venn.diagram(x = test_pep_subsets[c(1,2,3,7)],
                            filename = "Sig_Lab_Control.png",
                            col = c("yellow", "purple","lightblue", "grey"),
                            cex = 0.4,
                            cat.cex = 0.4,
                            cat.dist = 0.3,
                            cat.col =	c("yellow", "purple","lightblue", "grey"),
                            margin = 1)
  
  Con_SigLab <- Reduce(f = intersect, x = test_pep_subsets[c(1,2,3)])
  
  last_index <- length(test_pep_subsets)+1
  
  test_pep_subsets[[last_index]] <- Con_SigLab
  
  names(test_pep_subsets)[last_index] <- "Con_TRUE_Lab"
  
  #### Heatmaps (average trends) and output lists for next function
  
  pdf("Heatmaps.pdf")
  
  cat("\n Building Heatmaps of detected peptide subsets... \n")
  
  for (i in 1:length(test_pep_subsets)) {
    
    print(ProduceHeatmap(EnrichmentFileMat = EnrichmentFile,
                         ListOfPeptides = test_pep_subsets[[i]],
                         Autoscale = T, 
                         Km = 0,
                         Treatment = Treatment,
                         main = names(test_pep_subsets)[i]))
    
  }
  
  dev.off()
  
  return(test_pep_subsets)
  
}

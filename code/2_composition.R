# R code for Qin et al., Journal of Ecology, accepted Aug. 28, 2019

# Figure 3 and related statistics on community composition.
# Effects of global change treatments, soil properties, 
# and plant community properties on community composition
# of soil microbial communities.

# For this script to work, you must first run 0_setup.R
# in the same R session.

library(dplyr)
library(ggplot2)
library(tidyr)
library(gdm)

# Create directory to store intermediate results from GDM
results_output_dir <- "../results/"
if(!dir.exists(results_output_dir)) dir.create(results_output_dir)

# Add new columns to sampleData to represent the
# treatment interactions
sampleData %>%
  mutate(nc = n==2 & co2==2,
         nt = n==2 & heat==2,
         np = n==2 & precip==2,
         ct = co2==2 & heat==2,
         cp = co2==2 & precip==2,
         tp = heat==2 & precip==2,
         nct = n==2 & co2==2 & heat==2,
         ncp = n==2 & co2==2 & precip==2,
         ntp = n==2 & heat==2 & precip==2,
         ctp = co2==2 & heat==2 & precip==2,
         nctp = n==2 & co2==2 & heat==2 & precip==2) ->
  sampleData

# Run variable importance analysis.
# WARNING: Can take a long time to run.
# First, use backwards selection to find parsimonious model: choose the
# first model whose pval is below selection_pval_threshold.
# Then run gdm.varImp with fullModelOnly=TRUE to get model sig.,
# model dev. explained, variable sig., variable imp.
# (Apparently these seeds don't actually guarantee consistent results)
varimp_gdm_selection <- function(physeq, preds, outfile_name='gdm_outfile', 
                                 include_geo=TRUE, include_plant_dissim=FALSE, 
                                 nperm=50, selection=TRUE, test_var_imp=TRUE, 
                                 selection_pval_threshold=0.05,
                                 next_test_parsim_model=FALSE, next_nperm=50) {
  # If selection is FALSE, then next_test_parsim_model is irrelevant
  if(!selection & next_test_parsim_model) {
    print("Warning: If 'selection' is FALSE, then 'next_test_parsim_model' is irrelevant because there will be no post-selection model to test.")
  }
  if(selection & !test_var_imp) {
    print("Warning: If 'selection' is TRUE, then 'test_var_imp' must be TRUE. Setting to TRUE.")
  }
  
  # Ensure all preds are in sampleData, 
  # and all required predictors are included in preds
  preds <- intersect(preds, names(sampleData))
  preds <- union(preds, c("id","lat","lon"))
  
  # Load OTU table
  bioData <- as.data.frame(t(otu_table(physeq))@.Data)
  id <- sub("B_GCE", "", rownames(bioData))
  id <- sub("F_GCE", "", id)
  bioData$id <- id
  
  # Some variables (notably %n and %c) are missing values, which can raise errors,
  # so filter out the records where they are missing
  complete_case_vars <- c("perc.n","perc.c")
  for(v in complete_case_vars) {
    if(v %in% preds) {
      complete_id <- sampleData$id[which(!is.na(sampleData[[v]]))]
      bioData <- filter(bioData, id %in% complete_id)
    }
  }
  
  # If including plant dissimilarity, load plant community table
  if(include_plant_dissim) {
    plants <- Plants.identified[rownames(Plants.identified) %in% bioData$id,]
    plants.dist <- as.data.frame(as.matrix(vegdist(plants, method="bray")))
    plants.dist %>%
      mutate(id = rownames(plants)) %>%
      dplyr::select(id, num_range("",1:144)) ->
      plants.dist
  }
  
  # Make sample data to match OTU table
  sampleData %>% 
    filter(id %in% bioData$id) %>%
    dplyr::select(preds) %>%
    mutate_all(funs(as.numeric)) %>%
    mutate_if(function(x) all(range(x)==c(1,2)), function(x) x-1) -> # shift 1-2 factor to 0-1 factor
    sampleData.sel
  
  # Formate site-pair tables. Two options depending on whether to include plant dissim.
  if(include_plant_dissim) {
    sitepair.sel <- formatsitepair(bioData=bioData, bioFormat=1, dist="bray", 
                                   abundance=TRUE, siteColumn="id", predData=sampleData.sel,
                                   XColumn="lon", YColumn="lat", distPreds=list(plants.dist))
  } else {
    sitepair.sel <- formatsitepair(bioData=bioData, bioFormat=1, dist="bray", 
                                   abundance=TRUE, siteColumn="id", predData=sampleData.sel,
                                   XColumn="lon", YColumn="lat")
  }
  
  # If selection==FALSE, then actually just return the gdm model,
  # optionally saving the varimp results on the model if test_var_imp==TRUE
  if(!selection) {
    
    # If you aren't doing selection, then you don't necessarily need to run the 
    # variable importance tests. But if test_var_imp==TRUE, then run:
    if(test_var_imp) {
      # Varimp results (list of 4 tables) are saved as "[outfile_name]_parsim.Rdata"
      if(include_geo) {
        gdm.selection0 <- gdm.varImp(sitepair.sel, parallel=TRUE, geo=TRUE, 
                                     outFile=paste0(outfile_name,"_parsim"), nPerm=nperm,
                                     fullModelOnly=TRUE)
        
      } else {
        gdm.selection0 <- gdm.varImp(sitepair.sel, parallel=TRUE, geo=FALSE, 
                                     outFile=paste0(outfile_name,"_parsim"), nPerm=nperm,
                                     fullModelOnly=TRUE)
      }
      
      # Also save contents of selection process to csv's
      write.csv(gdm.selection0[[1]], file=paste0(results_output_dir,outfile_name,"_parsim_summ.csv"))
      write.csv(gdm.selection0[[2]], file=paste0(results_output_dir,outfile_name,"_parsim_imp.csv"))
      write.csv(gdm.selection0[[3]], file=paste0(results_output_dir,outfile_name,"_parsim_sig.csv"))
      
      print(paste0("Saved parsimonious model's summary, variable importance, and variable significance to ",
                   paste0(outfile_name, c("_parsim_summ.csv","_parsim_imp.csv","_parsim_sig.csv"), collapse=", ")))
    }
    
    # Return gdm model without selection
    if(include_geo) {
      return(gdm(sitepair.sel, geo=TRUE))
    } else {
      return(gdm(sitepair.sel, geo=FALSE))
    }
    
    # If selection==TRUE, then run backwards selection
  } else {
    
    # Contents of selection process (list of 4 tables) are saved as "[outfile_name]_sel.Rdata"
    if(include_geo) {
      gdm.selection <- gdm.varImp(sitepair.sel, parallel=TRUE, geo=TRUE, 
                                  outFile=paste0(outfile_name,"_sel"), nPerm=nperm)
    } else {
      gdm.selection <- gdm.varImp(sitepair.sel, parallel=TRUE, geo=FALSE, 
                                  outFile=paste0(outfile_name,"_sel"), nPerm=nperm)
    }
    
    # Also save contents of selection process to csv's
    write.csv(gdm.selection[[1]], file=paste0(results_output_dir,outfile_name,"_sel_summ.csv"))
    write.csv(gdm.selection[[2]], file=paste0(results_output_dir,outfile_name,"_sel_imp.csv"))
    write.csv(gdm.selection[[3]], file=paste0(results_output_dir,outfile_name,"_sel_sig.csv"))
    
    print(paste0("Saved model selection results to ",
                 paste0(outfile_name, c("_sel_summ.csv","_sel_imp.csv","_sel_sig.csv"), collapse=", ")))
    
    # Select the first model whose model p-value is less than selection threshold
    # (Called the "parsimonious model")
    model_ind <- min(which(gdm.selection[[1]]["Model p-value",] < selection_pval_threshold))
    
    # If no models meet selection threshold, return(NA)
    if(is.infinite(model_ind)) {
      print(paste0("Selection process failed to fit any model: ", outfile_name))
      return(NA)
    }
    
    # Get the list of variables in the first significant model ("parsimonious model")
    model_vars <- gdm.selection[[2]][,model_ind]
    preds_parsim <- names(model_vars[!is.na(model_vars)])
    include_geo_parsim <- FALSE
    include_plant_dissim_parsim <- FALSE
    if("Geographic" %in% preds_parsim) {
      include_geo_parsim <- TRUE
      preds_parsim <- preds_parsim[-which(preds_parsim=="Geographic")]
    }
    if("matrix_1" %in% preds_parsim) {
      include_plant_dissim_parsim <- TRUE
      preds_parsim <- preds_parsim[-which(preds_parsim=="matrix_1")]
    }
    preds_parsim <- union(preds, c("id","lat","lon"))
    
    # If a ("parsimonious") model did meet selection threshold, return some value depending
    # on whether to next_test_parsim_model:
    
    # If next_test_parsim_model==TRUE, then run this script one more time on parsim model's
    # variables, and without backwards-selection process
    if(next_test_parsim_model) {
      
      if(include_geo_parsim & include_plant_dissim_parsim) {
        return(varimp_gdm_selection(physeq, preds_parsim, outfile_name, include_geo=TRUE, 
                                    include_plant_dissim=TRUE, nperm=next_nperm,
                                    selection=FALSE, next_test_parsim_model=FALSE))
      } else if(include_geo_parsim & !include_plant_dissim_parsim) {
        return(varimp_gdm_selection(physeq, preds_parsim, outfile_name, include_geo=TRUE, 
                                    include_plant_dissim=FALSE, nperm=next_nperm,
                                    selection=FALSE, next_test_parsim_model=FALSE))
      } else if(!include_geo_parsim & include_plant_dissim_parsim) {
        return(varimp_gdm_selection(physeq, preds_parsim, outfile_name, include_geo=FALSE, 
                                    include_plant_dissim=TRUE, nperm=next_nperm,
                                    selection=FALSE, next_test_parsim_model=FALSE))
      } else if(!include_geo_parsim & !include_plant_dissim_parsim) {
        return(varimp_gdm_selection(physeq, preds_parsim, outfile_name, include_geo=FALSE, 
                                    include_plant_dissim=FALSE, nperm=next_nperm,
                                    selection=FALSE, next_test_parsim_model=FALSE))
      }
      
      # If next_test_parsim_model==FALSE, simply return the gdm object of the parsimonious model
    } else {
      
      sampleData %>% 
        filter(id %in% bioData$id) %>%
        dplyr::select(preds_parsim) %>%
        mutate_all(funs(as.numeric)) %>%
        mutate_if(function(x) all(range(x)==c(1,2)), function(x) x-1) -> # shift 1-2 factor to 0-1 factor
        sampleData.parsim
      
      if(include_plant_dissim_parsim) {
        sitepair.parsim <- formatsitepair(bioData=bioData, bioFormat=1, dist="bray", 
                                          abundance=TRUE, siteColumn="id", predData=sampleData.parsim,
                                          XColumn="lon", YColumn="lat", distPreds=list(plants.dist))
      } else {
        sitepair.parsim <- formatsitepair(bioData=bioData, bioFormat=1, dist="bray", 
                                          abundance=TRUE, siteColumn="id", predData=sampleData.parsim,
                                          XColumn="lon", YColumn="lat")
      }
      
      # Return parsimonious gdm
      if(include_geo_parsim) {
        return(gdm(sitepair.parsim, geo=TRUE))
      } else {
        return(gdm(sitepair.parsim, geo=FALSE))
      }
    }
  }
}



preds_main <- c("id","lon","lat","n","co2","heat","precip","burn.2011","wildfire.2003")
preds_interact <- c("id","lon","lat","n","co2","heat","precip","burn.2011","wildfire.2003",
                    "nc","nt","np","ct","cp","tp","nct","ncp","ntp","ctp","nctp")
preds_plantsoil <- c("id","lon","lat","n","co2","heat","precip","burn.2011","wildfire.2003",
                     "perc.n", "c.n", "ph", "npp", "water.content")

# WARNING: GDM model-fitting takes a long time. For example,
# gdm.BAC.main took 22 minutes on a MacBook Air with a 2.2 GHz processor.

# Fit main models
set.seed(1) # this model should fail
gdm.ITS.main <- varimp_gdm_selection(ITS.rare, preds_main, "gdmselection_ITS_main", 
                                     include_geo=TRUE, include_plant_dissim=FALSE,
                                     nperm=50, selection=TRUE, selection_pval_threshold=0.05,
                                     next_test_parsim_model=TRUE, next_nperm=200)
set.seed(1)
gdm.BAC.main <- varimp_gdm_selection(BAC.rare, preds_main, "gdmselection_BAC_main", 
                                     include_geo=TRUE, include_plant_dissim=FALSE,
                                     nperm=50, selection=TRUE, selection_pval_threshold=0.05,
                                     next_test_parsim_model=TRUE, next_nperm=200)

# Fit interaction models
set.seed(1) # this model should fail
gdm.ITS.interact <- varimp_gdm_selection(ITS.rare, preds_interact, "gdmselection_ITS_interact", 
                                         include_geo=TRUE, include_plant_dissim=FALSE,
                                         nperm=50, selection=TRUE, selection_pval_threshold=0.05,
                                         next_test_parsim_model=TRUE, next_nperm=200)
set.seed(1)
gdm.BAC.interact <- varimp_gdm_selection(BAC.rare, preds_interact, "gdmselection_BAC_interact", 
                                         include_geo=TRUE, include_plant_dissim=FALSE,
                                         nperm=50, selection=TRUE, selection_pval_threshold=0.05,
                                         next_test_parsim_model=TRUE, next_nperm=200)

# Fit plant-soil models
set.seed(1)
gdm.ITS.plantsoil <- varimp_gdm_selection(ITS.rare, preds_plantsoil, "gdmselection_ITS_plantsoil", 
                                          include_geo=TRUE, include_plant_dissim=TRUE,
                                          nperm=50, selection=TRUE, selection_pval_threshold=0.05,
                                          next_test_parsim_model=TRUE, next_nperm=200)
set.seed(1)
gdm.BAC.plantsoil <- varimp_gdm_selection(BAC.rare, preds_plantsoil, "gdmselection_BAC_plantsoil", 
                                          include_geo=TRUE, include_plant_dissim=TRUE,
                                          nperm=50, selection=TRUE, selection_pval_threshold=0.05,
                                          next_test_parsim_model=TRUE, next_nperm=200)

# Plot models that successfully fit
plot(gdm.ITS.plantsoil)
plot(gdm.BAC.main)
plot(gdm.BAC.interact)
plot(gdm.BAC.plantsoil)


# Didn't save the outputs from varimp_gdm_selection? If you know the model's predictors,
# manually fit the model. Then you can plot them using plot().
gdm.ITS.plantsoil.preds <- c("co2","burn.2011")
gdm.ITS.plantsoil <- varimp_gdm_selection(ITS.rare, preds=gdm.ITS.plantsoil.preds,
                                          include_geo=TRUE, include_plant_dissim=TRUE,
                                          selection=FALSE, test_var_imp=FALSE)

gdm.BAC.main.preds <- c("co2","precip","wildfire.2003","n","heat","burn.2011")
gdm.BAC.main <- varimp_gdm_selection(BAC.rare, preds=gdm.BAC.main.preds,
                                     include_geo=TRUE, include_plant_dissim=FALSE,
                                     selection=FALSE, test_var_imp=FALSE)

gdm.BAC.interact.preds <- c("cp","wildfire.2003")
gdm.BAC.interact <- varimp_gdm_selection(BAC.rare, preds=gdm.BAC.interact.preds,
                                         include_geo=TRUE, include_plant_dissim=FALSE,
                                         selection=FALSE, test_var_imp=FALSE)

gdm.BAC.plantsoil.preds <- c("co2","wildfire.2003")
gdm.BAC.plantsoil <- varimp_gdm_selection(BAC.rare, preds=gdm.BAC.plantsoil.preds,
                                          include_geo=FALSE, include_plant_dissim=TRUE,
                                          selection=FALSE, test_var_imp=FALSE)

###-------------AneufinderFileFilter-----------###

#----Thomas van Ravesteyn
#----Kops Group
#----Hubrecht Institute

###------ Aneufinder File Filter Function------###
###--------------------------------------------###

AneufinderFileFilter <- function(sampleIDs){

###-----Required packages----------------------###
  require("ggrepel")
  require("AneuFinder")
  require("colorspace")
  
###------loop over all samples to generate QC data and filters------###

	message("==> Start Aneufinder File Filter <==")
	message("Selected model: ", selected.model)
	message("== Filter settings ==")
	message("- Total read count: ", filter.total.read.count)
	message("- Spikiness: ", filter.spikiness)
	message("- Bhattacharyya: ", filter.bhattacharyya)
	message("- Number of segments: ", filter.num.segments)
	message("- Too low copy numbers: ", filter.min.copy.number)
	message("- Too high copy numbers: ", filter.max.copy.number)
	message("- Perfectly diploid cells: ", filter.diploid)
	confirmed <- readline(prompt="Settings correct? (y/n)..")
  
	if (confirmed == "n"){
		message("Settings were not approved, adjust and restart")
		stop(call. = FALSE)
    } else {
    
    #start timer
    total.time <- proc.time()
    
    #set random seed
    set.seed(1)
          
    #Create a shortcut and directory where summary data will be saved
    summary.shortcut <- c(paste0(project.name, "/", "Summary_Files_", selected.model, "/"))
    dir.create(summary.shortcut, recursive = TRUE)
      
    #create empty datacollection list
		datacollection <- list()
		selected.sampleIDs <- c()
		

		###------loop over all samples to generate QC data and apply filters------###
		
		for(sample in sampleIDs) {
			  
		    #create empty vectors
			  selected.files    <- c()
			  removed.files     <- c()
			  diploid.files     <- c()
			  selected.numbers   <- c()
			  selected.karyotype	 <- c()
			  unfiltered.karyotype <- c()	
		  
		    #create a directory to store QC files and data files after filtering
			  dir.create(paste0(project.name, "/", sample ,"/", "QC_FILES_", selected.model), recursive = TRUE )
		    shortcut <- paste0(project.name, "/", sample ,"/", "QC_FILES_", selected.model ,"/", sample)
			
		    #write filter settings to .txt file  
		    writeLines(c("=== Selected model ===", c(paste0("- ", selected.model)), "",
		                 "=== Filter settings ===","", 
		                 c(paste0("- Minimum Total read count = ", filter.total.read.count, ", ", min.total.read.count)),"",
		                 c(paste0("- Maximum Spikiness = ",filter.spikiness, ", ", max.spikiness)),"",
		                 c(paste0("- Minimum Bhattacharyya = ", filter.bhattacharyya, ", ", min.bhattacharyya)), "",
		                 c(paste0("- Minimum Number of segments = ", filter.num.segments, ", ", min.num.segments)), "",
		                 c(paste0("- Minimum weighted copy number = ", filter.min.copy.number, ", ", min.copy.number)), "", 
                     c(paste0("- Maximum weighted copy number = ", filter.max.copy.number, ", ", max.copy.number)), "", 
                     c(paste0("- Diploid cell removal, n set for X chromosome = ", filter.diploid, ", ", nXchr))), 
  		               c(paste0(shortcut,"_Filter parameters.txt")))  
		    
		    #load Files and models based on selected model
			  message("== Load files for ", sample, " ==")
			  if (selected.model == "dnaCopy"){
				my.files <- unlist(dnaFiles[[sample]])
			  } else if (selected.model == "hmm") {
				my.files <- unlist(hmmFiles[[sample]])
			  } else if (selected.model == "edivisive"){
				my.files <- unlist(edivisiveFiles[[sample]])
			  } else {
				stop("Incorrect model selection")
			  }
			  models <- loadFromFiles(my.files)
			  
			  #plot Genomewide heatmap prior to filtering and save to QC file directory
			  if (set.plot.genomewide == TRUE){
			    heatmapGenomewide(models, file = paste0(shortcut,"_unfiltered_genomewide.pdf"))
			  }
			    
			  #plot aneuploidy and heterogeneity per chromosome prior to filtering and save to QC file directory
			  if (set.plot.heterogeneity == TRUE){
  			  myplot <- plotHeterogeneity(models)
  			  pdf(file = paste0(shortcut, "_unfiltered_heterogeneity.pdf"))
  			  print(myplot)
  			  dev.off()
			  }
			  
			  #Generate karyotype measures prior to filtering
			  unfiltered.karyotype <- karyotypeMeasures(models)
			  
			  ###------------Apply filters------------###
			  
			  message("Apply filter settings for ",sample, "...")
			  
			  #obtain QC values for all cells within sample
			  myQC <- getQC(models)
			  
			  #Set all files as selected prior to filtering
			  selected.files <- rownames(myQC)
			  myQC$selected <- rownames(myQC) %in% selected.files 
			  
			  #obtain names of cells that fulfill - total readcount - criteria
			  if (filter.total.read.count == TRUE){
  				selected.files <- rownames(myQC[(myQC$selected == TRUE & myQC$total.read.count >= min.total.read.count),])
  				myQC$selected <- rownames(myQC) %in% selected.files
				}			  
			  
			  #obtain names of cells that fullfil - spikiness - criteria
			  if (filter.spikiness == TRUE){
  			  selected.files <- rownames(myQC[(myQC$selected == TRUE & myQC$spikiness <= max.spikiness),])
  			  myQC$selected <- rownames(myQC) %in% selected.files 
			  }
			  
			  #obtain names of cells that fullfil - bhattacharyya - criteria
			  if (filter.bhattacharyya == TRUE){
			    #convert NA values into zero
			    myQC$bhattacharyya[is.na(myQC$bhattacharyya)] <- c(0)
			    selected.files <- rownames(myQC[(myQC$selected == TRUE & myQC$bhattacharyya >= min.bhattacharyya),])
			    myQC$selected <- rownames(myQC) %in% selected.files 
			  }
			  
			  #obtain names of cells that fullfil - number of segments - criteria
			  if (filter.num.segments == TRUE){
  				selected.files <- rownames(myQC[(myQC$selected == TRUE & myQC$num.segments >= min.num.segments),])
  				myQC$selected <- rownames(myQC) %in% selected.files 
			  }

			  #obtain names of cells that fulfil - minimum weighted copy number - criteria
			  if (filter.min.copy.number == TRUE){
			    for (my.file in my.files){
			      
			      #set all cells to not too low copy number
			      myQC[my.file,'low.copy.number'] <- FALSE
			      
			      total_width <- c()
			      relative_width <- c()
			      relative_copy.number <- c()
			      
			      #determine if a cell has more than zero reads
			      if (models[[my.file]]$qualityInfo$total.read.count != 0){ 
			        
			        #determine if a cell has a too high copy number 
			        total_width <- sum(models[[my.file]]$segments@ranges@width)
			        relative_width <- (models[[my.file]]$segments@ranges@width)/total_width
			        relative_copy.number <- ((models[[my.file]]$segments$copy.number)*relative_width)
			        myQC[my.file,'weighted.copy.number'] <- sum(relative_copy.number)
			        if (myQC[my.file,'weighted.copy.number'] < min.copy.number){
			          myQC[my.file,'low.copy.number'] <- TRUE
			        }
			      }
			    }
			    
			    selected.files <- rownames(myQC[(myQC$selected == TRUE & myQC$low.copy.number == FALSE),])
			    myQC$selected <- rownames(myQC) %in% selected.files   
			    
			  }
			  
			  #obtain names of cells that fulfil - maximum weighted copy number - criteria
			  if (filter.max.copy.number == TRUE){
			    for (my.file in my.files){
			     
			      #set all cells to not too high copy number
			      myQC[my.file,'high.copy.number'] <- FALSE
			      
			      total_width <- c()
			      relative_width <- c()
			      relative_copy.number <- c()
			      
			      #determine if a cell has more than zero reads
			      if (models[[my.file]]$qualityInfo$total.read.count != 0){ 
			        
			        #determine if a cell has a too high copy number 
			        total_width <- sum(models[[my.file]]$segments@ranges@width)
			        relative_width <- (models[[my.file]]$segments@ranges@width)/total_width
			        relative_copy.number <- ((models[[my.file]]$segments$copy.number)*relative_width)
			        myQC[my.file,'weighted.copy.number'] <- sum(relative_copy.number)
			        if (myQC[my.file,'weighted.copy.number'] >= max.copy.number){
			          myQC[my.file,'high.copy.number'] <- TRUE
			        }
			      }
			    }
			    
			    selected.files <- rownames(myQC[(myQC$selected == TRUE & myQC$high.copy.number == FALSE),])
			    myQC$selected <- rownames(myQC) %in% selected.files   
			    
			  }
			  
			  #obtain names of cells that fulfil - diploid - criteria
			  if (filter.diploid == TRUE){
			      for (my.file in my.files){

			      #set all cells to non-diploid
  			    myQC[my.file,'diploid'] <- FALSE
  			    
  			    #determine if a cell has more than zero reads
  			    if (models[[my.file]]$qualityInfo$total.read.count != 0){ 
  			      #determine if a cell is perfectly diploid
  			      
  			      #logical test, are all segment lengths equal to 1? - there are no chromosome breaks)
  			      if (length(unique(models[[my.file]]$segments@seqnames@lengths)) == 1){
  			        #logical test, is the number of unique copy number states equal to 1? - have all chromosomes the same ploidy?
  			        if (length(unique(models[[my.file]]$segments$copy.number[1:22])) == 1){
  			          #logical test, is the single unique copy number state equal to 2? - is the common ploidy equal to 2?
  			          if (unique(models[[my.file]]$segments$copy.number[1:22]) == 2){
  			            #logical test, is the X chromosome copy number state equal to 1 (male) or 2 (female)?
  			            if (models[[my.file]]$segments$copy.number[23] == nXchr){
  			              myQC[my.file,'diploid'] <- TRUE
  			            }
  			          }
  			       }
  			     }
  			   }
			      }
        
			    selected.files <- rownames(myQC[(myQC$selected == TRUE & myQC$diploid == FALSE),])
			    myQC$selected <- rownames(myQC) %in% selected.files   
			    diploid.files <- rownames(myQC[myQC$diploid == TRUE,])
			  }
			  
			  #write QC values to a csv file
			  if (set.measures == TRUE){
			    write.csv(myQC, file = paste0(shortcut, "_QC_measurements.csv"))
			  }
			    
			  #determine which cells were removed
			  removed.files <- rownames(myQC[myQC$selected == FALSE,])
			  
			  #create dataframe with number of selected cells and total cells
			  selected.numbers <- data.frame(type = c("Total", "Selected"), number = c(length(my.files),length(selected.files)))
			  
			  #rapport filtering outcome and corresponding shortcut
			  message(length(selected.files)," out of ",length(my.files), " files were selected after filtering...")
			  selected.shortcut <- paste0(shortcut, "_", length(selected.files),"outof",length(my.files)) 
			  
			  if (length(selected.files) == 0){
					message("None of the cells in ", sample, " fullfills the filter settings..")
					message("Proceed to the next sample..")
			  } else {
					
					
					###--------Actions using selected files, executed per sample-----------###
					
			    #load selected models
			    selected.models <- loadFromFiles(selected.files)
			    
			    #store sampleID for which filtering results in more than zero cells
			    selected.sampleIDs <- append(selected.sampleIDs, sample)
					
					#plot Genomewide heatmap
					if (set.plot.genomewide == TRUE){
  			    message("Plot Genomewide heatmap...")
  			    heatmapGenomewide(selected.models, file = paste0(selected.shortcut,"_selected_genomewide.pdf"))
					}
			    
					#plot single cell profiles for selected and removed files
					if (set.plot.singles == TRUE){
			      message("Plot Single cell profiles...")
					  plotaCGH(selected.files, file = paste0(selected.shortcut,"_selected_singles.pdf"))
					  if (length(removed.files) > 0 ){
					    plotaCGH(removed.files, file = paste0(shortcut,"_removed_singles.pdf"))
					  }
				  }
			   
					#calculate aneuploidy and heterogeneity genomewide and per chromosome, and save as .csv
					selected.karyotype <- karyotypeMeasures(selected.models)
					if (set.measures == TRUE){
					write.csv(selected.karyotype$genomewide, file = paste0(selected.shortcut, "_selected_karyotype_measures_genomewide.csv"))
					write.csv(selected.karyotype$per.chromosome, file = paste0(selected.shortcut, "_selected_karyotype_measures_chromosome.csv"))
					}
					
					#plot aneuploidy and heterogeneity per chromosome and save to QC file directory
          if (set.plot.heterogeneity == TRUE){
  					if (length(selected.models) > 1){
              message("Plot Aneuploidy and Heterogeneity...")
    					myplot <- plotHeterogeneity(selected.models)
    					pdf(paste0(selected.shortcut,"_selected_heterogeneity.pdf")) 
    					print(myplot)
    					dev.off()
  					} else {
  					  message("No Aneuploidy and Heterogeneity Plot, less than 2 cells available...")
  					}
          }
				
					
					
				}
		
			  datacollection[[sample]]$myQC <- myQC
			  datacollection[[sample]]$selected.files <- selected.files
			  datacollection[[sample]]$removed.files  <- removed.files
			  datacollection[[sample]]$selected.numbers <- selected.numbers
			  datacollection[[sample]]$selected.karyotype <- selected.karyotype	
			  datacollection[[sample]]$unfiltered.karyotype <- unfiltered.karyotype	
			  
			  
			  #create directory to store data files after filtering of selected and removed models
			  if (set.move.model.files == TRUE){
  			 
  			  folder.selected.model.files <- c(paste0(project.name, "/", sample ,"/", "QC_FILES_", selected.model,"/selected_model-files"))
  			  dir.create(folder.selected.model.files)
  			  file.copy(selected.files, folder.selected.model.files)
  	
  			  folder.removed.model.files <- c(paste0(project.name, "/", sample ,"/", "QC_FILES_", selected.model,"/removed_model-files"))
  			  dir.create(folder.removed.model.files)
  			  file.copy(removed.files, folder.removed.model.files)
  			  }
			
			
		#determine which samples do not contain any successfull filtered cells
		nonselected.sampleIDs <- c()
		for (sample in sampleIDs){
		  if (!sample %in% selected.sampleIDs){
		    nonselected.sampleIDs <- append(nonselected.sampleIDs, sample)
		  }
		}
		
		#update datacollection with filtered and unfiltered sampleIDS
		datacollection$selected.sampleIDs <- selected.sampleIDs
		datacollection$nonselected.sampleIDs <- nonselected.sampleIDs
		
		#save datacollection to Rdata file for use by other functions
		save(datacollection, file = paste0(summary.shortcut, "datacollection.Rdata"))	
		
		}
		
		
		###------Generate Pdf with Filter Statistics ------###
		###------------------------------------------------###
		
	  ###------Create PDF to compare statistics of selected and removed files-------###
		  
		  message("==> Start Filter Statictics PDF <==")
		  
		  #find maximum values to make uniform plots
		  max.total.read.count <- c(0)
		  max.spikiness <- c(0)
		  max.bhattacharyya <- c(0)
		  max.num.segments <- c(0)
		  for (sample in sampleIDs){
		    max.sample.total.read.count <- max(datacollection[[sample]]$myQC[["total.read.count"]], na.rm = TRUE)
		    max.sample.spikiness <- max(datacollection[[sample]]$myQC[["spikiness"]], na.rm = TRUE)
		    max.sample.bhattacharyya <- max(datacollection[[sample]]$myQC[["bhattacharyya"]], na.rm = TRUE)
		    max.sample.num.segments <- max(datacollection[[sample]]$myQC[["num.segments"]], na.rm = TRUE)
		    if (max.sample.total.read.count > max.total.read.count){
		      max.total.read.count <- max.sample.total.read.count
		    }
		    if (max.sample.spikiness > max.spikiness){
		      max.spikiness <- max.sample.spikiness
		    }
		    if (max.sample.bhattacharyya > max.bhattacharyya){
		      max.bhattacharyya <- max.sample.bhattacharyya
		    }
		    if (max.sample.num.segments > max.num.segments){
		      max.num.segments <- max.sample.num.segments
		    }
		  }
		  
		  #find maximum values to make uniform Karyotype Measurement plots
		  max.aneuploidy <- c(0)
		  max.heterogeneity <- c(0)
		  for (sample in selected.sampleIDs){
		    max.sample.aneuploidy <- max(datacollection[[sample]]$selected.karyotype$genomewide$Aneuploidy)
		    max.sample.heterogeneity <- max(datacollection[[sample]]$selected.karyotype$genomewide$Heterogeneity)
		    if (max.sample.aneuploidy > max.aneuploidy){
		      max.aneuploidy <- max.sample.aneuploidy
		    }
		    if (max.sample.heterogeneity > max.heterogeneity){
		      max.heterogeneity <- max.sample.heterogeneity
		    }
		  }
		  
		  #set axis limits
		  total.read.count.limit <- max.total.read.count*1.05
		  spikiness.limit <- max.spikiness*1.05
		  bhattacharyya.limit <- max.bhattacharyya*1.05
		  num.segments.limit <- max.num.segments+2
		  aneuploidy.limit <- max.aneuploidy*1.05
		  heterogeneity.limit <- max.heterogeneity*1.05
		  
		  
		  ###----Generate plots to compare selected and removed data----###
		  
		  
		  #create sample plot list
		  sample.plots <- list()
		  
		  #Generate plots for samples that contain succesfully filtered cells
		  message("Generate filter plots of selected and removed cells...")
		  for (sample in selected.sampleIDs){
		      p.selected.numbers <- ggplot(datacollection[[sample]]$selected.numbers, aes( x=type, y=number)) +
		      geom_col(show.legend = FALSE) +
		      geom_text(aes(label=number), vjust=1.6, color="white", size = 4) +
		      labs(title=paste0(sample," Number of files selected"), x="Type", y = "Number of samples", size = 10) +
		      theme_classic()
		    
		    p.karyotype <- ggplot(datacollection[[sample]]$selected.karyotype$per.chromosome, aes(x=Aneuploidy, y=Heterogeneity)) +
		      geom_point(size=2) + 
		      xlim(0, aneuploidy.limit) + 
		      ylim(0, heterogeneity.limit) +
		      labs(title=paste0(sample," Karyotype measurements")) +
		      geom_text_repel(aes(label = rownames(datacollection[[sample]]$selected.karyotype$per.chromosome)), size = 3)+
		      theme_classic()
		    
		    p.readcount <- ggplot(datacollection[[sample]]$myQC, aes(x=selected, y=total.read.count, color=selected)) + 
		      geom_boxplot() + 
		      geom_point(shape=16, position = position_dodge2(width = 0.3)) +
		      ylim(0, total.read.count.limit) +
		      labs(title=paste0(sample," Total read count"),x="Selected", y = "Total read count") +
		      theme_classic() +
		      scale_color_discrete_diverging(palette = "Red-Green")
		    
		    
		    p.spikiness <- ggplot(datacollection[[sample]]$myQC, aes(x=selected, y=spikiness, color=selected)) + 
		      geom_boxplot() + 
		      ylim(0, spikiness.limit) +
		      geom_point(shape=16, position = position_dodge2(width = 0.3)) +
		      labs(title=paste0(sample," Spikiness"),x="Selected", y = "Spikiness") +
		      theme_classic() +
		      scale_color_discrete_diverging(palette = "Red-Green")
		    
		    p.bhattacharyya <- ggplot(datacollection[[sample]]$myQC, aes(x=selected, y=bhattacharyya, color=selected)) + 
		      geom_boxplot() + 
		      ylim(0, bhattacharyya.limit) +
		      geom_point(shape=16, position = position_dodge2(width = 0.3)) +
		      labs(title=paste0(sample," Bhattacharyya"),x="Selected", y = "Bhattacharyya") +
		      theme_classic() +
		      scale_color_discrete_diverging(palette = "Red-Green")
		    
		    
		    p.num.segments <- ggplot(datacollection[[sample]]$myQC, aes(x=selected, y=num.segments, color=selected)) + 
		      geom_boxplot() + 
		      ylim(0, num.segments.limit) +
		      geom_point(shape=16, position = position_dodge2(width = 0.3)) +
		      labs(title=paste0(sample," Number of segments"),x="Selected", y = "Number of segments")+
		      theme_classic() +
		      scale_color_discrete_diverging(palette = "Red-Green")
		  
		    #take all plots together in a single grid
		    sample.plots[[sample]] <- plot_grid(p.selected.numbers, p.karyotype, p.readcount, p.spikiness, p.bhattacharyya, p.num.segments,  
		                                        labels = c("A", "B", "C", "D", "E", "F"),
		                                        ncol = 2, nrow = 3, vjust = 1.5)
		    
		  }
		  
		  #Generate plots for samples that do not contain filtered cells
		  if (length(nonselected.sampleIDs)>0){
		    for (sample in nonselected.sampleIDs){
		      #message("Generate plots for samples that do not contain filered cells")
		      p.selected.numbers <- ggplot(datacollection[[sample]]$selected.numbers, aes( x=type, y=number)) +
		        geom_col(show.legend = FALSE) +
		        geom_text(aes(label=number), vjust=1.6, color="white", size = 4) +
		        labs(title=paste0(sample," Number of files selected"), x="Type", y = "Number of samples", size = 10) +
		        theme_classic()
		        
		      p.karyotype <- ggplot(datacollection[[sample]]$unfiltered.karyotype$per.chromosome, aes(x=Aneuploidy, y=Heterogeneity)) +
		        geom_point(size=2) + 
		        labs(title=paste0(sample," Karyotype unfiltered cells")) +
		        geom_text_repel(aes(label = rownames(datacollection[[sample]]$unfiltered.karyotype$per.chromosome)), size = 3) +
		        theme_classic()
		      
		      p.readcount <- ggplot(datacollection[[sample]]$myQC, aes(x=selected, y=total.read.count, color=selected)) + 
		        geom_boxplot(na.rm=TRUE) + 
		        ylim(0, total.read.count.limit) +
		        geom_point(shape=16, position = position_dodge2(width = 0.3), na.rm=TRUE) +
		        labs(title=paste0(sample," Total read count"),x="Selected", y = "Total read count")  +
		        theme_classic()
		      
		      p.spikiness <- ggplot(datacollection[[sample]]$myQC, aes(x=selected, y=spikiness, color=selected)) + 
		        geom_boxplot(na.rm=TRUE) + 
		        ylim(0, spikiness.limit) +
		        geom_point(shape=16, position = position_dodge2(width = 0.3), na.rm=TRUE) +
		        labs(title=paste0(sample," Spikiness"),x="Selected", y = "Spikiness")  +
		        theme_classic()
		      
		      p.bhattacharyya <- ggplot(datacollection[[sample]]$myQC, aes(x=selected, y=bhattacharyya, color=selected)) + 
		        geom_boxplot(na.rm=TRUE) + 
		        ylim(0, bhattacharyya.limit) +
		        geom_point(shape=16, position = position_dodge2(width = 0.3), na.rm=TRUE) +
		        labs(title=paste0(sample," Bhattacharyya"),x="Selected", y = "Bhattacharyya")  +
		        theme_classic()
		      
		      p.num.segments <- ggplot(datacollection[[sample]]$myQC, aes(x=selected, y=num.segments, color=selected)) + 
		        geom_boxplot(na.rm=TRUE) + 
		        ylim(0, num.segments.limit) +
		        geom_point(shape=16, position = position_dodge2(width = 0.3), na.rm=TRUE) +
		        labs(title=paste0(sample," Number of segments"),x="Selected", y = "Number of segments") +
		        theme_classic()
		      
		      #take all plots together in a single grid
		      sample.plots[[sample]] <- plot_grid(p.selected.numbers, p.karyotype, p.readcount, p.spikiness, p.bhattacharyya, p.num.segments,  
		                                          labels = c("A", "B", "C", "D", "E", "F"),
		                                          ncol = 2, nrow = 3, vjust = 1.5)
		    }
		  }

		  #combine plots to compare selected and removed files in a single pdf
		  message("Save Filtering statistics to PDF...")
		  grDevices::pdf(file=paste0(summary.shortcut,"Filtering statistics_overview.pdf"), paper="a4")
		  for (sample in sampleIDs) {
		    print(sample.plots[[sample]])
		  }
		  grDevices::dev.off()
		  
		
		  ###-----------------Generate PCA plot from all selected samples---------------------###
		  
		  message("==> Plot PCA and Karyotype Measurements to PDF <==")
		  
		  if (length(nonselected.sampleIDs)>0){
		    for (sample in nonselected.sampleIDs){
		      message(sample, " cannot be plotted. This sample does not contain any selected cells after filtering...")
		    }
		  }
		  
		  # create empty list and vectors
		  all.selected.models <- list()
		  classes <- c()
		  labels <-c()
		    
		  # collect data and assign classes and labels per sample
		  for(sample in selected.sampleIDs) {
		    all.selected.models[[sample]] <- datacollection[[sample]]$selected.files
		    classes <- append(classes, c(rep(sample, length(datacollection[[sample]]$selected.files))))
		    labels <- append(labels, c(paste(sample,1:length(datacollection[[sample]]$selected.files))))
		  }
		  
		  #use aneufinder package to calculate PCA plot
		  p.PCA <- plot_pca(unlist(all.selected.models), colorBy=classes, PC1=1, PC2=2)
		  
		  #Add aesthetics to PCA plots using ggplot2
		  PCA_plot <- ggplot(p.PCA$data, aes(x=PC1, y=PC2, color = color, label = labels))
		  PCA_plot_jitter <- PCA_plot + 
		    geom_point(position = position_jitter(width = 0.01, height = 0.01, seed = 1L)) +
		    labs(x=p.PCA[["labels"]][["x"]], y = p.PCA[["labels"]][["x"]]) +
		    theme_bw() +
		    scale_color_discrete_qualitative(palette = "dark 3")
		  PCA_plot_jitter_labels <- PCA_plot_jitter + 
		    geom_text_repel(size = 2, position = position_jitter(width = 0.01, height = 0.01, seed = 1L))
		  
		  #Save PCA plot to pdf
		  message("Save PCA plot to PDF...")
		  grDevices::pdf(file=paste0(summary.shortcut,"PCA_plot.pdf"), width = 10, height = 10)
		  
		  #print(p.PCA)
		  print(PCA_plot_jitter)
		  print(PCA_plot_jitter_labels)
		  grDevices::dev.off()
		  
		  
		  
		  ###------Create PDF with genomewide karyotype measurements and for each chromosome-------###
		  ###--------------------------------------------------------------------------------------###
		  
		  message("==> Plot Karyotype Measurements to PDF <==")
		  
		  #combine genomewide karyotype measurements for all samples in a single figure
		  combined.karyotype.genomewide <- data.frame()
		  
		  for (sample in selected.sampleIDs) {
		    my.karyotype.measurements <- datacollection[[sample]]$selected.karyotype$genomewide
		    combined.karyotype.genomewide <- rbind(combined.karyotype.genomewide, my.karyotype.measurements)
		  }
		  rownames(combined.karyotype.genomewide) <- c(selected.sampleIDs)
		  combined.karyotype.genomewide$color <-  rownames(combined.karyotype.genomewide)
		  
		  p.combined.karyotypes.genomewide <- ggplot(combined.karyotype.genomewide, aes(x=Aneuploidy, y=Heterogeneity, color = color)) +
		    geom_point(size=2) + 
		    xlim(0, NA) + 
		    ylim(0, NA) +
		    labs(title="Genomewide Karyotype Measurements") +
		    geom_text_repel(aes(label = rownames(combined.karyotype.genomewide)), size = 3) +
		    theme_bw() +
		    theme(legend.position="none")
		  
		  #calculate max aneuploidy and heterogeneity scores among all samples and chromosomes
		  max.aneuploidy <- c(0)
		  max.heterogeneity <- c(0)
		  for (sample in selected.sampleIDs){
		    max.sample.aneuploidy <- max(datacollection[[sample]]$selected.karyotype$per.chromosome$Aneuploidy)
		    max.sample.heterogeneity <- max(datacollection[[sample]]$selected.karyotype$per.chromosome$Heterogeneity)
		    if (max.sample.aneuploidy > max.aneuploidy){
		      max.aneuploidy <- max.sample.aneuploidy
		    }
		    if (max.sample.heterogeneity > max.heterogeneity){
		      max.heterogeneity <- max.sample.heterogeneity
		    }
		  }
		  
		  #set axis limits at max values plus 5%
		  aneuploidy.limit <- max.aneuploidy*1.05
		  heterogeneity.limit <- max.heterogeneity*1.05
		    
		  #create karyotype measurement plots for each chromosome
		  karyotype.chromosome.plots <- list()
		  my.chromosomes <- rownames(datacollection[[1]]$selected.karyotype$per.chromosome)
		  
		  for (my.chr in my.chromosomes){
		    combined.karyotype.chromosome <- data.frame()
		    for (sample in selected.sampleIDs){
		      sample.karyotype.measurements <- datacollection[[sample]]$selected.karyotype$per.chromosome 
		      combined.karyotype.chromosome <- rbind(combined.karyotype.chromosome, sample.karyotype.measurements[my.chr,])
		    }
		    rownames(combined.karyotype.chromosome) <- c(selected.sampleIDs)
		    combined.karyotype.chromosome$color <- rownames(combined.karyotype.chromosome)
		      
		    karyotype.chromosome.plots[[my.chr]] <- ggplot(combined.karyotype.chromosome, aes(x=Aneuploidy, y=Heterogeneity, color = color)) +
		      geom_point(size=2) + 
		      xlim(0, aneuploidy.limit) + 
		      ylim(0, heterogeneity.limit) +
		      labs(title=paste0("Chr ", my.chr)) +
		      geom_text_repel(aes(label = rownames(combined.karyotype.chromosome)), size = 3) +
		      theme_bw() +
		      theme(legend.position="none")
		    }
		    
		    #Save all karyotype measurement plots to a single pdf
		    message("Save karyotype measurements to PDF...")
		    grDevices::pdf(file=paste0(summary.shortcut,"Karyotype_measurements_overview.pdf"),paper="a4")
	      print(p.combined.karyotypes.genomewide)
		    for (my.chr in my.chromosomes){
		      print(karyotype.chromosome.plots[[my.chr]])  
		    }
		    grDevices::dev.off()
		    
		#end clock
		total.time <- proc.time() - total.time
		message("==> Total time spent: ", round(total.time[3]), "s <==")
	}

}


###---------Function to plot single cell profiles from a list of files-------------###

plotaCGH <- function(cells, file = "default.pdf") {
  
  pdf(file = file, width = 50/2.54, height = 15/2.54)
  for(i in 1:length(cells)) {
    print(plot(cells[i], type = 1))
    #message("Plotted cell ", i, "/", length(cells))
  }
  message("Plotted profiles to file ", file)
  dev.off()
  remove(i)
}


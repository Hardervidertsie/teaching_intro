AnalyzeAll <- function(imagePath = image_paths_dna,
                       adapt.width,
                       adapt.height,
                       nuc.offset, 
                       cyto.offset
                             )
{
  
  
  if( any(
    !is.numeric(
      c(adapt.width,
        adapt.height,
        nuc.offset,
        cyto.offset
         )
      )
    )
  )


  {
    stop("please define numeric arguments")
  }

  
  
  if(any(!grepl(
    ".tif$", imagePath
  )
  )
  )
  {
    stop("imageSet does not contain tif images")
  }
 
 image_nuclei_path <- imagePath
 image_GFP_path <- gsub("Channel2", "Channel1", imagePath)
 

    if(!file.exists("results")){
    dir.create("results")
  }
  nmpje <- gsub("_","" ,str_match( image_nuclei_path, '[A-H]{1}-[0-9]{2}.tif$'))
  

    nuclei.image <- readImage(image_nuclei_path)
    cytosol.image <- readImage(image_GFP_path)

    nuclei.image.f <- medianFilter(nuclei.image, size = 2)
    nmask <- thresh(nuclei.image.f, w=adapt.width, h=adapt.height, offset=nuc.offset)
    nmask <- opening(nmask, makeBrush(7, shape="disc"))
    img <- rgbImage(green= cytosol.image, blue=nuclei.image)
    nmask <- bwlabel(nmask) # label nuclei

    #Cytosol propagation:
    
    cytosolmask = opening(cytosol.image > cyto.offset, makeBrush(5, shape="disc")) # erode + dilate
    celmask <- propagate(cytosol.image, seeds=nmask, mask=cytosolmask)
    colorMode(celmask) <- 0
    colorMode(nmask) <- 0
    res <- paintObjects(celmask, img, col="#ff00ff")  
    res <- paintObjects(nmask, res, col= "#ffff00") 

    nuclei_intensity <- computeFeatures.basic(x = nmask[ , , 1], cytosol.image[ , ,1])[ , 1] # calculate nuclei intensities of channel 1
    cytosol_intensity <- computeFeatures.basic(x = celmask[ , , 1], cytosol.image[ , ,1])[ , 1 ] # calculate cell intensities of channel 1
    
    cytosol_intensity <- cytosol_intensity - nuclei_intensity # subtract nuclei intensity
    
    data_out <- list(nuclei_intensity = nuclei_intensity, cytosol_intensity = cytosol_intensity,
              ratio_cyto_nucl = cytosol_intensity / nuclei_intensity)
    
    
    data_out$cellN <- max(nmask)
    data_out$image.name <- as.character(nmpje)
    
    writeImage(res, paste("results/", nmpje, ".jpg", sep =""), quality = 85)
  
  print(paste( nmpje, " processed"))
  return(data_out)
} # end AnalyzeAllImages function
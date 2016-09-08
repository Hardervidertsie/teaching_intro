AnalyzeAll <- function(adapt.width,
                             adapt.height,
                             nuc.offset, 
                             cyto.offset,
                             imagePath,
                             is.cytosol)
{
  
  if( !is.numeric(c
                  (adapt.width,
                   adapt.height,
                   nuc.offset, 
                   cyto.offset
                  )
  )
  )
  
  
  {
    stop("please define numeric arguments")
  }
  
  
  imageSet <-dir(imagePath) 
  
  if(any(!grepl(
    ".tif$", imageSet
  )
  )
  )
  {
    stop("imageSet does not contain tif images")
  }
  
 image.list.nuclei <- imageSet[grep("c1.tif", imageSet)] 
 image.list.GFP <- imageSet[grep("c2.tif", imageSet)] 
  
  
  mean.int.data = list()
  if(!file.exists("output")){
    dir.create("output")
  }
  nmpje <- gsub("_","" ,str_match( imageSet[1], '[A-Z 1-9]{4,10}_'))
  dir.create(
    paste("output",
    nmpje
  , sep ="/")
  )
  
  for( i in seq_along(image.list.nuclei))
  {
  
    nuclei.image <- suppressWarnings(readImage(paste(imagePath, '/', image.list.nuclei[i], sep ="")))
    cytosol.image <- suppressWarnings(readImage(paste(imagePath,'/', image.list.GFP[i], sep ="")))
    
    nuclei.image.f <- medianFilter(nuclei.image, size = 2)
    nmask = thresh(nuclei.image.f, w=adapt.width, h=adapt.height, offset=nuc.offset)
    nmask = opening(nmask, makeBrush(7, shape="disc"))
    img <- rgbImage(green= cytosol.image, blue=nuclei.image)
    nmask = bwlabel(nmask) # label nuclei
    
    #Cytosol propagation:
    if(is.cytosol){
    cytosolmask = opening(cytosol.image > cyto.offset, makeBrush(5, shape="disc")) # erode + dilate
    celmask = propagate(cytosol.image, seeds=nmask, mask=cytosolmask)
    
    res <- paintObjects(celmask, 16*img, col="#ff00ff")  # increase intensity of original image to see better (10X in this case). Add cell mask outline.
    res <- paintObjects(nmask, res, col= "#ffff00") # add nuclei mask outline
    
    mean.int.data[[i]] <- as.data.frame(computeFeatures.basic(x = celmask, cytosol.image))
    image.name <- gsub("c1.tif", "", str_match( image.list.nuclei[i], "xy[0-9]{2}c1.tif$"))
    mean.int.data[[i]]$cellN <- max(nmask)
    mean.int.data[[i]]$image.name <- as.character(image.name)
    mean.int.data[[i]]$cell <- gsub("images/", "", imagePath)
    writeImage(res, paste("output/", nmpje, "/", image.name, ".jpg", sep =""), quality = 85)
    } else { # measure gfp in nuclei
        res <- paintObjects(nmask, 16*img, col= "#ffff00") # add nuclei mask outline
        mean.int.data[[i]] <- as.data.frame(computeFeatures.basic(x = nmask, cytosol.image))
        image.name <- gsub("c1.tif", "", str_match( image.list.nuclei[i], "xy[0-9]{2}c1.tif$"))
        mean.int.data[[i]]$cellN <- max(nmask)
        mean.int.data[[i]]$image.name <- as.character(image.name)
        mean.int.data[[i]]$cell <- gsub("images/", "", imagePath)
        writeImage(res, paste("output/", nmpje, "/", image.name, ".jpg", sep =""), quality = 85)
      
      }
    }
  
  all.mean.int.data <- do.call("rbind", mean.int.data)
  
  return(all.mean.int.data)
} # end AnalyzeAllImages function
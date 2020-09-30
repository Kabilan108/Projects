library(minfi)

## Needs to be uncompressed
manifestFile <- '../../data/manifest_mouse.csv'
stopifnot(file.exists(manifestFile))
maniTmp <- minfi:::read.manifest.Mammal(manifestFile)

## Checking
manifest <- maniTmp$manifest
address.all <- c(manifest$U, manifest$M)
sum(address.all == "")
sum(is.na(address.all))
address.all <- address.all[address.all != ""]
length(address.all)
stopifnot(!anyDuplicated(address.all))
library(illuminaio)
epic <- readIDAT("../../data/MouseArrayManifest.csv")
address.epic <- as.character(epic$MidBlock)
sum(!address.epic %in% address.all) ## Set of addresses in the IDAT file not part of the manifest
sum(! address.all %in% address.epic) ## set of addresses not in IDAT file.
nrow(manifest)
any(manifest$U != "" & !manifest$U %in% address.epic)
any(manifest$M != "" & !manifest$M %in% address.epic)
wh <- which(manifest$M != "" & !manifest$M %in% address.epic)
tmp <- manifest[wh,]
table(tmp$Infinium_Design_Type)
table(tmp$Color_Channel)
table(tmp$Methyl450_Loci)

## Controls ok
all(maniTmp$controls[,1] %in% address.epic)

## Manifest package
maniList <- maniTmp$manifestList
## Manually removing 1031 CpGs with missing addresses
## dropCpGs <- manifest$Name[manifest$AddressB != "" & !manifest$AddressB %in% address.epic]
## table(substr(dropCpGs, 1,2))
## maniList$TypeI <- maniList$TypeI[! maniList$TypeI$Name %in% dropCpGs,]

IlluminaHumanMethylationEPICmanifest <- IlluminaMethylationManifest(TypeI = maniList$TypeI,
                                                                    TypeII = maniList$TypeII,
                                                                    TypeControl = maniList$TypeControl,
                                                                    TypeSnpI = maniList$TypeSnpI,
                                                                    TypeSnpII = maniList$TypeSnpII,
                                                                    annotation = "Mouse Manifest")
stopifnot(validObject(IlluminaHumanMethylationEPICmanifest))
save(IlluminaHumanMethylationEPICmanifest, compress = "xz",
     file = "../../data/manifest.rda")
sessionInfo()
rm(list = ls())

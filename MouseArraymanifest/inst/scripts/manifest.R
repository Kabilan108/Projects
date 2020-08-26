library(minfi)

## Needs to be uncompressed
manifestFile <- "C:/Users/kabil/Desktop/manifest_epic_v10b2.zip"
stopifnot(file.exists(manifestFile))
maniTmp <- minfi:::read.manifest.EPIC(manifestFile)

## Checking
manifest <- maniTmp$manifest
address.all <- c(manifest$AddressA, manifest$AddressB)
sum(address.all == "")
sum(is.na(address.all))
address.all <- address.all[address.all != ""]
length(address.all)
stopifnot(!anyDuplicated(address.all))
library(illuminaio)
epic <- readIDAT("../../../data_files/Demo_Data_EPIC/200144450018/200144450018_R04C01_Grn.idat")
address.epic <- as.character(epic$MidBlock)
sum(!address.epic %in% address.all) ## Set of addresses in the IDAT file not part of the manifest
sum(! address.all %in% address.epic) ## set of addresses not in IDAT file.
nrow(manifest)
any(manifest$AddressA != "" & !manifest$AddressA %in% address.epic)
any(manifest$AddressB != "" & !manifest$AddressB %in% address.epic)
wh <- which(manifest$AddressB != "" & !manifest$AddressB %in% address.epic)
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
                                                                    annotation = "IlluminaHumanMethylationEPIC")
stopifnot(validObject(IlluminaHumanMethylationEPICmanifest))
save(IlluminaHumanMethylationEPICmanifest, compress = "xz",
     file = "../../data/IlluminaHumanMethylationEPICmanifest.rda")
sessionInfo()
rm(list = ls())

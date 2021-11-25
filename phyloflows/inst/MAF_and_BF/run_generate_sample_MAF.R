# Here I prepare the input to be 'sent' to generate_sample_MAF.R
require(data.table)

outdir <- '~/Documents/2021/MAFs'
outfilename <- paste0(gsub('-', '', Sys.Date()),'_MAF.csv')
outfilename <- file.path(outdir, outfilename)

dir.BF.locations <- '~/Documents/ratmann_pangea_deepsequencedata/'

BF.locations <- as.data.table(readRDS(file.path('BF_file_locations_updated.RDS') ))
BF.locations <- unique(BF.locations[, .(FILE, PANGEA_ID)])
BF.locations <- BF.locations[file.exists(FILE)]

if(0) # TODO: some PANGEA_IDs have multiple reads assigned to them. Are they the same? Seems like not
{
  idx <- BF.locations[,.N,by='PANGEA_ID'][N>1, PANGEA_ID]
  
  for(id in idx)
  {
    tmp <- BF.locations[PANGEA_ID == id]
    tmp1 <- tmp[, {z <- read.csv(FILE)} ,by=FILE]
    
  }
}else{
  BF.locations <- BF.locations[,list(FILE=FILE[1]),by=PANGEA_ID]
}
setcolorder(BF.locations, c('FILE', 'PANGEA_ID'))


path_file <- file.path(outdir, paste0(gsub('-','',Sys.Date()),'_bfpath.csv'))

# write.table instead of write.csv due to colnames still showing up.
write.table(BF.locations,  path_file, row.names = FALSE, col.names = FALSE, sep=',')
cmd <- paste0("Rscript generate_sample_MAF.R '", path_file, "' '", outfilename, "'")
system(cmd)

read.MAF.file <- function(file.path)
{
  df <- read.csv(file.path, header = F)
  colnames(df) <- 1:ncol(df)
  rownames(df) <- df[,1]; df <- df[,-1]
  colnames(df) <- paste0('pos',df[1,])
  df <- df[-1,]
  return(df)
}

MAF <- read.MAF.file(outfilename)
testMAF <- read.MAF.file('~/Documents/HIV-phyloTSI-main/ExampleInputs/testmaf.csv')

colSums(MAF)
colSums(testMAF)

colnames(testMAF) <- 1:ncol(testMAF)
rownames(testMAF) <- testMAF[,1]; testMAF <- MAF[,-1]





##############################################################
# I think I can try to run Tanya's program immediately after.
##############################################################
PATSTATS <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/patstats_02_05_30_min_read_100_max_read.csv'

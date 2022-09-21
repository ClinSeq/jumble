# Markus Mayrhofer 2022
# Dependencies ------------------------------------------------------------
{
    suppressPackageStartupMessages(library(optparse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
    #suppressPackageStartupMessages(library(csaw))
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    suppressPackageStartupMessages(library(BSgenome))
    suppressPackageStartupMessages(library(Repitools))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(VariantAnnotation))
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(MASS))
    suppressPackageStartupMessages(library(rtracklayer))
    suppressPackageStartupMessages(library(rospca))
}
# Options ------------------------------------------------------------

# Define options
option_list <- list(
    make_option(c("-i", "--input_folder"), action = "store", type = "character",default = '.', 
                help = "folder of count RDS files"),
    make_option(c("-o", "--output_folder"), action = "store", type = "character",default = '.', 
                help = "folder to put reference file"),
    make_option(c("-c", "--cores"), action = "store", type = "character",default = '4', 
                help = "cores to use")
)
opt <- parse_args(OptionParser(option_list = option_list))


if (!is.null(opt$cores)) registerDoParallel(cores=as.integer(opt$cores))

# opt <- list(input_folder='.',output_folder='.',cores=4)


# Read count files --------------------------------------------------------

files=dir(path = opt$input_folder,pattern = 'counts.RDS$',full.names = T)

# read files
{
    allcounts <- NULL
    ntargets <- 0
    bed_files <- ''
    for (i in 1:length(files)) {
        counts <- readRDS(files[i])
        # counts$ranges <- counts$ranges[-removals$V1]
        # counts$count <- counts$count[-removals$V1]
        # counts$count_short <- counts$count_short[-removals$V1]
        # saveRDS(counts,file = files[i])
        if (is.null(counts$input_bam_file))
            counts$input_bam_file <- files[i]
        ntargets[i] <- length(counts$count)
        if (!is.null(counts$target_bed_file)) 
            bed_files[i] <- str_remove(string = counts$target_bed_file,pattern = '^.*/') else bed_files[i] <- ''
        allcounts[[i]] <- counts
    }
    
    if (length(table(ntargets))>1) stop('Number of bins differs between samples.')
    if (length(table(bed_files))>1) stop('BED file differs between samples.')
}

wgs <- F
if (is.null(allcounts[[1]]$target_bed_file)) wgs <- T

# Make tables ------------------------------------------------------------

counts <- allcounts[[1]]

targets <- as.data.table(counts$ranges)
#targets[queryHits(overlaps),background:=subjectHits(overlaps)]
#background <- as.data.table(counts$background_ranges)


targets <- targets[,.(sample='',
                      bin=1:.N,
                      is_target=T,
                      chromosome=as.character(seqnames),
                      start,end,
                      mid=round((end+start)/2),
                      width)]
targets[width!=min(width),is_target:=F]
# background <- background[,.(sample='',
#                             background=1:.N,
#                             chromosome=as.character(seqnames),
#                             start,end,
#                             mid=round((end+start)/2),
#                             width)]




# Gene annotation ------------------------------------------------------------

mygenes <- c("ASXL1", "JAK2", "FANCA", "SDHB", "SDHC", "BMPR1A", "SDHD", 
  "POT1", "FANCM", "EPCAM", "GREM1", "PRSS1", "RECQL", "PHOX2B", 
  "FH", "EP300", "PAX5", "PDCD1", "CD274", "JAK1", "MLH3", "NKX3-1", 
  "SPEN", "CDC27", "RAD52", "CHD1", "EZH2", "MED12", "PIK3CG", 
  "PMS1", "SOX9", "HLA-C", "CDKN2C", "AR", "MSH6", "FGFR2", "MSH2", 
  "PALB2", "MTOR", "ATR", "BAP1", "CHEK2", "CDKN1B", "CTNNB1", 
  "DICER1", "ERCC2", "FBXW7", "GNAS", "VHL", "KEAP1", "ARID1A", 
  "KDM6A", "KMT2D", "MAP3K1", "MEN1", "NBN", "PMS2", "PRKAR1A", 
  "SMARCA4", "ARAF", "ATRX", "BARD1", "KMT2C", "RNF43", "FOXA1", 
  "MAP2K4", "PIK3R1", "AXIN2", "BRIP1", "CDK12", "MUTYH", "NCOR1", 
  "NF2", "SMAD2", "CTCF", "CYLD", "GATA3", "NOTCH1", "TGFBR2", 
  "ZFHX3", "B2M", "FAT1", "TBX3", "CREBBP", "MITF", "PBRM1", "SETD2", 
  "WT1", "IL6ST", "ARID2", "EPHA3", "RAD51D", "CIC", "KDM5C", "MGA", 
  "MRE11", "MSH3", "NOTCH2", "PPM1D", "SDHA", "CASP8", "KMT2A", 
  "RAD50", "RAD51C", "STAG2", "POLQ", "CHD4", "LZTR1", "AXIN1", 
  "FUBP1", "HLA-A", "KMT2B", "PLCG1", "PSIP1", "TCF12", "TLR4", 
  "USP9X", "ERCC4", "FANCC", "NTRK1", "ASXL2", "LATS1", "DDX3X", 
  "ZMYM2", "BLM", "TAF1", "GNA13", "GPS2", "HLA-B", "RASA1", "SMC3", 
  "TRAF3", "ACVR1B", "GABRA6", "IRF2", "KEL", "CHEK1", "RAD51B", 
  "NTRK3", "ARID5B", "CSDE1", "CYSLTR2", "ELF3", "H3F3C", "LATS2", 
  "ARHGAP35", "ALB", "HIST1H1E", "SMC1A", "SPTA1", "CD70", "BTG2", 
  "ATM", "BRCA1", "BRCA2", "PTEN", "NF1", "PDGFRA", "APC", "MET", 
  "PIK3CA", "RET", "PTCH1", "ALK", "CDH1", "ERBB2", "RB1", "SMAD4", 
  "TP53", "TSC1", "TSC2", "ESR1", "EGFR", "MLH1", "SMARCB1", "MAP2K1", 
  "FGFR1", "RAD51", "CDKN2A", "KIT", "FGFR3", "STK11", "AMER1", 
  "POLD1", "RBM10", "TCF7L2", "RAD21", "RECQL4", "CHD3", "SMARCA1", 
  "COL5A1", "ZMYM3", "CDKN1A", "HGF", "RPL5", "TBL1XR1", "PLCB4", 
  "ERCC3", "GNA11", "GNAQ", "U2AF1", "IDH2", "ERBB3", "SPOP", "ERBB4", 
  "PIK3CB", "MAX", "MYC", "CCND1", "STAT3", "ARID1B", "NFE2L2", 
  "PPP2R1A", "RAF1", "BCOR", "JUN", "KDR", "MAPK1", "MYCN", "RAC1", 
  "FGFR4", "SOS1", "ACVR1", "PPP6C", "RHOA", "CBL", "ERG", "HIST1H3B", 
  "SMAD3", "EPAS1", "GATA2", "SMO", "AKT2", "PREX2", "RHEB", "RIT1", 
  "BCL6", "CDK6", "DAXX", "FOXL2", "HIST1H3C", "SOCS1", "RRAS2", 
  "RXRA", "SOX17", "BCL2L11", "EIF1AX", "HIST1H1C", "CCND3", "IRF4", 
  "MAP2K2", "PTPRS", "RICTOR", "TP63", "AKT1", "BRAF", "HRAS", 
  "KRAS", "IDH1", "NRAS", "AKT3", "BCL10", "FOXO1", "MAP3K13", 
  "PTPRT", "SUZ12", "AXL", "IGF1R", "KLF4", "SDHAF2", "SMARCD1", 
  "TRAF7", "DIS3", "DNAJB1", "EIF4A2", "IRS2", "MYOD1", "PARP1", 
  "RPTOR", "ERRFI1", "GLI1", "GSK3B", "LYN", "PIK3CD", "PRKCI", 
  "CTLA4", "PAK7", "ANKRD11", "DNMT1", "DNMT3B", "FAM58A", "HIST1H3H", 
  "KNSTRN", "MDC1", "MST1", "NCOA3", "PDPK1", "PLK2", "PPP4R2", 
  "PRDM14", "RPS6KA4", "SESN2", "SPRED1", "STK19", "TGFBR1", "POLE", 
  "CUL3", "PIK3R2", "FOXP1", "H3F3A", "PIM1", "PTPRD", "EPHB1", 
  "GTF2I", "BRD4", "CRLF2", "ETV6", "TCF3", "TET1", "TNFRSF14", 
  "EPHA7", "MEF2B", "RPL22", "PGR", "INPPL1", "NUP93", "ACVR2A", 
  "TCEB1", "ELOC", "HNF1A", "CDK4", "SF3B1", "XPO1", "CBFB", "PTPN11", 
  "CARD11", "RUNX1", "NSD1", "BCL2", "CD79B", "CEBPA", "IKZF1", 
  "MYD88", "RARA", "TET2", "WHSC1", "NSD2", "SRSF2", "FLT3", "DNMT3A", 
  "EWSR1", "GPR126", "PDGFRB", "ROS1", "TERT", "TMPRSS2", "COMT", 
  "CYP2D6", "DPYD", "HOXB13", "HSD3B1", "NQO1", "TPMT", "UGT1A1", 
  "APOBEC3B", "AR_ENHANCER", "CCNE1", "CDKN2B", "CSMD1", "DROSHA", 
  "FANCD2", "FLT4", "GPC3", "INTS4", "LRP1B", "MCL1", "MDM2", "MDM4", 
  "MYCL", "NKX2", "NKX3", "PARK2", "PDE4D", "SKP2", "SOX2", "SUFU", 
  "WWOX", "ZBTB16", "MRE11A", "RAD54L", "ETV1")


# 
# d <- detailRanges(ucsc_ranges, orgdb=org.Hs.eg.db,
#                   txdb=TxDb.Hsapiens.UCSC.hg19.knownGene)
# 
# # targets can overlap multiple genes:
# t <- strsplit(d$overlap,',')
# t <- lapply(t,function(t) str_remove(t,pattern = ':.*'))
# t <- lapply(t,function(t) paste0('<',t,'>'))
# t <- lapply(t,function(t) paste0(t,collapse = ','))
# 
# targets[,gene:=''][targets$is_target==T,gene:=unlist(t)][gene=='<>',gene:='']
# targets[targets$is_target==F,gene:='Background']

mart <- fread('~/Analysis/2022-09_check_genes/mart_export.txt')
mart <- mart[`Gene type`=='protein_coding']
mart <- mart[`Chromosome/scaffold name` %in% c(1:22,'X','Y')]

# If targeted, will use only the pre-listed genes:
if (!wgs) mart <- mart[`Gene name` %in% mygenes | `Gene Synonym` %in% mygenes]

mart <- mart[,.(ensembl_id=`Gene stable ID`,chromosome=`Chromosome/scaffold name`,
                start=`Gene start (bp)`,end=`Gene end (bp)`,
                symbol=`Gene name`)]#,synonym=`Gene Synonym`
mart <- unique(mart)

# Add T_E_intergenic and AR_enhancer to gene table
#mart[symbol %in% c('TMPRSS2','ERG')]
te_int <- data.table(ensembl_id=NA,chromosome=21,start=40033704,end=42836478,symbol='T_E_intergenic')
mart <- rbind(mart,te_int)

# Add AR enhancer to gene table
#ggplot(targets[chromosome=='X' & is_target & start>60e6 & end<70e6]) + geom_point(aes(x=start,y=0,col=gene))
#targets[chromosome=='X' & is_target & start>65.5e6 & end<66.5e6]
ar_enh <- data.table(ensembl_id=NA,chromosome='X',start=66100404,end=66160987,symbol='AR_enhancer')
mart <- rbind(mart,ar_enh)

generanges <- makeGRangesFromDataFrame(mart)
binranges <- counts$ranges

# Compute overlap table
gene_overlap <- as.data.table(findOverlaps(binranges,generanges))
# Remove background matches
gene_overlap <- gene_overlap[queryHits %in% which(targets$is_target)][,bin:=queryHits][,gene:=subjectHits]
# Add symbol to overlap table
gene_overlap[,symbol:=mart[subjectHits]$symbol]

targets[,gene:='']
for (i in unique(gene_overlap$queryHits)) {
    targets[i,gene:=paste0(gene_overlap[queryHits==i]$symbol,collapse = ',')]
}


# # label some genes
# label_genes <- c('AR','ATM','BRCA2','PTEN','RB1','ERG','CDK12','TMPRSS2')
# targets[,label:=as.character(NA)]
# for (g in label_genes) targets[str_detect(gene,paste0('<',g,'>')),label:=g]
# 
# # remove the <>
# targets[,gene:=str_remove_all(gene,'[<>]')]

# GC content ------------------------------------------------------------

ucsc_ranges <- counts$ranges[targets$is_target==T]
seqlevelsStyle(ucsc_ranges) <- "UCSC"

targets[,gc:=as.double(NA)]
targets[is_target==T]$gc <- gcContentCalc(ucsc_ranges , organism=Hsapiens)

# Mappability ------------------------------------------------------------

targets[,map:=as.double(NA)]
targets[is_target==T]$map <- mappabilityCalc(ucsc_ranges , organism=Hsapiens)


# Backbone definition ------------------------------------------------------------

set.seed(25) # <------------------ To be reproducible.
max_backbone_in_gene <- 20 #  <--- applies to targeted


if (wgs) targets[,is_backbone:=chromosome %in% 1:22]

if (!wgs) {
    targets[,is_backbone:=chromosome %in% 1:22 & gene=='']
    for (g in unique(targets$gene)) {
        ix=targets[gene==g & chromosome %in% as.character(1:22),.I]
        n <- length(ix)
        if (n>max_backbone_in_gene) ix=ix[order(rnorm(n))][1:max_backbone_in_gene]
        targets[ix,is_backbone:=T]
    }
}
#background[,is_backbone:=chromosome %in% 1:22]


# Templates ready ------------------------------------------------------------
target_template <- copy(targets)
#background_template <- copy(background)


# Default SNP allele bias
# Cannot be used where there are no allele
# # if no VCF files are available, use these values for SNP allele bias.
# default_bias <- data.table(
#     type=c("C>T", "G>A", "T>A", "A>G", "G>T", "other",
#            "A>T", "C>G", "T>C", "A>C", "T>G", "C>A", "G>C"),
#     bias=c(0.057, 0.08,0.04, 0.006, 0.058, 0.078, 0.053, 
#            0.048, 0.008, 0.021, 0.001,0.073, 0.04))


# Read VCF files ------------------------------------------------------------

peakx <- function(data) {
    d <- density(data)
    maxy <- max(d$y)
    maxx <- d$x[d$y==maxy][1]
    return(maxx)
} 

vcf_files=dir(path = opt$input_folder,pattern = '.vcf',full.names = T)

if (length(vcf_files) > 0) {
    
    if (length(vcf_files) != length(files)) 
        warning('Not one VCF file per sample?')
    
    snp_tables <- foreach(i=1:length(vcf_files)) %dopar% try({
        vcf <- readVcf(vcf_files[i])
        
        # remove SNPs that do not have 2 alleles
        alleles <- as.data.table(table(as.data.table(alt(vcf))$group))$N
        vcf <- vcf[alleles==1 & sapply(geno(vcf)$AD,length)==2]
        
        this_list <- NULL
        for (j in 1:ncol(vcf)) {
            name <- colnames(vcf)[j]
            
            g <- geno(vcf[,j])
            
            # Prepare table of SNPs for this sample
            snp_table <- as.data.table(rowRanges(vcf))[,1:3][,chromosome:=seqnames][,c(4,2,3)]
            snp_table <- cbind(
                data.table(
                    sample=name,
                    id=names(vcf)),
                snp_table
            )
            
            # Get alleles
            snp_table$ref_allele <- as.character(ref(vcf))
            alt_allele <- as.data.table(alt(vcf))[, .(values = list(value)), by = group][,values]
            snp_table$n_alt_alleles <- sapply(alt_allele, length)
            snp_table$alt_allele <- sapply(alt_allele, "[[", 1)
            
            # Get SNP type
            # from <- rep('C/G',length(vcf)); from[snp_table$ref_allele %in% c('A','T')] <- 'A/T'
            # to <- rep('C/G',length(vcf)); to[snp_table$alt_allele %in% c('A','T')] <- 'A/T'
            snp_table[!str_detect(ref_allele,'^[ACGT]$'),ref_allele:='other']
            snp_table[!str_detect(alt_allele,'^[ACGT]$'),alt_allele:='other']
            snp_table[,type:=paste0(ref_allele,'>',alt_allele)]
            snp_table[str_length(ref_allele)!=1 | str_length(alt_allele)!=1]$type <- 'other'
            
            # Get read counts
            snp_table$AD <- sapply(g$AD, "[[", 2)
            snp_table$RD <- sapply(g$AD, "[[", 1)
            snp_table$DP <- snp_table$AD+snp_table$RD
            snp_table[,logDP:=log2(DP)]
            
            # Compute raw allele ratio
            raw_allele_ratio <- unname(round(snp_table$AD/snp_table$DP,4))
            raw_allele_ratio[is.nan(raw_allele_ratio)] <- 0
            snp_table$allele_ratio <- raw_allele_ratio
            snp_table[allele_ratio==0,type:='none']
            
            # Add to list of tables
            #fwrite(snp_table,paste0(name,'.snptable.csv'))
            this_list[[j]] <- snp_table
        }
        this_set <- rbindlist(this_list)
        return(this_set)
    }, silent=T)
    
    snp_table <- unique(rbindlist(snp_tables)) # in case one sample appears more than once
    
    snp_table[,snp:=paste(id,ref_allele,alt_allele)]
    
    # map snps to bins
    unique_table <- unique(snp_table[,.(id,chromosome,start,end)])
    overlap <- findOverlaps(makeGRangesFromDataFrame(unique_table),counts$ranges)
    unique_table[queryHits(overlap),bin:=subjectHits(overlap)]
    snp_table <- merge(snp_table,unique_table,by=c('id','chromosome','start','end'))
    
    # Drop SNPs not mapping to bins
    snp_table <- snp_table[bin %in% target_template[is_target==T]$bin]
    
    setkey(snp_table,'sample')
    rm(snp_tables)
    #saveRDS(snp_table,'snp_table.RDS')
    
    
    
    
    # SNP GC ------------------------------------------------------------
    
    
    gc_ranges <- makeGRangesFromDataFrame(unique(snp_table[,.(chromosome,start,end)]))
    seqlevelsStyle(gc_ranges) <- "UCSC"
    gc_ranges <- gc_ranges[width(gc_ranges)==1]
    start(gc_ranges) <- start(gc_ranges)-1
    end(gc_ranges) <- end(gc_ranges)+1
    gc_ranges$gc3 <- gcContentCalc(gc_ranges, organism=Hsapiens)
    start(gc_ranges) <- start(gc_ranges)-1
    end(gc_ranges) <- end(gc_ranges)+1
    gc_ranges$gc5 <- gcContentCalc(gc_ranges , organism=Hsapiens)
    start(gc_ranges) <- start(gc_ranges)-48
    end(gc_ranges) <- end(gc_ranges)+48
    gc_ranges$gc101 <- gcContentCalc(gc_ranges , organism=Hsapiens)
    start(gc_ranges) <- start(gc_ranges)+50
    end(gc_ranges) <- end(gc_ranges)-50
    seqlevelsStyle(gc_ranges) <- "NCBI"
    gc_table <- as.data.table(gc_ranges)[,-c('strand','width')][,chromosome:=seqnames][,-'seqnames']
    snp_table <- merge(snp_table,gc_table,all=T,by=c('chromosome','start','end'))
    #saveRDS(snp_table,'snp_table.RDS')
    
    
    # SNP training table ------------------------------------------------------------
    
    train_table <- snp_table[allele_ratio>.25]
    
    
    # Check numbers per type
    n_per_type <- as.data.table(table(train_table$type))[order(N,decreasing = T)]
    
    
    # Compute features
    train_table[,hom_alt:=F][chromosome %in% 1:22 & allele_ratio > .95,hom_alt:=T]
    snp_table[,hom_alt:=F][chromosome %in% 1:22 & allele_ratio > .95,hom_alt:=T]
    train_table[,het:=F][allele_ratio > .35 & allele_ratio < .6,het:=T]
    
    # Depth peaks by sample of log2 (DP of all het C>T and G>A, nonX)
    t <- c('C>T','G>A')
    train_table[,peak_ld:=peakx(log2(DP[type %in% t & hom_alt==T & chromosome %in% 1:22])),by=c('sample')]
    snp_table[,peak_ld:=peakx(log2(DP[type %in% t & hom_alt==T & chromosome %in% 1:22])),by=c('sample')]
    
    
    # SNP type specific homozygous peak depth by sample, as log2 (this)-log2(sample_median)
    c('C>T','G>A','A>G','T>C')
    train_table[,homC_T:=peakx(log2(DP[hom_alt==T & type=='C>T']))-peak_ld,by=c('sample')]
    train_table[,homT_C:=peakx(log2(DP[hom_alt==T & type=='T>C']))-peak_ld,by=c('sample')]
    train_table[,homA_G:=peakx(log2(DP[hom_alt==T & type=='A>G']))-peak_ld,by=c('sample')]
    train_table[,homG_A:=peakx(log2(DP[hom_alt==T & type=='G>A']))-peak_ld,by=c('sample')]
    
    snp_table[,homC_T:=peakx(log2(DP[hom_alt==T & type=='C>T']))-peak_ld,by=c('sample')]
    snp_table[,homT_C:=peakx(log2(DP[hom_alt==T & type=='T>C']))-peak_ld,by=c('sample')]
    snp_table[,homA_G:=peakx(log2(DP[hom_alt==T & type=='A>G']))-peak_ld,by=c('sample')]
    snp_table[,homG_A:=peakx(log2(DP[hom_alt==T & type=='G>A']))-peak_ld,by=c('sample')]
    
    # Look at how they correlate
    temp <- unique(train_table[,.(peak_ld,homC_T,homT_C,homA_G,homG_A)])
    #pairs(temp)
    
    # Heterozugous SNPs table
    het_table <- train_table[het==T]
    #plot(density(het_table$allele_ratio))
    
    
    # Individual reference bias (het, near-normal) as log2
    het_table[,bias:=round(log2(RD/AD),4)]
    #plot(density(het_table$bias))
    
    # Allele bias model, general ------------------------------------------------------------
    
    
    # build test model of bias from selected features, using RLM
    biasmodel <- rlm(bias ~ logDP+gc3+gc5+gc101+peak_ld+homC_T+homT_C+homA_G+homG_A, data=het_table)
    #summary(biasmodel)
    
    # same using lm
    lm_model <- lm(bias ~ logDP+gc3+gc5+gc101+peak_ld+homC_T+homT_C+homA_G+homG_A, data=het_table)
    #summary(lm_model)
    
    # But we really need separate models per variant type
    snp_rlm_model <- NULL
    for (ref in unique(het_table$ref_allele)) 
        for (alt in unique(het_table$alt_allele)) 
            if(ref!=alt & ref!='other') {
                data <- het_table[ref_allele==ref & alt_allele==alt]
                if (alt=='other') data <- het_table[ref_allele==ref]
                #if (ref=='other') data <- het_table[alt_allele==alt]
                biasmodel <- rlm(bias ~ logDP+gc3+gc5+gc101+peak_ld+homC_T+homT_C+homA_G+homG_A, data=data)
                cat(ref,'>',alt,':\n')
                print(summary(biasmodel))
                snp_rlm_model[[paste0(ref,'_',alt)]] <- biasmodel
            }
    
    # compute bias estimates
    snp_table[,bias_estimate:=peakx(het_table$bias)] # default
    het_table[,bias_estimate:=peakx(het_table$bias)] # default
    for (ref in unique(snp_table$ref_allele)) 
        for (alt in unique(snp_table$alt_allele)) 
            if(ref!=alt & ref != 'other') {
                cat(ref,'>',alt,'\n')
                thismodel <- snp_rlm_model[[paste0(ref,'_',alt)]]
                
                # in het snp table
                ix <- snp_table$ref_allele==ref & snp_table$alt_allele==alt
                cat(sum(ix),'total \n')
                predicted_bias <- predict(thismodel,snp_table[ix])
                cat(length(predicted_bias),'predicted \n')
                hist(predicted_bias)
                snp_table[ix,bias_estimate:=predicted_bias]
                
                # in snp table
                ix <- het_table$ref_allele==ref & het_table$alt_allele==alt
                cat(sum(ix),'total \n')
                predicted_bias <- predict(thismodel,het_table[ix])
                cat(length(predicted_bias),'predicted \n')
                hist(predicted_bias)
                het_table[ix,bias_estimate:=predicted_bias]
            }
    
    snp_table[,AD_corrected:=AD*(2^bias_estimate)]
    snp_table[,DP_corrected:=RD+AD_corrected]
    snp_table[,allele_ratio_corrected:=AD_corrected/DP_corrected]
    
    het_table[,AD_corrected:=AD*(2^bias_estimate)]
    het_table[,DP_corrected:=RD+AD_corrected]
    het_table[,allele_ratio_corrected:=AD_corrected/DP_corrected]
    
    # Allele bias model, common SNPs ------------------------------------------------------------
    
    # How many heterozygous in reference of each SNP?
    het_table[,het_in_ref:=.N,by=c('id','ref_allele','alt_allele')]
    
    hist(het_table$het_in_ref[1:100000],100)
    
    het_table[,residual_bias:=log2((1-allele_ratio_corrected)/allele_ratio_corrected)]
    
    
    het_table[,constant:=0]
    het_table[het_in_ref>=5,constant:=mean(residual_bias),by=snp]
    #het_table[het_in_ref>=20,constant:=peakx(residual_bias),by=snp]
    
    
    snp_coeff_table <- unique(het_table[het_in_ref>=5,.(snp=paste(id,ref_allele,alt_allele),het_in_ref,constant,logDP=0,peak_ld=0,homC_T=0,homT_C=0,homA_G=0,homG_A=0)])
    
    for (i in which(snp_coeff_table$het_in_ref>20)) {
        data <- het_table[snp==snp_coeff_table$snp[i]]
        model <- rlm(residual_bias ~ logDP+peak_ld+homC_T+homT_C+homA_G+homG_A, data=data)
        coeff <- model$coefficients
        snp_coeff_table[i,constant:=coeff['(Intercept)']]
        snp_coeff_table[i,logDP:=coeff['logDP']]
        snp_coeff_table[i,peak_ld:=coeff['peak_ld']]
        snp_coeff_table[i,homC_T:=coeff['homC_T']]
        snp_coeff_table[i,homT_C:=coeff['homT_C']]
        snp_coeff_table[i,homA_G:=coeff['homA_G']]
        snp_coeff_table[i,homG_A:=coeff['homG_A']]
    }
    
    # Now predict the residual bias for all SNPs in treaining data (for sanity)
    ix <- match(het_table$snp,snp_coeff_table$snp)
    het_table[,residual_bias_estimate:=0]
    het_table[,residual_bias_estimate:=
                  snp_coeff_table$constant[ix]+
                  snp_coeff_table$logDP[ix]*logDP+
                  snp_coeff_table$peak_ld[ix]*peak_ld+
                  snp_coeff_table$homC_T[ix]*homC_T+
                  snp_coeff_table$homT_C[ix]*homT_C+
                  snp_coeff_table$homA_G[ix]*homA_G+
                  snp_coeff_table$homG_A[ix]*homG_A]
    het_table[is.na(residual_bias_estimate),residual_bias_estimate:=0]
    
    # And in full data for use
    ix <- match(snp_table$snp,snp_coeff_table$snp)
    snp_table[,residual_bias_estimate:=0]
    snp_table[,residual_bias_estimate:=
                  snp_coeff_table$constant[ix]+
                  snp_coeff_table$logDP[ix]*logDP+
                  snp_coeff_table$peak_ld[ix]*peak_ld+
                  snp_coeff_table$homC_T[ix]*homC_T+
                  snp_coeff_table$homT_C[ix]*homT_C+
                  snp_coeff_table$homA_G[ix]*homA_G+
                  snp_coeff_table$homG_A[ix]*homG_A]
    snp_table[is.na(residual_bias_estimate),residual_bias_estimate:=0]
    
    
    
    # Correct for residual bias
    het_table[,AD_corrected2:=AD_corrected*(2^residual_bias_estimate)]
    het_table[,DP_corrected2:=RD+AD_corrected2]
    het_table[,allele_ratio_corrected2:=AD_corrected2/DP_corrected2]
    
    snp_table[,AD_corrected2:=AD_corrected*(2^residual_bias_estimate)]
    snp_table[,DP_corrected2:=RD+AD_corrected2]
    snp_table[,allele_ratio_corrected2:=AD_corrected2/DP_corrected2]
    
    # Compute the correction factor and cap it at 2/3 to 3/2. Note that only the model based, first correction factor is used for that.
    snp_table[,correct_factor:=2^(bias_estimate+residual_bias_estimate)]
    snp_table[correct_factor>1.5,correct_factor:=1.5]
    snp_table[correct_factor<2/3,correct_factor:=2/3]
    
    # Compute the correction number
    snp_table[,correct_number:=(correct_factor-1)*AD]
    snp_table[,correct_number_check:=AD_corrected2-AD] # should be same
    
}

# Iterate samples ------------------------------------------------------------

targetlist <- NULL
for (i in 1:length(allcounts)) {
    
    counts <- allcounts[[i]]
    targets <- copy(target_template)
    name <- basename(counts$input_bam_file)
    targets$sample <- name

    
    # Add counts
    targets[,count:=counts$count]
    targets[,count_short:=counts$count_short]

    targets[,snps:=0]
    targets[,allele_count_correction:=0]
    
    # Add SNP info if present
    if (exists('snp_table')) {
        
        # Match vcf names in SNP table to this BAM file name
        this_table <- snp_table[str_detect(name,sample)]
        # Proceed if one matching name
        if (nrow(this_table)>0) if (length(unique(this_table$sample))==1) {
            snps_by_bin <- this_table[,.N,by=bin]
            targets[snps_by_bin$bin,snps:=snps_by_bin$N]
            correction_by_bin <- this_table[,max(correct_number),by=bin]
            targets[correction_by_bin$bin,allele_count_correction:=correction_by_bin$V1]
        }
        
    }
    
    targetlist[[i]] <- targets
    #backgroundlist[[i]] <- background
}

targets <- rbindlist(targetlist)
#background <- rbindlist(backgroundlist)
rm(targetlist)


# Drop some bins ------------------------------------------------------------


# Low coverage threshold
threshold <- median(targets$count) * 0.05
keep_targets <- targets[,median(count),by=bin][V1 > threshold]
targets <- targets[bin %in% keep_targets$bin]

# High coverage threshold
threshold <- median(targets$count) / 0.05
keep_targets <- targets[,median(count),by=bin][V1 < threshold]
targets <- targets[bin %in% keep_targets$bin]

# Mappability 0 removed
targets <- targets[is_target==F | map>0]


min1 <- function(data) {
    data[data<1] <- 1
    data[is.na(data)] <- 1
    return(data)
}

# Quantify variability for "blacklisting" high-variability regions.
if (wgs) {
    # sdev along each sample
    targets[,rollsd:= frollapply(log2(min1(count)), 10, sd, align = 'center',na.rm=T), by = sample]
    # medain over samples, for each bin
    targets[,rollsd_median:= median(rollsd), by = bin]
    # standard deviation of that
    sdev <- sd(targets$rollsd_median,na.rm = T)
    # drop bins with more noise
    targets <- targets[!rollsd_median > sdev*4]
}

# SNP and Median correct ------------------------------------------------------------



# Basic logR, targets (with correction for SNPs if available)
targets[,rawLR:=log2(min1(count+allele_count_correction))]
targets[,rawLR_short:=log2(min1(count_short+allele_count_correction*(count_short/count)))] # allelic correction scaled for short-fragments

# Median correct
targets[,rawLR:=rawLR-median(rawLR[is_backbone]),by=c('sample','is_target')]
targets[,rawLR_short:=rawLR_short-median(rawLR_short[is_backbone]),by=c('sample','is_target')]



# X-Y chromosome correct ------------------------------------------------------------

# Double X values in samples where their median implies male
targets[,xmedian:=median(rawLR[chromosome=='X']),by=sample]
targets[,male:=2^xmedian < .75] # assign gender

targets[,nonPA:=chromosome %in% c('X') & end>2.70e6 & start<154.93e6] # hard coded PA

if (length(unique(targets[male==T]$sample))>0) { # if at least 1 male
    targets[chromosome=='X' & male==T & nonPA,rawLR:=rawLR+1]
    targets[chromosome=='X' & male==T & nonPA,rawLR_short:=rawLR_short+1]
}

# Double Y values where their median implies male
targets[,ymedian:=median(rawLR[chromosome=='Y']),by=sample] # median by sample
targets[,male:=2^ymedian > .25] # assign gender based on Y

if (length(unique(targets[male==T]$sample))>0) { # if at least 1 male
    targets[male==T & chromosome=='Y' & end<28.79e6,rawLR:=rawLR+1] # hard coded PA
    targets[male==T & chromosome=='Y' & end<28.79e6,rawLR_short:=rawLR_short+1]
}

# Y is NA where median implies female
targets[chromosome=='Y' & male==F,rawLR:=NA]
targets[chromosome=='Y' & male==F,rawLR_short:=NA]




# Reference set median correct ------------------------------------------------------------
targets[,refmedian:=median(rawLR,na.rm=T),by=bin][,rawLR:=rawLR-refmedian]
targets[,refmedian_short:=median(rawLR_short,na.rm=T),by=bin][,rawLR_short:=rawLR_short-refmedian_short]




# Impute missing ----------------------------------------------------------

# If any missing, replace with random value near 0 (which is already the median)
targets[is.na(rawLR),rawLR:=rnorm(n = .N,mean = 0,sd = .1)]
targets[is.na(rawLR_short),rawLR_short:=rnorm(n = .N,mean = 0,sd = .1)]


# Matrix form ------------------------------------------------------------
tmat <- dcast(data = targets[,.(bin,sample,rawLR)],formula = bin ~ sample, value.var = 'rawLR')
tmat_short <- dcast(data = targets[,.(bin,sample,rawLR_short)],formula = bin ~ sample, value.var = 'rawLR_short')




# Reference object ------------------------------------------------------------
if (wgs) allcounts[[1]]$target_bed_file <- 'wgs'
reference <- allcounts[[1]][c("target_bed_file","chromlength","ranges")]
reference$date <- date()
reference$samples <- unique(targets$sample)
reference$targets <- target_template
reference$keep <- unique(targets$bin)

reference$targets_ref <- tmat #tpca$x
reference$targets_ref_short <- tmat_short #tpca_short$x



reference$median <- targets[sample==sample[1]]$refmedian
reference$median_short <- targets[sample==sample[1]]$refmedian_short

if (exists('snp_table')) {
    reference$snp_rlm_model <- snp_rlm_model
    reference$snp_coeff_table <- snp_coeff_table
}

# Save ------------------------------------------------------------

name <- reference$target_bed_file
if (is.null(name)) name <- 'jumble.WGS'
saveRDS(reference,paste0(opt$output_folder,'/',
                         str_remove(name,'.*/'),
                         '.reference.RDS'))


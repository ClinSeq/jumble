





# Dependencies and arguments ----------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(RJSONIO))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(data.table))

#long, short(NA), argmask, datatype, desc
#argmask 0=no arg, 1=req, 2=optional
args <- rbind(
  c("tumor_cnr", NA, 1, "character", "tumor bin file from CNVkit"),
  c("tumor_cns", NA, 1, "character", "tumor segment file from CNVkit"),
  c("normal_cnr", NA, 1, "character", "normal bin file from CNVkit"),
  c("normal_cns", NA, 1, "character", "normal segment file from CNVkit"),
  c("het_snps_vcf", NA, 1, "character", "heterozygous SNPs .vcf file"),
  c("purecn_csv", NA, 1, "character", "PureCN result .csv file"),
  c("purecn_genes_csv", NA, 1, "character", "PureCN result _genes.csv"),
  c("purecn_loh_csv", NA, 1, "character", "PureCN result _loh.csv"),
  c("purecn_variants_csv", NA, 1, "character", "PureCN result _variants.csv"),
  c("svcaller_T_DEL", NA, 1, "character", "Tumor SV caller DEL-events.gtf"),
  c("svcaller_T_DUP", NA, 1, "character", "Tumor SV caller DUP-events.gtf"),
  c("svcaller_T_INV", NA, 1, "character", "Tumor SV caller INV-events.gtf"),
  c("svcaller_T_TRA", NA, 1, "character", "Tumor SV caller TRA-events.gtf"),
  c("svcaller_N_DEL", NA, 1, "character", "Normal SV caller DEL-events.gtf"),
  c("svcaller_N_DUP", NA, 1, "character", "Normal SV caller DUP-events.gtf"),
  c("svcaller_N_INV", NA, 1, "character", "Normal SV caller INV-events.gtf"),
  c("svcaller_N_TRA", NA, 1, "character", "Normal SV caller TRA-events.gtf"),
  c("germline_mut_vcf", NA, 1, "character", "germline mutation vcf file"),
  c("somatic_mut_vcf", NA, 1, "character", "somatic mutation vcf file"),
  c("plot_png", NA, 1, "character", "plot .png file name"),
  c("plot_png_normal", NA, 1, "character", "normal CNV plot .png file name"),
  c("cna_json", NA, 1, "character", "CNA output json file name"),
  c("purity_json", NA, 1, "character", "purity output json file name"),
  c("gene_track", NA, 1, "character", "gene_track csv file name")
)


opts <- getopt(args)

save.image('ws.Rdata')


chrsizes <- data.table(
    chr = c("1", "2", "3", "4", "5", "6", "7", "8",
            "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
            "20", "21", "22", "X", "Y", "MT"),
    size = c(249250621L, 243199373L,
             198022430L, 191154276L, 180915260L, 171115067L, 159138663L, 146364022L,
             141213431L, 135534747L, 135006516L, 133851895L, 115169878L, 107349540L,
             102531392L, 90354753L, 81195210L, 78077248L, 59128983L, 63025520L,
             48129895L, 51304566L, 155270560L, 59373566L, 16569L),
    cumstart = c(0,
                 249250621, 492449994, 690472424, 881626700, 1062541960, 1233657027,
                 1392795690, 1539159712, 1680373143, 1815907890, 1950914406, 2084766301,
                 2199936179, 2307285719, 2409817111, 2500171864, 2581367074, 2659444322,
                 2718573305, 2781598825, 2829728720, 2881033286, 3036303846, 3095677412
    )
  )



# Genes to plot -----------------------------------------------------------
genes <- data.table(label = c("AKT1", "APC", "AR", "ARID1A", "ARID2", "ATM", "ATR", "BARD1",
                           "BRAF", "BRCA1", "BRCA2", "BRIP1", "CCND1", "CDH1", "CDK12", "CDK4",
                           "CDK6", "CDKN1A", "CDKN1B", "CDKN2A", "CHD1", "CHEK2", "CTNNB1",
                           "CUL3", "DICER1", "DNMT3A", "FANCA", "FOXA1", "FOXO1", "HRAS", "IDH1",
                           "JAK1", "KDM6A", "KEAP1", "KMT2A", "KMT2C", "KMT2D", "KRAS", "MED12",
                           "MET", "MGA", "MLH1", "MLH3", "MRE11A", "MSH2", "MSH3", "MSH6", "MYC",
                           "NBN", "NCOR1", "NKX3-1", "NRAS", "PALB2", "PIK3CA", "PIK3CB", "PIK3CD",
                           "PIK3R1", "PIK3R2", "PMS1", "PMS2", "POLD1", "POLE", "PTEN", "RAD50",
                           "RAD51", "RAD51B", "RAD51C", "RAD51D", "RB1", "RNF43", "SETD2", "SF3B1",
                           "SPEN", "SPOP", "TMPRSS2", "TP53", "U2AF1", "XPO1", "ZBTB16", "ZFHX3",
                           "ZMYM3"),
                 chromosome = c("14", "5", "X", "1", "12", "11", "3", "2", "7", "17", "13", "17",
                                "11", "16", "17", "12", "7", "6", "12", "9", "5", "22", "3", "2",
                                "14", "2", "16", "14", "13", "11", "2", "1", "X", "19", "11", "7",
                                "12", "12", "X", "7", "15", "3", "14", "11", "2", "5", "2", "8",
                                "8", "17", "8", "1", "16", "3", "3", "1", "5", "19", "2", "7", "19",
                                "12", "10", "5", "15", "14", "17", "17", "13", "17", "3", "2", "1",
                                "17", "21", "17", "21", "2", "11", "16", "X"),
                 start = c(105235686, 112043195, 66764465, 27022524, 46123448, 108093211, 142168077,
                           215590370, 140419127, 41196312, 32889611, 59758627, 69455855, 68771128,
                           37617764, 58141510, 92234235, 36644305, 12867992, 21967751, 98190908,
                           29083731, 41236328, 225334867, 95552565, 25455845, 89803957, 38059189,
                           41129804, 532242, 209100951, 65298912, 44732757, 10596796, 118307205,
                           151832010, 49412758, 25357723, 70338406, 116312444, 41913422, 37034823,
                           75480467, 94152895, 47630108, 79950467, 47922669, 128747680, 90945564,
                           15932471, 23536206, 115247090, 23614488, 178865902, 138372860, 9711790,
                           67511548, 18263928, 190649107, 6012870, 50887461, 133200348, 89622870,
                           131891711, 40986972, 68286496, 56769934, 33426811, 48877887, 56429861,
                           47057919, 198254508, 16174359, 47676246, 42836478, 7565097, 44513066,
                           61704984, 113930315, 72816784, 70459474),
                 end = c(105262088, 112181936, 66950461, 27108595, 46301823, 108239829, 142297668,
                         215674428, 140624564, 41277500, 32973805, 59940882, 69469242, 68869451,
                         37721160, 58149796, 92465908, 36655116, 12875305, 21995300, 98262240,
                         29138410, 41301587, 225450110, 95624347, 25565459, 89883065, 38069245,
                         41240734, 537287, 209130798, 65432187, 44971847, 10614417, 118397539,
                         152133090, 49453557, 25403870, 70362303, 116438440, 42062141, 37107380,
                         75518235, 94227074, 47789450, 80172279, 48037240, 128753674, 91015456,
                         16121499, 23540440, 115259515, 23652631, 178957881, 138553780, 9789172,
                         67597649, 18281350, 190742355, 6048756, 50921273, 133263951, 89731687,
                         131980313, 41024354, 69196935, 56811703, 33448541, 49056122, 56494956,
                         47205457, 198299815, 16266955, 47755596, 42903043, 7590856, 44527697,
                         61765761, 114121398, 73093597, 70474996)
                 )
genes$cumstart <- genes$start + chrsizes$cumstart[match(genes$chromosome,chrsizes$chr)]
genes$cumend <- genes$end + chrsizes$cumstart[match(genes$chromosome,chrsizes$chr)]


# Read gene track ---------------------------------------------------------

data2 <- read.table(opts$gene_track, sep=",", na.strings=".", stringsAsFactors=FALSE, quote="", fill=FALSE , comment.char="#", header=TRUE)
data2$cumstart <- data2$Start + chrsizes$cumstart[match(data2$Chromosome,chrsizes$chr)]
data2$cumend <- data2$End + chrsizes$cumstart[match(data2$Chromosome,chrsizes$chr)]
data2_mRNA <- data2[data2$Feature == 'transcript',]
data2_3p <- data2[data2$Feature == 'three_prime_utr',]
data2_5p <- data2[data2$Feature == 'five_prime_utr',]
data2_exon <- data2[data2$Feature == 'exon',]







# Read SNP allele ratio ---------------------------------------------------

{
  vcf <- readVcf(opts$het_snps_vcf,genome = "GRCh37")
  g <- geno(vcf)
  r <- rowRanges(vcf)
  chr <- as.character(seqnames(r))
  pos <- data.frame(ranges(r))$start
  alf <- NULL
  smoothedAi <- NA
  if (length(pos)>0) {
    pos=as.numeric(pos)
    alf <- data.frame(chromosome=chr, start=pos, end=pos, stringsAsFactors = F, cumstart=NA, cumend=NA)
  }

  if (!is.null(alf)) if( ! all(is.na(alf$chromosome)) ) {

    for(chr in chrsizes$chr){
      ix <- which(alf$chromosome == chr)
      alf$cumstart[ix] <- alf$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
      alf$cumend[ix] <- alf$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    }

    alf$t <- as.numeric(sapply(g$AD[,2], "[", 2))/as.numeric(g$DP[,2])
    alf$t[is.nan(alf$t)]=NA # allele freq becomes NaN if cov=0. Then set to NA
    alf$n <- as.numeric(sapply(g$AD[,1], "[", 2))/as.numeric(g$DP[,1])
    alf$n[is.nan(alf$n)]=NA # allele freq becomes NaN if cov=0. Then set to NA
    alf$td <- as.numeric(g$DP[,2])
    alf$nd <- as.numeric(g$DP[,1])

    alf$ai=2*abs(alf$t-0.5)
    alf$ai_n=2*abs(alf$n-0.5)
  }
}



# Read somatic point mutations --------------------------------------------


{
  vcf <- readVcf(opts$somatic_mut_vcf,genome = "GRCh37")
  g <- geno(vcf)
  r=rowRanges(vcf)
  if (length(g)>0) { # if there are any somatic mutations...
    chr=as.character(seqnames(r))
    pos=data.frame(ranges(r))$start
    salf <- data.frame(N=1:length(chr),chromosome=chr,pos=pos,stringsAsFactors = F)
    rownames(salf)=names(r)
    salf$REF=as.data.frame(r$REF)[,1]
    salf$ALT=as.data.frame(r$ALT)[,3]
    salf$cumpos <- NA
    for(chr in chrsizes$chr){
      ix <- which(salf$chromosome == chr)
      salf$cumpos[ix] <- salf$pos[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    }

    salf$FILTER <- fixed(vcf)$FILTER

    salf$AF.T <- as.numeric(g$VAF[,2])
    salf$AO.T <- as.numeric(apply(g$DP4[,2,3:4],1,sum))  # sum alt forward and alt reverse
    salf$DP.T <- as.numeric(apply(g$DP4[,2,],1,sum))  # sum ref forward, ref reverse, alt forward and alt reverse
    salf$AO.N <- as.numeric(apply(g$DP4[,1,3:4],1,sum))  # sum alt forward and alt reverse
    salf$DP.N <- as.numeric(apply(g$DP4[,1,],1,sum))  # sum ref forward, ref reverse, alt forward and alt reverse

    salf$type='other'
    salf$type[isSNV(vcf)]='snv'
    salf$type[isDeletion(vcf)]='del'
    salf$type[isInsertion(vcf)]='ins'


    header=info(header(vcf))$Description
    ix=grep('Consequence annotations from Ensembl',header)
    header=strsplit(header[ix],'\\|')[[1]]
    header[1]='Allele'

    vep=info(vcf)$CSQ

    ## This loop creates a new "table" with all mutation effects.
    rowspermut=unlist(lapply(vep,length))
    table=#data.frame(
      matrix('',nrow = sum(rowspermut),ncol = length(header)+1)#,stringsAsFactors = F)
    colnames(table)=c('N',header)  # blir detta fel ibland?????
    for (i in 1:length(vep)) {
      #if (i %% 100 ==0) cat(i,'..')
      for (j in 1:length(vep[[i]])) { # for each effect
        t2=vep[[i]][j]
        t2=strsplit(t2,'[|]')[[1]] # pipe separated line with one effect of the mutation
        t2=c((i),t2)
        thisrow=sum(rowspermut[1:i])-rowspermut[i]+j
        table[thisrow,1:length(t2)]=t2
      }
    }

    salf=merge(salf,table,by='N',all=T)

    # Remove rejected variants
    salf = subset(salf, FILTER != "REJECT")

  } #end somatic mutations
  salf=(salf[,-1])

  # mark the type
  salf$pch=rep(22,nrow(salf))
  salf$pch[salf$type=='snv']=21
  salf$pch[salf$type=='del']=24
  salf$pch[salf$type=='ins']=25

  # Icke-NA/intron på konsekvens
  salf$hasConsequence=!salf$Consequence %in% c("intron_variant", "synonymous_variant",
                                               "splice_region_variant&intron_variant",
                                               "3_prime_UTR_variant", "intergenic_variant",
                                               "regulatory_region_variant", "upstream_gene_variant",
                                               "downstream_gene_variant",
                                               "intron_variant&non_coding_transcript_variant",
                                               "5_prime_UTR_variant", "splice_region_variant&synonymous_variant",
                                               "non_coding_transcript_exon_variant&non_coding_transcript_variant",
                                               "intron_variant&NMD_transcript_variant", "TF_binding_site_variant",
                                               "splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant",
                                               "mature_miRNA_variant", "coding_sequence_variant&5_prime_UTR_variant")
  salf$hasConsequence[is.na(salf$hasConsequence)]=F


  ## Vilka är tydligt deleterious? (inkl inframe)
  salf$deleterious=rep(F,nrow(salf))
  salf$deleterious[unique(c(
    grep('start',salf$Consequence),
    grep('stop',salf$Consequence),
    grep('frame',salf$Consequence),
    grep('acceptor',salf$Consequence),
    grep('donor',salf$Consequence)
  ))]=T
  salf$deleterious[!salf$hasConsequence]=F

}


# Read germline point mutations -------------------------------------------


{
  vcf <- readVcf(opts$germline_mut_vcf,genome = "GRCh37")
  vcf <- expand(vcf)
  g <- geno(vcf)
  r=rowRanges(vcf)
  f <- fixed(vcf)
  if (length(g)>0) { # if there are any somatic mutations...
    chr=as.character(seqnames(r))
    pos=data.frame(ranges(r))$start
    galf <- data.frame(N=1:length(chr),chromosome=chr,pos=pos,stringsAsFactors = F)
    rownames(galf)=names(r)
    galf$cumpos <- NA
    galf$REF <- as.character(f$REF)
    galf$ALT <- as.character(f$ALT)
    for(chr in chrsizes$chr){
      ix <- which(galf$chromosome == chr)
      galf$cumpos[ix] <- galf$pos[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    }

    galf$AF <- g$AD[,1,2]/g$DP[,1]
    #galf$AF_n <- g$AD[,2,2]/g$DP[,2]    # <----- this should be the normal sample allele ratio.
    galf$AO <- g$AD[,1,2]
    galf$DP <- g$DP[,1]

    galf$type='other'
    galf$type[isSNV(vcf)]='snv'
    galf$type[isDeletion(vcf)]='del'
    galf$type[isInsertion(vcf)]='ins'


    header=info(header(vcf))$Description
    ix=grep('Consequence annotations from Ensembl',header)
    header=strsplit(header[ix],'\\|')[[1]]
    header[1]='Allele'

    vep=info(vcf)$CSQ

    ## This loop creates a new "table" with all mutation effects.
    rowspermut=unlist(lapply(vep,length))
    table=#data.frame(
      matrix('',nrow = sum(rowspermut),ncol = length(header)+1)#,stringsAsFactors = F)
    colnames(table)=c('N',header)  # blir detta fel ibland?????
    for (i in 1:length(vep)) {
      for (j in 1:length(vep[[i]])) { # for each effect
        t2=vep[[i]][j]
        t2=strsplit(t2,'[|]')[[1]] # pipe separated line with one effect of the mutation
        t2=c((i),t2)
        thisrow=sum(rowspermut[1:i])-rowspermut[i]+j
        table[thisrow,1:length(t2)]=t2
      }
    }
    table=table[,-which(colnames(table)=="AF")]  # remove AF col from vep data to not confuse with AF from galf

    # Add alleles as vep names them, to use for merge
    galf$Allele = ifelse(galf$type=="snv", yes = galf$ALT, no = substr(galf$ALT, start=2, stop=nchar(galf$ALT)))
    galf$Allele[which(galf$Allele== "")] = "-"
    # Merge with vep table
    galf=merge(galf,table,by=c('N','Allele'),all.x=T)
  } #end germline mutations
  galf$gnomAD_AF=as.numeric(galf$gnomAD_AF)  # Make gnom_AD numerical so it can be used
  galf=(galf[which(galf$AO>=12 & galf$AF>=0.2 & (galf$gnomAD_AF < 0.05 | is.na(galf$gnomAD_AF))),-1])
  # mark the type
  galf$pch=rep(22,nrow(galf))
  galf$pch[galf$type=='snv']=21
  galf$pch[galf$type=='del']=24
  galf$pch[galf$type=='ins']=25

  # Icke-NA/intron på konsekvens
  galf$hasConsequence=!galf$Consequence %in% c("intron_variant", "synonymous_variant",
                                               "splice_region_variant&intron_variant",
                                               "3_prime_UTR_variant", "intergenic_variant",
                                               "regulatory_region_variant", "upstream_gene_variant",
                                               "downstream_gene_variant",
                                               "intron_variant&non_coding_transcript_variant",
                                               "5_prime_UTR_variant", "splice_region_variant&synonymous_variant",
                                               "non_coding_transcript_exon_variant", "non_coding_transcript_variant",
                                               "non_coding_transcript_exon_variant&non_coding_transcript_variant",
                                               "intron_variant&NMD_transcript_variant", "TF_binding_site_variant",
                                               "splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant",
                                               "mature_miRNA_variant", "coding_sequence_variant&5_prime_UTR_variant")
  galf$hasConsequence[is.na(galf$hasConsequence)]=F


  ## Vilka är tydligt deleterious? (inkl inframe)
  galf$deleterious=rep(F,nrow(galf))
  galf$deleterious[unique(c(
    grep('start',galf$Consequence),
    grep('stop',galf$Consequence),
    grep('frame',galf$Consequence),
    grep('acceptor',galf$Consequence),
    grep('donor',galf$Consequence)
  ))]=T
  galf$deleterious[!galf$hasConsequence]=F



}




# Read CNVkit copy number data --------------------------------------------


{ #tumor
  segments <- read.delim(opts$tumor_cns,stringsAsFactors = F)
  bins <- read.delim(opts$tumor_cnr,stringsAsFactors = F)

  ## Segment start and end pos
  segments$cumstart <- NA
  segments$cumend   <- NA
  for(chr in chrsizes$chr){
    ix <- which(segments$chromosome == chr)
    segments$cumstart[ix] <- segments$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    segments$cumend[ix] <- segments$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
  }
  segments$centerpos <- segments$cumstart+(segments$cumend-segments$cumstart)/2

  ## Bin start and end pos
  bins$cumstart <- NA
  bins$cumend   <- NA
  for(chr in chrsizes$chr){
    ix <- which(bins$chromosome == chr)
    bins$cumstart[ix] <- bins$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    bins$cumend[ix] <- bins$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
  }
  bins$centerpos <- bins$cumstart+(bins$cumend-bins$cumstart)/2
  bins=bins[order(bins$cumstart),]
}

{ # normal
  segments_n <- read.delim(opts$normal_cns,stringsAsFactors = F)
  bins_n <- read.delim(opts$normal_cnr,stringsAsFactors = F)

  ## Segment start and end pos
  segments_n$cumstart <- NA
  segments_n$cumend   <- NA
  for(chr in chrsizes$chr){
    ix <- which(segments_n$chromosome == chr)
    segments_n$cumstart[ix] <- segments_n$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    segments_n$cumend[ix] <- segments_n$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
  }
  segments_n$centerpos <- segments_n$cumstart+(segments_n$cumend-segments_n$cumstart)/2

  ## Bin start and end pos
  bins_n$cumstart <- NA
  bins_n$cumend   <- NA
  for(chr in chrsizes$chr){
    ix <- which(bins_n$chromosome == chr)
    bins_n$cumstart[ix] <- bins_n$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    bins_n$cumend[ix] <- bins_n$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
  }
  bins_n$centerpos <- bins_n$cumstart+(bins_n$cumend-bins_n$cumstart)/2
  bins_n=bins_n[order(bins_n$cumstart),]
}

### Read structural variant files
{ # tumor
  t_strvs=NULL
  try( {
    sv <- read.delim(opts$svcaller_T_DEL,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='DEL'
    t_strvs=rbind(t_strvs,sv)
  }, silent=T)
  try( {
    sv <- read.delim(opts$svcaller_T_DUP,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='DUP'
    t_strvs=rbind(t_strvs,sv)
  }, silent=T)
  try( {
    sv <- read.delim(opts$svcaller_T_INV,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='INV'
    t_strvs=rbind(t_strvs,sv)
  }, silent=T)
  try( {
    sv <- read.delim(opts$svcaller_T_TRA,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='TRA'
    t_strvs=rbind(t_strvs,sv)
  }, silent=T)
  t_strvs$cumstart=NA
  t_strvs$cumend=NA
  for(i in 1:nrow(chrsizes)){
    ix <- which(t_strvs$chr == chrsizes$chr[i])
    t_strvs$cumstart[ix] <- t_strvs$start[ix] + chrsizes$cumstart[i]
    t_strvs$cumend[ix] <- t_strvs$end[ix] + chrsizes$cumstart[i]
  }
}
{ # normal
  n_strvs=NULL
  try( {
    sv <- read.delim(opts$svcaller_N_DEL,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='DEL'
    n_strvs=rbind(n_strvs,sv)
  }, silent=T)
  try( {
    sv <- read.delim(opts$svcaller_N_DUP,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='DUP'
    n_strvs=rbind(n_strvs,sv)
  }, silent=T)
  try( {
    sv <- read.delim(opts$svcaller_N_INV,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='INV'
    n_strvs=rbind(n_strvs,sv)
  }, silent=T)
  try( {
    sv <- read.delim(opts$svcaller_N_TRA,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='TRA'
    n_strvs=rbind(n_strvs,sv)
  }, silent=T)
  n_strvs$cumstart=NA
  n_strvs$cumend=NA
  for(i in 1:nrow(chrsizes)){
    ix <- which(n_strvs$chr == chrsizes$chr[i])
    n_strvs$cumstart[ix] <- n_strvs$start[ix] + chrsizes$cumstart[i]
    n_strvs$cumend[ix] <- n_strvs$end[ix] + chrsizes$cumstart[i]
  }
}



# Read PureCN files -------------------------------------------------------


{
  # files to read (purity/ploidy, mutations and snps, gene copy number and LOH, segmented copy number and LOH)
  purecn_files = c(opts$purecn_csv, opts$purecn_variants_csv, opts$purecn_genes_csv, opts$purecn_loh_csv)
  # variables to assign to
  purecn_variables = c("purecn_stat", "purecn_vars", "purecn_genes", "purecn_loh")
  # colnames to use if no data avialable and mock df created
  purecn_colnames = list(
    c("Sampleid","Purity","Ploidy","Sex","Contamination","Flagged","Failed","Curated","Comment"),
    c("Sampleid", "chr", "start", "end", "ID", "REF", "ALT", "SOMATIC.M0", "SOMATIC.M1", "SOMATIC.M2",
      "SOMATIC.M3", "SOMATIC.M4", "SOMATIC.M5", "SOMATIC.M6", "SOMATIC.M7", "GERMLINE.M0", "GERMLINE.M1",
      "GERMLINE.M2", "GERMLINE.M3", "GERMLINE.M4", "GERMLINE.M5", "GERMLINE.M6", "GERMLINE.M7", "GERMLINE.CONTHIGH",
      "GERMLINE.CONTLOW", "GERMLINE.HOMOZYGOUS", "ML.SOMATIC", "POSTERIOR.SOMATIC", "ML.M", "ML.C", "ML.M.SEGMENT",
      "M.SEGMENT.POSTERIOR", "M.SEGMENT.FLAGGED", "ML.AR", "AR", "AR.ADJUSTED", "MAPPING.BIAS", "ML.LOH",
      "CN.SUBCLONAL", "CELLFRACTION", "FLAGGED", "log.ratio", "depth", "prior.somatic", "prior.contamination",
      "on.target", "seg.id", "gene.symbol"),
    c("Sampleid", "gene.symbol", "chr", "start", "end", "C", "seg.mean", "seg.id", "number.targets", "gene.mean",
      "gene.min", "gene.max", "focal", "breakpoints", "type", "num.snps.segment", "M", "M.flagged", "loh"),
    c("Sampleid", "chr", "start", "end", "arm", "C", "M", "type")
  )
  # sample id to use if no data avialable and mock df created
  sampid = sub(".csv","", basename(purecn_files[1]))  # NB: is this correct, or should there also be the "-nodups" ending to conform with samples with PureCN data?
  for (i in 1:length(purecn_files)) {
    if (file.size(purecn_files[i]) == 0) {  # if no valid PureCN output produced, create mock df
      assign(purecn_variables[i],
             data.frame(matrix(c(sampid, rep(NA, length(purecn_colnames[[i]])-1)), nrow=1,
                               dimnames=list(1, purecn_colnames[[i]])), stringsAsFactors = F))
    } else {  # if data available in the PureCN output, read it
      assign(purecn_variables[i], read.delim(purecn_files[i], sep = ',', stringsAsFactors = F))
    }
  }
  purecn_loh$cumstart=NA
  purecn_loh$cumend=NA
  for(i in 1:nrow(chrsizes)){
    ix <- which(purecn_loh$chr == chrsizes$chr[i])
    if (length(ix) == 0) next()  # skip if the chromosome is not present in the loh table
    purecn_loh$cumstart[ix] <- purecn_loh$start[ix] + chrsizes$cumstart[i]
    purecn_loh$cumend[ix] <- purecn_loh$end[ix] + chrsizes$cumstart[i]
  }
}




setwd('C:/Users/schon035/Documents/Projects/taxonomy/ath_bna/')

ath_bed = read.table('../../genome_resources/TAIR10.51.bed')
ath_order = unique(gsub('\\..*','',ath_bed[ath_bed[,5]=='protein_coding',4]))

ack = read.table('../../genome_resources/ACK_blocks_nogaps.tsv', header = T, sep='\t', comment.char = '')
ack_color_matrix = ack[,c("block","hex")]
ack_color_matrix = ack_color_matrix[!duplicated(ack_color_matrix),]
ack_colors = ack_color_matrix[,2]
rownames(ack) = ack[,1]
names(ack_colors) = ack_color_matrix[,1]
rm(ack_color_matrix)
ack_colors = ack_colors[sort(names(ack_colors))]

bna_gff = read.table('BnapusDarmor-bzh_annotation.gff', sep='\t')
bna_genes = bna_gff[bna_gff[,3]=='mRNA',]
bna_genes[,'ID'] = gsub('ID=([^;]+);.*$','\\1', bna_genes[,9])


orthologs = read.table('Bnapus_orthologs.tsv', sep='\t', header = T)
orthologs = orthologs[!duplicated(orthologs),]

diamond_forward = read.table('results/BnapusDarmor_bzh_arabidopsis_thaliana_UP000006548.tsv', sep='\t')
diamond_reverse = read.table('results/arabidopsis_thaliana_UP000006548_BnapusDarmor_bzh.tsv', sep='\t')
uniprot_agi = read.table('../ath_uniprot.tsv', header = F)
diamond_forward[,2] = unlist(lapply(strsplit(diamond_forward[,2],'\\|'), function(x)x[2]))
diamond_reverse[,1] = unlist(lapply(strsplit(diamond_reverse[,1],'\\|'), function(x)x[2]))
diamond_forward[,'agi'] = sapply(diamond_forward[,2], function(i)uniprot_agi[uniprot_agi[,"V1"]==i,2][1])
diamond_reverse[,'agi'] = sapply(diamond_reverse[,1], function(i)uniprot_agi[uniprot_agi[,"V1"]==i,2][1])
agi_ack = ack[,"block"]
names(agi_ack) = ack[,"genes"]

tophit = sapply(bna_genes[,"ID"], function(gene)diamond_forward[which(diamond_forward[,1]==gene)[1],"agi"])
syn_hit_match = agi_ack[unlist(tophit)] == syntelogs[names(tophit), "consensus_ACK"]


block_likelihoods_forward = matrix(0, nrow=length(unique(orthologs[,1])), ncol=length(ack_colors), dimnames = list(c(unique(orthologs[,1])), names(ack_colors)))
block_likelihoods_reverse = matrix(0, nrow=length(unique(orthologs[,1])), ncol=length(ack_colors), dimnames = list(c(unique(orthologs[,1])), names(ack_colors)))
block_likelihoods_reverse = block_likelihoods_reverse[nrow(block_likelihoods_reverse):1,]

chromosomes = gsub('p.*','',gsub('^BNapus_Darmor_BZH_','',rownames(block_likelihoods_forward)))
names(chromosomes) = rownames(block_likelihoods_forward)

decay = .9
for(i in 1:nrow(block_likelihoods_forward)){
  gene = rownames(block_likelihoods_forward)[i]
  chrom = chromosomes[gene]
  hits = diamond_forward[diamond_forward[,1] == gene,]
  hits = hits[!is.na(hits[,"agi"]),]
  if(nrow(hits) > 0){
    for(j in 1:nrow(hits)){
      block = agi_ack[hits[j,'agi']]
      if(!is.na(block)){
        block_likelihoods_forward[i,block] = max(block_likelihoods_forward[i,block], hits[j,12])
      }
    }
  }
  if(i>1){
    if(last_chrom == chrom){
      # apply a decayed bit score from the previous row
      block_likelihoods_forward[i,] = block_likelihoods_forward[i,] + decay * block_likelihoods_forward[i-1,]
    }
  }
  last_chrom = chrom
}


for(i in 1:nrow(block_likelihoods_reverse)){
  gene = rownames(block_likelihoods_reverse)[i]
  chrom = chromosomes[gene]
  hits = diamond_forward[diamond_forward[,1] == gene,]
  hits = hits[!is.na(hits[,"agi"]),]
  if(nrow(hits) > 0){
    for(j in 1:nrow(hits)){
      block = agi_ack[hits[j,'agi']]
      if(!is.na(block)){
        block_likelihoods_reverse[i,block] = max(block_likelihoods_reverse[i,block], hits[j,12])
      }
    }
  }
  if(i>1){
    if(last_chrom == chrom){
      # apply a decayed bit score from the previous row
      block_likelihoods_reverse[i,] = block_likelihoods_reverse[i,] + decay * block_likelihoods_reverse[i-1,]
    }
  }
  last_chrom = chrom
}

block_likelihoods_reverse = block_likelihoods_reverse[nrow(block_likelihoods_reverse):1,]

block_likelihoods = sqrt(block_likelihoods_forward * block_likelihoods_reverse)
block_consensus = apply(block_likelihoods, 1, function(i)colnames(block_likelihoods)[order(i, decreasing = T)[1]])
block_consensus_score = apply(block_likelihoods, 1, function(i)sort(i, decreasing = T)[1])

orthologs[,'consensus_ACK'] = block_consensus[orthologs[,"Bna"]]
orthologs[,'consensus_score'] = round(block_consensus_score[orthologs[,"Bna"]],2)

orthologs[,"Ath_position"] = as.numeric(sapply(orthologs[,"Ath_agi"], function(agi) which(ath_order==agi)[1]))
# plot(diff(orthologs[,"Ath_position"])[1:10000], type='l', ylim=c(-10,10))
# abline(h=0)


out_of_place = orthologs[,"ACK_block"] != orthologs[,"consensus_ACK"]

to_fix = orthologs[which(out_of_place),]
fixed_syntelogs = apply(to_fix, 1, function(x){
  best = x["Ath_agi"]
  hits = diamond_forward[diamond_forward[,1] == x[1], "agi"]
  if(x["consensus_ACK"] %in% agi_ack[hits]){
    best = hits[agi_ack[hits]==x["consensus_ACK"]][1]
    cat(x["Bna"], ' -> ', best, '(hit #', which(hits==best),')\n',sep='')
  }
  return(best)
})
names(fixed_syntelogs) = to_fix[,"Bna"]
sum(fixed_syntelogs != to_fix[,"Ath_agi"], na.rm = T)


syntelogs = orthologs
rownames(syntelogs) = make.names(syntelogs[,"Bna"], unique = T)
syntelogs[names(fixed_syntelogs),"Ath_agi"] = fixed_syntelogs
syntelogs[names(fixed_syntelogs),"ACK_block"] = agi_ack[fixed_syntelogs]
syntelogs[,"Ath_position"] = as.numeric(sapply(syntelogs[,"Ath_agi"], function(agi) which(ath_order==agi)[1]))
syntelogs[,"Ath_increment"] = c(0,diff(syntelogs[,"Ath_position"]))

write.table(syntelogs, 'Bnapus_syntelogs_20240606.tsv', quote = F, sep='\t', row.names = F)

orthologs[,"Ath_increment"] = c(0,diff(orthologs[,"Ath_position"]))


write.table(orthologs[order(orthologs[,"Ath_position"]),], 'Bnapus_orthologs_ath_order.tsv', quote = F, sep='\t', row.names = F)
write.table(orthologs, 'Bnapus_orthologs.tsv', quote = F, sep='\t', row.names = F)

write.table(round(block_likelihoods,0), 'Bnapus_ACK_block_scores.tsv', quote = F, sep='\t', row.names = T)
block_likelihoods = read.table("Bnapus_ACK_block_scores.tsv", sep='\t')
columns = data.frame(
  'ack'=names(ack_colors)
)
anno_col = data.frame(
  'ack'=ack_colors
)
pheatmap(block_likelihoods, cluster_rows = F, cluster_cols = F, show_rownames = F, scale = 'row')
dev.off()

show_blocks <- function(genelist){
  hits = sapply(genelist, function(g)which(orthologs[,"Ath_agi"]==g))
  plot(NA, ylim=c(0,nrow(orthologs)), xlim=c(0,length(genelist)), xlab="Query Genes", ylab="Target genes")
  abline(v=sapply(1:5, function(i)min(grep(paste('AT',i,'G',sep=''),ath_order))))
  abline(h=sapply(unique(chromosomes), function(i)min(which(orthologs[,1] %in% names(which(chromosomes==i))))))
  copynumber = rep(0,length(genelist))
  names(copynumber) = genelist
  for(x in 1:length(genelist)){
    y = hits[[x]]
    points(x=rep(x,length(y)), y=y, pch='.', cex=3, col=ack_colors[orthologs[y,"consensus_ACK"]])
    copynumber[x] = length(y)
  }
  copynumber[copynumber == 0] = NA
}

show_blocks_flip <- function(genelist){
  genetable = syntelogs
  hits = sapply(genelist, function(g)which(genetable[,"Ath_agi"]==g))
  plot(NA, xlim=c(0,nrow(genetable)), ylim=c(0,length(genelist)), ylab="A.thaliana genes", xlab="B.napus genes")
  abline(h=sapply(1:5, function(i)min(grep(paste('AT',i,'G',sep=''),ath_order))))
  abline(v=sapply(unique(chromosomes), function(i)min(which(genetable[,1] %in% names(which(chromosomes==i))))))
  copynumber = rep(0,length(genelist))
  names(copynumber) = genelist
  for(x in 1:length(genelist)){
    y = hits[[x]]
    points(y=rep(x,length(y)), x=y, pch='.', cex=3, col=ack_colors[genetable[y,"consensus_ACK"]])
    copynumber[x] = length(y)
  }
  copynumber[copynumber == 0] = NA
}


show_blocks_flip(ath_order[length(ath_order):1])
title(main = 'OrthoFinder gene order (synteny correction)')

### THE BLOCKS ###

blocks = list(
  ### A ###
  c('A10p00900','A10p00060'),
  c('A10p00930', 'A10p06130'),
  c('A06p04810','A06p08340'),
  c('A06p08360', 'A06p09570'),
  c('A06p09580','A06p15720'),
  c('A09p68870', 'A09p69380'),
  c('A09p68840','A09p59330'),
  c('A08p37750','A08p29780'),
  c('C05p00980','C05p00030'),
  c('C05p01000', 'C05p09890'),
  c('C05p09900','C05p11280'),
  c('C05p11300','C05p18040'),
  c('C03p00040','C03p00150'),
  c('C08p00100','C08p03620'),
  c('C08p11420','C08p19510'),
  c('C08p53410','C08p42590'),
  c('C08p53440','C08p54180'),
  ### B ###
  c('A06p15760','A06p17950'),
  c('A07p11990','A07p09380'),
  c('A07p12000','A07p13400'),
  c('A07p16250','A07p13410'),
  c('A05p25580','A05p23160'),
  c('A08p29770','A08p27850'),
  c('A08p26440','A08p27830'),
  c('A08p26430', 'A08p24410'),
  c('A08p11090','A08p10380'),
  c('A08p09440','A08p07630'),
  c('A09p42830','A09p40680'),
  c('A09p40670','A09p38660'),
  c('A09p38650','A09p31570'),
  c('C03p73640','C03p72280'),
  c('C03p73650','C03p76800'),
  c('C05p18070','C05p23380'),
  c('C05p25980','C05p23390'),
  c('C05p25990','C05p32560'),
  c('C05p39690','C05p36610'),
  c('C06p12880','C06p11060'),
  c('C07p23870','C07p20040'),
  c('C07p17930','C07p20030'),
  c('C07p17920','C07p14470'),
  c('C08p19520','C08p22370'),
  c('C08p22870','C08p22410')
  ### C ###
  
)
transitions = which(diff(as.numeric(as.factor(orthologs[,"consensus_ACK"])))!=0)


for(block in blocks){
  ortho_boundaries = which(orthologs[,"Bna"] %in% paste(block, '.1_BnaDAR', sep = ''))
  ortho_block = orthologs[min(ortho_boundaries):max(ortho_boundaries),]
  
  bna_range = which(bna_genes[,"ID"] %in% paste(block, '.1_BnaDAR', sep = ''))
  bna_genes_in_range = bna_genes[min(bna_range):max(bna_range),"ID"]
  length(bna_genes_in_range)
  
  ath_range = which(ath_order %in% orthologs[ortho_boundaries,"Ath_agi"])
  ath_genes_in_range = ath_order[min(ath_range):max(ath_range)]
  length(ath_genes_in_range)
  
  intersection = sum(unique(ortho_block[,"Ath_agi"]) %in% ath_genes_in_range)
  disjunction_bna = sum(!bna_genes_in_range %in% ortho_block[,"Bna"])
  disjunction_ath = sum(!ath_genes_in_range %in% ortho_block[,"Ath_agi"])
  total = sum(intersection, disjunction_ath, disjunction_bna)
  jaccard = intersection / total
  
  cat(
    names(sort(table(ortho_block[,"consensus_ACK"]), decreasing = T))[1],
    block[1], block[2],
    bna_range, ath_range,
    length(bna_genes_in_range),
    length(ath_genes_in_range),
    intersection,
    total,
    round(jaccard, 3), '\n'
  )
}


ack_bna = rep(NA,nrow(bna_genes))
names(ack_bna) = bna_genes[,"ID"]
descending = rep(0,nrow(bna_genes))
names(descending) = bna_genes[,"ID"]

ack_bna[bna_genes[,"ID"] %in% rownames(syntelogs)] = syntelogs[bna_genes[bna_genes[,"ID"] %in% rownames(syntelogs),'ID'],"consensus_ACK"]

picked = bna_genes[bna_genes[,"ID"] %in% rownames(syntelogs),'ID']
inc = syntelogs[picked,"Ath_increment"]
descending[picked][inc < 0 & inc > -20] = -1
descending[picked][inc > 0 & inc < 20] = 1


bna_color = ack_colors[ack_bna]
bna_color[is.na(bna_color)] = '#000000'
bna_rgb = apply(col2rgb(as.character(bna_color)),2, function(i)paste(i,collapse=','))
names(bna_rgb) = names(ack_bna)

descending_color = c('220,0,0','200,200,200','0,220,220')[2+descending]
names(descending_color) = names(descending)

descending_names = bna_genes[,'ID']
descending_names[descending_names %in% rownames(syntelogs)] = syntelogs[descending_names[descending_names %in% rownames(syntelogs)], "Ath_agi"]
descending_names = gsub('(.*)\\.1_BnaDAR','Bna_\\1',descending_names)

bna_bed = cbind(bna_genes[,c(1,4,5)],paste(gsub('\\.1_BnaDAR','',bna_genes[,10]),ack_bna, sep='.'),bna_genes[,c(6,7,4,5)], bna_rgb[bna_genes[,"ID"]], rep(1,nrow(bna_genes)), bna_genes[,5]-bna_genes[,4], rep(0,nrow(bna_genes)))
descending_bed = cbind(bna_genes[,c(1,4,5)],descending_names,bna_genes[,c(6,7,4,5)], descending_color[bna_genes[,"ID"]], rep(1,nrow(bna_genes)), bna_genes[,5]-bna_genes[,4], rep(0,nrow(bna_genes)))

write.table(bna_bed,'bna_ACK.bed', quote = F, row.names = F, col.names = F, sep='\t')
write.table(descending_bed,'bna_ascending.bed', quote = F, row.names = F, col.names = F, sep='\t')


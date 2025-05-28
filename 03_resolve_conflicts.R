setwd("C:/Users/schon035/Documents/Projects/taxonomy/ath_bna/bna_pangenome/")


forward = read.table('diamond_forward_darmor-v10.tsv', header = F)
hub_order = unique(forward[,1])
hub_position = 1:length(hub_order)
names(hub_position) = hub_order

reverse = read.table('diamond_reverse_darmor-v10.tsv', header = F)
bzh_order = unique(reverse[,1])
forward[,'hubchrom'] = gsub('scaffold.*','scaffold',gsub('HUBna','',gsub('[TG][0-9]+$','',forward[,1])))
forward[,'bzhchrom'] = gsub('p.*_BnaDAR','',gsub('.*scaffold.*','scaffold',forward[,2]))
reverse[,'hubchrom'] = gsub('scaffold.*','scaffold',gsub('HUBna','',gsub('[TG][0-9]+$','',reverse[,2])))
reverse[,'bzhchrom'] = gsub('p.*_BnaDAR','',gsub('.*scaffold.*','scaffold',reverse[,1]))

forward_samechrom = forward[forward[,"hubchrom"] == forward[,"bzhchrom"] | forward[,"bzhchrom"] %in% c('-','scaffold'),]
reverse_samechrom = reverse[reverse[,"hubchrom"] == reverse[,"bzhchrom"] | reverse[,"hubchrom"] %in% c('-','scaffold'),]

best_forward = forward_samechrom[!duplicated(forward_samechrom[,1]),]
nochrhit_forward = forward[!forward[,1] %in% best_forward[,1],]
best_forward = rbind(best_forward, nochrhit_forward[!duplicated(nochrhit_forward[,1]),])
rownames(best_forward) = best_forward[,1]
best_forward = best_forward[hub_order,]

best_reverse = reverse_samechrom[!duplicated(reverse_samechrom[,1]),]
nochrhit_reverse = reverse[!reverse[,1] %in% best_reverse[,1],]
best_reverse = rbind(best_reverse, nochrhit_reverse[!duplicated(nochrhit_reverse[,1]),])
rownames(best_reverse) = best_reverse[,1]
best_reverse = best_reverse[bzh_order,]


write.table(best_reverse[,c(1,2)],'darmor_to_hub.tsv', row.names = F, col.names = F, sep='\t', quote = F)

bbh = intersect(paste(best_forward[,2], best_forward[,1], sep='_'),
                paste(best_reverse[,1], best_reverse[,2], sep='_'))

bbhmat = matrix(unlist(strsplit(bbh, '_BnaDAR_')),ncol = 2, byrow = T)
bbhmat[,1] = paste(bbhmat[,1],'_BnaDAR',sep='')
write.table(bbhmat,'bbh_darmor-bzh.tsv', quote = F, sep='\t', row.names = F, col.names = F)


syntelogs = read.table('../Bnapus_syntelogs_20240606.tsv', sep='\t', header = T)
syntelogs = syntelogs[!duplicated(syntelogs[,1]),]
rownames(syntelogs) = syntelogs[,1]
orthoMCL = read.table('hub_to_ath.tsv', sep='\t', header = F)
rownames(orthoMCL) = orthoMCL[,1]

ack = read.table('../../../genome_resources/ACK_blocks_nogaps.tsv', header = T, sep='\t', comment.char = '')
rownames(ack) = ack[,"genes"]
bbh_orthologs = data.frame(
  'darmor-bzh'=best_reverse[,1],
  'hubID'=best_reverse[,2],
  'OrthoFinder'=syntelogs[best_reverse[,1],'Ath_agi'],
  'OrthoMCL'=orthoMCL[best_reverse[,2],2]
)
bbh_orthologs[,'OrthoFinder.ACK'] = ack[bbh_orthologs[,"OrthoFinder"],"block"]
bbh_orthologs[,'OrthoMCL.ACK'] = ack[bbh_orthologs[,"OrthoMCL"],"block"]
bbh_orthologs[,'agree'] = bbh_orthologs[,"OrthoFinder"] == bbh_orthologs[,"OrthoMCL"]

x = bbh_orthologs[!(is.na(bbh_orthologs[,"OrthoFinder"]) | is.na(bbh_orthologs[,"OrthoFinder.ACK"])),"OrthoFinder.ACK"]
mean(x[2:length(x)] == x[1:(length(x)-1)])
plot(as.integer(as.factor(x)), pch='.')

x = bbh_orthologs[!(is.na(bbh_orthologs[,"OrthoMCL"]) | is.na(bbh_orthologs[,"OrthoMCL.ACK"])),"OrthoMCL.ACK"]
mean(x[2:length(x)] == x[1:(length(x)-1)])

rownames(bbh_orthologs) = bbh_orthologs[,"darmor.bzh"]

write.table(bbh_orthologs, 'bzh_orthologs_full.tsv', sep='\t', quote = F, row.names = F)

ath_bed = read.table('../../../genome_resources/TAIR_ACK.bed')
ath_order = unique(gsub('\\..*','',ath_bed[ath_bed[,5]=='protein_coding',4]))
ath_position = 1:length(ath_order)
names(ath_position) = ath_order




mean(bbh_orthologs[,'agree'], na.rm = T)

disagree = which(bbh_orthologs[,"agree"] == F)

bbh_orthologs[disagree,"darmor.bzh"]

syntenic_OF = sapply(disagree, function(i){
  consensus = names(sort(table(bbh_orthologs[c(max((i-5),1):max((i-1),1),(i+1):(i+5)),"OrthoFinder.ACK"]), decreasing = T)[1])
  if(length(consensus) == 0){
    return(NA)
  }else{
    return(bbh_orthologs[i,"OrthoFinder.ACK"] == consensus)
  }
  }
)
names(syntenic_OF) = bbh_orthologs[disagree,"darmor.bzh"]


syntenic_OM = sapply(disagree, function(i){
  consensus = names(sort(table(bbh_orthologs[c(max((i-5),1):max((i-1),1),(i+1):(i+5)),"OrthoMCL.ACK"]), decreasing = T)[1])
  if(length(consensus) == 0){
    return(NA)
  }else{
    return(bbh_orthologs[i,"OrthoMCL.ACK"] == consensus)
  }
}
)
names(syntenic_OM) = bbh_orthologs[disagree,"darmor.bzh"]


is_syntenic <- function(i, n=5, row="OrthoFinder.ACK"){
  consensus = names(sort(table(bbh_orthologs[c(max((i-5),1):max((i-1),1),min((i+1),nrow(bbh_orthologs)):min(nrow(bbh_orthologs),(i+5))),"OrthoMCL.ACK"]), decreasing = T)[1])
}


## bestmatch Decision Tree ##

bbh_orthologs[,'bestmatch'] = NA
agree = which(bbh_orthologs[,"agree"])
bbh_orthologs[agree,"bestmatch"] = bbh_orthologs[agree,"OrthoFinder"]

onlyOF = !is.na(bbh_orthologs[,"OrthoFinder"]) & is.na(bbh_orthologs[,"OrthoMCL"])
bbh_orthologs[onlyOF,"bestmatch"] = bbh_orthologs[onlyOF,"OrthoFinder"]

onlyOM = is.na(bbh_orthologs[,"OrthoFinder"]) & !is.na(bbh_orthologs[,"OrthoMCL"])
bbh_orthologs[onlyOM,"bestmatch"] = bbh_orthologs[onlyOM,"OrthoMCL"]

sOF = names(which(!syntenic_OM & syntenic_OF))
bbh_orthologs[sOF,"bestmatch"] = bbh_orthologs[sOF,"OrthoFinder"]

sOM = names(which(syntenic_OM & !syntenic_OF))
bbh_orthologs[sOM,"bestmatch"] = bbh_orthologs[sOM,"OrthoMCL"]

rest = is.na(bbh_orthologs[,"bestmatch"])
bbh_orthologs[rest,"bestmatch"] = bbh_orthologs[rest,"OrthoFinder"]

write.table(bbh_orthologs, 'bzh_orthologs_full.20250204.tsv', sep='\t', quote = F, row.names = F)

ath_copies = sapply(ath_order, function(x)bbh_orthologs[which(bbh_orthologs[,"bestmatch"]==x),"darmor.bzh"])
ath_copytable = data.frame(
  'ath'=names(ath_copies),
  'bna_copies'=unlist(lapply(ath_copies, length)),
  'name'=geneName(names(ath_copies)),
  'bna_IDs'=unlist(lapply(ath_copies, function(i)paste(i,collapse = ',')))
)
write.table(ath_copytable, 'ath_copies_bna.tsv', sep='\t', quote=F, row.names = F)



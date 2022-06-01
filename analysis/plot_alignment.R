# Plot MSA output with domain annotated
library(tidyverse)
library(msa)
library(ggmsa)
library(RColorBrewer)

fp = '../tutorial/RASSF/aligned.fasta'       # path to DAMPSA output alignment.
tb_path = '../tutorial/RASSF/domain.txt'     # path to DAMPSA domain detection output

output = file.path(dirname(fp), 'vis',
                   paste(strsplit(basename(fp), '\\.')[[1]][1], '_vis.pdf', sep=''))

height = 8                  # output plot height
width = 10                  # output plot width
no_label_y = T              # label y (sequence name) or not.
no_legend = F               # include domain legend or not
master = c()                # master sequence ID, will plot on the top (y), if c() then not consider.
top = c()                   # sequences at the top (y), if c() then not consider.
yorder = c()                # force sequence order on the y axis.
no_id_match = F             # biopython usually split ID, so ID matching is necessary. In most case should set to F.

# For RASSF family order in the report
# yorder = c("sp|Q9NS23|RASF1_HUMAN/1-344 Ras association domain-containing protein 1 OS=Homo sapiens OX=9606 GN=RASSF1 PE=1 SV=1",  
#            "sp|P50749|RASF2_HUMAN/1-326 Ras association domain-containing protein 2 OS=Homo sapiens OX=9606 GN=RASSF2 PE=1 SV=1",  
#            "sp|Q86WH2|RASF3_HUMAN/1-238 Ras association domain-containing protein 3 OS=Homo sapiens OX=9606 GN=RASSF3 PE=1 SV=1",  
#            "sp|Q9H2L5|RASF4_HUMAN/1-321 Ras association domain-containing protein 4 OS=Homo sapiens OX=9606 GN=RASSF4 PE=1 SV=2",  
#            "sp|Q8WWW0|RASF5_HUMAN/1-418 Ras association domain-containing protein 5 OS=Homo sapiens OX=9606 GN=RASSF5 PE=1 SV=1",  
#            "sp|Q6ZTQ3|RASF6_HUMAN/1-369 Ras association domain-containing protein 6 OS=Homo sapiens OX=9606 GN=RASSF6 PE=1 SV=1",  
#            "sp|Q02833|RASF7_HUMAN/1-373 Ras association domain-containing protein 7 OS=Homo sapiens OX=9606 GN=RASSF7 PE=1 SV=1",  
#            "sp|Q8NHQ8|RASF8_HUMAN/1-419 Ras association domain-containing protein 8 OS=Homo sapiens OX=9606 GN=RASSF8 PE=1 SV=2", 
#            "sp|O75901|RASF9_HUMAN/1-435 Ras association domain-containing protein 9 OS=Homo sapiens OX=9606 GN=RASSF9 PE=2 SV=2", 
#            "sp|A6NK89|RASFA_HUMAN/1-507 Ras association domain-containing protein 10 OS=Homo sapiens OX=9606 GN=RASSF10 PE=1 SV=3")
yorder = yorder[length(yorder):1]

## read alignment ====
aln = readAAStringSet(fp)
aln_long = tidy_msa(aln)
aln = as.matrix(aln)
colnames(aln) = 1:ncol(aln)

## match seq ID ====
dom_tb = read.table(tb_path, sep='\t', header = T)
ids = unique(dom_tb$seq_id)
if (no_id_match){
  print('Not matching ID...')
} else {
  id_dic = as.character(unique(rownames(aln)))
  names(id_dic) = id_dic
  for (i in 1:length(id_dic)){
    for (j in ids){
      if (length(grep(gsub('\\|', '=', j), 
                      gsub('\\|', '=', names(id_dic)[i]))) > 0){
        names(id_dic)[i] = j
        break
      }
    }
  }
  dom_tb$seq_id = id_dic[dom_tb$seq_id]
}

## convert sequence coordinate to alignment corrdinate ====
seq2aln_loc = function(id, aln, x){
  # id: sequence ID (vector), aln: matrix, x: location(s) to convert
  rt = c()
  if (length(id) == 1 & length(x) != 1){
    id = rep(id, length(x))
  }
  for (j in 1:length(x)){
    i = x[j]
    while (i <= dim(aln)[2]){
      roi = aln[id[j],1:i]
      if (length(roi[roi!='-']) == x[j]){
        rt = c(rt, i)
        break
      }
      i = i + 1
    }
  }
  return(rt)
}
dom_tb$env_start_aln = NA
dom_tb$env_end_aln = NA

dom_tb$env_start_aln = seq2aln_loc(dom_tb$seq_id, 
                                   aln, dom_tb$env_start)
dom_tb$env_end_aln = seq2aln_loc(dom_tb$seq_id, 
                                 aln, dom_tb$env_end)

# generate domain color matrix ====
aln_col = matrix(NA, nrow=dim(aln)[1], ncol=dim(aln)[2])
rownames(aln_col) = rownames(aln)
colnames(aln_col) = colnames(aln)
for (i in 1:nrow(dom_tb)){
  aln_col[dom_tb$seq_id[i], 
          dom_tb$env_start_aln[i]:dom_tb$env_end_aln[i]] = dom_tb$clan[i]
}
aln_col = data.frame(aln_col)
colnames(aln_col) = colnames(aln)
aln_col['name'] = rownames(aln_col)
aln_col = gather(aln_col, key='position', value='col', 1:(ncol(aln_col)-1))
aln_col$position = as.double(aln_col$position)

aln_long = inner_join(aln_long, aln_col, by=c('name','position'))

aln_long$col[aln_long$character!='-' & is.na(aln_long$col)] = 'zlinker'
if (length(master) > 0){
  msub = aln_long[aln_long$name==master,]
  offset = sort(msub[msub$character!='-','position'])[1]
  aln_long$position = aln_long$position - offset
}
if (length(top) > 0){
  lvl = c(top, sort(setdiff(unique(aln_long$name), top)))
  aln_long$name = factor(aln_long$name, levels = lvl[length(lvl):1])
}
if (length(yorder) > 0){
  aln_long$name = factor(aln_long$name, levels = yorder)
}

n_clan = length(unique(aln_long$col))-2
g = ggplot(aln_long) +
  geom_tile(aes(x=position, y=name, fill=as.factor(col), color=as.factor(col)), 
            height=0.9) +
  scale_fill_manual(values=c(brewer.pal(n_clan,'Set2')[1:n_clan],'black'), 
                    na.value = NA) +
  scale_color_manual(values=c(brewer.pal(n_clan,'Set2')[1:n_clan],'black'), 
                    na.value = NA) +
  labs(x='', y='') +
  theme_classic()
if (no_label_y){
  g = g + theme(axis.line.y = element_blank(),
                axis.ticks = element_blank(), axis.text = element_blank())
}
if (no_legend){
  g = g + theme(legend.position = 'none')
}
g %>% ggsave(output, plot=., width=width,
             height=height,
             units='in',device='pdf',dpi=300, limitsize = F, bg = 'transparent')


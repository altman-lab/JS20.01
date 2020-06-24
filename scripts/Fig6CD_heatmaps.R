library(tidyverse)
library(ComplexHeatmap)
library(pvclust)
#Magma color scheme
library(viridis)
set.seed(546)

#### Data ####
counts <- read_csv("results/module_Shah_contrast_deepSplit3_minMod50/Shah_contrast_mod_voom_counts.csv") %>% 
  #Simplify module names
  mutate(module = gsub("module_Shah_contrast_","",module)) 

meta <- read_csv("data_clean/Shah.metadata.csv") %>% 
  mutate(cell = factor(cell, levels=c("WT","TKO")),
         status = factor(status, levels=c("Uninfected","Infected"))) %>% 
  #Sort and format to matrix
  select(sampID, status, cell) %>% 
  arrange(status, cell) %>% 
  column_to_rownames("sampID")

#### Calculate Z scores #####
counts.long <- counts %>% 
  #Remove module 0 
  filter(module != "00") %>% 
  #Add metadata
  pivot_longer(-module, names_to = "rowname") %>% 
  left_join(rownames_to_column(meta)) %>% 
  select(-rowname)

#Mean and SD for each module
mean.mod <- counts.long %>% 
  group_by(module) %>% 
  summarize(mod.mean = mean(value, na.rm=TRUE),
            mod.sd = sd(value, na.rm=TRUE)) %>% 
  ungroup() 

#Calculate Z-score
z.score <- full_join(counts.long, mean.mod,
                     by = c("module")) %>% 
  group_by(module, status, cell) %>% 
  mutate(z = (value - mod.mean) / mod.sd) %>% 
  ungroup() %>% 
  
  #Recode treatment groups
  mutate(group = paste(status,cell, sep="_")) %>% 
  
  #Mean z score for groups
  group_by(module, group) %>% 
  summarize(z.mean = mean(z, na.rm=TRUE)) %>% 
  #Wide format
  arrange(group) %>% 
  select(module, group, z.mean) %>% 
  pivot_wider(names_from = group, values_from = z.mean) %>% 
  #Remove leading 0 in module name
  mutate(module = sub("^0+", "", module)) %>% 
  #format to matrix
  column_to_rownames("module") %>% 
  as.matrix()
  
#### Format meta #####
meta.sub<- meta %>% 
  #Recode treatment groups
  mutate(group = paste(status,cell, sep="_")) %>% 
  arrange(group) %>% 
  distinct() %>% 
  column_to_rownames("group") %>% 
  as.matrix()

#### FIGURE 6C: MODULES ####
#### Tree ####
corr.pv <- pvclust(t(z.score), nboot=1000, 
                   method.hclust="average", method.dist="correlation")

#### column (module) annotation ####
col_annot.df <- data.frame(module = counts$module) %>% 
  mutate(Uninfected = 
           ifelse(module %in% c("14","07","16","09","13","15"),
                  "Up",
                  ifelse(module %in% c("05","11","17","10","02","06","12"),
                         "Down",
                         NA)),
         Infected = 
           ifelse(module %in% c("14","07","16","09","01","03","00"),
                  "Up",
                  ifelse(module %in% c("05","11","17","10","02","04","08"),
                         "Down",
                         NA))) %>% 
  #Remove module 0 
  filter(module != "00") %>% 
  #Add total genes in module
  left_join(count(read_csv(
    "results/module_Shah_contrast_deepSplit3_minMod50/Shah_contrast_genes_in_mod.csv"),
    module.char) , by=c("module"="module.char"))

col_annot <- HeatmapAnnotation(
  Uninfected =  anno_simple(col_annot.df$Uninfected,
                            col = c("Up"="#ca0020", "Down"="#0571b0"), 
                            na_col = "white"),
  Infected =  anno_simple(col_annot.df$Infected,
                          col = c("Up"="#ca0020", "Down"="#0571b0"), 
                          na_col = "white"))

col_legend <- Legend(title = "Significant\nfold change",
                     legend_gp = gpar(fill = c("#ca0020", "#0571b0")),
                     labels = c("Up","Down"))

#### heatmap ####
rowNames <- c(expression(Uninfected~italic(Tollip^"-/-")),
              expression(Uninfected~italic(Tollip^"+/+")),
              expression(Infected~italic(Tollip^"-/-")),
              expression(Infected~italic(Tollip^"+/+")))

mod_hm <- Heatmap(t(z.score), name = "Mean\nZ-score",
                  #Expression colors
                  col = magma(20),
                  #Sample annot
                  row_order = c("Uninfected_WT",
                                "Uninfected_TKO",
                                "Infected_WT",
                                "Infected_TKO"),
                  row_split = c(1,1,2,2),
                  row_gap = unit(5, "mm"),
                  show_row_dend = FALSE,
                  row_title = " ",
                  row_labels = rowNames,
                  #Module annot
                  column_names_side = "top",
                  top_annotation = col_annot,
                  cluster_columns = corr.pv$hclust,
                  column_split = 2, column_gap = unit(5, "mm"),
                  column_dend_height = unit(2, "cm"),
                  column_names_rot = 0,
                  column_title = " ",
                  column_names_centered = TRUE,
                  #Force square
                  heatmap_height = unit(10, "cm"),
                  heatmap_width = unit(24, "cm"))

#### Save ####
pdf(file = "figs/publication/heatmap_modules.4groups.Zscore.pdf", 
    height=6, width=12)

draw(mod_hm, annotation_legend_list=list(col_legend))
dev.off()

##### FIGURE 6D: HALLMARK TERMS #####
#### Data ####
hallmark <- read_csv("results/GSEA/GSEA_modules_H.csv")

FDR.cutoff <- .05

#### Format data: FDR####
#
hallmark_fdr <- hallmark %>% 
  #remove mod 0
  filter(group != "00") %>% 
  #Calculate proportion of genes in term
  select(group, Description, size.overlap.term, 
         size.group, p.adjust) %>% 
  mutate(fdr.scale = -log10(p.adjust)) %>% 
  #Format labels
  mutate(Description = gsub("HALLMARK_","",Description),
         Description = gsub("_", " ", Description) )

#list terms with at least 1 module FDR < 0.5
hallmark_summ <- hallmark_fdr %>% 
  group_by(Description) %>% 
  summarize(fdr.min = min(p.adjust)) %>% 
  arrange(-fdr.min)

terms.to.keep <- hallmark_summ %>% 
  filter(fdr.min<=FDR.cutoff) %>% 
  select(Description) %>% unlist(use.names = FALSE)

hallmark_sub <- hallmark_fdr %>% 
  filter(Description %in% terms.to.keep) %>% 
  #Wide format
  select(group, Description, fdr.scale) %>% 
  pivot_wider(names_from = Description, values_from = fdr.scale) %>% 
  arrange(group) %>% 
  #Fill NAs
  mutate_if(is.numeric, ~ifelse(is.na(.),0,.)) %>% 
  #To matrix
  column_to_rownames("group") %>% 
  as.matrix()


#### Format data ####
#Calculate percent of genes in term
hallmark_pct <- hallmark %>% 
  #remove mod 0
  filter(group != "00") %>% 
  #Calculate proportion of genes in term
  select(group, Description, size.overlap.term, 
         size.group, p.adjust) %>% 
  mutate(pct = size.overlap.term/size.group*100) %>% 
  #Format labels
  mutate(Description = gsub("HALLMARK_","",Description),
         Description = gsub("_", " ", Description) )

#list terms with at least 1 module FDR < 0.5
hallmark_summ <- hallmark_pct %>% 
  group_by(Description) %>% 
  summarize(pct.max=max(pct), fdr.min = min(p.adjust)) %>% 
  arrange(-fdr.min)

terms.to.keep <- hallmark_summ %>% 
  filter(fdr.min<=FDR.cutoff) %>% 
  select(Description) %>% unlist(use.names = FALSE)

hallmark_sub <- hallmark_pct %>% 
  filter(Description %in% terms.to.keep) %>% 
  #Wide format
  select(group, Description, pct) %>% 
  pivot_wider(names_from = Description, values_from = pct) %>% 
  arrange(group) %>% 
  #Fill NAs
  mutate_if(is.numeric, ~ifelse(is.na(.),0,.)) %>% 
  #To matrix
  column_to_rownames("group") %>% 
  as.matrix()

#### Column anno ####
col_annot <- HeatmapAnnotation(
  Uninfected =  anno_simple(col_annot.df$Uninfected,
                            col = c("Up"="#ca0020", "Down"="#0571b0"), 
                            na_col = "white"),
  Infected =  anno_simple(col_annot.df$Infected,
                          col = c("Up"="#ca0020", "Down"="#0571b0"), 
                          na_col = "white"),
  "Total genes" = anno_simple(col_annot.df$n,
                              col=circlize::colorRamp2(c(0, 600),
                                                       c("white", "darkgrey"))),
  "Total genes2" = anno_text(col_annot.df$n, rot=0, just="center",
                             location = unit(2, "npc")))

#### heatmap ####
hallmark_hm <- Heatmap(t(hallmark_sub), name = "-log10(FDR)",
                       #Expression colors
                       col = magma(20),
                       #Module annot
                       column_names_side = "top",
                       top_annotation = col_annot,
                       cluster_columns = corr.pv$hclust,
                       column_split = 2, column_gap = unit(5, "mm"),
                       column_dend_height = unit(2, "cm"),
                       column_names_rot = 0,
                       column_title = " ",
                       column_names_centered = TRUE,
                       #Rows
                       clustering_method_rows = "average",
                       #Force square
                       heatmap_height = unit(22, "cm"),
                       heatmap_width = unit(24, "cm"))

#### save
pdf(file = paste("figs/publication/heatmap_hallmark_FDR",
                 FDR.cutoff, "B.pdf", sep=""), 
    height=9, width=12)

draw(hallmark_hm, annotation_legend_list=list(col_legend))

dev.off()


library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(reshape2)
# install.packages("viridis")
library(viridis)
# install.packages('gghalves')
# library(gghalves)

# geom_split_violin()
# install.packages("devtools")
# devtools::install_github("iholzleitner/facefuns")
library(facefuns)

# ==== 1.1 S3: plot post filtered FSM and ISM isoform read length

# prepare read length data
if (T) {
  UMI_length <- vector("list",8)
  names(UMI_length) <- c('revio_sample1B','revio_sample1Bjumpcode','revio_sample2A','revio_sample2Ajumpcode',
                         'sequel_sample1B','sequel_sample1Bjumpcode','sequel_sample2A','sequel_sample2Ajumpcode')
  for(sample in names(UMI_length)){
    # sample='revio_sample1B'
    fai <- file.path('scratch_space/CB-MAS-Seq_figure_yl/Analysis/cromwell_data/',
                     sample,'filter_UMI_full-splice_match.fasta.fai')
    df <- read.table(fai, header = F, check.names = F, stringsAsFactors = F, sep = '\t')
    UMI_length[[sample]][['FSM']] <- data.frame(length=df$V2, type=rep("FSM",nrow(df)))
    
    fai <- file.path('scratch_space/CB-MAS-Seq_figure_yl/Analysis/cromwell_data/',
                     sample,'filter_UMI_incomplete-splice_match.fasta.fai')
    df <- read.table(fai, header = F, check.names = F, stringsAsFactors = F, sep = '\t')
    UMI_length[[sample]][['ISM']] <- data.frame(length=df$V2, type=rep("ISM",nrow(df)))
  }
  save(UMI_length, 
       file='data/UMI_length.RData')
}
# plot
if (T) {
  load('data/UMI_length.RData', verbose = T)
  
  # load read length of ONT
  sample='ONT_sample1B'
  df <- read.table('data/Earlham1PBMC_Iain.merged.sorted.bam.readlen.2.read_tags.FSM.dedup.tsv', 
                   header = F, check.names = F, stringsAsFactors = F, sep = '\t')
  UMI_length[[sample]][['FSM']] <- data.frame(length=df$V3, type=rep("FSM",nrow(df)))
  df <- read.table('data/Earlham1PBMC_Iain.merged.sorted.bam.readlen.2.read_tags.ISM.dedup.tsv', 
                   header = F, check.names = F, stringsAsFactors = F, sep = '\t')
  UMI_length[[sample]][['ISM']] <- data.frame(length=df$V3, type=rep("ISM",nrow(df)))

  # plot FSM using density plot
  for (type in c("FSM", "ISM")){
  
    df.plot <- data.frame(matrix(ncol = 4, nrow = 0))
    colnames(df.plot) <- c('length', 'type', 'sample', 'platform')
    for (sample in names(UMI_length)){
      df <- UMI_length[[sample]][[type]]
      df$sample <- sample
      platform <- unlist(str_split(sample, pattern = "_"))[1]
      df$platform <- platform
      df.plot <- rbind(df.plot, df)
    }
    p <- ggplot(data = df.plot, aes(x=length, y=..count.., color=sample))+
      geom_density(adjust=0.1)+xlim(0,2000)+scale_color_viridis_d()+
      theme_bw()+ggtitle(sprintf("post filtered %s isoform read length", type))
    p+facet_wrap(~platform, nrow=3, scales='free')
    # save as Sample1B.FSM_read_length.v6.2kb.pdf
    
    p <- ggplot(data = df.plot, aes(x=length, y=..count.., color=sample))+
      geom_density(adjust=0.1)+xlim(0,5000)+scale_color_viridis_d()+
      theme_bw()+ggtitle(sprintf("post filtered %s isoform read length", type))
    p+facet_wrap(~platform, nrow=3, scales='free')
    # save as Sample1B.FSM_read_length.v6.5kb.pdf
  }
}

# ==== 2. S4 post filtered FSM/ISM isoform length vs transcript reference 
if (T) {

  if (T) {
    # initialise data 
    df.plot <- data.frame(matrix(ncol = 4, nrow = 0))
    # type: pb_isoform, reference
    colnames(df.plot) <- c('length','catagory','type','sample')
    sample="ONT_sample1B"
    for (category in c('FSM', 'ISM')) {
      tsv <- read.table(file.path(sprintf("project_space/data/data/MASseq/nanopore/wf-single-cell-YL/Cribbs.run2/Earlham1BPBMC_Ian/sqanti.YL/%s_isoform.length.csv",
                                          category)), header = F, check.names = F, stringsAsFactors = F, sep = ',')
      # tail(tsv)
      #       iso.  ref.
      #         V2   V3
      # 32175  884 4640
      tsv <- tsv[,2:3]
      colnames(tsv) <- c('isoform', 'ref')
      tsv <- melt(tsv)
      # head(tsv)
      #      variable value
      # 1     isoform  1599

      # replace "length" -> "ISM_pb_isoform"
      # replace "ref_length" -> 'ISM_ref'
      tsv$variable <- paste(category, tsv$variable, sep = '_')
      tsv$sample <- sample
      tsv$catagoryt <- category
      colnames(tsv) <- c('type', 'length', 'sample','catagory')
      # head(tsv)
      #                 type length       sample catagory
      # 1 ONT_isoform_length   1599 ONT_sample1B      FSM
      df.plot <- rbind(df.plot, tsv)
    }  
    
    for(sample in c('revio_sample1B','revio_sample1Bjumpcode','sequel_sample1B','sequel_sample1Bjumpcode')){
      print(sample)
      annot_df <- read.table(file.path('project_space/CB-MAS-Seq_figure_yl/Analysis/cromwell_data',
                                       sample,"scisoseq_classification.filtered_lite_classification.deref.txt"),
                             header = T, check.names = F, stringsAsFactors = F, sep = '\t')
      # modify annot_df category for display
      annot_df_mod <- annot_df %>%
        mutate(structural_category = recode(structural_category, "novel_not_in_catalog" = "NNC", 
                                            "novel_in_catalog" = "NIC", 
                                            "incomplete-splice_match" = "ISM",
                                            "fusion" = "fusion",
                                            "full-splice_match" = "FSM",
                                            "genic" = "genic",
                                            "antisense" = "antisense",
                                            "intergenic" = "intergenic",
                                            "moreJunctions" = "more junctions"))
      
      for (category in c("FSM",'ISM')){
        
        # keep only ISM and FSM
        annot_df_mod.cat <- annot_df_mod[(annot_df_mod$structural_category==category),
                                         c('length','structural_category', "ref_length")]
        # colnames(annot_df_mod.cat) <- c('pb_isoform', 'category', 'reference')
        tmp <- melt(annot_df_mod.cat, 'structural_category')
        tmp$sample=sample
        tmp$variable <- as.character(tmp$variable)
        tail(tmp)
        # structural_category   variable value                  sample
        # 1               ISM   length   391 sequel_sample1Bjumpcode # this is PB isoform
        # 323797          ISM ref_length  2032 sequel_sample1Bjumpcode # this is ref
        
        # replace "length" -> "ISM_pb_isoform"
        # replace "ref_length" -> 'ISM_ref'
        tmp$variable[tmp$variable=='length'] <- 
          paste(tmp$structural_category[tmp$variable=='length'], 'isoform', sep='_')
        tmp$variable[tmp$variable=='ref_length'] <- 
          paste(tmp$structural_category[tmp$variable=='ref_length'], 'ref', sep='_')
        
        colnames(tmp) <- c('catagory', 'type', 'length', 'sample')
        # head(tmp)
        #              catagory        type length                   sample
        # 1                 ISM ISM_isoform    391  sequel_sample1Bjumpcode
        df.plot <- rbind(df.plot, tmp)
      }
    }
    
    save(isoform_length, 
         file='data/isoform_length.RData')
  }
  
  # density plot
  if (T){
    # load(file='data/isoform_length.RData', 
    #      verbose = T)
    head(df.plot)
    #          type length       sample catagory
    # 1 FSM_isoform   1599 ONT_sample1B      FSM
    # 2 FSM_isoform   2179 ONT_sample1B      FSM
    
    # plot density without color fill
    if(T){
      # p <- ggplot(data=df.plot, aes(x=length, y=..count.., group = type, color=type, linetype=type))+
      #   geom_density(adjust=0.1, alpha=0.4)+xlim(0,7500)+
      #   scale_linetype_manual(values=c("FSM_isoform"='solid',     
      #                                  "FSM_ref"='dotted',         
      #                                  "ISM_isoform"='solid',     
      #                                  "ISM_ref"='dotted'))+
      #   scale_color_manual(values=c("FSM_isoform"="#440154" ,    
      #                               "FSM_ref"="#440154",          
      #                               "ISM_isoform"="#31688e", 
      #                               "ISM_ref"="#31688e"))+
      #   theme_bw()
      # 
      # p
      
      p <- ggplot(data=df.plot, aes(x=length, y=..count.., group = type, color=type, linetype=type))+
        geom_density(adjust=0.1, alpha=0.4)+xlim(0,6000)+
        scale_linetype_manual(values=c("FSM_isoform"='solid',     
                                       "FSM_ref"='dashed',         
                                       "ISM_isoform"='solid',     
                                       "ISM_ref"='dashed'))+
        scale_color_manual(values=c("FSM_isoform"="#440154" ,    
                                    "FSM_ref"="#440154",          
                                    "ISM_isoform"="#31688e", 
                                    "ISM_ref"="#31688e"))+
        theme_bw()
      p+facet_wrap(~sample, nrow=5, scales='free')
      
      ggsave(file.path('data',
                       "pb_isoform_vs_ref.length.v4.pdf"), p)
    }
    
  }
  
}

# ==== 3. S5: Distribution of bases for post isoform filtered deduplicated reads that align to reference genome
if (T) {
  df.read_dist <- read.table('data/read_dist.incONT.bps.csv',
                             header = T, check.names = F, stringsAsFactors = F, sep = ',', row.names = 1)
  # df.read_dist<- df.read_dist[c(1,3,5,7,9,10),]
  df.p_pct <- data.frame(pct=c(df.read_dist$PCT_CODING_BASES,df.read_dist$PCT_UTR_BASES,df.read_dist$PCT_INTRONIC_BASES,
                               df.read_dist$PCT_INTERGENIC_BASES),
                         region=c(rep('CDS',nrow(df.read_dist)), rep('UTR',nrow(df.read_dist)),
                                  rep('intronic',nrow(df.read_dist)),rep('intergenic',nrow(df.read_dist))),
                         sample=c(rep(row.names(df.read_dist),4)))
  p <- ggplot(df.p_pct, aes(x=sample, y=pct, fill=region))+
    geom_bar(stat = 'identity')+scale_fill_viridis(discrete=T, option = 'H')+
    theme_bw()+ylab("%bases")+
    ggtitle('Distribution of aligned bases')+
    theme(plot.title = element_text(size=20), axis.title = element_text(size=15), axis.text = element_text(size = 12),
          legend.text = element_text(size=12), legend.title = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))
  p
  ggsave(file.path('data',
                   "read_distribution.pc_bases.sample1B.filtered.v2.pdf"), p, width = 20, height = 15, units = 'cm')
  
}
# plot CDS only transcript body coverage using output from rseqc geneBodyCoverage function
# ==== 1. check rseqc outputs =====
if (T) {
  txtype="FSM"
  len_range="1001-2000"
  # len_range="0-500"
  # len_range="501-1000"
  
  # prepare empty data matrix
  nb_sample=6
  tmp <- rep(1,nb_sample)
  dim(tmp) <- c(nb_sample,1)
  data_matrix <- data.frame(tmp%*%rep(0,100))
  if (nb_sample==6) {
    rownames(data_matrix) <- c("Illumina_sample1B", "Revio_sample1B", 
                               "Revio_sample1Bjumpcode",
                               "Sequel_sample1B",  
                               "Sequel_sample1Bjumpcode", "ONT_sample1B")
  } 
  rowLabel <- rownames(data_matrix)
  
  # load normalised coverage.
  for (sample in rownames(data_matrix)) {
    if (sample=="Illumina_sample1B"){
      df <- read.table(file=file.path(sprintf('data/%s.%s.geneBodyCoverage.txt',
                                              sample, len_range)),
                       header = T, check.names = F, sep = '\t')
      df <- df[2:101]
      # z-score norm
      df <- (df-min(df))/(max(df)-min(df))
      
    }else {
    df <- read.table(file=file.path(sprintf('data/%s.%s.range.geneBodyCoverage',
                                            sample, len_range)),
                            header = F, check.names = F, sep = '\t')
    df <- df[1:100]
    }
    
    data_matrix[sample,] <- df
  }
  rownames(data_matrix) <- rowLabel
}
# ==== 3. plot geneBodyCoverage using ggplt2: new color combination ====
if (T) {
  # install.packages("viridis") 
  library(viridis)
  library(ggplot2)
  library(stringr)
  
  combined_legend <- c()
  rowLabel <- rownames(data_matrix)
  for (sample in rowLabel) {
    combined_legend <- c(combined_legend, rep(sample, 100))
  }
  if (nb_sample==6){
    p5 <- ggplot(data=data.frame(coverage=c(t(data_matrix)), 
                                 percentile=rep(x=1:100,times=length(rowLabel)), 
                                 sample=combined_legend,
                                 #sample=sample, 
                                 technology=c(rep('Illumina',100), rep('Revio',200), 
                                              rep('Sequel', 200), rep('ONT',100)), 
                                 replicates=c(rep('1B',100), rep('1B',100), rep('1Bjumpcode',100), 
                                              rep('1B',100), rep('1Bjumpcode',100), 
                                              rep('ONT',100))
      ), 
    aes(x=percentile, y=coverage))+
      geom_line(aes(color=sample, linetype=sample), linewidth=1)+
      scale_linetype_manual(values=c("Illumina_sample1B"='dotted',     
                                     "Revio_sample1B"='dashed',         
                                     "Revio_sample1Bjumpcode"='dashed', 
                                     "Sequel_sample1B"='solid', 
                                     "Sequel_sample1Bjumpcode"='solid', "ONT_sample1B"='solid'))+
      scale_color_manual(values=c("Illumina_sample1B"="#009E73" ,    
                                  "Revio_sample1B"="#E69F00",          
                                  "Revio_sample1Bjumpcode"="#56B4E9", 
                                  "Sequel_sample1B"="#E69F00", 
                                  "Sequel_sample1Bjumpcode"="#56B4E9", "ONT_sample1B"="#000000"))+
      theme_bw()+xlab("Gene body percentile (5' -> 3')")+
      ggtitle(sprintf('Gene body coverage for transcript length %s bps',cur_len_range))+
      theme(plot.title = element_text(size=20), axis.title = element_text(size=15), axis.text = element_text(size = 12),
            legend.text = element_text(size=12), legend.title = element_text(size=15))
    # guides(linetype=guide_legend(override.aes = list(color=c(rep('red',3)))), 
    #        color=guide_legend(override.aes = list(linetype=c(rep('dotted',2), rep('dashed',4), rep('solid',4)))))
    
    ggsave(sprintf("nf_6samples.FSM.len%s.geneBodyCoverage.v4.pdf",cur_len_range), p5, width = 30, height = 25, units = 'cm')
    
  } 
  
}




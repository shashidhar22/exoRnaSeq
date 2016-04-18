This documents lists the pipeline used for the assembly of 8 mouse macrophage RNAseq samples:

1. Sample information:
  The study contained RNASeq samples from three conditions:
  1. control exosome: Contained three biological replicates with each replicate having two technical replicates.
  2. pIC exosome: Contained three biological replicates with each replicate having two technical replicates.
  3. PBS : Contained two biological replicates with each replicate having two technical replicates. However one of the biological replicate showed contamination in capture, and was left out of the analysis/
2. Reference information:
  The mouse mm10 reference genome and annotation files (UCSC) were downloaded from the Illumina igenomes, genome repository.
  [Mus musculus (mm10)](ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz)
3. TopHat:
   Tophat was run using the following commands:

  1. Control M1:
    ```{sh}
  tophat2 --no-coverage-search --b2-very-sensitive --no-novel-juncs -p 4 -G ~/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/Control_Exo_M1_MC1-31664696/Control_ExoM1_UCSC ~/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome ~/Control_Exo_M1_MC1-31664696/Control-Exo-M1-MC1_S5_L001_R1_001.fastq.gz,~/Control_Exo_M1_MC1-31664696/Control-Exo-M1-MC1_S5_L002_R1_001.fastq.gz ~/Control_Exo_M1_MC1-31664696/Control-Exo-M1-MC1_S5_L001_R2_001.fastq.gz,~/Control_Exo_M1_MC1-31664696/Control-Exo-M1-MC1_S5_L002_R2_001.fastq.gz
    ```

  2. Control M2:
    ```{sh}
  tophat2 --no-coverage-search --b2-very-sensitive --no-novel-juncs -p 4 -G ~/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/Control_Exo_M2_MC1-31646815/Control_ExoM2_UCSC ~/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome ~/Control_Exo_M2_MC1-31646815/Control-Exo-M2-MC1_S6_L001_R1_001.fastq.gz,~/Control_Exo_M2_MC1-31646815/Control-Exo-M2-MC1_S6_L002_R1_001.fastq.gz ~/Control_Exo_M2_MC1-31646815/Control-Exo-M2-MC1_S6_L001_R2_001.fastq.gz,~/Control_Exo_M2_MC1-31646815/Control-Exo-M2-MC1_S6_L002_R2_001.fastq.gz
  ```

  3. Control M3:
    ```{sh}
  tophat2 --no-coverage-search --b2-very-sensitive --no-novel-juncs -p 4 -G ~/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/Control_Exo_M3_MC1/Control_ExoM3_UCSC ~/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome ~/Control_Exo_M3_MC1/Control_Exo_M3_MC_S20_L001_R1_001.fastq.gz,~/Control_Exo_M3_MC1/Control_Exo_M3_MC_S20_L002_R1_001.fastq.gz ~/Control_Exo_M3_MC1/Control_Exo_M3_MC_S20_L001_R2_001.fastq.gz,~/Control_Exo_M3_MC1/Control_Exo_M3_MC_S20_L002_R2_001.fastq.gz
  ```

  4. pIC M1:
    ```{sh}
  tophat2 --no-coverage-search --b2-very-sensitive --no-novel-juncs -p 4 -G ~/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/pIC_exo_M1_MC1-31656710/pIC_exoM1_UCSC ~/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome ~/pIC_exo_M1_MC1-31656710/pIC-exo-M1-MC1_S7_L001_R1_001.fastq.gz,~/pIC_exo_M1_MC1-31656710/pIC-exo-M1-MC1_S7_L002_R1_001.fastq.gz ~/pIC_exo_M1_MC1-31656710/pIC-exo-M1-MC1_S7_L001_R2_001.fastq.gz,~/pIC_exo_M1_MC1-31656710/pIC-exo-M1-MC1_S7_L002_R2_001.fastq.gz
  ```

  5. pIC M2:
    ```{sh}
  tophat2 --no-coverage-search --b2-very-sensitive --no-novel-juncs -p 4 -G ~/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/pIC_exo_M2_MC1-31646822/pIC_exoM2_UCSC ~/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome ~/pIC_exo_M2_MC1-31646822/pIC-exo-M2-MC1_S8_L001_R1_001.fastq.gz,~/pIC_exo_M2_MC1-31646822/pIC-exo-M2-MC1_S8_L002_R1_001.fastq.gz ~/pIC_exo_M2_MC1-31646822/pIC-exo-M2-MC1_S8_L001_R2_001.fastq.gz,~/pIC_exo_M2_MC1-31646822/pIC-exo-M2-MC1_S8_L002_R2_001.fastq.gz
  ```

  6. pIC M3:
    ```{sh}
  tophat2 --no-coverage-search --b2-very-sensitive --no-novel-juncs -p 4 -G ~/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/pIC_exo_M3_MC1/pIC_exoM3_UCSC ~/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome ~/pIC_exo_M3_MC1/pIC-exo-M3-MC1_S17_L001_R1_001.fastq.gz,~/pIC_exo_M3_MC1/pIC-exo-M3-MC1_S17_L002_R1_001.fastq.gz ~/pIC_exo_M3_MC1/pIC-exo-M3-MC1_S17_L001_R2_001.fastq.gz,~/pIC_exo_M3_MC1/pIC-exo-M3-MC1_S17_L002_R2_001.fastq.gz
  ```

  7. PBS M1:
    ```{sh}
  tophat2 --no-coverage-search --b2-very-sensitive --no-novel-juncs -p 4 -G ~/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/PBS_M1_MC1-31658713/PBS_M1_UCSC ~/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome ~/PBS_M1_MC1-31658713/PBS-M1-MC1_S1_L001_R1_001.fastq.gz,~/PBS_M1_MC1-31658713/PBS-M1-MC1_S1_L002_R1_001.fastq.gz ~/PBS_M1_MC1-31658713/PBS-M1-MC1_S1_L001_R2_001.fastq.gz,~/PBS_M1_MC1-31658713/PBS-M1-MC1_S1_L002_R2_001.fastq.gz
  ```

  8. PBS M2:
    ```{sh}
  tophat2 --no-coverage-search --b2-very-sensitive --no-novel-juncs -p 4 -G ~/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/PBS_M2_MC1-31644830/PBS_M2_UCSC ~/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome ~/PBS_M2_MC1-31644830/PBS-M2-MC1_S2_L001_R1_001.fastq.gz,~/PBS_M2_MC1-31644830/PBS-M2-MC1_S2_L002_R1_001.fastq.gz ~/PBS_M2_MC1-31644830/PBS-M2-MC1_S2_L001_R2_001.fastq.gz,~/PBS_M2_MC1-31644830/PBS-M2-MC1_S2_L002_R2_001.fastq.gz
  ```

4. Cufflinks:
  1. Control M1:
    ```{sh}
  /projects/home/sravishankar9/tools/tuxedo_suite/cufflinks/cufflinks -p 4 -G /data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/Control_Exo_M1_MC1-31664696/Control_ExoM1_UCSC_cuff ~/Control_Exo_M1_MC1-31664696/Control_ExoM1_UCSC/accepted_hits.bam
  ```

  2. Control M2:
    ```{sh}
  /projects/home/sravishankar9/tools/tuxedo_suite/cufflinks/cufflinks -p 4 -G /data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/Control_Exo_M2_MC1-31646815/Control_ExoM2_UCSC_cuff ~/Control_Exo_M2_MC1-31646815/Control_ExoM2_UCSC/accepted_hits.bam
  ```

  3. Control M3:
    ```
  /projects/home/sravishankar9/tools/tuxedo_suite/cufflinks/cufflinks -p 4 -G /data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/Control_Exo_M3_MC1/Control_ExoM3_UCSC_cuff ~/Control_Exo_M3_MC1/Control_ExoM3_UCSC/accepted_hits.bam
  ```

  4. pIC M1:
    ```{sh}
  /projects/home/sravishankar9/tools/tuxedo_suite/cufflinks/cufflinks -p 4 -G /data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/pIC_exo_M1_MC1-31656710/pIC_exoM1_UCSC_cuff ~/pIC_exo_M1_MC1-31656710/pIC_exoM1_UCSC/accepted_hits.bam
  ```

  5. pIC M2:
    ```{sh}
  /projects/home/sravishankar9/tools/tuxedo_suite/cufflinks/cufflinks -p 4 -G /data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/pIC_exo_M2_MC1-31646822/pIC_exoM2_UCSC_cuff ~/pIC_exo_M2_MC1-31646822/pIC_exoM2_UCSC/accepted_hits.bam
  ```

  6. pIC M3:
    ```{sh}
  /projects/home/sravishankar9/tools/tuxedo_suite/cufflinks/cufflinks -p 4 -G /data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/pIC_exo_M3_MC1/pIC_exoM3_UCSC_cuff ~/pIC_exo_M3_MC1/pIC_exoM3_UCSC/accepted_hits.bam  
  ```

  7. PBS M1:
    ```{sh}
  /projects/home/sravishankar9/tools/tuxedo_suite/cufflinks/cufflinks -p 4 -G /data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/PBS_M1_MC1-31658713/PBS_M1_UCSC_cuff ~/PBS_M1_MC1-31658713/PBS_M1_UCSC/accepted_hits.bam  
  ```

  8. PBS M2:
    ```{sh}
  /projects/home/sravishankar9/tools/tuxedo_suite/cufflinks/cufflinks -p 4 -G /data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o ~/PBS_M2_MC1-31644830/PBS_M2_UCSC_cuff ~/PBS_M2_MC1-31644830/PBS_M2_UCSC/accepted_hits.bam  
  ```

5. Cuffmerge:
  1. Prepare text file with paths to cufflinks transcripts files:
    ```{sh}
  cat ~/Control_Exo_M1_MC1-31664696/Control_ExoM1_UCSC_cuff/transcripts.gtf ~/Control_Exo_M2_MC1-31646815/Control_ExoM2_UCSC_cuff/transcripts.gtf ~/Control_Exo_M3_MC1/Control_ExoM3_UCSC_cuff/transcripts.gtf ~/pIC_exo_M1_MC1-31656710/pIC_exoM1_UCSC_cuff/transcripts.gtf ~/pIC_exo_M2_MC1-31646822/pIC_exoM2_UCSC_cuff/transcripts.gtf ~/pIC_exo_M3_MC1/pIC_exoM3_UCSC_cuff/transcripts.gtf ~/PBS_M1_MC1-31658713/PBS_M1_UCSC_cuff/transcripts.gtf ~/PBS_M2_MC1-31644830/PBS_M2_UCSC_cuff/transcripts.gtf > ~/assemblies.txt
  ```

  2. Run cuffmerge to prepare merged transcript file:
    ```{sh}
  cuffmerge -g /data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -s /data/db/Mus_musculus/UCSC/mm10/Annotation/Genes/genome.fa -p 8 ~/assemblies.txt
  ```

6. Cuffdiff:
  1. Run cuffdiff between control exosome and pIC exosome:
    ```{sh}
  cuffdiff -o controlvspIC -b /data/db/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -p 4 -L Control,pIC -u ~/merged_asm/merged.gtf ~/Control_Exo_M1_MC1-31664696/Control_ExoM1_UCSC/accepted_hits.bam,~/Control_Exo_M2_MC1-31646815/Control_ExoM2_UCSC/accepted_hits.bam,~/Control_Exo_M3_MC1/Control_Exo_M3_outputs_UCSC/accepted_hits.bam ~/pIC_exo_M1_MC1-31656710/pIC_M1_UCSC/accepted_hits.bam,~/pIC_exo_M2_MC1-31646822/pIC_M2_UCSC/accepted_hits.bam,~/pIC_exo_M3_MC1/pIC-Exo-M3_outputs_UCSC/accepted_hits.bam
  ```

  2. Run cuffdiff between PBS and pIC exosome:
    ```{sh}
  cuffdiff -o PBSvspIC -b /data/db/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -p 4 -L PBS,pIC -u ~/merged_asm/merged.gtf ~/PBS_M1_MC1-31658713/PBSM1_UCSC/accepted_hits.bam,~/PBS_M2_MC1-31644830/PBSM2_UCSC/accepted_hits.bam ~/pIC_M1_UCSC/accepted_hits.bam,~/pIC_M2_UCSC/accepted_hits.bam
  ```

 7. Expression plots:
   The following R code plots the expression of the genes mentioned in "~/GenesList.txt" to create the expression plots
   ```
 library(ggplot2)
 library(Hmisc)
 library(gridExtra)
 library(grid)
 library(reshape2)
 genes <- read.table('~/GenesList.txt', header=TRUE,sep='\t')
 genel <- as.vector(genes$Gene.name)
 genel <- capitalize(tolower(genel))
 names(genes)[1] <- "gene"
 genes$gene <- capitalize(tolower(genes$gene))
 difexp <- read.table('~/controlvspIC/gene_exp.diff',header=TRUE,sep='\t')
 difexp <- subset(difexp,select = c(gene,value_1,value_2))
 difexp2 <- read.table('~/PBSvspICM2/gene_exp.diff', header=TRUE,sep='\t')
 difexp2 <- subset(difexp2, select=c(gene,value_1))
 names(difexp)[2] <- "Con_exo"
 names(difexp)[3] <- "pIC_exo"
 names(difexp2)[2] <- "PBS"
 difexp <- merge(difexp,difexp2)
 difexp$Con_exo <- log2(difexp$Con_exo+1)
 difexp$pIC_exo <- log2(difexp$pIC_exo+1)
 difexp$PBS <- log2(difexp$PBS+1)
 dmg <- difexp[is.element(difexp$gene,genel),]
 dmg <- merge(dmg, genes)
 dmg <- dmg[order(dmg[,5],dmg[,1]),]
 row.names(dmg) <- dmg$gene
 dmgl <- dmg[,2:5]
 dmgl <- dmgl[order(dmgl[,4],dmgl[,1],decreasing=TRUE),]
 dmgl <- dmgl[,1:3]
 dmgl$Gene <- row.names(dmgl)
 dmgl.m <- melt(dmgl)
 names(dmgl.m)[1] <- "Gene"
 names(dmgl.m)[2] <- "Sample"
 names(dmgl.m)[3] <- "log2_FPKM"
 dmgl.m$Gene <- factor(dmgl.m$Gene,levels=dmgl.m$Gene)
 df.m <- melt(difexp)
 names(df.m)[3] <- "log2_FPKM"
 names(df.m)[2] <- "Condition"
 scatter <- ggplot(data = dmgl.m, aes(x = Gene, y = log2_FPKM, shape=Sample,  fill=Sample, aplha=0.9, stroke=1 )) + geom_point(size=5) + coord_flip() + scale_shape_manual(values = c(21,21,21)) + scale_fill_manual(values=c("grey","black","white"))  + theme_bw() + theme(legend.position=c(1,0.8),legend.justification=c(1,1),axis.text.x=element_text(face="bold"), axis.text.y=element_text(face="italic"), axis.line=element_line(size=1, linetype = "solid"), panel.background=element_blank(), panel.border =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(y=bquote('Gene Expression('~Log[2]~')'),x=element_blank())
 density <- ggplot(df.m, aes(x=log2_FPKM, fill=Condition,stroke=12)) + geom_density(alpha=.6) + scale_fill_manual(values=c("grey","black","white")) +theme_bw() +theme(legend.position = "none", axis.title.y=element_blank(), axis.title.x=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +scale_x_continuous(limits = range(dmgl.m$log2_FPKM))
 scat <- ggplot_gtable(ggplot_build(scatter))
 dens <- ggplot_gtable(ggplot_build(density))
 maxWidth = unit.pmax(scat$widths[2:3], dens$widths[2:3])
 scat$widths[2:3] <- maxWidth
 dens$widths[2:3] <- maxWidth
 pdf('ExpressionPlots.pdf')
 grid.arrange(dens, scat, heights = c(3, 16))
 dev.off()
 ```
   Expected output:

   ![Expression plot](ExpressionPlots.png)

8. Gage:
  1. Aggregating the bamfiles:
    ```{sh}
  cp ~/Control_Exo_M1_MC1-31664696/Control_ExoM1_UCSC/accepted_hits.bam ~/bamfiles/controlm1.bam
  cp ~/Control_Exo_M2_MC1-31646815/Control_ExoM2_UCSC/accepted_hits.bam ~/bamfiles/controlm2.bam
  cp ~/Control_Exo_M3_MC1/Control_Exo_M3_accepted_hits.bam ~/bamfiles/controlm2.bam
  cp ~/Control_Exo_M3_MC1/Control_Exo_M3_outputs_UCSC/accepted_hits.bam ~/bamfiles/controlm3.bam
  cp ~/pIC_exo_M1_MC1-31656710/pIC_M1_UCSC/accepted_hits.bam ~/bamfiles/picm1.bam
  cp ~/pIC_exo_M2_MC1-31646822/pIC_M2_UCSC/accepted_hits.bam ~/bamfiles/picm2.bam
  cp ~/pIC_exo_M3_MC1/pIC-Exo-M3_outputs_UCSC/accepted_hits.bam ~/bamfiles/picm3.bam
  cp ~/PBS_M2_MC1-31644830/PBSM2_UCSC/accepted_hits.bam ~/bamfiles/pbsm2.bam
  ```

  2. Indexing all bamfiles:
    ```{sh}
  cd bamfiles
  samtools index *
  ```

  3. Count reads mapped to each gene:
    ```
  source("http://bioconductor.org/biocLite.R")
  biocLite(c("pathview", "gage", "gageData", "GenomicAlignments", "TxDb.Mmusculus.UCSC.mm10.knownGene"))
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  exByGn <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, "gene")
  library(GenomicAlignments)
  fls <- list.files("~/bamfiles/", pattern="bam$", full.names =T)
  bamfls <- BamFileList(fls)
  flag <- scanBamFlag(isSecondaryAlignment=FALSE, isProperPair=TRUE)
  param <- ScanBamParam(flag=flag)
  gnCnt <- summarizeOverlaps(exByGn, bamfls, mode="Union", ignore.strand=TRUE, singleEnd=FALSE, param=param)
  hnrnp.cnts=assay(gnCnt) #read counts
  ```

  4. Plot significant altered pathways:
    1. Control vs pIC Top ten pathways:
      ```
    library(gage)
    library(pathview)
    library(ggplot2)
    library(reshape2)
    library(plyr)
    cnts.norm <- hnrnp.cnts
    con_exo <- 1:3
    pic_exo <- 5:7
    kg.mouse<- kegg.gsets("mouse")
    kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
    cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = con_exo, samp = pic_exo, compare ="as.group")
    kegg_pathway <- as.data.frame(cnts.kegg.p$greater)
    kegg_pathway$pathway <- row.names(kegg_pathway)
    kegg_pathway <- subset(kegg_pathway, select = c(p.val, pathway))
    names(kegg_pathway)[1] <- "Con_exoVSpIC_exo"
    kp.m <- melt(kegg_pathway)
    kp.m$value <- -1 * log10(kp.m$value)
    kp.d <- arrange(kp.m, desc(value))
    kp.d <- kp.d[1:10,]
    kp.d <- arrange(kp.d, value)
    kp.d$pathway <- factor(kp.d$pathway,levels=kp.d$pathway)
    names(kp.d)[2] <- "Sample"
    plist <- c("mmu00190 Oxidative phosphorylation" = "Oxidative phosphorylation", "mmu04145 Phagosome" = "Phagosome", "mmu04380 Osteoclast differentiation" = "Oseteoclast differentiation", "mmu04142 Lysosome" = "Lysosome", "mmu03040 Spliceosome" = "Spliceosome", "mmu04064 NF−kappa B signaling pathway" = "NF-kappa B signaling pathway", "mmu03050 Proteasome" = "Proteasome", "mmu04062 Chemokine signaling pathway" = "Chemokine signaling pathway", "mmu04210 Apoptosis" = "Apoptosis", "mmu03010 Ribosome" = "Ribosome")
    pdf('Top10ControlvspIC1.pdf', width=10, height=8, )
    ggplot(data=kp.d, aes(x= pathway, y= value)) + geom_bar(stat="identity", width=0.5,position=position_dodge(0.1), colour="black", fill="black", size=1) + coord_flip() + theme(axis.title = element_text(size=20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size=18), panel.border =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank()) + labs(x=NULL, y='-Log10(p-value)', size=18) + scale_x_discrete(labels = plist)
    dev.off()
    ```
      Expected output:

      ![Top ten](cvptt.PNG)

    2. Control vs pIC pathways of interest:

      | Pathways                                  | Kegg ID     |
      | :---------------------------------------- | :---------- |
      | Chemokine signaling                       |   mmu04062  |
      | NFKB signaling pathway                    |   mmu04064  |
      | Toll like receptor signaling pathway      |   mmu04620  |
      | Antigen processing and presentation       |   mmu04612  |
      | TNF signaling pathway                     |   mmu04668  |
      | Cytosolic DNA sensing pathway             |   mmu04623  |

      ```
    library(gage)
    library(pathview)
    library(ggplot2)
    library(reshape2)
    library(plyr)
    cnts.norm <- hnrnp.cnts
    glist <- c("mmu04062 Chemokine signaling pathway", "mmu04064 NF-kappa B signaling pathway", "mmu04620 Toll-like receptor signaling pathway", "mmu04612 Antigen processing and presentation", "mmu04668 TNF signaling pathway", "mmu04623 Cytosolic DNA-sensing pathway")
    plist <- c("mmu04062 Chemokine signaling pathway" = "Chemokine signaling pathway", "mmu04064 NF-kappa B signaling pathway" = "NF-kappa B signaling pathway", "mmu04620 Toll-like receptor signaling pathway" = "Toll-like receptor signaling pathway", "mmu04612 Antigen processing and presentation" = "Antigen processing and presentation", "mmu04668 TNF signaling pathway" = "TNF signaling pathway", "mmu04623 Cytosolic DNA-sensing pathway" = "Cytosolic DNA-sensing pathway")
    con_exo <- 1:3
    pic_exo <- 5:7
    kg.mouse<- kegg.gsets("mouse")
    kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
    cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = con_exo, samp = pic_exo, compare ="as.group")
    kegg_pathway <- as.data.frame(cnts.kegg.p$greater)
    kegg_pathway$pathway <- row.names(kegg_pathway)
    kegg_pathway <- subset(kegg_pathway, select = c(p.val, pathway))
    names(kegg_pathway)[1] <- "Con_exoVSpIC_exo"
    kp.m <- melt(kegg_pathway)
    kp.m$value <- -1 * log10(kp.m$value)
    kp.d <- arrange(kp.m, desc(value))
    kp.d <- kp.d[is.element(kp.d$pathway, glist),]
    kp.d <- arrange(kp.d, value)
    kp.d$pathway <- factor(kp.d$pathway,levels=kp.d$pathway)
    names(kp.d)[2] <- "Sample"
    pdf('POIControlvspIC1.pdf', width=10, height=6)
    ggplot(data=kp.d, aes(x= pathway, y= value)) + geom_bar(stat="identity", width=0.5,position=position_dodge(0.1), colour="black", fill="black", size=1) + coord_flip() + theme(axis.title = element_text(size=20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size=18), panel.border =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank()) + labs(x=NULL, y='-Log10(p-value)', size=18) + scale_x_discrete(labels = plist)
    dev.off()
    ```
      Expected output:

      ![Pathways of interest](cvppoi.PNG)

    3. PBS vs pIC Top ten pathways:

      ```
    library(gage)
    library(pathview)
    library(ggplot2)
    library(reshape2)
    library(plyr)
    cnts.norm <- hnrnp.cnts
    pbs <- 4
    pic_exo <- 5:7
    kg.mouse<- kegg.gsets("mouse")
    kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
    cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = pbs, samp = pic_exo, compare ="as.group")
    kegg_pathway <- as.data.frame(cnts.kegg.p$greater)
    kegg_pathway$pathway <- row.names(kegg_pathway)
    kegg_pathway <- subset(kegg_pathway, select = c(p.val, pathway))
    names(kegg_pathway)[1] <- "PBSVSpIC_exo"
    kp.m <- melt(kegg_pathway)
    kp.m$value <- -1 * log10(kp.m$value)
    kp.d <- arrange(kp.m, desc(value))
    kp.d <- kp.d[1:10,]
    kp.d <- arrange(kp.d, value)
    kp.d$pathway <- factor(kp.d$pathway,levels=kp.d$pathway)
    names(kp.d)[2] <- "Sample"
    plist <- c("mmu00190 Oxidative phosphorylation" = "Oxidative phosphorylation", "mmu04145 Phagosome" = "Phagosome", "mmu04380 Osteoclast differentiation" = "Osteoclast differentiation", "mmu04142 Lysosome" = "Lysosome", "mmu03040 Spliceosome" = "Spliceosome", "mmu04064 NF−kappa B signaling pathway" = "NF-kappa B signaling pathway", "mmu04668 TNF signaling pathway" = "TNF signaling pathway", "mmu03050 Proteasome" = "Proteasome", "mmu03010 Ribosome" = "Ribsome", "mmu04062 Chemokine signaling pathway" = "Chemokine signaling pathway")
    pdf('Top10PBSvspIC1.pdf', width=10, height=8)
    ggplot(data=kp.d, aes(x= pathway, y= value)) + geom_bar(stat="identity", width=0.5,position=position_dodge(0.1), colour="black", fill="black", size=1) + coord_flip() + theme(axis.title = element_text(size=20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size=18), panel.border =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank()) + labs(x=NULL, y='-Log10(p-value)', size=18) + scale_x_discrete(labels = plist)
    dev.off()
    ```

      Expected output:

      ![Top ten](pvptt.PNG)

    4. PBS vs pIC pathways of interest:

      | Pathways                                  | Kegg ID     |
      | :---------------------------------------- | :---------- |
      | MAPK signaling pathway                    |   mmu04010  |
      | Cytokine-Cytokine receptor interactions   |   mmu04060  |
      | Chemokine signaling                       |   mmu04062  |
      | NFKB signaling pathway                    |   mmu04064  |
      | Toll like receptor signaling pathway      |   mmu04620  |
      | Antigen processing and presentation       |   mmu04612  |
      | NOD like receptor signaling pathway       |   mmu04621  |
      | Rig-1 like receptor signaling pathway     |   mmu04622  |
      | Snare interactions in vesicular transport |   mmu04130  |
      | Cytosolic DNA sensing pathway             |   mmu04623  |
      | JAK Stat sensing pathway                  |   mmu04630  |

      ```
    library(gage)
    library(pathview)
    library(ggplot2)
    library(reshape2)
    library(plyr)
    cnts.norm <- hnrnp.cnts
    glist <- c("mmu04010 MAPK signaling pathway", "mmu04060 Cytokine-cytokine receptor interaction", "mmu04062 Chemokine signaling pathway", "mmu04064 NF-kappa B signaling pathway", "mmu04620 Toll-like receptor signaling pathway", "mmu04612 Antigen processing and presentation", "mmu04621 NOD-like receptor signaling pathway", "mmu04622 RIG-I-like receptor signaling pathway", "mmu04130 SNARE interactions in vesicular transport", "mmu04623 Cytosolic DNA-sensing pathway", "mmu04630 Jak-STAT signaling pathway")
    pbs <- 4
    pic_exo <- 5:7
    kg.mouse<- kegg.gsets("mouse")
    kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
    cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = pbs, samp = pic_exo, compare ="as.group")
    kegg_pathway <- as.data.frame(cnts.kegg.p$greater)
    kegg_pathway$pathway <- row.names(kegg_pathway)
    kegg_pathway <- subset(kegg_pathway, select = c(p.val, pathway))
    names(kegg_pathway)[1] <- "PBSVSpIC_exo"
    kp.m <- melt(kegg_pathway)
    kp.m$value <- -1 * log10(kp.m$value)
    kp.d <- arrange(kp.m, desc(value))
    kp.d <- kp.d[is.element(kp.d$pathway, glist),]
    kp.d <- arrange(kp.d, value)
    kp.d$pathway <- factor(kp.d$pathway,levels=kp.d$pathway)
    names(kp.d)[2] <- "Sample"
    plist <- c("mmu04064 NF-kappa B signaling pathway" = "NF-kappa B signaling pathway", "mmu04062 Chemokine signaling pathway" = "Chemokine signaling pathway", "mmu04612 Antigen processing and presentation" = "Antigen processing and presentation", "mmu04620 Toll-like receptor signaling pathway" = "Toll-like receptor signaling pathway", "mmu04010 MAPK signaling pathway" = "MAPK signaling pathway", "mmu04621 NOD-like receptor signaling pathway" = "NOD-like receptor signaling pathway", "mmu04623 Cytosolic DNA-sensing pathway" = "Cytosolic DNA-sensing pathway", "mmu04060 Cytokine-cytokine receptor interaction" = "Cytokine-cytokine receptor interaction", "mmu04630 Jak-STAT signaling pathway" = "Jak-STAT signaling pathway", "mmu04130 SNARE interactions in vesicular transport" = "SNARE interactions in vesicular transport", "mmu04622 RIG-I-like receptor signaling pathway" = "RIG-I-like receptor signaling pathway")
    pdf('POIPBSvspIC1.pdf', width=10, height=8)
    ggplot(data=kp.d, aes(x= pathway, y= value)) + geom_bar(stat="identity", width=0.5,position=position_dodge(0.1), colour="black", fill="black", size=1) + coord_flip() + theme(axis.title = element_text(size=20), axis.text.y = element_text(size = 18), axis.text.x = element_text(size=18), panel.border =element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank()) + labs(x=NULL, y='-Log10(p-value)', size=18) + scale_x_discrete(labels = plist)
    dev.off()
    ```    

     Expected output:

      ![Pathways of interest](pvppoi.PNG)

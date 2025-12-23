# ============================================================
# 性能优化版本 - 添加缓存和并行处理
# ============================================================

# hg19和38数据库，需要建立索引
hg19 <- "/Users/fangy/fsdownload/blastdb/ucsc.hg19.fasta"
hg38 <- "/Users/fangy/fsdownload/hg38/hg38.fa"
#blastn工具，需要本地部署
bn <- "/Users/fangy/Downloads/ncbi-blast-2.15.0+/bin/blastn"
#balst使用的hg19数据库，需要本地部署
bdb <- "/Users/fangy/fsdownload/blastdb/hg19"

#绝对路径，需要指向当前目录
seqs_path <- "/Users/fangy/Desktop/sanger"

#bash中需要部署的第三方软件（bash调用）
# 1. samtools
# 2. Primer3

library(Biostrings)
library(memoise)
library(cachem)

# 创建缓存对象 - 最大缓存512MB，过期时间1小时
cache <- cache_mem(max_size = 512 * 1024^2, max_age = 3600)

# 缓存 samtools faidx 调用结果
get_sequence_cached <- memoise(function(hg, genome_loc) {
  cmd_samtools <- paste0("samtools faidx ", hg, " ", genome_loc)
  result <- system(cmd_samtools, intern = TRUE)
  if (length(result) > 1) {
    return(paste(result[2:length(result)], collapse = ""))
  }
  return(result[2])
}, cache = cache)
figheight <- function(obj, trim5=0, trim3=0, width=100, pixelsperrow=200, showtrim=FALSE) {
  traces <- obj@traceMatrix
  basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
  basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
  aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
  basecalls1 <- basecalls1[1:length(aveposition)] #####
  basecalls2 <- basecalls2[1:length(aveposition)] ######
  if(showtrim == FALSE) {
    if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
    else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
    if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
    else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
    aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
  }
  indexes <- 1:length(basecalls1)
  trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all false if not trimmed
  if (!is.null(trim3)) {
    traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, nrow(traces))),]
  }
  if (!is.null(trim5)) {
    offset <- max(c(1, aveposition[1] - 10))
    traces <- traces[offset:nrow(traces),]
    aveposition <- aveposition - (offset-1)
  }
 
  valuesperbase <- nrow(traces)/length(basecalls1)
  tracewidth <- width*valuesperbase
  breaks <- seq(1,nrow(traces), by=tracewidth) 
  
  numplots <- length(breaks)
  return(numplots*pixelsperrow)
}

wrap_fixed = function(string, width=80) {
  pattern = paste("(.{1,", width, "})", sep="")
  res = gsub(pattern, "\\1\n", string)
  return(res)
}

cleanstring <- function(string) {
  string <- gsub("^>.*?\n", "", string)
  string <- toupper(string)
  string <- gsub("[^ACGTRYSWKMBDHVN]", "", string, perl=TRUE)
  return(string)
}

alignchromatogram <- function(data, block.width=50, trim=FALSE, refseq, trim5, trim3) {
  if (is.null(data)) return(NULL)
  d <- setAllelePhase(data, refseq, trim5, trim3)
  refseq <- toString(d@primarySeq)
  altseq <- toString(d@secondarySeq)
  if (trim == TRUE) {
    altseq <- toString(d@secondarySeq[(trim5 + 1):(nchar(altseq) - trim3)])
  }
  names(altseq) <- "Alt Allele"
  names(refseq) <- d@primarySeqID
  if (trim == TRUE) {
    pa <- pairwiseAlignment(altseq, refseq, type="local", gapExtension=-2)
  } else {
    pa <- pairwiseAlignment(altseq, refseq, type="global", gapExtension=-2)
  }
  alignment <- paste(capture.output(writePairwiseAlignments(pa, block.width=block.width)), collapse="\n")
  results <- list(altseq=altseq, refseq=gsub("-", "", refseq), alignment=alignment)
  return(results)
}

pickchromatogram <- function(obj, trim5=0, trim3=0, width=100, pixelsperrow=200, showtrim=FALSE,inputloc,selectedGenomeVersion, custom_seq=NULL){
  extend_len <- 10
  if(selectedGenomeVersion == 'hg19'){
    hg <- hg19
  }
  if(selectedGenomeVersion == 'hg38'){
    hg <- hg38
  }
  
  
  
  # 坐标提取序列
  # s = readDNAStringSet(hg19)
  # 提取特定坐标的序列 13-20763485-AG-A
  # inputloc <- '1-11796321-G-A'
  inputloc_list <- strsplit(inputloc,"-")[[1]]
  locchr <- paste0('chr', inputloc_list[1])
  locid <- as.numeric(inputloc_list[2])
  # locchr <- "chr13"
  # locid  <- 20763485
  s1<- locid - extend_len
  e1 <- locid + extend_len + max(nchar(inputloc_list[[3]]),nchar(inputloc_list[[4]])) * 2
  # extracted_sequence <- subseq(s[[locchr]], start = s1, end = e1)
  #GATCAGCTGCAGGGCCCATAG
  genome_loc <- paste0(locchr,":",as.character(s1),"-",as.character(e1))
  
  if (selectedGenomeVersion == 'Custom') {
    if (is.null(custom_seq) || custom_seq == "") {
      stop("Custom sequence is required when 'Custom' genome is selected.")
    }
    # Clean up custom sequence (remove whitespace, newlines, numbers)
    custom_seq_clean <- gsub("[^a-zA-Z]", "", custom_seq)
    
    # Check bounds
    seq_len <- nchar(custom_seq_clean)
    
    # We need strictly extend_len padding for validation logic to work?
    # Actually check_del etc rely on offset from start of extracted_sequence.
    # If s1 < 1, we might need to pad?
    # For now, clamp and warn if too close to edge.
    
    start_pos <- max(1, s1)
    end_pos <- min(seq_len, e1)
    
    if (start_pos > end_pos) {
       stop("Calculated coordinates are outside the custom sequence range.")
    }
    
    seq_result <- substr(custom_seq_clean, start_pos, end_pos)
    extracted_sequence <- DNAString(seq_result)
    
  } else {
    # 使用缓存版本的序列获取
    seq_result <- get_sequence_cached(hg, genome_loc)
    extracted_sequence <- DNAString(seq_result)
  }
  
  
  
  #读取sanger测序结果
  # setwd(seqs_path)
  # setwd("/Users/fangy/Desktop/sanger/testupload2")
  # seq <- makeBaseCalls(readsangerseq("/Users/fangy/Desktop/sanger/test2/3990-66_66-R_TSS20231226-0371-02131_B01.ab1"), ratio=0.33)
  # seq <- obj
  seq <- makeBaseCalls(readsangerseq(obj), ratio = 0.33)
  # seq <- makeBaseCalls(obj, ratio = 0.33)
  
  # seq@primarySeq <- seq@primarySeq[1:200]
  # 10-letter DNAString object
  # seq: GTGACGTCAT
  seqstr <- primarySeq(seq,string=TRUE)
  # seqstr <- secondarySeq(seq,string=TRUE)
  # subseq <- substr(seqstr,226,nchar(seqstr))
  
  #比对目标序列
  # extracted_sequence <- "GTGATCGTAGCACACGTTCTTGCAGCCTGGCTGCAGGGTGTTGCAGACAA"
  # seqstr <- "GTGATCGTAGCTGGCTGCAGGGTGTTGCAGACAA"
  
  pa <- pairwiseAlignment(extracted_sequence,seqstr,
                          type = "global-local",gapOpening=20, gapExtension=2)
  writePairwiseAlignments(pa)
  
  start_value <- start(pa@subject@indel[[1]])
  width_value <- width(pa@subject@indel[[1]])
  
  #验证纯合缺失
  check_del_het <- function(start_value,width_value){
    
    if (start_value == extend_len + 2 & width_value == nchar(inputloc_list[[3]]) -1 ){
      # print("Del, Heterozygous")
      return("Homozygous")
    }else{
      return("Heterozygous")
    }
  }
  # check_del_het(start_value,width_value)
  # P1                1 GATCAGCTGCAGGGCCCATAG     21
  #                     ||||||||||||||||| || 
  #   S1            212 GATCAGCTGCAGGGCCCTTAC    232
  
  # subject_start <- pa@subject@range@start
  # subject_end <- pa@subject@range@start + pa@subject@range@width - 1
  # 
  # subject_loc <- pa@subject@range@start + extend_len
  # 
  # #subject DNA str
  # subject_str <- subseq(primarySeq(seq, string = TRUE), start = subject_start, width = pa@subject@range@width)
  # 
  
  
  subject_start <- pa@subject@range@start
  subject_end <- pa@subject@range@start + pa@subject@range@width - 1
  
  subject_loc <- pa@subject@range@start + extend_len
  
  subject_str <- subseq(primarySeq(seq, string = TRUE), start = subject_start, width = pa@subject@range@width)
  
  subject_str_secondary <- subseq(secondarySeq(seq, string = TRUE), start = subject_start, width = pa@subject@range@width)
  
  check_start_primary <- subseq(primarySeq(seq, string = TRUE), start = subject_loc, width = extend_len+1)
  
  check_start_secondary <- subseq(secondarySeq(seq, string = TRUE), start = subject_loc, width = extend_len+1)
  
  check_ref <- subseq(extracted_sequence, start = extend_len+1, end =nchar(extracted_sequence))
  
  # subject_str[subject_loc:nchar(subject_str)]
  
  check_del <- function(subject_str,subject_str_secondary,extracted_sequence){
    ref_str <- as.character(extracted_sequence)
    # del_loc <- extend_len + 2
    del_loc <- extend_len + 2
    # for (forward in 1:del_loc) {
    #   
    #   
    # }
    # ref_end_loc <- nchar(extracted_sequence)-1
    ref_end_loc <- del_loc + 5
    backward_vec <- vector()
    for (backward in del_loc:ref_end_loc) {
      
      # backward <-12
      p_str <- substr(subject_str,backward,backward)
      s_str <- substr(subject_str_secondary,backward,backward)
      ref_wark <- substr(ref_str,backward,backward+1)
      ps_wark <- paste0(p_str,s_str)
      if (sum(sort(strsplit(ref_wark, "")[[1]]) == sort(strsplit(ps_wark, "")[[1]]))==2){
        backward_vec<-c(backward_vec,TRUE)
      }else{
        backward_vec<-c(backward_vec,FALSE)
      }
      
    }
    print(backward_vec)
    if(sum(backward_vec) == length(del_loc:ref_end_loc)){
      print("Del")
      return("Del")
    }else{
        print("Not Del")
      return("Not Del")
    }
    
  }
  
  check_del_multiple <- function(subject_str,subject_str_secondary,extracted_sequence){
    ref_str <- as.character(extracted_sequence)
    # del_loc <- extend_len + 2
    del_loc <- extend_len + 2
    del_len <- nchar(inputloc_list[[3]]) -1
    # for (forward in 1:del_loc) {
    #   
    #   
    # }
    ref_end_loc <- nchar(extracted_sequence)-1
    backward_vec <- vector()
    del_loc_end <- del_loc+del_len
    
    for (backward in del_loc: del_loc_end) {
      
      # backward <-12
      p_str <- substr(subject_str,backward,backward)
      s_str <- substr(subject_str_secondary,backward,backward)
      ref_wark <- substr(ref_str,backward,backward)
      ref_wark2 <- substr(ref_str,backward + del_len ,backward + del_len)
      refref_wark <- paste0(ref_wark,ref_wark2)
      ps_wark <- paste0(p_str,s_str)
      if (sum(sort(strsplit(refref_wark, "")[[1]]) == sort(strsplit(ps_wark, "")[[1]]))==2){
        backward_vec<-c(backward_vec,TRUE)
      }else{
        backward_vec<-c(backward_vec,FALSE)
      }
      
    }
    print(backward_vec)
    if(sum(backward_vec) == length(del_loc:del_loc_end)){
      print("Del")
      return("Del")
    }else{
      print("Not Del")
      return("Not Del")
    }
    
  }
  
  check_ins <- function(subject_str,subject_str_secondary,extracted_sequence){
    ref_str <- as.character(extracted_sequence)
    ins_loc <- extend_len + 2
    # for (forward in 1:del_loc) {
    #   
    #   
    # }
    # ref_end_loc <- nchar(extracted_sequence)-1
    ref_end_loc <- ins_loc + 5
    backward_vec <- vector()
    for (backward in ins_loc:ref_end_loc) {
      
      # backward <-21
      p_str <- substr(subject_str,backward + 1,backward + 1)
      s_str <- substr(subject_str_secondary,backward + 1,backward + 1)
      ref_wark <- substr(ref_str,backward,backward+1)
      ps_wark <- paste0(p_str,s_str)
      if (sum(sort(strsplit(ref_wark, "")[[1]]) == sort(strsplit(ps_wark, "")[[1]]))==2){
        backward_vec<-c(backward_vec,TRUE)
      }else{
        backward_vec<-c(backward_vec,FALSE)
      }
      
    }
    if(sum(backward_vec) == length(ins_loc:ref_end_loc)){
      print("Ins")
      return("Ins")
    }else{
      print("Not Ins")
      return("Not Ins")
    }
    
  }
  
  check_ins_multiple <- function(subject_str,subject_str_secondary,extracted_sequence){
    ref_str <- as.character(extracted_sequence)
    ins_loc <- extend_len + 1
    ins_len <- nchar(inputloc_list[[3]])
    # for (forward in 1:del_loc) {
    #   
    #   
    # }
    ref_end_loc <- nchar(extracted_sequence)-1
     
    ins_loc_end <- ins_loc + ins_len
    backward_vec <- vector()
    for (backward in ins_loc_end:ref_end_loc) {
      
      # backward <-15
      p_str <- substr(subject_str,backward ,backward )
      s_str <- substr(subject_str_secondary,backward ,backward)
      ref_wark <- substr(ref_str,backward,backward)
      ref_wark2 <- substr(ref_str,backward - ins_len,backward - ins_len)
      ps_wark <- paste0(p_str,s_str)
      refref_wark <- paste0(ref_wark,ref_wark2)
      if (sum(sort(strsplit(refref_wark, "")[[1]]) == sort(strsplit(ps_wark, "")[[1]]))==2){
        backward_vec<-c(backward_vec,TRUE)
      }else{
        backward_vec<-c(backward_vec,FALSE)
      }
      
    }
    
    check_dup_vec <- vector()
    dup_end_loc <-ins_loc_end +nchar(inputloc_list[[3]]) -1
    for (backward in ins_loc_end:dup_end_loc) {
      # backward <-15
      # p_str <- substr(subject_str,backward ,backward )
      # s_str <- substr(subject_str_secondary,backward ,backward)
      ref_wark <- substr(ref_str,backward,backward)
      ref_wark2 <- substr(ref_str,backward - ins_len,backward - ins_len)
      # ps_wark <- paste0(p_str,s_str)
      refref_wark <- paste0(ref_wark,ref_wark2)
      forward_ref <- substr(ref_str,backward-4 ,backward-4 )
      check_dup_vec <- c(check_dup_vec,grepl(forward_ref, refref_wark))
    }
    
    
    if(sum(backward_vec) == length(ins_loc_end:ref_end_loc)){
      
      if(sum(check_dup_vec) == ins_len){
        print("Dup")
        return("Dup")
      } else{
        print("Ins")
        return("Ins")
        }
    }else{
      print("Not Ins")
      return("Not Ins")
    }
    
  }
  
  check_ins_het <- function(subject_str,subject_str_secondary){
    ins_loc <- extend_len + 2
    
    if (substr(subject_str,ins_loc,nchar(subject_str)) ==  substr(subject_str_secondary,ins_loc,nchar(subject_str_secondary))){
      # print("Del, Heterozygous")
      return("Homozygous")
    }else{
      return("Heterozygous")
    }
    
  }
  
  
  
  
  
  check_mis <- function(subject_str,subject_str_secondary,extracted_sequence){
    loc <- extend_len + 1
    ref_str <- as.character(extracted_sequence)
    if(substr(subject_str,loc,loc)  !=  substr(ref_str,loc,loc) | substr(subject_str_secondary,loc,loc)  !=  substr(ref_str,loc,loc)){
      loc_type <- "Substution"
    } else{
      loc_type <- "Not Substution"
    }
    
    if(substr(subject_str,loc,loc)  ==  substr(subject_str_secondary,loc,loc)){
      loc_peaks <- "Homozygous"
    } else{
      loc_peaks <- "Heterozygous"
    }
    res <- list(het=loc_peaks,mu=loc_type)
    return(res)
  }
  
  
  # check wild type or variant 
  # check_mu <- function(p1,s1,loc=11){
  #   if (substr(p1,loc,loc) ==  substr(s1,loc,loc)){
  #     return("Wild type")
  #   }else{
  #     return("Variant")
  #   }
  # }
  
  extracted_sequence_str <-as.character(extracted_sequence)
  
  # loc_type <- check_mu(extracted_sequence_str,subject_str,11)
  # print(loc_type)
  # 
  
  # 峰型分析函数
  # analyze_peaks_loc <- function(ab1_data,subject_loc,threshold=400) {
  #   # ab1_data <- seq
  #   # 获取当前位点的所有峰值
  #   peaks <- ab1_data@traceMatrix[ab1_data@peakPosMatrix[subject_loc],]
  #   # 判断峰值数量
  #   significant_peaks <- sum(peaks > threshold)  # 需要根据数据设置适当的阈值
  #   if (significant_peaks > 1) {
  #     return("heterozygous") #"杂合型"
  #   } else {
  #     return("homozygous") # #"纯合型"
  #   }
  #   
  # }
  
  #判断位点的纯合还是杂合
  # loc_peaks <- analyze_peaks_loc(seq,subject_loc)
  # print(loc_peaks)
  
  
  #截取关注位点前后10bp图像用于判断
  seq2 <- seq
  seq2@primarySeq <- DNAString(subject_str)
  seq2@secondarySeq <- DNAString(subject_str_secondary)
  seq2@peakPosMatrix <-  seq2@peakPosMatrix[subject_start:subject_end,]
  # seq@peakAmpMatrix <-  seq@peakAmpMatrix[subject_start:subject_end,]
  # chromatogram(seq2, width = 30, height = 3, trim5 = 0, trim3 = 0,showcalls = "both")
  #
  
  alignment <- paste(capture.output(writePairwiseAlignments(pa, block.width=50)), collapse="\n")
  
  #单碱基缺失
  if(nchar(inputloc_list[3]) > nchar(inputloc_list[4]) & nchar(inputloc_list[3]) == 2 ){
    # return('del')
    loc_type <- check_del(subject_str,subject_str_secondary,extracted_sequence)

    if (length(start_value) == 0) {
      loc_peaks <- "Heterozygous"
    } else {
      loc_peaks <- check_del_het(start_value,width_value)
    }

  }
  #长缺失
  if(nchar(inputloc_list[3]) > nchar(inputloc_list[4]) & nchar(inputloc_list[3]) > 2 ){
    # return('del')
    loc_type <- check_del_multiple(subject_str,subject_str_secondary,extracted_sequence)
    
    if (length(start_value) == 0) {
      loc_peaks <- "Heterozygous"
    } else {
      loc_peaks <- check_del_het(start_value,width_value)
    }
    
  }
  #单碱基ins
  if(nchar(inputloc_list[3]) < nchar(inputloc_list[4]) & nchar(inputloc_list[4]) == 2 ){
    # return('ins')
    loc_type <- check_ins(subject_str,subject_str_secondary,extracted_sequence)
    
    loc_peaks <- check_ins_het(subject_str,subject_str_secondary)
  }
  
  #长碱基ins
  if(nchar(inputloc_list[3]) < nchar(inputloc_list[4])  & nchar(inputloc_list[4]) > 2 ){
    # return('ins')
    loc_type <- check_ins_multiple(subject_str,subject_str_secondary,extracted_sequence)
    
    loc_peaks <- check_ins_het(subject_str,subject_str_secondary)
  }
  if(nchar(inputloc_list[3]) == nchar(inputloc_list[4])){
    # return('Substution')
    mis <- check_mis(subject_str,subject_str_secondary,extracted_sequence)
    loc_peaks <- mis$het
    loc_type <- mis$mu
  }

  
  
  res <- list(seq2=seq2,
              het = loc_peaks, 
              mu = loc_type,
              altseq=seqstr, 
              refseq=gsub("-", "", extracted_sequence),
              alignment=alignment,
              start_loc =s1,
              end_loc=e1)
  return(res)
  
  
}


pickchromatogram2 <- function(inputinfo,selectedGenomeVersion, custom_seq=NULL){
  # inputinfo <- data.frame(GenomicCoordinate=c("17-41245245-CT-C","13-32893229-G-A","13-20763485-AG-A"),
  #                         FilePath=c("/Users/fangy/Desktop/sanger/testupload/0037_32122081601256_(BRCA2252-1)_[BRCA1-R].ab1",
  #                                    "/Users/fangy/Desktop/sanger/testupload/BRCA225183_BRCA2251-83F_TSS20220804-029-00668_G01.ab1",
  #                                    "/Users/fangy/Desktop/sanger/testupload/0007_32123111000861_(ZHC1)_[GJB2-2F].ab1"))
  check_loc <- function(inputinfo,mu){
    # inputloc <- "13-20763485-AG-A"
    # mu <- mu_values
    
    
          perform_operation <- function(row) {
          inputloc_list <- strsplit(row["GenomicCoordinate"],"-")[[1]]
          # print(inputloc_list)
          if(nchar(inputloc_list[3]) > nchar(inputloc_list[4])){
            return('Del')
          }
          if(nchar(inputloc_list[3]) < nchar(inputloc_list[4]) & nchar(inputloc_list[3]) != 1){
            return('Dup')
          }
          if(nchar(inputloc_list[3]) < nchar(inputloc_list[4])){
            return('Ins')
          }

          if(nchar(inputloc_list[3]) == nchar(inputloc_list[4])){
            return('Substution')
          }
        }
      
        input_info <- apply(inputinfo, 1, perform_operation)
        
        
        res <- vector()
        for (each in 1:length(mu)) {
          if (mu[each] == input_info[each]){
            res<-c(res,TRUE)}
          else{
            res <- c(res,FALSE)
          }
          
        }
        return(res)
      

    
  }
  
  # pickchromatogram(obj, trim5=0, trim3=0, width=100, pixelsperrow=200, showtrim=FALSE,inputloc)
  
  # 使用 future_lapply 实现并行处理，显著提升多文件处理速度
  results <- future.apply::future_lapply(1:nrow(inputinfo), function(i) {
    obj <- inputinfo$FilePath[i] 
    inputloc <- inputinfo$GenomicCoordinate[i]
    genome_ver <- selectedGenomeVersion
    # 调用pickchromatogram函数，并传入其他需要的参数
    pickchromatogram(obj, trim5=0, trim3=0, width=100, pixelsperrow=200, showtrim=FALSE, inputloc=inputloc,selectedGenomeVersion=genome_ver, custom_seq=custom_seq)
  }, future.seed = TRUE)
  
  chromatogram_pick <- sapply(results, function(x) x$seq2)
  mu_values <- sapply(results, function(x) x$mu)
  het_values <- sapply(results, function(x) x$het)
  
  Genotype <- paste(mu_values,het_values,sep=",")
    
  start_loc_values <- sapply(results, function(x) x$start_loc)
  end_loc_values <- sapply(results, function(x) x$end_loc)
  
  altseq_values <- sapply(results, function(x) x$altseq)
  refseq_values <- sapply(results, function(x) x$refseq)
  alignment_values <- sapply(results, function(x) x$alignment)
  
  
  
  
  print(inputinfo)
  res_table <- data.frame(
    Location = inputinfo$GenomicCoordinate,
    Genotype = Genotype,
    Check = check_loc(inputinfo,mu_values),
    Start = start_loc_values,
    End = end_loc_values,
    Figure = seq_len(nrow(inputinfo))
  )
  
  res <- list(res_table=res_table,
              chromatogram_pick =chromatogram_pick,
              altseq = altseq_values,
              refseq = refseq_values,
              alignment = alignment_values)
  return(res)
 # data.frame(
 #    Location = c("Location1", "Location2", "Location3"),
 #    Genotype = c("TypeA", "TypeB", "TypeC"),
 #    Check = c(TRUE, FALSE, TRUE),
 #    Start = c(100, 200, 300),
 #    End = c(150, 250, 350),
 #    Figure = seq_len(3)
 #  )

}

generate <- function(input) {
  # 这里可以根据input参数进行复杂的数据生成逻辑
  # 以下是临时示例数据
  data.frame(
    GenomeLocation = c("Location1", "Location2", "Location3"),
    Genotype = c("TypeA", "TypeB", "TypeC"),
    Start = c(100, 200, 300),
    End = c(150, 250, 350),
    Check = c(TRUE, FALSE, TRUE), # 示例的布尔值
    Figure = seq_len(3) # 用于标识每个按钮
  )
}
primerRes <- function(inputloc,selectedGenomeVersion, custom_seq=NULL,
                      primer_min_size=18, primer_opt_size=20, primer_max_size=27,
                      primer_min_tm=58, primer_opt_tm=60, primer_max_tm=64,
                      primer_min_gc=35, primer_opt_gc=50, primer_max_gc=65,
                      product_min_size=300, product_max_size=1000,
                      region_upstream=500, region_downstream=500){

  # s = readDNAStringSet(hg19)
  # 提取特定坐标的序列 13-20763485-AG-A
  # inputloc <- "13-20763485-AG-A"
  inputloc_list <- strsplit(inputloc,"-")[[1]]
  locchr <- paste0('chr', inputloc_list[1])
  locid <- as.numeric(inputloc_list[2])

  # locid  <- 20763485
  s1<- locid - region_upstream
  e1 <- locid + region_downstream
  # seqs <- subseq(s[[locchr]], start = s1, end = e1)
  # names(seqs) <- "1000bp seq"
  # seq.ids <- inputloc
  # hg <- hg19
  if(selectedGenomeVersion == 'hg19'){
    hg <- hg19
  }
  if(selectedGenomeVersion == 'hg38'){
    hg <- hg38
  }
  
  if (selectedGenomeVersion == 'Custom') {
    if (is.null(custom_seq) || custom_seq == "") {
      stop("Custom sequence is required when 'Custom' genome is selected.")
    }
    custom_seq_clean <- gsub("[^a-zA-Z]", "", custom_seq)
    seq_len <- nchar(custom_seq_clean)
    start_pos <- max(1, s1)
    end_pos <- min(seq_len, e1)
    
    if (start_pos > end_pos) {
       stop("Calculated coordinates are outside the custom sequence range.")
    }
    seq_result <- substr(custom_seq_clean, start_pos, end_pos)
  } else {
    genome_loc <- paste0(locchr,":",as.character(s1),"-",as.character(e1))
    
    # 使用缓存版本的序列获取
    seq_result <- get_sequence_cached(hg, genome_loc)
  }
  seqs <- DNAString(seq_result)

  seqs.path <- seqs_path
  # seqs.path <- dirname(seqs.fasta.file)
  #* 设置Primer3参数文件 p3_settings_file
  p3.paras <- list(
    PRIMER_TASK="generic",
    PRIMER_NUM_RETURN=1000,
    PRIMER_DNA_CONC=500,
    PRIMER_DNTP_CONC=0.8,
    PRIMER_SALT_CORRECTIONS=1,
    PRIMER_SALT_MONOVALENT=50,
    PRIMER_SALT_DIVALENT=1.5,
    PRIMER_MAX_END_GC=2,
    PRIMER_MIN_GC=primer_min_gc,
    PRIMER_MAX_GC=primer_max_gc,
    PRIMER_OPT_GC_PERCENT=primer_opt_gc,
    PRIMER_MIN_TM=primer_min_tm,
    PRIMER_MAX_TM=primer_max_tm,
    PRIMER_OPT_TM=primer_opt_tm,
    PRIMER_MIN_SIZE=primer_min_size,
    PRIMER_MAX_SIZE=primer_max_size,
    PRIMER_OPT_SIZE=primer_opt_size,
    PRIMER_PAIR_MAX_DIFF_TM=2,
    PRIMER_PRODUCT_SIZE_RANGE=paste0(product_min_size, "-", product_max_size),
    PRIMER_PICK_ANYWAY=1
  )
  p3.paras <- paste(names(p3.paras), '=', p3.paras, sep='')
  p3.paras <- c(
    "Primer3 File - http://primer3.sourceforge.net",
    "P3_FILE_TYPE=settings",
    "",
    p3.paras,
    "="
  )
  p3.settings.file <- file.path(seqs.path, "p3.settings.file")
  writeLines(p3.paras, p3.settings.file)
  tmpfiles <- p3.settings.file
  #* 设置Primer3参数文件 input_file
  seq.ids <- paste("SEQUENCE_ID=", inputloc, sep='')
  seq.templates <- paste("SEQUENCE_TEMPLATE=", as.character(seqs), sep='')
  content.input.file <- paste(seq.ids, seq.templates, '=', sep="\n")
  input.file <- file.path(seqs.path, "p3.input.file")
  writeLines(content.input.file, input.file)
  tmpfiles <- c(tmpfiles, input.file)
  #* 运行 Primer3 获取引物
  output.file <- file.path(seqs.path, "p3.temp1")
  p3.settings <- paste("-p3_settings_file=", p3.settings.file, sep='')
  p3.output <- paste("-output=", output.file, sep='')
  cmd <- paste("primer3_core", p3.settings, p3.output, input.file,"-default_version=1")
  system(cmd)
  tmpfiles <- c(tmpfiles, output.file)
  #* 解析 Primer3 输出文件
  #
  # p3.results <- readLines("/Users/fangy/Desktop/sanger/p3.temp2.temp1")
  p3.results <- readLines(output.file)
  group.start <- grep("SEQUENCE_ID", p3.results)
  group.end <- c(group.start[-1]-1, length(p3.results))
  # seq.ids <- names(seqs)
  for(i in 1:length(seq.ids)){
    sel <- group.start[i]:group.end[i]
    p3.results[sel] <- paste(seq.ids[i], p3.results[sel], sep="_")
  }
  writeLines(p3.results, output.file)
  primers.seq <- p3.results[grep("(LEFT|RIGHT)_[0-9]+_SEQUENCE=", p3.results)]
  primers.name <- gsub("(.+)_PRIMER(_[^=]+)_SEQUENCE.*", "\\1\\2", primers.seq)
  primers.name <- gsub("LEFT", "L", primers.name)
  primers.name <- gsub("RIGHT", "R", primers.name)
  primers.name <- primers.name[1:20]
  primers.seq <- gsub(".+=(.+)", "\\1", primers.seq)
  primers.seq <- primers.seq[1:20]
  primers.tm <- p3.results[grep("(LEFT|RIGHT)_[0-9]+_TM=", p3.results)]
  primers.tm <- gsub(".+=(.+)", "\\1", primers.tm)
  primers.tm <- primers.tm[1:20]
  primers.start <- p3.results[grep("(LEFT|RIGHT)_[0-9]=", p3.results)]
  primers.start <- gsub(".+=(.+)", "\\1", primers.start)
  primers.start <- sapply(strsplit(primers.start, ","), function(x) x[1])
  primers.start <- primers.start[1:20]
  primers.gc <- p3.results[grep("(LEFT|RIGHT)_[0-9]+_GC_PERCENT=", p3.results)]
  primers.gc <- gsub(".+=(.+)", "\\1", primers.gc)
  primers.gc <- primers.gc[1:20]
  primers.any <- p3.results[grep("(LEFT|RIGHT)_[0-9]+_SELF_ANY=", p3.results)]
  primers.any <- gsub(".+=(.+)", "\\1", primers.any)
  primers.any <- primers.any[1:20]
  primers.end <- p3.results[grep("(LEFT|RIGHT)_[0-9]+_SELF_END=", p3.results)]
  primers.end <- gsub(".+=(.+)", "\\1", primers.end)
  primers.end <- primers.end[1:20]
 
  primers.len <- p3.results[grep("_[0-9]+_PRODUCT_SIZE=", p3.results)]
  primers.len <- gsub(".+=(.+)", "\\1", primers.len)
  primers.len <- rep(primers.len, each =2)
  primers.len <- primers.len[1:20]
  
  primers.row <- rep(c("Forward primer","Reverse primer"),6)
  primers.row <- primers.row[1:20]
  
  primers <- data.frame( "Forward/Reverse"= primers.row,
                         Name=primers.name,
                         "Sequence (5'->3')"=primers.seq, 
                        Start=primers.start,
                        TM=primers.tm,
                        "%GC" =primers.gc,
                        "Self complementarity" =primers.any,
                        "Self 3' complementarity"=primers.end,
                        "Product length" = primers.len,
                        check.names = FALSE
                       )
  
  primers <- primers[1:20,]
  #* BLAST 分析，输出便于程序解析的 m8 格式
  if (selectedGenomeVersion != 'Custom') {
      blast.in <- file.path(seqs.path, "blast.in")
      xxx <- paste(">", primers.name, sep='')
      xxx <- paste(xxx, primers.seq, sep="\n")
      writeLines(xxx, blast.in)
      rm(xxx)
      tmpfiles <- c(tmpfiles, blast.in)
      # blast.db <- tk_choose.files(default = "~//", multi = FALSE, caption = "BLAST database")
      # blast.db <- sub("(.+)\\.[^\\.]+", "\\1", blast.db)
      blast.db <- bdb
      blast.out <- file.path(seqs.path, "blast.out")
      cmd <- paste(bn," -evalue 1e-1 -outfmt 6 -num_threads 6 -task blastn-short -query",
                   blast.in, "-out", blast.out, "-db", blast.db)
      # cmd <- paste("blastall -p blastn -e 10000 -F F -m 8 -a 4 -i",
      #              blast.in, "-o", blast.out, "-d", blast.db)
      system(cmd)
      tmpfiles <- c(tmpfiles, blast.out)
      #* 解析 BALST 输出结果。 m8 结果共有12列，分别为：
      # 1. Query id
      # 2. Subject id
      # 3. % identity
      # 4. alignment length
      # 5. mismatches
      # 6. gap openings
      # 7. q.start
      # 8. q.end
      # 9. s.start
      # 10. s.end
      # 11. e-value
      # 12. bit score
      # 这里我们仅要求 q.start=1, q.end=引物长度 的比对结果有且仅有一个，即目标序列的匹配
      if (file.exists(blast.out) && file.size(blast.out) > 0) {
        blast.result <- read.table(blast.out, stringsAsFactors = FALSE)[,c(1,7,8)]
        sel <- blast.result[,2]==1
        blast.result <- blast.result[sel,]
        primers.n <- length(primers.name)
        sel <- rep(FALSE, primers.n)
        for(i in 1:primers.n){
            sel.sub <- blast.result[,1]==primers.name[i]
            blast.sub <- blast.result[sel.sub,3]
            if (length(blast.sub) > 0) {
                max.qend <- max(blast.sub)
                blast.sub <- blast.sub[blast.sub==max.qend]
                if(length(blast.sub)==1 & max.qend==nchar(primers.seq[i]))
                sel[i] <- TRUE
            }
        }
        primers <- primers[sel,]
      } else {
        primers <- primers[FALSE,]
      }
  }
  primers$Name <- sub("SEQUENCE_ID=([^_]+)_.*", "\\1", primers$Name)
  # 修改列名
  names(primers)[names(primers) == "Name"] <- "Genome location"
  #* 删除临时文件，输出结果
  # file.remove(tmpfiles)
  result.file <- file.path(seqs.path, "primer3.results.csv")
  write.csv(primers, result.file, quote=FALSE, row.names=FALSE)
  return(primers)


}
primerRes2 <- function(inputloc){
  tem <- data.frame(A=c(1,2,3),B=c("A","B","C"))
  return(tem)
}

# # 自定义直线的开始和结尾
# start <- 0
# end <- 10
# # 自定义刻度的间隔
# interval <- 2
# # 创建一个空白的绘图窗口，只有x轴刻度
# plot(0, type = "n", xlim = c(start, end), ylim = c(0, 1), xlab = "X", ylab = "", yaxt = "n",)
# # 绘制直线
# lines(c(start, end), c(0.5, 0.5))
# # 添加x轴刻度
# axis(side = 1, at = seq(start, end, interval))
# # 去掉外边框
# box(bty = "n")
# 
# 
# # primerRes <- function(inputloc){
# start <- 0
# end <- 10
# datax<-1:10
# datay <- rep(c(1,2),5)
# plot(0,
#      xlim = c(start, end),
#      axes=FALSE,
#      pch=0,
#      col="red",ylim=c(-1,3))
# plot(Na,
#      xlim = c(start, end),
#      axes=FALSE)
# 
# axis(3,at=1:10)
# # }
#
loc_axis <- function(x=1){
  plot(runif(2), runif(2),
       xlim=c(189586, 189596), ylim=c(0,1000),
       axes=FALSE, #Don't plot the axis
       type="n",  #hide the points
       ylab="", xlab="") #No axis labels
  axis(1, seq(189586, 189596, 1))
}



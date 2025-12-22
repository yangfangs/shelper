
# Package environment to store configuration
.shelper_env <- new.env()
.shelper_env$hg19 <- "/Users/fangy/fsdownload/blastdb/ucsc.hg19.fasta"
.shelper_env$hg38 <- "/Users/fangy/fsdownload/hg38/hg38.fa"
.shelper_env$blastn <- "/Users/fangy/Downloads/ncbi-blast-2.15.0+/bin/blastn"
.shelper_env$blast_db <- "/Users/fangy/fsdownload/blastdb/hg19"
.shelper_env$temp_path <- tempdir()

#' Set Shelper Configuration
#' 
#' @param hg19 Path to hg19 fasta
#' @param hg38 Path to hg38 fasta
#' @param blastn Path to blastn executable
#' @param blast_db Path to blast database
#' @param temp_path Path for temporary files
#' @export
set_shelper_config <- function(hg19 = NULL, hg38 = NULL, blastn = NULL, blast_db = NULL, temp_path = NULL) {
  if (!is.null(hg19)) .shelper_env$hg19 <- hg19
  if (!is.null(hg38)) .shelper_env$hg38 <- hg38
  if (!is.null(blastn)) .shelper_env$blastn <- blastn
  if (!is.null(blast_db)) .shelper_env$blast_db <- blast_db
  if (!is.null(temp_path)) .shelper_env$temp_path <- temp_path
}

# Cache setup
# We'll use a local cache for the package session
# Ideally this should be initialized in .onLoad or similar, but for simplicity:
.shelper_cache <- cachem::cache_mem(max_size = 512 * 1024^2, max_age = 3600)

#' Get Sequence Cached
#' @export
get_sequence_cached <- memoise::memoise(function(hg, genome_loc) {
  cmd_samtools <- paste0("samtools faidx ", hg, " ", genome_loc)
  # Check if samtools is available
  if (Sys.which("samtools") == "") {
    warning("samtools not found in PATH. Make sure samtools is installed.")
  }
  result <- system(cmd_samtools, intern = TRUE)
  if (length(result) > 1) {
    return(paste(result[2:length(result)], collapse = ""))
  }
  return(result[2])
}, cache = .shelper_cache)

#' Calculate Figure Height
#' @export
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

#' Wrap Fixed Width String
#' @export
wrap_fixed = function(string, width=80) {
  pattern = paste("(.{1,", width, "})", sep="")
  res = gsub(pattern, "\\1\n", string)
  return(res)
}

#' Clean String
#' @export
cleanstring <- function(string) {
  string <- gsub("^>.*?\n", "", string)
  string <- toupper(string)
  string <- gsub("[^ACGTRYSWKMBDHVN]", "", string, perl=TRUE)
  return(string)
}

#' Align Chromatogram
#' @export
alignchromatogram <- function(data, block.width=50, trim=FALSE, refseq, trim5, trim3) {
  if (is.null(data)) return(NULL)
  d <- sangerseqR::setAllelePhase(data, refseq, trim5, trim3)
  refseq <- toString(d@primarySeq)
  altseq <- toString(d@secondarySeq)
  if (trim == TRUE) {
    altseq <- toString(d@secondarySeq[(trim5 + 1):(nchar(altseq) - trim3)])
  }
  names(altseq) <- "Alt Allele"
  names(refseq) <- d@primarySeqID
  if (trim == TRUE) {
    pa <- Biostrings::pairwiseAlignment(altseq, refseq, type="local", gapExtension=-2)
  } else {
    pa <- Biostrings::pairwiseAlignment(altseq, refseq, type="global", gapExtension=-2)
  }
  alignment <- paste(capture.output(Biostrings::writePairwiseAlignments(pa, block.width=block.width)), collapse="\n")
  results <- list(altseq=altseq, refseq=gsub("-", "", refseq), alignment=alignment)
  return(results)
}

#' Pick Chromatogram
#' @export
pickchromatogram <- function(obj, trim5=0, trim3=0, width=100, pixelsperrow=200, showtrim=FALSE,inputloc,selectedGenomeVersion){
  extend_len <- 10
  if(selectedGenomeVersion == 'hg19'){
    hg <- .shelper_env$hg19
  }
  if(selectedGenomeVersion == 'hg38'){
    hg <- .shelper_env$hg38
  }
  
  # 坐标提取序列
  inputloc_list <- strsplit(inputloc,"-")[[1]]
  locchr <- paste0('chr', inputloc_list[1])
  locid <- as.numeric(inputloc_list[2])

  s1<- locid - extend_len
  e1 <- locid + extend_len + max(nchar(inputloc_list[[3]]),nchar(inputloc_list[[4]])) * 2

  genome_loc <- paste0(locchr,":",as.character(s1),"-",as.character(e1))
  
  # 使用缓存版本的序列获取
  seq_result <- get_sequence_cached(hg, genome_loc)
  extracted_sequence <- Biostrings::DNAString(seq_result)
  
  #读取sanger测序结果
  seq <- sangerseqR::makeBaseCalls(sangerseqR::readsangerseq(obj), ratio = 0.33)

  seqstr <- sangerseqR::primarySeq(seq,string=TRUE)

  pa <- Biostrings::pairwiseAlignment(extracted_sequence,seqstr,
                          type = "global-local",gapOpening=20, gapExtension=2)
  # writePairwiseAlignments(pa)
  
  start_value <- Biostrings::start(pa@subject@indel[[1]])
  width_value <- Biostrings::width(pa@subject@indel[[1]])
  
  #验证纯合缺失
  check_del_het <- function(start_value,width_value){
    
    if (start_value == extend_len + 2 & width_value == nchar(inputloc_list[[3]]) -1 ){
      # print("Del, Heterozygous")
      return("Homozygous")
    }else{
      return("Heterozygous")
    }
  }

  subject_start <- pa@subject@range@start
  subject_end <- pa@subject@range@start + pa@subject@range@width - 1
  
  subject_loc <- pa@subject@range@start + extend_len
  
  subject_str <- Biostrings::subseq(sangerseqR::primarySeq(seq, string = TRUE), start = subject_start, width = pa@subject@range@width)
  
  subject_str_secondary <- Biostrings::subseq(sangerseqR::secondarySeq(seq, string = TRUE), start = subject_start, width = pa@subject@range@width)
  
  check_start_primary <- Biostrings::subseq(sangerseqR::primarySeq(seq, string = TRUE), start = subject_loc, width = extend_len+1)
  
  check_start_secondary <- Biostrings::subseq(sangerseqR::secondarySeq(seq, string = TRUE), start = subject_loc, width = extend_len+1)
  
  check_ref <- Biostrings::subseq(extracted_sequence, start = extend_len+1, end =nchar(extracted_sequence))
  
  check_del <- function(subject_str,subject_str_secondary,extracted_sequence){
    ref_str <- as.character(extracted_sequence)
    del_loc <- extend_len + 2
    ref_end_loc <- del_loc + 5
    backward_vec <- vector()
    for (backward in del_loc:ref_end_loc) {
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
    # print(backward_vec)
    if(sum(backward_vec) == length(del_loc:ref_end_loc)){
      # print("Del")
      return("Del")
    }else{
      # print("Not Del")
      return("Not Del")
    }
  }
  
  check_del_multiple <- function(subject_str,subject_str_secondary,extracted_sequence){
    ref_str <- as.character(extracted_sequence)
    del_loc <- extend_len + 2
    del_len <- nchar(inputloc_list[[3]]) -1
    ref_end_loc <- nchar(extracted_sequence)-1
    backward_vec <- vector()
    del_loc_end <- del_loc+del_len
    
    for (backward in del_loc: del_loc_end) {
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
    # print(backward_vec)
    if(sum(backward_vec) == length(del_loc:del_loc_end)){
      # print("Del")
      return("Del")
    }else{
      # print("Not Del")
      return("Not Del")
    }
  }
  
  check_ins <- function(subject_str,subject_str_secondary,extracted_sequence){
    ref_str <- as.character(extracted_sequence)
    ins_loc <- extend_len + 2
    ref_end_loc <- ins_loc + 5
    backward_vec <- vector()
    for (backward in ins_loc:ref_end_loc) {
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
      # print("Ins")
      return("Ins")
    }else{
      # print("Not Ins")
      return("Not Ins")
    }
  }
  
  check_ins_multiple <- function(subject_str,subject_str_secondary,extracted_sequence){
    ref_str <- as.character(extracted_sequence)
    ins_loc <- extend_len + 1
    ins_len <- nchar(inputloc_list[[3]])
    ref_end_loc <- nchar(extracted_sequence)-1
     
    ins_loc_end <- ins_loc + ins_len
    backward_vec <- vector()
    for (backward in ins_loc_end:ref_end_loc) {
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
      ref_wark <- substr(ref_str,backward,backward)
      ref_wark2 <- substr(ref_str,backward - ins_len,backward - ins_len)
      refref_wark <- paste0(ref_wark,ref_wark2)
      forward_ref <- substr(ref_str,backward-4 ,backward-4 )
      check_dup_vec <- c(check_dup_vec,grepl(forward_ref, refref_wark))
    }
    
    if(sum(backward_vec) == length(ins_loc_end:ref_end_loc)){
      if(sum(check_dup_vec) == ins_len){
        # print("Dup")
        return("Dup")
      } else{
        # print("Ins")
        return("Ins")
        }
    }else{
      # print("Not Ins")
      return("Not Ins")
    }
  }
  
  check_ins_het <- function(subject_str,subject_str_secondary){
    ins_loc <- extend_len + 2
    if (substr(subject_str,ins_loc,nchar(subject_str)) ==  substr(subject_str_secondary,ins_loc,nchar(subject_str_secondary))){
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
  
  extracted_sequence_str <-as.character(extracted_sequence)
  
  #截取关注位点前后10bp图像用于判断
  seq2 <- seq
  seq2@primarySeq <- Biostrings::DNAString(subject_str)
  seq2@secondarySeq <- Biostrings::DNAString(subject_str_secondary)
  seq2@peakPosMatrix <-  seq2@peakPosMatrix[subject_start:subject_end,]
  
  alignment <- paste(capture.output(Biostrings::writePairwiseAlignments(pa, block.width=50)), collapse="\n")
  
  # Initialize loc_type and loc_peaks
  loc_type <- "Unknown"
  loc_peaks <- "Unknown"

  #单碱基缺失
  if(nchar(inputloc_list[3]) > nchar(inputloc_list[4]) & nchar(inputloc_list[3]) == 2 ){
    loc_type <- check_del(subject_str,subject_str_secondary,extracted_sequence)

    if (length(start_value) == 0) {
      loc_peaks <- "Heterozygous"
    } else {
      loc_peaks <- check_del_het(start_value,width_value)
    }
  }
  #长缺失
  if(nchar(inputloc_list[3]) > nchar(inputloc_list[4]) & nchar(inputloc_list[3]) > 2 ){
    loc_type <- check_del_multiple(subject_str,subject_str_secondary,extracted_sequence)
    
    if (length(start_value) == 0) {
      loc_peaks <- "Heterozygous"
    } else {
      loc_peaks <- check_del_het(start_value,width_value)
    }
  }
  #单碱基ins
  if(nchar(inputloc_list[3]) < nchar(inputloc_list[4]) & nchar(inputloc_list[4]) == 2 ){
    loc_type <- check_ins(subject_str,subject_str_secondary,extracted_sequence)
    loc_peaks <- check_ins_het(subject_str,subject_str_secondary)
  }
  
  #长碱基ins
  if(nchar(inputloc_list[3]) < nchar(inputloc_list[4])  & nchar(inputloc_list[4]) > 2 ){
    loc_type <- check_ins_multiple(subject_str,subject_str_secondary,extracted_sequence)
    loc_peaks <- check_ins_het(subject_str,subject_str_secondary)
  }
  if(nchar(inputloc_list[3]) == nchar(inputloc_list[4])){
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

#' Pick Chromatogram 2 (Parallel)
#' @export
pickchromatogram2 <- function(inputinfo,selectedGenomeVersion){
  check_loc <- function(inputinfo,mu){
          perform_operation <- function(row) {
          inputloc_list <- strsplit(row["GenomicCoordinate"],"-")[[1]]
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
  
  # 使用 future_lapply 实现并行处理，显著提升多文件处理速度
  # Note: ensure plan() is called before this function if parallelism is desired
  results <- future.apply::future_lapply(1:nrow(inputinfo), function(i) {
    obj <- inputinfo$FilePath[i] 
    inputloc <- inputinfo$GenomicCoordinate[i]
    genome_ver <- selectedGenomeVersion
    # 调用pickchromatogram函数，并传入其他需要的参数
    pickchromatogram(obj, trim5=0, trim3=0, width=100, pixelsperrow=200, showtrim=FALSE, inputloc=inputloc,selectedGenomeVersion=genome_ver)
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
}

#' Generate Sample Data
#' @export
generate <- function(input) {
  data.frame(
    GenomeLocation = c("Location1", "Location2", "Location3"),
    Genotype = c("TypeA", "TypeB", "TypeC"),
    Start = c(100, 200, 300),
    End = c(150, 250, 350),
    Check = c(TRUE, FALSE, TRUE),
    Figure = seq_len(3)
  )
}

#' Primer Result
#' @export
primerRes <- function(inputloc,selectedGenomeVersion){

  inputloc_list <- strsplit(inputloc,"-")[[1]]
  locchr <- paste0('chr', inputloc_list[1])
  locid <- as.numeric(inputloc_list[2])

  s1<- locid - 500
  e1 <- locid + 500
  
  hg <- .shelper_env$hg19
  if(selectedGenomeVersion == 'hg19'){
    hg <- .shelper_env$hg19
  }
  if(selectedGenomeVersion == 'hg38'){
    hg <- .shelper_env$hg38
  }
  genome_loc <- paste0(locchr,":",as.character(s1),"-",as.character(e1))
  
  # 使用缓存版本的序列获取
  seq_result <- get_sequence_cached(hg, genome_loc)
  seqs <- Biostrings::DNAString(seq_result)

  # Use temp path from config
  seqs.path <- .shelper_env$temp_path
  if (!dir.exists(seqs.path)) dir.create(seqs.path, recursive = TRUE)

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
    PRIMER_MIN_GC=35,
    PRIMER_MAX_GC=65,
    PRIMER_MIN_TM=58,
    PRIMER_MAX_TM=64,
    PRIMER_OPT_TM=60,
    PRIMER_MIN_SIZE=18,
    PRIMER_MAX_SIZE=26,
    PRIMER_OPT_SIZE=20,
    PRIMER_PAIR_MAX_DIFF_TM=2,
    PRIMER_PRODUCT_SIZE_RANGE="300-1000",
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
  # Check if primer3_core is installed
  if (Sys.which("primer3_core") == "") {
      warning("primer3_core not found in PATH.")
  }
  system(cmd)
  tmpfiles <- c(tmpfiles, output.file)
  #* 解析 Primer3 输出文件
  if (!file.exists(output.file)) {
      stop("Primer3 output file not created. Check if primer3_core is installed and working.")
  }
  p3.results <- readLines(output.file)
  group.start <- grep("SEQUENCE_ID", p3.results)
  group.end <- c(group.start[-1]-1, length(p3.results))
  # seq.ids <- names(seqs)
  # Note: seq.ids here is a single string "SEQUENCE_ID=...", but the loop implies multiple?
  # The original code loop seems to iterate over seq.ids, but here seq.ids is length 1.
  # Assuming inputloc is single location.
  
  # for(i in 1:length(seq.ids)){ # This was weird in original code if seq.ids is not vector
  #   sel <- group.start[i]:group.end[i]
  #   p3.results[sel] <- paste(seq.ids[i], p3.results[sel], sep="_")
  # }
  # Simplification for single input:
  # p3.results is already for one sequence.
  
  # Replicating original logic roughly but safely
  # In original code: seq.ids <- paste("SEQUENCE_ID=", inputloc, sep='') (length 1)
  # for(i in 1:length(seq.ids)) -> i=1
  # sel <- group.start[1]:group.end[1] -> all lines
  # p3.results[sel] <- paste(seq.ids[1], p3.results[sel], sep="_")
  # This prepends SEQUENCE_ID=..._ to every line.
  
  # Let's keep it as is
  for(i in 1:1){
    if (length(group.start) > 0) {
        sel <- group.start[i]:group.end[i]
        p3.results[sel] <- paste(seq.ids[i], p3.results[sel], sep="_")
    }
  }

  writeLines(p3.results, output.file)
  primers.seq <- p3.results[grep("(LEFT|RIGHT)_[0-9]+_SEQUENCE=", p3.results)]
  primers.name <- gsub("(.+)_PRIMER(_[^=]+)_SEQUENCE.*", "\\1\\2", primers.seq)
  primers.name <- gsub("LEFT", "L", primers.name)
  primers.name <- gsub("RIGHT", "R", primers.name)
  
  limit <- min(20, length(primers.name))
  primers.name <- primers.name[1:limit]
  primers.seq <- gsub(".+=(.+)", "\\1", primers.seq)
  primers.seq <- primers.seq[1:limit]
  primers.tm <- p3.results[grep("(LEFT|RIGHT)_[0-9]+_TM=", p3.results)]
  primers.tm <- gsub(".+=(.+)", "\\1", primers.tm)
  primers.tm <- primers.tm[1:limit]
  primers.start <- p3.results[grep("(LEFT|RIGHT)_[0-9]=", p3.results)]
  primers.start <- gsub(".+=(.+)", "\\1", primers.start)
  primers.start <- sapply(strsplit(primers.start, ","), function(x) x[1])
  primers.start <- primers.start[1:limit]
  primers.gc <- p3.results[grep("(LEFT|RIGHT)_[0-9]+_GC_PERCENT=", p3.results)]
  primers.gc <- gsub(".+=(.+)", "\\1", primers.gc)
  primers.gc <- primers.gc[1:limit]
  primers.any <- p3.results[grep("(LEFT|RIGHT)_[0-9]+_SELF_ANY=", p3.results)]
  primers.any <- gsub(".+=(.+)", "\\1", primers.any)
  primers.any <- primers.any[1:limit]
  primers.end <- p3.results[grep("(LEFT|RIGHT)_[0-9]+_SELF_END=", p3.results)]
  primers.end <- gsub(".+=(.+)", "\\1", primers.end)
  primers.end <- primers.end[1:limit]
 
  primers.len <- p3.results[grep("_[0-9]+_PRODUCT_SIZE=", p3.results)]
  primers.len <- gsub(".+=(.+)", "\\1", primers.len)
  primers.len <- rep(primers.len, each =2)
  primers.len <- primers.len[1:limit]
  
  primers.row <- rep(c("Forward primer","Reverse primer"), ceiling(limit/2))
  primers.row <- primers.row[1:limit]
  
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
  
  # primers <- primers[1:20,] # handled by limit
  
  #* BLAST 分析，输出便于程序解析的 m8 格式
  blast.in <- file.path(seqs.path, "blast.in")
  xxx <- paste(">", primers.name, sep='')
  xxx <- paste(xxx, primers.seq, sep="\n")
  writeLines(xxx, blast.in)
  rm(xxx)
  tmpfiles <- c(tmpfiles, blast.in)
  
  blast.db <- .shelper_env$blast_db
  blast.out <- file.path(seqs.path, "blast.out")
  
  bn <- .shelper_env$blastn
  
  cmd <- paste(bn," -evalue 1e-1 -outfmt 6 -num_threads 6 -task blastn-short -query",
               blast.in, "-out", blast.out, "-db", blast.db)
  
  if (Sys.which("blastn") == "" && !file.exists(bn)) {
      warning("blastn not found.")
  }
  
  system(cmd)
  tmpfiles <- c(tmpfiles, blast.out)
  
  #* 解析 BALST 输出结果。
  if (!file.exists(blast.out)) {
      # Return primers without blast filtering if blast fails
      return(primers)
  }
  
  blast.result <- tryCatch({
      read.table(blast.out, stringsAsFactors = FALSE)[,c(1,7,8)]
  }, error = function(e) NULL)
  
  if (!is.null(blast.result) && nrow(blast.result) > 0) {
      sel <- blast.result[,2]==1
      blast.result <- blast.result[sel,]
      primers.n <- length(primers.name)
      sel <- rep(FALSE, primers.n)
      for(i in 1:primers.n){
        sel.sub <- blast.result[,1]==primers.name[i]
        blast.sub <- blast.result[sel.sub,3]
        if(length(blast.sub) > 0) {
            max.qend <- max(blast.sub)
            blast.sub <- blast.sub[blast.sub==max.qend]
            if(length(blast.sub)==1 & max.qend==nchar(primers.seq[i]))
              sel[i] <- TRUE
        }
      }
      primers <- primers[sel,]
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

#' Primer Result 2 (Dummy)
#' @export
primerRes2 <- function(inputloc){
  tem <- data.frame(A=c(1,2,3),B=c("A","B","C"))
  return(tem)
}

#' Local Axis Plot
#' @export
loc_axis <- function(x=1){
  plot(runif(2), runif(2),
       xlim=c(189586, 189596), ylim=c(0,1000),
       axes=FALSE, #Don't plot the axis
       type="n",  #hide the points
       ylab="", xlab="") #No axis labels
  axis(1, seq(189586, 189596, 1))
}

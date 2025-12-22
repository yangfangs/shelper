# ============================================================
# SangerPriVar 配置文件
# 根据服务器环境调整以下参数以获得最佳性能
# ============================================================

# ----------------------------------------------------------
# 并行处理配置
# ----------------------------------------------------------

# 并行工作进程数量
# 建议值: CPU核心数 - 2 (保留系统资源)
# 设置为 NULL 时自动检测
PARALLEL_WORKERS <- NULL

# 并行策略: "multisession" (推荐用于服务器) 或 "multicore" (仅Linux)
PARALLEL_PLAN <- "multisession"

# ----------------------------------------------------------
# 缓存配置
# ----------------------------------------------------------

# 内存缓存最大大小 (字节)
# 默认 512MB
CACHE_MAX_SIZE <- 512 * 1024^2

# 缓存过期时间 (秒)
# 默认 1小时
CACHE_MAX_AGE <- 3600

# 是否启用序列缓存
ENABLE_SEQUENCE_CACHE <- TRUE

# ----------------------------------------------------------
# 外部工具路径配置
# ----------------------------------------------------------

# hg19 参考基因组路径
HG19_PATH <- "/Users/fangy/fsdownload/blastdb/ucsc.hg19.fasta"

# hg38 参考基因组路径
HG38_PATH <- "/Users/fangy/fsdownload/hg38/hg38.fa"

# BLASTN 工具路径
BLASTN_PATH <- "/Users/fangy/Downloads/ncbi-blast-2.15.0+/bin/blastn"

# BLAST 数据库路径
BLAST_DB_PATH <- "/Users/fangy/fsdownload/blastdb/hg19"

# 临时文件路径
TEMP_PATH <- "/Users/fangy/Desktop/sanger"

# ----------------------------------------------------------
# Shiny 服务器配置
# ----------------------------------------------------------

# 最大上传文件大小 (MB)
MAX_UPLOAD_SIZE_MB <- 100

# 会话超时时间 (秒)
SESSION_TIMEOUT <- 600

# ----------------------------------------------------------
# 日志配置
# ----------------------------------------------------------

# 是否启用详细日志
VERBOSE_LOGGING <- FALSE

# 日志文件路径 (NULL 表示不写入文件)
LOG_FILE_PATH <- NULL


# ============================================================
# 加载配置的辅助函数
# ============================================================

load_config <- function() {
  # 设置上传限制
  options(shiny.maxRequestSize = MAX_UPLOAD_SIZE_MB * 1024^2)
  
  # 设置并行workers
  if (is.null(PARALLEL_WORKERS)) {
    n_workers <- max(1, parallel::detectCores() - 2)
  } else {
    n_workers <- PARALLEL_WORKERS
  }
  
  # 配置future计划
  if (PARALLEL_PLAN == "multicore" && .Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = n_workers)
  } else {
    future::plan(future::multisession, workers = n_workers)
  }
  
  message(sprintf("SangerPriVar 配置已加载:"))
  message(sprintf("  - 并行工作进程: %d", n_workers))
  message(sprintf("  - 缓存大小: %.0f MB", CACHE_MAX_SIZE / 1024^2))
  message(sprintf("  - 最大上传: %.0f MB", MAX_UPLOAD_SIZE_MB))
}

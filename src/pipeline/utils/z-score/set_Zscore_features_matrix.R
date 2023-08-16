usage <- "
USAGE:

Rscript this_script.R --viral=path/2/file --nonviral=path/2/file --nohit=path/2/file --dir=path/2/dir

- Output directory (--dir) is required;
- All other arguments are optional but at least one must be provided;
"

# =============================================================
# -- Auxiliary functions declarations -------------------------
# =============================================================

bind_df <- function(matrix, path, similarity_label) {
	
	cat('\n-- bind_df --\n')
	cat('similarity_label:', similarity_label, "\n")
	cat('matrix.shape [prev]', dim(matrix), "\n")
	
	df <- read.delim(path, sep="#", header=F, row.names=NULL, quote=NULL)
	cat('df.shape [1]', dim(df), "\n")
	
	df <- t(as.matrix(df))
	cat('df.shape [2]', dim(df), "\n")
	
	# df <- as.data.frame(df[-1,])
	df <- as.data.frame(df)
	cat('df.shape [3]', dim(df), "\n")

	df$similarity_label <- similarity_label
	cat('df.shape [4]', dim(df), "\n")
	
	matrix <- rbind(matrix, df)
	cat('matrix.shape [next]', dim(matrix), "\n")
	
	return(matrix)
}

validate_arg_path <- function(arg, path, is_dir) {
	if (!is.character(path)
		|| (is_dir && !dir.exists(path)) # Test valid directory
		|| (!is_dir && (!file.exists(path) || !grepl("\\.tab$", path))) # Test valid tab file
	) {
		path_type <- ifelse(is_dir, "directory", "file")
		stop(paste("ERROR: '", arg,"' must be a valid ", path_type, ". (Given value: '", path, "')\n", usage, sep=""))
	}
}

# =============================================================
# -- Parse input args -----------------------------------------
# =============================================================

args <- commandArgs(trailingOnly=TRUE)

matrix <- data.frame()
has_input <- FALSE
output_dir <- ""

for (arg in args) {
  
    parts <- strsplit(arg, "=")[[1]]
    if (length(parts) != 2) {
        next
    }

	arg_name <- parts[1]
	arg_value <- parts[2]

	# Parse output directory
	if (arg_name == "--dir") {
		validate_arg_path(arg_name, arg_value, TRUE)
		output_dir <- arg_value
		next
	}

	# Parse input files
    similarity_type_args <- c("--viral", "--nonviral", "--nohit")

    if (arg_name %in% similarity_type_args) {
		validate_arg_path(arg_name, arg_value, FALSE)
        similarity_label <- gsub("(--)(.+)", "\\2", arg_name)
        matrix <- bind_df(matrix, arg_value, similarity_label)
        has_input <- TRUE
    }
}

# Validate minimum conditions to forward
if (!has_input) {
	stop(paste("ERROR: There ain't no input data to be processed\n", usage))
}

has_out_dir <- nchar(trimws(output_dir)) > 0
if (!has_out_dir) {
	stop(paste("ERROR: Argument '--dir' is required.\n", usage))
}

# =============================================================
# -- Main -----------------------------------------------------
# =============================================================

cat('\n-- main --\n')
cat('matrix.shape [1]', dim(matrix), "\n")

# Parse numeric values
matrix[,2:50]<- as.data.frame(lapply(matrix[,2:50], function(x) {
	as.numeric(as.character(x))
}))
cat('matrix.shape [2]', dim(matrix), "\n")


# Normalize numeric values
matrix[,44:46] <- matrix[,44:46] / matrix[,ncol(matrix)-1]
matrix[,ncol(matrix) - 2] <- matrix[,ncol(matrix) - 2] / matrix[,ncol(matrix)-1]
cat('matrix.shape [3]', dim(matrix), "\n")

# Add column names
n<-seq(15,35,by=1)
n2<- seq(-15,-35,by=-1)
nd <- c("dens15to18","dens20to22","dens25to29","ratiosi_pi","ratio_si","dens18to35","length")

n3<-c("Contigs_ID", c(n, n2, nd), "Similarity_label")
colnames(matrix) <- n3
cat('matrix.shape [4]', dim(matrix), "\n")

# Who knows...
matrix[is.na(matrix)] <- 0
matrix <- matrix[rowSums(matrix[2:50]) > 0,]
cat('matrix.shape [5]', dim(matrix), "\n")

matrix[, 44:49][matrix[, 44:49] == 0] <- matrix[, 44:49][matrix[, 44:49] == 0] + 0.00001
matrix[, 44:49] <- log2(matrix[,44:49])
matrix <- matrix[, c(1, ncol(matrix), 2:(ncol(matrix) - 1))]
cat('matrix.shape [6]', dim(matrix), "\n")

# Write output table
write.table(matrix, paste0(output_dir, "/Zscore_and_features_matrix.tab"), sep="\t", col.names=TRUE, row.names=F, quote=F)
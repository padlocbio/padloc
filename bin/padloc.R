#!/usr/bin/env Rscript
#
#         ___       ___                           ___       ___   
#        /  /\     /  /\     _____               /  /\     /  /\  
#       /  /::\   /  /::\   /  /::\   ___  __   /  /::\   /  /:/_
#      /  /:/:/  /  /:/::\ /  /:/\:\ /  /\/ /\ /  /:/\:\ /  /:/ /\
#      \  \::/   \  \::/\/ \  \:\/:/ \  \:\/:/ \  \:\/:/ \  \:\/:/
#       \  \:\    \  \:\    \  \::/   \  \::/   \  \::/   \  \::/  
#        \__\/     \__\/     \__\/     \__\/     \__\/     \__\/                 
# 
# padloc: Locate antiviral defence systems in prokaryotic genomes

# LOAD PACKAGES ---------------------------------------------------------------

library(yaml)
library(plyr)
suppressMessages(library(dplyr))
library(tidyr)
library(stringr)
library(rlang)
library(readr)
library(parallel)
library(getopt)
library(readxl)

# SCRIPT UTILITIES ------------------------------------------------------------

# debug_msg(message)
# print message when using debug
debug_msg <- function(msg) {
  if ( DEBUG_COUNTER > 0 ) {
    write(paste0("(", format(Sys.time(), "%X"), 
                 ") DEBUG  >>  ", msg), stdout())
  }
}

# warning_msg(message)
# print warning message
warning_msg <- function(msg) {
  if ( QUIET < 1 ) {
    write(paste0("(", format(Sys.time(), "%X"), 
                 ") WARNING  >>  ", msg), stdout())
  }
}

# die(message)
# print error message and quit
die <- function(msg) {
  write(paste0("\n(", format(Sys.time(), "%X"), 
               ") ERROR  >>  ", msg, "\n"), stdout())
  stop("Fatal error encountered, stopping padloc ...\n\n", call. = FALSE)
}

# ARGUMENT PARSING ------------------------------------------------------------

# set spec
spec = matrix(c(
  'domtbl_path' , 'd', 1, "character", "path to domain table",
  'featbl_path' , 'f', 1, "character", "path to feature table",
  'alias_path'  , 'a', 1, "character", "path to hmm alias table",
  'summary_path', 's', 1, "character", "path to system summary table",
  'yaml_dir'    , 'y', 1, "character", "path to yaml directory",
  'output_dir'  , 'o', 1, "character", "path to output directory",
  'debug'       , 'b', 1, "integer"  , "print debug messages",
  'quiet'       , 'q', 1, "integer"  , "suppress warnings"
), byrow = TRUE, ncol = 5)

# process options
opt = getopt::getopt(spec)

DOMTBL_PATH   <- opt$domtbl_path
FEATBL_PATH   <- opt$featbl_path
ALIAS_PATH    <- opt$alias_path
SUMMARY_PATH  <- opt$summary_path
YAML_DIR      <- opt$yaml_dir
OUTPUT_DIR    <- opt$output_dir
DEBUG_COUNTER <- opt$debug
QUIET         <- opt$quiet

### DEBUG ###
# DOMTBL_PATH   <- "~/tools/padloc/tests/domtblout/GCA_000830595.1_ASM83059v1.domtblout"
# FEATBL_PATH   <- "~/tools/padloc/tests/GCA_000830595.1_ASM83059v1_feature_table.txt"
# ALIAS_PATH    <- "~/tools/padloc/HMM_META.xlsx"
# SUMMARY_PATH  <- "~/tools/padloc/SYS_META.xlsx"
# YAML_DIR      <- "~/tools/padloc/SYS"
# OUTPUT_DIR    <- "~/tools/padloc/tests/"
# DEBUG_COUNTER <- 1
# QUIET         <- 0

# CHECK ARGUMENTS -------------------------------------------------------------

# list required parameters
required_parameters <- 
  list("domain table" = DOMTBL_PATH, "feature table" = FEATBL_PATH, 
       "alias table" = ALIAS_PATH, "system summary" = SUMMARY_PATH, 
       "yaml directory" = YAML_DIR, "output directory" = OUTPUT_DIR)

# list required files
required_files <- 
  c("domain table", "feature table", "alias table", "system summary")

# list required directories
required_directories <- 
  c("yaml directory", "output directory")

# check parameters
for ( i in 1:length(required_parameters) ) {
  
  parameter_name <- names(required_parameters[i])
  parameter_value <- required_parameters[i]
  
  # check that required parameters were assigned values
  if (is.null(parameter_value)) {
    die(paste0(parameter_name, " is a required parameter"))
  } else {
    debug_msg(paste0(parameter_name, ": ", parameter_value))
  }
  
  # check that files exist
  if ( parameter_name %in% required_files ) {
    if ( file.exists(as.character(parameter_value)) == FALSE ) {
      die(paste0(parameter_name, " file does not exist"))
    }
  }
  
  # check that directories exist and append "/" to path if needed
  if ( parameter_name %in% required_directories ) {

    if ( ! str_ends(parameter_value, "/") ) {
      required_parameters[i] <- paste0(eval(parse(text = "parameter_value")), "/")
    }
    
    if ( dir.exists(as.character(parameter_value)) == FALSE ) {
      die(paste0(parameter_name, " directory does not exist"))
    }
      
  }
  
}

# update directory paths to append "/"
YAML_DIR   <- unlist(unname(required_parameters["yaml directory"]))
OUTPUT_DIR <- unlist(unname(required_parameters["output directory"]))

# FUNCTIONS -------------------------------------------------------------------

# read_aliastbl(path to alias file)
# read in aliastbl file (i.e. hmm_meta.xlsx)
read_aliastbl <- function(aliastbl_path) {
  
  # read in the HMM alias table
  hmm_aliases <- read_xlsx(
    aliastbl_path, 
    skip = 1,
    col_names = c(
      "hmm.acc", "hmm.name", "hmm.description", 
      "protein.name", "system.definition.shortcut", "author",
      "originator", "number.seq", "length.hmm", 
      "e.value.threshold", "hmm.coverage.threshold", 
      "target.coverage.threshold", "system", "literature.ref",
      "database.ref", "comments")
    ) 
  
  # allow ambiguous protein name assignments
  hmm_aliases <- hmm_aliases %>% 
    separate_rows(protein.name, sep = "\\|") %>%
    separate_rows(system.definition.shortcut, sep = "\\|")
  
  # check for empty HMM parameters
  empty_params <- hmm_aliases %>%
    filter(is.na(hmm.name) == FALSE) %>%
    filter(is.na(e.value.threshold) | 
             is.na(hmm.coverage.threshold) | 
             is.na(target.coverage.threshold) == TRUE) %>%
    select(hmm.name, e.value.threshold, hmm.coverage.threshold, 
           target.coverage.threshold) %>%
    mutate_at(c("e.value.threshold", 
                "hmm.coverage.threshold", 
                "target.coverage.threshold"), 
              ~replace(., is.na(.), "MISSING"))
  
  if ( nrow(empty_params) != 0 ) {
    
    # warn of empty parameters
    warning_msg("HMMs missing required parameters\n")
    
    if ( QUIET < 1 ) {
      print.data.frame(empty_params, row.names = FALSE)
    }
    
    warning_msg("\nFilling missing parameters with default values")
    warning_msg("e.value.threshold = 1.00E-05")
    warning_msg("hmm.coverage.threshold = 0.5")
    warning_msg("target.coverage.threshold = 0.5")
    
    # fill with defaults
    hmm_aliases <- hmm_aliases %>%
      mutate_at("e.value.threshold", ~replace(., is.na(.), 1.00E-05)) %>%
      mutate_at("hmm.coverage.threshold", ~replace(., is.na(.), 0.5)) %>%
      mutate_at("target.coverage.threshold", ~replace(., is.na(.), 0.5))
    
  }
  
  return(hmm_aliases)
  
}

# read_systbl(path to systems summary file)
# read in the systbl file (i.e. sys_meta.xlsx) and get the names of systems
read_systbl <- function(systbl_path) {
  
  # read in the systems summary
  system_definitions <<- read_xlsx(
    SUMMARY_PATH,
    skip = 1,
    col_names = c(
      "system", "type", "yaml.name", "search", "within", 
      "around", "notes")) 
     
  
  # grab the names of systems that are marked TRUE for searching
  system_names <- system_definitions %>% 
    filter(search == "T") %>% 
    select(yaml.name) %>% 
    as.matrix()
  
}

# read_domtbl(path to domain table file)
# used for parsing domain table into a dataframe
read_domtbl <- function(domtbl_path) {
  
  # generate the assembly name from the domtbl path
  assembly_name <<- domtbl_path %>% 
    sub('.*\\/', '', .) %>% 
    sub('.domtblout', '', .)

  # adapted from: 
  # Zebulun Arendsee (2017). rhmmer: Utilities Parsing 'HMMER' Results. 
  # R package version 0.1.0. https://CRAN.R-project.org/package=rhmmer

  # specify column types
  col_types <- readr::cols(
    target.name = col_character(), 
    target.accession = col_character(),
    target.length = col_integer(), 
    hmm.name = col_character(),
    hmm.accession = col_character(), 
    hmm.length = col_integer(),
    full.seq.E.value = col_double(), 
    full.seq.score = col_double(),
    full.seq_bias = col_double(), 
    domain.number = col_integer(),
    total.domains = col_integer(), 
    domain.cE.value = col_double(),
    domain.iE.value = col_double(), 
    domain.score = col_double(),
    domain.bias = col_double(), 
    hmm.coord.from = col_double(),
    hmm.coord.to = col_double(), 
    alignment.coord.from = col_double(),
    alignment.coord.to = col_double(), 
    envelope.coord.from = col_double(),
    envelope.coord.to = col_double(), 
    accuracy = col_double(),
    target.description = col_character())
  
  readr::read_lines(domtbl_path) %>%
    # substitute the whitespace between the 'acc' and 'desctription' 
    # columns with '\t'
    sub(pattern = sprintf(
      "(%s) *(.*)", paste0(rep('\\S+', 22), collapse = " +")), 
      replacement = '\\1\t\\2') %>%
    # collapse everything to a single line
    paste0(collapse = "\n") %>%
    # re-parse the table as tsv
    readr::read_tsv(
      col_names = c('temp', 'target.description'), comment = '#', na = '-') %>%
    # separate the temp column into actual columns
    tidyr::separate(
      .data$temp, head(names(col_types$cols), -1), sep = ' +') %>%
    # apply colum types
    readr::type_convert(col_types = col_types)
    
    # TODO: COULD JUST SELECT FOR RELEVANT COLUMNS
    # select(
    #   target_name, 
    #   hmm_name, 
    #   full_seq_E_value, 
    #   domain_iE_value, 
    #   target_description
    # )
  
}

# read_featbl(path to feature table file)
# used for parsing feature table into a dataframe
read_featbl <- function(featbl_path) {

  read.table(
    file = featbl_path, 
    as.is = TRUE, 
    sep = '\t', 
    quote = "\"",
    col.names = c(
      "feature", "class", "assembly", "assembly.unit", "seq.type",
      "chromosome", "genomic.accession", "start", "end", "strand", 
      "target.name", "non-redundant.refseq", "related.accession", 
      "name", "symbol", "GeneID", "locus.tag", 
      "feature.interval.length", "product.length", 
      "attributes")) %>%
    # filter for proteins
    filter(target.name != "") %>%
    # arrange by contig and position in contig
    arrange(genomic.accession, start) %>%
    # TODO: COULD PROCESS STRANDS SEPARATELY?
    # dplyr::group_by(strand) %>%
    # add a column for the relative position of the gene in the genome
    dplyr::mutate(relative.position = row_number()) %>%
    # select relevant columns
    select(
      assembly,
      genomic.accession,
      target.name, 
      start, 
      end, 
      strand, 
      locus.tag, 
      relative.position
    )
  
}

# merge_tbls(domtbl, featbl, aliastbl)
# used for merging the domtbl, featbl, and aliastbl
merge_tbls <- function(domtbl, featbl, aliastbl) {
  
  # merge the domain table with the feature table and alias table, and arrange 
  # hits by genomic location
  merged <- domtbl %>% 
    left_join(featbl, by = "target.name") %>%
    left_join(aliastbl, by = "hmm.name") %>%
    arrange(relative.position)
  
  # check for unknown HMMs
  unknown_aliases <- merged %>%
    filter(is.na(protein.name) == TRUE) %>% 
    group_by(hmm.name) %>% 
    dplyr::summarise(counts = n())
  
  if ( nrow(unknown_aliases) != 0 ) {
    
    # warn of unknown HMMs
    warning_msg(paste0(nrow(unknown_aliases), 
                       " HMMs not listed in the aliases table"))
    warning_msg("These HMMs will not be included in analysis")
    warning_msg("Please add HMM metadata to hmm_aliases.xlsx\n")
    
    if ( QUIET < 1 ) {
      print.data.frame(unknown_aliases, row.names = FALSE)
      message()
    }
    
  }
  
  # ignore anything not in the alias list
  merged <- merged %>% filter(is.na(protein.name) != TRUE )  
  
  # calculate hit coverage
  merged <- merged %>%
    dplyr::mutate(
      hmm.coverage = (hmm.coord.to - hmm.coord.from) / hmm.length,
      target.coverage = (alignment.coord.to - alignment.coord.from) / 
        target.length, combined.coverage = hmm.coverage + target.coverage)
  
}

# search_system(name of system, merged domtbl/featbl dataframe)
# used for identifying complete defence systems in the domain table data
search_system <- function(system_type, merged_tbls) {
  
  # generate path to YAML file
  yaml_path <- paste0(YAML_DIR, system_type, ".yaml")
  
  # check that the YAML file exists
  if( ! file.exists(yaml_path) ) {
    
    warning_msg(paste0("YAML does not exist for system: ", system_type))
    return(NULL)
    
  }
  
  # read in the yaml file
  system_param <- yaml::read_yaml(yaml_path)
  
  max_space        <- system_param$maximum_intergenic_space
  min_core         <- system_param$minimum_core
  min_total        <- system_param$minimum_total
  core_genes       <- system_param$core_genes
  other_genes      <- system_param$other_genes
  prohibited_genes <- system_param$prohibited_genes
  
  gene_types <- list(core_genes = core_genes, other_genes = other_genes, 
                     prohibited_genes = prohibited_genes)
  
  for ( i in 1:length(gene_types) ) {
    
    gene_type_name <- names(gene_types[i])
    gene_type_values <- gene_types[i]
    genes <- unlist(gene_type_values)
    
    for ( gene in genes ) {
      
      # assign protein names based on system.definition.shortcut
      if ( gene %in% aliastbl$system.definition.shortcut ) {
        additional_genes <- 
          unlist(aliastbl %>% filter(system.definition.shortcut == gene) %>% 
                   select(protein.name), use.names = FALSE)
        do.call("<-", list(eval(parse(text = "gene_type_name")), 
                           unique(c(gene_type_values, additional_genes))))
      }
    }
  }
  
  # remove any genes from 'other' that are also listed as 'core' or 'prohibited'
  other_genes <- unlist(other_genes %>% setdiff(c(core_genes, prohibited_genes)))
  # remove any genes from 'prohibited' that are also listed as 'core'
  prohibited_genes <- unlist(prohibited_genes %>% setdiff(core_genes))
  core_genes <- unlist(core_genes)
  
  # filter for relevant genes
  relevance_check <- merged_tbls %>%
    dplyr::mutate(relevant = ifelse(
      protein.name %in% c(core_genes, other_genes, prohibited_genes), 
      TRUE, FALSE))
  
  relevance_checked <- relevance_check %>%
    filter(relevant == TRUE)
  
  # specify gene classifications
  genes_classified <- relevance_checked %>% 
    dplyr::mutate(
      is.core = ifelse(protein.name %in% core_genes, TRUE, FALSE),
      is.accessory = ifelse(protein.name %in% other_genes, TRUE, FALSE),
      is.prohibited = ifelse(protein.name %in% prohibited_genes, TRUE, FALSE))
  
  # add column listing all domains for a particular hmm hiting a particular 
  # target
  collapsed_domains <- genes_classified %>%
    group_by(assembly, locus.tag) %>% 
    dplyr::mutate(all.domains = paste0(
      hmm.name, ", ", domain.number, ", ", domain.iE.value, ", ", 
      target.coverage, ", ", hmm.coverage, collapse = " | ")) %>%
    ungroup()
  
  # filter for top scoring domain per hmm
  top_domain <- collapsed_domains %>%
    group_by(assembly, locus.tag, hmm.name) %>% 
    top_n(-1, domain.iE.value) %>% 
    ungroup()
  
  # add column listing all hits for a particular target
  collapsed_hits <- top_domain %>%
    group_by(assembly, target.name) %>% 
    dplyr::mutate(all.hits = paste0(
      hmm.name, ", ", hmm.accession, ", ", full.seq.E.value, 
      collapse = " | ")) %>%
    ungroup()
  
  # filter for acceptable coverage
  coverage_check <- collapsed_hits %>% 
    dplyr::mutate(coverage.check = ifelse(
      hmm.coverage >= hmm.coverage.threshold & target.coverage >= 
        target.coverage.threshold, 
      TRUE, FALSE))
  
  coverage_checked <- coverage_check %>%
    filter(coverage.check == TRUE)
  
  # filter for acceptable e-value
  e_value_check <- coverage_checked %>% 
    dplyr::mutate(e.value.check = ifelse(
      domain.iE.value <= e.value.threshold, TRUE, FALSE))
  
  e_value_checked <- e_value_check %>%
    filter(e.value.check == TRUE)
  
  # for each target, order each hit by type (core > accessory > prohibited)
  # then by e-value (low > high), then assign rank (1 is best)
  ranked_hits <- e_value_checked %>%
    group_by(assembly, locus.tag) %>% 
    dplyr::arrange(
      desc(is.core),
      desc(is.accessory),
      desc(is.prohibited),
      domain.iE.value,
      desc(combined.coverage),
      .by_group = TRUE) %>%  
    dplyr::mutate(rank = row_number()) %>%
    ungroup()
  
  # filter for top ranked hits for each target
  top_hits <- ranked_hits %>%
    filter(rank == 1)
  
  # add a column that groups hits into clusters
  clusters <- top_hits %>% 
    arrange(relative.position) %>%
    dplyr::mutate(cluster = cumsum(c(1, abs(diff(relative.position)) > max_space)))
  
  # count the number of unique hits in each cluster
  clusters_unique <- clusters %>%
    group_by(cluster, is.core, is.accessory, is.prohibited) %>%
    dplyr::mutate(unq.counts = length(unique(protein.name))) %>%
    ungroup()
  
  # determine the actual number of hits per cluster
  clusters_counts <- clusters_unique %>%
    dplyr::mutate(
      unq.core = is.core * unq.counts,
      unq.accessory = is.accessory * unq.counts,
      unq.prohibited = is.prohibited * unq.counts) %>%
    group_by(cluster) %>%
    dplyr::mutate(
      unq.core = max(unq.core),
      unq.accessory = max(unq.accessory),
      unq.prohibited = max(unq.prohibited)) %>%
    dplyr::mutate(unq.total = unq.core + unq.accessory) %>%
    ungroup()
  
  # identify candidate systems
  candidates_check <- clusters_counts %>% 
    dplyr::mutate(candidate = ifelse(
      # all core genes, at least the minimum number of total genes, and no 
      # prohibited genes must be present
      unq.core >= min_core & unq.total >= min_total & 
        unq.prohibited == 0, TRUE, FALSE))
  
  candidates_checked <- candidates_check %>% 
    filter(candidate == TRUE)
  
  ### return genes within/around the defence system ###

  # check that a system was actually found
  if ( nrow(candidates_checked) != 0 ) {

    # check whether we want genes within cluster
    within_check <- system_definitions %>%
      filter(yaml.name == system_type) %>%
      select(within)

    ifelse(is.na(within_check) == TRUE,
           within_check <- FALSE,
           within_check <- as.logical(within_check))

    if ( within_check == TRUE ) {

      # grab the limits of the system (for each instance of the system)
      within_limits <- candidates_checked %>%
        group_by(cluster) %>%
        do(data.frame(lower = min(.$relative.position),
                      upper = max(.$relative.position))) %>%
        ungroup()

      # pull out the genes from the feature table within those limits
      genes_within <- apply(X = within_limits, MARGIN = 1, FUN = pull_features,
                            featbl = featbl, range = 0)

      genes_within <- ldply(genes_within, rbind, .id = NULL) %>%
        dplyr::mutate(system = system_type)

      genes_within_out <<- genes_within_out %>% rbind(genes_within)

    }

    # check whether we want genes around the cluster
    around_check <-
      system_definitions %>%
      filter(yaml.name == system_type) %>%
      select(around)

    ifelse(is.na(around_check) == TRUE,
           around_check <- FALSE,
           around_check <- as.numeric(around_check))

    if ( around_check > 0 ) {

      around_limits <- candidates_checked %>%
        group_by(cluster) %>%
        do(data.frame(lower = min(.$relative.position),
                      upper = max(.$relative.position))) %>%
        ungroup()

      genes_around <- apply(
        X = around_limits,
        MARGIN = 1,
        FUN = pull_features,
        featbl = featbl,
        range = as.numeric(around_check))
      
      genes_around <- ldply(genes_around, rbind, .id = NULL) %>%
        dplyr::mutate(system = system_type)
      
      genes_around_out <<- genes_around_out %>% rbind(genes_around)

    }

  }
  
  # select relevant columns for final output
  results <- candidates_checked %>%
    dplyr::mutate(system = system_type) %>%
    select(assembly, genomic.accession, system, target.name, hmm.name, protein.name, 
           full.seq.E.value, domain.iE.value, target.coverage, hmm.coverage, 
           locus.tag, start, end, strand, target.description, relative.position, 
           all.domains, all.hits)
    
}

# pull_features(dataframe of relative positions, feature table, range)
# used for pulling out the information of genes, specified by relative position
pull_features <- function(df, featbl, range) {

  # grab limits and modify by range
  lower_limit <- df[2] - range
  upper_limit <- df[3] + range

  # filter for genes within limits
  featbl <- featbl %>% 
    filter(relative.position %in% lower_limit:upper_limit)
  
}

# generate_gff(padloc output dataframe)
# used to generate an annotation file for the padloc output table
generate_gff <- function(padloc_out) {
  
# generate an annotation file for Geneious (.gff variant)
gff <- padloc_out %>%
  dplyr::mutate(
    seqid = gsub("\\..*", "", genomic.accession),
    source = "padloc",
    type = system,
    score = as.numeric(domain.iE.value),
    phase = ".",
    attributes = paste0(
      "Name=", protein.name, ";HMM=", hmm.name,
      ";Target.coverage=", target.coverage, 
      ";HMM.coverage=", hmm.coverage)) %>%
  select(seqid, source, type, start, end, score, strand, phase, attributes)

}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# MAIN
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

start_time <- Sys.time()
debug_msg(paste0("Start time: ", start_time))

# read in the hmm alias table
debug_msg(paste0("Reading ", basename(ALIAS_PATH)))
aliastbl <- read_aliastbl(ALIAS_PATH)

# read in the systems summary
debug_msg(paste0("Reading ", basename(SUMMARY_PATH)))
systbl <- read_systbl(SUMMARY_PATH) 

# read in domain table
debug_msg(paste0("Reading ", basename(DOMTBL_PATH)))
domtbl <- read_domtbl(DOMTBL_PATH)

# read in feature table
debug_msg(paste0("Reading ", basename(FEATBL_PATH)))
featbl <- read_featbl(FEATBL_PATH)

# merge the domain, feature, and alias tables
debug_msg("Merging domain, alias, and feature tables")
merged <- merge_tbls(domtbl, featbl, aliastbl)

# set the genes_within_out and genes_around_out dataframe
genes_within_out <- genes_around_out <- data.frame(
    "assembly" = character(), "genomic_accession" = character(),
    "target_name" = character(), "start" = character(), "end" = character(), 
    "strand" = character(), "locus_tag" = character(), 
    "relative_position" = character(), stringsAsFactors = FALSE)

# search for systems
debug_msg("Searching for defence systems")
systems <- lapply(X = systbl, FUN = search_system, merged_tbls = merged)
# rename the dataframes generated above
names(systems) <- systbl
# merge into a single table
padloc_out <- ldply(systems, rbind, .id = NULL)

# TODO: POST-PROCESSING STEP WILL PROBABLY GO HERE? 

# output table of defence systems and annotation file
if ( nrow(padloc_out) > 0 ) {
  
  padloc_out %>%
    write_csv(path = paste0(OUTPUT_DIR, "/", assembly_name, ".csv"))

  gff <- generate_gff(padloc_out)
  gff %>% write_delim(path = paste0(OUTPUT_DIR, "/", assembly_name, ".gff"),
                      delim = "\t", col_names = FALSE)
    
}

# output table of genes within defence systems
if ( nrow(genes_within_out) > 0 ) {
  
  genes_within_out %>% 
    write_csv(path = paste0(OUTPUT_DIR, "/", assembly_name, "_within.csv"))
  
}

# output table of genes around defence systems
if ( nrow(genes_around_out) > 0 ) {
  
  genes_around_out %>% 
    write_csv(path = paste0(OUTPUT_DIR, "/", assembly_name, "_around.csv"))
  
}

end_time <- Sys.time()
debug_msg(paste0("End time: ", end_time))
run_time <- end_time - start_time
debug_msg(paste0("Run time: ", run_time))

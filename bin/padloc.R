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
# PADLOC: Locate antiviral defence systems in prokaryotic genomes

# LOAD PACKAGES ----------------------------------------------------------------

suppressMessages(library(tidyverse))
library(getopt)
library(yaml)
options(dplyr.summarise.inform = FALSE)

# SCRIPT UTILITIES -------------------------------------------------------------

# Print message.
msg <- function(msg) {
  if ( QUIET < 1 ) {
    write(paste0("[", format(Sys.time(), "%X"), "] >> ", msg), stdout())
  }
}

# Print message when using debug.
debug_msg <- function(msg) {
  if ( DEBUG_COUNTER > 0 ) {
    write(paste0("[", format(Sys.time(), "%X"), "] DEBUG >> ", msg), stdout())
  }
}

# Print warning message and issue warning.
warning_msg <- function(msg) {
  if ( QUIET < 1 ) {
    write(paste0("[", format(Sys.time(), "%X"), "] WARNING >> ", msg), stdout())
  }
}

# Print message and exit.
die <- function(msg) {
  write(paste0("[", format(Sys.time(), "%X"), "] ERROR >> ", msg), stdout())
  invokeRestart("abort")
}

# ARGUMENT PARSING -------------------------------------------------------------

# set spec
spec = matrix(c(
  'domtbl_path'  , 'd', 1, "character", "path to domain table",
  'gff_file_path', 'f', 1, "character", "path to feature table",
  'crispr_file_path', 'c', 1, "character", "optional path to crispr input",
  'hmm_meta_path', 'h', 1, "character", "path to hmm_meta table",
  'sys_meta_path', 's', 1, "character", "path to system summary table",
  'yaml_dir'     , 'y', 1, "character", "path to yaml directory",
  'output_dir'   , 'o', 1, "character", "path to output directory",
  'debug'        , 'b', 1, "integer"  , "print debug messages",
  'quiet'        , 'q', 1, "integer"  , "suppress warnings",
  'prodigal'     , 'p', 1, "integer"  , "prodigal was run"
), byrow = TRUE, ncol = 5)

# process options
opt = getopt::getopt(spec)

DOMTBL_PATH   <- opt$domtbl_path
GFF_PATH      <- opt$gff_file_path
CRISPR_PATH   <- opt$crispr_file_path
HMM_META_PATH <- opt$hmm_meta_path
SYS_META_PATH <- opt$sys_meta_path
YAML_DIR      <- opt$yaml_dir
OUTPUT_DIR    <- str_remove(opt$output_dir, "\\/$")
DEBUG_COUNTER <- opt$debug
QUIET         <- opt$quiet
PRODIGAL      <- opt$prodigal

### DEBUG ###
DOMTBL_PATH   <- "~/tools/padloc/test/output/GCF_003182315.1_faa.domtblout"
GFF_PATH      <- "~/tools/padloc/test/input/GCF_003182315.1_faa.gff"
HMM_META_PATH <- "~/tools/padloc/data/hmm_meta.txt"
SYS_META_PATH <- "~/tools/padloc/data/sys_meta.txt"
YAML_DIR      <- "~/tools/padloc/data/sys/"
OUTPUT_DIR    <- "~/tools/padloc/debug"
DEBUG_COUNTER <- 1
QUIET         <- 0
PRODIGAL      <- 0
# 
# DOMTBL_PATH   <- "D:/git_clone/padloc/test/0954a438-c024-455d-9e51-5d1f9133c931.domtblout"
# GFF_PATH      <- "D:/git_clone/padloc/test/0954a438-c024-455d-9e51-5d1f9133c931_prodigal.gff"
# CRISPR_PATH   <- "D:/git_clone/padloc/test/0954a438-c024-455d-9e51-5d1f9133c931_crispr.txt.gff"
# HMM_META_PATH <- "D:/git_clone/padloc/data/hmm_meta.txt"
# SYS_META_PATH <- "D:/git_clone/padloc/data/sys_meta.txt"
# YAML_DIR      <- "D:/git_clone/padloc/data/sys/"
# OUTPUT_DIR    <- "D:/git_clone/padloc/test/"
# DEBUG_COUNTER <- 1
# QUIET         <- 0
# PRODIGAL      <- 0



# FUNCTIONS --------------------------------------------------------------------

# read_hmm_meta(path to hmm_meta)
# Read in hmm_meta file.
read_hmm_meta <- function(file) {
  
  cols <- cols(
    hmm.acc = col_character(),
    hmm.name = col_character(),
    hmm.description = col_character(),
    protein.name = col_character(),
    system.definition.shortcut = col_character(),
    author = col_character(),
    number.seq = col_double(),
    length.hmm = col_double(),
    e.value.threshold = col_double(),
    hmm.coverage.threshold = col_double(),
    target.coverage.threshold = col_double(),
    system = col_character(),
    literature.ref = col_character(),
    database.ref = col_character(),
    comments = col_character()
  )
  
  # read in hmm_meta
  raw <- read_tsv(
    file,
    skip = 1,
    col_names = names(cols$cols),
    col_types = cols
   )
  
  fail <- raw %>%
    filter(is.na(hmm.acc) | is.na(hmm.name) | is.na(protein.name))
  
  if(nrow(fail) > 0) {
    warning_msg("Some rows of hmm_meta.txt are missing values in required columns (hmm.accession, hmm.name, protein.name)\n")
    print.data.frame(fail)
    die("Failed to read hmm_meta.txt")
  }
  
  warn <- raw %>%
    filter(is.na(e.value.threshold) | is.na(hmm.coverage.threshold) | is.na(target.coverage.threshold))
  
  if(nrow(warn) > 0 & QUIET == 0) {
    warning_msg("Some rows of hmm_meta.txt are missing values in required columns (e.val.threshold, hmm.coverage.threshold, target.coverage.threshold)\n")
    print.data.frame(warn)
    warning_msg("These columns will be filled with default values, respectively: 1E-05, 0.3, 0.3")
  }
  
  out <- raw %>%
    # remove space characters from ambiguous protein assignments (e.g. AbcD | EfgH should be AbcD|EfgH)
    mutate(protein.name = str_remove_all(protein.name," "),
           system.definition.shortcut = str_remove_all(system.definition.shortcut," "))%>%
    # allow ambiguous protein name assignments
    separate_rows(protein.name, sep = "\\|") %>%
    separate_rows(system.definition.shortcut, sep = "\\|") %>%
    # Fill empty fields with defaults
    mutate(
      e.value.threshold = ifelse(is.na(e.value.threshold), 1.00E-05, e.value.threshold),
      hmm.coverage.threshold = ifelse(is.na(hmm.coverage.threshold), 0.3, hmm.coverage.threshold),
      target.coverage.threshold = ifelse(is.na(target.coverage.threshold), 0.3, target.coverage.threshold),
    )
  
}

# read_sys_meta(path to sys_meta)
# Read in the sys_meta file.
read_sys_meta <- function(file) {
  
  cols <- cols(
    system = col_character(),
    type = col_character(),
    yaml.name = col_character(),
    search = col_logical(),
    notes = col_character()
  )
  
  # read in sys_meta
  out <- read_tsv(
    file, 
    skip = 1,
    col_names = names(cols$cols),
    col_types = cols
  )

}

# read_gff(path to gff file)
# Read in a GFF file.
read_gff <- function(gff_path) {
  
  cols <- cols(
    seqid = col_character(),
    source = col_character(),
    type = col_character(),
    start = col_double(),
    end = col_double(),
    score = col_character(),
    strand = col_character(),
    phase = col_character(),
    attributes = col_character()
  )
  
  out <- read_tsv(
    gff_path,
    col_names = names(cols$cols),
    col_types = cols,
    comment = "#"
  )
  
}

# separate_attributes(gff)
# Separate the attributes column of a GFF file.
separate_attributes <- function(gff) {
  
  out <- gff %>%
    separate_rows(attributes, sep = ";") %>%
    filter(attributes != "") %>%
    separate(attributes, into = c("key", "value"), sep = "=") %>% 
    spread(key = key, value = value, fill = NA)
  
}

# gather_attributes(gff)
# Consolidate the attributes column of a GFF file.
gather_attributes <- function(gff) {
  
  # WARNING: Attributes are not necessarily gathered back into original order
  
  fields <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase")
  attributes <- names(gff) %>% setdiff(fields)
  
  out <- gff_separated %>% 
    pivot_longer(cols = attributes, names_to = "attribute", values_to = "value") %>%
    filter(is.na(value) == FALSE) %>%
    mutate(attribute.value = paste(attribute, value, sep = "=")) %>%
    group_by(seqid, start, end) %>%
    mutate(element = paste(seqid, start, end, sep = ";")) %>%
    mutate(attributes = paste(attribute.value, collapse = ";")) %>%
    distinct(element, .keep_all = TRUE) %>%
    select(all_of(fields), attributes)
  
}

# read_domtbl(path to domtbl)
# Read in domtbl.
read_domtbl <- function(domtbl_path) {

  # adapted from: 
  # Zebulun Arendsee (2017). rhmmer: Utilities Parsing 'HMMER' Results. 
  # R package version 0.1.0. https://CRAN.R-project.org/package=rhmmer
  
  cols <- cols(
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
  
  out <- read_lines(domtbl_path) %>%
    # substitute the whitespace between the 'acc' and 'description' 
    # columns with '\t'
    sub(pattern = sprintf("(%s) *(.*)", paste0(rep('\\S+', 22), collapse = " +")), replacement = '\\1\t\\2') %>%
    # remove '#'s that aren't comments, stops breaking w/ subsequent read_tsv()
    str_remove_all(pattern = '(?<!^)#') %>%
    # collapse everything to a single line
    paste0(collapse = "\n") %>%
    # re-parse the table as tsv
    read_tsv(col_names = c('temp', 'target.description'), comment = "#", na = '-', show_col_types = FALSE) %>%
    # separate the temp column into actual columns
    separate(.data$temp, head(names(cols$cols), -1), sep = ' +') %>%
    # apply colum types
    type_convert(col_types = cols)
  
}

# Read and process CRISPRDetect.gff outputs to add CRISPR arrays. Requires pre-loading the protein-coding gene table.
process_crisprdetect_gff<-function(CRISPR_PATH,gff.genes){
  
  crispr.in <- read_gff(CRISPR_PATH)
  
  if(nrow(crispr.in) == 0) {
    return(tibble(NULL))
  }
  
  crispr <- crispr.in %>%
                    filter(type=="repeat_region") %>%
                    # CRISPRDetect doesn't continue CRISPR numbering across contigs... 
                    arrange(seqid,start) %>%
                    mutate(array.id = paste0("CRISPR",sprintf("%03d",row_number()))) %>%
                    group_by(array.id) %>%
                    separate_attributes() %>%
                    mutate(type = "CRISPR_array",
                           target.name = array.id,
                           crispr.repeat = Note,
                           domain.iE.value=as.numeric(Array_quality_score),
                           target.description = paste0(array.id,"; repeat=",crispr.repeat,"; score=",Array_quality_score)) %>%
                    # reformat seqid for CRISPRDetect output
                    mutate(seqid = str_remove_all(seqid,"-.*")) %>%
                    # assing temp relative.position in genome
                    group_by(seqid) %>%
                    arrange(seqid,start)%>%
                    mutate(relative.position.offset = row_number()/100,  # assign CRISPR array positions as decimals (note, edge case if n>=100)
                           relative.position = 0)
  
  # need to update relative.position
  crispr.processed <- gff.genes %>% mutate(relative.position.offset = 0) %>%
            bind_rows(crispr) %>%
            arrange(seqid,start) %>%
            group_by(seqid) %>%
            mutate(contig.end = max(relative.position)) %>%
            mutate(relative.position = cummax(relative.position)) %>%
            mutate(relative.position = relative.position + relative.position.offset) %>%
            mutate(protein.name = "CRISPR_array") %>%
            ungroup() %>%
            filter(type == "CRISPR_array") %>%
            select(seqid,
                   type,
                   protein.name,
                   start,
                   end,
                   strand,
                   domain.iE.value, # gets converted to 'score' later
                   relative.position,
                   contig.end,
                   target.name,
                   target.description)
    
  return(crispr.processed)
}

# merge_tbls(domtbl, GFF, hmm_meta)
# Merge the domtbl, GFF, and hmm_meta
merge_tbls <- function(domtbl, gff, hmm_meta) {
  
  # Join domtbl with GFF, then with hmm_meta
  merged <- domtbl %>% 
    left_join(gff, by = "target.name") %>%
    left_join(hmm_meta, by = "hmm.name") %>%
    arrange(seqid,relative.position)
  
  # Catch HMMs that are missing from hmm_meta
  warn <- merged %>% 
    filter(is.na(protein.name))
  
  if(nrow(warn) > 0 & QUIET == 0) {
    warning_msg("hmm_meta.tsv - Found HMMs that aren't listed in hmm_meta - see below")
    message()
    print.data.frame(warn)
    message()
    warning_msg("hmm_meta.tsv - These HMMs will not be included in analysis")
  }
  
  out <- merged %>% 
    # Ignore anything not in the alias list
    filter(!is.na(protein.name)) %>%
    # Calculate hit coverage
    mutate(
      hmm.coverage = (hmm.coord.to - hmm.coord.from) / hmm.length,
      target.coverage = (alignment.coord.to - alignment.coord.from) / 
        target.length, combined.coverage = hmm.coverage + target.coverage) %>%
    mutate(
      hmm.coverage = round(hmm.coverage, 3),
      target.coverage = round(target.coverage, 3)
    ) %>%
    # remove the system.definition.shortcut
    select(-system.definition.shortcut) %>%
    distinct()
  
}

# search_system(system name, merged tbls)
# Main system identification logic.
search_system <- function(system_type, merged_tbls) {
  
  # system_type<-"cas_type_arrays"
  # system_type<-"septu_other"
  # merged_tbls <- merged
  
  # generate path to YAML file
  yaml_path <- paste0(YAML_DIR, system_type, ".yaml")
  
  # check that the YAML file exists
  if( ! file.exists(yaml_path) ) {
    warning_msg(paste0("YAML does not exist for system: ", system_type))
    return(NULL)
  }
  
  # read in the yaml file
  system_param <- yaml::read_yaml(yaml_path)
  max_space        <- system_param$maximum_separation
  min_core         <- system_param$minimum_core
  min_total        <- system_param$minimum_total
  core_genes       <- system_param$core_genes
  optional_genes   <- system_param$optional_genes
  prohibited_genes <- system_param$prohibited_genes
  
  # check boolean parameters
  force_strand<-ifelse("force_strand" %in% names(system_param),system_param$force_strand,F)
  
  expand_secondary_gene_assignments<-function(primary_gene_list){
  
    additional_genes<-hmm_meta %>% filter(system.definition.shortcut %in% primary_gene_list) %>% 
                                    select(protein.name) %>%
                                    distinct()
                    
    return(c(primary_gene_list,unlist(additional_genes,use.names = F)))
  }
  
  core_genes<-expand_secondary_gene_assignments(core_genes)
  optional_genes<-expand_secondary_gene_assignments(optional_genes)
  prohibited_genes<-expand_secondary_gene_assignments(prohibited_genes)
  
  # remove any genes from 'other' that are also listed as 'core' or 'prohibited'
  optional_genes <- unlist(optional_genes %>% setdiff(c(core_genes, prohibited_genes)))
  # remove any genes from 'prohibited' that are also listed as 'core'
  prohibited_genes <- unlist(prohibited_genes %>% setdiff(core_genes))
  core_genes <- unlist(core_genes)
  
  # check whether to add CRISPR arrays
  if (include.crispr.arrays==T & "CRISPR_array" %in% c(core_genes,optional_genes,prohibited_genes)){
    crispr.add <- T} else {
    crispr.add <- F}

  
  # filter for relevant genes
  relevance_check <- merged_tbls %>%
    dplyr::mutate(relevant = ifelse(
      protein.name %in% c(core_genes, optional_genes, prohibited_genes), 
      TRUE, FALSE))
  
  relevance_checked <- relevance_check %>%
    filter(relevant == TRUE)
  
  # specify gene classifications
  genes_classified <- relevance_checked %>% 
    dplyr::mutate(
      is.core = ifelse(protein.name %in% core_genes, TRUE, FALSE),
      is.accessory = ifelse(protein.name %in% optional_genes, TRUE, FALSE),
      is.prohibited = ifelse(protein.name %in% prohibited_genes, TRUE, FALSE))
  
  # add column listing all domains for a particular hmm hiting a particular 
  # target
  collapsed_domains <- genes_classified %>%
    group_by(seqid, start, end) %>% 
    dplyr::mutate(all.domains = paste0(
      hmm.accession, ",", hmm.name, ", ", domain.number, ", ", 
      domain.iE.value, ", ", target.coverage, ", ", hmm.coverage, 
      collapse = " | ")) %>%
    ungroup()
  
  # filter for top scoring domain per hmm
  top_domain <- collapsed_domains %>%
    group_by(seqid, start, end, hmm.name) %>% 
    top_n(-1, domain.iE.value) %>% 
    ungroup()
  
  # add column listing the best domain of the top 5 hits for a particular target
  collapsed_hits <- top_domain %>%
    group_by(seqid, start, end) %>%
    dplyr::arrange(domain.iE.value, .by_group = TRUE) %>%  
    top_n(-5, domain.iE.value) %>%
    dplyr::mutate(best.hits = paste0(
      hmm.accession, ", ", hmm.name, ", ", domain.iE.value, ", ", 
      target.coverage, ", ", hmm.coverage, 
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
    group_by(seqid, start, end) %>% 
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
    filter(rank == 1) %>% 
    # avoid type incompatibility when binding CRISPR arrays and systems later
    mutate(target.description=as.character(target.description))
  
  # add CRISPR array data (if required)
  
  if(crispr.add == T){
    to.add <- crispr.arrays %>%
      mutate(
      is.core = ifelse("CRISPR_array" %in% core_genes, TRUE, FALSE),
      is.accessory = ifelse("CRISPR_array" %in% optional_genes, TRUE, FALSE),
      is.prohibited = ifelse("CRISPR_array" %in% prohibited_genes, TRUE, FALSE)
      )
    top_hits <- top_hits %>% bind_rows(to.add)
  }
  
  # add a column that groups hits into clusters
    # first, set grouping by strand if force_strand == T
    # also make sure every cluster id is unique
      if (force_strand == T) {
        clusters <- top_hits %>% 
          arrange(seqid,strand,relative.position) %>% 
          group_by(seqid,strand) %>% 
          dplyr::mutate(cluster.tmp = cumsum(c(1, abs(diff(relative.position)) > max_space + 1))) %>%
          group_by(seqid,strand,cluster.tmp) %>% 
          mutate(cluster = cur_group_id()) %>%
          ungroup() %>%
          select(-cluster.tmp)
      } else {
        clusters <- top_hits %>% 
          arrange(seqid,relative.position) %>% 
          group_by(seqid) %>% 
          dplyr::mutate(cluster.tmp = cumsum(c(1, abs(diff(relative.position)) > max_space + 1))) %>%
          group_by(seqid,cluster.tmp) %>% 
          mutate(cluster = cur_group_id()) %>%
          ungroup() %>%
          select(-cluster.tmp)
      }
  
  # count the number of unique hits in each cluster
  clusters_unique <- clusters %>%
    group_by(seqid, cluster, is.core, is.accessory, is.prohibited) %>%
    dplyr::mutate(unq.counts = length(unique(protein.name))) %>%
    ungroup()
  
  # determine the actual number of hits per cluster
  clusters_counts <- clusters_unique %>%
    dplyr::mutate(
      unq.core = is.core * unq.counts,
      unq.accessory = is.accessory * unq.counts,
      unq.prohibited = is.prohibited * unq.counts) %>%
    group_by(seqid, cluster) %>%
    dplyr::mutate(
      unq.core = suppressWarnings(max(unq.core)),
      unq.accessory = suppressWarnings(max(unq.accessory)),
      unq.prohibited = suppressWarnings(max(unq.prohibited))) %>%
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
  
  # select relevant columns for final output
  results <- candidates_checked %>%
    dplyr::mutate(system = system_type) %>%
    select(seqid, cluster, system, target.name, hmm.accession, 
           hmm.name, protein.name, full.seq.E.value, domain.iE.value, 
           target.coverage, hmm.coverage, start, end, strand, 
           target.description, relative.position, contig.end, all.domains, best.hits)

}

# generate_gff(PADLOC output dataframe)
# Generate a GFF file from the PADLOC output table
generate_gff <- function(padloc_out) {
  
  # generate a .gff variant
  gff <- padloc_out %>%
    mutate(
      source = "padloc",
      type = system,
      score = as.numeric(domain.iE.value),
      phase = ".",
      attributes = paste0(
        "Name=", protein.name, ";HMM.accession=", hmm.accession,
        ";HMM.name=", hmm.name, ";Target.coverage=", target.coverage, 
        ";HMM.coverage=", hmm.coverage, ";Alt=", best.hits)) %>%
    select(seqid, source, type, start, end, score, strand, phase, attributes)
  
}

# MAIN -------------------------------------------------------------------------

# Record start time.
start_time <- Sys.time()
debug_msg(paste0("Start time: ", start_time))

# Read in hmm_meta.
debug_msg(paste0("Reading ", basename(HMM_META_PATH)))
hmm_meta <- read_hmm_meta(HMM_META_PATH)

# Read in sys_meta.
debug_msg(paste0("Reading ", basename(SYS_META_PATH)))
sys_meta <- read_sys_meta(SYS_META_PATH) 
# Get names of systems that are marked TRUE for searching.
system_names <- sys_meta %>% 
  filter(search == TRUE) %>% 
  select(yaml.name) %>% 
  as.matrix()

# Generate the assembly name from the domtbl path.
assembly_name <- DOMTBL_PATH %>% 
  sub('.*\\/', '', .) %>% 
  sub('.domtblout', '', .)
# Read in domtbl.
debug_msg(paste0("Reading ", basename(DOMTBL_PATH)))
domtbl <- read_domtbl(DOMTBL_PATH)

# Read in GFF that specifies protein-coding genes.
debug_msg(paste0("Reading ", basename(GFF_PATH)))
gff <- read_gff(GFF_PATH)

# Process GFF
gff <- gff %>%
  filter(type == "CDS") %>%
  separate_attributes()
#if ("pseudo" %in% names(gff)) { gff <- gff %>% filter(is.na(pseudo)) }
if ("pseudo" %in% names(gff)) { gff <- gff %>% mutate(ID = ifelse(is.na(pseudo),ID,Name)) }
if (!"protein_id" %in% names(gff)){gff <- gff %>% mutate(protein_id=ID)}

gff <- gff %>%
  arrange(seqid, start) %>%
  group_by(seqid) %>%
  mutate(relative.position = row_number()) %>%
  mutate(contig.end = max(relative.position)) %>%
  group_by(ID) %>%
  mutate(
    target.name = case_when(
                            # Fix prodigal-generated GFF
                            PRODIGAL == 1 ~  paste0(seqid, str_remove(ID, "^[0-9]*")),
                            # RefSeq & GenBank 'ID' fields have 'cds-' in front of the protein name 
                            # that needs to be removed, prokka GFFs are fine.
                            grepl("cds-",ID) ~ str_remove(ID, "cds-"),
                            # ensure pseudogenes are updated
                            grepl("pseudo_sub",ID) ~ ID,
                            # old genbank files don't always have the same ID format
                            T ~ protein_id)
  ) %>%
  ungroup()

# Check that all proteins listed in the .faa are also in the .gff

# Merge tables.
debug_msg("Merging domain, alias, and feature tables")
merged <- merge_tbls(domtbl, gff, hmm_meta)

# Warn and exit if there are proteins not found in the .GFF
unknown_protein_count <- nrow(merged %>%
                                filter(is.na(start)) %>%
                                select(target.name) %>%
                                distinct()
                               )

if (unknown_protein_count > 0) {
  die(paste0(unknown_protein_count, " protein sequence IDs are missing from GFF file"))
}

# Read in the optional CRISPR file (if present)

if(CRISPR_PATH != ""){
    debug_msg(paste0("Include CRISPR arrays from ", basename(CRISPR_PATH)))
    crispr.arrays <- process_crisprdetect_gff(CRISPR_PATH,gff)

    if(nrow(crispr.arrays) == 0){      
      debug_msg("CRISPR array file is empty, arrays will not be included in output")
      include.crispr.arrays <- F
    } else {include.crispr.arrays <- T}
} else {include.crispr.arrays <- F}


# Search for systems.
debug_msg("Searching for defence systems")
systems <- lapply(X = system_names, FUN = search_system, merged_tbls = merged)
# Rename the data frames generated above.
names(systems) <- system_names
# Merge into a single table
padloc_out <- bind_rows(systems)

# Formatting.
format_output<-function(to_format){
  
  # Formatting.
  formatted <- to_format %>% 
    mutate(
      full.seq.E.value = signif(full.seq.E.value, 3),
      domain.iE.value = signif(domain.iE.value, 3),
      target.description = ifelse(is.na(target.description) == T, target.name, target.description)
    ) %>%
    mutate(target.description = str_remove(target.description, "MULTISPECIES: "))
  
    # if checking for CRISPR_arrays, remove the system type (to prevent users misinterpreting array types)
  
  if(include.crispr.arrays == T){
    formatted <- formatted %>% 
      mutate(
        system = ifelse(protein.name == "CRISPR_array", "CRISPR_array", system),
        cluster = ifelse(protein.name == "CRISPR_array", 0, cluster)
      ) %>% #TODO: depreciate system.number? (requires updating operon display script and webserver)
    distinct()
  }
  

  # Remove "other" systems that overlap canonical systems
  formatted <- formatted %>% 
    mutate(
      system.class = gsub("_.*", "", system),
      is.other = ifelse(grepl("_other", system) == T, 1, 0)
    ) %>%
    group_by(seqid, target.name, system.class) %>% 
    mutate(remove = ifelse(is.other == 1 & min(is.other) == 0, 1, 0)) %>%
    filter(remove == 0) %>%
    group_by(seqid, system, cluster) %>% #TODO: depreciate system.number? (requires updating operon display script and webserver)
    mutate(system.number = cur_group_id()) %>%
    ungroup() %>%
    select(-c(cluster, system.class, is.other, remove)) %>%
    select(system.number, everything()) %>%
    arrange(seqid, relative.position)

  return(formatted)

}

# Output table of defence systems and annotation file
if ( nrow(padloc_out) > 0 ) {
  padloc_out <- format_output(padloc_out)
  msg(paste0("Writing output to '", OUTPUT_DIR, "/", assembly_name, "_padloc.csv'"))
  
  padloc_out %>%
    write_csv(file = paste0(OUTPUT_DIR, "/", assembly_name, "_padloc.csv"))
  
  gff <- generate_gff(padloc_out)
  gff %>% write_delim(file = paste0(OUTPUT_DIR, "/", assembly_name, "_padloc.gff"),
                      delim = "\t", col_names = FALSE)
  
} else {
  msg(paste0("Nothing found for ", assembly_name))
}

# Record end time.
end_time <- Sys.time()
debug_msg(paste0("End time: ", end_time))
run_time <- end_time - start_time
debug_msg(paste0("Run time: ", run_time))




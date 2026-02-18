library(fst)
library(dplyr)
library(httr)
library(jsonlite)
library(purrr)
library(webchem)

# Read the original data
data_path <- "results/combined_dose_data_summarized_no_okl_LFQ_only.fst"
df <- read_fst(data_path)

# Get unique CIDs
unique_cids <- unique(df$compound_pubchem_cid) %>% 
  as.character() %>% 
  na.omit()

cat("Total unique CIDs to query:", length(unique_cids), "\n")

# Function to get PubMed ID counts from PubChem in batches
get_pubmed_counts <- function(cids) {
  # Batch size for PubChem is usually fine up to 400-500
  cid_string <- paste(cids, collapse = ",")
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid_string, "/xrefs/PubMedID/JSON")
  
  res <- tryCatch({
    GET(url)
  }, error = function(e) {
    message("Error querying PubChem: ", e$message)
    return(NULL)
  })
  
  if (is.null(res) || status_code(res) != 200) {
    # If a batch fails, it might be due to a single bad CID or just server error
    return(data.frame(compound_pubchem_cid = as.numeric(cids), pubmed_count = 0))
  }
  
  content_raw <- content(res, "text")
  data_list <- fromJSON(content_raw, simplifyVector = FALSE)
  
  info <- data_list$InformationList$Information
  
  result <- map_dfr(info, function(x) {
    data.frame(
      compound_pubchem_cid = as.numeric(x$CID),
      pubmed_count = length(x$PubMedID)
    )
  })
  
  return(result)
}

# Function to get clinical data from MyChem.info in batches
get_mychem_info <- function(cids) {
  # POST to MyChem.info (up to 1000 per batch)
  # Using cleaner, more quantitative fields from the schema
  fields <- "chembl.max_phase,chembl.withdrawn_flag,aeolus.indications.count,aeolus.no_of_outcomes,gtopdb.approved,unichem.clinicaltrials"
  res <- tryCatch({
    POST("http://mychem.info/v1/query", 
         body = list(q = paste(cids, collapse = ","), 
                    scopes = "pubchem.cid", 
                    fields = fields), 
         encode = "form")
  }, error = function(e) {
    message("Error querying MyChem: ", e$message)
    return(NULL)
  })
  
  if (is.null(res) || status_code(res) != 200) return(NULL)
  
  data_list <- fromJSON(content(res, "text"), simplifyVector = FALSE)
  
  # Process each result
  map_dfr(data_list, function(x) {
    # If not found, x will have 'notfound': TRUE
    cid <- as.numeric(x$query)
    if (!is.null(x$notfound) && x$notfound) {
      return(data.frame(
        compound_pubchem_cid = cid,
        max_phase = 0,
        withdrawn = FALSE,
        indications_count = 0,
        outcomes_count = 0,
        gtopdb_approved = FALSE,
        clinical_trials_count = 0
      ))
    }
    
    # Extract annotations using the specific schema paths
    max_phase <- x$chembl$max_phase %||% 0
    withdrawn <- x$chembl$withdrawn_flag %||% FALSE
    
    # aeolus.indications.count is safer than length() of the indications object
    indications_count <- x$aeolus$indications$count %||% 0
    outcomes_count <- x$aeolus$no_of_outcomes %||% 0
    
    # gtopdb.approved is a very clean boolean marker
    gtopdb_approved <- x$gtopdb$approved %||% FALSE
    
    # Extract clinical trials count (from unichem)
    clinical_trials_count <- 0
    if (!is.null(x$unichem)) {
      trial_fields <- grep("clinicaltrials", names(x$unichem), value = TRUE)
      clinical_trials_count <- length(trial_fields)
    }
    
    data.frame(
      compound_pubchem_cid = cid,
      max_phase = as.numeric(max_phase),
      withdrawn = isTRUE(withdrawn),
      indications_count = as.numeric(indications_count),
      outcomes_count = as.numeric(outcomes_count),
      gtopdb_approved = isTRUE(gtopdb_approved),
      clinical_trials_count = clinical_trials_count
    )
  })
}

# Chunking for batches
chunk_size <- 200
chunks <- split(unique_cids, ceiling(seq_along(unique_cids) / chunk_size))

cat("Querying PubChem for PubMed ID counts...\n")
pubmed_results <- list()
for (i in seq_along(chunks)) {
  cat("PubChem chunk", i, "of", length(chunks), "...\n")
  res <- get_pubmed_counts(chunks[[i]])
  if (!is.null(res)) {
    pubmed_results[[i]] <- res
  }
  Sys.sleep(0.5)
}
pubmed_df <- bind_rows(pubmed_results) %>% 
  distinct(compound_pubchem_cid, .keep_all = TRUE)

cat("Querying MyChem.info for clinical annotations...\n")
# Smaller batches for MyChem can be faster/safer
mychem_results <- list()
# Chunks for MyChem (can be larger, up to 1000)
mychem_chunks <- split(unique_cids, ceiling(seq_along(unique_cids) / 500))
for (i in seq_along(mychem_chunks)) {
  cat("MyChem chunk", i, "of", length(mychem_chunks), "...\n")
  res <- get_mychem_info(mychem_chunks[[i]])
  if (!is.null(res)) {
    mychem_results[[i]] <- res
  }
  Sys.sleep(0.5)
}
mychem_df <- bind_rows(mychem_results) %>% 
  distinct(compound_pubchem_cid, .keep_all = TRUE)

# Join it all together
cat("Joining data...\n")
final_annotations <- full_join(pubmed_df, mychem_df, by = "compound_pubchem_cid")

# Also add the synonym count and clinical suffix logic from before
cat("Retrieving synonyms and mapping clinical suffixes...\n")
# This part is slower if we do it for everyone. We can batch it too.
# But we already have a previous synonym count? Let's re-calculate if needed or join.
# Given time, I'll batch synonym counts as well.
get_synonym_info <- function(cids) {
  syns <- tryCatch({ pc_synonyms(cids, from = "cid") }, error = function(e) return(NULL))
  if (is.null(syns)) return(NULL)
  map_dfr(names(syns), function(cid) {
    s <- syns[[cid]]
    if (length(s) == 1 && is.na(s)) return(data.frame(compound_pubchem_cid = as.numeric(cid), n_synonyms = 0, has_clinical_suffix = FALSE))
    clinical_suffixes <- c("nib$", "mib$", "mab$", "tinib$", "ciclib$", "parib$", "lisib$", "staurin$", "sertib$", "ertib$")
    has_clinical_suffix <- any(sapply(clinical_suffixes, function(p) any(grepl(p, s, ignore.case = TRUE))))
    data.frame(compound_pubchem_cid = as.numeric(cid), n_synonyms = length(unique(s)), has_clinical_suffix = has_clinical_suffix)
  })
}

synonym_results <- list()
for (i in seq_along(chunks)) {
    cat("Synonym chunk", i, "of", length(chunks), "...\n")
    res <- get_synonym_info(chunks[[i]])
    if (!is.null(res)) synonym_results[[i]] <- res
    Sys.sleep(0.5)
}
synonym_df <- bind_rows(synonym_results) %>% distinct(compound_pubchem_cid, .keep_all = TRUE)

final_df <- df %>%
  mutate(compound_pubchem_cid = as.numeric(compound_pubchem_cid)) %>%
  left_join(final_annotations, by = "compound_pubchem_cid") %>%
  left_join(synonym_df, by = "compound_pubchem_cid") %>%
  mutate(
    # Synthesize a final robust "clinical" flag
    is_clinical = (!withdrawn) & (has_clinical_suffix | (max_phase >= 2) | gtopdb_approved | (indications_count > 0))
  )

cat("Saving result...\n")
write_fst(final_df, "results/combined_dose_data_summarized_no_okl_LFQ_only_annotated_v2.fst")
cat("Done! Output saved to results/combined_dose_data_summarized_no_okl_LFQ_only_annotated_v2.fst\n")

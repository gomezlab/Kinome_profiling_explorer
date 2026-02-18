library(fst)
library(dplyr)
library(webchem)
library(purrr)

# Read the original data
data_path <- "results/combined_dose_data_summarized_LFQ_only.fst"
df <- read_fst(data_path)

# Get unique CIDs
unique_cids <- unique(df$compound_pubchem_cid) %>% 
  as.character() %>% 
  na.omit()

# For testing, just use first 5
# unique_cids <- unique_cids[1:5]

cat("Total unique CIDs to query:", length(unique_cids), "\n")

# Function to get synonym info for a batch of CIDs
get_batch_synonyms <- function(cids) {
  # pc_synonyms can handle a vector of CIDs
  syns <- tryCatch({
    pc_synonyms(cids, from = "cid")
  }, error = function(e) {
    message("Error in batch: ", e$message)
    return(NULL)
  })
  
  if (is.null(syns)) return(NULL)
  
  # Process results
  lapply(names(syns), function(cid) {
    s <- syns[[cid]]
    if (length(s) == 1 && is.na(s)) {
      return(data.frame(
        compound_pubchem_cid = as.numeric(cid),
        is_clinical = FALSE,
        n_synonyms = 0
      ))
    }
    
    # Check for DrugBank ID (DB followed by 5 digits)
    has_drugbank <- any(grepl("^DB[0-9]{5}$", s))
    
    # Check for clinical suffixes in any of the synonyms
    # (nib, mib, mab, tinib, ciclib, parib, lisib, etc.)
    clinical_suffixes <- c("nib$", "mib$", "mab$", "tinib$", "ciclib$", "parib$", "lisib$", "staurin$", "sertib$", "ertib$")
    has_clinical_suffix <- any(sapply(clinical_suffixes, function(p) any(grepl(p, s, ignore.case = TRUE))))
    
    data.frame(
      compound_pubchem_cid = as.numeric(cid),
      is_clinical = has_drugbank | has_clinical_suffix,
      n_synonyms = length(unique(s)) # Count unique synonyms
    )
  }) %>% 
    bind_rows() %>% 
    distinct() # Ensure unique CIDs in output
}

# Run in chunks to be safe
chunk_size <- 100
chunks <- split(unique_cids, ceiling(seq_along(unique_cids) / chunk_size))

all_synonym_info <- list()

for (i in seq_along(chunks)) {
  cat("Querying chunk", i, "of", length(chunks), "...\n")
  info <- get_batch_synonyms(chunks[[i]])
  if (!is.null(info)) {
    all_synonym_info[[i]] <- info
  }
  Sys.sleep(0.5) # Avoid hitting API limits
}

synonym_df <- bind_rows(all_synonym_info) %>%
  distinct(compound_pubchem_cid, .keep_all = TRUE)

# Add synonyms info back to the original df
# Ensure same type for join
df_annotated <- df %>% 
  mutate(compound_pubchem_cid = as.numeric(compound_pubchem_cid)) %>%
  left_join(synonym_df, by = "compound_pubchem_cid")

# Save the result
write_fst(df_annotated, "results/combined_dose_data_summarized_LFQ_only_annotated.fst")

cat("Done! Annotated data saved to results/combined_dose_data_summarized_LFQ_only_annotated.fst\n")

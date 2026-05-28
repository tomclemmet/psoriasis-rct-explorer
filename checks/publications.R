#!/usr/bin/env Rscript
# Publication-table quality checks. Run from the project root:
#   Rscript checks/publications.R
#
# Surfaces problems for the curator to fix at source (Access DB):
#
#   1. The same publication recorded twice against one trial. Three signals
#      are checked in order of confidence: identical DOI, identical
#      year/journal/volume/issue/page_start, and identical normalised title.
#      The same publication legitimately appearing under multiple
#      study_ids (a paper that reports more than one trial) is NOT a
#      duplicate -- it is listed separately as informational so each
#      occurrence can be verified.
#
#   2. Suspected text-encoding errors in title / authors / journal:
#      mojibake byte sequences (Windows-1252 read as UTF-8 or vice versa),
#      the Unicode replacement character U+FFFD, undecoded HTML entities,
#      and stray ASCII control characters.

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

find_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(fa)) return(dirname(normalizePath(fa[1], mustWork = FALSE)))
  for (i in seq_along(sys.frames())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of, mustWork = FALSE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                  error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p, mustWork = FALSE)))
  }
  NULL
}

here <- find_script_dir()
sqlite_p <- NULL
for (cand in c(
  if (!is.null(here)) file.path(here, "..", "app", "psoriasis-rcts.sqlite"),
  "app/psoriasis-rcts.sqlite",
  "psoriasis-rcts.sqlite"
)) {
  if (file.exists(cand)) { sqlite_p <- normalizePath(cand); break }
}
if (is.null(sqlite_p))
  stop("psoriasis-rcts.sqlite not found - run app/convert.R first.")

con <- dbConnect(SQLite(), sqlite_p, flags = SQLITE_RO)

options(width = 200)

q <- function(sql) dbGetQuery(con, sql)

trunc_chars <- function(df, max_chars = 70) {
  for (nm in names(df)) {
    if (is.character(df[[nm]])) {
      long <- !is.na(df[[nm]]) & nchar(df[[nm]]) > max_chars
      df[[nm]][long] <- paste0(substr(df[[nm]][long], 1, max_chars - 1), "...")
    }
  }
  df
}

section <- function(title, df, note = NULL, max_rows = Inf) {
  cat("\n--- ", title, " ---\n", sep = "")
  if (!is.null(note)) cat(note, "\n", sep = "")
  if (nrow(df) == 0) {
    cat("OK (0 problems)\n")
  } else {
    cat(sprintf("%d problem row(s):\n", nrow(df)))
    shown <- trunc_chars(if (nrow(df) > max_rows) head(df, max_rows) else df)
    print(shown, row.names = FALSE, right = FALSE)
    if (nrow(df) > max_rows) {
      cat(sprintf("  ... %d more rows not shown.\n", nrow(df) - max_rows))
    }
  }
}

cat("Checking", sqlite_p, "\n")

# ---------------------------------------------------------------------------
# 1. Duplicate references within a single trial
# ---------------------------------------------------------------------------

# 1a. Same DOI used twice (or more) within the same study_id. Highest-
# confidence signal: a DOI uniquely identifies a publication, so any repeat
# inside one study is a duplicate row.
section(
  "Duplicate DOI within a trial",
  q("
    WITH dup AS (
      SELECT study_id, lower(trim(doi)) AS doi_n
      FROM   publications
      WHERE  doi IS NOT NULL AND trim(doi) != ''
      GROUP  BY study_id, lower(trim(doi))
      HAVING COUNT(*) > 1
    )
    SELECT p.study_id, p.publication_id, p.is_primary, p.year,
           p.journal, p.doi, p.title
    FROM   publications p
    JOIN   dup d
      ON   d.study_id = p.study_id
      AND  lower(trim(p.doi)) = d.doi_n
    ORDER BY p.study_id, p.doi, p.is_primary DESC, p.publication_id
  "),
  note = "Each block of rows shares one (study_id, DOI). Keep one; flag the rest."
)

# 1b. Same year/journal/volume/issue/page_start within a study. Catches
# duplicates where the DOI is missing or differently formatted (e.g. with
# or without https://doi.org/ prefix).
section(
  "Duplicate citation tuple (year+journal+volume+issue+page_start) within a trial",
  q("
    WITH dup AS (
      SELECT study_id, year, lower(trim(journal)) AS journal_n,
             trim(volume) AS vol_n, trim(issue) AS iss_n,
             trim(page_start) AS p_n
      FROM   publications
      WHERE  year     IS NOT NULL
        AND  journal  IS NOT NULL AND trim(journal)    != ''
        AND  volume   IS NOT NULL AND trim(volume)     != ''
        AND  page_start IS NOT NULL AND trim(page_start) != ''
      GROUP  BY study_id, year, journal_n, vol_n, iss_n, p_n
      HAVING COUNT(*) > 1
    )
    SELECT p.study_id, p.publication_id, p.is_primary, p.year, p.journal,
           p.volume, p.issue, p.page_start, p.doi, p.title
    FROM   publications p
    JOIN   dup d
      ON   d.study_id = p.study_id
      AND  d.year = p.year
      AND  d.journal_n = lower(trim(p.journal))
      AND  d.vol_n = trim(p.volume)
      AND  d.iss_n = trim(p.issue)
      AND  d.p_n   = trim(p.page_start)
    ORDER BY p.study_id, p.year, p.journal, p.publication_id
  "),
  note = "Catches duplicates where DOI is missing or inconsistently formatted."
)

# 1c. Same normalised title within a study. Lower confidence than 1a/1b
# because near-identical titles can be a preprint + final publication, an
# erratum, or a supplement; eyeball before flagging in Access.
pubs <- q("SELECT publication_id, study_id, is_primary, year, journal,
                  title, doi FROM publications")

norm_title <- function(s) {
  s <- ifelse(is.na(s), "", s)
  s <- tolower(s)
  s <- gsub("[^a-z0-9]+", " ", s, perl = TRUE)
  trimws(gsub("\\s+", " ", s))
}
pubs$title_n <- norm_title(pubs$title)
has_title <- nzchar(pubs$title_n)
key <- paste(pubs$study_id, pubs$title_n, sep = "|")
dup_keys <- unique(key[has_title & duplicated(key[has_title])])
title_dups <- pubs[has_title & key %in% dup_keys,
                   c("study_id", "publication_id", "is_primary",
                     "year", "journal", "doi", "title")]
title_dups <- title_dups[order(title_dups$study_id, title_dups$title,
                               -title_dups$is_primary,
                               title_dups$publication_id), ]
section(
  "Duplicate normalised title within a trial",
  title_dups,
  note = paste(
    "Titles compared after lowercasing and stripping punctuation/whitespace.",
    "Lower confidence: review each group -- could be erratum, supplement,",
    "or preprint vs final paper."
  )
)

# 1d. Informational: same DOI appears under multiple study_ids. Expected
# when one paper reports more than one trial; printed so each occurrence
# can be confirmed intentional rather than a study_id mix-up.
section(
  "[info] Same DOI shared across multiple trials (likely legitimate)",
  q("
    WITH shared AS (
      SELECT lower(trim(doi)) AS doi_n
      FROM   publications
      WHERE  doi IS NOT NULL AND trim(doi) != ''
      GROUP  BY lower(trim(doi))
      HAVING COUNT(DISTINCT study_id) > 1
    )
    SELECT p.study_id, p.publication_id, p.is_primary, p.year,
           p.journal, p.doi, p.title
    FROM   publications p
    JOIN   shared s ON lower(trim(p.doi)) = s.doi_n
    ORDER BY p.doi, p.study_id, p.publication_id
  "),
  note = "Not flagged as an error; verify each block reflects a paper that genuinely reports multiple trials."
)

# ---------------------------------------------------------------------------
# 2. Suspected encoding errors
# ---------------------------------------------------------------------------

# Mojibake byte signature: when Windows-1252 bytes are mis-decoded as UTF-8
# (the dominant cause in Access -> CSV -> SQLite pipelines), the original
# byte becomes one of {0xC3, 0xC2, 0xE2} -- which print as the Latin-1
# characters U+00C3 (Ã), U+00C2 (Â), U+00E2 (â) -- followed by a byte in
# U+0080..U+00BF. We match that pair using explicit Unicode escapes so the
# script's own encoding cannot mangle the pattern.
mojibake_re <- "[\u00C3\u00C2\u00E2][\u0080-\u00BF]|\u00E2\u20AC"
suspect_re <- paste(
  mojibake_re,
  "�",                        # Unicode REPLACEMENT CHARACTER (lost byte)
  "&[a-zA-Z]{2,8};",               # named HTML entity (&amp;, &nbsp;, ...)
  "&#[0-9]{1,5};",                 # numeric HTML entity
  "&#x[0-9a-fA-F]{1,5};",          # hex HTML entity
  "[\x01-\x08\x0B\x0C\x0E-\x1F]",  # control chars excluding TAB/LF/CR
  sep = "|"
)

# Return the first regex match per element, NA where there is none.
# (regmatches() drops non-matching elements -- can't use it directly here
# because we need to keep alignment with the input vector.)
first_match <- function(x, re) {
  m <- regexpr(re, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok  <- m > 0
  out[ok] <- substring(x[ok], m[ok],
                       m[ok] + attr(m, "match.length")[ok] - 1L)
  out
}

pubs2 <- q("SELECT publication_id, study_id, is_primary, year,
                   title, authors, journal FROM publications")

find_encoding_issues <- function(df, fields) {
  out_rows <- list()
  for (fld in fields) {
    val <- df[[fld]] %||% rep(NA_character_, nrow(df))
    val[is.na(val)] <- ""
    hits <- first_match(val, suspect_re)
    flag <- !is.na(hits)
    if (any(flag)) {
      out_rows[[fld]] <- data.frame(
        publication_id = df$publication_id[flag],
        study_id       = df$study_id[flag],
        field          = fld,
        pattern_hit    = hits[flag],
        value          = val[flag],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out_rows)) return(data.frame(
    publication_id = integer(), study_id = integer(), field = character(),
    pattern_hit = character(), value = character(),
    stringsAsFactors = FALSE))
  do.call(rbind, out_rows)
}

enc_issues <- find_encoding_issues(pubs2, c("title", "authors", "journal"))
if (nrow(enc_issues)) {
  enc_issues <- enc_issues[order(enc_issues$field, enc_issues$study_id,
                                 enc_issues$publication_id), ]
}
section(
  "Suspected encoding errors in title / authors / journal",
  enc_issues,
  note = paste(
    "`pattern_hit` shows the first suspicious sequence found in `value`.",
    "Mojibake pairs, the U+FFFD replacement character, undecoded HTML",
    "entities, and stray control characters are all flagged."
  )
)

cat("\nDone.\n")

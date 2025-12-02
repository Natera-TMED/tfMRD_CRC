# ---- Packages ----
# install.packages(c("readxl","dplyr","lubridate","survival","broom"))
library(readxl)
library(dplyr)
library(lubridate)
library(survival)
library(broom)

rm(list=ls()) # start clean
setwd("~/Documents/Clinical Data Mining/CRC/Galaxy/tfMRD/")

# ---- Inputs ----
in_path <- "PostDefinitiveTherapy_HR.xlsx"

# If TRUE: time between Post.Definitive.Start.Date and first Collection.date is treated as ctDNA-negative.
# If FALSE: that pre-test time is excluded from the model.
assume_negative_before_first <- TRUE

# ---- Helpers ----
clean_names_like_janitor <- function(x) {
  y <- tolower(trimws(x))
  y <- gsub("[^a-z0-9]+", ".", y)
  y <- gsub("\\.+", ".", y)
  y <- gsub("^\\.|\\.$", "", y)
  y
}

pick_first_present <- function(nms, candidates) {
  cand <- intersect(candidates, nms)
  if (length(cand) > 0) cand[1] else NA_character_
}

# Map cleaned column names to required roles (now uses post-definitive start as origin)
infer_schema <- function(nms) {
  list(
    id        = pick_first_present(nms, c("id","patient.id","patientid","mrn","subject","patient")),
    # NEW ORIGIN:
    postdef_start = pick_first_present(
      nms,
      c(
        "post.definitive.start.date","postdefinitive.start.date","postdefinitive.start",
        "post.definitive.therapy.start.date","postdef.tx.start.date","postdef.start",
        "postdefinitive.tx.start","post.definitive.tx.start","post.tx.start.date"
      )
    ),
    # keep surgery optional (unused) just for diagnostics if present
    surgery   = pick_first_present(nms, c("surgery.date","surgerydate","date.of.surgery","dos")),
    collect   = pick_first_present(nms, c("collection.date","collect.date","draw.date","ctdna.collection.date","ctdna.date")),
    ctdna     = pick_first_present(nms, c("ctdna","ctdna.result","ctdna.status","ctdna.binary","ctdna.posneg","ct.dna")),
    rfs_date  = pick_first_present(nms, c("rfs.date","rfsdate","relapse.free.survival.date","recurrence.free.survival.date")),
    # recognize plain "rfs" as the event/result column
    rfs_result= pick_first_present(
      nms,
      c("rfs","rfs.result","rfs.event","rfs.status","rfs.outcome","rfs.binary",
        "relapse","relapse.result","relapse.event","relapse.status",
        "event","status"))
  )
}

# ---- Read & standardize columns ----
raw <- read_excel(in_path)
orig_names <- names(raw)
cleaned <- clean_names_like_janitor(orig_names)
names(raw) <- cleaned

schema <- infer_schema(cleaned)
missing <- names(schema)[vapply(schema, function(x) is.na(x) || !nzchar(x), logical(1))]
# Require id, postdef_start, collect, ctdna, rfs_date, rfs_result (collect/ctdna not strictly required to run, but helpful)
required <- c("id","postdef_start","rfs_date","rfs_result")
missing_required <- intersect(required, missing)
if (length(missing_required) > 0) {
  stop(sprintf(
    paste0(
      "Missing required column(s): %s\n",
      "Cleaned names available:\n  %s\n",
      "Tip: Ensure columns like ID, Post.Definitive.Start.Date, Collection.date, ctDNA, RFS.Date, and RFS are present."
    ),
    paste(missing_required, collapse = ", "),
    paste(cleaned, collapse = ", ")
  ))
}

# ---- Build base dataset ----
dat <- raw %>%
  transmute(
    id            = as.character(.data[[schema$id]]),
    origin        = as_date(.data[[schema$postdef_start]]),   # CHANGED: origin is post-definitive start
    collect       = suppressWarnings(as_date(.data[[schema$collect]])),
    ctdna_raw     = .data[[schema$ctdna]],
    rfs_date      = as_date(.data[[schema$rfs_date]]),
    rfs_event_raw = .data[[schema$rfs_result]])

# Robust recode for ctDNA to 0/1 (1 = Positive/Detected)
to01_ctdna <- function(x) {
  if (is.numeric(x) || is.logical(x)) {
    return(ifelse(is.na(x), NA_integer_, as.integer(as.numeric(x) > 0)))
  }
  xs <- tolower(trimws(as.character(x)))
  xs <- gsub("\\s+", "", xs)
  pos <- grepl("^(pos|positive|det|detected|true|yes|y|1)$", xs)
  neg <- grepl("^(neg|negative|und|undetected|false|no|n|0)$", xs)
  out <- ifelse(pos, 1L, ifelse(neg, 0L, NA_integer_))
  if (any(is.na(out))) {
    suppressWarnings({
      xn <- as.numeric(xs)
      out[is.na(out) & !is.na(xn)] <- as.integer(xn[is.na(out) & !is.na(xn)] > 0)
    })
  }
  out
}

# Robust recode for RFS event to 0/1 (1 = relapsed/event)
to01_event <- function(x) {
  if (is.numeric(x) || is.logical(x)) {
    return(ifelse(is.na(x), NA_integer_, as.integer(as.numeric(x) > 0)))
  }
  xs <- tolower(trimws(as.character(x)))
  xs <- gsub("\\s+", "", xs)
  ev1 <- grepl("(relaps|recur|event|progress|yes|true|1)", xs)
  ev0 <- grepl("^(no|none|censor|free|false|0)$", xs)
  out <- ifelse(ev1 & !ev0, 1L, ifelse(ev0 & !ev1, 0L, NA_integer_))
  if (any(is.na(out))) {
    suppressWarnings({
      xn <- as.numeric(xs)
      out[is.na(out) & !is.na(xn)] <- as.integer(xn[is.na(out) & !is.na(xn)] > 0)
    })
  }
  out
}

dat <- dat %>%
  mutate(
    ctdna     = to01_ctdna(ctdna_raw),
    rfs_event = to01_event(rfs_event_raw)
  )

# Keep essential dates/fields and ensure draws inside [origin, rfs_date]
dat <- dat %>%
  filter(!is.na(id), !is.na(origin), !is.na(rfs_date)) %>%
  filter(is.na(collect) | (collect >= origin & collect <= rfs_date))

# Choose one origin per patient (earliest post-definitive start)
origin_ref <- dat %>%
  group_by(id) %>%
  summarize(origin = min(origin, na.rm = TRUE), .groups = "drop")

# One RFS per patient: latest date; event = max (prefers 1 if any row shows event)
rfs_ref <- dat %>%
  group_by(id) %>%
  summarize(
    rfs_date  = max(rfs_date, na.rm = TRUE),
    rfs_event = max(rfs_event, na.rm = TRUE),
    .groups = "drop"
  )

# All ctDNA tests per patient (dedupe same-day draws; keep last)
tests <- dat %>%
  select(id, collect, ctdna) %>%
  filter(!is.na(collect), !is.na(ctdna)) %>%
  arrange(id, collect) %>%
  group_by(id, collect) %>%
  slice_tail(n = 1) %>%
  ungroup()

# Master frame for interval construction
frame <- origin_ref %>%
  inner_join(rfs_ref, by = "id") %>%
  left_join(tests, by = "id", relationship = "one-to-many") %>%
  filter(is.na(collect) | (collect >= origin & collect <= rfs_date)) %>%
  arrange(id, collect)

# Build patient-specific intervals carrying ctDNA forward (origin = Post.Definitive.Start.Date)
build_intervals_one <- function(df) {
  id        <- df$id[1]
  origin    <- df$origin[1]
  rfs_date  <- df$rfs_date[1]
  rfs_event <- df$rfs_event[1]
  
  timeline <- df %>%
    select(collect, ctdna) %>%
    distinct() %>%
    arrange(collect)
  
  rows <- list()
  
  # Pre-first-test interval (from ORIGIN to first collect)
  if (nrow(timeline) > 0) {
    first_collect <- timeline$collect[1]
    if (!is.na(first_collect) && first_collect > origin && assume_negative_before_first) {
      rows[[length(rows) + 1]] <- tibble(
        id = id,
        start_date = origin,
        stop_date  = min(first_collect, rfs_date),
        ctdna      = 0L,
        event      = as.integer(rfs_event == 1 && min(first_collect, rfs_date) == rfs_date)
      )
    }
  } else {
    # No tests at all
    if (assume_negative_before_first) {
      rows[[length(rows) + 1]] <- tibble(
        id = id,
        start_date = origin,
        stop_date  = rfs_date,
        ctdna      = 0L,
        event      = as.integer(rfs_event == 1)
      )
    } else {
      return(tibble(id = character(), start_date = as_date(character()),
                    stop_date = as_date(character()), ctdna = integer(), event = integer()))
    }
  }
  
  # Intervals from each test date to next test (or RFS), bounded by origin/rfs_date
  if (nrow(timeline) > 0) {
    for (i in seq_len(nrow(timeline))) {
      this_date <- timeline$collect[i]
      this_val  <- timeline$ctdna[i]
      next_date <- if (i < nrow(timeline)) timeline$collect[i + 1] else rfs_date
      
      s <- max(this_date, origin)
      e <- min(next_date, rfs_date)
      if (!is.na(s) && !is.na(e) && e > s) {
        rows[[length(rows) + 1]] <- tibble(
          id = id,
          start_date = s,
          stop_date  = e,
          ctdna      = as.integer(this_val),
          event      = as.integer(rfs_event == 1 && e == rfs_date)
        )
      }
    }
  }
  
  bind_rows(rows)
}

intervals <- frame %>%
  group_by(id) %>%
  group_modify(~build_intervals_one(.x)) %>%
  ungroup()

# Counting-process times relative to Post.Definitive.Start.Date
intervals <- intervals %>%
  left_join(origin_ref, by = "id") %>%
  mutate(
    start = as.numeric(difftime(start_date, origin, units = "days")),
    stop  = as.numeric(difftime(stop_date , origin, units = "days"))
  ) %>%
  filter(stop > start) %>%
  select(id, start, stop, event, ctdna)

# ---- Fit Cox model ----
fit <- coxph(Surv(start, stop, event) ~ ctdna + cluster(id),
             data = intervals, ties = "efron")

cat("\nCox model with time-dependent ctDNA (origin = Post.Definitive.Start.Date, end = RFS.Date)\n")
print(summary(fit))

hr <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::filter(term == "ctdna") %>%
  select(term, estimate, conf.low, conf.high, p.value)

cat("\nHazard Ratio for ctDNA-positive vs ctDNA-negative:\n")
print(hr)

# ---- Diagnostics ----
cat("\n[Auto-mapped columns]\n")
cat(sprintf("ID                         : %s\n", schema$id))
cat(sprintf("Post.Definitive.Start.Date: %s\n", schema$postdef_start))
cat(sprintf("Surgery.Date (unused)     : %s\n", schema$surgery))
cat(sprintf("Collection                : %s\n", schema$collect))
cat(sprintf("ctDNA                     : %s\n", schema$ctdna))
cat(sprintf("RFS.Date                  : %s\n", schema$rfs_date))
cat(sprintf("RFS (event)               : %s\n\n", schema$rfs_result))

cat("[Interval counts per patient]\n")
print(intervals %>% count(id, name = "n_intervals") %>% arrange(desc(n_intervals)) %>% head(10))

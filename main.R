library(seqinr)
library(ggplot2)
library(dplyr)
source("hmm_for_r.txt")

logSumExp <- function(log_probs) {
  max_log <- max(log_probs)
  if (is.infinite(max_log)) return(-Inf)
  return(max_log + log(sum(exp(log_probs - max_log))))
}

check_domain_detection <- function(results_df) {
  results_df <- results_df %>%
    group_by(Sequence) %>%
    mutate(ExpectedDomains = strsplit(Sequence, "_"),
           Rank = rank(-LogOdds, ties.method = "first")) %>%
    ungroup()
  
  # Napravi novu tablicu gdje je TRUE ako je predviđena domena među top očekivanima
  accuracy_check <- results_df %>%
    rowwise() %>%
    mutate(IsExpected = Domain %in% unlist(ExpectedDomains),
           CorrectTop = IsExpected && Rank <= length(unlist(ExpectedDomains))) %>%
    ungroup()
  
  # Rezime po sekvenci
  summary_per_seq <- accuracy_check %>%
    group_by(Sequence) %>%
    summarise(
      Expected = paste(sort(unique(unlist(ExpectedDomains))), collapse = ","),
      TopPredicted = paste(Domain[order(-LogOdds)][1:length(unlist(ExpectedDomains))], collapse = ","),
      Matches = sum(CorrectTop),
      TotalExpected = length(unlist(ExpectedDomains)),
      AllCorrect = Matches == TotalExpected
    )
  
  print(summary_per_seq)
  cat("\n✅ Točnost po sekvenci:", sum(summary_per_seq$AllCorrect), "/", nrow(summary_per_seq), "sekvenci su imale sve točne HMM domene u top skoru.\n")
  
  return(summary_per_seq)
}

simulate_until_end <- function(hmm, start_state = NULL, end_states = NULL) {
  if (is.null(start_state)) {
    start_state <- sample(hmm$States, size = 1, prob = hmm$startProbs)
  }
  if (is.null(end_states)) {
    end_states <- c("E")
  }
  
  states <- character()
  observations <- character()
  current_state <- start_state
  
  repeat {
    states <- c(states, current_state)
    
    obs_probs <- hmm$emissionProbs[current_state, ]
    if (all(obs_probs == 0)) {
      # Stanje bez emisija - ne emitiraj simbol, samo nastavi
      observations <- c(observations, "")  # možeš staviti NA ili ""
    } else {
      obs <- sample(hmm$Symbols, size = 1, prob = obs_probs)
      observations <- c(observations, obs)
    }
    
    if (current_state %in% end_states) {
      break
    }
    
    trans_probs <- hmm$transProbs[current_state, ]
    if (all(trans_probs == 0)) {
      warning("Nema prijelaznih vjerojatnosti iz stanja ", current_state, ". Simulacija završena ranije.")
      break
    }
    
    current_state <- sample(hmm$States, size = 1, prob = trans_probs)
  }
  
  list(states = states, observations = observations)
}

forward_custom <- function(hmm, obs_seq) {
  States <- hmm$States
  Symbols <- hmm$Symbols
  transProbs <- hmm$transProbs
  emissionProbs <- hmm$emissionProbs
  
  nStates <- length(States)
  T_len <- length(obs_seq)
  
  state_idx <- setNames(seq_along(States), States)
  symbol_idx <- setNames(seq_along(Symbols), Symbols)
  obs_idx <- sapply(obs_seq, function(s) symbol_idx[[s]])
  
  # Smoothing: zamijeni 0 s vrlo malim vrijednostima, pa logiraj
  emissionProbs[emissionProbs == 0] <- 1e-10
  emissionProbs <- emissionProbs / rowSums(emissionProbs)
  transProbs[transProbs == 0] <- 1e-10
  transProbs <- transProbs / rowSums(transProbs)

  logTrans <- log(transProbs)
  logEmiss <- log(emissionProbs)
  
  logAlpha <- matrix(-Inf, nrow = nStates, ncol = T_len + 1)
  rownames(logAlpha) <- States
  
  # Inicijalizacija
  start_state <- "S"
  start_idx <- state_idx[[start_state]]
  logAlpha[start_idx, 1] <- 0  # log(1)
  
  # Pomocna funkcija za propagaciju "null" stanja (bez emisije)
  propagate_null_states <- function(logAlpha_col) {
    converged <- FALSE
    iter <- 0
    max_iter <- 1000
    
    while(!converged) {
      iter <- iter + 1
      prev <- logAlpha_col
      
      for (j in seq_len(nStates)) {
        state_j <- States[j]
        emits <- any(emissionProbs[j, ] > 1e-9)  # sigurno ne emitira

        if (!emits) {
          incoming <- logAlpha_col + logTrans[, j]
          incoming[is.nan(incoming)] <- -Inf
          logAlpha_col[j] <- logSumExp(incoming)
        }
      }
      
      diff <- abs(logAlpha_col - prev)
      diff[is.na(diff)] <- 0
      converged <- all(diff < 1e-10) || iter >= max_iter
      
      if (iter >= max_iter) {
        warning("propagate_null_states nije konvergirao unutar ", max_iter, " iteracija")
        break
      }
    }
    
    return(logAlpha_col)
  }
  
  # Propagacija kroz null-stanja na početku
  logAlpha[, 1] <- propagate_null_states(logAlpha[, 1])
  
  # Glavni koraci kroz promatranu sekvencu
  for (t in seq_len(T_len)) {
    logAlpha_t <- rep(-Inf, nStates)
    
    for (j in seq_len(nStates)) {
      emits <- any(emissionProbs[j, ] > 1e-9)
      
      if (emits) {
        e_prob <- logEmiss[j, obs_idx[t]]
        incoming <- logAlpha[, t] + logTrans[, j]
        incoming[is.nan(incoming)] <- -Inf
        logAlpha_t[j] <- logSumExp(incoming) + e_prob
      }
    }
    
    logAlpha[, t + 1] <- propagate_null_states(logAlpha_t)
  }
  
  # Vraćamo cijelu matricu logAlpha ako korisniku treba detaljnije
  return(list(logAlpha = logAlpha))
}

forward_bg <- function(hmm, obs_seq) {
  n_states <- length(hmm$States)
  T_len <- length(obs_seq)
  
  # log-forward matrica (redovi: stanja, stupci: vrijeme)
  logAlpha <- matrix(-Inf, nrow = n_states, ncol = T_len)
  rownames(logAlpha) <- hmm$States
  
  # Početne log vjerojatnosti (startProbs + emisija prve promatrane vrijednosti)
  for (s in 1:n_states) {
    emisija_prob <- hmm$emissionProbs[s, obs_seq[1]]
    logAlpha[s, 1] <- log(hmm$startProbs[s]) + log(emisija_prob)
  }
  
  # Forward iteracija
  for (t in 2:T_len) {
    for (s in 1:n_states) {
      emisija_prob <- hmm$emissionProbs[s, obs_seq[t]]
      logSum <- -Inf
      for (sp in 1:n_states) {
        trans_prob <- hmm$transProbs[sp, s]
        logSum <- logSumExp(c(logSum, logAlpha[sp, t-1] + log(trans_prob)))
      }
      logAlpha[s, t] <- log(emisija_prob) + logSum
    }
  }
  
  return(logAlpha)
}

# Funkcija za generiranje sekvenci i određivanje duljine prozora
determine_window_size <- function(hmm, n = 50) {
  lengths <- numeric(n)
  for (i in 1:n) {
    sim <- simulate_until_end(hmm, start_state = "S", end_states = c("E"))
    obs <- sim$observations[sim$observations != ""]
    lengths[i] <- length(obs)
  }
  round(mean(lengths))
}

# Učitaj sve modele iz .GlobalEnv
all_domains <- grep("^transProbs_", ls(), value = TRUE)
all_domains <- sub("transProbs_", "", all_domains)

# Napravi HMM-ove kao liste
hmm_list <- list()
window_sizes <- numeric()
symbols_global <- get(paste0("symbols_", all_domains[1]))

for (domain in all_domains) {
  hmm <- list(
    States = get(paste0("states_", domain)),
    Symbols = get(paste0("symbols_", domain)),
    startProbs = get(paste0("startProbs_", domain)),
    transProbs = get(paste0("transProbs_", domain)),
    emissionProbs = get(paste0("emissionProbs_", domain))
  )
    # ✅ Postavi dimnames da radi indeksiranje po imenima
  rownames(hmm$transProbs) <- hmm$States
  colnames(hmm$transProbs) <- hmm$States
  rownames(hmm$emissionProbs) <- hmm$States
  colnames(hmm$emissionProbs) <- hmm$Symbols
  hmm_list[[domain]] <- hmm
  window_sizes[domain] <- determine_window_size(hmm)
}
#definicija pozadinskog HMM-a
states_bg <- c("B")

emissionProbs_bg <- matrix(1 / length(symbols_global), 
                           nrow = 1, ncol = length(symbols_global),
                           dimnames = list(states_bg, symbols_global))

transProbs_bg <- matrix(1, nrow = 1, ncol = 1,
                        dimnames = list(states_bg, states_bg))

hmm_bg <- list(
  States = states_bg,
  Symbols = symbols_global,
  startProbs = c(B = 1),
  transProbs = transProbs_bg,
  emissionProbs = emissionProbs_bg
)

# Kalkulacija log-odds
calc_log_prob <- function(hmm, obs_seq, use_custom = TRUE) {
  if (use_custom) {
    res <- forward_custom(hmm, obs_seq)
    if ("E" %in% hmm$States) {
      return(res$logAlpha["E", length(obs_seq) + 1])
    } else {
      return(logSumExp(res$logAlpha[, length(obs_seq) + 1]))
    }
  } else {
    res <- forward_bg(hmm, obs_seq)
    return(logSumExp(res[, ncol(res)]))
  }
}

# Procesiraj sve fasta datoteke
fasta_files <- list.files("test_sequences", pattern = "\\.fa$", full.names = TRUE)
results <- data.frame()

for (fasta_file in fasta_files) {
  seq_name <- tools::file_path_sans_ext(basename(fasta_file))
  cat("▶️ Obrada sekvence:", seq_name, "\n")
  
  seqs <- read.fasta(fasta_file, seqtype = "AA", as.string = FALSE)
  
  for (i in seq_along(seqs)) {
    aa_seq <- toupper(seqs[[i]])
    
    for (domain in all_domains) {
      hmm <- hmm_list[[domain]]
      w_size <- window_sizes[domain]
      max_log_odds <- -Inf
      
      cat("  ➡️ HMM:", domain, "- prozor:", w_size, "- duljina sekvence:", length(aa_seq), "\n")
      
      if (length(aa_seq) < w_size) {
        cat(sprintf("    ⚠️ Sekvenca prekratka (%d) za prozor (%d) → koristi se cijela\n", length(aa_seq), w_size))
        
        window_seq <- aa_seq
        if (!all(window_seq %in% hmm$Symbols)) next
        
        log_p_hmm <- calc_log_prob(hmm, window_seq, TRUE)
        log_p_bg <- calc_log_prob(hmm_bg, window_seq, FALSE)
        log_odds <- log_p_hmm - log_p_bg
        
        if (!is.na(log_odds) && is.finite(log_odds)) {
          max_log_odds <- log_odds
        }
        
      } else {
        for (start in seq(1, length(aa_seq) - w_size + 1, by = floor(w_size / 5))) {
          cat(sprintf("    📍 Prozor na poziciji %d\n", start))
          
          window_seq <- aa_seq[start:(start + w_size - 1)]
          if (!all(window_seq %in% hmm$Symbols)) next
          
          log_p_hmm <- calc_log_prob(hmm, window_seq, TRUE)
          log_p_bg <- calc_log_prob(hmm_bg, window_seq, FALSE)
          log_odds <- log_p_hmm - log_p_bg
          
          if (!is.na(log_odds) && is.finite(log_odds)) {
            max_log_odds <- max(max_log_odds, log_odds)
          }
        }
      }
      
      results <- rbind(results, data.frame(
        Sequence = seq_name,
        Domain = domain,
        LogOdds = max_log_odds
      ))
    }
  }
}

# Pozovi funkciju
domain_eval_summary <- check_domain_detection(results)
# Vizualizacija
# Definiraj do 12 različitih oblika (možeš ih proširiti ako ima više domena)
available_domains <- unique(results$Domain)
n_domains <- length(available_domains)
shape_vals <- c(0:5, 15:20, 21:25)[1:n_domains]  # do 12 oblika

ggplot(results, aes(x = Sequence, y = LogOdds, color = Domain, shape = Domain)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 2) +
  scale_shape_manual(values = shape_vals) +
  labs(
    title = "Log-Odds vrijednosti po sekvenci i HMM-u",
    x = "Naziv FASTA sekvence", y = "MAX log-odds score",
    color = "HMM domena", shape = "HMM domena"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),  # makni horizontalne
    panel.grid.major.x = element_line(linetype = "solid", color = "grey60", size = 0.3)  # tanke sive vertikale
  )




#pokretanje simulacija
#Rscript main.R
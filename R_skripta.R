
logSumExp <- function(log_probs) {
  max_log <- max(log_probs)
  if (is.infinite(max_log)) return(-Inf)
  return(max_log + log(sum(exp(log_probs - max_log))))
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

source("hmm_for_r.txt")

# Smoothing emisijske matrice
emissionProbs_smoothed <- emissionProbs  # kopija originala

for (i in 1:nrow(emissionProbs_smoothed)) {
  if (all(emissionProbs_smoothed[i, ] == 0)) {
    # Neemitirajuće stanje (npr. D, S, E) – preskoči
    next
  }
  emissionProbs_smoothed[i, emissionProbs_smoothed[i, ] == 0] <- 1e-10
  emissionProbs_smoothed[i, ] <- emissionProbs_smoothed[i, ] / sum(emissionProbs_smoothed[i, ])
}

# Smoothing prijelazne matrice
transProbs_smoothed <- transProbs  # kopija originala

for (i in 1:nrow(transProbs_smoothed)) {
  if (all(transProbs_smoothed[i, ] == 0)) {
    # Stanje bez prijelaza – preskoči
    next
  }
  transProbs_smoothed[i, transProbs_smoothed[i, ] == 0] <- 1e-10
  transProbs_smoothed[i, ] <- transProbs_smoothed[i, ] / sum(transProbs_smoothed[i, ])
}

library(HMM)

# Inicijaliziraj HMM
hmm <- initHMM(States = states, Symbols = symbols,
               startProbs = c(ifelse(states == "S", 1, 0)),
               transProbs = transProbs,
               emissionProbs = emissionProbs)

# Simuliraj sekvencu dok ne dođeš do stanja "E"
gen <- simulate_until_end(hmm, start_state = "S", end_states = c("E"))

# Ispiši sekvencu stanja
cat("Stanja:\n")
print(gen$states)

# Filtriraj promatrane simbole
obs_seq <- gen$observations[!is.na(gen$observations) & gen$observations != ""]
cat("Simboli:\n")
cat(paste(obs_seq, collapse = ""), "\n")

# Izračunaj log-vjerojatnost sekvence
res <- forward_custom(hmm, obs_seq)
log_prob <- res$logAlpha["E", length(obs_seq) + 1]
cat("Log vjerojatnost:\n")
print(log_prob)

cat("Obična vjerojatnost:\n")
print(exp(log_prob))

print("Razultati za netočnu sekvencu:")
#ubaci netočnu sekvencu MAIMFVPVGEQIKPIYEGFKYVKNIEKVYLIASNKTERYAHDIKKRIEYIYDTEIVLVDPEKLDDIMEKLIDVVSENKDR
#EIISNITGGTKVMSHACYILCSYLGGDAFYIFKRDDGSMEYVDMPMLKIKLNSVIEDRSTRERILEKLMESEYDSMTSLA
#RELRIKDSTLSVILDQLKEQGLVILERSGRNLRIKISKTGKILLRLRKLKK
wrong_seq <- c("M", "A", "I", "M", "F", "V", "P", "V", "G", "E", "Q", "I", "K", "P", "I", "Y", "E", "G", "F", "K",
              "Y", "V", "K", "N", "I", "E", "K", "V", "Y", "L", "I", "A", "S", "N", "K", "T", "E",
              "R", "Y", "A", "H", "D", "I", "K")

# Izračunaj log-vjerojatnost sekvence
res <- forward_custom(hmm, wrong_seq)
log_prob <- res$logAlpha["E", length(wrong_seq) + 1]
cat("Log vjerojatnost netočne sekvence:\n")
print(log_prob)

cat("Obična vjerojatnost netočne sekvence:\n")
print(exp(log_prob))
#----------------------------------------------------------------
#---------------------------------------------------------------
# Definiraj pozadinski trivijalni model: jedno stanje, emisije prema frekvencijama simbola u svim sekvencama ili uniformno

# Naivni pozadinski model: jedno stanje B, emitira simbole s uniformnom vjerojatnošću
states_bg <- c("B")
symbols_bg <- symbols

transProbs_bg <- matrix(1, nrow=1, ncol=1)  # uvijek ostaje u B

emissionProbs_bg <- matrix(1/length(symbols_bg), nrow=1, ncol=length(symbols_bg))
colnames(emissionProbs_bg) <- symbols_bg
rownames(emissionProbs_bg) <- states_bg

# Postavljanje dimnames je OBAVEZNO za ugrađenu forward funkciju iz HMM paketa
rownames(transProbs_bg) <- states_bg
colnames(transProbs_bg) <- states_bg

rownames(emissionProbs_bg) <- states_bg
colnames(emissionProbs_bg) <- symbols_bg
#print("Pozadinski model ima matrice emisije i prijelaza:")
#print(emissionProbs_bg)
#print("------------------")
#print(transProbs_bg)

hmm_bg <- list(
  States = states_bg,
  Symbols = symbols_bg,
  startProbs = c(B = 1),
  transProbs = transProbs_bg,
  emissionProbs = emissionProbs_bg
)


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

# Funkcija koja računa log-vjerojatnost sekvence za neki HMM
calc_log_prob <- function(hmm, obs_seq, use_builtin_bg = FALSE) {
  if (use_builtin_bg) {
    # koristi forward_bg
    logAlpha <- forward_bg(hmm, obs_seq)
    # Zadnji stupac sadrži log-vjerojatnosti po stanjima za cijelu sekvencu
    # Suma preko svih stanja u log prostoru daje ukupnu vjerojatnost sekvence
    return(logSumExp(logAlpha[, ncol(logAlpha)]))
  } else {
    # koristi vlastiti custom forward (koji vraća logAlpha s +1 dužinom za "E" stanje)
    res <- forward_custom(hmm, obs_seq)
    if ("E" %in% hmm$States) {
      return(res$logAlpha["E", length(obs_seq) + 1])
    } else {
      return(logSumExp(res$logAlpha[, length(obs_seq) + 1]))
    }
  }
}

# Izračun log-vjerojatnosti za generiranu sekvencu
log_prob_hmm <- calc_log_prob(hmm, obs_seq, use_builtin_bg = FALSE)
#print("do ovde")
log_prob_bg <- calc_log_prob(hmm_bg, obs_seq, use_builtin_bg = TRUE)
log_odds <- log_prob_hmm - log_prob_bg

cat("Log-vjerojatnost modela:\n")
print(log_prob_hmm)
cat("Log-vjerojatnost pozadinskog modela:\n")
print(log_prob_bg)
cat("Log-odds score za generiranu sekvencu:\n")
print(log_odds)


# Isto za netočnu sekvencu
wrong_seq <- c("M", "A", "I", "M", "F", "V", "P", "V", "G", "E", "Q", "I", "K", "P", "I", "Y", "E", "G", "F", "K",
              "Y", "V", "K", "N", "I", "E", "K", "V", "Y", "L", "I", "A", "S", "N", "K", "T", "E",
              "R", "Y", "A", "H", "D", "I", "K")
#wrong_seq <- c("M", "A", "I", "M")
log_prob_wrong_hmm <- calc_log_prob(hmm, wrong_seq, use_builtin_bg = FALSE)
log_prob_wrong_bg <- calc_log_prob(hmm_bg, wrong_seq, use_builtin_bg = TRUE)
log_odds_wrong <- log_prob_wrong_hmm - log_prob_wrong_bg

cat("\nNetocna sekvenca:\n")
cat("Log-vjerojatnost modela za netočnu sekvencu:\n")
print(log_prob_wrong_hmm)
cat("Log-vjerojatnost pozadinskog modela za netočnu sekvencu:\n")
print(log_prob_wrong_bg)
cat("Log-odds score za netočnu sekvencu:\n")
print(log_odds_wrong)

# Zaključi klasifikaciju po log-odds
if (log_odds > log_odds_wrong) {
  cat("\nModel bolje objašnjava generiranu sekvencu nego netočnu.\n")
} else {
  cat("\nNetocna sekvenca ima viši log-odds score, treba dodatno istražiti!\n")
}
#-----------------------------------------------------------------
#----------------------------------------------------------------
# Ako je sekvenca nemoguća, istraži zašto
if (!is.finite(log_prob)) {
  cat("⚠️ Sekvenca NIJE moguća prema HMM-u (log-vjerojatnost = -Inf).\n")

  # Provjera 1: može li se doći do E iz S
  cat("\n🔍 Provjera prijelaza od S do E:\n")
  reachable <- function(transMat, states, from, to, visited = character()) {
    if (from == to) return(TRUE)
    if (from %in% visited) return(FALSE)
    visited <- c(visited, from)
    idx_from <- which(states == from)
    for (i in seq_along(states)) {
      if (transMat[idx_from, i] > 0) {
        if (reachable(transMat, states, states[i], to, visited)) return(TRUE)
      }
    }
    return(FALSE)
  }

  path_exists <- reachable(hmm$transProbs, hmm$States, "S", "E")
  if (path_exists) {
    cat("✅ Put od 'S' do 'E' postoji.\n")
  } else {
    cat("❌ Put od 'S' do 'E' NE postoji! Provjeri prijelaznu matricu.\n")
  }

  # Provjera 2: simboli imaju barem jedno moguće stanje koje ih emitira
  cat("\n🔍 Provjera emisija za svaki simbol u sekvenci:\n")
  for (sym in obs_seq) {
    emitters <- which(hmm$emissionProbs[, sym] > 0)
    if (length(emitters) == 0) {
      cat("❌ Nijedno stanje ne emitira simbol '", sym, "'\n", sep = "")
    } else {
      cat("✅ Simbol '", sym, "' emitiraju stanja: ", paste(hmm$States[emitters], collapse = ", "), "\n", sep = "")
    }
  }
}



#Rscript R_skripta.R
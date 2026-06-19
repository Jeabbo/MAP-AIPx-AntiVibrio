# =========================
# KD + helixvis wheel + Eisenberg bubble + left footer + AMP verdict
# =========================

# 
pkgs <- c("ggplot2","patchwork","zoo","helixvis","grid","Peptides")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
if (length(to_install)) install.packages(to_install)
library(ggplot2); library(patchwork); library(zoo); library(helixvis); library(grid); library(Peptides)

# ---- User input ----
SEQ     <- "SLLWRWNSMKPVGASCRDHLECGTKYCRRSICTF"  
WINDOW  <- 9                              


# --------
residue_class <- function(a){
  if(a %in% c('D','E')) return('Acidic (D,E)')
  if(a %in% c('R','K','H')) return('Basic (R,K,H)')
  if(a %in% c('F','W','Y')) return('Aromatic (F,W,Y)')
  if(a %in% c('A','I','L','V','M')) return('Aliphatic (A,I,L,V,M)')
  if(a %in% c('S','T','N','Q')) return('Polar (S,T,N,Q)')
  if(a %in% c('C')) return('Cysteine (C)')
  if(a %in% c('G','P')) return('Others (G,P)')
  'Unknown'
}
pal <- c('Acidic (D,E)'='#F4A7B9','Basic (R,K,H)'='#4E9DE7','Aromatic (F,W,Y)'='#98D27F',
         'Aliphatic (A,I,L,V,M)'='#CFCFCF','Polar (S,T,N,Q)'='#A7E0B2',
         'Cysteine (C)'='#FFD93D','Others (G,P)'='#333333','Unknown'='#999999')

# ---- Scales ----
KD_raw <- c(I=4.5,V=4.2,L=3.8,F=2.8,C=2.5,M=1.9,A=1.8,G=-0.4,T=-0.7,S=-0.8,
            W=-0.9,Y=-1.3,P=-1.6,H=-3.2,E=-3.5,Q=-3.5,D=-3.5,N=-3.5,K=-3.9,R=-4.5)
kd_rescale <- function(x) (x+4.5)/9
KD <- kd_rescale(KD_raw)

# Eisenberg hydrophobicity for μH
EIS <- c(A=0.62,R=-2.53,N=-0.78,D=-0.90,C=0.29,Q=-0.85,E=-0.74,G=0.48,
         H=-0.40,I=1.38,L=1.06,K=-1.50,M=0.64,F=1.19,P=0.12,S=-0.18,
         T=-0.05,W=0.81,Y=0.26,V=1.08)

# Wimley–White interfacial hydrophobicity (kcal/mol)
WWIHP <- c(
  A = 0.17, R = 0.81, N = 0.42, D = 1.23, C = -0.24, Q = 0.58, E = 2.02, G = 0.01,
  H = 0.96, I = -1.12, L = -1.25, K = 0.99, M = -0.67, F = -1.71, P = 0.45, S = 0.13,
  T = 0.14, W = -2.09, Y = -0.71, V = -0.46
)

# ---- Chou–Fasman (compact %SS) ----
CF <- data.frame(
  aa=c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"),
  Pa=c(1.42,0.98,0.67,1.01,0.70,1.11,1.51,0.57,1.00,1.08,1.21,1.16,1.45,1.13,0.57,0.77,0.83,1.08,0.69,1.06),
  Pb=c(0.83,0.93,0.89,0.54,1.19,1.10,0.37,0.75,0.87,1.60,1.30,0.74,1.05,1.38,0.55,0.75,1.19,1.37,1.47,1.70),
  Pc=c(0.66,0.95,1.56,1.46,1.19,0.98,1.01,1.64,0.95,0.47,0.59,1.01,0.60,0.60,1.52,1.43,0.96,0.77,0.88,0.50),
  row.names=c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
)

# ---- Utils ----
clean_seq <- function(seq){
  s <- toupper(gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", seq))
  if(nchar(s)==0) stop("Sequence is empty after cleaning.")
  strsplit(s, "", fixed=TRUE)[[1]]
}
sliding_mean_centered <- function(x, w){
  stopifnot(w %% 2 == 1, w >= 1)
  n <- length(x); h <- (w-1)/2; out <- rep(NA_real_, n)
  for(i in (1+h):(n-h)) out[i] <- mean(x[(i-h):(i+h)])
  out
}
build_df <- function(seq, scale_vec){
  aa <- clean_seq(seq); vals <- as.numeric(scale_vec[aa])
  data.frame(position=seq_along(aa), residue=aa, value=vals,
             class=vapply(aa, residue_class, character(1)))
}
cf_percent <- function(seq){
  aa <- clean_seq(seq); P <- CF[aa, c("Pa","Pb","Pc")]
  pr <- 100 * colSums(P) / sum(colSums(P))
  setNames(round(pr,1), c("Helix","Sheet","Coil"))
}
kd_scores <- function(seq){
  aa <- clean_seq(seq); vals_raw <- KD_raw[aa]
  list(mean_raw = mean(vals_raw),
       mean_norm_0_1 = kd_rescale(mean(vals_raw)))
}
wwihp_mean <- function(seq){
  aa <- clean_seq(seq); round(mean(WWIHP[aa], na.rm=TRUE), 2)
}

# ---- KD panel (y-axis smaller & not bold) ----
kd_panel <- function(seq, window=WINDOW){
  df <- build_df(seq, KD); df$slide <- sliding_mean_centered(df$value, window)
  ggplot(df, aes(factor(position), value, fill=class)) +
    geom_col(color="grey25", width=0.9) +
    geom_line(aes(y=slide, group=1), color="black", linewidth=0.7, na.rm=TRUE) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_fill_manual(values=pal) +
    scale_x_discrete(labels=df$residue) +
    labs(title=sprintf("Hydropathy Profile (Kyte–Doolittle, window=%d)", window),
         y="Hydropathy Index (KD, normalized 0–1)", x=NULL, fill="Residue Class") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 9, face = "plain"),
          axis.text.y = element_text(size = 7,  face = "plain"),
          axis.title.y = element_text(size = 8, face = "plain"))
}

# ---- μH sliding & best window ----
muH_sliding <- function(seq, window=11, angle_deg=100){
  aa <- clean_seq(seq); h <- as.numeric(EIS[aa]); n <- length(h)
  if (window > n || window %% 2 == 0) return(data.frame(pos=1:n, muH=rep(NA_real_,n)))
  half <- floor(window/2); theta <- angle_deg*pi/180; out <- rep(NA_real_, n)
  for (i in (1+half):(n-half)){
    w <- h[(i-half):(i+half)]; k <- -half:half
    mx <- sum(w*cos(k*theta)); my <- sum(w*sin(k*theta))
    out[i] <- sqrt(mx^2 + my^2)/window
  }
  data.frame(pos=1:n, muH=out)
}
auto_muH_best <- function(seq){
  aa <- clean_seq(seq); n <- length(aa)
  cand <- c(9,11,13); cand <- cand[cand <= n]
  if (!length(cand)) cand <- ifelse(n %% 2 == 1, n, n-1)
  best <- list(win=NA, pos=NA, max=-Inf)
  for (w in cand){
    df <- muH_sliding(seq, window=w); m <- max(df$muH, na.rm=TRUE); p <- which.max(df$muH)
    if (is.finite(m) && m > best$max) best <- list(win=w, pos=p, max=round(m,3))
  }
  best
}

# ---- Eisenberg bubble wheel ----
wedge_poly <- function(phi, wdeg, r_in=0, r_out=1.45, n=120){
  start <- phi-(wdeg*pi/180)/2; end <- phi+(wdeg*pi/180)/2
  th_outer <- seq(start,end,length.out=n); th_inner <- rev(th_outer)
  xo <- r_out*cos(th_outer); yo <- r_out*sin(th_outer)
  xi <- r_in*cos(th_inner);  yi <- r_in*sin(th_inner)
  data.frame(x=c(xo,xi), y=c(yo,yi))
}
bubble_wheel <- function(seq, start, end, angle_deg=100, face_width_deg=100){
  aa <- clean_seq(seq); h <- as.numeric(EIS[aa])
  w_aa <- aa[start:end]; w_h <- h[start:end]; k <- 0:(length(w_aa)-1)
  theta <- angle_deg*pi/180; ang <- (k*theta)%%(2*pi)
  x <- cos(ang); y <- sin(ang)
  mx <- sum(w_h*cos(ang)); my <- sum(w_h*sin(ang))
  muH <- sqrt(mx^2+my^2)/length(w_aa); dir <- atan2(my,mx)
  ax <- 1.25*muH*cos(dir); ay <- 1.25*muH*sin(dir)
  df <- data.frame(x,y,aa=w_aa,h=w_h,abs_h=abs(w_h))
  wedge_hydrophob <- wedge_poly(dir,face_width_deg)
  wedge_hydrophil <- wedge_poly(dir+pi,face_width_deg)

  ggplot(df,aes(x,y))+
    geom_polygon(data=wedge_hydrophob,aes(x,y),inherit.aes=FALSE, fill="#ffa500",alpha=.12)+
    geom_polygon(data=wedge_hydrophil,aes(x,y),inherit.aes=FALSE, fill="#2c7bb6",alpha=.08)+
    geom_hline(yintercept=0,linetype="dotted")+
    geom_vline(xintercept=0,linetype="dotted")+
    geom_point(aes(size=abs_h,color=h),alpha=.95) +
    scale_color_gradient2(low="#2c7bb6", mid="grey85", high="#d7191c", midpoint=0,
      name="Hydrophobicity h\n(+ hydrophobic, − hydrophilic)") +
    scale_size(range=c(3,8), guide="none") +
    geom_text(aes(label=aa),vjust=-1.1,size=3.6,fontface="bold")+
    geom_segment(x=0,y=0,xend=ax,yend=ay, arrow=arrow(length=grid::unit(.25,"cm")),
                 linewidth=1.1,color="black")+
    coord_fixed(xlim=c(-1.65,1.65),ylim=c(-1.75,1.75))+
    theme_void()+
    theme(plot.title   = element_text(size=9,  face="plain"),
          legend.title = element_text(size=8,  face="plain"),
          legend.text  = element_text(size=8,  face="plain"),
          legend.position="right",
          plot.margin=margin(0,5,0,0))+
    ggtitle(sprintf("Eisenberg bubble wheel (%d–%d) | μH = %.3f", start,end,muH))
}

# ---- helixvis wheel ----
helixvis_wheel <- function(seq, start, end){
  subseq <- paste0(clean_seq(seq)[start:end], collapse = "")
  helixvis::draw_wheel(sequence = subseq, labels = TRUE) +
    theme(
      plot.title   = element_text(size = 9, face = "plain"),
      plot.margin  = margin(0,10,0,0),
      legend.title = element_text(size = 9, face = "plain"),
      legend.text  = element_text(size = 8, face = "plain")
    ) +
    ggtitle(sprintf("helixvis wheel (%d–%d)", start, end))
}

# --------
HRC <- function(seq){
  hyd <- c("A","V","L","I","M","F","W","P")
  aa <- clean_seq(seq); round(100*mean(aa %in% hyd),2)
}
SAI <- function(seq){
  aa <- clean_seq(seq); h <- mean(KD_raw[aa])
  net <- sum(aa %in% c("R","K","H")) - sum(aa %in% c("D","E"))
  round(h/(abs(net)+1),3)
}
HelixProp <- function(seq){
  tab <- c(A=1.42,R=0.98,N=0.67,D=1.01,C=0.70,E=1.51,Q=1.11,G=0.57,H=1.00,
           I=1.08,L=1.21,K=1.16,M=1.45,F=1.13,P=0.57,S=0.77,T=0.83,W=1.08,
           Y=0.69,V=1.06)
  aa <- clean_seq(seq); round(mean(tab[aa]),2)
}

# ---- AMP verdict helper ----
amp_verdict <- function(charge, pI, muH, KD_raw_mean, instability,
                        helix_propensity, HRC, Boman, WW_mean = NA) {

  passed <- c(
    charge >= 2,
    pI >= 8.5,
    muH >= 0.45,
    KD_raw_mean >= 0.0 & KD_raw_mean <= 1.2,
    instability < 40,
    helix_propensity >= 1
  )

  # treat |WW| <= 0.2 as neutral; only flag clearly positive WW
  risk_flags <- c(
    Boman_high = !is.na(Boman) && Boman > 3.0,
    HRC_high   = !is.na(HRC)   && HRC > 65,
    WW_pos     = !is.na(WW_mean) && WW_mean > 0.2
  )

  pass_count <- sum(passed, na.rm = TRUE)

  if (pass_count >= 4 && !any(risk_flags)) {
    decision <- sprintf("MAP Verdict: Membrane-Interaction Plausible ✅ (passed %d/6 criteria)", pass_count)
    color <- "darkgreen"
  } else if (pass_count >= 4 && any(risk_flags)) {
    reason <- paste(names(risk_flags)[risk_flags], collapse = ", ")
    decision <- sprintf("MAP Verdict: Membrane-Interaction Borderline ⚠️ (passed %d/6 criteria; risk: %s)", pass_count, reason)
    color <- "orange"
  } else if (pass_count >= 3) {
    decision <- sprintf("MAP Verdict: Membrane-Interaction Borderline ⚠️ (passed %d/6 criteria)", pass_count)
    color <- "orange"
  } else {
    decision <- sprintf("MAP Verdict: Low Membrane-Interaction Plausibility ❌ (passed %d/6 criteria)", pass_count)
    color <- "red3"
  }

  list(text = decision, color = color)
}


aipx_score <- function(seq){

  # --- Reference means and SDs (adjust if needed) ---
  ref_mean <- list(muH=0.50, Charge=0.15, GRAVY=0.10, HelixProp=1.05, HRC=45.0, SAI=0.20)
  ref_sd   <- list(muH=0.15, Charge=0.10, GRAVY=0.25, HelixProp=0.20, HRC=15.0, SAI=0.15)

 
  aa <- clean_seq(seq)
  ks <- kd_scores(seq)
  muH <- auto_muH_best(seq)$max
  charge_density <- (sum(aa %in% c("R","K","H")) - sum(aa %in% c("D","E"))) / length(aa)
  GRAVY <- ks$mean_raw
  HelixP <- HelixProp(seq)
  HRCv <- HRC(seq)
  SAIv <- SAI(seq)

  
  z <- function(x, m, s) (x - m) / (ifelse(s == 0, 1e-6, s))
  z_muH  <- z(muH, ref_mean$muH, ref_sd$muH)
  z_chg  <- z(charge_density, ref_mean$Charge, ref_sd$Charge)
  z_grav <- z(GRAVY, ref_mean$GRAVY, ref_sd$GRAVY)
  z_hx   <- z(HelixP, ref_mean$HelixProp, ref_sd$HelixProp)
  z_hrc  <- z(HRCv, ref_mean$HRC, ref_sd$HRC)
  z_sai  <- z(SAIv, ref_mean$SAI, ref_sd$SAI)

  
  coefs <- c(muH=0.25, Charge=0.20, GRAVY=0.15, HelixProp=0.15, HRC=0.15, SAI=0.10)
  zsum <- sum(coefs * c(z_muH, z_chg, z_grav, z_hx, z_hrc, z_sai))

 
  AIPx <- 1 / (1 + exp(-zsum))

  
  class <- if (AIPx > 0.60) {
    "Strong AMP-like ✅"
  } else if (AIPx >= 0.45) {
    "Borderline ⚠️"
  } else {
    "Unlikely AMP ❌"
  }

  # --- Timestamp for logging ---
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  
  cat(sprintf(
    "\n[%s] Raw Z-sum for %s\n  μH_z=%.2f  Charge_z=%.2f  GRAVY_z=%.2f  Helix_z=%.2f  HRC_z=%.2f  SAI_z=%.2f\n  → Z-sum=%.3f   AIPx=%.3f   → %s\n",
    timestamp, seq, z_muH, z_chg, z_grav, z_hx, z_hrc, z_sai, zsum, AIPx, class
  ))

  # --- Log to file with timestamp ---
  log_line <- sprintf("[%s]\t%s\tZ-sum=%.3f\tAIPx=%.3f\t%s\n",
                      timestamp, seq, zsum, AIPx, class)
  write(log_line, file="AIPx_log.txt", append=TRUE)

  
  round(AIPx, 3)
}



AIPx_value <- aipx_score(SEQ)
cat(sprintf("\nAIPx score (KD uses KD mean 0–1): %.3f\n", AIPx_value))



# ---- Footer (left-aligned) ----
footer_panel <- function(seq, mu_best){
  aa <- clean_seq(seq)
  pct <- cf_percent(seq)
  ks  <- kd_scores(seq)
  nB  <- sum(aa %in% c("R","K","H"))
  nA  <- sum(aa %in% c("D","E"))
  net <- nB - nA

  
  boman <- tryCatch(Peptides::boman(seq),      error=function(e) NA)
  inst  <- tryCatch(Peptides::instaIndex(seq), error=function(e) NA)
  chg7  <- tryCatch(Peptides::charge(seq, pH=7), error=function(e) NA)
  pIval <- tryCatch(Peptides::pI(seq),           error=function(e) NA)
  hmom  <- tryCatch(Peptides::hmoment(seq, angle=100, window=mu_best$win),
                    error=function(e) NA)

  ww    <- wwihp_mean(seq)

  msg1 <- sprintf("Residues  Basic: %d  Acidic: %d   Net: %+d    Charge(pH7): %.2f    pI: %.2f ", nB, nA, net, chg7, pIval)
  msg2 <- sprintf(
  "KD (GRAVY, raw): %.2f   KD mean (0–1): %.2f   μH max: %.2f (window %d)    AIPx: %.3f",
  ks$mean_raw, ks$mean_norm_0_1, mu_best$max, mu_best$win, aipx_score(seq))
  msg3 <- sprintf("Helix: %.1f%%   Sheet: %.1f%%   Coil: %.1f%%   |   HRC: %.2f%%   SAI: %.3f   Helix_Propensity: %.2f",
                  pct['Helix'], pct['Sheet'], pct['Coil'], HRC(seq), SAI(seq), HelixProp(seq))
  msg4 <- sprintf("Boman: %.2f   Instability: %.2f   Hmoment: %.3f   WWIHP(mean): %.2f",
                  boman, inst, hmom, ww)

  verdict <- amp_verdict(
    charge = net,
    pI = pIval,
    muH = mu_best$max,
    KD_raw_mean = ks$mean_raw,
    instability = inst,
    helix_propensity = HelixProp(seq),
    HRC = HRC(seq),
    Boman = boman,
    WW_mean = ww
  )

  ggplot() +
    annotate("text", x=0, y=1.00, hjust=0, label=msg1, size=3.2, fontface="bold", colour="#1F77B4") +
    annotate("text", x=0, y=0.75, hjust=0, label=msg2, size=3.2, fontface="bold", colour="#1F77B4") +
    annotate("text", x=0, y=0.51, hjust=0, label=msg3, size=3.2, fontface="bold", colour="#1F77B4") +
    annotate("text", x=0, y=0.265, hjust=0, label=msg4, size=3.2, fontface="bold", colour="#1F77B4") +
    annotate("text", x=0, y=0.025, hjust=0, label=verdict$text, size=3, fontface="plain", colour=verdict$color) +
    xlim(0,1) + ylim(0,1) + theme_void()
}

# ---- Compose ----
plot_amp <- function(seq=SEQ, window=WINDOW){
  p_kd <- kd_panel(seq, window)
  best <- auto_muH_best(seq)

  n  <- length(clean_seq(seq)); h <- (best$win-1)/2
  st <- max(1, best$pos - h); en <- min(n, best$pos + h)

  p_left  <- helixvis_wheel(seq, st, en)
  p_right <- bubble_wheel(seq, st, en)

  wheels_row <- (p_left | p_right) + plot_layout(widths = c(1,1))
  p_footer   <- footer_panel(seq, best)

  (p_kd / wheels_row / p_footer) + plot_layout(heights = c(1.0, 1.8, 0.68))
}

# ---- Console summary ----

cat(sprintf("\nAIPx score (z-weighted, sigmoid 0–1): %.3f\n", aipx_score(SEQ)))

# ---- Run ----
p <- plot_amp()
print(p)
# ggsave("AMP_KD_wheels_footer_verdict.png", p, width=12, height=8, dpi=300)

#
##   AIPx (0–1)	Meaning 
##    > 0.75	    Strong AMP
###   0.50 – 0.75	  Moderate 
###   < 0.45	    Weak AMP 

##  Z-sum ≈ 0 → Balanced / borderline
##  Z-sum > 0.8 → Strong AMP-like
##  Z-sum < -0.3 → Weak or non-AMP






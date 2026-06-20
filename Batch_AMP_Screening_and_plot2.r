# ============================================================
# AMP screens from FASTA + Figure per AMP-like sequence (≥6/9)
# KD bars (top) + helixvis wheel (left) + μH bubble (right)
# Footer: KD/μH/SS + HRC/SAI/HelixProp + Boman/Instab/pI/Charge/Hmoment/WWIHP
# Output: AMP_summary.csv + AMP_full_*.png  (12 x 8 in @ 300 dpi)
# Clean, academic styling
# ============================================================

# ---------- User: set FASTA path ----------
FASTA_FILE <- "peptides.fasta"   # <-- put your multi-sequence FASTA here

# ---------- Packages ----------
need <- c("ggplot2","patchwork","zoo","helixvis","Peptides","seqinr")
to_install <- need[!sapply(need, requireNamespace, quietly=TRUE)]
if(length(to_install)) install.packages(to_install)
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(zoo)
  library(helixvis); library(Peptides); library(seqinr)
})

# ---------- Constants & helper tables ----------
# KD (raw) -> for GRAVY; KD_rescaled 0..1 for bars to match Peptalyzer-style
KD_raw <- c(I=4.5,V=4.2,L=3.8,F=2.8,C=2.5,M=1.9,A=1.8,G=-0.4,T=-0.7,S=-0.8,
            W=-0.9,Y=-1.3,P=-1.6,H=-3.2,E=-3.5,Q=-3.5,D=-3.5,N=-3.5,K=-3.9,R=-4.5)
kd_rescale <- function(x) (x + 4.5) / 9
KD <- kd_rescale(KD_raw)

# Eisenberg scale (h) for μH/bubble
EIS <- c(A=0.62,R=-2.53,N=-0.78,D=-0.90,C=0.29,Q=-0.85,E=-0.74,G=0.48,
         H=-0.40,I=1.38,L=1.06,K=-1.50,M=0.64,F=1.19,P=0.12,S=-0.18,
         T=-0.05,W=0.81,Y=0.26,V=1.08)

# Wimley–White interfacial (kcal/mol), peptide→interface (helical/neutral)
# (used to form a simple "WWIHP" sum; more negative ≈ better membrane transfer)
WWIHP <- c(
  A=0.17, R=2.58, N=0.42, D=1.23, C=-0.24, Q=0.58, E=2.02, G=0.01, H=0.96,
  I=-0.31, L=-0.56, K=0.99, M=-0.23, F=-1.13, P=0.45, S=0.13, T=0.14,
  W=-1.85, Y=-0.94, V=-0.07
)

# Chou–Fasman compact (%)
CF <- data.frame(
  aa=c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"),
  Pa=c(1.42,0.98,0.67,1.01,0.70,1.11,1.51,0.57,1.00,1.08,1.21,1.16,1.45,1.13,0.57,0.77,0.83,1.08,0.69,1.06),
  Pb=c(0.83,0.93,0.89,0.54,1.19,1.10,0.37,0.75,0.87,1.60,1.30,0.74,1.05,1.38,0.55,0.75,1.19,1.37,1.47,1.70),
  Pc=c(0.66,0.95,1.56,1.46,1.19,0.98,1.01,1.64,0.95,0.47,0.59,1.01,0.60,0.60,1.52,1.43,0.96,0.77,0.88,0.50),
  row.names=c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
)

# KD bar classes
residue_class <- function(a){
  if(a %in% c('D','E')) return('Acidic (D,E)')
  if(a %in% c('R','K','H')) return('Basic (R,K,H)')
  if(a %in% c('F','W','Y')) return('Aromatic (F,W,Y)')
  if(a %in% c('A','I','L','V','M')) return('Aliphatic (A,I,L,V,M)')
  if(a %in% c('S','T','N','Q')) return('Polar (S,T,N,Q)')
  if(a == 'C') return('Cysteine (C)')
  if(a %in% c('G','P')) return('Others (G,P)')
  'Unknown'
}
pal <- c('Acidic (D,E)'='#F4A7B9','Basic (R,K,H)'='#4E9DE7','Aromatic (F,W,Y)'='#98D27F',
         'Aliphatic (A,I,L,V,M)'='#CFCFCF','Polar (S,T,N,Q)'='#A7E0B2',
         'Cysteine (C)'='#FFD93D','Others (G,P)'='#333333','Unknown'='#999999')

# ---------- Utilities ----------
clean_seq <- function(s){
  x <- toupper(gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", s))
  if(nchar(x)==0) stop("Empty sequence after cleaning.")
  unlist(strsplit(x,""))
}
sliding_mean_centered <- function(x, w){
  stopifnot(w %% 2 == 1, w>0); n <- length(x); h <- (w-1)/2
  out <- rep(NA_real_, n)
  for(i in (1+h):(n-h)) out[i] <- mean(x[(i-h):(i+h)])
  out
}
cf_percent <- function(seq){
  aa <- clean_seq(seq); P <- CF[aa, c("Pa","Pb","Pc"), drop=FALSE]
  pr <- 100 * colSums(P) / sum(colSums(P))
  setNames(round(pr,1), c("Helix","Sheet","Coil"))
}

# μH scan (Eisenberg) and best window among {9,11,13}
muH_sliding <- function(seq, window=11, angle_deg=100){
  aa <- clean_seq(seq); h <- as.numeric(EIS[aa]); n <- length(h)
  if(window > n || window %% 2 == 0) return(rep(NA_real_,n))
  half <- floor(window/2); theta <- angle_deg*pi/180; out <- rep(NA_real_, n)
  for(i in (1+half):(n-half)){
    w <- h[(i-half):(i+half)]; k <- -half:half
    mx <- sum(w*cos(k*theta)); my <- sum(w*sin(k*theta))
    out[i] <- sqrt(mx^2+my^2)/window
  }
  out
}
best_muH <- function(seq){
  aa <- clean_seq(seq); n <- length(aa)
  wins <- c(9,11,13); wins <- wins[wins<=n]
  if(!length(wins)) wins <- ifelse(n %% 2==1, n, n-1)
  best <- list(win=wins[1], pos=NA, muH=-Inf)
  for(w in wins){
    x <- muH_sliding(seq, w)
    m <- suppressWarnings(max(x, na.rm=TRUE))
    if(is.finite(m) && m>best$muH) best <- list(win=w, pos=which.max(x), muH=m)
  }
  best$muH <- round(best$muH,3); best
}

# ---------- Panels ----------
kd_panel <- function(seq, window=9){
  aa <- clean_seq(seq); vals <- KD[aa]; df <- data.frame(
    pos=seq_along(aa), aa=aa, val=vals,
    class=vapply(aa, residue_class, character(1))
  )
  df$trend <- sliding_mean_centered(df$val, window)
  ggplot(df, aes(factor(pos), val, fill=class))+
    geom_col(color="grey25", width=0.9)+
    geom_line(aes(y=trend, group=1), color="black", linewidth=0.7, na.rm=TRUE)+
    geom_hline(yintercept=0, linetype="dashed")+
    scale_fill_manual(values=pal)+
    scale_x_discrete(labels=df$aa)+
    labs(title=sprintf("Hydropathy Profile (Kyte–Doolittle, window=%d)", window),
         y="Hydropathy Index (KD, normalized 0–1)", x=NULL, fill="Residue Class")+
    theme_bw()+
    theme(axis.text.x=element_text(size=10,face="bold"),
          axis.title.y=element_text(size=9,face="plain"),
          axis.text.y =element_text(size=8,face="plain"),
          legend.title=element_text(size=9,face="bold"),
          legend.text =element_text(size=8,face="plain"))
}

# Bubble wheel (uses Eisenberg h; color legend only)
wedge_poly <- function(phi, wdeg, r_in=0, r_out=1.45, n=120){
  start <- phi-(wdeg*pi/180)/2; end <- phi+(wdeg*pi/180)/2
  th <- seq(start,end,length.out=n)
  data.frame(x=c(r_out*cos(th), r_in*cos(rev(th))),
             y=c(r_out*sin(th), r_in*sin(rev(th))))
}
bubble_wheel <- function(seq, start, end, angle_deg=100, face_width_deg=100, title=""){
  aa <- clean_seq(seq); seg <- aa[start:end]; h <- as.numeric(EIS[seg])
  k <- 0:(length(seg)-1); theta <- angle_deg*pi/180; ang <- (k*theta)%%(2*pi)
  x <- cos(ang); y <- sin(ang)
  mx <- sum(h*cos(ang)); my <- sum(h*sin(ang))
  muH <- sqrt(mx^2+my^2)/length(seg); dir <- atan2(my,mx)
  ax <- 1.25*muH*cos(dir); ay <- 1.25*muH*sin(dir)
  df <- data.frame(x,y,aa=seg,h=h,abs_h=abs(h))
  hyd <- wedge_poly(dir,face_width_deg); phil <- wedge_poly(dir+pi, face_width_deg)

  ggplot(df,aes(x,y))+
    geom_polygon(data=hyd,aes(x,y),inherit.aes=FALSE, fill="#ffa500",alpha=.12)+
    geom_polygon(data=phil,aes(x,y),inherit.aes=FALSE, fill="#2c7bb6",alpha=.08)+
    geom_hline(yintercept=0,linetype="dotted")+
    geom_vline(xintercept=0,linetype="dotted")+
    geom_point(aes(color=h, size=abs_h), alpha=.95)+
    geom_text(aes(label=aa), vjust=-1.05, size=3.6, fontface="bold")+
    geom_segment(x=0,y=0,xend=ax,yend=ay,
                 arrow=arrow(length=grid::unit(.25,"cm")), linewidth=1.1)+
    scale_color_gradient2(low="#2c7bb6", mid="grey85", high="#d7191c",
                          midpoint=0, name="Hydrophobicity h\n(+ hydrophobic, − hydrophilic)",
                          guide=guide_colorbar(title.position = "top", title.hjust=0.5))+
    scale_size(range=c(3,8), guide="none")+
    coord_fixed(xlim=c(-1.65,1.65), ylim=c(-1.75,1.75))+
    labs(title=title)+
    theme_void()+
    theme(plot.title=element_text(size=10,face="bold"),
          legend.text =element_text(size=8,face="plain"),
          legend.title=element_text(size=9,face="plain"))
}

helixvis_wheel <- function(seq, start, end, title=""){
  seg <- paste0(clean_seq(seq)[start:end], collapse="")
  helixvis::draw_wheel(seg, labels = TRUE)+
    labs(title=title)+
    theme(plot.title=element_text(size=10,face="bold"))
}

# ---------- Descriptors ----------
HRC <- function(seq){ hyd <- c("A","V","L","I","M","F","W","P")
  aa <- clean_seq(seq); round(100*mean(aa %in% hyd),2) }
SAI <- function(seq){ aa <- clean_seq(seq); h <- mean(KD_raw[aa])
  net <- sum(aa %in% c("R","K","H")) - sum(aa %in% c("D","E"))
  round(h/(abs(net)+1),3) }
HelixProp <- function(seq){
  tab <- c(A=1.42,R=0.98,N=0.67,D=1.01,C=0.70,E=1.51,Q=1.11,G=0.57,H=1.00,
           I=1.08,L=1.21,K=1.16,M=1.45,F=1.13,P=0.57,S=0.77,T=0.83,W=1.08,
           Y=0.69,V=1.06)
  aa <- clean_seq(seq); round(mean(tab[aa]),2)
}
WW_sum <- function(seq){ aa <- clean_seq(seq); round(sum(WWIHP[aa]),3) }

# ---------- AMP scoring (9 rules) ----------
score_rules <- function(seq){
  aa <- paste0(clean_seq(seq), collapse="")
  len <- nchar(aa)
  # core metrics
  mu <- best_muH(aa)
  gravy <- mean(KD_raw[clean_seq(aa)])
  pct <- cf_percent(aa)
  hr  <- HRC(aa)
  ww  <- WW_sum(aa)
  net <- Peptides::charge(aa, pH=7)
  hm  <- tryCatch(Peptides::hmoment(aa, angle=100, window=11), error=function(e) NA_real_)
  bom <- Peptides::boman(aa)
  inst<- Peptides::instaIndex(aa)
  pi  <- Peptides::pI(aa)

  rules <- c(
    between_len      = len >= 12 & len <= 40,
    charge_pos       = net >= 2,
    muH_min          = is.finite(mu$muH) & mu$muH >= 0.40,
    HRC_mid          = hr >= 35 & hr <= 65,
    kd_moderate      = gravy >= -1.5 & gravy <= 1.5,
    ww_favorable     = ww <= -1.0,
    helix_bias       = pct["Helix"] >= 25,
    boman_not_low    = bom >= 0.5,
    hmoment_ok       = is.finite(hm) & hm >= 0.30
  )

  list(
    len=len, mu_best=mu, gravy=gravy, pct=pct, HRC=hr, WW=ww, net=net,
    hmoment=hm, boman=bom, instab=inst, pI=pi,
    rules=rules, passed=sum(rules), verdict=ifelse(sum(rules)>=6,"AMP-like","Unlikely")
  )
}

# ---------- Footer panel ----------
footer_panel <- function(seq, met){
  pct <- met$pct; mu <- met$mu_best
  txt1 <- sprintf("Residues  Basic: %d   Acidic: %d   Net: %+d",
                  sum(clean_seq(seq) %in% c("R","K","H")),
                  sum(clean_seq(seq) %in% c("D","E")),
                  met$net)
  txt2 <- sprintf("KD (GRAVY, raw): %.2f   KD mean (0–1): %.2f   μH max: %.2f (window %d)",
                  met$gravy, kd_rescale(met$gravy), mu$muH, mu$win)
  txt3 <- sprintf("Helix: %.1f%%   Sheet: %.1f%%   Coil: %.1f%%   |   HRC: %.2f%%   SAI: %.3f   Helix_Propensity: %.2f",
                  pct["Helix"], pct["Sheet"], pct["Coil"], met$HRC, SAI(seq), HelixProp(seq))
  txt4 <- sprintf("Boman: %.2f   Instability: %.2f   Charge@pH7: %.2f   pI: %.2f   Hmoment: %.3f   WWIHP: %.3f   |   Verdict: %s (%d/9)",
                  met$boman, met$instab, met$net, met$pI, met$hmoment, met$WW,
                  met$verdict, met$passed)

  ggplot()+
    annotate("text", x=0, y=1.00, hjust=0, label=txt1, size=3.2, fontface="bold", colour="#1F77B4")+
    annotate("text", x=0, y=0.75, hjust=0, label=txt2, size=3.2, fontface="bold", colour="#1F77B4")+
    annotate("text", x=0, y=0.50, hjust=0, label=txt3, size=3.2, fontface="bold", colour="#1F77B4")+
    annotate("text", x=0, y=0.25, hjust=0, label=txt4, size=3.2, fontface="bold", colour="#1F77B4")+
    xlim(0,1)+ylim(0,1)+theme_void()
}

# ---------- FASTA read ----------
fa <- read.fasta(FASTA_FILE, seqtype="AA", as.string=TRUE, forceDNAtolower=FALSE)
ids <- names(fa); seqs <- vapply(fa, as.character, character(1))
if(length(seqs)==0) stop("No sequences found in FASTA.")

# ---------- Evaluate & save summary ----------
res_list <- lapply(seqs, score_rules)
summary_df <- data.frame(
  id = ids,
  sequence = seqs,
  length = vapply(res_list, `[[`, numeric(1), "len"),
  gravy  = vapply(res_list, `[[`, numeric(1), "gravy"),
  muH    = vapply(res_list, function(x) x$mu_best$muH, numeric(1)),
  muH_win= vapply(res_list, function(x) x$mu_best$win, numeric(1)),
  HRC    = vapply(res_list, `[[`, numeric(1), "HRC"),
  WWIHP  = vapply(res_list, `[[`, numeric(1), "WW"),
  charge = vapply(res_list, `[[`, numeric(1), "net"),
  hmoment= vapply(res_list, `[[`, numeric(1), "hmoment"),
  boman  = vapply(res_list, `[[`, numeric(1), "boman"),
  instability = vapply(res_list, `[[`, numeric(1), "instab"),
  pI     = vapply(res_list, `[[`, numeric(1), "pI"),
  helix_pct = vapply(res_list, function(x) x$pct["Helix"], numeric(1)),
  sheet_pct = vapply(res_list, function(x) x$pct["Sheet"], numeric(1)),
  coil_pct  = vapply(res_list, function(x) x$pct["Coil"],  numeric(1)),
  passed = vapply(res_list, `[[`, numeric(1), "passed"),
  verdict= vapply(res_list, `[[`, character(1), "verdict"),
  stringsAsFactors = FALSE
)
write.csv(summary_df, "AMP_summary.csv", row.names=FALSE)
message("Wrote: AMP_summary.csv")

# ---------- Make full figures for all AMP-like ----------
amp_idx <- which(summary_df$verdict=="AMP-like")
if(length(amp_idx)==0){
  message("No AMP-like (≥6/9) sequences. Nothing to plot.")
} else {
  for(i in seq_along(amp_idx)){
    idx <- amp_idx[i]; id <- summary_df$id[idx]; seq <- summary_df$sequence[idx]
    met <- res_list[[idx]]
    best <- met$mu_best
    n <- nchar(seq); h <- (best$win-1)/2
    st <- max(1, best$pos - h); en <- min(n, best$pos + h)

    p_top <- kd_panel(seq, window=9)
    p_left  <- helixvis_wheel(seq, st, en, title=sprintf("helixvis wheel (%d–%d)", st, en))
    p_right <- bubble_wheel(seq, st, en,
                            title=sprintf("Eisenberg bubble wheel (%d–%d) | μH = %.3f",
                                          st,en,best$muH))
    wheels <- (p_left | p_right) + plot_layout(widths=c(1,1))
    foot   <- footer_panel(seq, met)

    p_all <- (p_top / wheels / foot) + plot_layout(heights=c(1.0,1.8,0.55))
    outfile <- sprintf("AMP_full_%02d.png", i)
    ggsave(outfile, p_all, width=12, height=8, dpi=300)
    message(sprintf("Saved: %s  [%s]", outfile, id))
  }
}

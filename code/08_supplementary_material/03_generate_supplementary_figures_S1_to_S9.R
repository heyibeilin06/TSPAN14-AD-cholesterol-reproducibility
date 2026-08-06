#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(scales)
})

root <- normalizePath(ifelse(length(commandArgs(trailingOnly=TRUE)), commandArgs(trailingOnly=TRUE)[1], "."), winslash="/")
out <- file.path(root, "outputs/supplement_v19/figures")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
read_tsv2 <- function(x) readr::read_tsv(file.path(root, x), show_col_types=FALSE, progress=FALSE)
blue <- "#1769AA"; teal <- "#128C7E"; orange <- "#D98200"; purple <- "#7A5195"; red <- "#C44E52"; grey <- "#7A8791"; light <- "#EAF0F3"
theme_pub <- theme_classic(base_family="Arial", base_size=8.2) + theme(
  plot.title=element_text(size=9.2, face="bold", color="#172B3A", margin=margin(b=5,l=14)),
  axis.title=element_text(size=8.2), axis.text=element_text(size=7.2, color="#263640"),
  legend.title=element_text(size=7.5, face="bold"), legend.text=element_text(size=7),
  legend.position="bottom", plot.margin=margin(6,8,6,8), strip.background=element_blank(),
  strip.text=element_text(face="bold", size=8), panel.grid.major.y=element_line(color="#E7ECEF", linewidth=.25)
)
tag_theme <- theme(plot.tag=element_text(family="Arial", face="bold", size=12), plot.tag.position=c(.01,.99))
save_fig <- function(p, n, w=190, h=150) {
  stem <- file.path(out, paste0("Supplementary_Figure_S", n))
  ggsave(paste0(stem,".png"), p, width=w, height=h, units="mm", dpi=300, bg="white")
  ggsave(paste0(stem,".pdf"), p, width=w, height=h, units="mm", device=cairo_pdf, bg="white")
  ggsave(paste0(stem,".svg"), p, width=w, height=h, units="mm", bg="white")
  ggsave(paste0(stem,".tiff"), p, width=w, height=h, units="mm", dpi=600, compression="lzw", bg="white")
}

# S1: APOE sensitivity
ldall <- read_tsv2("outputs/main_figures_v9/source_data/Figure_1_ldsc_apoe_conditioning.tsv")
ld <- ldall %>% filter(model=="baseline") %>% distinct(trait_label, rg, se, p, lo, hi)
p1a <- ggplot(ld, aes(rg, reorder(trait_label, rg))) + geom_vline(xintercept=0, linetype=2, color="#9AA5AC") +
  geom_errorbarh(aes(xmin=lo,xmax=hi), height=.16, color=grey) + geom_point(aes(color=p<.05), size=2.5) +
  scale_color_manual(values=c(`TRUE`=orange,`FALSE`=grey), guide="none") + labs(x="Genetic correlation (95% CI)",y=NULL,title="Baseline AD-lipid genetic correlation") + theme_pub
ap <- ldall %>% filter(analysis_group=="Extended-APOE sensitivity" | (model=="baseline" & trait=="HDL")) %>% transmute(model,rg,se,p,analysis=ifelse(grepl("condition",analysis_class,ignore.case=TRUE),"LD-conditioned","Physical-window sensitivity"),lo,hi,label=gsub("_"," ",model))
p1b <- ggplot(ap, aes(rg, reorder(label, rg), color=analysis)) + geom_vline(xintercept=0,linetype=2,color="#9AA5AC") + geom_errorbarh(aes(xmin=lo,xmax=hi),height=.14) + geom_point(size=2.2) +
  scale_color_manual(values=c("LD-conditioned"=blue,"Physical-window sensitivity"=teal)) + labs(x="AD-HDL-C genetic correlation (95% CI)",y=NULL,title="Extended-APOE sensitivity",color=NULL) + theme_pub
save_fig((p1a|p1b)+plot_annotation(tag_levels="A")&tag_theme,1,190,105)

# S2: complete regional screen
reg <- read_tsv2("outputs/main_figures_v9/source_data/Figure_1_regional_screen.tsv") %>% mutate(logit=log10((PP.H4+1e-4)/(1-PP.H4+1e-4)), focus=ifelse(grepl("TSPAN14",locus),"TSPAN14","Other loci"))
p2a <- ggplot(reg,aes(trait_label,reorder(locus,midpoint),size=PP.H4,color=PP.H4)) + geom_point(alpha=.95) + scale_size(range=c(1.5,7),limits=c(0,1)) + scale_color_gradient(low="#DDE6EA",high=blue,limits=c(0,1)) + labs(x=NULL,y=NULL,size="PP.H4",color="PP.H4",title="All screened locus-trait pairs") + theme_pub + theme(axis.text.x=element_text(angle=35,hjust=1))
p2b <- reg %>% arrange(PP.H4) %>% mutate(pair=paste(locus,trait_label,sep=" | ")) %>% ggplot(aes(PP.H4,reorder(pair,PP.H4),color=focus)) + geom_segment(aes(x=0,xend=PP.H4,yend=reorder(pair,PP.H4)),color="#C9D2D8") + geom_point(size=2.5) + geom_vline(xintercept=.8,linetype=2,color=orange) + scale_color_manual(values=c("TSPAN14"=red,"Other loci"=grey)) + labs(x="Posterior probability of a shared signal",y=NULL,title="Regional evidence ranking",color=NULL) + theme_pub
save_fig((p2a|p2b)+plot_annotation(tag_levels="A")&tag_theme,2,190,135)

# S3: colocalization and fine mapping
vf <- read_tsv2("outputs/main_figures_v9/source_data/Figure_1_variant_fingerprint.tsv") %>% filter(rank_by_product_pip<=15)
p3a <- ggplot(vf,aes(trait_label,reorder(snp,mean_pip),size=pip,color=pip)) + geom_point() + scale_size(range=c(1,6)) + scale_color_gradient(low="#DCE7EC",high=purple) + labs(x=NULL,y=NULL,size="PIP",color="PIP",title="Trait-specific posterior inclusion probabilities") + theme_pub + theme(axis.text.x=element_text(angle=30,hjust=1))
cs <- read_tsv2("outputs/supplementary_material/supplementary_table_sources/Table_S5_SuSiE_finemapping_signals.tsv")
p3b <- ggplot(cs,aes(n_variants,reorder(paste(trait,credible_set,sep=" | "),n_variants),color=top_pip)) + geom_segment(aes(x=0,xend=n_variants,yend=reorder(paste(trait,credible_set,sep=" | "),n_variants)),color="#CED6DB") + geom_point(size=3) + scale_color_gradient(low=teal,high=red) + labs(x="Variants in credible set",y=NULL,color="Top PIP",title="Multiple-signal credible sets") + theme_pub
ec <- read_tsv2("outputs/main_figures_v9/source_data/Figure_2_exact_event_coloc.tsv")
p3c <- ggplot(ec,aes(pph4,reorder(trait,pph4))) + geom_col(width=.62,fill=blue) + geom_vline(xintercept=.8,linetype=2,color=orange) + coord_cartesian(xlim=c(.8,1)) + labs(x="Default-prior exact-event PP.H4",y=NULL,title="Exact exon5-6 colocalization") + theme_pub
save_fig(((p3a|p3b)/p3c)+plot_annotation(tag_levels="A")&tag_theme,3,190,170)

# S4: raw IER magnitude and robustness
raw <- read_tsv2("outputs/mentor_revision/leafcutter_delta_ier/02_ba24_sample_raw_ier_and_genotype.tsv")
gcol <- intersect(names(raw),c("genotype","alt_allele_dosage","risk_allele_dosage"))[1]; icol <- intersect(names(raw),c("ier","raw_ier","target_ier"))[1]
raw$geno <- factor(raw[[gcol]]); raw$IER <- raw[[icol]]*ifelse(max(raw[[icol]],na.rm=TRUE)<=1.2,100,1)
p4a <- ggplot(raw,aes(geno,IER,fill=geno)) + geom_violin(width=.75,alpha=.28,color=NA,trim=TRUE) + geom_boxplot(width=.28,outlier.shape=NA,alpha=.72) + geom_jitter(width=.08,size=.8,alpha=.55,color="#243746") + scale_fill_manual(values=c("0"="#BCD8E8","1"=blue,"2"=purple),guide="none") + labs(x="Risk-allele dosage",y="Exon5-6 IER (%)",title="Donor-level effect magnitude") + theme_pub
dep <- read_tsv2("outputs/mentor_revision/leafcutter_delta_ier/04_depth_sensitivity_delta_ier.tsv")
p4b <- ggplot(dep,aes(minimum_cluster_reads,per_alt_allele_delta_ier_percentage_points)) + geom_ribbon(aes(ymin=bootstrap_95ci_low,ymax=bootstrap_95ci_high),fill="#B9D9D3",alpha=.6) + geom_line(color=teal,linewidth=.8)+geom_point(color=teal,size=2) + labs(x="Minimum cluster reads",y="Per-allele change (percentage points)",title="Read-depth sensitivity") + theme_pub
cc <- read_tsv2("outputs/mentor_revision/leafcutter_delta_ier/06_recount3_official_target_count_concordance.tsv")
num <- names(cc)[vapply(cc,is.numeric,logical(1))]; xcol <- num[1]; ycol <- num[2]
p4c <- ggplot(cc,aes(.data[[xcol]],.data[[ycol]])) + geom_point(size=.9,alpha=.5,color=blue) + geom_smooth(method="lm",se=FALSE,color=red,linewidth=.7) + labs(x=gsub("_"," ",xcol),y=gsub("_"," ",ycol),title="Independent count-source concordance") + theme_pub
save_fig(((p4a|p4b)/p4c)+plot_annotation(tag_levels="A")&tag_theme,4,190,165)

# S5: exact replication
rep <- read_tsv2("outputs/main_figures_v9/source_data/Figure_3_exact_replication_matrix.tsv") %>% mutate(signed=-log10(p_value)*sign(nes))
p5a <- ggplot(rep,aes(tissue_label,reorder(snp_label,position),fill=signed)) + geom_tile(color="white",linewidth=.4) + scale_fill_gradient2(low=purple,mid="white",high=orange,midpoint=0,name="Signed\n-log10(P)") + labs(x=NULL,y=NULL,title="Coordinate-identical exon5-6 sQTL") + theme_pub + theme(axis.text.x=element_text(angle=30,hjust=1))
p5b <- rep %>% group_by(tissue_label) %>% summarise(median_nes=median(nes),min_p=min(p_value),n_variants=n(),.groups="drop") %>% ggplot(aes(median_nes,reorder(tissue_label,median_nes))) + geom_vline(xintercept=0,linetype=2,color="#A8B0B5") + geom_point(aes(size=-log10(min_p)),color=blue) + labs(x="Median risk-aligned NES",y=NULL,size="Strongest\n-log10(P)",title="Direction across neural tissues") + theme_pub
save_fig((p5a|p5b)+plot_annotation(tag_levels="A")&tag_theme,5,190,110)

# S6: adjacent-junction co-usage
co <- read_tsv2("outputs/main_figures_v9/source_data/Figure_3_brain_cousage_summary.tsv") %>% mutate(tissue_label=ifelse(is.na(tissue_label),tissue,tissue_label))
p6a <- ggplot(co,aes(spearman_rho,reorder(tissue_label,spearman_rho),size=n_both_nonzero,color=spearman_rho)) + geom_vline(xintercept=0,linetype=2,color="#A8B0B5") + geom_point(alpha=.9) + scale_color_gradient(low="#BCD8E8",high=teal) + coord_cartesian(xlim=c(0,1)) + labs(x="Spearman correlation",y=NULL,size="Samples with\nboth junctions",color="rho",title="Exon5-6 and exon6-7 co-usage") + theme_pub
p6b <- co %>% select(tissue_label,n_samples,n_donors,n_both_nonzero) %>% pivot_longer(-tissue_label,names_to="metric",values_to="count") %>% ggplot(aes(count,reorder(tissue_label,count),color=metric)) + geom_point(size=2) + facet_wrap(~metric,scales="free_x",nrow=1) + scale_color_manual(values=c(n_samples=blue,n_donors=orange,n_both_nonzero=teal),guide="none") + labs(x="Count",y=NULL,title="Coverage supporting the correlation estimate") + theme_pub
save_fig((p6a/p6b)+plot_annotation(tag_levels="A")&tag_theme,6,190,145)

# S7: cis-MR diagnostics
mr <- read_tsv2("outputs/main_figures_v9/source_data/Figure_4_ld_aware_cis_mr.tsv")
p7a <- ggplot(mr,aes(estimate,reorder(paste(outcome_label,method_label,sep=" | "),estimate),color=method_label)) + geom_vline(xintercept=0,linetype=2,color="#A8B0B5") + geom_errorbarh(aes(xmin=lo,xmax=hi),height=.15) + geom_point(size=2.2) + scale_color_manual(values=c("LD-aware generalized IVW"=blue,"Lead-instrument Wald"=orange)) + labs(x="Risk-aligned estimate (95% CI)",y=NULL,color=NULL,title="Exact-event cis-MR") + theme_pub
dg <- read_tsv2("outputs/main_figures_v9/source_data/Figure_4_cis_mr_diagnostics.tsv")
p7b <- ggplot(dg,aes(diagnostic_label,outcome_label,fill=pass)) + geom_tile(color="white") + geom_text(aes(label=ifelse(pass,"PASS","FLAG")),size=2.4) + scale_fill_manual(values=c(`TRUE`="#B8DDD4",`FALSE`="#F4C7C3"),guide="none") + labs(x=NULL,y=NULL,title="Instrument and model diagnostics") + theme_pub + theme(axis.text.x=element_text(angle=30,hjust=1))
ov <- read_tsv2("outputs/mentor_revision/cis_mr_sample_overlap/04_overlap_sensitivity_summary.tsv")
p7c <- ggplot(ov,aes(baseline_estimate,reorder(outcome,baseline_estimate))) + geom_errorbarh(aes(xmin=overlap_sensitivity_min_estimate,xmax=overlap_sensitivity_max_estimate),height=.18,color=teal,linewidth=.8) + geom_point(size=2.3,color=blue) + labs(x="Estimate under overlap scenarios",y=NULL,title="Maximal participant-overlap sensitivity") + theme_pub
save_fig(((p7a|p7b)/p7c)+plot_annotation(tag_levels="A")&tag_theme,7,190,165)

# S8: causal scope
gm <- read_tsv2("outputs/mentor_revision/complete_mr/03_genomewide_bidirectional_mr.tsv") %>% mutate(lo=estimate-1.96*se,hi=estimate+1.96*se,pair=paste(exposure,"to",outcome))
p8a <- ggplot(gm,aes(estimate,reorder(pair,estimate))) + geom_vline(xintercept=0,linetype=2,color="#A8B0B5") + geom_errorbarh(aes(xmin=lo,xmax=hi),height=.15,color=grey) + geom_point(aes(color=pvalue<.05),size=2.2) + scale_color_manual(values=c(`TRUE`=red,`FALSE`=grey),guide="none") + labs(x="Genome-wide MR estimate (95% CI)",y=NULL,title="Bidirectional systemic effects") + theme_pub
pc <- read_tsv2("outputs/mentor_revision/mediation_rescue/14_pc_dimension_mediation_sensitivity.tsv")
p8b <- ggplot(pc,aes(n_pcs,indirect_estimate,color=lipid)) + geom_hline(yintercept=0,linetype=2,color="#A8B0B5") + geom_line(linewidth=.65,na.rm=TRUE)+geom_point(aes(shape=all_strength_F_ge_10),size=1.8,na.rm=TRUE) + scale_color_manual(values=c(TC=orange,LDL=blue,nonHDL=teal),labels=c(TC="TC",LDL="LDL-C",nonHDL="non-HDL-C")) + labs(x="Retained LD principal components",y="Estimated indirect effect",color=NULL,shape="All F >= 10",title="PC-GMM dimension sensitivity") + guides(shape=guide_legend(order=1,nrow=1),color=guide_legend(order=2,nrow=1)) + theme_pub + theme(legend.box="vertical",legend.key.width=grid::unit(4,"mm"),legend.margin=margin(0,0,0,0))
p8c <- ggplot(pc,aes(n_pcs,pmin(first_step_instrument_F,splice_conditional_F,lipid_conditional_F,na.rm=TRUE),color=lipid)) + geom_hline(yintercept=10,linetype=2,color=red) + geom_line(linewidth=.65,na.rm=TRUE)+geom_point(size=1.6,na.rm=TRUE) + scale_color_manual(values=c(TC=orange,LDL=blue,nonHDL=teal),guide="none") + labs(x="Retained LD principal components",y="Minimum instrument-strength statistic",title="Identification strength") + theme_pub
save_fig(((p8a|p8b)/p8c)+plot_annotation(tag_levels="A")&tag_theme,8,190,165)

# S9: cell context and transcript-to-structure support
atlas <- read_tsv2("outputs/main_figures_v9/source_data/Figure_5_cell_context_atlas.tsv") %>% mutate(score=as.numeric(evidence_strength))
p9a <- ggplot(atlas,aes(evidence_layer,reorder(context,score),size=score,color=evidence_class)) + geom_point(alpha=.9) + scale_size(range=c(2,6),breaks=2:4,labels=c("Contextual","Moderate","Direct")) + labs(x=NULL,y=NULL,size="Evidence",color=NULL,title="Neural cell-context evidence") + theme_pub + theme(axis.text.x=element_text(angle=28,hjust=1))
dis <- read_tsv2("outputs/main_figures_v9/source_data/Figure_5_disease_state_rna.tsv") %>% mutate(label=paste(source_label,cell_label,sep=" | "))
p9b <- ggplot(dis,aes(estimate,reorder(label,estimate),color=estimate>0)) + geom_vline(xintercept=0,linetype=2,color="#A8B0B5") + geom_errorbarh(aes(xmin=lo,xmax=hi),height=.15,na.rm=TRUE) + geom_point(size=2.2) + scale_color_manual(values=c(`TRUE`=red,`FALSE`=blue),guide="none") + labs(x="Disease-state estimate (95% CI)",y=NULL,title="Disease-state RNA sensitivity") + theme_pub
st <- read_tsv2("outputs/main_figures_v9/source_data/Figure_5_ec2_structure.tsv")
p9c <- ggplot(st,aes(residue,pLDDT)) + annotate("rect",xmin=114,xmax=232,ymin=-Inf,ymax=Inf,fill="#DDF0EC",alpha=.6) + geom_line(linewidth=.65,color=grey) + geom_vline(xintercept=150.5,color=purple,linewidth=.8) + annotate("text",x=153,y=max(st$pLDDT,na.rm=TRUE),label="AA150/151",hjust=0,size=2.7,color=purple) + labs(x="TSPAN14 residue",y="AlphaFold pLDDT",title="Exact splice boundary within the EC2 region") + theme_pub
save_fig(((p9a|p9b)/p9c)+plot_annotation(tag_levels="A")&tag_theme,9,190,170)

cat("Generated 9 supplementary figures in", out, "\n")

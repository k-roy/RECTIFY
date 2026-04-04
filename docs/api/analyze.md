# rectify.core.analyze

Downstream analysis modules: clustering, differential expression, APA detection, GO enrichment, and motif discovery.

---

## clustering

Valley-based CPA site clustering and count matrix building.

::: rectify.core.analyze.clustering
    options:
      members:
        - cluster_cpa_sites_adaptive
        - cluster_cpa_sites
        - build_cluster_count_matrix

---

## deseq2

DESeq2 differential expression at gene and cluster resolution.

::: rectify.core.analyze.deseq2
    options:
      members:
        - run_deseq2_gene_level
        - run_deseq2_cluster_level
        - detect_control_samples

---

## apa_detection

Alternative polyadenylation isoform detection and usage quantification.

::: rectify.core.analyze.apa_detection
    options:
      members:
        - APAIsoform
        - GeneAPAProfile
        - detect_apa_isoforms
        - build_gene_apa_profiles
        - quantify_apa_usage
        - identify_proximal_distal_tes

---

## shift_analysis

APA shift analysis between conditions using Jensen-Shannon divergence.

::: rectify.core.analyze.shift_analysis
    options:
      members:
        - analyze_cluster_shifts
        - get_top_shifted_genes
        - compute_js_divergence

---

## gene_attribution

Read-to-gene assignment based on 5' and 3' end positions.

::: rectify.core.analyze.gene_attribution
    options:
      members:
        - compute_read_gene_attribution
        - aggregate_attributions_for_3prime_end

---

## go_enrichment

GO term enrichment using hypergeometric test with BH FDR correction.

::: rectify.core.analyze.go_enrichment
    options:
      members:
        - run_go_enrichment
        - plot_go_enrichment

---

## motif_discovery

De novo motif discovery using STREME (MEME Suite).

::: rectify.core.analyze.motif_discovery
    options:
      members:
        - run_motif_discovery
        - extract_sequences_around_clusters

---

## deconvolution

NET-seq point-spread function fitting and NNLS deconvolution.

::: rectify.core.analyze.deconvolution
    options:
      members:
        - build_convolution_matrix
        - deconvolve_signal
        - fit_psf_from_controls

---

## genomic_distribution

Classify CPA positions into genomic regions (3' UTR, CDS, intron, intergenic).

::: rectify.core.analyze.genomic_distribution
    options:
      members:
        - classify_positions_by_region
        - run_genomic_distribution_analysis

---

## pca

Principal component analysis of samples by CPA cluster usage.

::: rectify.core.analyze.pca
    options:
      members:
        - run_pca_analysis
        - plot_pca

---

## heatmap

Heatmap visualization of clusters and samples.

::: rectify.core.analyze.heatmap
    options:
      members:
        - plot_sample_heatmap
        - plot_cluster_heatmap

---

## summary

HTML report generation.

::: rectify.core.analyze.summary
    options:
      members:
        - generate_analysis_summary
        - generate_html_report

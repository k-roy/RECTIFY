#!/usr/bin/env python3
"""
Summary Report Generation Module

Generates summary tables and HTML reports for RECTIFY analysis.

Author: Kevin R. Roy
Date: 2026-03-17
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd


def generate_summary_tables(
    deseq2_gene_results: Dict[str, pd.DataFrame],
    deseq2_cluster_results: Dict[str, pd.DataFrame],
    shift_results: Optional[pd.DataFrame] = None,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
    output_dir: str = '.',
) -> Dict[str, str]:
    """
    Generate summary tables from analysis results.

    Creates the following tables:
    1. Significantly upregulated genes
    2. Significantly downregulated genes
    3. Significantly changed clusters
    4. Top shifted genes (if shift analysis provided)

    Args:
        deseq2_gene_results: Gene-level DESeq2 results
        deseq2_cluster_results: Cluster-level DESeq2 results
        shift_results: Optional shift analysis results
        padj_threshold: Adjusted p-value threshold
        lfc_threshold: Log2 fold-change threshold
        output_dir: Output directory for tables

    Returns:
        Dict mapping table_name -> output_path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files = {}

    # Gene-level tables
    for condition, results_df in deseq2_gene_results.items():
        # Significant genes
        sig = results_df['padj'] < padj_threshold

        # Upregulated
        up = sig & (results_df['log2FoldChange'] > lfc_threshold)
        up_df = results_df[up].sort_values('log2FoldChange', ascending=False)

        up_path = output_dir / f'genes_upregulated_{condition}.tsv'
        up_df.to_csv(up_path, sep='\t')
        output_files[f'genes_up_{condition}'] = str(up_path)

        # Downregulated
        down = sig & (results_df['log2FoldChange'] < -lfc_threshold)
        down_df = results_df[down].sort_values('log2FoldChange')

        down_path = output_dir / f'genes_downregulated_{condition}.tsv'
        down_df.to_csv(down_path, sep='\t')
        output_files[f'genes_down_{condition}'] = str(down_path)

    # Cluster-level tables
    for condition, results_df in deseq2_cluster_results.items():
        sig = results_df['padj'] < padj_threshold

        # Upregulated clusters
        up = sig & (results_df['log2FoldChange'] > lfc_threshold)
        up_df = results_df[up].sort_values('log2FoldChange', ascending=False)

        up_path = output_dir / f'clusters_upregulated_{condition}.tsv'
        up_df.to_csv(up_path, sep='\t')
        output_files[f'clusters_up_{condition}'] = str(up_path)

        # Downregulated clusters
        down = sig & (results_df['log2FoldChange'] < -lfc_threshold)
        down_df = results_df[down].sort_values('log2FoldChange')

        down_path = output_dir / f'clusters_downregulated_{condition}.tsv'
        down_df.to_csv(down_path, sep='\t')
        output_files[f'clusters_down_{condition}'] = str(down_path)

    # Shift analysis
    if shift_results is not None and not shift_results.empty:
        # Top shifted genes
        top_shifted = shift_results.sort_values(
            'distribution_divergence', ascending=False
        ).head(100)

        shift_path = output_dir / 'top_shifted_genes.tsv'
        top_shifted.to_csv(shift_path, sep='\t', index=False)
        output_files['top_shifted'] = str(shift_path)

    return output_files


def generate_analysis_summary(
    n_samples: int,
    n_clusters: int,
    n_genes: int,
    deseq2_gene_results: Dict[str, pd.DataFrame],
    deseq2_cluster_results: Dict[str, pd.DataFrame],
    reference_condition: str,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
) -> pd.DataFrame:
    """
    Generate overall analysis summary statistics.

    Args:
        n_samples: Number of samples
        n_clusters: Number of CPA clusters
        n_genes: Number of genes with clusters
        deseq2_gene_results: Gene-level DESeq2 results
        deseq2_cluster_results: Cluster-level DESeq2 results
        reference_condition: Reference condition name
        padj_threshold: Significance threshold
        lfc_threshold: Fold-change threshold

    Returns:
        Summary DataFrame
    """
    rows = []

    # Overall stats
    rows.append({'category': 'Input', 'metric': 'Samples', 'value': n_samples})
    rows.append({'category': 'Input', 'metric': 'CPA Clusters', 'value': n_clusters})
    rows.append({'category': 'Input', 'metric': 'Genes with Clusters', 'value': n_genes})
    rows.append({'category': 'Input', 'metric': 'Reference Condition', 'value': reference_condition})

    # Gene-level results
    for condition, df in deseq2_gene_results.items():
        sig = df['padj'] < padj_threshold
        up = sig & (df['log2FoldChange'] > lfc_threshold)
        down = sig & (df['log2FoldChange'] < -lfc_threshold)

        rows.append({'category': f'Gene-level ({condition})', 'metric': 'Tested', 'value': len(df)})
        rows.append({'category': f'Gene-level ({condition})', 'metric': 'Significant (padj<0.05)', 'value': sig.sum()})
        rows.append({'category': f'Gene-level ({condition})', 'metric': 'Upregulated (LFC>1)', 'value': up.sum()})
        rows.append({'category': f'Gene-level ({condition})', 'metric': 'Downregulated (LFC<-1)', 'value': down.sum()})

    # Cluster-level results
    for condition, df in deseq2_cluster_results.items():
        sig = df['padj'] < padj_threshold
        up = sig & (df['log2FoldChange'] > lfc_threshold)
        down = sig & (df['log2FoldChange'] < -lfc_threshold)

        rows.append({'category': f'Cluster-level ({condition})', 'metric': 'Tested', 'value': len(df)})
        rows.append({'category': f'Cluster-level ({condition})', 'metric': 'Significant (padj<0.05)', 'value': sig.sum()})
        rows.append({'category': f'Cluster-level ({condition})', 'metric': 'Upregulated (LFC>1)', 'value': up.sum()})
        rows.append({'category': f'Cluster-level ({condition})', 'metric': 'Downregulated (LFC<-1)', 'value': down.sum()})

    return pd.DataFrame(rows)


def generate_html_report(
    summary_df: pd.DataFrame,
    plots: Dict[str, str],
    tables: Dict[str, str],
    output_path: str,
    title: str = 'RECTIFY Analysis Report',
) -> str:
    """
    Generate HTML report with embedded plots and table links.

    Args:
        summary_df: Summary statistics DataFrame
        plots: Dict mapping plot_name -> plot_path
        tables: Dict mapping table_name -> table_path
        output_path: Output HTML file path
        title: Report title

    Returns:
        Path to generated HTML file
    """
    import base64

    html_parts = [
        '<!DOCTYPE html>',
        '<html lang="en">',
        '<head>',
        '<meta charset="UTF-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
        f'<title>{title}</title>',
        '<style>',
        '''
        body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }
        .container { max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .plot-container { margin: 20px 0; text-align: center; }
        .plot-container img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; }
        .links { background-color: #ecf0f1; padding: 15px; border-radius: 4px; margin: 20px 0; }
        .links a { display: inline-block; margin: 5px 10px; color: #2980b9; }
        .timestamp { color: #7f8c8d; font-size: 12px; text-align: right; }
        ''',
        '</style>',
        '</head>',
        '<body>',
        '<div class="container">',
        f'<h1>{title}</h1>',
        f'<p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>',
    ]

    # Summary table
    html_parts.append('<h2>Summary Statistics</h2>')
    html_parts.append(summary_df.to_html(index=False, classes='summary-table'))

    # Plots
    if plots:
        html_parts.append('<h2>Visualizations</h2>')
        for name, path in plots.items():
            if Path(path).exists():
                # Embed plot as base64
                with open(path, 'rb') as f:
                    img_data = base64.b64encode(f.read()).decode()
                ext = Path(path).suffix.lower()
                mime = 'image/png' if ext == '.png' else 'image/svg+xml' if ext == '.svg' else 'image/jpeg'

                html_parts.append(f'<div class="plot-container">')
                html_parts.append(f'<h3>{name.replace("_", " ").title()}</h3>')
                html_parts.append(f'<img src="data:{mime};base64,{img_data}" alt="{name}">')
                html_parts.append('</div>')

    # Table links
    if tables:
        html_parts.append('<h2>Output Tables</h2>')
        html_parts.append('<div class="links">')
        for name, path in tables.items():
            rel_path = Path(path).name
            html_parts.append(f'<a href="{rel_path}">{name.replace("_", " ").title()}</a>')
        html_parts.append('</div>')

    # Close HTML
    html_parts.extend([
        '</div>',
        '</body>',
        '</html>',
    ])

    # Write HTML
    with open(output_path, 'w') as f:
        f.write('\n'.join(html_parts))

    return output_path


def generate_consensus_html_report(
    stats: dict,
    output_path: str,
    sample_name: str = '',
) -> str:
    """Generate a standalone HTML report for consensus aligner selection stats.

    Args:
        stats:       Stats dict returned by run_consensus_selection().
        output_path: Path to write the HTML file.
        sample_name: Optional sample / prefix label shown in the title.

    Returns:
        Path to the generated HTML file.
    """
    total     = stats.get('total_reads', 0)
    by_aligner = dict(stats.get('by_aligner', {}))
    high      = stats.get('consensus_high', 0)
    medium    = stats.get('consensus_medium', 0)
    low       = stats.get('consensus_low', 0)
    rescued   = stats.get('5prime_rescued', 0)
    tied      = stats.get('tied_score', 0)
    chimeric  = stats.get('chimeric_reads', 0)

    def _pct(n: int) -> str:
        return f"{100.0 * n / total:.1f}%" if total else "—"

    title = f"RECTIFY Consensus Report"
    if sample_name:
        title += f" — {sample_name}"

    # Build aligner rows
    aligner_rows = ""
    for aligner, count in sorted(by_aligner.items(), key=lambda x: -x[1]):
        aligner_rows += (
            f"<tr><td>{aligner}</td><td>{count:,}</td>"
            f"<td>{_pct(count)}</td></tr>\n"
        )

    # Build aligner-combo rows. Keys are frozensets of the aligner panel that
    # was available/compared for each read (diagnoses which combos dominate and
    # where panel coverage is uneven). May be absent on older stats dicts.
    by_aligner_combo = stats.get('by_aligner_combo', {})
    combo_rows = ""
    for combo, count in sorted(
        by_aligner_combo.items(), key=lambda x: -x[1]
    ):
        label = '+'.join(sorted(combo)) if combo else '(none)'
        combo_rows += (
            f"<tr><td>{label}</td><td>{count:,}</td>"
            f"<td>{_pct(count)}</td></tr>\n"
        )

    # Build summary rows
    summary_rows = (
        f"<tr><td>Total reads</td><td>{total:,}</td><td>100%</td></tr>\n"
        f"<tr><td>High confidence</td><td>{high:,}</td><td>{_pct(high)}</td></tr>\n"
        f"<tr><td>Medium confidence</td><td>{medium:,}</td><td>{_pct(medium)}</td></tr>\n"
        f"<tr><td>Low confidence</td><td>{low:,}</td><td>{_pct(low)}</td></tr>\n"
        f"<tr><td>5′ splice rescued</td><td>{rescued:,}</td><td>{_pct(rescued)}</td></tr>\n"
        f"<tr><td>Tiebreaker used</td><td>{tied:,}</td><td>{_pct(tied)}</td></tr>\n"
    )
    if chimeric:
        summary_rows += (
            f"<tr><td>Chimeric reads</td><td>{chimeric:,}</td>"
            f"<td>{_pct(chimeric)}</td></tr>\n"
        )

    combo_section = ""
    if combo_rows:
        combo_section = f"""
<h2>Aligner Combination Breakdown</h2>
<p style="color:#7f8c8d;font-size:13px;">Panel of aligners available/compared per read.</p>
<table>
  <tr><th>Aligner Combination</th><th>Reads</th><th>Percent of Total</th></tr>
  {combo_rows}
</table>
"""

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }}
  .container {{ max-width: 900px; margin: 0 auto; background: white; padding: 20px;
               border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
  h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
  h2 {{ color: #34495e; margin-top: 30px; }}
  table {{ border-collapse: collapse; width: 100%; margin: 16px 0; }}
  th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
  th {{ background: #3498db; color: white; }}
  tr:nth-child(even) {{ background: #f2f2f2; }}
  .timestamp {{ color: #7f8c8d; font-size: 12px; text-align: right; }}
  .winner {{ font-weight: bold; color: #27ae60; }}
</style>
</head>
<body>
<div class="container">
<h1>{title}</h1>
<p class="timestamp">Generated: {ts}</p>

<h2>Aligner Win Distribution</h2>
<table>
  <tr><th>Aligner</th><th>Reads Won</th><th>Percent</th></tr>
  {aligner_rows}
</table>

<h2>Selection Summary</h2>
<table>
  <tr><th>Metric</th><th>Count</th><th>Percent of Total</th></tr>
  {summary_rows}
</table>
{combo_section}
</div>
</body>
</html>
"""

    with open(output_path, 'w') as fh:
        fh.write(html)

    return output_path


def create_excel_report(
    summary_df: pd.DataFrame,
    deseq2_gene_results: Dict[str, pd.DataFrame],
    deseq2_cluster_results: Dict[str, pd.DataFrame],
    shift_results: Optional[pd.DataFrame],
    output_path: str,
) -> str:
    """
    Create Excel workbook with multiple sheets.

    Args:
        summary_df: Summary statistics
        deseq2_gene_results: Gene-level results
        deseq2_cluster_results: Cluster-level results
        shift_results: Shift analysis results
        output_path: Output Excel file path

    Returns:
        Path to generated Excel file
    """
    try:
        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            # Summary sheet
            summary_df.to_excel(writer, sheet_name='Summary', index=False)

            # Gene-level results
            for condition, df in deseq2_gene_results.items():
                sheet_name = f'Genes_{condition[:25]}'  # Excel sheet name limit (31 chars total)
                df.to_excel(writer, sheet_name=sheet_name)

            # Cluster-level results
            for condition, df in deseq2_cluster_results.items():
                sheet_name = f'Clusters_{condition[:22]}'  # Excel sheet name limit (31 chars total)
                df.to_excel(writer, sheet_name=sheet_name)

            # Shift results
            if shift_results is not None and not shift_results.empty:
                shift_results.to_excel(writer, sheet_name='Shift_Analysis', index=False)

        return output_path

    except ImportError:
        print("Warning: openpyxl required for Excel output. Install with: pip install openpyxl")
        return ''

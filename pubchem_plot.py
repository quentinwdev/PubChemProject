import argparse
import os
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

# -------------------------
# Helpers
# -------------------------


def detect_delimiter(filepath, nbytes=4096):
    """Try to detect delimiter using csv.Sniffer. Fallback to comma."""
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        sample = f.read(nbytes)
    try:
        dialect = csv.Sniffer().sniff(sample)
        return dialect.delimiter
    except Exception:
        return ','


def clean_count_series(s):
    """
    Convert a series that may contain commas, whitespace, or non-numeric markers to numeric.
    Non-convertible entries become NaN.
    """
    s2 = s.astype(str).str.strip()
    # Remove thousands separators like "1,234"
    s2 = s2.str.replace(',', '', regex=False)
    # Convert to numeric
    return pd.to_numeric(s2, errors='coerce')


def ecdf_values(arr):
    arr_sorted = np.sort(arr)
    n = arr_sorted.size
    y = np.arange(1, n + 1) / n
    return arr_sorted, y


def small_pmf_plot(series, name, outpath):
    vc = series.value_counts(normalize=True).sort_index()
    plt.figure(figsize=(8, 4))
    vc.plot(kind='bar')
    plt.xlabel(name)
    plt.ylabel('Probability (relative frequency)')
    plt.title(f'PMF (value counts) of {name}')
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()


def save_summary_row(rows, filename):
    df = pd.DataFrame(rows)
    df.to_csv(filename, index=False)

# -------------------------
# Main plotting function
# -------------------------


def plot_distributions(series, name, outdir, summary_rows):
    s = series.dropna()
    if s.size == 0:
        print(
            f"[warning] {name} has no numeric values after cleaning. Skipping.")
        return

    # Basic stats
    stats = OrderedDict()
    stats['column'] = name
    stats['count_nonnull'] = int(s.size)
    stats['count_missing'] = int(series.size - s.size)
    stats['mean'] = float(s.mean())
    stats['median'] = float(s.median())
    stats['std'] = float(s.std())
    stats['min'] = float(s.min())
    stats['25%'] = float(s.quantile(0.25))
    stats['75%'] = float(s.quantile(0.75))
    stats['max'] = float(s.max())
    stats['unique_count'] = int(s.nunique())
    stats['zero_count'] = int((s == 0).sum())
    stats['num_positive'] = int((s > 0).sum())
    summary_rows.append(stats)

    # Histogram (normalized as probability density)
    plt.figure(figsize=(7, 4))
    plt.hist(s, bins='auto', density=True, alpha=0.75)
    plt.xlabel(name)
    plt.ylabel('Probability density')
    plt.title(f'Histogram (density) of {name}')
    plt.tight_layout()
    outfile = os.path.join(outdir, f"{name}_hist.png")
    plt.savefig(outfile, dpi=150)
    plt.close()

    # Histogram on log1p scale (handles zeros)
    plt.figure(figsize=(7, 4))
    plt.hist(np.log1p(s), bins='auto', density=True, alpha=0.75)
    plt.xlabel(f'log1p({name})')
    plt.ylabel('Probability density')
    plt.title(f'Histogram of log1p({name})')
    plt.tight_layout()
    outfile = os.path.join(outdir, f"{name}_hist_log1p.png")
    plt.savefig(outfile, dpi=150)
    plt.close()

    # ECDF
    x_ecdf, y_ecdf = ecdf_values(s.values)
    plt.figure(figsize=(7, 4))
    plt.step(x_ecdf, y_ecdf, where='post')
    plt.xlabel(name)
    plt.ylabel('ECDF (proportion)')
    plt.title(f'ECDF of {name}')
    plt.tight_layout()
    outfile = os.path.join(outdir, f"{name}_ecdf.png")
    plt.savefig(outfile, dpi=150)
    plt.close()

    # Boxplot (horizontal)
    plt.figure(figsize=(7, 2.5))
    plt.boxplot(s, vert=False, notch=True, showfliers=True)
    plt.xlabel(name)
    plt.title(f'Boxplot of {name}')
    plt.tight_layout()
    outfile = os.path.join(outdir, f"{name}_boxplot.png")
    plt.savefig(outfile, dpi=150)
    plt.close()

    # PMF (value counts) only if not too many unique values
    if stats['unique_count'] <= 120:
        outfile = os.path.join(outdir, f"{name}_pmf.png")
        small_pmf_plot(s, name, outfile)

    print(f"[info] Saved plots for {name} in {outdir}")

# -------------------------
# CLI & main
# -------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Plot distribution of PubMed_Count and Patent_Count")
    parser.add_argument("csvfile", help="Path to PubChemLite CSV file")
    parser.add_argument("--outdir", default="results/plots",
                        help="Output directory to save plots")
    parser.add_argument("--columns", nargs='*', default=["PubMed_Count", "Patent_Count"],
                        help="Columns to analyze (default: PubMed_Count Patent_Count)")
    args = parser.parse_args()

    infile = args.csvfile
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    delim = detect_delimiter(infile)
    print(f"[info] Using delimiter: '{delim}' to read {infile}")

    # Read CSV
    try:
        df = pd.read_csv(infile, sep=delim, engine='c',
                         dtype=str, encoding='utf-8', low_memory=False)
    except Exception as e:
        # Fallback to python engine if c engine fails
        try:
            df = pd.read_csv(infile, sep=delim, engine='python',
                             dtype=str, encoding='utf-8')
        except Exception as e2:
            raise SystemExit(f"Failed to read CSV: {e2}")

    # Clean and plot each target column
    summary_rows = []
    for col in args.columns:
        if col not in df.columns:
            print(
                f"[warning] Column '{col}' not found in CSV. Available columns: {list(df.columns)[:10]} ...")
            continue
        cleaned = clean_count_series(df[col])
        plot_distributions(cleaned, col, outdir, summary_rows)

    # Save summary stats
    if summary_rows:
        summary_file = os.path.join(outdir, "summary_stats.csv")
        save_summary_row(summary_rows, summary_file)
        print(f"[info] Wrote summary stats to {summary_file}")

    print("[done] All finished.")


if __name__ == "__main__":
    main()

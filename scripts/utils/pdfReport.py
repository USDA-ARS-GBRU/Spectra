#!/usr/bin/env python3
import logging
import os
import argparse
import re

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet

### Prepares a ReportLab Image object with restricted size.
def image_prep(path, width=6.5 * inch, height=9 * inch):
    try:
        img = Image(path)
        img._restrictSize(width, height)
        return img
    except Exception as e:
        logging.error(f"Failed to load image {path}: {e}")
        return None

### Adds an image to the story if it exists, otherwise adds an error message.
def add_safe_image(story, path, width, height, styles, spacer=0.1):
    if os.path.exists(path):
        if spacer:
            story.append(Spacer(1, spacer * inch))
        img = image_prep(path, width, height)
        if img:
            story.append(img)
            return True

    story.append(Paragraph(f"<b>ERROR:</b> Could not find file {os.path.basename(path)}.", styles["Normal"]))
    return False

### Adds a title page with run information.
def add_title_page(story, styles, logo_path, prefix, raw_files, assembly_file, counter, mer_size):
    if logo_path and os.path.exists(logo_path):
        logo = image_prep(logo_path, 4 * inch, 4 * inch)
        if logo:
            story.append(Spacer(1, 1 * inch))
            story.append(logo)
            story.append(Spacer(1, 0.5 * inch))

    story.append(Paragraph(f"<font size=24><b>Spectra Analysis Report</b></font>", styles["Title"]))
    story.append(Spacer(1, 0.5 * inch))
    story.append(Paragraph(f"<font size=14><b>Project Prefix:</b> {prefix}</font>", styles["Normal"]))
    story.append(Spacer(1, 0.2 * inch))

    # Parameters table
    data = [
        [Paragraph("<b>Parameter</b>", styles["Normal"]), Paragraph("<b>Value</b>", styles["Normal"])],
        ["K-mer size", str(mer_size)],
        ["Counter", counter],
        [Paragraph("Assembly file", styles["Normal"]), Paragraph(os.path.basename(assembly_file) if assembly_file else "N/A", styles["Normal"])]
    ]

    if raw_files:
        for i, raw in enumerate(raw_files):
            label = "Raw file(s)" if i == 0 else ""
            data.append([label, Paragraph(os.path.basename(raw), styles["Normal"])])

    table = Table(data, colWidths=[1.5 * inch, 4.5 * inch])
    table.setStyle([
        ('GRID', (0, 0), (-1, -1), 0.5, 'grey'),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('PADDING', (0, 0), (-1, -1), 6)
    ])
    story.append(table)
    story.append(PageBreak())

### Adds a row of images using a Table.
def add_image_row(story, paths, widths, heights, styles):
    row = []
    all_exist = True
    for path, w, h in zip(paths, widths, heights):
        if os.path.exists(path):
            img = image_prep(path, w, h)
            if img:
                row.append(img)
            else:
                all_exist = False
        else:
            all_exist = False

    if all_exist and len(row) == len(paths):
        table = Table([row], colWidths=widths)
        story.append(table)
    else:
        missing = [os.path.basename(p) for p in paths if not os.path.exists(p)]
        story.append(Paragraph(f"<b>ERROR:</b> Could not find one or more files: {', '.join(missing)}.", styles["Normal"]))

### Discovers sequence names from images in the directory.
def get_sequence_names(image_dir, prefix):

    # Pattern: {prefix}_spectra_{sequence}.png, excluding spectra_gff
    pattern = re.compile(rf"^{re.escape(prefix)}_spectra_(.+)\.png$")
    sequences = []
    for f in os.listdir(image_dir):
        match = pattern.match(f)
        if match:
            seq_name = match.group(1)
            if not seq_name.startswith("gff_"):
                sequences.append(seq_name)
    return sorted(sequences)

### Adds the K-mer distribution plots section.
def add_kmer_distribution_section(story, image_dir, prefix, mer, styles, percentile, percentile_low=None, percentile_high=None):
    if percentile_low is None: percentile_low = percentile
    if percentile_high is None: percentile_high = percentile

    story.append(Paragraph(
        f"<b>K={mer} distributions:</b> K-mer prevalence (left) in raw data [x-axis, log-scale] against prevalence in assembled data [y-axis, log-scale]. "
        f"K-mer prevalence (right) when filtered for the bottom {percentile_low:.2f}% and top {percentile_high:.2f}% of k-mers by shift in abundance between datasets. "
        f"The scatter plot provides a global overview of how k-mer frequencies in the assembly match the raw sequencing reads. "
        f"Ideally, k-mers should cluster along the diagonal, with peaks representing the expected sequencing coverage.",
        styles["Normal"]))

    extreme_suffix = f"{percentile_low}pct" if percentile_low == percentile_high else f"L{percentile_low}_H{percentile_high}pct"
    paths = [
        os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_scatter.png"),
        os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_scatter_extreme_{extreme_suffix}.png")
    ]
    add_image_row(story, paths, [3.75 * inch] * 2, [4.5 * inch] * 2, styles)

    story.append(Spacer(1, 0.1 * inch))
    story.append(Paragraph(
        f"<b>K={mer} abundance shift:</b> log-fold change in k-mer representation between raw and assembled data. "
        f"Peaks in change should roughly corroborate the sequencing coverage of the genome. "
        f"Positive shifts indicate k-mers over-represented in the assembly, while negative shifts indicate k-mers that are more "
        f"abundant in the raw data than in the final assembly (potentially collapsed or missing regions).",
        styles["Normal"]))

    density_path = os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_density.png")
    add_safe_image(story, density_path, 6 * inch, 3 * inch, styles)

### Adds the abundance density and ECDF section.
def add_abundance_density_section(story, image_dir, prefix, mer, styles):
    story.append(PageBreak())
    story.append(Paragraph(
        f"<b>K={mer} abundance density:</b> Kernel density estimation (left) and violin plots (right) of k-mers in raw [blue] and assembled [orange] data. "
        f"Graphical estimations might not be smoothed depending on the data's composition. "
        f"These plots compare the overall distribution of k-mer multiplicities. A well-assembled genome should "
        f"closely mirror the distribution of the raw data, particularly at the primary coverage peak.",
        styles["Normal"]))

    paths = [
        os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_back2back_density.png"),
        os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_violin.png")
    ]
    add_image_row(story, paths, [3.75 * inch] * 2, [4.5 * inch] * 2, styles)

    story.append(Spacer(1, 0.1 * inch))
    story.append(Paragraph(
        f"<b>K={mer} empirical cumulative distribution (ECDF):</b> Measure of how many k-mers (and their cumulative probability) "
        f"are observed at each sequential log-fold change in frequency. "
        f"The ECDF helps identify the proportion of k-mers that fall within certain shift ranges, "
        f"providing a quantitative measure of assembly completeness and consistency.",
        styles["Normal"]))

    ecdf_path = os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_ecdf.png")
    add_safe_image(story, ecdf_path, 6 * inch, 3 * inch, styles)

### Adds the sequence-specific breakdown pages.
def add_sequence_breakdown_section(story, image_dir, prefix, mer, sequence_names, max_output, ngaps, bins, styles, spectra_dir=None, canonical=False, percentile_low=5, percentile_high=5):
    story.append(PageBreak())

    paragraph_text = (f"<b>Sequence-specific spectra breakdowns:</b> the following pages are a breakdown of spectra (K=3 mer distribution) "
                      f"and the K={mer} localization of exact kmer matches for highest {percentile_high:.2f}% and lowest {percentile_low:.2f}% in abundance change. "
                      f"Spectra plots show the 64 K=3 mers. Each page will have: high-abundance kmers(top), spectra (middle), "
                      f"and low abundance (bottom). K={mer} abundance plots use consistent scaling across sequences.")

    if ngaps:
        paragraph_text += f" Gaps in the sequence are denoted by vertical black bars at their positions. Bars are wider than actual gap size for visibility, please refer to output file {prefix}_ngaps.gff for precise sizes."
    if bins:
        paragraph_text += " Predicted shifts in sequence identity are labeled below Spectra 3-mers."
    if canonical:
        paragraph_text += " Secondary spectra plots represent distribution of canonical 3-mers."

    if len(sequence_names) > max_output:
        logging.warning(f"Too many contigs ({len(sequence_names)}) to tabulate. Only the first {max_output} will be output.")
        paragraph_text += f" There were too many sequences to reliably construct the report. Only the first {max_output} alphabetically are reported here. Rerun with '--max-out {len(sequence_names)} or find images in the output directory."

    story.append(Paragraph(paragraph_text, styles["Normal"]))
    story.append(Spacer(1, 0.5 * inch))

    circular_path = os.path.join(image_dir, f"{prefix}_standard_circular.png")
    add_safe_image(story, circular_path, 6.5 * inch, 6.5 * inch, styles, spacer=0)

    # Add legend below circular plot
    if spectra_dir:
        legend_path = os.path.join(spectra_dir, "includes", "Spectra-legend.png")
        if legend_path and os.path.exists(legend_path):
            story.append(Spacer(1, 0.2 * inch))
            legend_img = image_prep(legend_path, 5 * inch, 5 * inch)
            if legend_img:
                story.append(legend_img)

            legend_desc = ("<b>Spectra Legend:</b> The color keys above represent the 64 possible 3-mers (trinucleotides). "
                           "In the following sequence plots, these colors indicate the local composition and shifts in 3-mer distributions "
                           "across the assembly. Each color corresponds to a specific 3-mer as shown in the grid.")
            story.append(Spacer(1, 0.1 * inch))
            story.append(Paragraph(legend_desc, styles["Normal"]))

    for sequence in sequence_names[:max_output]:
        story.append(PageBreak())
        story.append(Paragraph(f"<b>Sequence {sequence}:</b>", styles["Normal"]))
        story.append(Spacer(1, 0.2 * inch))

        # 1. High abundance plot
        high_path = os.path.join(image_dir, f"{prefix}_mass_{sequence}_high.png")
        add_safe_image(story, high_path, 6.5 * inch, 4 * inch, styles, spacer=0)

        # 2. Non-canonical spectra plot
        story.append(Spacer(1, 0.1 * inch))
        spectra_path = os.path.join(image_dir, f"{prefix}_spectra_standard_{sequence}.png")
        add_safe_image(story, spectra_path, 6.5 * inch, 4 * inch, styles, spacer=0)

        # 3. Bins for non-canonical spectra
        if bins:
            gff_path = os.path.join(image_dir, f"{prefix}_spectra_standard_gff_{sequence}.png")
            if os.path.exists(gff_path):
                img = image_prep(gff_path, 6.5 * inch, 4 * inch)
                if img:
                    story.append(img)
            else:
                story.append(Paragraph(f"<b>ERROR:</b> Could not find file {os.path.basename(gff_path)}.", styles["Normal"]))

        # 4. Canonical spectra plot (if requested)
        if canonical:
            story.append(Spacer(1, 0.1 * inch))
            canon_path = os.path.join(image_dir, f"{prefix}_spectra_canonical_{sequence}.png")
            add_safe_image(story, canon_path, 6.5 * inch, 4 * inch, styles, spacer=0)

            # 5. Bins for canonical spectra
            if bins:
                canon_gff_path = os.path.join(image_dir, f"{prefix}_spectra_canonical_gff_{sequence}.png")
                if os.path.exists(canon_gff_path):
                    img = image_prep(canon_gff_path, 6.5 * inch, 4 * inch)
                    if img:
                        story.append(img)
                else:
                    story.append(Paragraph(f"<b>ERROR:</b> Could not find file {os.path.basename(canon_gff_path)}.", styles["Normal"]))

        # 6. Low abundance plot
        story.append(Spacer(1, 0.1 * inch))
        low_path = os.path.join(image_dir, f"{prefix}_mass_{sequence}_low.png")
        add_safe_image(story, low_path, 6.5 * inch, 4 * inch, styles, spacer=0)

### Main function to construct the PDF report.
def make_report(output_pdf, image_dir, mer, prefix, bins=False, ngaps=False, max_output=50, percentile=None, percentile_low=5, percentile_high=5,
                raw_files=None, assembly_file=None, counter=None, spectra_dir=None, canonical=False):
    if percentile is not None:
        percentile_low = percentile
        percentile_high = percentile
    doc = SimpleDocTemplate(output_pdf, pagesize=letter)
    doc.title = f'Spectra output report: {prefix}'
    story = []
    styles = getSampleStyleSheet()

    # Title Page
    logo_path = os.path.join(spectra_dir, "includes", "Spectra-Logo.png") if spectra_dir else None
    add_title_page(story, styles, logo_path, prefix, raw_files, assembly_file, counter, mer)

    # Introduction
    intro_text = ("<b>Spectra pipeline output report:</b> The following figures were auto-generated by the Spectra pipeline. "
                  "These figures show the relationship between k-mers in raw sequence data and in a genome assembly. "
                  "These figures were generated with a random sampling of k-mers.")
    story.append(Paragraph(intro_text, styles["Normal"]))
    story.append(Spacer(1, 0.1 * inch))

    # Sections
    add_kmer_distribution_section(story, image_dir, prefix, mer, styles, percentile=percentile_low, percentile_low=percentile_low, percentile_high=percentile_high)
    add_abundance_density_section(story, image_dir, prefix, mer, styles)

    sequence_names = get_sequence_names(image_dir, prefix)
    add_sequence_breakdown_section(story, image_dir, prefix, mer, sequence_names, max_output, ngaps, bins, styles, spectra_dir=spectra_dir, canonical=canonical, percentile_low=percentile_low, percentile_high=percentile_high)

    # Build PDF
    try:
        doc.build(story)
        logging.info(f"Report successfully generated: {output_pdf}")
    except Exception as e:
        logging.error(f"Failed to build PDF: {e}")

def main():
    parser = argparse.ArgumentParser(description='Spectra pipeline report writer')
    parser.add_argument('-i', '--image-directory', dest='directory', type=str, help='Input image directory', required=True)
    parser.add_argument('-o', '--output', dest='output', type=str, help='Output pdf filename', default='spectra_report.pdf')
    parser.add_argument('-m', '--mer-size', dest='mer_size', type=int, help='kmer size ran.', default=20)
    parser.add_argument('-n', '--n-gaps', dest='ngaps', action='store_true', help='Label gaps in the assembly in the final report', default=False)
    parser.add_argument('-b', '--bin-identify', dest='bins', action='store_true', help='Label bin regions in the genome assembly', default=False)
    parser.add_argument('--canonical', dest='canonical', action='store_true', help='Canonical spectra was generated', default=False)
    parser.add_argument('-p', '--prefix', dest='prefix', type=str, required=True)
    parser.add_argument('-e', '--percentile', dest='percentile', type=float, default=None)
    parser.add_argument('--percentile-low', dest='percentile_low', type=float, default=5)
    parser.add_argument('--percentile-high', dest='percentile_high', type=float, default=5)
    parser.add_argument('-x', '--max-output', dest='to_output', type=int, help='Contigs to include individual plots for, taken alphabetically.', default=50)
    parser.add_argument('-r', '--raw', dest='raw', nargs='+', help='Input raw fasta/fastq read file(s).', default=None)
    parser.add_argument('-a', '--assembled', dest='assembled', help='Input fasta assembly file.', default=None)
    parser.add_argument('-c', '--counter', dest='counter', help='K-mer counter used.', default=None)
    parser.add_argument('-s', '--spectra-dir', dest='spectra_dir', help='Path to Spectra directory for logo and legend.', default=None)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    make_report(
        output_pdf=args.output,
        image_dir=args.directory,
        mer=args.mer_size,
        prefix=args.prefix,
        bins=args.bins,
        ngaps=args.ngaps,
        max_output=args.to_output,
        percentile=args.percentile,
        percentile_low=args.percentile_low,
        percentile_high=args.percentile_high,
        raw_files=args.raw,
        assembly_file=args.assembled,
        counter=args.counter,
        spectra_dir=args.spectra_dir,
        canonical=args.canonical
    )

if __name__ == "__main__":
    main()

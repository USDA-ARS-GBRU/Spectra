#!/usr/bin/env python3
import logging
import os
import argparse
import re

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table
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
def add_kmer_distribution_section(story, image_dir, prefix, mer, styles, percentile):
    story.append(Paragraph(
        f"<b>K={mer} distributions:</b> Kmer prevalence (left) in raw data [x-axis, log-scale] against prevalence in assembled data [y-axis, log-scale]. Kmer prevalence (right) when filtered for the top and bottom {percentile}% of kmers by shift in abundance between datasets.",
        styles["Normal"]))

    paths = [
        os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_scatter.png"),
        os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_scatter_extreme_{percentile}pct.png")
    ]
    add_image_row(story, paths, [3.75 * inch] * 2, [4.5 * inch] * 2, styles)

    story.append(Spacer(1, 0.1 * inch))
    story.append(Paragraph(
        f"<b>K={mer} abundance shift:</b> log-fold change in kmer representation between raw and assembled data. Peaks in change should roughly corroborate the sequencing coverage of the genome.",
        styles["Normal"]))

    density_path = os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_density.png")
    add_safe_image(story, density_path, 6 * inch, 3 * inch, styles)

### Adds the abundance density and ECDF section.
def add_abundance_density_section(story, image_dir, prefix, mer, styles):
    story.append(PageBreak())
    story.append(Paragraph(
        f"<b>K={mer} abundance density:</b> Kernal density estimation (left) and violin plots (right) of kmers in raw [blue] and assembled [orange] data. Graphical estimations might not be smoothed depending on the data's composition.",
        styles["Normal"]))

    paths = [
        os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_back2back_density.png"),
        os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_violin.png")
    ]
    add_image_row(story, paths, [3.75 * inch] * 2, [4.5 * inch] * 2, styles)

    story.append(Spacer(1, 0.1 * inch))
    story.append(Paragraph(
        f"<b>K{mer} empirical cumulative distribution:</b> Measure of how many kmers (and their cumulative probability) are observed at each sequential log-fold change in frequency.",
        styles["Normal"]))

    ecdf_path = os.path.join(image_dir, f"{prefix}_kmer_comp_k{mer}_ecdf.png")
    add_safe_image(story, ecdf_path, 6 * inch, 3 * inch, styles)

### Adds the sequence-specific breakdown pages.
def add_sequence_breakdown_section(story, image_dir, prefix, mer, sequence_names, max_output, ngaps, bins, styles):
    story.append(PageBreak())

    paragraph_text = (f"<b>Sequence-specific spectra breakdowns:</b> the following pages are a breakdown of spectra (K=3 mer distribution) "
                      f"and the K={mer} localization of exact kmer matches for highest 5% and lowest 5% in abundance change. "
                      f"Spectra plots show the 64 K=3 mers. Each page will have: high-abundance kmers(top), spectra (middle), "
                      f"and low abundance (bottom). K={mer} abundance plots are not to scale.")

    if ngaps:
        paragraph_text += " Gaps in the sequence are denoted by solid black bars at their positions."
    if bins:
        paragraph_text += " Predicted shifts in sequence identity are labeled below Spectra 3-mers"

    if len(sequence_names) > max_output:
        logging.warning(f"Too many contigs ({len(sequence_names)}) to tabulate. Only the first {max_output} will be output.")
        paragraph_text += f" There were too many sequences to reliably construct the report. Only the first {max_output} alphabetically are reported here."

    story.append(Paragraph(paragraph_text, styles["Normal"]))
    story.append(Spacer(1, 0.5 * inch))

    circular_path = os.path.join(image_dir, f"{prefix}_circular.png")
    add_safe_image(story, circular_path, 6.5 * inch, 6.5 * inch, styles, spacer=0)

    for sequence in sequence_names[:max_output]:
        story.append(PageBreak())
        story.append(Paragraph(f"<b>Sequence {sequence}:</b>", styles["Normal"]))
        story.append(Spacer(1, 0.3 * inch))

        # High abundance plot
        high_path = os.path.join(image_dir, f"{prefix}_mass_{sequence}_high.png")
        add_safe_image(story, high_path, 6.5 * inch, 4 * inch, styles, spacer=0)

        # Spectra plot
        story.append(Spacer(1, 0.1 * inch))
        spectra_path = os.path.join(image_dir, f"{prefix}_spectra_{sequence}.png")
        add_safe_image(story, spectra_path, 6.5 * inch, 4 * inch, styles, spacer=0)

        # GFF/Bins plot
        if bins:
            story.append(Spacer(1, 0.1 * inch))
            gff_path = os.path.join(image_dir, f"{prefix}_spectra_gff_{sequence}.png")
            if os.path.exists(gff_path):
                # Using a table to match the original padding logic if needed,
                # but simplified for now.
                img = image_prep(gff_path, 6.6 * inch, 4 * inch)
                if img:
                    row = Table([[Spacer(0.003 * inch, 0.01 * inch), img]])
                    story.append(row)
            else:
                story.append(Paragraph(f"<b>ERROR:</b> Could not find file {os.path.basename(gff_path)}.", styles["Normal"]))

        # Low abundance plot
        story.append(Spacer(1, 0.1 * inch))
        low_path = os.path.join(image_dir, f"{prefix}_mass_{sequence}_low.png")
        add_safe_image(story, low_path, 6.5 * inch, 4 * inch, styles, spacer=0)

### Main function to construct the PDF report.
def make_report(output_pdf, image_dir, mer, prefix, bins=False, ngaps=False, max_output=50, percentile=5):
    doc = SimpleDocTemplate(output_pdf, pagesize=letter)
    doc.title = f'Spectra output report: {prefix}'
    story = []
    styles = getSampleStyleSheet()

    # Introduction
    intro_text = ("<b>Spectra pipeline output report:</b> The following figures were auto-generated by the Spectra pipeline. "
                  "These figures show the relationship between kmers in raw sequence data and in a genome assembly. "
                  "These figures were generated with a random sampling of kmers.")
    story.append(Paragraph(intro_text, styles["Normal"]))
    story.append(Spacer(1, 0.1 * inch))

    # Sections
    add_kmer_distribution_section(story, image_dir, prefix, mer, styles, percentile)
    add_abundance_density_section(story, image_dir, prefix, mer, styles)

    sequence_names = get_sequence_names(image_dir, prefix)
    add_sequence_breakdown_section(story, image_dir, prefix, mer, sequence_names, max_output, ngaps, bins, styles)

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
    parser.add_argument('-p', '--prefix', dest='prefix', type=str, required=True)
    parser.add_argument('-e', '--percentile', dest='percentile', type=int, default=5)
    parser.add_argument('-x', '--max_output', dest='to_output', type=int, help='Contigs to include individual plots for, taken alphabetically.', default=50)
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
        percentile=args.percentile
    )

if __name__ == "__main__":
    main()

from pathlib import Path
from io import BytesIO
from pypdf import PdfReader, PdfWriter, PageObject, Transformation
from pypdf.generic import RectangleObject, ContentStream, NameObject, NumberObject
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4

base_symbol = Path("/home/global_EOCRC_microbiology")
output_dir = (
    base_symbol / "downstream_analysis/output/submission_version_261102/metabolomics_v35"
)
source_a = (
    base_symbol / "FUSCC_metabonomics_analysis/v35/references/functional_A_original.pdf"
)
source_be = output_dir / "plots/Fig4_B_E_A4.pdf"

# The original functional PDF is composed here; its raw data are unavailable.
a_reader = PdfReader(source_a)
be_reader = PdfReader(source_be)
assert len(a_reader.pages) == len(be_reader.pages) == 1
a_page = a_reader.pages[0]
be_page = be_reader.pages[0]

a_content = ContentStream(a_page.get_contents(), a_reader)
operations = a_content.operations
legend_rects = [
    i
    for i, (args, op) in enumerate(operations)
    if op == b"re"
    and len(args) == 4
    and abs(float(args[0]) - 502.785) < 0.01
    and abs(float(args[1]) - 500.840) < 0.01
]
title_matrices = [
    i
    for i, (args, op) in enumerate(operations)
    if op == b"Tm"
    and len(args) == 6
    and abs(float(args[4]) - 300.410156) < 0.01
    and abs(float(args[5]) - 15.102539) < 0.01
]
assert len(legend_rects) == len(title_matrices) == 1
legend_start, legend_end = legend_rects[0] - 1, title_matrices[0]
assert operations[legend_start][1] == b"rg" and legend_start < legend_end
a_content.operations = (
    operations[:legend_start]
    + [
        ([], b"BT"),
        ([NameObject("/f-1-0"), NumberObject(1)], b"Tf"),
    ]
    + operations[legend_end:]
)
a_page[NameObject("/Contents")] = a_content

cohort_items = [("EO", (0.933333, 0, 0)), ("LO", (0.231373, 0.286275, 0.572549))]
legend_items = [
    ("Significant (q < 0.1)", (0.988235, 0.552941, 0.349020)),
    ("Non-significant", (0.650980, 0.337255, 0.156863)),
    ("NA", (0.850980, 0.850980, 0.850980)),
    ("Consistent (QC_pass)", (0.596078, 0.305882, 0.639216)),
]
a_page.cropbox = RectangleObject([0, 612 - 496, 916, 612])

# Compose the original panel on the final A4 canvas.
width, height = A4
mm = 72 / 25.4
a_target_width = 160 * mm
a_scale = a_target_width / float(a_page.mediabox.width)
a_height = float(a_page.mediabox.height) * a_scale
page = PageObject.create_blank_page(width=width, height=height)
page.merge_page(be_page)
page.merge_transformed_page(
    a_page,
    Transformation().scale(a_scale).translate(
        4 * mm,
        height - 5 * mm - a_height,
    ),
)

overlay = BytesIO()
c = canvas.Canvas(overlay, pagesize=(width, height))
c.setFont("Helvetica-Bold", 12)
c.drawString(1 * mm, height - 6 * mm, "A")
c.setFont("Helvetica-Bold", 6.8)
c.drawString(168 * mm, height - 31 * mm, "Cohort Group")
c.setFont("Helvetica", 6.5)

for i, (label, color) in enumerate(cohort_items):
    legend_y = height - (36 + i * 5) * mm
    c.setFillColorRGB(*color)
    c.rect(168 * mm, legend_y - 0.45 * mm, 2.8 * mm, 2.8 * mm, fill=1, stroke=0)
    c.setFillColorRGB(0, 0, 0)
    c.drawString(174 * mm, legend_y, label)

c.setFont("Helvetica-Bold", 6.8)
c.drawString(168 * mm, height - 52 * mm, "Statistical Metrics")
c.setFont("Helvetica", 6.5)

for i, (label, color) in enumerate(legend_items):
    legend_y = height - (57 + i * 5) * mm
    c.setFillColorRGB(*color)
    c.rect(168.8 * mm, legend_y + 0.2 * mm, 1.25 * mm, 1.25 * mm, fill=1, stroke=0)
    c.setFillColorRGB(0, 0, 0)
    c.drawString(174 * mm, legend_y, label)

c.save()
overlay.seek(0)
page.merge_page(PdfReader(overlay).pages[0])
writer = PdfWriter()
writer.add_page(page)
writer.add_metadata(
    {
        "/Title": (
            "Fig4: original microbial functional panel and revised FUSCC metabolomics"
        ),
        "/Subject": (
            "A plot results preserved; both legends vertical on right; "
            "B-D taller, C wider, D narrower."
        ),
    }
)
target = output_dir / "plots/Fig4_A_E_v35.pdf"
with target.open("wb") as stream:
    writer.write(stream)

# Match the final supplementary PDF page box to A4.
supplement_target = output_dir / "plots/supplementary_figure_composite_v35.pdf"
supplement_reader = PdfReader(supplement_target)
assert len(supplement_reader.pages) == 1
needs_supplement_normalization = any(
    abs(float(actual) - expected) > 1e-6
    for actual, expected in zip(
        supplement_reader.pages[0].mediabox,
        [0, 0, width, height],
    )
)
supplement_writer = PdfWriter()
supplement_page = supplement_reader.pages[0]
supplement_page.mediabox = RectangleObject([0, 0, width, height])
supplement_page.cropbox = RectangleObject([0, 0, width, height])
supplement_writer.add_page(supplement_page)
supplement_writer.add_metadata(
    {"/Title": "FUSCC supplementary composite: A4 portrait; panel captions removed"}
)

if needs_supplement_normalization:
    with supplement_target.open("wb") as stream:
        supplement_writer.write(stream)

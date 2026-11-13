# Figure 1 artwork and README image

The README displays [Figure_1.png](Figure_1.png), rendered at 2592 × 1080 pixels from the user-supplied [PowerPoint PDF export](Figure_1.pdf). The editable [source/Fig1.pptx](source/Fig1.pptx) is preserved unchanged. The PDF render was visually inspected: Greek alpha/beta, MSA braces, labels, contact maps and structure panels are present and legible.

The [provenance record](provenance.json) contains hashes for the PDF, PNG and editable source, renderer versions and the prior manuscript-preview record. The supplied PDF is retained byte-for-byte. No scientific values, artwork content or crop boundaries were changed.

## Recreate the PNG

From the repository root, with `pypdfium2==5.11.0` and `Pillow==12.2.0` installed (optional artwork-rendering dependencies, separate from the analysis requirements):

```bash
python - <<'PYCODE'
from pathlib import Path
import pypdfium2 as pdfium
root = Path('figures/Figure_1')
with pdfium.PdfDocument(root / 'Figure_1.pdf') as pdf:
    assert len(pdf) == 1
    pdf[0].render(scale=2).to_pil().save(root / 'Figure_1.png')
PYCODE
```

The earlier `Figure_1_manuscript_preview.png` is retained in the ignored local archive, with its original hash recorded in provenance. It is no longer embedded in the root README. The failed Quick Look render of the PPT was not added; the PDF export resolves its Greek-letter and brace rendering problems.

The editable layout contains embedded contact and structure artwork. Its availability does not by itself reproduce the illustrative calculations or AlphaFold inference.

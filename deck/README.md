# Study-plan deck

`build.js` generates `pan-cancer-spatial-single-cell-plan.pptx` (31 slides): the
data-acquisition, analysis and validation plan for the pan-cancer paired
spatial + single-cell study this repository implements.

```bash
npm install pptxgenjs        # once
node build.js                # -> pan-cancer-spatial-single-cell-plan.pptx
python3 qa.py pan-cancer-spatial-single-cell-plan.pptx
```

The deck is generated rather than hand-authored so that every figure and number
can be traced to its source and regenerated when the underlying facts change.

**Provenance of the numbers.** Dataset counts on the "live GEO survey" slide were
queried from NCBI E-utilities (`db=gds`) on 9 September 2026 and are keyword-search
upper bounds, not a curated cohort — the slide says so. Method-performance claims
are sourced to the published benchmarks cited in each slide's footer.

**QA.** `qa.py` checks slide-bound violations, text overflow (average-advance-width
model), text-on-text overlap, and edge margins. It substitutes for render-based
visual QA, which needs LibreOffice; install it and use the pptx skill's
`soffice.py` + `pdftoppm` path if you want image-level inspection too.

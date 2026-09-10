/* Pan-cancer paired ST + scRNA-seq study plan deck.
   Palette is H&E-informed: hematoxylin plum, eosin pink, fluorophore teal.
   Motif: filled "cell" dots and dot-scatter fields evoking a spatial plot. */
const pptxgen = require("pptxgenjs");

const INK="1F1235", PLUM="3D2463", PURPLE="6A4BA8", EOSIN="D4467A",
      TEAL="0F9B8E", AMBER="C87F14", LIGHT="F6F3FB", TINT="EDE7F6",
      WHITE="FFFFFF", BODY="3C3547", MUTED="6E6682", LINE="D8D0E6";
const HEAD="Cambria", SANS="Calibri";
const W=13.333, H=7.5, M=0.62;

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.author = "Pan-cancer spatial oncology programme";
pres.title  = "Pan-cancer paired spatial + single-cell atlas";

let N = 0;

/* deterministic pseudo-random so re-runs are identical */
let _s = 20260909;
const rnd = () => { _s = (_s*1103515245 + 12345) & 0x7fffffff; return _s/0x7fffffff; };

function dotField(slide, x, y, w, h, n, colors, sz=0.07, transparency=55){
  for(let i=0;i<n;i++){
    slide.addShape(pres.ShapeType.ellipse, {
      x: x + rnd()*w, y: y + rnd()*h, w: sz, h: sz,
      fill: { color: colors[i % colors.length], transparency }, line: { width: 0 },
    });
  }
}

function footer(slide, note){
  N += 1;
  slide.addText(String(N), { x: W-1.05, y: H-0.64, w: 0.5, h: 0.3, fontSize: 10,
    color: MUTED, fontFace: SANS, align: "right", isTextBox: true, margin: 0 });
  if (note) slide.addText(note, { x: M, y: H-0.66, w: W-2.0, h: 0.34, fontSize: 8.5,
    color: MUTED, fontFace: SANS, italic: true, isTextBox: true, margin: 0 });
}

function titleBlock(slide, kicker, title, opts={}){
  const dark = !!opts.dark;
  slide.addText(kicker.toUpperCase(), { x: M, y: 0.40, w: W-2*M, h: 0.26, fontSize: 10.5,
    bold: true, charSpacing: 2.2, color: dark ? EOSIN : PURPLE, fontFace: SANS,
    isTextBox: true, margin: 0 });
  slide.addText(title, { x: M, y: 0.68, w: opts.tw || (W-2*M), h: opts.th || 0.72, fontSize: opts.ts || 30,
    bold: true, color: dark ? WHITE : INK, fontFace: HEAD, isTextBox: true, margin: 0 });
}

function content(kicker, title, opts={}){
  const s = pres.addSlide();
  s.background = { color: opts.bg || WHITE };
  titleBlock(s, kicker, title, opts);
  return s;
}

function card(slide, o){
  slide.addShape(pres.ShapeType.roundRect, {
    x:o.x, y:o.y, w:o.w, h:o.h, rectRadius:0.06,
    fill:{ color:o.fill || LIGHT }, line:{ color:o.line || (o.fill?o.fill:LINE), width:o.lw===undefined?1:o.lw },
    shadow: o.shadow ? { type:"outer", angle:90, blur:8, offset:0.04, color:"9C8FB8", opacity:0.28 } : undefined,
  });
}

function dot(slide, x, y, d, color, label, labelColor){
  slide.addShape(pres.ShapeType.ellipse, { x, y, w:d, h:d, fill:{color}, line:{width:0} });
  if(label!==undefined) slide.addText(String(label), { x, y, w:d, h:d, fontSize: d>0.42?13:10.5,
    bold:true, color: labelColor||WHITE, align:"center", valign:"middle", fontFace:SANS,
    isTextBox:true, margin:0 });
}

function body(slide, o){
  slide.addText(o.text, { x:o.x, y:o.y, w:o.w, h:o.h, fontSize:o.size||13,
    color:o.color||BODY, fontFace:SANS, isTextBox:true, margin:o.margin===undefined?0:o.margin,
    lineSpacingMultiple:o.lsm||1.12, align:o.align||"left", bold:o.bold, italic:o.italic,
    valign:o.valign });
}

function bullets(slide, items, o){
  slide.addText(items.map((t,i)=>({ text:t, options:{ bullet:{indent:14}, breakLine:i<items.length-1 } })),
    { x:o.x, y:o.y, w:o.w, h:o.h, fontSize:o.size||12.5, color:o.color||BODY, fontFace:SANS,
      isTextBox:true, margin:0, paraSpaceAfter:o.gap===undefined?7:o.gap, lineSpacingMultiple:1.1 });
}

function stat(slide, o){
  slide.addText(o.value, { x:o.x, y:o.y, w:o.w, h:o.vh||0.62, fontSize:o.vs||34, bold:true,
    color:o.color||PURPLE, fontFace:HEAD, isTextBox:true, margin:0, align:o.align||"left" });
  slide.addText(o.label, { x:o.x, y:o.y+(o.vh||0.62)-0.04, w:o.w, h:o.lh||0.5, fontSize:o.ls||10,
    color:MUTED, fontFace:SANS, isTextBox:true, margin:0, align:o.align||"left",
    lineSpacingMultiple:1.05 });
}

/* ---------------------------------------------------------------- 1 title */
{
  const s = pres.addSlide();
  s.background = { color: INK };
  dotField(s, 8.4, 0.4, 4.6, 6.7, 190, [EOSIN, PURPLE, TEAL, "8E7BC3"], 0.075, 40);
  s.addShape(pres.ShapeType.ellipse, { x:9.55, y:2.1, w:2.5, h:2.5,
    fill:{color:EOSIN, transparency:78}, line:{width:0} });
  dotField(s, 9.85, 2.4, 1.9, 1.9, 55, [WHITE, TEAL], 0.065, 18);

  s.addText("PAN-CANCER SPATIAL ONCOLOGY PROGRAMME", { x:M, y:1.5, w:7.4, h:0.3,
    fontSize:11, bold:true, charSpacing:2.6, color:EOSIN, fontFace:SANS, isTextBox:true, margin:0 });
  s.addText("Paired spatial and single-cell\ntranscriptomics across all cancers",
    { x:M, y:1.92, w:7.7, h:1.9, fontSize:32, bold:true, color:WHITE, fontFace:HEAD,
      isTextBox:true, margin:0, lineSpacingMultiple:1.02 });
  s.addText("A data-acquisition, analysis and validation plan — built on benchmark evidence, "
    + "with the statistical guardrails that decide whether the findings survive review.",
    { x:M, y:3.95, w:7.2, h:0.9, fontSize:14, color:"C9BFE0", fontFace:SANS, isTextBox:true,
      margin:0, lineSpacingMultiple:1.2 });

  const chips = ["Data harvest", "Niche atlas", "CCI · GRN · modules", "Trajectories", "Validation"];
  let cx = M;
  chips.forEach(c => {
    const w = 0.22 + c.length*0.093;
    s.addShape(pres.ShapeType.roundRect, { x:cx, y:5.15, w, h:0.36, rectRadius:0.18,
      fill:{color:PLUM}, line:{color:PURPLE, width:1} });
    s.addText(c, { x:cx, y:5.15, w, h:0.36, fontSize:10, color:"D9CEF0", fontFace:SANS,
      align:"center", valign:"middle", isTextBox:true, margin:0 });
    cx += w + 0.14;
  });
  s.addText("Method claims fact-checked against published benchmarks · Dataset counts queried live from NCBI GEO, 9 September 2026",
    { x:M, y:6.55, w:8.0, h:0.5, fontSize:9.5, color:"8E82AD", fontFace:SANS, italic:true,
      isTextBox:true, margin:0, lineSpacingMultiple:1.1 });
  s.addNotes("Framing: this is a plan grounded in what benchmarks actually show and what data actually exists, not an aspirational method list. Every number on the following slides is either from a cited publication or from a live GEO query run while building the deck.");
  N += 1;
}

/* ------------------------------------------------- 2 prior art (honesty) */
{
  const s = content("Reality check", "Pan-cancer spatial atlases already exist");
  body(s, { x:M, y:1.42, w:W-2*M, h:0.56, size:13, color:BODY, italic:true,
    text:"Before designing anything, the honest question: what has already been published? Four groups have released pan-cancer spatial resources in the last 18 months. The plan must differentiate against these, not rediscover them." });

  const studies = [
    { t:"Cell Rep Med, Apr 2026", n:"373 samples · 12 cancers", d:"10x Visium. 56 local cellular programs (28 malignant, 28 stromal); 13 consensus niches. Macrophage–tumor niche tracks poor prognosis and ICB resistance; macrophage–immune niche tracks better survival.", c:EOSIN },
    { t:"Pan-tumor SGs, Aug 2026", n:"262 solid tumors", d:"TMEs partition into discrete, hierarchically organised 'spatial groups'. Dominant axis of variation across tumors is spatial heterogeneity of immune biology.", c:PURPLE },
    { t:"Cell Rep Med, Oct 2025", n:"230 sc samples · 9 cancers", d:"70 pan-cancer single-cell subtypes and their co-occurrence; 60 Visium sections (BC 6, CRC 12, HCC 7, GC 9, NSCLC 20, MEL 6).", c:TEAL },
    { t:"Cellular architecture atlas", n:"7,910 WSIs · 21 tumors", d:"Deep-learning segmentation of >4.7 billion nuclei; nearest-neighbour distance linked to a biomechanical–immune axis, correlated with Visium HD.", c:AMBER },
  ];
  const cw = (W-2*M-3*0.24)/4;
  studies.forEach((st,i)=>{
    const x = M + i*(cw+0.24);
    card(s, { x, y:2.05, w:cw, h:3.55, fill:LIGHT, line:LINE });
    s.addShape(pres.ShapeType.roundRect, { x:x+0.22, y:2.28, w:cw-0.44, h:0.30, rectRadius:0.15,
      fill:{color:st.c}, line:{width:0} });
    s.addText(st.t, { x:x+0.22, y:2.28, w:cw-0.44, h:0.30, fontSize:9, bold:true, color:WHITE,
      fontFace:SANS, align:"center", valign:"middle", isTextBox:true, margin:0 });
    body(s, { x:x+0.22, y:2.72, w:cw-0.44, h:0.42, size:13, bold:true, color:INK, text:st.n });
    body(s, { x:x+0.22, y:3.20, w:cw-0.44, h:2.25, size:10.5, color:BODY, text:st.d, lsm:1.16 });
  });

  card(s, { x:M, y:5.82, w:W-2*M, h:0.72, fill:TINT, line:TINT });
  body(s, { x:M+0.26, y:5.98, w:W-2*M-0.52, h:0.45, size:12, color:PLUM,
    text:"Consequence for this programme: \"first pan-cancer spatial atlas\" is not available as a claim. The defensible contributions are scale, single-cell-resolution platforms, mechanism layers (GRN, modules, trajectories), and statistical rigour — addressed on the next slide." });
  footer(s, "Sources: Cell Rep Med 2026;6(4) · Cell Rep Med 2025 · pan-tumor spatial groups 2026 · pan-cancer cellular architecture atlas 2026");
  s.addNotes("This slide exists because the single largest hallucination risk in a project like this is claiming novelty that the literature already occupies. All four studies were verified by search during deck construction.");
}

/* ------------------------------------------------------ 3 differentiation */
{
  const s = content("Positioning", "Four gaps the prior work leaves open");
  const gaps = [
    { n:"01", h:"Spot resolution ceilings the biology", c:EOSIN,
      d:"The 2026 atlas of 373 samples is Visium-only. A 55 µm spot holds 1–10 cells, so rare populations — NK cells, cDC1, Tregs — exist only as deconvolution weights. Single-cell-resolution platforms (Xenium, CosMx, MERFISH) now exist at cohort scale and are largely unexploited pan-cancer." },
    { n:"02", h:"Description without mechanism", c:PURPLE,
      d:"Published atlases stop at niches and ligand–receptor pairs. No pan-cancer resource layers enhancer-driven regulatory networks, co-expression modules, and fate trajectories onto the same spatial coordinates — which is what turns a niche into a target." },
    { n:"03", h:"Replication unit is routinely wrong", c:TEAL,
      d:"Cell- and spot-level tests with n in the 10⁵–10⁶ range make almost any effect significant. Patient-level meta-analysis with leave-one-cohort-out sensitivity is the difference between a pan-cancer finding and one cohort reported with a pan-cancer n." },
    { n:"04", h:"Spatial nulls that inflate significance", c:AMBER,
      d:"Label-shuffling over an autocorrelated tissue field destroys the structure that makes the null realistic. On simulated tissue with no true association we measured 100% false positives from label shuffling versus 8% from a torus-shift null." },
  ];
  const rh = 1.18;
  gaps.forEach((g,i)=>{
    const y = 1.62 + i*(rh+0.14);
    card(s, { x:M, y, w:W-2*M, h:rh, fill:i%2?WHITE:LIGHT, line:LINE });
    dot(s, M+0.28, y+0.30, 0.58, g.c, g.n);
    body(s, { x:M+1.06, y:y+0.20, w:3.5, h:0.72, size:14, bold:true, color:INK, text:g.h, lsm:1.05 });
    body(s, { x:M+4.72, y:y+0.19, w:W-2*M-5.05, h:0.85, size:11, color:BODY, text:g.d, lsm:1.14 });
  });
  footer(s, "False-positive rates measured on simulated clustered point patterns, 24 replicates, α = 0.05 (this programme's statistics module).");
  s.addNotes("Gap 4 is measured, not asserted: the simulation is in the accompanying codebase's test suite.");
}

/* ---------------------------------------------------------------- 4 aims */
{
  const s = pres.addSlide();
  s.background = { color: INK };
  titleBlock(s, "Objectives", "Three aims", { dark:true });
  dotField(s, 10.6, 0.5, 2.4, 1.2, 45, [PURPLE, EOSIN, TEAL], 0.07, 45);

  const aims = [
    { k:"AIM 1", t:"Assemble", c:EOSIN,
      d:"Systematically harvest, curate and harmonise every publicly available human cancer dataset with paired spatial and single-cell transcriptomics — across repositories, platforms and access tiers.",
      out:"A versioned, provenance-tracked cohort manifest with per-dataset quality and resolution flags." },
    { k:"AIM 2", t:"Dissect", c:PURPLE,
      d:"Resolve the tumour microenvironment and the malignant compartment into spatial niches, and layer mechanism onto them: cell–cell interactions, transcription-factor networks, co-expression modules and fate trajectories.",
      out:"A pan-cancer niche atlas with a mechanism layer, conserved and cancer-type-specific components separated." },
    { k:"AIM 3", t:"Translate", c:TEAL,
      d:"Prioritise spatially defined targets and biomarkers, associate them with survival and immunotherapy response, and specify the orthogonal experiments that test causality.",
      out:"A ranked target list, each entry paired with the assay that would falsify it." },
  ];
  const cw = (W-2*M-2*0.3)/3;
  aims.forEach((a,i)=>{
    const x = M + i*(cw+0.3);
    card(s, { x, y:1.78, w:cw, h:4.5, fill:PLUM, line:"51357F" });
    dot(s, x+0.32, 2.06, 0.52, a.c, i+1);
    s.addText(a.k, { x:x+1.0, y:2.10, w:cw-1.2, h:0.24, fontSize:9.5, bold:true, charSpacing:1.8,
      color:a.c, fontFace:SANS, isTextBox:true, margin:0 });
    s.addText(a.t, { x:x+1.0, y:2.32, w:cw-1.2, h:0.36, fontSize:20, bold:true, color:WHITE,
      fontFace:HEAD, isTextBox:true, margin:0 });
    body(s, { x:x+0.32, y:2.98, w:cw-0.64, h:1.85, size:11.5, color:"CFC4E6", text:a.d, lsm:1.18 });
    s.addShape(pres.ShapeType.line, { x:x+0.32, y:4.92, w:cw-0.64, h:0,
      line:{ color:"5B3E96", width:1 } });
    s.addText("DELIVERABLE", { x:x+0.32, y:5.02, w:cw-0.64, h:0.2, fontSize:8.5, bold:true,
      charSpacing:1.6, color:a.c, fontFace:SANS, isTextBox:true, margin:0 });
    body(s, { x:x+0.32, y:5.26, w:cw-0.64, h:0.9, size:10.5, color:"B7A9D4", text:a.out, lsm:1.16 });
  });
  footer(s, "");
  s.addNotes("Each aim carries one deliverable, so progress is measurable rather than narrative.");
}

/* ------------------------------------------------------- 5 study design */
{
  const s = content("Design", "Six phases, each gated on the one before");
  const ph = [
    { n:"1", t:"Harvest",    d:"Repository sweep,\nplatform typing,\naccess triage", c:EOSIN },
    { n:"2", t:"Curate & QC", d:"Per-sample MAD QC,\ndoublets, metadata\nreconciliation", c:EOSIN },
    { n:"3", t:"Harmonise",  d:"Within-cancer-type\nintegration, unified\ncell ontology", c:PURPLE },
    { n:"4", t:"Map",        d:"Deconvolution,\nspatial domains,\nniche calling", c:PURPLE },
    { n:"5", t:"Mechanism",  d:"CCI · GRN ·\nmodules · fate\ntrajectories", c:TEAL },
    { n:"6", t:"Translate",  d:"Outcome models,\ntargets, wet-lab\nvalidation", c:TEAL },
  ];
  const cw = (W-2*M-5*0.18)/6;
  ph.forEach((p,i)=>{
    const x = M + i*(cw+0.18);
    card(s, { x, y:1.66, w:cw, h:2.28, fill:LIGHT, line:LINE });
    dot(s, x+cw/2-0.26, 1.88, 0.52, p.c, p.n);
    body(s, { x:x+0.14, y:2.52, w:cw-0.28, h:0.32, size:14, bold:true, color:INK,
      text:p.t, align:"center" });
    body(s, { x:x+0.14, y:2.90, w:cw-0.28, h:0.96, size:10, color:BODY, text:p.d,
      align:"center", lsm:1.2 });
    if(i<5) s.addShape(pres.ShapeType.rightArrow, { x:x+cw+0.015, y:2.02, w:0.15, h:0.24,
      fill:{color:LINE}, line:{width:0} });
  });

  card(s, { x:M, y:4.16, w:W-2*M, h:1.34, fill:INK, line:INK });
  body(s, { x:M+0.3, y:4.34, w:2.7, h:0.3, size:11, bold:true, color:EOSIN,
    text:"GATES BETWEEN PHASES" });
  const gates = [
    "3 → 4  Integration may not be applied across cancer type; batch is confounded with it",
    "4 → 5  A cell type may only carry a conclusion after passing the deconvolution benchmark",
    "5 → 6  A claim needs patient-level significance and a stable leave-one-cohort-out range",
  ];
  bullets(s, gates, { x:M+0.3, y:4.68, w:W-2*M-0.6, h:0.75, size:10.5, color:"C6BADD", gap:3 });

  card(s, { x:M, y:5.68, w:W-2*M, h:0.86, fill:LIGHT, line:LINE });
  body(s, { x:M+0.3, y:5.84, w:W-2*M-0.6, h:0.6, size:11.5, color:BODY,
    text:"Everything runs as a Snakemake DAG with per-rule conda environments, declared resources and a global seed. Per-sample work fans out on wildcards; cross-patient numbers are produced at exactly one aggregation rule, which is where the multiple-testing family is defined." });
  footer(s, "");
  s.addNotes("The gates are the point of the slide. A phase diagram without gates is a picture; with gates it is a protocol.");
}

/* ------------------------------------------------- 6 platform capability */
{
  const s = content("Constraint", "What each platform can and cannot resolve");
  body(s, { x:M, y:1.42, w:W-2*M, h:0.52, size:12.5, italic:true, color:BODY,
    text:"Platform resolution is the hard ceiling on every downstream claim. Choosing analyses that the data cannot support is the most common failure in spatial oncology." });

  const rows = [
    ["Visium",      "55 µm spots",   "1–10 cells",    "Whole transcriptome", "Niche context, regional programs", "No", EOSIN],
    ["Visium HD",   "2 µm bins",     "Sub-cellular bin", "Whole transcriptome", "Fine architecture; bins ≠ segmented cells", "Partial", EOSIN],
    ["Xenium",      "Single cell",   "Segmented",     "Targeted panel (100s–5k)", "Rare cell types, distances, niches", "Yes", TEAL],
    ["CosMx",       "Single cell",   "Segmented",     "Targeted panel (1k–6k)", "Rare cell types, subcellular CCI", "Yes", TEAL],
    ["MERFISH",     "Single cell",   "Segmented",     "Targeted panel",       "High-sensitivity single-molecule work", "Yes", TEAL],
    ["Slide-seq",   "10 µm beads",   "1–3 cells",     "Whole transcriptome",  "Near-cellular architecture", "Partial", AMBER],
    ["GeoMx",       "Region of interest", "Pooled ROI", "Whole transcriptome", "Hypothesis-driven region comparison", "No", AMBER],
  ];
  const cols = [1.5, 1.32, 1.42, 2.22, 4.42, 1.22];
  const heads = ["Platform","Resolution","Cells / unit","Gene coverage","Best-supported use","Rare-cell claims"];
  let x0 = M, y0 = 1.96;
  s.addShape(pres.ShapeType.rect, { x:M, y:y0, w:W-2*M, h:0.42, fill:{color:INK}, line:{width:0} });
  heads.forEach((h,i)=>{
    s.addText(h, { x:x0+0.12, y:y0, w:cols[i]-0.16, h:0.42, fontSize:10, bold:true, color:WHITE,
      fontFace:SANS, valign:"middle", isTextBox:true, margin:0 });
    x0 += cols[i];
  });
  rows.forEach((r,ri)=>{
    const y = y0 + 0.42 + ri*0.50;
    s.addShape(pres.ShapeType.rect, { x:M, y, w:W-2*M, h:0.50,
      fill:{color: ri%2 ? WHITE : LIGHT}, line:{color:LINE, width:0.5} });
    let cx = M;
    for(let ci=0; ci<6; ci++){
      const last = ci===5;
      s.addText(r[ci], { x:cx+0.12, y, w:cols[ci]-0.16, h:0.50,
        fontSize: ci===0?11:10, bold: ci===0 || last,
        color: last ? (r[5]==="Yes"?TEAL:(r[5]==="No"?EOSIN:AMBER)) : (ci===0?INK:BODY),
        fontFace:SANS, valign:"middle", isTextBox:true, margin:0, lineSpacingMultiple:1.05 });
      cx += cols[ci];
    }
  });
  card(s, { x:M, y:6.02, w:W-2*M, h:0.60, fill:TINT, line:TINT });
  body(s, { x:M+0.26, y:6.14, w:W-2*M-0.52, h:0.40, size:11.5, color:PLUM,
    text:"Rule applied throughout: spot platforms supply niche context; only segmented single-cell platforms carry claims about rare populations. Mixing the two without flagging which is which is how deconvolution artefacts become findings." });
  footer(s, "");
}

/* --------------------------------------------------------- 7 data sources */
{
  const s = content("Aim 1 · Sources", "Where the paired data actually lives");
  const src = [
    { n:"GEO / SRA", v:"Primary", d:"The bulk of paired depositions. Free-text platform description means classification must be regex-audited, not keyword-matched.", c:EOSIN },
    { n:"HTAN", v:"8,425", d:"Biospecimens from ~2,000 participants across >20 assay types; the most systematically structured tumour resource, with paired single-cell and spatial by design.", c:PURPLE },
    { n:"STOmicsDB", v:"218", d:"Manually curated spatial datasets across 17 species, with cell types, spatial regions and CCI pre-computed.", c:TEAL },
    { n:"CROST", v:"182", d:"Spatial datasets comprising 1,033 sub-datasets and 48,043 tumour-related spatially variable genes.", c:TEAL },
    { n:"HEST-1k", v:"1,108", d:"Paired ST and H&E whole-slide images: 1.5M spots, >60M cells — the substrate for histology-based scale-out.", c:AMBER },
    { n:"dbGaP / EGA", v:"Controlled", d:"Raw-access tiers holding much of the clinically annotated material. Access applications are on the critical path and start in month 1.", c:MUTED },
  ];
  const cw = (W-2*M-2*0.26)/3, ch = 2.16;
  src.forEach((o,i)=>{
    const x = M + (i%3)*(cw+0.26), y = 1.62 + Math.floor(i/3)*(ch+0.26);
    card(s, { x, y, w:cw, h:ch, fill:WHITE, line:LINE, shadow:true });
    s.addShape(pres.ShapeType.roundRect, { x:x+0.22, y:y+0.24, w:0.09, h:0.42, rectRadius:0.045,
      fill:{color:o.c}, line:{width:0} });
    body(s, { x:x+0.44, y:y+0.24, w:cw-0.66, h:0.3, size:14, bold:true, color:INK, text:o.n });
    body(s, { x:x+0.44, y:y+0.54, w:cw-0.66, h:0.26, size:11, bold:true, color:o.c, text:o.v });
    body(s, { x:x+0.22, y:y+0.90, w:cw-0.44, h:1.1, size:10.5, color:BODY, text:o.d, lsm:1.16 });
  });
  footer(s, "Counts from: HTAN data-sharing report (Nat Methods 2025) · STOmicsDB (NAR 2024) · CROST (NAR 2024) · HEST-1k (NeurIPS 2024).");
}

/* --------------------------------------------- 8 live GEO survey (chart) */
{
  const s = content("Aim 1 · Evidence", "What a live GEO survey returns today");
  body(s, { x:M, y:1.44, w:6.0, h:0.4, size:12, italic:true, color:BODY,
    text:"Queried against NCBI E-utilities on 9 September 2026, not estimated." });

  const funnel = [
    { v:"67,774", l:"Human cancer GSE series", c:MUTED },
    { v:"3,392",  l:"…with single-cell RNA-seq", c:PURPLE },
    { v:"647",    l:"…with any spatial platform", c:EOSIN },
    { v:"188",    l:"…with BOTH in one series", c:TEAL },
  ];
  funnel.forEach((f,i)=>{
    const y = 1.94 + i*0.78;
    card(s, { x:M, y, w:4.5, h:0.68, fill: i===3?TINT:LIGHT, line: i===3?PURPLE:LINE });
    stat(s, { x:M+0.22, y:y+0.06, w:1.5, value:f.v, vs:20, vh:0.34, label:"", ls:1, color:f.c });
    body(s, { x:M+1.8, y:y+0.16, w:2.5, h:0.4, size:10.5, color:BODY, text:f.l, lsm:1.05 });
  });

  s.addChart(pres.ChartType.bar, [{
    name: "Paired series",
    labels: ["BRCA","GBM","PDAC","SKCM","CRC","LUAD/LUSC","HCC","PRAD","HNSCC","RCC","OV","STAD"],
    values: [26,24,20,17,10,10,10,10,6,6,6,4],
  }], {
    x:5.45, y:1.72, w:7.3, h:4.3,
    barDir:"bar", showTitle:true, title:"Paired ST + scRNA-seq series per tumour type",
    titleFontSize:12, titleColor:INK, titleFontFace:SANS,
    chartColors:[PURPLE], showLegend:false,
    showValue:true, dataLabelPosition:"outEnd", dataLabelFontSize:9, dataLabelColor:BODY,
    catAxisLabelColor:BODY, catAxisLabelFontSize:9.5, catAxisLabelFontFace:SANS,
    valAxisLabelColor:MUTED, valAxisLabelFontSize:9,
    valGridLine:{ color:"E8E2F0", size:1 }, catGridLine:{ style:"none" },
    barGapWidthPct:45,
  });

  card(s, { x:M, y:5.22, w:4.5, h:1.28, fill:INK, line:INK });
  body(s, { x:M+0.22, y:5.36, w:4.06, h:1.0, size:10.5, color:"C6BADD", lsm:1.16,
    text:"These are keyword-search upper bounds, not a curated cohort. Free-text platform mentions produce false positives (a paper citing Visium is not a Visium deposit) and false negatives (unnamed platforms). Curation is Phase 2 — see the funnel on the next slide." });
  footer(s, "NCBI GEO DataSets (db=gds), Homo sapiens, GSE entry type, 9 Sep 2026. Platform split of the paired set: Visium 54 · Xenium 27 · CosMx 15 · MERFISH 2 · Stereo-seq 1.");
  s.addNotes("The disclaimer matters: reporting 188 as if it were a curated cohort would be the exact kind of number inflation this programme is designed to avoid.");
}

/* ------------------------------------------------------ 9 curation funnel */
{
  const s = content("Aim 1 · Curation", "From search hits to an analysable cohort");
  const steps = [
    { t:"Automated classification", d:"Regex platform typing with the matched substring retained for audit, and explicit negation handling so \"not Visium\" never classifies as Visium.", c:EOSIN },
    { t:"Structure verification", d:"List each series' supplementary directory over HTTPS; confirm a count matrix and spatial coordinates exist. Record what is missing rather than dropping the series.", c:EOSIN },
    { t:"Pairing level assignment", d:"Distinguish paired (both modalities, one series), matched (same patient IDs recoverable across series) and unpaired. These are different strengths of evidence and are never collapsed.", c:PURPLE },
    { t:"Metadata reconciliation", d:"Tumour type, treatment status, stage, and patient identity harmonised to a controlled vocabulary. Treatment status is unstated in most depositions and must be recovered from the source publication.", c:PURPLE },
    { t:"Access triage", d:"Open versus dbGaP/EGA-controlled. Controlled-access applications start immediately; the open subset defines what analysis can begin in month 1.", c:TEAL },
    { t:"Manual review", d:"Every retained series read by a curator against its publication. Automated classification proposes; a human disposes.", c:TEAL },
  ];
  const cw = (W-2*M-0.28)/2, ch = 1.42;
  steps.forEach((o,i)=>{
    const x = M + (i%2)*(cw+0.28), y = 1.58 + Math.floor(i/2)*(ch+0.2);
    card(s, { x, y, w:cw, h:ch, fill: i%2 ? LIGHT : WHITE, line:LINE });
    dot(s, x+0.24, y+0.26, 0.48, o.c, i+1);
    body(s, { x:x+0.86, y:y+0.22, w:cw-1.1, h:0.3, size:13, bold:true, color:INK, text:o.t });
    body(s, { x:x+0.86, y:y+0.56, w:cw-1.1, h:0.76, size:10.5, color:BODY, text:o.d, lsm:1.15 });
  });
  card(s, { x:M, y:6.0, w:W-2*M, h:0.66, fill:TINT, line:TINT });
  body(s, { x:M+0.26, y:6.14, w:W-2*M-0.52, h:0.42, size:11.5, color:PLUM,
    text:"Realistic expectation: of ~188 candidate series, the analysable paired cohort after curation and access resolution is a substantially smaller number. That number is an output of Phase 2, not an assumption of the proposal." });
  footer(s, "");
}

/* ------------------------------------------------------------- 10 QC */
{
  const s = content("Aim 1 · Quality", "QC that adapts to tissue instead of imposing on it");
  const left = [
    { h:"Per-sample MAD thresholds", d:"Mitochondrial content, feature count and UMI depth are cut at median ± 3 MAD within each sample. A flat 10% mitochondrial filter deletes real cells in kidney, liver and heart, and retains debris elsewhere." },
    { h:"Doublet removal", d:"DoubletFinder or scDblFinder per sample, before integration. Undetected doublets are the main source of spurious NK/T and macrophage/tumour \"intermediate\" states." },
    { h:"Ambient RNA", d:"SoupX or CellBender. Ambient contamination systematically inflates marker genes of abundant cell types, which directly biases deconvolution references." },
  ];
  const right = [
    { h:"Spatial-specific QC", d:"Spots off-tissue, sections with staining artefacts, and segmentation failures on imaging platforms. Segmentation error propagates into every distance measurement downstream." },
    { h:"Coordinate normalisation", d:"All coordinates converted to microns at load. Mixing pixel coordinates from one platform with micron coordinates from another silently produces meaningless distances." },
    { h:"Audit output", d:"Every threshold and every cell count lost is written to a per-sample table that becomes a supplementary file. Silent filtering is unreproducible filtering." },
  ];
  [[left, M], [right, M + (W-2*M)/2 + 0.16]].forEach(([col, x])=>{
    const cw = (W-2*M)/2 - 0.16;
    col.forEach((o,i)=>{
      const y = 1.62 + i*1.62;
      card(s, { x, y, w:cw, h:1.44, fill:WHITE, line:LINE, shadow:true });
      s.addShape(pres.ShapeType.ellipse, { x:x+0.24, y:y+0.28, w:0.16, h:0.16,
        fill:{color: x===M ? EOSIN : TEAL}, line:{width:0} });
      body(s, { x:x+0.52, y:y+0.22, w:cw-0.76, h:0.3, size:13, bold:true, color:INK, text:o.h });
      body(s, { x:x+0.52, y:y+0.56, w:cw-0.76, h:0.8, size:10.5, color:BODY, text:o.d, lsm:1.16 });
    });
  });
  footer(s, "");
}

/* --------------------------------------------------- 11 integration trap */
{
  const s = content("Aim 1 · Harmonisation", "The confound that breaks pan-cancer integration");
  card(s, { x:M, y:1.52, w:5.9, h:1.5, fill:INK, line:INK });
  body(s, { x:M+0.3, y:1.72, w:5.3, h:0.3, size:12, bold:true, color:EOSIN,
    text:"THE PROBLEM" });
  body(s, { x:M+0.3, y:2.04, w:5.3, h:0.86, size:12, color:"CFC4E6", lsm:1.16,
    text:"In a cohort assembled from public depositions, \"batch\" and \"cancer type\" are very nearly the same variable. Each cancer type arrives from different labs, chemistries and years. Batch correction that removes this removes the biology the study is about." });

  card(s, { x:M+6.18, y:1.52, w:W-2*M-6.18, h:1.5, fill:TINT, line:PURPLE });
  body(s, { x:M+6.46, y:1.72, w:W-2*M-6.74, h:0.3, size:12, bold:true, color:PURPLE,
    text:"THE APPROACH" });
  body(s, { x:M+6.46, y:2.04, w:W-2*M-6.74, h:0.86, size:12, color:PLUM, lsm:1.16,
    text:"Integrate within cancer type, over patient. Compare across cancer types at the level of shared cell-type labels and per-sample summary statistics — replication across cohorts, not a shared latent space that may have erased the difference being tested." });

  const tools = [
    { n:"Harmony", d:"Fast linear correction; strong default within a cancer type, and the baseline that single-cell foundation models have repeatedly failed to beat in zero-shot settings.", c:PURPLE },
    { n:"scVI / scANVI", d:"Probabilistic latent space; scANVI uses partial labels. Preferred where nonlinear batch effects dominate or where uncertainty is needed downstream.", c:TEAL },
    { n:"Seurat v5 layers", d:"Layer-based objects with BPCells on-disk backing; the practical route to cohort-scale objects without exhausting node memory.", c:EOSIN },
    { n:"scArches", d:"Reference mapping rather than de novo integration — the right tool for projecting new cohorts onto a frozen reference without refitting it.", c:AMBER },
  ];
  const cw = (W-2*M-3*0.22)/4;
  tools.forEach((t,i)=>{
    const x = M + i*(cw+0.22);
    card(s, { x, y:3.28, w:cw, h:2.02, fill:LIGHT, line:LINE });
    body(s, { x:x+0.22, y:3.48, w:cw-0.44, h:0.3, size:13.5, bold:true, color:t.c, text:t.n });
    body(s, { x:x+0.22, y:3.84, w:cw-0.44, h:1.32, size:10.5, color:BODY, text:t.d, lsm:1.16 });
  });
  card(s, { x:M, y:5.52, w:W-2*M, h:1.0, fill:WHITE, line:LINE });
  body(s, { x:M+0.28, y:5.66, w:W-2*M-0.56, h:0.76, size:11.5, color:BODY, lsm:1.16,
    text:"Integration is assessed with both batch-mixing and biological-conservation metrics (scIB), reported per cancer type. An integration that maximises mixing while collapsing known cell-type structure has failed, however good the UMAP looks — and the choice of method is documented per cohort rather than applied globally." });
  footer(s, "");
}

/* ------------------------------------------- 12 annotation + foundation */
{
  const s = content("Aim 1 · Annotation", "Cell typing, and an honest read on foundation models");
  const std = [
    "Marker-based scoring against a curated pan-cancer ontology, with an explicit unassigned class",
    "Reference mapping: Azimuth, SingleR, CellTypist, scArches onto frozen references",
    "Rare and ambiguous populations reported, not absorbed — NKT, γδ T and ILC1 sit between NK and T",
    "Every label carries a confidence and the evidence that produced it",
  ];
  card(s, { x:M, y:1.58, w:6.05, h:2.52, fill:LIGHT, line:LINE });
  body(s, { x:M+0.28, y:1.78, w:5.5, h:0.3, size:13.5, bold:true, color:PURPLE,
    text:"Established methods carry the annotation" });
  bullets(s, std, { x:M+0.28, y:2.16, w:5.5, h:1.8, size:11, gap:8 });

  card(s, { x:M+6.31, y:1.58, w:W-2*M-6.31, h:2.52, fill:INK, line:INK });
  body(s, { x:M+6.59, y:1.78, w:W-2*M-6.87, h:0.3, size:13.5, bold:true, color:EOSIN,
    text:"What the benchmarks actually say" });
  body(s, { x:M+6.59, y:2.16, w:W-2*M-6.87, h:1.8, size:11, color:"CFC4E6", lsm:1.2,
    text:"A 2025 Genome Biology zero-shot evaluation found Geneformer and scGPT performed poorly relative to scVI and Harmony for cell-type clustering — in some cases worse than simply using 2,000 highly variable genes, and worse than the same architectures initialised to random weights. Both also showed limited ability to predict held-out gene expression." });

  card(s, { x:M, y:4.3, w:W-2*M, h:2.16, fill:WHITE, line:PURPLE });
  body(s, { x:M+0.3, y:4.48, w:W-2*M-0.6, h:0.3, size:13.5, bold:true, color:INK,
    text:"Where foundation models do earn their place: spatial pretraining" });
  body(s, { x:M+0.3, y:4.84, w:W-2*M-0.6, h:1.4, size:11.5, color:BODY, lsm:1.2,
    text:"Nicheformer (Nature Methods, 2025) was pretrained on SpatialCorpus-110M — over 57 million dissociated and 53 million spatially resolved cells across 73 tissues — and outperforms scGPT, Geneformer, UCE and CellPLM on spatial composition and spatial label prediction under linear probing and fine-tuning. Notably, a model trained on 1% of the spatial data beat models trained on three times as much dissociated data, indicating that spatial context is not recoverable from dissociated profiles alone.\n\nPosition adopted: foundation models are used fine-tuned, for spatially defined tasks, and always benchmarked against scVI, Harmony and an HVG baseline on our own held-out patients. A foundation model that does not beat those baselines here does not enter the pipeline." });
  footer(s, "Sources: Genome Biology 2025;26:101 (zero-shot evaluation) · Nature Methods 2025 (Nicheformer).");
}

/* ------------------------------------------------------ 13 deconvolution */
{
  const s = content("Aim 2 · Mapping", "Deconvolution: the benchmarks disagree, so run a consensus");
  card(s, { x:M, y:1.5, w:6.05, h:1.66, fill:LIGHT, line:LINE });
  body(s, { x:M+0.26, y:1.66, w:5.53, h:0.28, size:11.5, bold:true, color:PURPLE,
    text:"Nature Communications, 2023" });
  body(s, { x:M+0.26, y:1.98, w:5.53, h:1.02, size:11, color:BODY, lsm:1.16,
    text:"18 methods across 50 real and simulated datasets. CARD, cell2location and Tangram ranked best overall for cellular deconvolution." });

  card(s, { x:M+6.31, y:1.5, w:W-2*M-6.31, h:1.66, fill:LIGHT, line:LINE });
  body(s, { x:M+6.57, y:1.66, w:W-2*M-6.83, h:0.28, size:11.5, bold:true, color:EOSIN,
    text:"Bioinformatics, 2023" });
  body(s, { x:M+6.57, y:1.98, w:W-2*M-6.83, h:1.02, size:11, color:BODY, lsm:1.16,
    text:"RCTD and cell2location were the top performers on both AUPR and Jensen–Shannon divergence. A separate 2026 benchmark put cell2location, RCTD and Tangram first on brain tissue." });

  card(s, { x:M, y:3.34, w:W-2*M, h:1.16, fill:INK, line:INK });
  body(s, { x:M+0.3, y:3.5, w:W-2*M-0.6, h:0.3, size:12.5, bold:true, color:TEAL,
    text:"Only cell2location appears at the top of every list. That is the actual signal." });
  body(s, { x:M+0.3, y:3.84, w:W-2*M-0.6, h:0.52, size:11, color:"C6BADD", lsm:1.16,
    text:"Rankings shift with tissue, platform, reference quality and metric. Choosing one winner from one benchmark and applying it pan-cancer imports that benchmark's tissue bias into every result." });

  const plan = [
    { n:"Run four", d:"cell2location, RCTD, CARD and Tangram on every spot-based section — the union of the top performers across benchmarks.", c:PURPLE },
    { n:"Consensus", d:"Report the per-cell-type median weight and the inter-method concordance. Disagreement between methods is itself a reliability signal and is retained, not averaged away.", c:EOSIN },
    { n:"Both RCTD modes", d:"Full mode returns a weight for every cell type; doublet mode commits to at most two. Where they disagree, the full-mode weight is largely reference similarity rather than evidence of presence.", c:TEAL },
    { n:"Benchmark locally", d:"Simulated spot mixtures built from our own reference give exact ground-truth proportions, so reliability is measured on this cohort rather than assumed from a paper.", c:AMBER },
  ];
  const cw = (W-2*M-3*0.22)/4;
  plan.forEach((p,i)=>{
    const x = M + i*(cw+0.22);
    card(s, { x, y:4.66, w:cw, h:1.82, fill:WHITE, line:LINE, shadow:true });
    body(s, { x:x+0.22, y:4.84, w:cw-0.44, h:0.3, size:13, bold:true, color:p.c, text:p.n });
    body(s, { x:x+0.22, y:5.18, w:cw-0.44, h:1.18, size:10.5, color:BODY, text:p.d, lsm:1.16 });
  });
  footer(s, "Sources: Nat Commun 2023;14:1548 · Bioinformatics 2023;39:btac805 · comparative benchmarking 2026.");
}

/* --------------------------------------------------- 14 reliability gate */
{
  const s = content("Aim 2 · Gate", "Rare cell types are where deconvolution fails");
  body(s, { x:M, y:1.46, w:W-2*M, h:0.36, size:12.5, italic:true, color:BODY,
    text:"The weakest link in every spot-based pan-cancer study, and the one most often left unstated." });

  const facts = [
    { v:"0.5–3%", l:"NK cells as a fraction\nof a typical solid tumour", c:EOSIN },
    { v:"1–10", l:"Cells inside a single\n55 µm Visium spot", c:PURPLE },
    { v:"≈0", l:"Separation between NK and CD8\neffector reference profiles at depth", c:TEAL },
  ];
  const sw = 2.60;
  facts.forEach((f,i)=>{
    const x = M + i*(sw+0.3);
    card(s, { x, y:1.96, w:sw, h:1.5, fill:INK, line:INK });
    stat(s, { x:x+0.24, y:2.16, w:sw-0.48, value:f.v, vs:30, vh:0.56, label:f.l, ls:10,
      color:f.c, lh:0.7 });
  });

  card(s, { x:M+3*(sw+0.3), y:1.96, w:W-M-(M+3*(sw+0.3)), h:1.5, fill:TINT, line:PURPLE });
  body(s, { x:M+3*(sw+0.3)+0.26, y:2.14, w:W-M-(M+3*(sw+0.3))-0.52, h:1.14, size:11,
    color:PLUM, lsm:1.16,
    text:"An \"NK cell\" on Visium is a similarity score against a profile whose nearest competitor is CD8 T. The resulting map often just reproduces the CD8 map." });

  card(s, { x:M, y:3.66, w:W-2*M, h:2.4, fill:WHITE, line:LINE, shadow:true });
  body(s, { x:M+0.3, y:3.84, w:W-2*M-0.6, h:0.3, size:13.5, bold:true, color:INK,
    text:"The gate, applied before any cell type carries a conclusion" });
  const gate = [
    "Simulated spot mixtures are built by summing counts from a known number of reference cells, giving exact ground-truth proportions",
    "Per-cell-type Pearson r and RMSE against that truth become an explicit reliability flag written to the results table",
    "Cell types below threshold are marked deconvolution-limited: their weights may describe niche context, never abundance or a target claim",
    "Rare-population claims are restricted to Xenium, CosMx and MERFISH, where cells are segmented rather than inferred",
    "Every figure states which platform class supports it, so a reader can see at a glance whether a claim rests on measurement or on inference",
  ];
  bullets(s, gate, { x:M+0.3, y:4.22, w:W-2*M-0.6, h:1.7, size:11, gap:6 });
  footer(s, "");
  s.addNotes("If one slide from this deck survives into the methods section, it should be this one.");
}

/* -------------------------------------------------- 15 spatial domains */
{
  const s = content("Aim 2 · Domains", "Spatial domain detection, chosen on benchmark evidence");
  s.addChart(pres.ChartType.bar, [{
    name:"Mean ARI (10x Visium)",
    labels:["GraphST","STAGATE","CCST"],
    values:[0.552, 0.515, 0.481],
  }], {
    x:M, y:1.62, w:6.0, h:2.7,
    barDir:"col", showTitle:true, title:"Visium benchmark, mean adjusted Rand index",
    titleFontSize:11.5, titleColor:INK, titleFontFace:SANS,
    chartColors:[PURPLE], showLegend:false, showValue:true, dataLabelPosition:"outEnd",
    dataLabelFontSize:9.5, dataLabelColor:BODY, dataLabelFormatCode:"0.000",
    catAxisLabelColor:BODY, catAxisLabelFontSize:9.5,
    valAxisLabelColor:MUTED, valAxisLabelFontSize:9, valAxisMaxVal:0.65,
    valGridLine:{ color:"E8E2F0", size:1 }, catGridLine:{ style:"none" },
    barGapWidthPct:50,
  });

  card(s, { x:M+6.26, y:1.62, w:W-2*M-6.26, h:2.7, fill:LIGHT, line:LINE });
  body(s, { x:M+6.52, y:1.82, w:W-2*M-6.78, h:0.3, size:12.5, bold:true, color:PURPLE,
    text:"Nucleic Acids Research, 2025" });
  body(s, { x:M+6.52, y:2.16, w:W-2*M-6.78, h:2.0, size:11, color:BODY, lsm:1.2,
    text:"Across a systematic benchmark, GraphST ranked first overall, closely followed by BayesSpace, SpaGCN and STAGATE. Platform matters: on legacy ST data BayesSpace led decisively (ARI 0.642), ahead of Leiden (0.562).\n\nNo method wins everywhere. Method selection is therefore made per platform, and the choice is recorded with the result." });

  const meth = [
    { n:"GraphST", d:"Graph contrastive self-supervised learning on HVGs. Default for Visium.", c:PURPLE },
    { n:"BayesSpace", d:"Bayesian spatial prior with sub-spot resolution. Default for legacy ST and low-spot-count sections.", c:EOSIN },
    { n:"STAGATE", d:"Adaptive graph attention auto-encoder; strong and stable second opinion.", c:TEAL },
    { n:"SpaGCN", d:"Graph convolution integrating histology; useful where the image carries real signal.", c:AMBER },
  ];
  const cw = (W-2*M-3*0.22)/4;
  meth.forEach((m,i)=>{
    const x = M + i*(cw+0.22);
    card(s, { x, y:4.52, w:cw, h:1.34, fill:WHITE, line:LINE });
    body(s, { x:x+0.2, y:4.70, w:cw-0.4, h:0.28, size:13, bold:true, color:m.c, text:m.n });
    body(s, { x:x+0.2, y:5.02, w:cw-0.4, h:0.76, size:10.5, color:BODY, text:m.d, lsm:1.16 });
  });
  card(s, { x:M, y:6.04, w:W-2*M, h:0.6, fill:TINT, line:TINT });
  body(s, { x:M+0.26, y:6.16, w:W-2*M-0.52, h:0.38, size:11, color:PLUM,
    text:"Domains are called with at least two methods per section; regions where they disagree are reported as boundary-uncertain rather than assigned." });
  footer(s, "Source: Nucleic Acids Research 2025;53(7):gkaf303. Bars show mean ARI on 10x Visium for the three methods with values reported; BayesSpace and SpaGCN ranked highly overall, with BayesSpace leading on legacy ST (ARI 0.642).");
}

/* ------------------------------------------------------------ 16 niches */
{
  const s = content("Aim 2 · Niches", "Niche calling at single-cell resolution");
  const rows = [
    { n:"NicheCompass", d:"Graph deep learning that models cellular communication to learn interpretable embeddings encoding signalling events. Uniquely recovered spatially contiguous niches and led on spatial consistency and niche coherence against BANKSY, GraphST and CellCharter.", tag:"Primary", c:PURPLE },
    { n:"BANKSY", d:"Augments each cell's profile with its neighbours', with a tunable weight on the microenvironment. Fast, well-understood, scales past 3 million cells.", tag:"Scale", c:TEAL },
    { n:"CellCharter", d:"Neighbourhood aggregation with explicit cluster-number selection; strong on imaging-based platforms.", tag:"Second opinion", c:EOSIN },
    { n:"UTAG / scNiche", d:"Linear neighbour weighting and multi-view learning respectively; both scale to >3M cells, useful as fast baselines on the largest sections.", tag:"Baseline", c:AMBER },
  ];
  rows.forEach((r,i)=>{
    const y = 1.6 + i*1.18;
    card(s, { x:M, y, w:W-2*M, h:1.06, fill: i%2 ? LIGHT : WHITE, line:LINE });
    body(s, { x:M+0.28, y:y+0.2, w:2.3, h:0.32, size:14, bold:true, color:INK, text:r.n });
    s.addShape(pres.ShapeType.roundRect, { x:M+0.28, y:y+0.58, w:1.34, h:0.26, rectRadius:0.13,
      fill:{color:r.c}, line:{width:0} });
    body(s, { x:M+0.28, y:y+0.58, w:1.34, h:0.26, size:8.5, bold:true, color:WHITE,
      text:r.tag, align:"center", valign:"middle" });
    body(s, { x:M+2.78, y:y+0.2, w:W-2*M-3.06, h:0.72, size:11, color:BODY, text:r.d, lsm:1.16 });
  });
  body(s, { x:M, y:6.30, w:W-2*M, h:0.42, size:11, color:PLUM, italic:true,
    text:"Niches are called per section, matched across sections by composition, then tested for recurrence across patients and cancer types — a niche seen in one tumour is an observation; a niche recurring across cohorts is a finding." });
  footer(s, "Sources: NicheCompass (Nature Genetics 2025) · niche-identification benchmark 2026.");
}

/* ------------------------------------------------- 17 malignant / CNV */
{
  const s = content("Aim 2 · Tumour cells", "Identifying the malignant compartment");
  body(s, { x:M, y:1.42, w:W-2*M, h:0.52, size:12.5, italic:true, color:BODY,
    text:"Every tumour-intrinsic claim depends on correctly separating malignant from normal cells. Marker genes cannot do this; inferred copy number can." });

  const tools = [
    { n:"Numbat", d:"Uses expression and B-allele frequency together. Best overall performance for distinguishing tumour from normal in benchmark, and the default wherever allele information is recoverable.", note:"Highest CPU cost of the tools tested", c:PURPLE },
    { n:"CopyKAT", d:"Expression-only. Best of the expression-only tools, and notably robust to shallow sequencing — still separated tumour from normal at a median of 1,000 UMIs per cell.", note:"Default when only a count matrix exists", c:TEAL },
    { n:"inferCNV", d:"The field standard and the most widely reported. Retained for comparability with published atlases rather than as the primary caller.", note:"Moderate runtime", c:EOSIN },
    { n:"SCEVAN / CaSpER", d:"SCEVAN is among the fastest; CaSpER also uses allele information. Both run as cross-checks on ambiguous samples.", note:"Cross-check tier", c:AMBER },
  ];
  const cw = (W-2*M-3*0.22)/4;
  tools.forEach((t,i)=>{
    const x = M + i*(cw+0.22);
    card(s, { x, y:1.94, w:cw, h:2.52, fill:WHITE, line:LINE, shadow:true });
    body(s, { x:x+0.22, y:2.12, w:cw-0.44, h:0.3, size:14, bold:true, color:t.c, text:t.n });
    body(s, { x:x+0.22, y:2.48, w:cw-0.44, h:1.42, size:10.5, color:BODY, text:t.d, lsm:1.16 });
    body(s, { x:x+0.22, y:3.98, w:cw-0.44, h:0.4, size:9.5, color:MUTED, italic:true, text:t.note, lsm:1.1 });
  });

  card(s, { x:M, y:4.68, w:6.05, h:1.82, fill:LIGHT, line:LINE });
  body(s, { x:M+0.28, y:4.86, w:5.49, h:0.3, size:13, bold:true, color:INK,
    text:"Beyond the binary call" });
  bullets(s, [
    "Subclone reconstruction from CNV profiles, then mapped back onto tissue coordinates",
    "Malignant transcriptional programs derived per patient by cNMF, then matched across patients into recurrent meta-programs",
    "Clone-to-niche assignment: which subclone occupies which microenvironment",
  ], { x:M+0.28, y:5.24, w:5.49, h:1.16, size:10.5, gap:6 });

  card(s, { x:M+6.31, y:4.68, w:W-2*M-6.31, h:1.82, fill:INK, line:INK });
  body(s, { x:M+6.59, y:4.86, w:W-2*M-6.87, h:0.3, size:13, bold:true, color:EOSIN,
    text:"The caveat that must travel with these results" });
  body(s, { x:M+6.59, y:5.24, w:W-2*M-6.87, h:1.16, size:11, color:"C6BADD", lsm:1.18,
    text:"Classification accuracy degrades with sequencing depth for every tool, sharply for Numbat. Depth is therefore reported alongside every malignant-fraction estimate, and cross-cohort comparisons of tumour purity are made only between samples of comparable depth." });
  footer(s, "Sources: Brief Bioinform 2025;26(2):bbaf076 · Nat Commun 2025 (scRNA-seq CNV caller benchmark).");
}

/* --------------------------------------------------------- 18 CCI */
{
  const s = content("Aim 2 · Mechanism", "Cell–cell interaction, constrained by physical distance");
  card(s, { x:M, y:1.5, w:W-2*M, h:1.04, fill:INK, line:INK });
  body(s, { x:M+0.3, y:1.64, w:W-2*M-0.6, h:0.3, size:12.5, bold:true, color:EOSIN,
    text:"The problem with ligand–receptor inference from dissociated data" });
  body(s, { x:M+0.3, y:1.96, w:W-2*M-0.6, h:0.46, size:11, color:"C6BADD", lsm:1.14,
    text:"The LIANA comparison of seven methods and sixteen resources found low overlap among top-ranked interactions — the tools disagree because their scoring strategies differ, not because one is right. Spatial data resolves this by adding the constraint that interacting cells must actually be adjacent." });

  const tools = [
    { n:"CellPhoneDB v5", d:"Permutation-based specificity scoring; robust to noisy input and to errors in the underlying resource. Widely reported, so it anchors comparability with prior atlases.", c:PURPLE },
    { n:"CellChat v2", d:"Pathway-level aggregation with spatial visualisation directly on tissue architecture; strong balance of spatial consistency and agreement with commonly found interactions.", c:PURPLE },
    { n:"NicheNet", d:"Links ligands to downstream target-gene expression in the receiver. Ranked the most consistent with spatial information in benchmark — the tool that comes closest to testing consequence, not just co-expression.", c:TEAL },
    { n:"COMMOT", d:"Collective optimal transport over spatial coordinates, accounting for competition between ligands and receptors; assigns signalling activity to individual cell positions.", c:TEAL },
    { n:"LIANA+", d:"Consensus framework across methods and resources; used to produce a rank-aggregated call rather than trusting any single scorer.", c:EOSIN },
    { n:"SpaCCI", d:"Spatially informed scoring that reached the highest mean normalised F1 in a recent spatial CCC benchmark (92.7%), ahead of CellPhoneDB v3 (77.4%) and CellChat v2 (69.6%).", c:AMBER },
  ];
  const cw = (W-2*M-2*0.24)/3, ch = 1.5;
  tools.forEach((t,i)=>{
    const x = M + (i%3)*(cw+0.24), y = 2.72 + Math.floor(i/3)*(ch+0.22);
    card(s, { x, y, w:cw, h:ch, fill: LIGHT, line:LINE });
    body(s, { x:x+0.22, y:y+0.18, w:cw-0.44, h:0.28, size:13, bold:true, color:t.c, text:t.n });
    body(s, { x:x+0.22, y:y+0.5, w:cw-0.44, h:0.9, size:10, color:BODY, text:t.d, lsm:1.14 });
  });
  card(s, { x:M, y:6.02, w:W-2*M, h:0.66, fill:TINT, line:TINT });
  body(s, { x:M+0.26, y:6.16, w:W-2*M-0.52, h:0.42, size:11.5, color:PLUM,
    text:"Interpretation discipline: a significant ligand–receptor score is co-expression in adjacent neighbourhoods. It is a hypothesis about signalling, never evidence of it. Language in every output distinguishes co-localisation from communication." });
  footer(s, "Sources: Nat Commun 2022;13:3224 (LIANA) · spatial CCC benchmark 2025 · COMMOT (Nat Methods 2023).");
}

/* --------------------------------------------------------- 19 GRN */
{
  const s = content("Aim 2 · Mechanism", "Transcription-factor networks, spatially resolved");
  const steps = [
    { n:"1", t:"Co-expression modules", d:"GRNBoost2 / GENIE3 link candidate transcription factors to targets across the cohort.", c:EOSIN },
    { n:"2", t:"Motif pruning", d:"cisTarget retains only targets whose promoters and enhancers carry the TF's motif — the step that separates regulons from correlation.", c:PURPLE },
    { n:"3", t:"Enhancer evidence", d:"SCENIC+ predicts genomic enhancers, links them to upstream TFs and target genes, where paired scATAC exists.", c:TEAL },
    { n:"4", t:"Spatial projection", d:"Regulon activity scored per cell or spot, then mapped onto tissue coordinates and tested for niche enrichment.", c:AMBER },
  ];
  const cw = (W-2*M-3*0.2)/4;
  steps.forEach((o,i)=>{
    const x = M + i*(cw+0.2);
    card(s, { x, y:1.6, w:cw, h:2.2, fill:WHITE, line:LINE, shadow:true });
    dot(s, x+cw/2-0.24, 1.82, 0.48, o.c, o.n);
    body(s, { x:x+0.18, y:2.42, w:cw-0.36, h:0.34, size:12.5, bold:true, color:INK,
      text:o.t, align:"center", lsm:1.05 });
    body(s, { x:x+0.18, y:2.84, w:cw-0.36, h:0.86, size:10, color:BODY, text:o.d,
      align:"center", lsm:1.16 });
    if(i<3) s.addShape(pres.ShapeType.rightArrow, { x:x+cw+0.02, y:1.94, w:0.16, h:0.24,
      fill:{color:LINE}, line:{width:0} });
  });

  card(s, { x:M, y:4.0, w:6.05, h:2.4, fill:INK, line:INK });
  body(s, { x:M+0.28, y:4.18, w:5.49, h:0.3, size:13, bold:true, color:EOSIN,
    text:"What the GRN benchmark says" });
  body(s, { x:M+0.28, y:4.54, w:5.49, h:1.74, size:11, color:"C6BADD", lsm:1.2,
    text:"BEELINE (Nature Methods, 2020) evaluated GRN inference on synthetic and experimental single-cell data and found heterogeneous performance, with moderate area under the precision–recall curve and moderate early precision across algorithms. Methods that do not require pseudotime-ordered cells were generally more accurate.\n\nRegulons are therefore treated as prioritised hypotheses for perturbation, not as a validated network." });

  card(s, { x:M+6.31, y:4.0, w:W-2*M-6.31, h:2.4, fill:LIGHT, line:LINE });
  body(s, { x:M+6.59, y:4.18, w:W-2*M-6.87, h:0.3, size:13, bold:true, color:PURPLE,
    text:"How this cohort adds evidence a single study cannot" });
  bullets(s, [
    "A regulon is called conserved only if it recurs across cancer types, tested at patient level",
    "Regulon activity is required to co-localise with the niche it is claimed to drive",
    "Where paired scATAC exists, enhancer support is required rather than motif presence alone",
    "Predicted TF–target edges are cross-referenced against perturbation resources before any is nominated as a target",
  ], { x:M+6.59, y:4.56, w:W-2*M-6.87, h:1.72, size:10.5, gap:7 });
  footer(s, "Sources: Nature Methods 2020;17:147 (BEELINE) · SCENIC+ (Nature Methods 2023).");
}

/* ------------------------------------------------------- 20 co-expression */
{
  const s = content("Aim 2 · Mechanism", "Co-expression modules across cells and space");
  card(s, { x:M, y:1.56, w:6.05, h:2.34, fill:LIGHT, line:LINE });
  body(s, { x:M+0.28, y:1.74, w:5.49, h:0.3, size:13.5, bold:true, color:PURPLE,
    text:"hdWGCNA as the backbone" });
  body(s, { x:M+0.28, y:2.1, w:5.49, h:1.66, size:11, color:BODY, lsm:1.2,
    text:"Weighted co-expression network analysis adapted to high-dimensional data. Highly similar cells are collapsed into metacells to counter single-cell sparsity while retaining heterogeneity, and separate networks are constructed per cell population — the modular design that makes context-specific networks possible across a cellular and spatial hierarchy. Seurat-compatible and demonstrated on datasets approaching one million cells." });

  card(s, { x:M+6.31, y:1.56, w:W-2*M-6.31, h:2.34, fill:WHITE, line:LINE, shadow:true });
  body(s, { x:M+6.59, y:1.74, w:W-2*M-6.87, h:0.3, size:13.5, bold:true, color:TEAL,
    text:"Complementary decompositions" });
  bullets(s, [
    "cNMF for malignant transcriptional programs, matched into recurrent meta-programs across patients",
    "MOFA+ for variance decomposition across modalities and cancer types",
    "Non-negative matrix factorisation on spatial spots to recover regional programs directly",
    "Spatially variable gene detection (SPARK-X, nnSVG, SpatialDE) to seed spatially anchored modules",
  ], { x:M+6.59, y:2.12, w:W-2*M-6.87, h:1.64, size:10.5, gap:6 });

  card(s, { x:M, y:4.1, w:W-2*M, h:2.32, fill:WHITE, line:PURPLE });
  body(s, { x:M+0.3, y:4.28, w:W-2*M-0.6, h:0.3, size:13.5, bold:true, color:INK,
    text:"Turning modules into something testable" });
  const cw2 = (W-2*M-0.6-2*0.3)/3;
  [
    { h:"Anchor in space", d:"Each module's eigengene is scored per spot and tested for enrichment in specific niches, so a module becomes a property of a location rather than of a cluster." },
    { h:"Test for conservation", d:"Module preservation statistics across cancer types separate genuinely pan-cancer programs from ones driven by a single cohort." },
    { h:"Connect to the network", d:"Module hub genes are intersected with regulon targets and with ligand–receptor partners, so modules, TF networks and interactions describe one mechanism rather than three parallel results." },
  ].forEach((o,i)=>{
    const x = M+0.3 + i*(cw2+0.3);
    body(s, { x, y:4.68, w:cw2, h:0.3, size:12, bold:true, color:[EOSIN,PURPLE,TEAL][i], text:o.h });
    body(s, { x, y:5.0, w:cw2, h:1.24, size:10.5, color:BODY, text:o.d, lsm:1.18 });
  });
  footer(s, "Source: hdWGCNA, Cell Reports Methods 2023;3:100498.");
}

/* ------------------------------------------------------ 21 trajectories */
{
  const s = content("Aim 2 · Dynamics", "Trajectories and fate, with velocity's limits stated");
  const items = [
    { n:"CellRank 2", d:"Unified fate mapping over multiple data views — pseudotime, RNA velocity, real time points, similarity — scaling to millions of cells. The primary framework, because it does not require velocity to work.", c:PURPLE, tag:"Primary" },
    { n:"Monocle 3 / Slingshot / PAGA", d:"Graph and lineage-based pseudotime without a kinetic model. Robust, interpretable, and the reference against which velocity-derived directions are checked.", c:TEAL, tag:"Baseline" },
    { n:"scVelo / UniTVelo / veloVI", d:"Dynamical RNA velocity. Used where its assumptions can be argued, and never as the sole evidence for a direction of change.", c:AMBER, tag:"Conditional" },
    { n:"Spatially informed velocity", d:"Emerging methods that condition dynamics on tissue position, testing whether fate transitions align with spatial gradients rather than assuming they do.", c:EOSIN, tag:"Exploratory" },
  ];
  items.forEach((o,i)=>{
    const y = 1.58 + i*1.06;
    card(s, { x:M, y, w:7.35, h:0.94, fill: i%2 ? LIGHT : WHITE, line:LINE });
    body(s, { x:M+0.24, y:y+0.14, w:2.6, h:0.36, size:12.5, bold:true, color:INK, text:o.n, lsm:1.0 });
    s.addShape(pres.ShapeType.roundRect, { x:M+0.24, y:y+0.56, w:1.16, h:0.24, rectRadius:0.12,
      fill:{color:o.c}, line:{width:0} });
    body(s, { x:M+0.24, y:y+0.56, w:1.16, h:0.24, size:8, bold:true, color:WHITE, text:o.tag,
      align:"center", valign:"middle" });
    body(s, { x:M+3.02, y:y+0.14, w:4.5, h:0.68, size:10, color:BODY, text:o.d, lsm:1.14 });
  });

  card(s, { x:M+7.63, y:1.58, w:W-2*M-7.63, h:4.24, fill:INK, line:INK });
  body(s, { x:M+7.91, y:1.76, w:W-2*M-8.19, h:0.56, size:13, bold:true, color:EOSIN,
    text:"Why velocity is conditional,\nnot default" });
  body(s, { x:M+7.91, y:2.34, w:W-2*M-8.19, h:3.32, size:10.5, color:"C6BADD", lsm:1.22,
    text:"Two assumptions underpinning RNA velocity are readily violated in practice: a common splicing rate across genes, and that the system's equilibria are actually observed during the experiment. Where they fail, inference yields incorrect results.\n\nVelocity methods also fit genes independently, so they do not preserve gene–gene coherence: the same cell can be placed at the start of one gene's trajectory and the end of another's. And in systems dominated by post-transcriptional regulation, splicing kinetics need not reflect cellular dynamics at all.\n\nTumour tissue is exactly such a system. Directions of change are therefore corroborated by at least two independent views before being reported." });
  footer(s, "Sources: CellRank 2 (Nature Methods 2024) · PLOS Comput Biol 2022;18:e1010492 (RNA velocity unraveled).");
}

/* --------------------------------------------- 22 alignment + histology */
{
  const s = content("Aim 2 · Scale", "Multi-section alignment and histology-based extension");
  card(s, { x:M, y:1.56, w:6.05, h:2.48, fill:LIGHT, line:LINE });
  body(s, { x:M+0.28, y:1.74, w:5.49, h:0.3, size:13.5, bold:true, color:PURPLE,
    text:"Aligning sections into a common frame" });
  bullets(s, [
    "PASTE / PASTE2 — partial fused Gromov–Wasserstein optimal transport for slices that only partly overlap, with the overlap fraction estimated rather than assumed",
    "STalign — diffeomorphic metric mapping, handling local non-linear distortion across sections, samples and technologies",
    "SLAT — graph-based alignment, notable for working across distinct technologies and modalities",
    "In benchmark, SPACEL and PASTE gave the highest layer-wise alignment accuracy; STAligner led among embedding-based approaches",
  ], { x:M+0.28, y:2.12, w:5.49, h:1.8, size:10.5, gap:6 });

  card(s, { x:M+6.31, y:1.56, w:W-2*M-6.31, h:2.48, fill:WHITE, line:LINE, shadow:true });
  body(s, { x:M+6.59, y:1.74, w:W-2*M-6.87, h:0.3, size:13.5, bold:true, color:TEAL,
    text:"Extending reach with histology" });
  body(s, { x:M+6.59, y:2.12, w:W-2*M-6.87, h:1.8, size:10.5, color:BODY, lsm:1.2,
    text:"H&E slides vastly outnumber spatial assays. HEST-1k pairs 1,108 ST samples with whole-slide images — 1.5 million spots and over 60 million cells — and its benchmark evaluates 11 pathology foundation models on gene-expression prediction across nine organs and eight cancer types.\n\nModels such as UNI, Virchow, CONCH and Prov-GigaPath, and generative approaches like STPath and OmiCLIP, make it feasible to project spatial programs onto archival cohorts with survival data that no one will ever sequence spatially." });

  card(s, { x:M, y:4.24, w:W-2*M, h:2.2, fill:INK, line:INK });
  body(s, { x:M+0.3, y:4.42, w:W-2*M-0.6, h:0.3, size:13, bold:true, color:EOSIN,
    text:"The discipline that keeps this from becoming circular" });
  const cw3 = (W-2*M-0.6-2*0.3)/3;
  [
    { h:"Predicted ≠ measured", d:"Expression predicted from morphology is an imputation. It may generate hypotheses and extend cohorts; it may never serve as the measurement that validates the model that produced it." },
    { h:"Held-out patients only", d:"Prediction performance is reported on patients whose slides never entered training, per organ and per cancer type — pooled metrics hide the organs where the model fails." },
    { h:"Report the failures", d:"HEST-Benchmark shows foundation-model performance varies widely by task and organ. Organs where prediction is poor are named in the results, not omitted." },
  ].forEach((o,i)=>{
    const x = M+0.3 + i*(cw3+0.3);
    body(s, { x, y:4.8, w:cw3, h:0.28, size:12, bold:true, color:[EOSIN,PURPLE,TEAL][i], text:o.h });
    body(s, { x, y:5.12, w:cw3, h:1.2, size:10.5, color:"C6BADD", text:o.d, lsm:1.18 });
  });
  footer(s, "Sources: Genome Biology 2024;25:212 (alignment benchmark) · HEST-1k (NeurIPS 2024) · SLAT (Nat Commun 2023).");
}

/* -------------------------------------------------- 23 predictive models */
{
  const s = content("Aim 3 · Prediction", "Outcome models, and the leakage that invalidates them");
  card(s, { x:M, y:1.52, w:6.05, h:2.6, fill:LIGHT, line:LINE });
  body(s, { x:M+0.28, y:1.7, w:5.49, h:0.3, size:13.5, bold:true, color:PURPLE,
    text:"Models to be fitted" });
  bullets(s, [
    "Graph neural networks over the tissue graph, predicting response and survival from niche composition and spatial arrangement",
    "Multiple-instance learning over sections, so a patient rather than a spot is the unit the model predicts",
    "Cox and competing-risks models on niche abundance, adjusted for stage, age and treatment",
    "Interpretability by attention attribution and SHAP, validated against known biology before being reported as discovery",
  ], { x:M+0.28, y:2.08, w:5.49, h:1.9, size:10.5, gap:6 });

  card(s, { x:M+6.31, y:1.52, w:W-2*M-6.31, h:2.6, fill:INK, line:INK });
  body(s, { x:M+6.59, y:1.7, w:W-2*M-6.87, h:0.3, size:13.5, bold:true, color:EOSIN,
    text:"Why spatial ML scores are so often meaningless" });
  body(s, { x:M+6.59, y:2.08, w:W-2*M-6.87, h:1.9, size:11, color:"C6BADD", lsm:1.2,
    text:"Two spots 100 µm apart share microenvironment, patient, batch and — on spot platforms — sometimes the same cells. Put one in training and one in test and the model memorises rather than generalises. An AUC of 0.95 obtained this way means nothing, and it is the default outcome of a random split." });

  card(s, { x:M, y:4.32, w:W-2*M, h:2.12, fill:WHITE, line:PURPLE });
  body(s, { x:M+0.3, y:4.5, w:W-2*M-0.6, h:0.3, size:13.5, bold:true, color:INK,
    text:"The split protocol, asserted in code inside the cross-validation loop" });
  const cw4 = (W-2*M-0.6-3*0.26)/4;
  [
    { n:"01", h:"Patient-level folds", d:"Whole patients move between folds together; no patient appears on both sides." },
    { n:"02", h:"Spatial buffer", d:"Where a section is deliberately split by region, training points within several spot pitches of any test point are removed." },
    { n:"03", h:"Fit inside the fold", d:"HVG selection, scaling and integration are fitted on training data only. Fitting them cohort-wide leaks the test set." },
    { n:"04", h:"External cohort", d:"Final performance reported on a cohort held out entirely, never touched during development." },
  ].forEach((o,i)=>{
    const x = M+0.3 + i*(cw4+0.26);
    dot(s, x, 4.9, 0.42, [EOSIN,PURPLE,TEAL,AMBER][i], o.n);
    body(s, { x:x+0.54, y:4.94, w:cw4-0.54, h:0.3, size:11.5, bold:true, color:INK, text:o.h, lsm:1.0 });
    body(s, { x, y:5.44, w:cw4, h:0.92, size:10, color:BODY, text:o.d, lsm:1.16 });
  });
  footer(s, "");
}

/* -------------------------------------------------- 24 statistics gates */
{
  const s = pres.addSlide();
  s.background = { color: INK };
  titleBlock(s, "Rigour", "The statistical guardrails", { dark:true });
  dotField(s, 10.9, 0.42, 2.1, 1.1, 38, [PURPLE, EOSIN, TEAL], 0.07, 45);
  body(s, { x:M, y:1.40, w:11.4, h:0.5, size:12.5, color:"A99CC6", italic:true,
    text:"Each is enforced in code, because every one is routinely violated in published spatial work." });

  const g = [
    { h:"The replication unit is the patient", d:"Compute per section, pool sections to patients, meta-analyse across patients with a random-effects model. Report I² and the leave-one-cohort-out range. An effect that collapses when one cohort is dropped is a single-cohort finding reported with a pan-cancer n.", c:EOSIN },
    { h:"Spatial nulls preserve autocorrelation", d:"Torus-shift and block-bootstrap nulls translate the focal pattern rigidly, keeping its own clustering intact and randomising only its registration against the anchor. Tissue masks are mandatory — without them, shifted points land off-tissue and the null manufactures attraction.", c:PURPLE },
    { h:"One multiple-testing family", d:"FDR is corrected across everything tested — all cell-type pairs, all anchors, all cancer types — with the family size stated in the figure legend. Correcting each cancer type separately and reporting the union has not corrected at all.", c:TEAL },
    { h:"Effect sizes, not just p-values", d:"Every association is reported as an effect size with a confidence interval and the number of patients contributing. With cohorts this size, significance is cheap and magnitude is the informative quantity.", c:AMBER },
  ];
  const cw = (W-2*M-0.28)/2;
  g.forEach((o,i)=>{
    const x = M + (i%2)*(cw+0.28), y = 1.96 + Math.floor(i/2)*1.62;
    card(s, { x, y, w:cw, h:1.5, fill:PLUM, line:"51357F" });
    s.addShape(pres.ShapeType.ellipse, { x:x+0.24, y:y+0.28, w:0.17, h:0.17,
      fill:{color:o.c}, line:{width:0} });
    body(s, { x:x+0.54, y:y+0.2, w:cw-0.8, h:0.3, size:13, bold:true, color:WHITE, text:o.h });
    body(s, { x:x+0.54, y:y+0.54, w:cw-0.8, h:0.86, size:10.5, color:"C6BADD", text:o.d, lsm:1.16 });
  });

  card(s, { x:M, y:5.24, w:W-2*M, h:1.22, fill:"2E1C4F", line:EOSIN });
  body(s, { x:M+0.3, y:5.4, w:5.4, h:0.3, size:12, bold:true, color:EOSIN,
    text:"MEASURED, NOT ASSERTED" });
  body(s, { x:M+0.3, y:5.72, w:W-2*M-0.6, h:0.6, size:11.5, color:"D3C8E8", lsm:1.16,
    text:"On simulated tissue containing two independently clustered populations with no true association, a naive label shuffle returned a significant result in 24 of 24 replicates at α = 0.05. The torus-shift null returned 2 of 24. The choice of null is not a technicality — it is the difference between a finding and an artefact." });
  footer(s, "");
}

/* ------------------------------------------------------- 25 translation */
{
  const s = content("Aim 3 · Translation", "From spatial niche to prioritised target");
  const flow = [
    { n:"1", t:"Signature", d:"Each recurrent niche is reduced to a compact expression signature that can be scored in bulk data.", c:EOSIN },
    { n:"2", t:"Deconvolve bulk", d:"Score the signature across TCGA, ICGC and immunotherapy trial cohorts using CIBERSORTx, BayesPrism or MCP-counter.", c:EOSIN },
    { n:"3", t:"Associate", d:"Survival and response models adjusted for stage, age, treatment and tumour purity — never a univariate hazard ratio alone.", c:PURPLE },
    { n:"4", t:"Prioritise", d:"Intersect niche-driving ligands, receptors, regulon TFs and module hubs with druggability and dependency resources.", c:TEAL },
    { n:"5", t:"Falsify", d:"Every nominated target ships with the experiment that would refute it.", c:AMBER },
  ];
  const cw = (W-2*M-4*0.18)/5;
  flow.forEach((f,i)=>{
    const x = M + i*(cw+0.18);
    card(s, { x, y:1.6, w:cw, h:2.34, fill:WHITE, line:LINE, shadow:true });
    dot(s, x+cw/2-0.24, 1.82, 0.48, f.c, f.n);
    body(s, { x:x+0.16, y:2.42, w:cw-0.32, h:0.32, size:12.5, bold:true, color:INK,
      text:f.t, align:"center" });
    body(s, { x:x+0.16, y:2.8, w:cw-0.32, h:1.02, size:10, color:BODY, text:f.d,
      align:"center", lsm:1.16 });
    if(i<4) s.addShape(pres.ShapeType.rightArrow, { x:x+cw+0.01, y:1.94, w:0.16, h:0.24,
      fill:{color:LINE}, line:{width:0} });
  });

  card(s, { x:M, y:4.14, w:6.05, h:2.3, fill:LIGHT, line:LINE });
  body(s, { x:M+0.28, y:4.32, w:5.49, h:0.3, size:13, bold:true, color:PURPLE,
    text:"Resources for prioritisation" });
  bullets(s, [
    "DepMap dependency scores — is the target essential in the right lineage?",
    "Open Targets and DGIdb — existing chemical matter and clinical precedent",
    "CRISPR screens in matched models — genetic evidence independent of expression",
    "Structural feasibility for surface receptors and secreted ligands",
  ], { x:M+0.28, y:4.7, w:5.49, h:1.6, size:10.5, gap:7 });

  card(s, { x:M+6.31, y:4.14, w:W-2*M-6.31, h:2.3, fill:INK, line:INK });
  body(s, { x:M+6.59, y:4.32, w:W-2*M-6.87, h:0.3, size:13, bold:true, color:EOSIN,
    text:"The evidence ceiling, stated plainly" });
  body(s, { x:M+6.59, y:4.7, w:W-2*M-6.87, h:1.6, size:11, color:"C6BADD", lsm:1.2,
    text:"Everything upstream of the validation experiments is observational. A hazard ratio derived from bulk deconvolution of a different cohort is a separate and weaker line of evidence than the spatial observation it is invoked to support — not a validation of it.\n\nVerbs are matched to evidence throughout: associated with for observation, drives or mediates only where a perturbation exists." });
  footer(s, "");
}

/* -------------------------------------------------------- 26 validation */
{
  const s = content("Aim 3 · Validation", "The experiments that make it a paper");
  const tiers = [
    { tag:"TIER 1 — Confirm the observation", c:EOSIN, items:[
      "Multiplexed immunofluorescence or imaging mass cytometry on an independent tissue microarray, with the niche's defining markers",
      "Spatial proteomics (CODEX, IMC) to confirm the niche exists at protein level, not only in transcript",
      "Independent spatial cohort at single-cell resolution, analysed with the frozen pipeline and no parameter refitting",
    ]},
    { tag:"TIER 2 — Test the mechanism", c:PURPLE, items:[
      "Patient-derived organoid co-cultures reconstituting the niche's cell pairs, with the predicted signalling readout",
      "Ligand or receptor blockade in co-culture, measuring the downstream target genes the interaction model predicts",
      "Regulon perturbation — CRISPRi or overexpression of the nominated transcription factor in the relevant lineage",
    ]},
    { tag:"TIER 3 — Test causality in vivo", c:TEAL, items:[
      "Syngeneic or humanised models with the pathway perturbed, read out spatially rather than only by tumour volume",
      "Combination with checkpoint blockade where the niche predicts resistance, testing whether disrupting it restores response",
    ]},
  ];
  let y = 1.56;
  tiers.forEach(t=>{
    const h = 0.44 + t.items.length*0.38;
    card(s, { x:M, y, w:W-2*M, h, fill:WHITE, line:LINE, shadow:true });
    s.addShape(pres.ShapeType.roundRect, { x:M+0.24, y:y+0.18, w:2.9, h:0.3, rectRadius:0.15,
      fill:{color:t.c}, line:{width:0} });
    body(s, { x:M+0.24, y:y+0.18, w:2.9, h:0.3, size:9, bold:true, color:WHITE, text:t.tag,
      align:"center", valign:"middle" });
    bullets(s, t.items, { x:M+3.34, y:y+0.16, w:W-2*M-3.62, h:h-0.3, size:10.5, gap:5 });
    y += h + 0.14;
  });
  card(s, { x:M, y:6.24, w:W-2*M, h:0.54, fill:TINT, line:TINT });
  body(s, { x:M+0.26, y:6.34, w:W-2*M-0.52, h:0.38, size:11, color:PLUM,
    text:"Tier 1 confirms that the thing we described is real. Tier 2 and 3 test whether it does what we claim. Conflating the two is the most common reason a computational atlas reads as descriptive to reviewers." });
  footer(s, "");
}

/* ----------------------------------------------------- 27 reproducibility */
{
  const s = content("Execution", "Compute, reproducibility and provenance");
  const cols = [
    { h:"Workflow", c:EOSIN, items:[
      "Snakemake DAG; per-rule conda environments and declared resources",
      "Per-sample fan-out on wildcards; one aggregation rule for cross-patient numbers",
      "Global seed; every stochastic step (UMAP, Leiden, NMF, permutation) seeded from it",
      "Restartable from intermediates — a cohort-scale run must survive a node failure",
    ]},
    { h:"Scale", c:PURPLE, items:[
      "On-disk backing (BPCells, h5ad backed mode) rather than in-memory cohort objects",
      "GPU nodes for graph deep learning and foundation-model fine-tuning",
      "Deconvolution is the memory-bound step and is resourced separately",
      "Cost and runtime profiled per rule, so the pipeline is portable to another cluster",
    ]},
    { h:"Provenance", c:TEAL, items:[
      "Container digests and package versions pinned and recorded per run",
      "Every analysis parameter versioned in configuration, never passed ad hoc",
      "Curation decisions logged with the evidence that produced them",
      "Public release of code, manifest and derived objects at submission",
    ]},
  ];
  const cw = (W-2*M-2*0.28)/3;
  cols.forEach((c,i)=>{
    const x = M + i*(cw+0.28);
    card(s, { x, y:1.6, w:cw, h:3.5, fill:LIGHT, line:LINE });
    s.addShape(pres.ShapeType.ellipse, { x:x+0.24, y:1.82, w:0.18, h:0.18, fill:{color:c.c}, line:{width:0} });
    body(s, { x:x+0.54, y:1.76, w:cw-0.78, h:0.3, size:14, bold:true, color:INK, text:c.h });
    bullets(s, c.items, { x:x+0.24, y:2.2, w:cw-0.48, h:2.7, size:10.5, gap:9 });
  });
  card(s, { x:M, y:5.3, w:W-2*M, h:1.16, fill:INK, line:INK });
  body(s, { x:M+0.3, y:5.46, w:W-2*M-0.6, h:0.3, size:12.5, bold:true, color:EOSIN,
    text:"Where LLM assistance is used, and where it is not" });
  body(s, { x:M+0.3, y:5.8, w:W-2*M-0.6, h:0.56, size:11, color:"C6BADD", lsm:1.16,
    text:"Prompts are versioned files with declared variables and provenance fingerprints; every call writes a record of its inputs, model and output. Manuscript-drafting prompts are gated: they render only against a verified results table built from executed analyses, so no number can enter a draft that was not computed." });
  footer(s, "");
}

/* --------------------------------------------------------- 28 timeline */
{
  const s = content("Execution", "Phasing over 24 months");
  const phases = [
    { q:"Months 1–4", t:"Harvest and curate", c:EOSIN, d:"Repository sweep, classification, structure verification, controlled-access applications submitted. Deliverable: the cohort manifest." },
    { q:"Months 3–9", t:"QC, harmonise, annotate", c:EOSIN, d:"Per-sample QC, within-cancer-type integration, unified cell ontology, deconvolution reliability benchmark. Deliverable: the harmonised object and reliability table." },
    { q:"Months 8–15", t:"Map niches", c:PURPLE, d:"Spatial domains, niche calling, recurrence testing across patients and cancer types. Deliverable: the pan-cancer niche atlas." },
    { q:"Months 12–19", t:"Mechanism layer", c:TEAL, d:"Cell–cell interaction, regulatory networks, co-expression modules, trajectories, all projected onto niches. Deliverable: the mechanism atlas." },
    { q:"Months 16–22", t:"Translate", c:AMBER, d:"Outcome models on external cohorts, target prioritisation, histology-based extension. Deliverable: the ranked target list." },
    { q:"Months 15–24", t:"Validate and write", c:AMBER, d:"Tier 1–3 experiments running in parallel from month 15. Deliverable: submission." },
  ];
  const rh = 0.68;
  phases.forEach((p,i)=>{
    const y = 1.6 + i*(rh+0.10);
    card(s, { x:M, y, w:W-2*M, h:rh, fill: i%2 ? WHITE : LIGHT, line:LINE });
    s.addShape(pres.ShapeType.roundRect, { x:M+0.22, y:y+0.17, w:1.42, h:0.34, rectRadius:0.17,
      fill:{color:p.c}, line:{width:0} });
    body(s, { x:M+0.22, y:y+0.17, w:1.42, h:0.34, size:9.5, bold:true, color:WHITE, text:p.q,
      align:"center", valign:"middle" });
    body(s, { x:M+1.84, y:y+0.19, w:2.5, h:0.32, size:13, bold:true, color:INK, text:p.t });
    body(s, { x:M+4.5, y:y+0.14, w:W-2*M-4.78, h:0.42, size:10, color:BODY, text:p.d, lsm:1.08 });
  });
  body(s, { x:M, y:6.28, w:W-2*M, h:0.44, size:11, color:PLUM, italic:true,
    text:"Phases overlap deliberately: validation begins while the mechanism layer is still being built, because the wet-lab timeline — not the compute — is the critical path." });
  footer(s, "");
}

/* -------------------------------------------------------- 29 limitations */
{
  const s = content("Rigour", "What this study will not be able to claim");
  body(s, { x:M, y:1.42, w:W-2*M, h:0.52, size:12.5, italic:true, color:BODY,
    text:"Stated at the outset, because a limitation discovered by a reviewer costs a year and one declared up front costs a paragraph." });

  const lims = [
    { h:"Not the first pan-cancer spatial atlas", d:"Several exist. The contribution is resolution, mechanism layering and statistical rigour — the framing must reflect that.", c:EOSIN },
    { h:"Public data carries selection bias", d:"Deposited cohorts over-represent resectable, treatment-naive, well-funded tumour types. Breast, glioma and pancreas dominate; gastric, oesophageal and rare cancers are thin. Pan-cancer means the cancers people have deposited.", c:EOSIN },
    { h:"Causality is out of reach computationally", d:"Spatial co-localisation and ligand–receptor scores generate hypotheses. Only the Tier 2–3 experiments test mechanism, and only for the handful of axes they cover.", c:PURPLE },
    { h:"Targeted panels constrain discovery", d:"Xenium, CosMx and MERFISH resolve single cells but measure a chosen gene set. Genes outside the panel are invisible, so discovery and resolution trade against each other.", c:PURPLE },
    { h:"Batch and biology remain partly entangled", d:"Within-cancer-type integration mitigates the confound; it does not remove it. Cross-cancer-type comparisons stay at the level of labels and summary statistics.", c:TEAL },
    { h:"Clinical annotation is sparse and inconsistent", d:"Treatment status, stage and outcome are unstated in most depositions and must be recovered from publications. Survival analyses will run on a fraction of the cohort.", c:TEAL },
  ];
  const cw = (W-2*M-0.28)/2, ch = 1.44;
  lims.forEach((o,i)=>{
    const x = M + (i%2)*(cw+0.28), y = 1.94 + Math.floor(i/2)*(ch+0.18);
    card(s, { x, y, w:cw, h:ch, fill: WHITE, line:LINE });
    s.addShape(pres.ShapeType.ellipse, { x:x+0.24, y:y+0.28, w:0.17, h:0.17, fill:{color:o.c}, line:{width:0} });
    body(s, { x:x+0.54, y:y+0.2, w:cw-0.8, h:0.34, size:12.5, bold:true, color:INK, text:o.h, lsm:1.0 });
    body(s, { x:x+0.54, y:y+0.62, w:cw-0.8, h:0.72, size:10.5, color:BODY, text:o.d, lsm:1.16 });
  });
  footer(s, "");
}

/* ------------------------------------------------------ 30 method stack */
{
  const s = content("Reference", "The method stack, one line per layer");
  const rows = [
    ["Ingest & curation", "Entrez Direct, BioPython, custom regex typing", "Auditable classification with retained evidence spans"],
    ["QC", "scanpy, Seurat v5, DoubletFinder / scDblFinder, SoupX / CellBender", "Per-sample MAD thresholds, not global cuts"],
    ["Integration", "Harmony, scVI / scANVI, Seurat v5 + BPCells, scArches", "Within cancer type, over patient; scIB metrics reported"],
    ["Annotation", "CellTypist, SingleR, Azimuth; Nicheformer fine-tuned", "Foundation models must beat scVI / Harmony / HVG on held-out patients"],
    ["Deconvolution", "cell2location, RCTD, CARD, Tangram", "Consensus of four; per-cell-type reliability flag from simulated mixtures"],
    ["Spatial domains", "GraphST, BayesSpace, STAGATE, SpaGCN", "Method chosen per platform; disagreement reported as uncertain"],
    ["Niches", "NicheCompass, BANKSY, CellCharter, UTAG / scNiche", "Recurrence across patients required before a niche is a finding"],
    ["Malignant cells", "Numbat, CopyKAT, inferCNV, SCEVAN; cNMF", "Numbat with alleles, CopyKAT expression-only; depth reported"],
    ["Interactions", "CellPhoneDB v5, CellChat v2, NicheNet, COMMOT, LIANA+", "Spatially constrained; co-localisation never called communication"],
    ["Regulatory networks", "pySCENIC, SCENIC+, cisTarget", "Enhancer support where scATAC exists; regulons are hypotheses"],
    ["Modules", "hdWGCNA, cNMF, MOFA+, SPARK-X / nnSVG", "Metacell aggregation; preservation tested across cancer types"],
    ["Dynamics", "CellRank 2, Monocle 3, Slingshot, PAGA, scVelo", "Velocity conditional; directions corroborated by two views"],
    ["Alignment", "PASTE / PASTE2, STalign, SLAT, STAligner", "Overlap fraction estimated, not assumed"],
    ["Histology", "UNI, Virchow, CONCH, Prov-GigaPath; HEST-Benchmark", "Predicted expression is imputation, never validation"],
    ["Statistics", "Torus-shift nulls, random-effects meta-analysis, BH-FDR", "Patient-level replication; one family; effect sizes with CIs"],
  ];
  const cols = [2.5, 4.6, 5.01];
  let y0 = 1.5;
  s.addShape(pres.ShapeType.rect, { x:M, y:y0, w:W-2*M, h:0.34, fill:{color:INK}, line:{width:0} });
  ["Layer","Primary tools","Governing rule"].forEach((h,i)=>{
    const x = M + cols.slice(0,i).reduce((a,b)=>a+b,0);
    s.addText(h, { x:x+0.12, y:y0, w:cols[i]-0.16, h:0.34, fontSize:9.5, bold:true, color:WHITE,
      fontFace:SANS, valign:"middle", isTextBox:true, margin:0 });
  });
  rows.forEach((r,ri)=>{
    const y = y0 + 0.34 + ri*0.328;
    s.addShape(pres.ShapeType.rect, { x:M, y, w:W-2*M, h:0.328,
      fill:{color: ri%2 ? WHITE : LIGHT}, line:{color:LINE, width:0.4} });
    r.forEach((c,ci)=>{
      const x = M + cols.slice(0,ci).reduce((a,b)=>a+b,0);
      s.addText(c, { x:x+0.12, y, w:cols[ci]-0.16, h:0.328, fontSize:8.6,
        bold: ci===0, color: ci===0 ? INK : (ci===2 ? PLUM : BODY), italic: ci===2,
        fontFace:SANS, valign:"middle", isTextBox:true, margin:0 });
    });
  });
  footer(s, "");
}

/* ------------------------------------------------------------ 31 close */
{
  const s = pres.addSlide();
  s.background = { color: INK };
  dotField(s, 0.3, 0.3, 4.2, 6.9, 150, [EOSIN, PURPLE, TEAL, "8E7BC3"], 0.072, 45);
  s.addShape(pres.ShapeType.ellipse, { x:1.2, y:2.5, w:2.2, h:2.2,
    fill:{color:TEAL, transparency:80}, line:{width:0} });

  s.addText("THE COMMITMENT", { x:5.1, y:1.5, w:7.4, h:0.3, fontSize:11, bold:true,
    charSpacing:2.6, color:EOSIN, fontFace:SANS, isTextBox:true, margin:0 });
  s.addText("Every claim traceable\nto a measurement", { x:5.1, y:1.9, w:7.4, h:1.5,
    fontSize:34, bold:true, color:WHITE, fontFace:HEAD, isTextBox:true, margin:0,
    lineSpacingMultiple:1.04 });

  const pts = [
    ["Scope", "All human cancers with paired spatial and single-cell data, curated rather than keyword-matched"],
    ["Method", "Benchmark-selected tools run in consensus, with disagreement reported instead of hidden"],
    ["Rigour", "Patient-level replication, autocorrelation-preserving nulls, one multiple-testing family"],
    ["Honesty", "Resolution limits, selection bias and the observational ceiling declared up front"],
  ];
  pts.forEach((p,i)=>{
    const y = 3.68 + i*0.67;
    s.addShape(pres.ShapeType.ellipse, { x:5.1, y:y+0.08, w:0.16, h:0.16,
      fill:{color:[EOSIN,PURPLE,TEAL,AMBER][i]}, line:{width:0} });
    s.addText(p[0], { x:5.42, y, w:1.15, h:0.3, fontSize:12, bold:true, color:WHITE,
      fontFace:SANS, isTextBox:true, margin:0 });
    s.addText(p[1], { x:6.62, y:y-0.02, w:5.9, h:0.56, fontSize:11, color:"B7A9D4",
      fontFace:SANS, isTextBox:true, margin:0, lineSpacingMultiple:1.12 });
  });
  s.addText("Dataset counts queried live from NCBI GEO on 9 September 2026. Method performance claims sourced to the published benchmarks cited on each slide.",
    { x:5.1, y:6.36, w:6.9, h:0.46, fontSize:9, color:"7A6E96", fontFace:SANS, italic:true,
      isTextBox:true, margin:0, lineSpacingMultiple:1.1 });
  footer(s, "");
}

pres.writeFile({ fileName: "pan-cancer-spatial-single-cell-plan.pptx" })
  .then(f => console.log("wrote", f, "|", N, "slides"));

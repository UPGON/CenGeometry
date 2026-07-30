const pptxgen = require("pptxgenjs");
const P = "/Users/hatzopou/Dropbox/MyDocs/Scripts/CenGeom/presentations/";
const F = P + "figs/";

// Palette taken from the model's own rendering colours, so the decks and the
// tool read as one thing.
const DARK = "1B1035", PURPLE = "6B3FA0", CYAN = "2FA8A8", GREEN = "5C7A45",
      BLUE = "5B8DD9", INK = "1A1A1A", MUTE = "6B7280", TINT = "F4F1F9",
      OK = "2E7D32", WARN = "D97706", BAD = "C62828", WHITE = "FFFFFF";
const HEAD = "Cambria", BODY = "Calibri";

const t = (s, o) => Object.assign({ fontFace: BODY, color: INK }, o);

function titleSlide(p, kicker, title, sub, foot) {
  const s = p.addSlide();
  s.background = { color: DARK };
  s.addText(kicker, t(kicker, { x: 0.9, y: 1.5, w: 11.5, h: 0.4, fontSize: 15,
    color: CYAN, bold: true, charSpacing: 3 }));
  s.addText(title, { x: 0.9, y: 1.95, w: 11.9, h: 2.15, fontSize: 41, bold: true,
    color: WHITE, fontFace: HEAD });
  s.addText(sub, t(sub, { x: 0.9, y: 4.25, w: 10.9, h: 1.1, fontSize: 18, color: "C9C4DB" }));
  if (foot) s.addText(foot, t(foot, { x: 0.9, y: 6.3, w: 11.5, h: 0.6, fontSize: 12,
    color: "9A93B5", italic: true }));
  return s;
}

function head(s, title, kicker) {
  if (kicker) s.addText(kicker, t(kicker, { x: 0.7, y: 0.34, w: 11.9, h: 0.3,
    fontSize: 12, color: PURPLE, bold: true, charSpacing: 2 }));
  s.addText(title, { x: 0.7, y: kicker ? 0.62 : 0.5, w: 11.9, h: 0.75, fontSize: 34,
    bold: true, color: INK, fontFace: HEAD });
}

function card(s, x, y, w, h, fill) {
  s.addShape("roundRect", { x, y, w, h, fill: { color: fill || TINT },
    line: { color: fill || TINT }, rectRadius: 0.08 });
}

function dot(s, x, y, label, colour, d) {
  const dd = d || 0.42;
  s.addShape("ellipse", { x, y, w: dd, h: dd, fill: { color: colour }, line: { color: colour } });
  s.addText(label, { x, y, w: dd, h: dd, fontSize: 13, bold: true, color: WHITE,
    align: "center", valign: "middle", fontFace: BODY, margin: 0 });
}

/* =====================================================================
   DECK 1 — Introduction for biologists
   ===================================================================== */
function deckOne() {
  const p = new pptxgen();
  p.layout = "LAYOUT_WIDE";
  p.author = "CenGeometry";

  titleSlide(p, "CENGEOMETRY", "Predicting how centriole\narchitecture responds to perturbation",
    "A geometric and mechanical model of the centriole cross-section, built from cryo-ET measurements",
    "A hypothesis-generating model — not a physics simulation. See the final slide for what it is not.");

  // --- the question
  let s = p.addSlide();
  head(s, "The question", "MOTIVATION");
  s.addText("A centriole is a ring of microtubule triplets held around a central cartwheel. " +
    "Every part is connected to every other part, so nothing can change alone.",
    t("", { x: 0.7, y: 1.5, w: 11.9, h: 0.8, fontSize: 17, color: MUTE }));
  const qs = [
    ["SAS-6 mutation", "changes cartwheel symmetry from 9-fold to 8-fold. What happens to the triplets?", PURPLE],
    ["Longer coiled-coil", "pushes the ring outward. Does the structure still close, and at what cost?", CYAN],
    ["C-tubule lost", "leaves doublets instead of triplets. Is that geometrically viable?", GREEN],
  ];
  qs.forEach(([h1, b, c], i) => {
    const y = 2.5 + i * 1.45;
    card(s, 0.7, y, 11.9, 1.25);
    dot(s, 1.0, y + 0.42, String(i + 1), c);
    s.addText(h1, t("", { x: 1.62, y: y + 0.2, w: 3.1, h: 0.42, fontSize: 17, bold: true, color: c, margin: 0 }));
    s.addText(b, t("", { x: 4.75, y: y + 0.2, w: 7.6, h: 0.85, fontSize: 15, color: INK, margin: 0 }));
  });
  s.addText("Answering these by intuition is unreliable — the constraints are geometric and they interact.",
    t("", { x: 0.7, y: 6.55, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: MUTE }));

  // --- the idea
  s = p.addSlide();
  head(s, "The idea: a linkage, not a drawing", "HOW IT WORKS");
  s.addImage({ path: F + "d_wt.png", x: 7.5, y: 1.35, w: 5.1, h: 5.1 });
  s.addText("One repeating unit is a set of rigid bodies joined by connections that can rotate:",
    t("", { x: 0.7, y: 1.4, w: 6.5, h: 0.6, fontSize: 16 }));
  const chain = [
    ["SAS-6 head", "sits on the hub ring", BLUE],
    ["Coiled-coil spoke", "radial at rest", BLUE],
    ["Pinhead", "grips A-tubule protofilaments 3 & 4", PURPLE],
    ["MT triplet", "A + B + C, one rigid body", "8A8F98"],
    ["Triplet base", "branches to the linker", GREEN],
    ["A-C linker", "reaches the next triplet", CYAN],
  ];
  chain.forEach(([n, d, c], i) => {
    const y = 2.15 + i * 0.66;
    s.addShape("ellipse", { x: 0.78, y: y + 0.07, w: 0.2, h: 0.2, fill: { color: c }, line: { color: c } });
    s.addText(n, t("", { x: 1.15, y, w: 2.7, h: 0.34, fontSize: 14, bold: true, margin: 0 }));
    s.addText(d, t("", { x: 3.75, y, w: 3.5, h: 0.34, fontSize: 13, color: MUTE, margin: 0 }));
  });
  card(s, 0.7, 6.25, 6.55, 0.85, "EDE7F6");
  s.addText("Because the units must close into a ring, changing one part forces the others to adapt. " +
    "The model finds the least-strained arrangement that still fits together.",
    t("", { x: 0.95, y: 6.38, w: 6.1, h: 0.65, fontSize: 13, bold: true, color: "4A2B7A", margin: 0 }));

  // --- measured, not estimated
  s = p.addSlide();
  head(s, "Every dimension is measured, not estimated", "PROVENANCE");
  s.addText("Read off a cryo-ET-derived schematic by a re-runnable script, calibrated against its 500 Å scale bar (4.132 Å per unit).",
    t("", { x: 0.7, y: 1.42, w: 11.9, h: 0.6, fontSize: 15, color: MUTE }));
  const rows = [
    [{ text: "Element", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "Measured value", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "Element", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "Measured value", options: { bold: true, color: WHITE, fill: { color: PURPLE } } }],
    ["SAS-6 dimer (total)", "50.0 nm", "Tubule radius", "11.9 nm"],
    ["  · head", "4.98 nm", "A→B spacing", "18.8 nm"],
    ["  · coiled coil", "45.0 nm", "B→C spacing", "19.6 nm"],
    ["Pinhead span", "21.3 nm", "Triplet axis vs spoke", "−56.9°"],
    ["Triplet base", "34.7 nm", "Protofilament pitch", "27.7°"],
    ["A-C linker arms", "14.6 / 11.6 nm", "Linker vertex angle", "138.7°"],
  ];
  s.addTable(rows, { x: 0.7, y: 2.0, w: 11.9, colW: [3.3, 2.65, 3.3, 2.65], fontSize: 14,
    fontFace: BODY, border: { type: "solid", color: "E3E0EC", pt: 1 }, rowH: 0.42, valign: "middle" });
  card(s, 0.7, 5.55, 11.9, 1.1, "FFF4E5");
  s.addText("One exception, stated plainly: the SAS-6 head-head spacing cannot be read off the drawing " +
    "(its centre is ambiguous), so it is fitted as the value that lets the measured unit close with least strain.",
    t("", { x: 0.95, y: 5.72, w: 11.4, h: 0.8, fontSize: 14, color: "8A5A00", margin: 0 }));

  // --- what measurement revealed
  s = p.addSlide();
  head(s, "Two things the measurement revealed", "FINDINGS FROM CALIBRATION");
  [["B and C tubules share the A-tubule's lattice",
    "They carry 9 protofilaments each — but at the same ~27.7° pitch as A's 13, i.e. 360/13, not 360/9. " +
    "All three tubules sit on one lattice; B and C simply stop short, leaving ~138° open where each fuses inward.", CYAN],
   ["The A-C linker is bent, not straight",
    "Two arms of 14.6 nm (to the C-tubule) and 11.6 nm (to the neighbouring A-tubule), meeting at 138.7°. " +
    "The triplet base attaches at that vertex. A straight linker would not reproduce the geometry.", GREEN]
  ].forEach(([h1, b, c], i) => {
    const y = 1.6 + i * 2.5;
    card(s, 0.7, y, 11.9, 2.15);
    dot(s, 1.05, y + 0.35, String(i + 1), c, 0.5);
    s.addText(h1, t("", { x: 1.8, y: y + 0.32, w: 10.4, h: 0.5, fontSize: 20, bold: true, color: c, fontFace: HEAD, margin: 0 }));
    s.addText(b, t("", { x: 1.8, y: y + 0.92, w: 10.4, h: 1.1, fontSize: 15, margin: 0 }));
  });
  s.addText("Neither was assumed — both fell out of measuring the schematic properly.",
    t("", { x: 0.7, y: 6.7, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: MUTE }));

  // --- mechanics
  s = p.addSlide();
  head(s, "The mechanical rules", "HOW IT DECIDES");
  [["Nothing stretches", "Segments have fixed length. Under compression they may buckle — bow so the ends come closer — but never extend.", PURPLE],
   ["Bonds have strengths", "Each connection is a spring stiffened by its bond strength, so weak bonds yield first and load distribution is a result, not an assumption.", CYAN],
   ["Joints have individual limits", "A rotation means different things at a 45 nm spoke and a 13 nm linker arm. Each is graded against its own thresholds.", GREEN],
   ["Connections cannot separate", "Each body is built onto the one that feeds it, so a connection point holds exactly \u2014 by construction, not by penalty.", BLUE]
  ].forEach(([h1, b, c], i) => {
    const x = 0.7 + (i % 2) * 6.15, y = 1.6 + Math.floor(i / 2) * 2.5;
    card(s, x, y, 5.75, 2.2);
    dot(s, x + 0.32, y + 0.32, "✓", c, 0.46);
    s.addText(h1, t("", { x: x + 0.95, y: y + 0.3, w: 4.5, h: 0.45, fontSize: 18, bold: true, color: c, fontFace: HEAD, margin: 0 }));
    s.addText(b, t("", { x: x + 0.35, y: y + 0.92, w: 5.1, h: 1.15, fontSize: 14, margin: 0 }));
  });
  s.addText("Grades: OK  ·  HARD  ·  SEVERE — assumed tolerances, not feasibility verdicts",
    t("", { x: 0.7, y: 6.65, w: 11.9, h: 0.4, fontSize: 14, bold: true, color: MUTE }));

  // --- per-joint limits
  s = p.addSlide();
  head(s, "Every joint has its own tolerance", "HOW STRAIN IS GRADED");
  s.addText("A rotation means different things at different connections. 15\u00b0 at the base of the " +
    "45 nm spoke moves its tip 11.6 nm \u2014 a whole tubule radius \u2014 while the same angle at a " +
    "13 nm linker arm moves it 3.4 nm. So each is graded against its own limits.",
    t("", { x: 0.7, y: 1.42, w: 11.9, h: 0.75, fontSize: 15, color: MUTE }));
  const jb = [
    [{ text: "Joint / contact", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "OK", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "HARD", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "SEVERE", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "Why", options: { bold: true, color: WHITE, fill: { color: PURPLE } } }],
    ["Linker to A-tubule", "8\u00b0", "8\u201320\u00b0", "> 20\u00b0", "grips a rigid, ordered lattice"],
    ["Linker to C-tubule", "8\u00b0", "8\u201320\u00b0", "> 20\u00b0", "same lattice, slightly weaker bond"],
    ["Pinhead to A-tubule", "10\u00b0", "10\u201322\u00b0", "> 22\u00b0", "lattice contact, larger footprint"],
    ["Spoke vs radial", "15\u00b0", "15\u201335\u00b0", "> 35\u00b0", "many SAS-6 rings stack axially"],
    ["Pinhead vs spoke", "15\u00b0", "15\u201330\u00b0", "> 30\u00b0", "multi-protein, three contacts"],
    ["Triplet axis vs spoke", "20\u00b0", "20\u201340\u00b0", "> 40\u00b0", "a composite, not one interface"],
    ["Base vs spoke", "20\u00b0", "20\u201340\u00b0", "> 40\u00b0", "designed to reorient"],
  ];
  s.addTable(jb, { x: 0.7, y: 2.3, w: 11.9, colW: [2.9, 1.0, 1.3, 1.2, 5.5], fontSize: 12.5,
    fontFace: BODY, border: { type: "solid", color: "E3E0EC", pt: 1 }, rowH: 0.45, valign: "middle" });
  card(s, 0.7, 5.9, 11.9, 1.3, "FFEBEE");
  s.addText("These are assumed tolerances, NOT verdicts on feasibility. The values were reasoned, " +
    "not measured. A SAS-6 spoke shortened by 17 nm builds real centrioles, yet the model grades its " +
    "linker contacts SEVERE while reporting zero clashes and a perfectly buildable geometry. " +
    "Read them as \u201chow far from wild type\u201d and judge from the geometry itself.",
    t("", { x: 0.95, y: 6.05, w: 11.4, h: 1.05, fontSize: 13.5, color: "8E1B1B", margin: 0 }));

  // --- validation
  s = p.addSlide();
  head(s, "The model was never told 9-fold is special", "VALIDATION");
  s.addImage({ path: F + "d_validation.png", x: 0.7, y: 1.5, w: 6.5, h: 4.1 });
  s.addText("Sweeping each parameter independently, the measured wild-type value sits at a strain minimum every time:",
    t("", { x: 7.5, y: 1.45, w: 5.2, h: 0.85, fontSize: 14, color: MUTE }));
  const v = [
    [{ text: "Parameter", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "Measured", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "Best", options: { bold: true, color: WHITE, fill: { color: PURPLE } } }],
    ["Symmetry", "9", "9"],
    ["SAS-6 coiled coil", "45.0 nm", "45 nm"],
    ["Triplet base", "34.7 nm", "34.7 nm"],
    ["A-tubule protofilaments", "13", "13"],
    ["Tubules per blade", "3", "3"],
  ];
  s.addTable(v, { x: 7.5, y: 2.4, w: 5.15, colW: [2.55, 1.35, 1.25], fontSize: 13, fontFace: BODY,
    border: { type: "solid", color: "E3E0EC", pt: 1 }, rowH: 0.4, valign: "middle" });
  card(s, 7.5, 5.15, 5.15, 1.5, "E8F5E9");
  s.addText("Wild-type diameter comes out at 255 nm, against roughly 250 nm measured in real centrioles — " +
    "and diameter was never fitted.",
    t("", { x: 7.72, y: 5.32, w: 4.75, h: 1.2, fontSize: 14, bold: true, color: "1B5E20", margin: 0 }));
  s.addText("Strain rises steeply either side of 9: 11.1° at 8-fold, 9.1° at 10-fold, against 1.3° at 9.",
    t("", { x: 0.7, y: 5.75, w: 6.5, h: 0.5, fontSize: 14, italic: true, color: MUTE }));

  // --- outputs
  s = p.addSlide();
  head(s, "What it tells you", "OUTPUTS");
  const o = [
    [{ text: "Readout", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "Biological meaning", options: { bold: true, color: WHITE, fill: { color: PURPLE } } }],
    ["Joint rotation", "How far each connection had to turn, graded OK / HARD / SEVERE against its own limits (assumed tolerances, not feasibility)"],
    ["Bond load", "Which connection carries most strain — i.e. which would rupture first"],
    ["Buckling", "Which segments were forced to bow because space became tight"],
    ["Clashes", "Microtubules overlapping in space. Anything above zero is physically impossible"],
    ["Clearance", "How close the linker, base or spoke passes to a microtubule; negative means passing through one"],
    ["Diameter", "Outer, A-tubule ring, and lumen, in nm — directly comparable to EM measurements"],
    ["Unattached triplets", "Triplets the cartwheel could not reach — expected when the two symmetries differ"],
  ];
  s.addTable(o, { x: 0.7, y: 1.5, w: 11.9, colW: [2.9, 9.0], fontSize: 14, fontFace: BODY,
    border: { type: "solid", color: "E3E0EC", pt: 1 }, rowH: 0.55, valign: "middle" });
  s.addText("Read grades comparatively — mutant against wild type — rather than as absolute verdicts.",
    t("", { x: 0.7, y: 6.5, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: MUTE }));

  // --- worked prediction
  s = p.addSlide();
  head(s, "A worked prediction: symmetry mismatch", "WHAT IT IS FOR");
  s.addImage({ path: F + "d_mismatch.png", x: 0.55, y: 1.35, w: 12.2, h: 4.3 });
  card(s, 0.7, 5.85, 11.9, 1.25, "EDE7F6");
  s.addText("If SAS-6 makes an 8-fold cartwheel but nine triplets are still built, the model predicts the " +
    "triplet ring keeps wild-type spacing and diameter (255 nm, unchanged) while the cartwheel absorbs the " +
    "entire mismatch — spokes strained to HARD, one triplet left with no pinhead at all. A fully 8-fold " +
    "centriole instead shrinks to 233 nm.",
    t("", { x: 0.95, y: 6.0, w: 11.4, h: 1.0, fontSize: 14, color: "4A2B7A", margin: 0 }));

  // --- how to use
  s = p.addSlide();
  head(s, "Using it takes no coding", "PRACTICALITIES");
  [["Double-click the launcher", "Sets everything up the first time and opens in your browser.", "1"],
   ["Cross-section tab", "Change any parameter; the geometry redraws with every metric.", "2"],
   ["Parameter scan tab", "Vary one parameter; all metrics plotted, with cross-sections and a wild-type reference.", "3"],
   ["2-parameter scan tab", "Cross any two — e.g. coiled-coil length against symmetry — as heatmaps or curves.", "4"]
  ].forEach(([h1, b, n], i) => {
    const y = 1.6 + i * 1.28;
    card(s, 0.7, y, 11.9, 1.1);
    dot(s, 1.02, y + 0.32, n, PURPLE, 0.46);
    s.addText(h1, t("", { x: 1.68, y: y + 0.16, w: 3.6, h: 0.4, fontSize: 17, bold: true, margin: 0 }));
    s.addText(b, t("", { x: 5.3, y: y + 0.16, w: 7.0, h: 0.8, fontSize: 14, color: MUTE, margin: 0 }));
  });
  s.addText("Sliders for symmetry, typed boxes for every distance, and CSV download on every scan.",
    t("", { x: 0.7, y: 6.85, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: MUTE }));

  // --- caveats
  s = p.addSlide();
  s.background = { color: DARK };
  s.addText("What this model is not", { x: 0.9, y: 0.75, w: 11.5, h: 0.8, fontSize: 34,
    bold: true, color: WHITE, fontFace: HEAD });
  [["It is not a physics simulation", "It is a geometric and mechanical model for generating hypotheses. Bond strengths are a rank ordering, not measured energies."],
   ["The joint thresholds are reasoned, not measured", "No sub-tomogram angular variance was available to calibrate OK / HARD / SEVERE. Their rank ordering is defensible; the absolute numbers are not."],
   ["It is a single 2D cross-section", "The real centriole twists along its length. Anything longitudinal is outside its scope."],
   ["Grades are not feasibility verdicts", "The tolerances were reasoned, not measured. Structures assemble in conditions graded SEVERE."],
   ["Calibration rests on one schematic", "Protofilament positions carry about 2° (~0.4 nm) of reading noise, which sets the precision floor."]
  ].forEach(([h1, b], i) => {
    const y = 1.72 + i * 0.99;
    s.addShape("roundRect", { x: 0.9, y, w: 11.5, h: 0.86, fill: { color: "2A1B4D" },
      line: { color: "2A1B4D" }, rectRadius: 0.06 });
    s.addText(h1, t("", { x: 1.2, y: y + 0.06, w: 4.55, h: 0.74, fontSize: 14.5, bold: true, color: CYAN, margin: 0 }));
    s.addText(b, t("", { x: 5.85, y: y + 0.06, w: 6.3, h: 0.74, fontSize: 12, color: "C9C4DB", margin: 0 }));
  });
  s.addText("Used within those limits, it answers a question intuition cannot: given everything is connected, what has to give?",
    t("", { x: 0.9, y: 6.85, w: 11.5, h: 0.45, fontSize: 15, italic: true, color: CYAN }));

  return p.writeFile({ fileName: P + "01_CenGeometry_Introduction.pptx" });
}

/* =====================================================================
   DECK 2 — Critical evaluation
   ===================================================================== */
function deckTwo() {
  const p = new pptxgen();
  p.layout = "LAYOUT_WIDE";

  titleSlide(p, "CRITICAL EVALUATION", "What this model cannot\ncurrently be trusted to tell you",
    "An assessment of CenGeometry's limitations, ranked by severity",
    "Caveat on this review: it was written by the same process that built the tool. Self-assessment is weaker than independent review — an outside critic should repeat it.");

  // --- fair summary
  let s = p.addSlide();
  head(s, "First, what it genuinely achieves", "A FAIR STARTING POINT");
  [["Traceable calibration", "Every dimension but one is derived from a schematic by a re-runnable script. That is better practice than most models of this kind.", OK],
   ["Real structural findings", "Measurement revealed the shared protofilament lattice and the bent A-C linker. Neither was assumed.", OK],
   ["An independent number that matches", "Wild-type diameter of 255 nm against ~250 nm measured. Diameter was never fitted.", OK],
   ["Honest internal documentation", "Failed analyses are parked with their blockers recorded rather than quietly dropped.", OK]
  ].forEach(([h1, b, c], i) => {
    const x = 0.7 + (i % 2) * 6.15, y = 1.6 + Math.floor(i / 2) * 2.45;
    card(s, x, y, 5.75, 2.15, "E8F5E9");
    s.addText(h1, t("", { x: x + 0.35, y: y + 0.25, w: 5.1, h: 0.45, fontSize: 18, bold: true, color: "1B5E20", fontFace: HEAD, margin: 0 }));
    s.addText(b, t("", { x: x + 0.35, y: y + 0.82, w: 5.1, h: 1.15, fontSize: 14, margin: 0 }));
  });
  s.addText("The problems below are real, but they are problems of interpretation and calibration — not of workmanship.",
    t("", { x: 0.7, y: 6.6, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: MUTE }));

  // --- THE CRUX
  s = p.addSlide();
  head(s, "The headline validation is largely circular", "FUNDAMENTAL · THE CRUX");
  card(s, 0.7, 1.5, 11.9, 1.5, "FFEBEE");
  s.addText("The claim: “the measured wild-type value sits at a strain minimum for five parameters independently.”\n" +
    "The problem: the model's energy is the sum of bond gaps and deviations from rest angles — and those rest " +
    "angles were themselves measured from the wild-type geometry.",
    t("", { x: 0.95, y: 1.68, w: 11.4, h: 1.2, fontSize: 15, color: "8E1B1B", margin: 0 }));
  s.addText("Why that makes the minimum near-automatic", t("", { x: 0.7, y: 3.2, w: 11.9, h: 0.4,
    fontSize: 17, bold: true, color: INK, fontFace: HEAD }));
  [["At the measured configuration, every angle deviation is zero by definition.",
    "The energy is a sum of squares, so zero deviation with near-zero closure gap is its global minimum. Any departure can only raise it."],
   ["So “wild type is optimal” restates how the objective was defined.",
    "It is not an independent prediction. A model built this way would place a minimum at whatever geometry it was calibrated on."]
  ].forEach(([h1, b], i) => {
    const y = 3.7 + i * 1.05;
    dot(s, 0.7, y, String(i + 1), BAD, 0.4);
    s.addText(h1, t("", { x: 1.28, y: y - 0.04, w: 11.3, h: 0.4, fontSize: 15, bold: true, margin: 0 }));
    s.addText(b, t("", { x: 1.28, y: y + 0.36, w: 11.3, h: 0.5, fontSize: 14, color: MUTE, margin: 0 }));
  });
  card(s, 0.7, 5.85, 11.9, 1.3, "E8F5E9");
  s.addText("What survives: two things are genuinely non-trivial. That a geometry measured from TWO units closes " +
    "into a NINE-fold ring at all (to 0.022 nm) — though one parameter was fitted to help. And that the resulting " +
    "diameter, 255 nm, matches the ~250 nm measured independently. The diameter agreement is the real validation; " +
    "the strain minima are largely bookkeeping.",
    t("", { x: 0.95, y: 6.0, w: 11.4, h: 1.05, fontSize: 14, color: "1B5E20", margin: 0 }));

  // --- calibration
  s = p.addSlide();
  head(s, "Everything rests on one hand-drawn schematic", "FUNDAMENTAL");
  [["Sample size of one", "All dimensions come from a single drawing of two units. There is no biological replicate, no error bar on any input, and no way to propagate uncertainty to the outputs.", BAD],
   ["The source contradicts itself", "The schematic's two triplets sit 41.4° apart while its two spokes sit exactly 40° apart. A 9-fold structure cannot be both. The model resolves this by enforcing exact symmetry and absorbing the discrepancy elsewhere.", BAD],
   ["Reading noise sets the floor", "Protofilament positions were recovered from text-label anchors, carrying ~2° (~0.4 nm) of error. Numbers are reported to 0.01 nm — three orders of magnitude finer than the calibration supports.", WARN],
   ["One input is fitted, not measured", "The SAS-6 head-head spacing cannot be read off the drawing and is tuned so wild type closes. That is the same parameter the closure validation then relies on.", WARN]
  ].forEach(([h1, b, c], i) => {
    const y = 1.55 + i * 1.34;
    card(s, 0.7, y, 11.9, 1.18);
    s.addShape("ellipse", { x: 1.0, y: y + 0.42, w: 0.34, h: 0.34, fill: { color: c }, line: { color: c } });
    s.addText(h1, t("", { x: 1.55, y: y + 0.1, w: 3.45, h: 0.95, fontSize: 16, bold: true, color: c, margin: 0 }));
    s.addText(b, t("", { x: 5.1, y: y + 0.14, w: 7.3, h: 0.95, fontSize: 13.5, margin: 0 }));
  });

  // --- 2D
  s = p.addSlide();
  head(s, "It is one flat slice of a helical object", "FUNDAMENTAL");
  s.addText("The centriole twists along its length. A single cross-section cannot represent that, and the " +
    "consequence is not hypothetical — it already blocked a real investigation.",
    t("", { x: 0.7, y: 1.45, w: 11.9, h: 0.7, fontSize: 16, color: MUTE }));
  card(s, 0.7, 2.35, 11.9, 2.0, "FFEBEE");
  s.addText("Worked example: the blooming question", t("", { x: 1.0, y: 2.55, w: 11.3, h: 0.4,
    fontSize: 18, bold: true, color: "8E1B1B", fontFace: HEAD, margin: 0 }));
  s.addText("A published iris-like motion is described as arising from A-C linker twist in the LONGITUDINAL " +
    "direction. The model was run and returned a confident negative — no in-plane iris mode exists. That answer " +
    "is irrelevant to the claim, because the claim is about a motion the model cannot represent. The investigation " +
    "was parked for exactly this reason.",
    t("", { x: 1.0, y: 3.05, w: 11.3, h: 1.2, fontSize: 14.5, margin: 0 }));
  s.addText("The general risk", t("", { x: 0.7, y: 4.55, w: 11.9, h: 0.4, fontSize: 17, bold: true, fontFace: HEAD }));
  s.addText("A 2D model still returns a confident-looking number when asked a 3D question. Nothing in the output " +
    "flags that the question was out of scope. Any result about twist, axial gradients, procentriole elongation " +
    "or stacked cartwheel layers should be treated as unanswerable rather than answered.",
    t("", { x: 0.7, y: 5.0, w: 11.9, h: 1.0, fontSize: 15 }));
  s.addText("Out-of-plane tilt also gives real interfaces accommodation the model lacks — so strain is probably overestimated.",
    t("", { x: 0.7, y: 6.2, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: MUTE }));

  // --- fixable
  s = p.addSlide();
  head(s, "Fixable: the numbers behind the verdicts", "SECOND TIER");
  [["Joint thresholds are invented", "OK / HARD / SEVERE limits were reasoned, not measured. Tested: scaling all bands 0.5×–2× leaves the qualitative conclusion intact, but flips the spoke's label between OK and HARD. So the ranking is robust; the labels are not.", WARN],
   ["Bond strengths are ordinal", "Six bonds carry values from 1.00 to 0.25 encoding only a rank order. The model then reports a specific “closest to rupture” bond and a numeric load. That precision is not supported by the input — only the ordering is.", WARN],
   ["No absolute energy scale", "Without bond energies in kT, nothing converts to a predicted fluctuation amplitude comparable with data. This is exactly what made the soft-mode analysis uninterpretable and led to it being abandoned.", WARN]
  ].forEach(([h1, b, c], i) => {
    const y = 1.55 + i * 1.72;
    card(s, 0.7, y, 11.9, 1.55, "FFF8E1");
    s.addText(h1, t("", { x: 1.05, y: y + 0.2, w: 11.0, h: 0.42, fontSize: 18, bold: true, color: "8A5A00", fontFace: HEAD, margin: 0 }));
    s.addText(b, t("", { x: 1.05, y: y + 0.7, w: 11.0, h: 0.8, fontSize: 14, margin: 0 }));
  });
  s.addText("A band_sensitivity() helper exists and should be run on any conclusion before it is reported.",
    t("", { x: 0.7, y: 6.85, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: MUTE }));

  // --- robustness
  s = p.addSlide();
  head(s, "Fixable: the solver is not dependable at the edges", "SECOND TIER");
  s.addText("Solve cost is driven by convergence difficulty, not problem size — so it is unpredictable:",
    t("", { x: 0.7, y: 1.45, w: 11.9, h: 0.4, fontSize: 15, color: MUTE }));
  s.addChart(p.ChartType.bar, [{ name: "Solve time (s)", labels: ["typical", "n_pf = 9", "n_pf = 18"],
    values: [2, 28.5, 257] }], { x: 0.7, y: 2.0, w: 5.6, h: 3.2, barDir: "col",
    showTitle: true, title: "Wall-clock per configuration", titleFontSize: 13,
    chartColors: [WARN], showValue: true, dataLabelPosition: "outEnd", dataLabelFontSize: 11,
    showLegend: false, catAxisLabelColor: MUTE, valAxisLabelColor: MUTE,
    valGridLine: { color: "EAEAEA", size: 1 }, catGridLine: { style: "none" } });
  [["Some inputs never converge", "n_pf_A = 20 fails outright. Non-convergence is now surfaced, but the geometry is still drawn and every metric still reported."],
   ["Accuracy was traded for responsiveness", "Solver effort is capped in the interface so it cannot stall. Extreme configurations therefore return lower-confidence answers by design."],
   ["An earlier expectation was simply wrong", "Clashes were predicted at low protofilament counts. Tested 9 through 18: zero clashes at every value. The claim came from a superseded architecture and had not been rechecked."]
  ].forEach(([h1, b], i) => {
    const y = 2.0 + i * 1.15;
    s.addText(h1, t("", { x: 6.6, y, w: 6.0, h: 0.38, fontSize: 15, bold: true, color: WARN, margin: 0 }));
    s.addText(b, t("", { x: 6.6, y: y + 0.36, w: 6.0, h: 0.75, fontSize: 13, margin: 0 }));
  });
  s.addText("The last point is the important one: a stated result survived an architecture change without being re-tested.",
    t("", { x: 0.7, y: 5.6, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: BAD }));

  // --- what must change
  s = p.addSlide();
  head(s, "What would have to change", "BEFORE PUBLICATION");
  const w = [
    [{ text: "Requirement", options: { bold: true, color: WHITE, fill: { color: BAD } } },
     { text: "Why it is blocking", options: { bold: true, color: WHITE, fill: { color: BAD } } },
     { text: "Effort", options: { bold: true, color: WHITE, fill: { color: BAD } } }],
    ["Validate against data NOT used in calibration", "Removes the circularity. Predict a mutant or a species the model has never seen, then measure it.", "High"],
    ["Calibrate joint bands on sub-tomogram variance", "Converts invented thresholds into measured ones; makes grades meaningful in absolute terms.", "Medium"],
    ["Propagate input uncertainty to outputs", "Every number is currently a point estimate from one noisy drawing. Error bars are needed.", "Medium"],
    ["Put bond strengths on an absolute scale", "Required before any force, rupture or fluctuation claim can be quantitative.", "High"],
    ["State the scope boundary in the output", "So a 3D question cannot silently receive a 2D answer.", "Low"],
  ];
  s.addTable(w, { x: 0.7, y: 1.5, w: 11.9, colW: [3.7, 6.6, 1.6], fontSize: 13.5, fontFace: BODY,
    border: { type: "solid", color: "E3E0EC", pt: 1 }, rowH: 0.72, valign: "middle" });
  s.addText("None of these invalidate the tool as a hypothesis generator. They are what stand between it and a quantitative claim.",
    t("", { x: 0.7, y: 6.4, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: MUTE }));

  // --- safe vs not
  s = p.addSlide();
  s.background = { color: DARK };
  s.addText("What you can and cannot conclude today", { x: 0.9, y: 0.7, w: 11.5, h: 0.8,
    fontSize: 32, bold: true, color: WHITE, fontFace: HEAD });
  s.addShape("roundRect", { x: 0.9, y: 1.75, w: 5.6, h: 4.4, fill: { color: "16351F" }, line: { color: OK }, rectRadius: 0.08 });
  s.addText("SAFE", { x: 1.2, y: 1.95, w: 5.0, h: 0.4, fontSize: 17, bold: true, color: "7BC67B", fontFace: BODY });
  s.addText([
    { text: "Relative comparisons — mutant against wild type, one perturbation against another", options: { bullet: true, breakLine: true } },
    { text: "The ordering of which connection is most loaded", options: { bullet: true, breakLine: true } },
    { text: "Whether a proposed architecture is geometrically possible at all", options: { bullet: true, breakLine: true } },
    { text: "Diameter predictions, the one output with independent support", options: { bullet: true, breakLine: true } },
    { text: "Generating hypotheses worth testing at the bench", options: { bullet: true } },
  ], t("", { x: 1.2, y: 2.45, w: 5.0, h: 3.5, fontSize: 13.5, color: "D8E8D8", paraSpaceAfter: 8 }));

  s.addShape("roundRect", { x: 6.85, y: 1.75, w: 5.55, h: 4.4, fill: { color: "3A1620" }, line: { color: BAD }, rectRadius: 0.08 });
  s.addText("NOT SAFE", { x: 7.15, y: 1.95, w: 5.0, h: 0.4, fontSize: 17, bold: true, color: "F09A9A", fontFace: BODY });
  s.addText([
    { text: "“Wild type is optimal” as an independent finding — it is largely built in", options: { bullet: true, breakLine: true } },
    { text: "Absolute strain values, forces, or rupture thresholds", options: { bullet: true, breakLine: true } },
    { text: "Any claim about twist, axial gradients or assembly over time", options: { bullet: true, breakLine: true } },
    { text: "Precision beyond about 1 nm and a few degrees", options: { bullet: true, breakLine: true } },
    { text: "Negative results on questions the 2D geometry cannot represent", options: { bullet: true } },
  ], t("", { x: 7.15, y: 2.45, w: 5.0, h: 3.5, fontSize: 13.5, color: "F0D8D8", paraSpaceAfter: 8 }));

  s.addText("Bottom line: a well-built hypothesis generator that has not yet been independently validated. Use it to decide what to measure, not to conclude what is true.",
    t("", { x: 0.9, y: 6.45, w: 11.5, h: 0.7, fontSize: 15, italic: true, color: CYAN }));

  return p.writeFile({ fileName: P + "02_CenGeometry_Critical_Evaluation.pptx" });
}

/* =====================================================================
   DECK 3 — Expert recommendations
   ===================================================================== */
function deckThree() {
  const p = new pptxgen();
  p.layout = "LAYOUT_WIDE";

  titleSlide(p, "EXPERT REVIEW", "Where to take this model,\nand what to ask it",
    "Prioritised improvements, published knowledge worth building in, and five testable predictions",
    "Citations below were drawn from the project's own literature notes. Each should be checked against the source before use in a manuscript.");

  // --- assessment
  let s = p.addSlide();
  head(s, "Assessment in one slide", "STARTING POSITION");
  card(s, 0.7, 1.5, 11.9, 1.35, "EDE7F6");
  s.addText("A geometrically careful, honestly documented model whose main weakness is not its mechanics but its " +
    "evidence base: one calibration source, invented thresholds, ordinal bond strengths, and a validation that is " +
    "largely internal. The mechanics are ready for more than the calibration currently supports.",
    t("", { x: 0.98, y: 1.68, w: 11.35, h: 1.05, fontSize: 15, color: "4A2B7A", margin: 0 }));
  [["Ready to use now", "Comparative predictions, feasibility screening, hypothesis generation", OK],
   ["Needs work first", "Absolute forces, fluctuation amplitudes, anything longitudinal", WARN],
   ["The single highest-value fix", "External validation against data not used in calibration", BAD]
  ].forEach(([h1, b, c], i) => {
    const x = 0.7 + i * 4.07;
    card(s, x, 3.1, 3.75, 2.0);
    dot(s, x + 0.28, 3.35, String(i + 1), c, 0.45);
    s.addText(h1, t("", { x: x + 0.28, y: 3.95, w: 3.2, h: 0.5, fontSize: 16, bold: true, color: c, fontFace: HEAD, margin: 0 }));
    s.addText(b, t("", { x: x + 0.28, y: 4.5, w: 3.2, h: 0.5, fontSize: 13, margin: 0 }));
  });
  s.addText("The tool's own documentation already identifies most of these. The recommendation is to act on them in order of leverage.",
    t("", { x: 0.7, y: 5.4, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: MUTE }));

  // --- Part A
  s = p.addSlide();
  head(s, "Part A — Improvements, by value for effort", "ROADMAP");
  const a = [
    [{ text: "#", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "Improvement", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "What it unlocks", options: { bold: true, color: WHITE, fill: { color: PURPLE } } },
     { text: "Effort", options: { bold: true, color: WHITE, fill: { color: PURPLE } } }],
    ["1", "Hold-out validation", "Breaks the circularity. Predict a system never used in calibration, then measure it. Without this, no result is independently supported.", "High"],
    ["2", "Bands from sub-tomogram variance", "Turns OK / HARD / SEVERE from reasoned guesses into measured tolerances. The observed particle-to-particle spread IS the OK band.", "Medium"],
    ["3", "Bond energies in kT", "Enables forces, rupture predictions and fluctuation amplitudes. Explicitly the blocker that killed the soft-mode analysis.", "High"],
    ["4", "2.5D stack of coupled slices", "Several cross-sections at fixed axial spacing with a twist increment. Reaches longitudinal questions at a fraction of full-3D cost.", "Medium"],
    ["5", "Uncertainty propagation", "Resample inputs within their ~0.4 nm noise and report output ranges instead of point estimates.", "Low"],
    ["6", "Per-state geometries", "Procentriole, distal cartwheel-free region, Chlamydomonas and Trichonympha as their own calibrations, not perturbations of one mature human slice.", "Medium"],
  ];
  s.addTable(a, { x: 0.7, y: 1.5, w: 11.9, colW: [0.5, 3.0, 6.9, 1.5], fontSize: 12.5, fontFace: BODY,
    border: { type: "solid", color: "E3E0EC", pt: 1 }, rowH: 0.62, valign: "middle" });
  s.addText("Do 5 first — it is cheap and immediately makes every other number more honest. Then 2, then 1.",
    t("", { x: 0.7, y: 6.55, w: 11.9, h: 0.4, fontSize: 14, italic: true, color: MUTE }));

  // --- Part B
  s = p.addSlide();
  head(s, "Part B — Published knowledge worth building in", "GROUNDING");
  [["Stage-resolved ring diameters", "U-ExM measurements give PLK4 77, SAS-6 75, STIL 89, CPAP 120 and CEP135 85 nm. These are direct calibration targets for assembly-stage geometries, and CEP135 spanning 85→182 nm across assembly bounds how far the triplet base can reach.", CYAN],
   ["Protein identity per element", "Assigning real proteins to each modelled body — CEP135 to the triplet base, SAS-6 to hub and spoke — lets mutant data enter as constraints and makes every prediction addressable by a specific reagent.", PURPLE],
   ["A-C linker as an integrity determinant", "The linker is reported to control centriole structural integrity and duplication. The model independently makes it the load-bearing element, so published linker mutants are a ready-made validation set.", GREEN],
   ["Known mutant phenotypes", "Any documented symmetry or length mutant with a measured diameter is a free hold-out test. This is the cheapest available route to non-circular validation.", BLUE]
  ].forEach(([h1, b, c], i) => {
    const y = 1.5 + i * 1.36;
    card(s, 0.7, y, 11.9, 1.2);
    s.addShape("ellipse", { x: 1.0, y: y + 0.43, w: 0.34, h: 0.34, fill: { color: c }, line: { color: c } });
    s.addText(h1, t("", { x: 1.55, y: y + 0.1, w: 3.35, h: 0.95, fontSize: 15.5, bold: true, color: c, margin: 0 }));
    s.addText(b, t("", { x: 5.05, y: y + 0.12, w: 7.35, h: 1.0, fontSize: 13, margin: 0 }));
  });
  s.addText("Verify each source directly before citing — these were taken from the project's working notes, not from the papers.",
    t("", { x: 0.7, y: 6.95, w: 11.9, h: 0.4, fontSize: 13, italic: true, color: BAD }));

  // --- Part C intro + Q1,2
  s = p.addSlide();
  head(s, "Part C — Five testable questions", "PREDICTIONS");
  s.addText("Chosen for making non-obvious predictions that could plausibly be wrong. A prediction that restates a known fact is not a test.",
    t("", { x: 0.7, y: 1.4, w: 11.9, h: 0.55, fontSize: 14, color: MUTE }));
  [["Q1", "Does the triplet ring or the cartwheel absorb a symmetry mismatch?",
    "PREDICTION  An 8-fold cartwheel with 9 triplets keeps wild-type diameter (255 nm, versus 233 nm if the whole centriole were 8-fold), leaves one triplet with no pinhead, and strains the spoke to 18.5° while triplets stay under 8°.",
    "TEST  A SAS-6 mutant with altered cartwheel symmetry; cryo-ET; measure diameter and count spoke-to-triplet attachments.",
    "FALSIFIED IF  diameter falls toward 233 nm, i.e. the triplet ring follows the cartwheel instead of resisting it.", PURPLE],
   ["Q2", "Which connection fails first under stress?",
    "PREDICTION  The triplet-base-to-A-C-linker junction carries the highest load in every perturbation tested — the weakest bond yields first, ahead of the linker-to-tubule contacts.",
    "TEST  Graded depletion or point mutation of each junction; look for which linkage is lost first, or which is most disordered in sub-tomogram averages of stressed centrioles.",
    "FALSIFIED IF  a tubule contact fails before the base-linker junction.", CYAN]
  ].forEach(([n, q, pr, te, fa, c], i) => {
    const y = 1.95 + i * 2.45;
    card(s, 0.7, y, 11.9, 2.25);
    dot(s, 1.0, y + 0.22, n, c, 0.52);
    s.addText(q, t("", { x: 1.72, y: y + 0.18, w: 10.6, h: 0.42, fontSize: 17, bold: true, color: c, fontFace: HEAD, margin: 0 }));
    s.addText(pr, t("", { x: 1.72, y: y + 0.66, w: 10.6, h: 0.62, fontSize: 12.5, margin: 0 }));
    s.addText(te, t("", { x: 1.72, y: y + 1.3, w: 10.6, h: 0.48, fontSize: 12.5, color: MUTE, margin: 0 }));
    s.addText(fa, t("", { x: 1.72, y: y + 1.78, w: 10.6, h: 0.36, fontSize: 12.5, color: BAD, margin: 0 }));
  });

  // --- Q3,4,5
  s = p.addSlide();
  head(s, "Five testable questions, continued", "PREDICTIONS");
  [["Q3", "Is a doublet centriole a shrunken triplet centriole?",
    "PREDICTION  Removing the C-tubule collapses diameter from 255 to 213 nm — a 42 nm drop — and raises joint strain from 1.3° to 36°, i.e. doublets are strongly disfavoured rather than merely smaller.",
    "TEST  Measure diameter in systems with native doublets, or after C-tubule loss.  FALSIFIED IF  doublet diameter is close to triplet diameter.", GREEN],
   ["Q4", "Does protofilament number set centriole diameter?",
    "PREDICTION  Diameter scales smoothly with A-tubule protofilament count — 169 nm at 9, 255 nm at 13, 299 nm at 18 — with no steric barrier anywhere in that range.",
    "TEST  Compare species or mutants with non-canonical protofilament numbers.  FALSIFIED IF  diameter is insensitive to protofilament count, or clashes appear where the model predicts none.", BLUE],
   ["Q5", "Does a longer SAS-6 coiled coil widen the centriole proportionally?",
    "PREDICTION  No — the response saturates. 40→45 nm adds 4.3 nm of diameter; 45→50 nm adds only 2.4 nm, while strain rises fivefold. Length buys progressively less width and progressively more strain.",
    "TEST  Engineered SAS-6 coiled-coil length variants; measure diameter.  FALSIFIED IF  diameter tracks coiled-coil length linearly.", PURPLE]
  ].forEach(([n, q, pr, te, c], i) => {
    const y = 1.5 + i * 1.82;
    card(s, 0.7, y, 11.9, 1.62);
    dot(s, 1.0, y + 0.2, n, c, 0.5);
    s.addText(q, t("", { x: 1.7, y: y + 0.16, w: 10.6, h: 0.4, fontSize: 16.5, bold: true, color: c, fontFace: HEAD, margin: 0 }));
    s.addText(pr, t("", { x: 1.7, y: y + 0.6, w: 10.6, h: 0.55, fontSize: 12.5, margin: 0 }));
    s.addText(te, t("", { x: 1.7, y: y + 1.15, w: 10.6, h: 0.4, fontSize: 12.5, color: MUTE, margin: 0 }));
  });
  s.addText("Q1 is the strongest: the prediction is counter-intuitive, quantitative, and distinguishable by a single diameter measurement.",
    t("", { x: 0.7, y: 6.95, w: 11.9, h: 0.4, fontSize: 13.5, italic: true, color: MUTE }));

  // --- closing
  s = p.addSlide();
  s.background = { color: DARK };
  s.addText("Recommended order of work", { x: 0.9, y: 0.75, w: 11.5, h: 0.8, fontSize: 34,
    bold: true, color: WHITE, fontFace: HEAD });
  [["Now", "Run Q1 as a hold-out test. It needs no model changes and would be the first non-circular validation.", CYAN],
   ["Next", "Add uncertainty propagation, then calibrate the joint bands against whatever angular variance the data can give.", "9C8FD4"],
   ["Then", "Bond energies in kT, which unlocks force and fluctuation claims and revives the parked soft-mode analysis.", "9C8FD4"],
   ["Later", "The 2.5D stack, once there is a longitudinal question worth the cost of answering.", "8A82A8"]
  ].forEach(([k, b, c], i) => {
    const y = 1.85 + i * 1.2;
    s.addShape("roundRect", { x: 0.9, y, w: 11.5, h: 1.02, fill: { color: "2A1B4D" },
      line: { color: "2A1B4D" }, rectRadius: 0.06 });
    s.addText(k, t("", { x: 1.2, y: y + 0.28, w: 1.5, h: 0.45, fontSize: 18, bold: true, color: c, fontFace: HEAD, margin: 0 }));
    s.addText(b, t("", { x: 2.8, y: y + 0.2, w: 9.3, h: 0.7, fontSize: 14, color: "C9C4DB", margin: 0 }));
  });
  s.addText("The mechanics are further ahead than the evidence. Closing that gap is worth more than any new feature.",
    t("", { x: 0.9, y: 6.75, w: 11.5, h: 0.5, fontSize: 15, italic: true, color: CYAN }));

  return p.writeFile({ fileName: P + "03_CenGeometry_Expert_Recommendations.pptx" });
}

deckOne()
  .then(() => deckTwo())
  .then(() => deckThree())
  .then(() => console.log("all three decks written"))
  .catch(e => { console.error("FAILED:", e); process.exit(1); });

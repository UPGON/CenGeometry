/* Layout check for the decks.
 *
 * Everything here is absolutely positioned, so the two ways a slide breaks are
 * an element running off the canvas and two elements sitting on top of each
 * other. Neither shows up in the build -- pptxgenjs writes the file happily
 * either way -- and there is no LibreOffice on this machine to render and
 * eyeball. So intercept every placement call and check the geometry directly.
 *
 *     node check_layout.js
 *
 * Overlap is reported only between things that actually occlude: text over
 * text, text over an image or table. Cards and shapes are deliberately drawn
 * underneath text, so those pairs are ignored.
 */
const pptxgen = require("pptxgenjs");
const fs = require("fs");

const W = 13.333, H = 7.5;      // LAYOUT_WIDE, inches
const EPS = 0.02;               // ignore hairline contact
const rec = [];
let slideNo = 0;

const origAddSlide = pptxgen.prototype.addSlide;
pptxgen.prototype.addSlide = function (...a) {
  const s = origAddSlide.apply(this, a);
  const no = ++slideNo;
  for (const m of ["addText", "addImage", "addShape", "addTable", "addChart"]) {
    const orig = s[m].bind(s);
    s[m] = (...args) => {
      const o = args[args.length - 1] || {};
      let h = +o.h;
      // addTable is usually given rowH and no height at all, so it slipped
      // through both checks unmeasured. Reconstruct it from the row count.
      if (m === "addTable" && !Number.isFinite(h) && Array.isArray(args[0]))
        h = args[0].length * (o.rowH || 0.4);
      let label = "";
      if (typeof args[0] === "string") label = args[0].replace(/\s+/g, " ").slice(0, 46);
      else if (m === "addImage" && o.path) label = o.path.split("/").pop();
      else if (m === "addShape") label = String(args[0]);
      rec.push({ slide: no, kind: m, x: +o.x, y: +o.y, w: +o.w, h, label });
      return orig(...args);
    };
  }
  return s;
};

const OCCLUDING = new Set(["addText", "addImage", "addTable", "addChart"]);

function overlap(a, b) {
  const ox = Math.min(a.x + a.w, b.x + b.w) - Math.max(a.x, b.x);
  const oy = Math.min(a.y + a.h, b.y + b.h) - Math.max(a.y, b.y);
  return ox > EPS && oy > EPS ? ox * oy : 0;
}

process.on("exit", () => {
  fs.writeFileSync(__dirname + "/figs/_layout.json", JSON.stringify(rec));
  const bad = [];
  const overflow = rec.filter(r =>
    Number.isFinite(r.x) && Number.isFinite(r.w) &&
    (r.x < -EPS || r.y < -EPS || r.x + r.w > W + EPS || r.y + r.h > H + EPS));

  const bySlide = {};
  rec.forEach(r => { (bySlide[r.slide] = bySlide[r.slide] || []).push(r); });
  const collisions = [];
  Object.entries(bySlide).forEach(([no, items]) => {
    const occ = items.filter(r => OCCLUDING.has(r.kind) &&
      Number.isFinite(r.x) && Number.isFinite(r.w) && Number.isFinite(r.h));
    for (let i = 0; i < occ.length; i++)
      for (let j = i + 1; j < occ.length; j++) {
        const a = occ[i], b = occ[j], area = overlap(a, b);
        // text boxes are routinely taller than their content; only flag a
        // genuine collision, i.e. more than a third of the smaller element
        const small = Math.min(a.w * a.h, b.w * b.h);
        if (area > 0.34 * small)
          collisions.push({ slide: +no, a, b, area: area.toFixed(2) });
      }
  });

  console.log(`\n${slideNo} slides, ${rec.length} placed elements`);
  if (overflow.length) {
    console.log(`\nOFF-CANVAS (${overflow.length}):`);
    overflow.forEach(r => console.log(
      `  slide ${r.slide}  ${r.kind}  x=${r.x} y=${r.y} w=${r.w} h=${r.h}` +
      `  -> right ${(r.x + r.w).toFixed(2)} bottom ${(r.y + r.h).toFixed(2)}   "${r.label}"`));
    bad.push("overflow");
  } else console.log("  no element runs off the canvas");

  if (collisions.length) {
    console.log(`\nOVERLAPS (${collisions.length}):`);
    collisions.forEach(c => console.log(
      `  slide ${c.slide}  ${c.area} in²  "${c.a.label}" (${c.a.kind}) x "${c.b.label}" (${c.b.kind})`));
    bad.push("overlap");
  } else console.log("  no occluding elements collide");

  if (bad.length) process.exitCode = 1;
});

require("./build_decks.js");

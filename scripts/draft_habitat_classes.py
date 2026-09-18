#!/usr/bin/env python3
"""Draft environment classes for habitat values of phenotype_data.tsv that config/habitat_classes.tsv does not list yet.

Prints rows in the table's format (habitat, class, class_label, n_organisms, basis) for review; add the reviewed rows
to config/habitat_classes.tsv. Keywords assign classes; a value matching several classes gets a combined class when the
combination is a known one, otherwise basis "unresolved" and class "other" so that it is decided by hand.

    python3 scripts/draft_habitat_classes.py [--phenotype data/electron/phenotype_data.tsv] [--table config/habitat_classes.tsv]
"""

from __future__ import annotations

import argparse
import collections
import csv
import re
import sys
from pathlib import Path
from typing import List, Tuple

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ko_detector.phenotype import clean_habitat, is_shifted  # noqa: E402

# (class, pattern) in display order; the patterns of the reviewed 2026-09-18 draft
CLASSES = [
    ("vent", r"hydrothermal|\bvents?\b|chimney|black smoker|\bsmokers?\b|pacific rise|diffuse-flow|alvinella|tube ?worm"),
    ("hotspring", r"hot springs?|hot pools?|geotherm|solfatar|fumarol|volcan|caldera|\bthermal\b|mud pots?|steamers?|hot water|"
                  r"boiling|hot (mud|sediment)|obsidian pool|yellowstone"),
    ("marine", r"marine|sea ?water|\bseas?\b|ocean|coast|estuar|tidal|salt marsh|mangrove|beach|shore|\bbay\b|sea ?floor|deep-sea|"
               r"submarine|\bgulf\b|lagoon|fjord|\binlet\b|pelagic|seaweed|brackish|harbou?r"),
    ("saline", r"\bsoda\b|hypersaline|saline|\bsalt\b|brine|salar|saltern|alkaline lake|halite|salty"),
    ("freshwater", r"fresh ?water|\blakes?\b|rivers?|ponds?|\bstreams?\b|groundwater|ground water|aquifer|\bwells?\b|creek|wetland|"
                   r"\bmarsh|swamp|\bbogs?\b|\bpeat|\bfen\b|springs?\b|drinking water|tap water|chemocline"),
    ("subsurface", r"\boil\b|petroleum|subsurface|borehole|gas field|oil ?field|oil reservoir|deep (granitic|biosphere)|deep aquifer|"
                   r"drill|underground|crude"),
    ("mine", r"\bmines?\b|mining|acid mine|pyrit|\bores?\b|tailings?|\bcoal\b|bioleach|heap leach|gold|mineral heap"),
    ("engineered", r"sewage|sludge|wastewater|waste ?water|bioreactor|reactor|digest[eo]r|treatment|biofilter|landfill|compost|sewer|"
                   r"factory|industrial|concrete|pipe|tyre|biogas|enrichment|\bculture\b|wwtp|fuel cell|silage|sour mash|piggery|"
                   r"waste\b|aquaculture|aquarium|buffer"),
    ("soil", r"soil|rhizosphere|paddy|farmland|agricultur|rice field|\broots?\b|forest|straw|\bseeds?\b|plants? (roots?|tissue|material)|"
             r"medicago|wetwood|tree|permafrost|tundra"),
    ("host", r"rumen|\bgut\b|hindgut|intestin|feces|faeces|termite|animal|insect|human|clinical|blood|\boral\b|manure|\bgills?\b|"
             r"symbio|\bworms?\b|mussel|\bclams?\b|sponge|snail|shrimp|\bfish\b|cattle|sheep|\bpigs?\b|coral|wombat|glider|"
             r"millipede|synovial|bronch|patholog|tuberculosis|\bmilk\b|chicken|poultry"),
]
FALLBACK = [("sediment", r"sediment|\bmud\b|\bsand\b|\bmats?\b|clay|rock"), ("water", r"\bwater\b")]
COMBINED = {
    "hotspring+soil": "geothermal_soil",
    "marine+saline": "salt_marsh",
    "engineered+host": "livestock_waste",
    "freshwater+engineered": "polluted_freshwater",
    "marine+subsurface": "marine_subsurface",
}
TERRESTRIAL = r"terrestrial|springs?\b|\bponds?\b|\bpools?\b|caldera|yellowstone|\bland\b"
EXPLICIT_FRESH = r"fresh ?water|\brivers?\b"
ORDER = [key for key, _ in CLASSES] + [key for key, _ in FALLBACK]


def classify(value: str) -> List[str]:
    low = value.lower()
    hits = [key for key, pattern in CLASSES if re.search(pattern, low)]
    if "vent" in hits and re.search(TERRESTRIAL, low) and "marine" not in hits:
        hits.remove("vent")                                               # terrestrial hydrothermal -> hot spring
        if "hotspring" not in hits:
            hits.append("hotspring")
    if "vent" in hits:
        hits = [h for h in hits if h not in ("marine", "hotspring")]      # marine is implied for vents
    if "hotspring" in hits and "marine" in hits:
        hits = [h for h in hits if h not in ("hotspring", "marine")] + ["vent"]   # marine geothermal
    if "freshwater" in hits and any(k in hits for k in ("marine", "saline", "hotspring", "vent", "subsurface")) \
            and not re.search(EXPLICIT_FRESH, low):
        hits.remove("freshwater")                                         # "lake", "spring" taken by the stronger class
    if not hits:
        hits = [key for key, pattern in FALLBACK if re.search(pattern, low)][:1]
    return sorted(set(hits), key=ORDER.index) or ["other"]


def propose(keys: List[str]) -> Tuple[str, str]:
    if len(keys) == 1:
        return keys[0], "unclassified" if keys[0] == "other" else "single"
    combo = "+".join(keys)
    if combo in COMBINED:
        return COMBINED[combo], "combined"
    return "other", f"unresolved ({combo})"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--phenotype", type=Path, default=ROOT / "data" / "electron" / "phenotype_data.tsv")
    parser.add_argument("--table", type=Path, default=ROOT / "config" / "habitat_classes.tsv")
    args = parser.parse_args()

    with args.table.open(encoding="utf-8") as handle:
        known = {row["habitat"] for row in csv.DictReader((l for l in handle if not l.startswith("#")), delimiter="\t")}
    with args.phenotype.open(encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if not is_shifted(row)]
    values = collections.Counter(clean_habitat(row.get("habitat", "")) for row in rows)
    new = {value: n for value, n in values.items() if value and value not in known}
    if not new:
        print(f"every habitat value is in {args.table}", file=sys.stderr)
        return 0
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(["habitat", "class", "class_label", "n_organisms", "basis"])
    for value, n in sorted(new.items(), key=lambda kv: kv[0].lower()):
        klass, basis = propose(classify(value))
        writer.writerow([value, klass, "", n, basis])
    print(f"{len(new)} new habitat value(s); review the class column before adding the rows to {args.table}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

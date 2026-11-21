"""Space-Track fetches for the CZ-8A Y13 model (4 live queries; cached by st.py).

data/hist_2026-243.json   gp_history of CZ-8A Y12 (2026-243) payloads + stage, launch -> now
data/gp_2025-295.json     current gp of the target plane group 2025-295 (CZ-12 Y4)
data/gp_2026-176_221.json current gp of the Y10/Y11 groups and their partners (for plane checks)
Inputs already on disk and reused, not refetched: ../launch_traj_12/data/hist_joiners.json (Y10/Y11 history),
../launch_traj_12/data/gp_hulianwang_current.json (all Guowang planes, 2026-10-25).
"""
import json
import os
import sys

sys.path.insert(0, os.path.expanduser("~/.claude/skills/space-track"))
from st import query  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

sat = query(["class", "satcat", "INTLDES", "~~2026-243", "orderby", "NORAD_CAT_ID asc"])
json.dump(sat, open(f"{OUT}/satcat_2026-243.json", "w"))
ids = [int(r["NORAD_CAT_ID"]) for r in sat]
print("2026-243 objects:", ids)
rows = query(["class", "gp_history", "NORAD_CAT_ID", ",".join(map(str, ids)),
              "EPOCH", "2026-10-25--2026-11-21", "orderby", "NORAD_CAT_ID asc,EPOCH asc"])
json.dump(rows, open(f"{OUT}/hist_2026-243.json", "w"))
print(len(rows), "history rows")

g295 = query(["class", "gp", "NORAD_CAT_ID", ",".join(str(n) for n in range(67059, 67069))], ttl_h=1.0)
json.dump(g295, open(f"{OUT}/gp_2025-295.json", "w"))
print(len(g295), "2025-295 gp rows")

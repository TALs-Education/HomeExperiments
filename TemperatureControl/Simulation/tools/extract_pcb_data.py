"""
Extract heated-PCB-beam geometry from the Fusion/EAGLE gerber zips and the
experiment logs, and write them as JS data files for pcb_heater_sim.html.

    python Simulation/tools/extract_pcb_data.py

Outputs (next to the HTML so it works from file://):
    Simulation/pcb_geometry.js     window.PCB_GEOMETRY
    Simulation/experiment_data.js  window.EXPERIMENT_DATA

All coordinates are board millimetres, origin at the gerber origin
(Y runs along the beam, heater 1 at low Y).
"""
import csv
import json
import math
import os
import re
import zipfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
OUT_DIR = os.path.join(ROOT, 'Simulation')

BOARDS = {
    # V01: heater copper area at both ends of the beam
    'dual': os.path.join(ROOT, 'PCB', 'Heating Element', 'V01', 'TemperatureControlHeater_PCB_V01.zip'),
    # V02: single heater, beam extends to a free end
    'single': os.path.join(ROOT, 'PCB', 'Heating Element', 'V02', 'TemperatureControlHeater_PCB_V02.zip'),
}

EXPERIMENTS = [
    # (id, file, board, resistor ohms (each of 2 in series), description)
    ('single_10ohm', 'Setup2_10Ohm.txt', 'single', 10.0,
     'Single heater, 2 x 10 ohm, 100% step then cool-down'),
    ('dual_15ohm', 'Setup1.txt', 'dual', 15.0,
     'Dual heater, 2 x 15 ohm per heater: H1, H2, then both'),
]

R3 = lambda v: round(v, 4)


# ---------------------------------------------------------------- gerber parsing
def parse_gerber(text):
    """Return ordered primitives: ('R', pol, [contour...]) | ('L', pol, x0,y0,x1,y1,w) | ('F', pol, x,y,shape,w,h)."""
    m = re.search(r'%FSLAX(\d)(\d)Y\d\d\*%', text)
    dec = int(m.group(2)) if m else 4
    scale = 10.0 ** -dec
    if '%MOIN*%' in text:
        scale *= 25.4
    aps = {}
    for a in re.finditer(r'%ADD(\d+)([CRO]),([\d.]+)(?:X([\d.]+))?\*%', text):
        w = float(a.group(3)); h = float(a.group(4) or w)
        if '%MOIN*%' in text:
            w *= 25.4; h *= 25.4
        aps[int(a.group(1))] = (a.group(2), w, h)

    prims = []
    pol = 'D'
    x = y = 0.0
    ap = None
    interp = 1  # 1 linear, 2 cw, 3 ccw
    region = None
    contour = None

    def arc_points(x0, y0, x1, y1, i, j, cw):
        cx, cy = x0 + i, y0 + j
        r = math.hypot(i, j)
        a0 = math.atan2(y0 - cy, x0 - cx)
        a1 = math.atan2(y1 - cy, x1 - cx)
        if cw:
            if a1 >= a0: a1 -= 2 * math.pi
        else:
            if a1 <= a0: a1 += 2 * math.pi
        n = max(2, int(abs(a1 - a0) / (math.pi / 32)))
        return [(cx + r * math.cos(a0 + (a1 - a0) * k / n), cy + r * math.sin(a0 + (a1 - a0) * k / n)) for k in range(1, n + 1)]

    for tok in re.split(r'(%[^%]*%)', text):
        if tok.startswith('%'):
            if tok.startswith('%LP'):
                pol = tok[3]
            continue
        for c in tok.split('*'):
            c = c.strip()
            if not c or c.startswith('G04'):
                continue
            if c == 'G36':
                region = []; contour = None; continue
            if c == 'G37':
                if contour and len(contour) > 2: region.append(contour)
                prims.append(('R', pol, region)); region = None; contour = None; continue
            g = re.match(r'^G0?([123])(.*)$', c)
            if g:
                interp = int(g.group(1)); c = g.group(2)
                if not c: continue
            d = re.match(r'^D(\d+)$', c)
            if d and int(d.group(1)) >= 10:
                ap = int(d.group(1)); continue
            m = re.match(r'^(?:X(-?\d+))?(?:Y(-?\d+))?(?:I(-?\d+))?(?:J(-?\d+))?D0?([123])$', c)
            if not m:
                continue
            nx = int(m.group(1)) * scale if m.group(1) else x
            ny = int(m.group(2)) * scale if m.group(2) else y
            i = int(m.group(3)) * scale if m.group(3) else 0.0
            j = int(m.group(4)) * scale if m.group(4) else 0.0
            op = m.group(5)
            if op == '1':
                pts = [(nx, ny)] if interp == 1 else arc_points(x, y, nx, ny, i, j, interp == 2)
                if region is not None:
                    if contour is None: contour = [(x, y)]
                    contour.extend(pts)
                else:
                    a = aps.get(ap, ('C', 0.0, 0.0))
                    px, py = x, y
                    for qx, qy in pts:
                        prims.append(('L', pol, px, py, qx, qy, a[1])); px, py = qx, qy
            elif op == '2':
                if region is not None:
                    if contour and len(contour) > 2: region.append(contour)
                    contour = [(nx, ny)]
            elif op == '3':
                a = aps.get(ap, ('C', 0.0, 0.0))
                prims.append(('F', pol, nx, ny, a[0], a[1], a[2]))
            x, y = nx, ny
    return prims


def parse_drill(text):
    tools, cur, holes = {}, None, []
    div = 1000.0
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('METRIC'):
            fmt = re.search(r'(\d+)\.(\d+)', line)
            if fmt: div = 10.0 ** len(fmt.group(2))
        t = re.match(r'^T(\d+)C([\d.]+)$', line)
        if t: tools[t.group(1)] = float(t.group(2)); continue
        t = re.match(r'^T(\d+)$', line)
        if t: cur = t.group(1); continue
        p = re.match(r'^X(-?\d+)Y(-?\d+)$', line)
        if p: holes.append([R3(int(p.group(1)) / div), R3(int(p.group(2)) / div), tools[cur]])
    return holes


def chain_loops(lines):
    """Chain consecutive line primitives into closed loops."""
    loops, cur = [], None
    for _, _, x0, y0, x1, y1, _ in lines:
        if cur and abs(cur[-1][0] - x0) < 1e-3 and abs(cur[-1][1] - y0) < 1e-3:
            cur.append((x1, y1))
        else:
            if cur: loops.append(cur)
            cur = [(x0, y0), (x1, y1)]
    if cur: loops.append(cur)
    return [l for l in loops if len(l) > 3 and math.hypot(l[0][0] - l[-1][0], l[0][1] - l[-1][1]) < 1e-3]


def bbox(pts):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    return [R3(min(xs)), R3(min(ys)), R3(max(xs)), R3(max(ys))]


def cluster(points, dist):
    groups = []
    for p in points:
        for g in groups:
            if any(math.hypot(p[0] - q[0], p[1] - q[1]) < dist for q in g):
                g.append(p); break
        else:
            groups.append([p])
    return groups


def pack_layer(prims):
    out = []
    for p in prims:
        if p[0] == 'R':
            out.append(['R', p[1], [[[R3(a), R3(b)] for a, b in c] for c in p[2]]])
        elif p[0] == 'L':
            out.append(['L', p[1]] + [R3(v) for v in p[2:7]])
        else:
            out.append(['F', p[1], R3(p[2]), R3(p[3]), p[4], R3(p[5]), R3(p[6])])
    return out


def extract_board(zip_path):
    z = zipfile.ZipFile(zip_path)
    files = {os.path.splitext(n)[1].upper().lstrip('.'): n for n in z.namelist()}
    read = lambda ext: z.read(files[ext]).decode('latin-1')

    gtl, gbl = parse_gerber(read('GTL')), parse_gerber(read('GBL'))
    gtp, gbp = parse_gerber(read('GTP')), parse_gerber(read('GBP'))
    gml = parse_gerber(read('GML'))
    gto = parse_gerber(read('GTO'))
    gts, gbs = parse_gerber(read('GTS')), parse_gerber(read('GBS'))

    # Outline: zero-width contour is the board edge; closed loops of the
    # wider milling contour (excluding the board edge copy) are cut-outs.
    edge = chain_loops([p for p in gml if p[0] == 'L' and p[6] == 0.0])
    board = max(edge, key=lambda l: (bbox(l)[2] - bbox(l)[0]) * (bbox(l)[3] - bbox(l)[1]))
    bb = bbox(board)
    cutouts = []
    widths = sorted({p[6] for p in gml if p[0] == 'L' and p[6] > 0})
    for w in widths:
        for l in chain_loops([p for p in gml if p[0] == 'L' and p[6] == w]):
            b = bbox(l)
            if abs(b[0] - bb[0]) < 0.01 and abs(b[2] - bb[2]) < 0.01:
                continue  # board edge duplicate
            cutouts.append(l)
        if cutouts:
            break

    # Heater resistors: bottom paste pads in pairs along Y (2512 footprint).
    rpads = [(p[2], p[3], p[5], p[6]) for p in gbp if p[0] == 'F']
    resistors = []
    used = set()
    for i, a in enumerate(rpads):
        if i in used: continue
        for j, b in enumerate(rpads):
            if j <= i or j in used: continue
            if abs(a[0] - b[0]) < 0.01 and 3.0 < abs(a[1] - b[1]) < 8.0:
                used |= {i, j}
                y0, y1 = sorted([a[1], b[1]])
                resistors.append({
                    'center': [R3(a[0]), R3((y0 + y1) / 2)],
                    'pads': [[R3(a[0]), R3(y0)], [R3(b[0]), R3(y1)]],
                    'pad_size': [a[2], a[3]],
                    # body + terminations footprint (heat input area)
                    'rect': [R3(a[0] - a[2] / 2), R3(y0 - a[3] / 2), R3(a[0] + a[2] / 2), R3(y1 + a[3] / 2)],
                })
    heaters = []
    for g in cluster([r['center'] for r in resistors], 8.0):
        rs = [r for r in resistors if r['center'] in g]
        heaters.append({'center': [R3(sum(r['center'][0] for r in rs) / len(rs)),
                                   R3(sum(r['center'][1] for r in rs) / len(rs))],
                        'resistors': rs})
    heaters.sort(key=lambda h: h['center'][1])

    # Temperature sensors: SOT-23 (3 pad) clusters of the 1.32 x 0.61 top paste pads.
    spads = [(p[2], p[3]) for p in gtp if p[0] == 'F' and abs(p[5] - 1.3208) < 0.01]
    sensors = []
    for g in cluster(spads, 2.5):
        b = bbox(g)
        sensors.append({'center': [R3((b[0] + b[2]) / 2), R3((b[1] + b[3]) / 2)], 'pads': [[R3(a), R3(c)] for a, c in g]})
    sensors.sort(key=lambda s: s['center'][1])

    # Heater copper area and beam strip (largest top-copper pours).
    pours = [bbox([pt for c in p[2] for pt in c]) for p in gtl if p[0] == 'R' and p[1] == 'D']

    return {
        'board': [[R3(a), R3(b)] for a, b in board],
        'bbox': bb,
        'cutouts': [[[R3(a), R3(b)] for a, b in l] for l in cutouts],
        'top': pack_layer(gtl),
        'bottom': pack_layer(gbl),
        'silk_top': pack_layer([p for p in gto if p[0] == 'L']),
        # solder-mask openings (exposed metal); everything else is mask-covered
        'mask_top': pack_layer(gts),
        'mask_bottom': pack_layer(gbs),
        'drills': parse_drill(read('TXT')),
        'pours': pours,
        'heaters': heaters,
        'sensors': sensors,
        'source': os.path.relpath(zip_path, ROOT).replace('\\', '/'),
    }


def extract_experiment(path):
    rows = list(csv.DictReader(open(path)))
    t0 = float(rows[0]['t_sec'])
    cols = ['ambient_c', 'beam1_c', 'beam2_c', 'beam3_c', 'beam4_c', 'beam5_c', 'heater1_pwm', 'heater2_pwm', 'fan_pwm']
    data = [[round(float(r['t_sec']) - t0, 3)] + [float(r[c]) for c in cols] for r in rows]
    return {'columns': ['t'] + cols, 'rows': data}


def main():
    geo = {k: extract_board(v) for k, v in BOARDS.items()}
    for k, g in geo.items():
        print(f"{k}: bbox={g['bbox']} cutouts={len(g['cutouts'])} heaters={[h['center'] for h in g['heaters']]} "
              f"sensors={[s['center'] for s in g['sensors']]} pours={g['pours']}")
    with open(os.path.join(OUT_DIR, 'pcb_geometry.js'), 'w') as f:
        f.write('// Generated by tools/extract_pcb_data.py from the heater PCB gerbers. Do not edit.\n')
        f.write('window.PCB_GEOMETRY = ' + json.dumps(geo, separators=(',', ':')) + ';\n')

    exps = {}
    for eid, fname, board, ohm, desc in EXPERIMENTS:
        e = extract_experiment(os.path.join(ROOT, 'Experiments Results', fname))
        e.update({'board': board, 'resistor_ohm': ohm, 'description': desc, 'file': fname})
        exps[eid] = e
        print(f"{eid}: {len(e['rows'])} rows, {e['rows'][-1][0]:.0f} s")
    with open(os.path.join(OUT_DIR, 'experiment_data.js'), 'w') as f:
        f.write('// Generated by tools/extract_pcb_data.py from Experiments Results logs. Do not edit.\n')
        f.write('window.EXPERIMENT_DATA = ' + json.dumps(exps, separators=(',', ':')) + ';\n')


if __name__ == '__main__':
    main()

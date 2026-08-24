import glob, os

CASES = {
    'v1': 'v_1_-_Haffner',
    'v1b': 'v_1b_-_Haffner_con_flujo_termico',
    'v2': 'v_2_-_DECAY__2007_',
    'v3': 'v_3_-_DECAY__2025_',
    'v4': 'v_4_-_DECAY-XSECTIONS__2025_',
}


def parse(path):
    """{'atoms_irr','atoms_cool','act_irr','act_cool'} -> {nuc: {t: val}}, tablas BY ZONE."""
    txt = open(path, 'rb').read().decode('latin-1').replace('\r', '')
    lines = txt.split('\n')
    out = {k: {} for k in ('atoms_irr', 'atoms_cool', 'act_irr', 'act_cool')}
    for i, L in enumerate(lines):
        if 'NUMBER OF ATOMS' in L:
            mag = 'atoms'
        elif 'NUCLIDE RADIOACTIVITY' in L:
            mag = 'act'
        else:
            continue
        ctx = '\n'.join(lines[max(0, i - 8):i])
        if 'BY ZONE' not in ctx:
            continue
        fase = 'irr' if 'DURING IRRADIATION' in ctx else 'cool'
        key = f'{mag}_{fase}'
        j = i + 1
        while j < min(i + 8, len(lines)) and 'INITIAL' not in lines[j]:
            j += 1
        if j >= len(lines) or 'INITIAL' not in lines[j]:
            continue
        hdr = lines[j].split()
        times = [h if h in ('SHUTDOWN', 'RESTART', 'INITIAL') else float(h) for h in hdr]
        # C.15: en enfriamiento, el primer RESTART es el apagado real; los de los
        # conjuntos siguientes duplican el ultimo punto del anterior.
        if fase == 'cool' and 'RESTART' in times:
            primero = not any('SHUTDOWN' in d for d in out[key].values())
            times = ['SHUTDOWN' if (t == 'RESTART' and primero) else
                     ('DUP' if t == 'RESTART' else t) for t in times]
        j += 1
        while j < len(lines):
            r = lines[j]
            if not r.strip() or r[0] in '01':
                break
            nuc, vals = r[:7].strip(), r[7:].split()
            if nuc and len(vals) == len(times):
                d = out[key].setdefault(nuc, {})
                for t, v in zip(times, vals):
                    if t == 'DUP':
                        continue
                    d[t] = float(v)
            j += 1
    return out


def load_all(root='/home/claude/sim'):
    res = {}
    for k, d in CASES.items():
        f = glob.glob(os.path.join(root, d, '*', 'fort.6'))[0]
        res[k] = parse(f)
    return res


def cool_times(d):
    return sorted([t for t in d if isinstance(t, float)])


if __name__ == '__main__':
    R = load_all()
    for k, r in R.items():
        a = r['act_cool']['I131']
        print(k, 'nucleidos:', len(r['act_cool']),
              '| instantes:', len(cool_times(a)),
              '| SHUTDOWN:', a.get('SHUTDOWN'))
        print('   t:', cool_times(a))

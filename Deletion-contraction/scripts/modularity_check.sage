load("qarr.sage")

def build_flats_with_content(Delta, n, q):
    """Same construction as fast_characteristic_polynomial, but returns the
    flats dict (key -> (rank, content_bitmask)) instead of the polynomial."""
    normals = generate_normals(Delta, n, q)
    m = len(normals)
    F = GF(q)
    G = matrix(F, normals)
    bottom_key = ()
    flats = {bottom_key: (0, 0)}
    frontier = [bottom_key]
    while frontier:
        next_frontier = []
        for key in frontier:
            rank, content = flats[key]
            if rank == n:
                continue
            B = matrix(F, list(key)) if key else matrix(F, 0, n)
            R = G.__copy__()
            for row in B.rows():
                piv = next(i for i in range(n) if row[i] != 0)
                col = R.column(piv)
                if not col.is_zero():
                    R = R - col.column() * row.row()
            groups = {}
            for idx in range(m):
                if (content >> idx) & 1:
                    continue
                r = R.row(idx)
                if r.is_zero():
                    continue
                piv = next(i for i in range(n) if r[i] != 0)
                inv = r[piv] ** (-1)
                canon = tuple((inv * r[i]) for i in range(n))
                groups.setdefault(canon, []).append(idx)
            for canon, idxs in groups.items():
                rep = normals[idxs[0]]
                new_key = _echelon_key(list(key) + [rep], q, n)
                new_content = content
                for idx in idxs:
                    new_content |= (1 << idx)
                if new_key not in flats:
                    flats[new_key] = (rank + 1, new_content)
                    next_frontier.append(new_key)
        frontier = next_frontier
    return flats, normals

def check_modularity(Delta, F_face, n, q):
    flats, normals = build_flats_with_content(Delta, n, q)
    m = len(normals)
    F_face = frozenset(F_face)
    # content of X_F: all ground elements (hyperplanes) whose support is a subset of F_face
    XF_content = 0
    for idx, v in enumerate(normals):
        support = frozenset(i for i in range(n) if v[i] != 0)
        if support <= F_face:
            XF_content |= (1 << idx)
    # find XF's rank among enumerated flats
    rank_of = {}
    for key, (rank, content) in flats.items():
        rank_of[content] = rank
    if XF_content not in rank_of:
        return None, "X_F content not found as an enumerated flat (unexpected)"
    rk_XF = rank_of[XF_content]

    failures = []
    all_contents = list(rank_of.keys())
    for Ycontent in all_contents:
        rk_Y = rank_of[Ycontent]
        meet_content = XF_content & Ycontent
        if meet_content not in rank_of:
            failures.append(("meet not a flat?!", Ycontent))
            continue
        rk_meet = rank_of[meet_content]
        # join = smallest enumerated flat whose content contains XF_content | Ycontent
        union_c = XF_content | Ycontent
        candidates = [c for c in all_contents if (c & union_c) == union_c]
        rk_join = min(rank_of[c] for c in candidates)
        lhs = rk_XF + rk_Y
        rhs = rk_join + rk_meet
        if lhs != rhs:
            failures.append((Ycontent, rk_Y, rk_join, rk_meet, lhs, rhs))
    return (len(failures) == 0), (rk_XF, len(all_contents), len(failures))

# Test 1: book_3pages -- where the isolating-edge theorem HOLDS
Delta1 = complex_from_maximal([(0,1,2),(0,1,3),(0,1,4)])
ok1, info1 = check_modularity(Delta1, {0,1,2}, 5, 2)
print("book_3pages, X_F for F={0,1,2}: modular?", ok1, info1)

# Test 2: prism (sides open) -- where NO edge works at all
Delta2 = complex_from_maximal([(0,1,2),(3,4,5),(0,3),(1,4),(2,5)])
ok2, info2 = check_modularity(Delta2, {0,1,2}, 6, 2)
print("prism, X_F for F={0,1,2}: modular?", ok2, info2)

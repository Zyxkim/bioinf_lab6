import itertools
from typing import Callable
from tqdm import tqdm


class Fasta:
    def __init__(self):
        self.name = ""
        self.sequence = ""


def read_fasta(file):
    name, seq = None, []
    for line in file:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield name, ''.join(seq)


def parse(file_name):
    res = []
    with open(file_name) as f:
        for name, seq in read_fasta(f):
            record = Fasta()
            record.name = name
            record.sequence = seq
            res.append(record)
    return res


def score_fun(a: str,
              b: str,
              match_score: int = 5,
              mismatch_score: int = -4) -> int:
    return match_score if a == b else mismatch_score


def SW(seq1: str, seq2: str, score_fun: Callable = score_fun, gap_score: int = -5):
    m, n = len(seq1) + 1, len(seq2) + 1

    matrix = [[0] * n for _ in range(m)]

    for i in range(m):
        matrix[i][0] = i * gap_score
    for j in range(n):
        matrix[0][j] = j * gap_score

    for i, j in itertools.product(range(1, m), range(1, n)):
        matrix[i][j] = max(matrix[i - 1][j - 1] + score_fun(seq1[i - 1], seq2[j - 1]),
                           matrix[i - 1][j] + gap_score,
                           matrix[i][j - 1] + gap_score)
    score = gap_score * max(n, m)
    i, j = -1, -1
    for k in range(1, n):
        for l in range(1, m):
            if matrix[k][l] > score:
                score = matrix[k][l]
                i, j = k, l
    aln1 = ""
    aln2 = ""
    while i > 0 or j > 0:
        a, b = '-', '-'
        # (A, B)
        if i > 0 and j > 0 and matrix[i][j] == matrix[i - 1][j - 1] + score_fun(seq1[i - 1], seq2[j - 1]):
            a = seq1[i - 1]
            b = seq2[j - 1]
            i -= 1
            j -= 1

        # (A, -)
        elif i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_score:
            a = seq1[i - 1]
            i -= 1

        # (-, A)
        elif j > 0 and matrix[i][j] == matrix[i][j - 1] + gap_score:
            b = seq2[j - 1]
            j -= 1

        aln1 += a
        aln2 += b
    return aln1[::-1], aln2[::-1], score


def get_blosum62(q, t):
    blosum62 = {
        'WF': 1, 'LR': -2, 'SP': -1, 'VT': 0,
        'QQ': 5, 'NA': -2, 'ZY': -2, 'WR': -3,
        'QA': -1, 'SD': 0, 'HH': 8, 'SH': -1,
        'HD': -1, 'LN': -3, 'WA': -3, 'YM': -1,
        'GR': -2, 'YI': -1, 'YE': -2, 'BY': -3,
        'YA': -2, 'VD': -3, 'BS': 0, 'YY': 7,
        'GN': 0, 'EC': -4, 'YQ': -1, 'ZZ': 4,
        'VA': 0, 'CC': 9, 'MR': -1, 'VE': -2,
        'TN': 0, 'PP': 7, 'VI': 3, 'VS': -2,
        'ZP': -1, 'VM': 1, 'TF': -2, 'VQ': -2,
        'KK': 5, 'PD': -1, 'IH': -3, 'ID': -3,
        'TR': -1, 'PL': -3, 'KG': -2, 'MN': -2,
        'PH': -2, 'FQ': -3, 'ZG': -2, 'XL': -1,
        'TM': -1, 'ZC': -3, 'XH': -1, 'DR': -2,
        'BW': -4, 'XD': -1, 'ZK': 1, 'FA': -2,
        'ZW': -3, 'FE': -3, 'DN': 1, 'BK': 0,
        'XX': -1, 'FI': 0, 'BG': -1, 'XT': 0,
        'FM': 0, 'BC': -3, 'ZI': -3, 'ZV': -2,
        'SS': 4, 'LQ': -2, 'WE': -3, 'QR': 1,
        'NN': 6, 'WM': -1, 'QC': -3, 'WI': -3,
        'SC': -1, 'LA': -1, 'SG': 0, 'LE': -3,
        'WQ': -2, 'HG': -2, 'SK': 0, 'QN': 0,
        'NR': 0, 'HC': -3, 'YN': -2, 'GQ': -2,
        'YF': 3, 'CA': 0, 'VL': 1, 'GE': -2,
        'GA': 0, 'KR': 2, 'ED': 2, 'YR': -2,
        'MQ': 0, 'TI': -1, 'CD': -3, 'VF': -1,
        'TA': 0, 'TP': -1, 'BP': -2, 'TE': -1,
        'VN': -3, 'PG': -2, 'MA': -1, 'KH': -1,
        'VR': -3, 'PC': -3, 'ME': -2, 'KL': -2,
        'VV': 4, 'MI': 1, 'TQ': -1, 'IG': -4,
        'PK': -1, 'MM': 5, 'KD': -1, 'IC': -1,
        'ZD': 1, 'FR': -3, 'XK': -1, 'QD': 0,
        'XG': -1, 'ZL': -3, 'XC': -2, 'ZH': 0,
        'BL': -4, 'BH': 0, 'FF': 6, 'XW': -2,
        'BD': 4, 'DA': -2, 'SL': -2, 'XS': 0,
        'FN': -3, 'SR': -1, 'WD': -4, 'VY': -1,
        'WL': -2, 'HR': 0, 'WH': -2, 'HN': 1,
        'WT': -2, 'TT': 5, 'SF': -2, 'WP': -4,
        'LD': -4, 'BI': -3, 'LH': -3, 'SN': 1,
        'BT': -1, 'LL': 4, 'YK': -2, 'EQ': 2,
        'YG': -3, 'ZS': 0, 'YC': -2, 'GD': -1,
        'BV': -3, 'EA': -1, 'YW': 2, 'EE': 5,
        'YS': -2, 'CN': -3, 'VC': -1, 'TH': -2,
        'PR': -2, 'VG': -3, 'TL': -1, 'VK': -2,
        'KQ': 1, 'RA': -1, 'IR': -3, 'TD': -1,
        'PF': -4, 'IN': -3, 'KI': -3, 'MD': -3,
        'VW': -3, 'WW': 11, 'MH': -2, 'PN': -2,
        'KA': -1, 'ML': 2, 'KE': 1, 'ZE': 4,
        'XN': -1, 'ZA': -1, 'ZM': -1, 'XF': -1,
        'KC': -3, 'BQ': 0, 'XB': -1, 'BM': -3,
        'FC': -2, 'ZQ': 3, 'XZ': -1, 'FG': -3,
        'BE': 1, 'XV': -1, 'FK': -3, 'BA': -2,
        'XR': -1, 'DD': 6, 'WG': -2, 'ZF': -3,
        'SQ': 0, 'WC': -2, 'WK': -3, 'HQ': 0,
        'LC': -1, 'WN': -4, 'SA': 1, 'LG': -4,
        'WS': -3, 'SE': 0, 'HE': 0, 'SI': -2,
        'HA': -2, 'SM': -1, 'YL': -1, 'YH': 2,
        'YD': -3, 'ER': 0, 'XP': -2, 'GG': 6,
        'GC': -3, 'EN': 0, 'YT': -2, 'YP': -3,
        'TK': -1, 'AA': 4, 'PQ': -1, 'TC': -1,
        'VH': -3, 'TG': -2, 'IQ': -3, 'ZT': -1,
        'CR': -3, 'VP': -2, 'PE': -1, 'MC': -1,
        'KN': 0, 'II': 4, 'PA': -1, 'MG': -3,
        'TS': 1, 'IE': -3, 'PM': -2, 'MK': -1,
        'IA': -1, 'PI': -3, 'RR': 5, 'XM': -1,
        'LI': 2, 'XI': -1, 'ZB': 1, 'XE': -1,
        'ZN': 0, 'XA': 0, 'BR': -1, 'BN': 3,
        'FD': -3, 'XY': -1, 'ZR': 0, 'FH': -1,
        'BF': -3, 'FL': 0, 'XQ': -1, 'BB': 4
    }

    return blosum62[q + t]


def build_diag(a, b):
    diag = {}
    for item_a in a:
        if item_a in b:
            for start_a in a[item_a]:
                for start_b in b[item_a]:
                    d = start_b - start_a
                    if d not in diag:
                        diag[d] = {
                            'matches': 1,
                            'intervals': [[start_a, start_a]]
                        }
                    else:
                        diag[d]['matches'] += 1

                        tmp = 0
                        for elem in diag[d]['intervals']:
                            if elem[0] == start_a + 1:
                                elem[0] -= 1
                                tmp += 1
                            elif elem[1] == start_a - 1:
                                elem[1] += 1
                                tmp += 1
                        if tmp == 0:
                            diag[d]['intervals'].append([start_a, start_a])

    return {
        key: diag[key] for key, value in diag.items() if value['matches'] >= 10
    }


def recount_blossum62(a, b, diagonals):
    new = {}
    for key in diagonals:
        for subdiag in diagonals[key]['intervals']:
            length = subdiag[1] - subdiag[0] + 1

            new[(subdiag[0], key + subdiag[0])] = {
                'length': length,
                'score': sum(
                    get_blosum62(a[subdiag[0] + i], b[key + subdiag[0] + i])
                    for i in range(subdiag[1] - subdiag[0] + 1)
                )
            }
    return new


def group_up(a, n):
    result = {}
    for i in range(len(a) - n + 1):
        if a[i:i + n] in result:
            result[a[i:i + n]].append(i)
        else:
            result[a[i:i + n]] = [i]
    return result


def merge_diag(diagonals, gap_score):
    res = {}
    for diag_1 in list(diagonals.keys()):
        count = 0
        for diag_2 in list(diagonals.keys()):
            if diag_1 != diag_2:
                dist_1 = diag_2[0] - diag_1[0] + diagonals[diag_1]['length'] - 1
                dist_2 = diag_2[1] - diag_1[1] + diagonals[diag_1]['length'] - 1

                if dist_1 >= 0 and dist_2 >= 0:
                    score = (dist_1 + dist_2) * gap_score + diagonals[diag_1]['score'] + diagonals[diag_2]['score']
                    if score > diagonals[diag_1]['score']:
                        diagonals[(diag_1[0], diag_1[1])] = {
                            'length': diagonals[diag_1]['length'] + dist_1 + diagonals[diag_2]['length'],
                            'score': score
                        }
                        count += 1
        if count == 0: res[diag_1] = diagonals[diag_1]
        diagonals = {key: diagonals[key] for key in diagonals if key != diag_1}
    return res


if __name__ == "__main__":
    fasta_path = input('Path to .fasta file: ')
    search = int(input('Index of element to search: '))
    group_size = int(input('Input group size: '))
    gap_score = int(input('Input gap score: '))
    filter = int(input('Input score size filter: '))

    records = parse(fasta_path)
    search_elem = records[search]
    print(search_elem.name, search_elem.sequence)
    search_groups = group_up(search_elem.sequence, group_size)

    for i in tqdm(range(len(records))):
        record = records[i]
        if record.name == search_elem.name:
            continue

        diagonals = build_diag(search_groups, group_up(record.sequence, group_size))
        diagonals = recount_blossum62(search_elem.sequence, record.sequence, diagonals)
        diagonals = {key: diagonals[key] for key in diagonals if diagonals[key]['score'] >= filter}
        diagonals = merge_diag(diagonals, gap_score)

        result = [
            SW(search_elem.sequence[seq[0]: (seq[0] + diagonals[seq]['length'])],
               record.sequence[seq[1]: (seq[1] + diagonals[seq]['length'])],
               score_fun, gap_score)
            for seq in diagonals
        ]
        score = max(s[2] for s in result) if [s[2] for s in result] else -1

        if score > 0:
            print(search_elem.name, record.name, score)

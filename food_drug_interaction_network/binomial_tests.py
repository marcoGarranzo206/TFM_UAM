from collections import Counter
from collections import defaultdict
from scipy.stats import rankdata

def build_ref(ref_file):

    """
    Build reference from a file where the first column is
    compound and the rest are annotations for it
    """

    classification = defaultdict(set)
    counter = defaultdict(int)
    N = 0

    with open(ref_file, "r") as f:

        for line in f:

            N += 1
            splitted = line[:-1].split("\t")
            cp = splitted[0]
            labels = set(splitted[1:])
            classification[cp] = labels
            for label in labels:

                counter[label] += 1

    return counter, N, classification


def binom(query, ref_counts, N, adj = "BH", classification = None):

    """
    query: iterable of nodes
    classification: mapping of node name to classification
    ref counts: counts of each class in reference
    N: total number of nodes
    """

    counts = Counter()
    n = 0
    pvals = dict()

    for node in query:

        if classification is not None:

            if node in classification:

                counts.update(classification[node])

        else:

            counts.update([node])

        n += 1

    ret = {}
    classes = ref_counts.keys()
    pvals = np.zeros(len(classes))
    signs = np.zeros(len(classes))
    for i, class_ in enumerate(classes):

        p = ref_counts[class_]/N
        E = p*n
        pvals[i] = stats.binom_test(counts[class_], n = n, p = p)
        signs[i] = np.sign(counts[class_]-E)

    if adj == "BH":

        pvals = fdr(pvals)


    return {class_ : (p,s) for class_,p,s in zip(classes, pvals, signs)}


def fdr(p_vals):

    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr

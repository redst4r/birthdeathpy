from scipy import stats
import numpy as np
import pandas as pd
import seaborn as sns
import sys
import matplotlib.pyplot as plt
sys.path.append('/home/michi/postdoc_seattle/scBarcoding_RNASeqWorkflow/')
from bc_extract.upset_plot import upset_wrapper
from anytree import Node, RenderTree
import tqdm
import collections
def discrete_birth_death_process(initial_vector, p_div):
    """
    initial_vector: 1 x n_barcodes
    """

    assert initial_vector.shape[0] == 1 and initial_vector.ndim == 2
    timecourse = [initial_vector]
    current_vector = initial_vector
    gen_count = 0
    while np.sum(current_vector) < 1e6 and np.sum(current_vector) > 0:
        print()
        gen_count+=1
        """
        for each barcode, each cell throws a dice to survive or die
        if theres 100 cells, the number of cells in the next gen will be binomial: k = Binom(N=100, p)
        """
        dividing_cells = np.random.binomial(current_vector, p=p_div) # the number of cells that will divide
        next_gen = 2 * dividing_cells

        print(f'Gen: {gen_count}, cells {np.sum(current_vector)}')
        timecourse.append(next_gen)
        current_vector = next_gen

    A = np.concatenate(timecourse)

    B = A[:, A.sum(0)>0] # onle barcodes >0

#     plt.figure()
#     plt.subplot(221)
#     plt.plot(B)
#     plt.subplot(222)
#     plt.plot(B+1)
#     plt.yscale('log')
#     plt.subplot(223)
#     plt.scatter(B[0]+1, B[-1]+1)
#     plt.xscale('log')
#     plt.yscale('log')

#     plt.subplot(224)
#     p_ini = B[0] / B[0].sum()
#     p_final = B[-1] / B[-1].sum()
#     plt.scatter(p_ini, p_final, alpha=0.5)
#     plt.xscale('log')
#     plt.yscale('log')

    return A


def PCR_amplification_with_errors(initial_seq, rounds, per_base_error_rate):
    "simulating the tree exatly // slow"
    root = Node(initial_seq)
    current_level = [root]
    for cycle in range(rounds):
        print(cycle)
        next_level = []
        for n in current_level:
            # duplicate, potentially with an error
            s1 = [_ if np.random.rand() > per_base_error_rate else np.random.choice(['A','T','G','C'])for _ in n.name]
            s2 = [_ if np.random.rand() > per_base_error_rate else np.random.choice(['A','T','G','C'])for _ in n.name]
            s1 = "".join(s1)
            s2 = "".join(s2)
            next_level.append(Node(s1, parent=n))
            next_level.append(Node(s2, parent=n))
        current_level = next_level


    for pre, fill, node in RenderTree(root):
        print("%s%s" % (pre, node.name))

def print_tree(root):
    for pre, fill, node in RenderTree(root):
        print("%s%s" % (pre, node.name))

def mutate_1BP(input_seq):
    error_position = np.random.choice(len(input_seq))
    bases = [_ for _ in ['A','T','G','C'] if _!= input_seq[error_position]] #make sure that we mutate intoa different base
    error_nucleotide = np.random.choice(bases)
    s1 = [input_seq[i] if i!= error_position else error_nucleotide for i in range(len(input_seq))]
    s1 = "".join(s1)

#     print(s1)
#     print(seq)
#     print('--------------')
    assert s1!=input_seq, f'asked to mutate, but no mutation {error_position}, {error_nucleotide}, {s1}, {input_seq}, {bases}'
    return s1


def PCR_amplification_with_errors_faster(initial_seq_freq:dict, rounds, per_base_error_rate, prob_death, VERBOSE=False, track_tree=True):
    "simulating the tree exactly, but storing only the errors made // faster"
    "molecules can also disapear from the pool by chance"
    from anytree import Node, RenderTree
    bplength = len(initial_seq_freq[0]['seq'])
    n_mutations = 0
    if track_tree:
        root = Node({'seq': 'XXXXX', 'freq': 1})
        current_level = [Node({'seq': n['seq'], 'freq': n['freq']}, parent=root) for n in initial_seq_freq]
    else:
        root = None
        current_level = [{'seq': n['seq'], 'freq': n['freq']} for n in initial_seq_freq]

    for cycle in range(rounds):
        next_level = []
        total_entities = sum([n.name['freq'] for n in current_level]) if track_tree else sum([n['freq'] for n in current_level])
        print(f'cycle {cycle}:\t{total_entities:.3e} molecules total')
        for n in tqdm.tqdm(current_level, desc=f'Amplifiying// {len(current_level)} subclones'):

            if track_tree:
                seq, f = n.name['seq'], n.name['freq']
            else:
                seq, f = n['seq'], n['freq']

            # a small fraction will disapear
            died = np.random.binomial(f, p=prob_death)
#             print(f'{died}/{f} died')
            if f == died and VERBOSE:
                print(f'\t\t{n.name["seq"]} exterminated')
            f = f-died

            # most molecules wont have an error:
            error_free_prob = (1- per_base_error_rate)**bplength
            error_free = np.random.binomial(n=f, p=error_free_prob)
            if error_free > 0:
                if track_tree:
                    next_level.append(Node({'seq': seq, 'freq': 2*error_free}, parent=n))
                else:
                    next_level.append({'seq': seq, 'freq': 2*error_free})

            # prob_one_error = (1- per_base_error_rate)**(bplength -1) * per_base_error_rate * bplength
            one_error = f - error_free
            if one_error > 0 and VERBOSE:
                print('\t\tMutation')
            n_mutations += one_error
            for _ in range(one_error):
                s1 = mutate_1BP(seq)
                if track_tree:
                    next_level.append(Node({'seq': s1, 'freq': 2}, parent=n))
                else:
                    next_level.append({'seq': s1, 'freq': 2})

        current_level = next_level

    total_entities = sum([n.name['freq'] for n in current_level]) if track_tree else sum([n['freq'] for n in current_level])
    print(f'Final:\t{total_entities:.3e} molecules total, {n_mutations} errors, {len(next_level)} subclones')


    import collections
    final_freqs = collections.defaultdict(int)
    for n in current_level:
        if track_tree:
            final_freqs[n.name['seq']] = final_freqs[n.name['seq']] + n.name['freq']
        else:
            final_freqs[n['seq']] = final_freqs[n['seq']] + n['freq']

    return dict(final_freqs), root


def PCR_amplification_with_errors_faster_DICT(initial_seq_freq:dict, rounds, per_base_error_rate, prob_death, VERBOSE=False, track_tree=True):
    "same as above, but acting on a DICT instead of a list of dicts"

    assert isinstance(initial_seq_freq, dict), 'first argument MUST be a dict'

    bplength = len(list(initial_seq_freq.keys())[0])
    n_mutations = 0
    root = None
    current_level = initial_seq_freq.copy()

    for cycle in range(rounds):
        next_level = collections.defaultdict(int)
        total_entities = sum(list(current_level.values()))
        # total_entities = sum([n.name['freq'] for n in current_level]) if track_tree else sum([n['freq'] for n in current_level])

        print(f'cycle {cycle}:\t{total_entities:.3e} molecules total')
        for seq, f in tqdm.tqdm(current_level.items(), desc=f'Amplifiying// {len(current_level)} subclones'):

            # a small fraction will disapear
            died = np.random.binomial(f, p=prob_death)
#             print(f'{died}/{f} died')
            if f == died and VERBOSE:
                print(f'\t\t{seq} exterminated')
            f = f-died

            # most molecules wont have an error:
            error_free_prob = (1- per_base_error_rate)**bplength
            error_free = np.random.binomial(n=f, p=error_free_prob)
            if error_free > 0:
                next_level[seq]+=2*error_free

            # prob_one_error = (1- per_base_error_rate)**(bplength -1) * per_base_error_rate * bplength
            one_error = f - error_free
            if one_error > 0 and VERBOSE:
                print('\t\tMutation')
            n_mutations += one_error
            for _ in range(one_error):
                s1 = mutate_1BP(seq)
                next_level[s1] += 2

        current_level = next_level

    total_entities = sum(current_level.values())
    print(f'Final:\t{total_entities:.3e} molecules total, {n_mutations} errors, {len(next_level)} subclones')


    final_freqs = current_level.copy()
    return final_freqs, root



def main():
    import seaborn as sns
    import toolz

    nstart = 250
    initial_seq_freq = [
        {'seq': "".join(np.random.choice(['A','G','T','C'], 20)), 'freq':1} for _ in range(nstart)
    ]
    rounds  =np.ceil(np.log2(1e6) - np.log2(nstart))
    """
    lets simulate the expansion of a population of barcoded cells
    - barcodes might mutate
    - cells might die, ie, barcodes could go extinct
    """
    final_bcs, tree = PCR_amplification_with_errors_faster(initial_seq_freq, rounds=int(rounds),
                                         per_base_error_rate=1e-7,
                                         prob_death=0.01)

    for pre, fill, node in RenderTree(tree):
        print("%s%s" % (pre, node.name))

    # look at the clone-size distribution of the non mutated barcodes
    sns.distplot(list(final_bcs.values()))

    sns.distplot(list(toolz.dicttoolz.valfilter(lambda v: v>1000, final_bcs).values()))


def main_old():
    N = 1_000_000 # cells in the population
    n_barcodes = 100_000

    alphas = 0.005*np.ones(n_barcodes)  # relatively uniform
    # randomly increase the frequency of a hanful of barcode_Sequence
    # for i in np.random.choice(len(alphas), 20):
    #     alphas[i] += 30000

    p_barcode = stats.dirichlet(alphas).rvs().flatten()
    plt.figure()
    sns.distplot(p_barcode)
    cell_initial = [250, 250, 500, 500, 1000, 1000, 2000,2000, 4000, 4000, 8000, 8000]

    X = []
    neutral_drift = []
    p_div = 0.9
    for c in cell_initial:
        X_tmp = np.random.multinomial(c, p_barcode, size=1)
        X.append(X_tmp)

        # apply the neutral drift to the cells
        drift = discrete_birth_death_process(X_tmp, p_div)

        plt.figure()
        plt.subplot(131)
        plt.plot(drift[:,X_tmp.flatten()>0 ])
        plt.subplot(132)
        sns.distplot(drift[-1,X_tmp.flatten()>0 ])
        plt.subplot(133)
        sns.distplot(drift[0,X_tmp.flatten()>0 ])
        neutral_drift.append(drift[-1:])

    X = np.concatenate(X)
    neutral_drift = np.concatenate(neutral_drift)

    # compare inital vs drifted proportions
    plt.figure()
    for i in range(len(X)):
        ix_existing = X[i]>0 # only plot barcodes that got picked in the first place!
        x_ = X[i][ix_existing]
        y_ = neutral_drift[i][ix_existing]
    #     plt.scatter(x_, y_, s=5, alpha=0.5, label=i)
        plt.scatter(x_/ x_.sum(), y_/y_.sum(), s=5, alpha=0.5,  label=i)
    plt.legend()

    # compare clone size accross experiments
    plt.figure()
    n_exp = (neutral_drift>0).sum(0)
    for i in range(neutral_drift.shape[0]):
        for j in range(i+1, neutral_drift.shape[0]):
            print(i,j)
            plt.subplot2grid([neutral_drift.shape[0], neutral_drift.shape[0]], (i,j))
            plt.scatter(np.log10(neutral_drift[i][n_exp>1]),
                        np.log10(neutral_drift[j][n_exp>1]),
                        alpha=0.5, s=5)


    # do the set overlap either on the inital population or on the drifted
    DO_DRIFTED = True

    # its only about seeing the barcode vs not seeing it, the quantities are not meaninful (PCR bias)
    # per experiment, annotate if that barcode was seen
    Y = neutral_drift > 0 if DO_DRIFTED else X > 0

    # kick out all barcodes never seeing
    Y = Y[:,Y.sum(0)>0]

    df = pd.DataFrame(Y.T.astype('int'), columns=['experiment_%d'%i for i in range(1,len(cell_initial)+1)])
    # df['barcode'] = np.arange(len(df))
    upset_wrapper(df, '/tmp/upset.pdf')

if __name__ == '__main__':
    main()

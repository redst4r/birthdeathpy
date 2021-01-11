from scipy import stats
import numpy as np
import pandas as pd
import seaborn as sns
import tqdm
import collections
# import upsetplot
import fire
import toolz



def mutate_1BP(input_seq):
    error_position = np.random.choice(len(input_seq))
    bases = [_ for _ in ['A','T','G','C'] if _!= input_seq[error_position]] #make sure that we mutate intoa different base
    error_nucleotide = np.random.choice(bases)
    s1 = [input_seq[i] if i!= error_position else error_nucleotide for i in range(len(input_seq))]
    s1 = "".join(s1)

    assert s1!=input_seq, f'asked to mutate, but no mutation {error_position}, {error_nucleotide}, {s1}, {input_seq}, {bases}'
    return s1


def amplify_one_round(initial_seq_freq:dict, per_base_error_rate:float, prob_death:float, VERBOSE=False):
    current_level = initial_seq_freq.copy()
    next_level = collections.defaultdict(int)
    bplength = len(list(initial_seq_freq.keys())[0])

    for seq, f in tqdm.tqdm(current_level.items(), desc=f'Amplifiying// {len(current_level)} subclones'):

        # a small fraction will disapear
        died = np.random.binomial(f, p=prob_death)
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
        for _ in range(one_error):
            s1 = mutate_1BP(seq)
            next_level[s1] += 2

    return next_level


def PCR_amplification_with_errors_faster_DICT(initial_seq_freq:dict, rounds, per_base_error_rate, prob_death, VERBOSE=False):
    "same as above, but acting on a DICT instead of a list of dicts"

    assert isinstance(initial_seq_freq, dict), 'first argument MUST be a dict'
    current_level = initial_seq_freq.copy()

    for cycle in range(rounds):
        total_entities = sum(list(current_level.values()))
        print(f'cycle {cycle}:\t{total_entities:.3e} molecules total')
        next_level = amplify_one_round(current_level, per_base_error_rate, prob_death, VERBOSE=False)
        current_level = next_level

    total_entities = sum(current_level.values())
    print(f'Final:\t{total_entities:.3e} molecules total, {len(next_level)} subclones')

    final_freqs = current_level.copy()
    return final_freqs


def main():

    nstart = 2500
    initial_seq_freq = {
        "".join(np.random.choice(['A','G','T','C'], 20)):1 for _ in range(nstart)
    }
    rounds  = 50
    """
    lets simulate the expansion of a population of barcoded cells
    - barcodes might mutate
    - cells might die, ie, barcodes could go extinct
    """
    final_bcs = PCR_amplification_with_errors_faster_DICT(
                    initial_seq_freq, rounds=rounds,
                    per_base_error_rate=0,
                    prob_death=0.1)


    # look at the clone-size distribution of the non mutated barcodes
    sns.distplot(list(final_bcs.values()))

    sns.distplot(list(toolz.dicttoolz.valfilter(lambda v: v>1000, final_bcs).values()))

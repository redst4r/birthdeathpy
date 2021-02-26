import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import fire
from birthdeathpy.simulation import PCR_amplification_with_errors

#
# examples:
# python cli.py simulate_uniform_abundance --n_species 1000 --abundance 1 --rounds 20 --mutation_rate 0.001 --death_probability 0.2 --output_prefix test_
# python cli.py simulate_dirichlet_abundance --n_species 1000 --n_individuals 10000 --alpha 0.5 --rounds 20 --mutation_rate 0 --death_probability 0.2 --output_prefix test_


def timecourse_to_dataframe(timecourse):

    bcs = set(timecourse[0].keys())
    for t in timecourse[1:]:
        bcs = bcs | t.keys()
    bcs = sorted(list(bcs))

    df = []
    for cycle, t in enumerate(timecourse):
        q = pd.Series(t)
        q.name= cycle
        df.append(q)

    df = pd.DataFrame(df).T
    df[np.isnan(df)] = 0
    return df.astype(int)


def plot_results(output_prefix, timecourse):
    """
    clone size distribution at t=0 and t-final
    also a scatter of clone size evolution

    :param output_prefix: plots will be saved in a filename with this prefix
    :param initial_seq_freq: dict of the initial frequencies
    :param final_bcs: dict of the final frequencies
    """
    X0 = list(timecourse[0].values())
    Xt = list(timecourse[-1].values())
    # final clone size
    sns.displot(X0)
    plt.xlabel('Initial clone size')
    plt.ylabel('Frequency')
    plt.title(f'Initial clones: {len(X0)}')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}initial_clonesize.png', dpi=300)

    # final clone size
    sns.displot(Xt)
    plt.xlabel('Final clone size')
    plt.ylabel('Frequency')
    plt.title(f'Final clones: {len(Xt)}')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}final_clonesize.png', dpi=300)


    # the entries dont neccesearily match between t0 and tn
    # luckily, we can turn them into defautldicts
    all_clones = list(set(timecourse[0].keys()) | set(timecourse[-1].keys()))

    initial_seq_freq = timecourse[0]
    final_bcs = timecourse[-1]
    X0 = [initial_seq_freq[clone] if clone in initial_seq_freq else 0 for clone in all_clones]
    Xt = [final_bcs[clone] if clone in final_bcs else 0 for clone in all_clones]
    plt.figure()
    plt.scatter(X0, Xt, s=2, alpha=0.5)
    plt.xlabel('Initial clone size')
    plt.ylabel('Final clone size')
    plt.tight_layout()
    plt.savefig(f'{output_prefix}clonesize_scatter.png', dpi=300)


    df = timecourse_to_dataframe(timecourse)
    plt.figure()
    df = df.T
    if len(df) > 100:
        df = df.sample(100)
    plt.plot(df, alpha=0.4)
    plt.yscale('log')
    plt.xlabel('Cycles')
    plt.ylabel('Clone abundance')
    plt.savefig(f'{output_prefix}clonesize_plot.png', dpi=300)

    # X0 = [initial_seq_freq[clone] if clone in initial_seq_freq else 0 for clone in all_clones]
    # Xt = [final_bcs[clone] if clone in final_bcs else 0 for clone in all_clones]
    # plt.figure()
    # sns.jointplot(x=X0, y=Xt, marginal_kws=dict(bins=25))
    # plt.savefig(f'{output_prefix}2d_density.png', dpi=300)


def export_clonesizes(the_dict, fname):
    """
    turn the clonesize dictionary into a dataframe, and export it into a file

    :param the_dict: frequency dictionary
    :param fname: filename
    """
    s = pd.Series(the_dict, name='frequency')
    df = pd.DataFrame(s)
    df.to_csv(fname)


def simulate_uniform_abundance(n_species, abundance, rounds, mutation_rate, death_probability, output_prefix, no_plot=False):
    """
    each clone starts with the same abundance

    :param n_species: number of different clones
    :param abundance: the initial abundace of each clone, i.e. each clone has exactly that many cells
    :param rounds: Number of amplification/division rounds
    :param mutation_rate: per-base mutation rate, which could change the barcodes, i.e. spawn a new clone
    :param death_probability: probability for a single individual to die (before replication)
    :param output_prefix: all output files will be prefixed with this
    :param no_plot: if True, skip plotting
    """
    initial_seq_freq = {
            "".join(np.random.choice(['A','G','T','C'], 20)): abundance for _ in range(n_species)
        }
    timecourses = PCR_amplification_with_errors(
                    initial_seq_freq, rounds=rounds,
                    per_base_error_rate=mutation_rate,
                    prob_death=death_probability)
    final_bcs = timecourses[-1]

    if not no_plot:
        plot_results(output_prefix, timecourses)

    export_clonesizes(final_bcs, f'{output_prefix}clonesizes.csv')


def simulate_dirichlet_abundance(n_species, n_individuals, alpha, rounds, mutation_rate, death_probability, output_prefix,no_plot=False):
    """
    the initial abundace comes from a Dirichlet/Multinomial model

    :param n_species: number of different clones
    :param n_individuals: the total amount of cells/individuals at the beginning (distributed across the species according to the Dirichlet distribution)
    :param alpha: alpha parameter of the Dirichlet distribution that determines the initial distribution
    :param rounds: Number of amplification/division rounds
    :param mutation_rate: per-base mutation rate, which could change the barcodes, i.e. spawn a new clone
    :param death_probability: probability for a single individual to die (before replication)
    :param no_plot: if True, skip plotting
    :param output_prefix: all output files will be prefixed with this
    """
    p = np.random.dirichlet(np.full(n_species, alpha), size=1).flatten()
    x0 = np.random.multinomial(n_individuals, p, size=1).flatten()


    initial_seq_freq = {
            "".join(np.random.choice(['A','G','T','C'], 20)): x0[i] for i in range(n_species)
        }
    timecourses = PCR_amplification_with_errors(
                    initial_seq_freq, rounds=rounds,
                    per_base_error_rate=mutation_rate,
                    prob_death=death_probability)
    final_bcs = timecourses[-1]

    if not no_plot:
        print('plotting')
        plot_results(output_prefix, timecourses)
    export_clonesizes(final_bcs, f'{output_prefix}clonesizes.csv')


def main():
    fire.Fire({
        'simulate_uniform_abundance': simulate_uniform_abundance,
        'simulate_dirichlet_abundance': simulate_dirichlet_abundance,
    })
if __name__ == '__main__':
    main()

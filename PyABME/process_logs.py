import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import re

def process_logfile(filename):
    non_decimal = re.compile(r'[^\d.]+')

    with open(filename) as f:
        lengths, lengths_b, lengths_i = [], [], []
        num_inds, killed, born, died = [], [], [], []
        ages, ages_rep, lifespans = [], [], []
        mrf_b, mrf_i = [], []
        mri_b, mri_i = [], []
        mrd_b, mrd_i = [], []
        mrt_b, mrt_i = [], []
        mrm = []
        mrf_p = []

        for n, line in enumerate(f):
            search_string_chr_length = "[Overall] Avg. chromosome length: "
            search_string_chr_length_b = "[Behaviour] Avg. chromosome length: "
            search_string_chr_length_i = "[Interaction] Avg. chromosome length: "
            search_string_num_ind = "Num. individuals = "
            search_string_age = "Avg. age: "
            search_string_rep = "Avg. reproductive age: "
            search_string_lifespan = "Avg. lifespan: "
            search_string_mrf_b = "[Behaviour] Avg. mut. rate (flip): "
            search_string_mri_b = "[Behaviour] Avg. mut. rate (ins.): "
            search_string_mrd_b = "[Behaviour] Avg. mut. rate (del.): "
            search_string_mrt_b = "[Behaviour] Avg. mut. rate (trans.): "
            search_string_mrf_i = "[Interaction] Avg. mut. rate (flip): "
            search_string_mri_i = "[Interaction] Avg. mut. rate (ins.): "
            search_string_mrd_i = "[Interaction] Avg. mut. rate (del.): "
            search_string_mrt_i = "[Interaction] Avg. mut. rate (trans.): "
            search_string_mrf_p = "[Params] Avg. mut. rate (flip): "
            search_string_mrm = "Avg. mut. rate (meta): "

            
            if search_string_chr_length in line:
                lengths.append(float(line.replace(search_string_chr_length, "")))
            if search_string_chr_length_b in line:
                lengths_b.append(float(line.replace(search_string_chr_length_b, "")))
            if search_string_chr_length_i in line:
                lengths_i.append(float(line.replace(search_string_chr_length_i, "")))
            if search_string_num_ind in line:
                first = line.find(search_string_num_ind) + len(search_string_num_ind)
                first_par = line.find("(")
                born_sep = line.find(" born")
                cycle_sep = line.find("cycle, ")
                killed_sep = line.find(" killed, ")
                died_sep = line.find(" died")
                num_inds.append(int(non_decimal.sub('', line[first:first_par])))
                born.append(int(non_decimal.sub('', line[first_par:born_sep])))
                killed.append(int(non_decimal.sub('', line[cycle_sep:killed_sep])))
                died.append(int(non_decimal.sub('', line[killed_sep:died_sep])))
            if search_string_age in line:
                ages.append(float(line.replace(search_string_age, "")))
            if search_string_rep in line:
                ages_rep.append(float(line.replace(search_string_rep, "")))
            if search_string_lifespan in line:
                lifespans.append(float(line.replace(search_string_lifespan, "")))
            if search_string_mrf_b in line:
                mrf_b.append(float(line.replace(search_string_mrf_b, "")))
            if search_string_mri_b in line:
                mri_b.append(float(line.replace(search_string_mri_b, "")))
            if search_string_mrd_b in line:
                mrd_b.append(float(line.replace(search_string_mrd_b, "")))
            if search_string_mrt_b in line:
                mrt_b.append(float(line.replace(search_string_mrt_b, "")))
            if search_string_mrf_i in line:
                mrf_i.append(float(line.replace(search_string_mrf_i, "")))
            if search_string_mri_i in line:
                mri_i.append(float(line.replace(search_string_mri_i, "")))
            if search_string_mrd_i in line:
                mrd_i.append(float(line.replace(search_string_mrd_i, "")))
            if search_string_mrt_i in line:
                mrt_i.append(float(line.replace(search_string_mrt_i, "")))            
            if search_string_mrm in line:
                mrm.append(float(line.replace(search_string_mrm, "")))
            if search_string_mrf_p in line:
                mrf_p.append(float(line.replace(search_string_mrf_p, "")))

    fig, ax = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False)
    fig.set_size_inches(30, 20)
    
    ax[0,0].plot(lengths, label='Overall')
    ax[0,0].plot(lengths_b, label='Behaviour')
    ax[0,0].plot(lengths_i, label='Interaction')

    ax[0,1].plot(num_inds, label='Live individuals')
    ax[0,1].plot(born, label='Born')
    ax[0,1].plot(killed, label='Killed')
    ax[0,1].plot(died, label='Died')

    ax[1,0].plot(mrf_b, label="Flip")
    ax[1,0].plot(mri_b, label="Insertion")
    ax[1,0].plot(mrd_b, label="Deletion")
    ax[1,0].plot(mrt_b, label="Trans.")

    ax[1,1].plot(mrf_i, label="Flip")
    ax[1,1].plot(mri_i, label="Insertion")
    ax[1,1].plot(mrd_i, label="Deletion")
    ax[1,1].plot(mrt_i, label="Trans.")

    ax[2,0].plot(ages, label='Age')
    ax[2,0].plot(ages_rep, label='Repr. age')
    ax[2,0].plot(lifespans, label='Lifespans')

    ax[2,1].plot(mrm, label="Meta")
    ax[2,1].plot(mrf_p, label="Flip (param)")


    ax[2,0].set_xlabel("Iteration Number (x100)")
    ax[2,1].set_xlabel("Iteration Number (x100)")
    
    ax[0,0].set_ylabel("Chr. length")
    ax[0,1].set_ylabel("Cycle metrics")
    ax[2,0].set_ylabel("Ages")
    ax[1,0].set_ylabel("Beh. mut. rate")
    ax[1,1].set_ylabel("Int. mut. rate")
    ax[2,1].set_ylabel("Meta rates")

    ax[0,0].set_ylim(ymin=0)
    ax[0,1].set_ylim(ymin=0)
    ax[2,0].set_ylim(ymin=0)
    ax[1,0].set_ylim(ymin=0)
    ax[1,1].set_ylim(ymin=0)
    ax[2,1].set_ylim(ymin=0)

    ax[0,0].legend(loc="upper center", shadow=False, prop={'size': 6})
    ax[0,1].legend(loc="upper center", shadow=False, prop={'size': 6})
    ax[2,0].legend(loc="upper center", shadow=False, prop={'size': 6})
    ax[1,0].legend(loc="upper center", shadow=False, prop={'size': 6})
    ax[1,1].legend(loc="upper center", shadow=False, prop={'size': 6})
    ax[2,1].legend(loc="upper center", shadow=False, prop={'size': 6})
    
    plt.tight_layout()
    fig.savefig("{}_plots.png".format(filename))
    plt.close()
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--logfile', required=True, help='Log file to process')
    
    args = parser.parse_args()

    if args.logfile: process_logfile(filename=args.logfile)

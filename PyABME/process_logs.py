import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def process_logfile(filename):
    with open(filename) as f:
        lengths = []
        num_inds = []
        ages = []
        mrf = []
        mri = []
        mrd = []
        mrm = []

        for n, line in enumerate(f):
            search_string_chr_length = "Avg. chromosome length: "
            search_string_num_ind = "Num. individuals = "
            search_string_age = "Avg. age: "
            search_string_mrf = "Avg. mut. rate (flip): "
            search_string_mri = "Avg. mut. rate (ins.): "
            search_string_mrd = "Avg. mut. rate (del.): "
            search_string_mrm = "Avg. mut. rate (meta): "

            
            if search_string_chr_length in line:
                lengths.append(float(line.replace(search_string_chr_length, "")))
            if search_string_num_ind in line:
                first = line.find(search_string_num_ind) + len(search_string_num_ind)
                last = line.find("(")
                num_inds.append(int(line[first:last]))
            if search_string_age in line:
                ages.append(float(line.replace(search_string_age, "")))
            if search_string_mrf in line:
                mrf.append(float(line.replace(search_string_mrf, "")))
            if search_string_mri in line:
                mri.append(float(line.replace(search_string_mri, "")))
            if search_string_mrd in line:
                mrd.append(float(line.replace(search_string_mrd, "")))
            if search_string_mrm in line:
                mrm.append(float(line.replace(search_string_mrm, "")))

    fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True, sharey=False)
    ax[0].plot(lengths)
    ax[1].plot(num_inds)
    ax[2].plot(ages)
    ax[3].plot(mrf, label="Flip")
    ax[3].plot(mri, label="Insertion")
    ax[3].plot(mrd, label="Deletion")
    ax[3].plot(mrm, label="Meta")

    ax[3].set_xlabel("Iteration Number (x100)")
    ax[0].set_ylabel("Avg. chr. length")
    ax[1].set_ylabel("Num. individuals")
    ax[2].set_ylabel("Avg. age")
    ax[3].set_ylabel("Avg. mut. rate")
    ax[0].set_ylim(ymin=0)
    ax[1].set_ylim(ymin=0)
    ax[2].set_ylim(ymin=0)
    ax[3].set_ylim(ymin=0)
    ax[0].xaxis.set_major_locator(MaxNLocator(integer=True))
    ax[1].xaxis.set_major_locator(MaxNLocator(integer=True))
    ax[2].xaxis.set_major_locator(MaxNLocator(integer=True))
    ax[3].xaxis.set_major_locator(MaxNLocator(integer=True))
    ax[3].legend(loc="upper right", shadow=False, prop={'size': 6})
    plt.tight_layout()
    fig.savefig("{}_plots.png".format(filename))
    plt.close()
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--logfile', required=True, help='Log file to process')
    
    args = parser.parse_args()

    if args.logfile: process_logfile(filename=args.logfile)
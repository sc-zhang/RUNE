#!/usr/bin/env python3
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import time

mpl.use("Agg")


def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument("-i", "--input", help="Input unique kmer file", required=True)
    group.add_argument(
        "-g",
        "--group",
        help="Group file for filtering unique kmers and draw bars by group",
        required=True,
    )
    group.add_argument("-o", "--output", help="Output bar plot", required=True)
    return group.parse_args()


def time_print(info):
    print(
        "\033[32m%s\033[0m %s"
        % (time.strftime("[%H:%M:%S]", time.localtime(time.time())), info)
    )


def main():
    opts = get_opts()
    in_uniq_kmer_file = opts.input
    in_grp_file = opts.group
    out_pic = opts.output

    time_print("Loading group file")
    grp_order = []
    grp_db = {}
    smp_set = set()
    with open(in_grp_file, "r") as fin:
        for line in fin:
            data = line.strip().split()
            grp = data[0]
            smp = data[1]
            if grp not in grp_db:
                grp_db[grp] = []
                grp_order.append(grp)
            grp_db[grp].append(smp)
            smp_set.add(smp)

    time_print("Loading kmer file")
    cnt_db = {}
    with open(in_uniq_kmer_file, "r") as fin:
        for line in fin:
            data = line.strip().split()
            if data[1] not in smp_set:
                continue
            if data[1] not in cnt_db:
                cnt_db[data[1]] = 0
            cnt_db[data[1]] += 1

    time_print("Plotting bar plot")
    scnt = len(cnt_db)
    fig_width = scnt // 2
    fig_width = 10 if fig_width < 10 else fig_width
    fig_height = 8
    plt.figure(figsize=(fig_width, fig_height), dpi=100)
    idx = 0
    x_ticks = []
    for grp in grp_order:
        plt.bar(
            x=[_ for _ in range(idx, idx + len(grp_db[grp]))],
            height=[cnt_db[_] for _ in grp_db[grp]],
            width=0.8,
        )
        idx += len(grp_db[grp])
        for sid in grp_db[grp]:
            x_ticks.append(sid)
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xticks([_ for _ in range(len(x_ticks))], x_ticks, fontsize=15, rotation=-90)
    plt.xlim(-.5, len(x_ticks) - .5)
    plt.ylabel("Counts", fontsize=20)
    plt.yticks(fontsize=15)
    plt.savefig(out_pic, bbox_inches="tight")

    time_print("Finished")


if __name__ == "__main__":
    main()

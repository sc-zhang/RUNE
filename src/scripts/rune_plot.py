#!/usr/bin/env python3
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import argparse
import time

mpl.use("Agg")


def get_opts():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="sub commands")

    parser_bar_plot = subparsers.add_parser("bar", help="Draw bar plot")
    parser_bar_plot.add_argument("-i", "--input", help="Input unique kmer file", required=True)
    parser_bar_plot.add_argument(
        "-g",
        "--group",
        help="Group file for filtering unique kmers and draw bars by group",
        required=True,
    )
    parser_bar_plot.add_argument("-l", "--length", help="Input sequence length file", required=True)
    parser_bar_plot.add_argument("-o", "--output", help="Output bar plot", required=True)
    parser_bar_plot.set_defaults(func=draw_bar_plot)

    parser_line_plot = subparsers.add_parser("line", help="Draw line plot")
    parser_line_plot.add_argument("-i", "--input", help="Input unique kmer file", required=True)
    parser_line_plot.add_argument("-l", "--length", help="Input sequence length file", required=True)
    parser_line_plot.add_argument("-w", "--window", help="Window size, default=1e5", default="1e5")
    parser_line_plot.add_argument("-o", "--output", help="Output line plot", required=True)
    parser_line_plot.set_defaults(func=draw_line_plot)

    return parser


def time_print(info):
    print(
        "\033[32m%s\033[0m %s"
        % (time.strftime("[%H:%M:%S]", time.localtime(time.time())), info)
    )


def draw_bar_plot(args):
    in_uniq_kmer_file = args.input
    in_grp_file = args.group
    in_len_file = args.length
    out_pic = args.output
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

    time_print("Loading length file")
    len_db = {}
    with open(in_len_file, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            len_db[data[0]] = int(data[1])

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
            height=[cnt_db[_] * 1. / len_db[_] for _ in grp_db[grp]],
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
    plt.ylabel("Counts/Length", fontsize=20)
    plt.yticks(fontsize=15)
    plt.savefig(out_pic, bbox_inches="tight")

    time_print("Finished")


def draw_line_plot(args):
    in_uniq_kmer_file = args.input
    in_len_file = args.length
    win_size = float(args.window)
    out_pic = args.output

    time_print("Loading length file")
    len_db = {}
    chr_order = []
    with open(in_len_file, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            len_db[data[0]] = int(data[1])
            chr_order.append(data[0])

    time_print("Loading kmer file")
    region_cnt_db = {chrn: [0 for _ in range(math.ceil(len_db[chrn] * 1. / win_size))] for chrn in len_db}
    with open(in_uniq_kmer_file, "r") as fin:
        for line in fin:
            data = line.strip().split()
            if data[1] not in len_db:
                continue
            pos = int(data[2])
            idx = int(pos * 1. / win_size)
            region_cnt_db[data[1]][idx] += 1

    time_print("Plotting line plot")
    cnt = len(region_cnt_db)
    fig_width = cnt // 2
    fig_width = 10 if fig_width < 10 else fig_width
    fig_height = 8
    plt.figure(figsize=(fig_width, fig_height), dpi=100)
    offset = 0
    offset_list = []
    max_y = 0
    for chrn in chr_order:
        x = [offset + _ for _ in range(len(region_cnt_db[chrn]))]
        offset += len(region_cnt_db[chrn])
        max_y = max(max_y, max(region_cnt_db[chrn]))
        plt.plot(x, region_cnt_db[chrn], color="steelblue", lw=1)
        offset_list.append(offset)

    for offset in offset_list[:-1]:
        plt.plot([offset, offset], [0, max_y], lw=3, color='lightgrey', ls=":")

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    x_ticks = []
    for idx in range(len(offset_list)):
        if idx == 0:
            x_ticks.append(offset_list[idx] // 2)
        else:
            x_ticks.append((offset_list[idx] - offset_list[idx - 1]) // 2 + offset_list[idx - 1])

    plt.xticks(x_ticks, chr_order, fontsize=15, rotation=-90)
    plt.xlim(0, offset_list[-1])
    plt.xlabel("Bins", fontsize=20)
    plt.ylabel("Counts", fontsize=20)
    plt.yticks(fontsize=15)
    plt.ylim(0, max_y * 1.01)
    plt.savefig(out_pic, bbox_inches="tight")

    time_print("Finished")


def main():
    parser = get_opts()
    try:
        args = parser.parse_args()
        args.func(args)
    except AttributeError:
        parser.print_help()


if __name__ == "__main__":
    main()

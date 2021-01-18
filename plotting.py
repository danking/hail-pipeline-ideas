def plot_agg_hist(hist: hl.Struct):
    binsp = hist.bin_edges
    weights = hist.bin_freq
    bin_width = 1.0 / (len(bins) - 1)
    plt.figure()
    plt.bar(bins[:-1], weights, width = bin_width)

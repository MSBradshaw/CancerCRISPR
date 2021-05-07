import pandas as pd
import matplotlib.pyplot as plt
import pysam
import statistics
import os
import ntpath
import argparse


def get_top_x_regions(data, x=10):
    """
    :param data: dataframe with columns ['chr', 'start', 'end', 'count', 'z']
    :param x: number of top sites to return, default 10
    :return: pandas data frame with the x sites with the highest z score
    """
    data = data.sort_values('z', ascending=False)
    if data.shape[0] > 5000:
        data = data.iloc[0:5000, :]
    sub_data = data
    # remove ones that over lap with any windows above them
    for i in range(x):
        # start fall in side current range
        sub_data = sub_data.reset_index(drop=True)
        cond1 = (sub_data['start'] >= sub_data.iloc[i, 1]) & (sub_data['start'] <= sub_data.iloc[i, 2])
        # start fall in side current range
        cond2 = (sub_data['end'] >= sub_data.iloc[i, 1]) & (sub_data['end'] <= sub_data.iloc[i, 2])
        # on current chr
        cond3 = (sub_data['chr'] <= sub_data.iloc[i, 0])
        # is not self (+1 because reset index makes it
        cond4 = (sub_data.index != i)
        sub_data = sub_data[~(((cond1 | cond2) & cond3) & cond4)]
    return sub_data.iloc[0:x, :]


def plot_sliding_window(data, output, threshold):
    """
    Saves a scatter of the SNP counts and z scores
    :param data: dataframe with columns ['chr', 'start', 'end', 'count', 'z']
    :param output: filename for outptu
    :param threshold: int the threshold used
    :return: None
    """
    fig, axes = plt.subplots(23, 2)
    for ax, col in zip(axes[0], ['DP > 10', 'Z Score # DP > 10']):
        ax.set_title(col)
    # for each of the chromosomes
    for i, c in enumerate([str(x) for x in range(1, 22)] + ['X', 'Y']):
        print('plotting chromosome ' + c)
        sub = data[data['chr'] == c]
        axes[i, 0].scatter(sub['start'], sub['count'], s=5)
        axes[i, 1].scatter(sub['start'], sub['z'], s=5)
        axes[i, 0].set_xticks([])
        axes[i, 1].set_xticks([])
        axes[i, 0].set_ylim([-5, 150])
        axes[i, 1].set_ylim([-5, 150])
        axes[i, 1].set_xlabel('Chromosome ' + c)
        axes[i, 0].set_xlabel('Chromosome ' + c)

    fig.set_size_inches(10, 46)
    plt.ylabel('Freq DP > ' + str(threshold))
    plt.tight_layout()
    plt.savefig(output)
    plt.clf()


def avg(l):
    l = [x for x in l if x != 0]
    return sum(l) / float(len(l))


def sd(l):
    l = [x for x in l if x != 0]
    return statistics.stdev(l)


def zscore(row, averages, st_devs):
    return row['count'] - float(averages[row['chr']]) / st_devs[row['chr']]


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf',
                        dest='vcf',
                        help='BGZIPED and TABIXED vcf file')

    parser.add_argument('--plot_output',
                        dest='plot_output',
                        help='path for where to save plot, if not used plot is not generated')

    parser.add_argument('--top_x',
                        dest='top_x',
                        default=10,
                        help='number of highest z score windows to report on (default 10)')

    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help='path for output file to be saved')

    parser.add_argument('--dp',
                        dest='dp',
                        default=10,
                        help='DP treshold for counting a SNP in the vcf as a real mutation (default 10)')

    parser.add_argument('--window_size',
                        dest='window',
                        default=1000,
                        help='window size (default 1000)')

    parser.add_argument('--shift_size',
                        dest='shift',
                        default=100,
                        help='amount to slide the window by (default 100), should be <= than the window size')

    args = parser.parse_args()
    return args


def save_top_snps(data, vcf_obj, output):
    # get the top x SNPs in each region, ranked based on allele depth abs(difference) / number of alleles
    for i in range(data.shape[0]):
        row = data.iloc[i, :]
        print('----------------------------------')
        print(row['chr'])
        stuff = list(vcf_obj.fetch(row['chr'], row['start'], row['end']))
        info = {'chr': [],
                  'start': [],
                  'end': [],
                  'ref': [],
                  'alt': [],
                  'AD': [],
                  'ratio': []}
        for j in range(len(stuff)):
            print('\t' + str(stuff[j].contig))
            tup = stuff[j].samples[0]['AD']
            # do the sum of absolute differences
            numerator = 0
            for k in range(len(tup)):
                biggest = max(tup)
                numerator += abs(tup[k] - biggest)
            ratio = numerator / float(sum(tup))
            info['chr'].append(stuff[j].contig)
            info['start'].append(stuff[j].start)
            info['end'].append(stuff[j].stop)
            info['AD'].append(','.join([str(x) for x in stuff[j].samples[0]['AD']]))
            info['ratio'].append(ratio)
            info['ref'].append(stuff[j].ref)
            info['alt'].append(','.join(stuff[j].alts))
        temp = pd.DataFrame(info)
        temp = temp.sort_values('ratio', ascending=False)
        with open(output, 'a') as output_file:
            # with open('del.txt', 'a') as output_file:
            output_file.write('index\t' + '\t'.join(list(data.columns)) + '\n')
            row = [str(x) for x in row]
            output_file.write(str(i) + '\t' + '\t'.join(row) + '\n')
        temp.to_csv(output, mode='a', sep='\t', index=False)
    # temp.to_csv('del.txt', mode='a', sep='\t', index=False)


args = get_args()

# get all the DP's
vcf = pysam.VariantFile(args.vcf)
DP_threshold = int(args.dp)
window_size = int(args.window)
shift_size = int(args.shift)

# get just the name of the input file
input_file = ntpath.basename(args.vcf)
file_path = 'sliding_window_' + input_file + '_' + '_'.join(
    [str(x) for x in [DP_threshold, window_size, shift_size]]) + '.tsv'

# check if file already exists
if not os.path.exists(file_path):
    print('Generating Count and Z score information,\nThis may take a while.')
    starts = []
    ends = []
    chrs = []
    counts = []
    for contig in vcf.header.contigs:
        print(contig)
        for i in range(0, vcf.header.contigs[contig].length, shift_size):
            sites = list(vcf.fetch(contig, i, i + window_size))
            over_tresh = [x.info.get('DP') > DP_threshold for x in sites]
            starts.append(i)
            ends.append(i + window_size)
            chrs.append(contig)
            counts.append(len(over_tresh))

    # make a df of the data
    df = pd.DataFrame({'chr': chrs, 'start': starts, 'end': ends, 'count': counts})
    # plot one for each chrom
    df = df[df['chr'].isin([str(x) for x in range(22)] + ['X', 'Y'])]

    avgs = {c: avg(list(df[df['chr'] == c]['count'])) for c in set(df['chr'])}
    sds = {c: sd(list(df[df['chr'] == c]['count'])) for c in set(df['chr'])}
    df['z'] = [zscore(df.iloc[i, :], avgs, sds) for i in range(df.shape[0])]

    df['chr'] = df['chr'].astype(str)
    df.to_csv(file_path, sep='\t', index=False)
else:
    print('Loading count and z score information from cached file')
    df = pd.read_csv(file_path, sep='\t')

# ensure that all chromosomes are strings
df['chr'] = df['chr'].astype(str)

# remove the sites with 0 mutations
nzdf = df[df['count'] != 0]

# remove the huge unnecessary df
del df

# get the x windows with the highest z score
top_x = get_top_x_regions(nzdf, int(args.top_x))

# make the figure
if args.plot_output is not None:
    plot_sliding_window(nzdf, args.plot_output, DP_threshold)

# remove the big unnecessary df
del nzdf

# pull the SNPs from the VCF in the top regions, get their AD ratios sort and save them
save_top_snps(top_x, vcf, args.output)

"""
Use gnomad pop gen information
can we simulate the mutation process and the ability of this tool to pluck out the best mutations
"""

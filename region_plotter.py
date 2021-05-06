import pandas as pd
import matplotlib.pyplot as plt
import pysam
import statistics

# get all the DP's

vcf = pysam.VariantFile('output.g.vcf.gz')
ads = []
dps = []
dps_count = {}
BIN_SIZE = 25000
DP_threshold = 10

for contig in vcf.header.contigs:
    for i in range(0, vcf.header.contigs[contig].length, BIN_SIZE):
        if contig == '1' and i < 1000000:
            print(str(int(i / BIN_SIZE)))
        key = str(contig) + '_' + str(int(i / BIN_SIZE))
        dps_count[key] = 0

for r in vcf.fetch():
    num_samples = len(r.samples)
    for i in range(num_samples):
        sample = r.samples[i]
        ad = sample['AD']
        dps.append(sample['DP'])
        bin_num = int(r.start / BIN_SIZE)
        key = str(r.contig) + '_' + str(bin_num)
        if sample['DP'] > DP_threshold:
            dps_count[key] += 1
        try:
            ratio = ad[1] / (float(ad[0]) + ad[1])
        except ZeroDivisionError:
            ratio = 0
        ads.append(ratio)

plt.hist(ads)
# plt.yscale('log')
plt.xlabel('ALT freq / (ALT freq + REF freq)')
plt.ylabel('Occurrence Count')
plt.show()
plt.savefig('Figures/AD_ratio_hist.png')

plt.hist(dps)
plt.yscale('log')
plt.xlabel('DP')
plt.ylabel('Occurrence Count')
plt.show()
plt.savefig('Figures/DP_hist.png')

dps_count

keys = list(dps_count.keys())
dps_count_df = pd.DataFrame({'bin_id': keys,
                             'count': [dps_count[x] for x in keys],
                             'chromosome': [x.split('_')[0] for x in keys]})

g = dps_count_df.groupby('chromosome')
fig, axes = plt.subplots(23, 2)
for ax, col in zip(axes[0], ['DP > 10', 'Z Score # DP > 10']):
    ax.set_title(col)
for i, contig in enumerate(dps_count_df.chromosome.unique()):
    if i > 22:
        break
    sub = g.get_group(contig).reset_index()
    mean = sub['count'].sum() / float(sub.shape[0])
    std = statistics.stdev(list(sub['count']))
    z_scores = [(x - mean) / std for x in list(sub['count'])]
    sub['z_scores'] = z_scores
    sub['start'] = [int(x.split('_')[1]) * BIN_SIZE for x in sub['bin_id']]
    sub['end'] = [(int(x.split('_')[1]) * BIN_SIZE) + BIN_SIZE for x in sub['bin_id']]
    ss = sub[['chromosome','start','end','count','z_scores']]
    ss[ss['z_scores'] >= 5].to_csv('z_score_over_size.bed', mode='a', header=False,index=False)
    axes[i, 0].bar(sub['bin_id'], sub['count'])
    axes[i, 1].bar(sub['bin_id'], z_scores)
    axes[i, 0].set_xticks([])
    axes[i, 1].set_xticks([])
    axes[i, 1].set_ylim([-5, 20])
    axes[i, 1].set_xlabel('Chromosome ' + contig)
    axes[i, 0].set_xlabel('Chromosome ' + contig)

fig.set_size_inches(10, 46)
plt.ylabel('Freq DP > ' + str(DP_threshold))
plt.tight_layout()
plt.savefig('Figures/DP.png')
plt.show()


# TODO
"""
In the regions output to the bed file, which SNPs in these regions are present in virtually all reads 
Future and which are not in the patient's germline vcf 
"""



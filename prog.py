import pybedtools 

a = pybedtools.BedTool("GSM439776_p63_1_peaks.bed")
b = pybedtools.BedTool("GSM439777_p63_2_peaks.bed")
c = pybedtools.BedTool("GSM538930_p63_3_peaks.bed") 
refgenes = pybedtools.BedTool("refGene.txt")

peak1 = a.intersect(b, f = 0.2, r = True)
peak2 = b.intersect(c, f = 0.2,r = True)
peak3 = c.intersect(a, f = 0.2,r = True)
list1 = peak1.window(peak2,w = 1000)
list = list1.window(peak3,w = 1000)
print list

import pybedtools

for i in list:
    feat1 = pybedtools.cbedtools.create_interval_from_list([i[0],i[2],i[2]])
    mid1 = pybedtools.featurefuncs.midpoint(feat1)
    feat2 = pybedtools.cbedtools.create_interval_from_list([i[4],i[5],i[6]])
    mid2 = pybedtools.featurefuncs.center(feat2)
    feat3 = pybedtools.cbedtools.create_interval_from_list([i[8],i[9],i[10]])
    mid3 = pybedtools.featurefuncs.center(feat3)
    average = int(mid1[1] + mid2[1] + mid3[1])/3

print average 

filer1 = peak1.filter(b.start < 10000)
filter2 = peak2.filter(b.start < 10000)
filter3 = peak3.filter(b.start < 10000)

a1 = peak1.annotate(files = refgenes)
a2 = peak2.annotate(files = refgenes)
a3 = peak3.annotate(files = refgenes)

b1 = peak1.closest('refgenes', s = True)
b2 = peak2.closest('refgenes', s = True)
b3 = peak3.closest('refgenes', s = True)


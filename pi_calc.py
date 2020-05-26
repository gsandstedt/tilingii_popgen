import vcf
import random

v = vcf.Reader(filename='4g_2nas_alltil_high_varfilt_PASS_DP.vcf')
s = v.samples

def rd_bases1(x):
    g = list(x.split('/'))
    a1 = g[0]
    return a1

def rd_bases2(x):
    g = list(x.split('/'))
    a2 = g[1]
    return a2

# these two functions grab the bases read at a given site, and if read, splits them into a list
# rd_bases1 returns the first base, and rd_bases2 returns the second base; this works b/c items in a list have a set order

out = open('MIM.scaffold_1.txt', 'w')   # open a text file where output will be written
b = sorted(s)                                 # b is the column headings, this line addes a column heading for each sample line
b.insert(0,'pos')                               # insert column heading for site position
b.insert(0,'scaffold')                      # insert column heading for scaffold
#b.append('N.alleles')
out.write('\t'.join(b)+'\n')            # this line writes the first line to the output file, which is the column headings
curr_scaff = 'scaffold_1'

for rec in v:
    scaff = rec.CHROM
    if scaff != curr_scaff:
        out.close()
        out = open('MIM.'+scaff+'.txt', 'w')
        b = sorted(s)
        b.insert(0,'pos')                               # insert column heading for site position
        b.insert(0,'scaffold')                      # insert column heading for scaffold
#        b.append('N.alleles')
        out.write('\t'.join(b)+'\n')
        curr_scaff = scaff
    pos = rec.POS
    m = [scaff, str(pos)]   # this line creates a list that will be the data for each site in the outputed text file, pos 0 and 1 are given to 'scaff' and 'position'
    alleleA = 0
    alleleC = 0
    alleleG = 0
    alleleT = 0
    allelecount = 0
    n = 0
    for line in sorted(s):
        base = rec.genotype(line).gt_bases
        if base: # line genotyped

            if rec.is_monomorphic == True: #if site is monomorphic
                if rec.genotype(line).data[1] >= 10:         # miminum depth requirement coded in
                    if rec.genotype(line).data[1] <= 500:
                        n += 1
                        gen = rd_bases1(base)
                        m.append(gen)
                        if gen == 'A':  #counting number of alleles of each base in MIM samples
                            alleleA += 1
                        if gen == 'C':
                            alleleC += 1
                        if gen == 'G':
                            alleleG += 1
                        if gen == 'T':
                            alleleT += 1
                    else:   # line depth greather than max threshold, do not use
                        m.append('N')
                else:       # line depth less than min requirement
                    m.append('N')

            else: #site is not monomorphic, genotype calls contain more data
                if rec.genotype(line).data[2] >= 10:         # miminum depth requirement coded in
                    if rec.genotype(line).data[2] <= 500:
                        n += 1
                        gen1 = rd_bases1(base)
                        gen2 = rd_bases2(base)
                        gen = random.choice([gen1,gen2]) #random choice so that don't lose hets from data set! #asign gen so grab random choice once, and can go to same call
                        m.append(gen)   #m.append(rd_bases3(base))
                        if gen == 'A':  #counting number of alleles of each base in MIM samples
                            alleleA += 1
                        if gen == 'C':
                            alleleC += 1
                        if gen == 'G':
                            alleleG += 1
                        if gen == 'T':
                            alleleT += 1
                    else:   # line depth greather than max threshold, do not use
                        m.append('N')
                else:       # line depth less than min requirement
                    m.append('N')

        else:   # line not genotyped
            m.append('N')

    if alleleA != 0:     #count number of alleles per site in samples
        allelecount += 1
    if alleleC != 0:
        allelecount += 1
    if alleleG != 0:
        allelecount += 1
    if alleleT != 0:
        allelecount += 1

#    m.append(str(allelecount)) #append number of alleles

    if allelecount <= 2:             # code in limit on #of alleles to biallelic
        if n >= 2:                          # need at least two lines genotyped
            out.write('\t'.join(m)+'\n')
    if not pos % 1000 :
        print scaff, pos, line, base

out.close()

from collections import Counter
import re

geneindex = {}
reverseindex = {}
indexfile = open("BnaGeneindex.v0.txt")
fields = indexfile.readline().rstrip().split('\t')
for line in indexfile:
    l = [re.sub(r'Bna(.*)G0',r'Bna\1T0',x) for x in line.rstrip().split('\t')]
    geneindex[l[0]] = {}
    for k,v in zip(fields[1:],l[1:]):
        geneindex[l[0]][k] = v
        reverseindex[v] = l[0]

cultivars = {'GG':'ganganF73','NO':'no2127','QA':'quintaA','SL':'shengli3','TA':'tapidor3','WE':'westar','ZY':'zheyou73','ZS':'zs11'}

darmorfile = open('../BnapusDarmor_bzh.fasta')
darmor = [line[1:].rstrip() for line in darmorfile if line[0]=='>']

diamondforward = {}
diamondreverse = {}
for abbr,cul in cultivars.items():
    diamondforward[abbr] = {}
    diamondreverse[abbr] = {}
    print(abbr)
    file = open(f'darmor-bzh_{cul}.all.v0.pep.dmnd.tsv')
    for line in file:
        l = line.rstrip().split('\t')
        query = l[0]
        target = l[1]
        bitscore = float(l[11])
        if target not in diamondforward[abbr]:
            diamondforward[abbr][target] = Counter()
        
        diamondforward[abbr][target][query] += bitscore
    
    file = open(f'{cul}.all.v0.pep_darmor-bzh.tsv')
    for line in file:
        l = line.rstrip().split('\t')
        query = l[0]
        target = l[1]
        bitscore = float(l[11])
        if target not in diamondreverse[abbr]:
            diamondreverse[abbr][target] = Counter()
        
        diamondreverse[abbr][target][query] += bitscore

def getgene(hubID):
    g = geneindex[hubID]
    out = Counter()
    for k,v in g.items():
        if k in diamondforward:
            if v in diamondforward[k]:
                out += diamondforward[k][v]
    return out


def getreverse(bnaID):
    out = Counter()
    for cul in cultivars.keys():
        if cul in diamondreverse:
            if bnaID in diamondreverse[cul]:
                for hit,score in diamondreverse[cul][bnaID].most_common():
                    if hit in reverseindex:
                        out[reverseindex[hit]] += score
    return out


if __name__ == '__main__':
    outfile = open('diamond_forward_darmor-v10.tsv','w')
    for id in geneindex.keys():
        hits = getgene(id)
        if len(hits) == 0:
            trash = outfile.write(f'{id}\t-\t-\n')
        else:
            for k,v in hits.most_common():
                trash = outfile.write(f'{id}\t{k}\t{round(float(v),1)}\n')

    outfile.close()
    
    outfile = open('diamond_reverse_darmor-v10.tsv','w')
    for id in darmor:
        hits = getreverse(id)
        if len(hits) == 0:
            trash = outfile.write(f'{id}\t-\t-\n')
        else:
            for k,v in hits.most_common():
                trash = outfile.write(f'{id}\t{k}\t{round(float(v),1)}\n')

    outfile.close()


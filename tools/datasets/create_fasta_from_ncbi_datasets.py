import os
import json
import sys
import time
from Bio import SeqIO

suffix = 1
size = 0
TOTAL_SIZE = float(sys.argv[1]) * 1e+9
databases = sys.argv[2]

databases = databases.split(',')

fstream = open('contamination_{}.fsa'.format(suffix), 'w')
fstream_tax = open('contamination_{}_taxid'.format(suffix), 'w')
fstream_idx = open('contamination_{}.idx'.format(suffix), 'w')

ids = {}
for db in databases:
    print('\nParsing : ' + db)
    assemblies = {}
    with open('{}/ncbi_dataset/data/assembly_data_report.jsonl'.format(db)) as fjson:
        for line in fjson.readlines():
            d = json.loads(line)
            assemblies[d['assemblyInfo']['assemblyAccession']] = d['taxId']
    count = 0
    total = len(assemblies)
    start = time.time()
    for s in assemblies:
        count += 1
        files = [f for ds, dr, files in os.walk('{}/ncbi_dataset/data/{}'.format(db, s)) for f in files if f.endswith('.fna')]
        for f in files:
            for r in SeqIO.parse('{}/ncbi_dataset/data/{}/{}'.format(db, s, f), 'fasta'):
                if not r.id.startswith('NW_') and not r.id.startswith('NZ_'):
                    v = ids.setdefault(r.id, False)
                    if not v:
                        l = len(r.seq)
                        if size + l > TOTAL_SIZE:
                            ids = {}
                            suffix += 1
                            size = l
                            fstream.close()
                            fstream_tax.close()
                            fstream_idx.close()
                            fstream = open('contamination_{}.fsa'.format(suffix), 'w')
                            fstream_tax = open('contamination_{}_taxid'.format(suffix), 'w')
                            fstream_idx = open('contamination_{}.idx'.format(suffix), 'w')
                        else:
                            size += l
                        ids[r.id] = True
                        fstream_tax.write('{}\t{}\n'.format(r.id, assemblies[s]))
                        fstream_idx.write('{}\t{}\t{}\n'.format(r.id, fstream.tell(), assemblies[s]))
                        fstream.write(r.format('fasta'))
        end = time.time()
        print('{}/{} {:.1f}%. Assembly: {:16s} Taxa: {:7d} Output: contamination_{}.fsa ({:.2f}/{:.2f} GB) '.format(
            count, total, count * 100/total,
            s, assemblies[s], suffix, size/1e+9, TOTAL_SIZE/1e+9),
            'Elapsed: {:.0f} sec Remaining {:.0f} sec'.format( (end - start), ((end-start)/count) * (total - count)), end='\r')
print()
fstream.close()
fstream_tax.close()
fstream_idx.close()

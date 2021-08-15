#!/usr/bin/python3
"""
Anna Sidera 2021.
"""



import json
import operator



# Read OAS File #
# ============= #
def read_oas_file(filename):
    imgt_length = 106
    with open(filename, 'r') as f:
        lines = f.read().splitlines()
    lines.pop(0)
    output_dictionary = {}
    for line in lines:
        data_level_a = json.loads(line)
        s = data_level_a['v']
        (gene, allele) = s.split('*')
        data_level_b = data_level_a['data']
        data_level_b = json.loads(data_level_b)
        sequence = []
        for i in range(0, imgt_length):
            sequence.append('.')
        for k in data_level_b:
            data_level_c = data_level_b[k]
            for k1 in data_level_c:
                if k != 'cdrh3' or k1 == '105' or k1 == '106':
                    index = int(k1) - 1
                    if index < imgt_length:
                        sequence[index] = data_level_c[k1]
        if sequence[9] == '.' and sequence[72] == '.':
            s = ''
            sequence = s.join(sequence)
            if gene not in output_dictionary:
                output_dictionary[gene] = {}
            if allele not in output_dictionary[gene]:
                output_dictionary[gene][allele] = []
            output_dictionary[gene][allele].append(sequence)
    return output_dictionary



# Read IMGT File #
# ============== #
def read_imgt_file():
    imgt_length = 106

    with open('codon_translation.txt', 'r') as f:
        a = f.read().splitlines()
    codon_translation = {}
    for s in a:
        b = s.split(',')
        codon_translation[b[0]] = b[1]
    codon_translation['...'] = '.'
    with open('imgt_fasta.txt', 'r') as f:
        s = f.read()
    a = s.split('>')
    a.pop(0)
    a = [s.replace('\n', '') for s in a]
    genes = {}
    for s in a:
        b = s.split('|')
        if b[3] == 'F' and b[13] == ' ':
            s0 = b[1]
            (s1, s2) = s0.split('*')
            s3 = b[15].upper()
            s3 = s3[0:imgt_length * 3]
            i = 0
            j = 0
            s4 = ''
            while i < imgt_length * 3 and j == 0:
                s5 = s3[i:i + 3]
                if s5 in codon_translation:
                    s4 += codon_translation[s5]
                    i += 3
                else:
                    j = 1
            if j == 0:
                if s1 not in genes:
                    genes[s1] = {}
                genes[s1][s2] = s4
    return genes



# Many Genes Plot #
# =============== #
def many_genes_plot(filename, min_number, min_ratio, IgM_filename):
    imgt_d = read_imgt_file()
    imgt_length = 106
    stat = []
    for i in range(0, imgt_length):
        a = []
        for s1 in imgt_d:
            for s2 in imgt_d[s1]:
                s = imgt_d[s1][s2]
                if s[i] != '.' and s[i] not in a:
                    a.append(s[i])
        j = len(a)
        stat.append(j)
    s = ''
    for n in stat:
        s += str(round(float(n) / 20, 2))
        s += ' '
    s += '\n'
    oas_d = read_oas_file(IgM_filename)
    allele_d = {}
    for gene in oas_d:
        alleles = []
        a = []
        for allele in oas_d[gene]:
            a.append((len(oas_d[gene][allele]), allele))
        a.sort(reverse = True)
        if a[0][0] >= min_number:
            if len(a) == 1 or a[0][0] >= a[1][0] * min_ratio:
                alleles.append(a[0][1]) 
            elif a[1][0] >= min_number:
                if len(a) == 2 or a[1][0] >= a[2][0] * min_ratio:
                    alleles.append(a[0][1])
                    alleles.append(a[1][1])
        if len(alleles) > 0:
            allele_d[gene] = alleles
    oas_d = read_oas_file(filename)
    plot_d = {}
    for gene in oas_d:
        for allele in oas_d[gene]:
            if gene in allele_d:
                if allele in allele_d[gene]:
                    if gene not in plot_d:
                        plot_d[gene] = {}
                    plot_d[gene][allele] = oas_d[gene][allele]
    stat = []
    n = 0
    for i in range(0, imgt_length):
        stat.append(0)
    for gene in plot_d:
        for allele in plot_d[gene]:
            if gene in imgt_d:
                if allele in imgt_d[gene]:
                    s1 = imgt_d[gene][allele]
                    for s2 in plot_d[gene][allele]:
                        for i in range(0, imgt_length):
                            if s1[i] != s2[i]:
                                stat[i] += 1
                    n += len(plot_d[gene][allele])
    for i in range(0, imgt_length):
        s += str(round(float(stat[i]) / n, 2))
        s += ' '
    s += '\n'
    f = open('many_genes_plot.txt', 'w')
    f.write(s)
    f.close()



# Many Genes Heatmap #
# ================== #
def many_genes_heatmap(filename, min_number, min_ratio, IgM_filename, plot_thr_s, plot_thr_a):
    imgt_d = read_imgt_file()
    imgt_length = 106
    oas_d = read_oas_file(IgM_filename)
    alleles = []
    for gene in oas_d:
        a = []
        for allele in oas_d[gene]:
            a.append((len(oas_d[gene][allele]), allele))
        a.sort(reverse = True)
        if a[0][0] >= min_number:
            if len(a) == 1 or a[0][0] >= a[1][0] * min_ratio:
                alleles.append((gene, a[0][1]))
            elif a[1][0] >= min_number:
                if len(a) == 2 or a[1][0] >= a[2][0] * min_ratio:
                    alleles.append((gene, a[0][1]))
                    alleles.append((gene, a[1][1]))
    alleles.sort()
    oas_d = read_oas_file(filename)
    s = ''
    s3 = ''
    count = 0
    b = []
    for i in range(0, imgt_length):
        b.append(0)
    m = 0
    for (gene, allele) in alleles:
        if gene in oas_d:
            if allele in oas_d[gene]:
                if gene in imgt_d:
                    if allele in imgt_d[gene]:
                        a = []
                        for i in range(0, imgt_length):
                            a.append(0)
                        s1 = imgt_d[gene][allele]
                        for s2 in oas_d[gene][allele]:
                            for i in range(0, imgt_length):
                                if s1[i] != s2[i]:
                                    a[i] += 1
                        n = len(oas_d[gene][allele])
                        if n >= plot_thr_s and count < plot_thr_a:
                            for i in range(0, imgt_length):
                                s += str(round(float(a[i]) / n, 2))
                                s += ' '
                            s += '\n'
                            s3 += gene[4:]
                            s3 += '*'
                            s3 += allele
                            s3 += '('
                            s3 += str(n)
                            s3 += ')'
                            s3 += ' '
                            count += 1
                        for i in range(0, imgt_length):
                            b[i] += a[i]
                        m += n
    for i in range(0, imgt_length):
        s += str(round(float(b[i]) / m, 2))
        s += ' '
    s += '\n'
    s3 += 'all_oas('
    s3 += str(m)
    s3 += ')'
    s3 += ' '
    b = []
    for i in range(0, imgt_length):
        a = []
        for s1 in imgt_d:
            for s2 in imgt_d[s1]:
                s0 = imgt_d[s1][s2]
                if s0[i] != '.' and s0[i] not in a:
                    a.append(s0[i])
        j = len(a)
        b.append(j)
    for i in range(0, imgt_length):
        s += str(round(float(b[i]) / 20, 2))
        s += ' '
    s += '\n'
    s3 += 'germline'
    s3 += ' '
    s3 += '\n'
    f = open('many_genes_heatmap.txt', 'w')
    f.write(s)
    f.write(s3)
    f.close()



# Single Gene Heatmap #
# =================== #
def single_gene_heatmap(filenames, gene, min_number, min_ratio, IgM_filenames, reassign = 0):
    imgt_d = read_imgt_file()
    imgt_length = 106
    s = ''
    s3 = ''
    b = []
    for i in range(0, imgt_length):
        b.append(0)
    m = 0
    for index in range(0, len(filenames)):
        filename = IgM_filenames[index]
        oas_d = read_oas_file(filename)
        alleles = []
        a = []
        for allele in oas_d[gene]:
            a.append((len(oas_d[gene][allele]), allele))
        a.sort(reverse = True)
        if a[0][0] >= min_number:
            if len(a) == 1 or a[0][0] >= a[1][0] * min_ratio:
                alleles.append(a[0][1]) 
            elif a[1][0] >= min_number:
                if len(a) == 2 or a[1][0] >= a[2][0] * min_ratio:
                    alleles.append(a[0][1])
                    alleles.append(a[1][1])
        filename = filenames[index]
        oas_d = read_oas_file(filename)
        for allele in alleles:
            oas_alleles = []
            if reassign == 0 or len(alleles) > 1:
                oas_alleles.append(allele)
            else:
                for oas_allele in oas_d[gene]:
                    oas_alleles.append(oas_allele)
            a = []
            for i in range(0, imgt_length):
                a.append(0)
            n = 0
            s1 = imgt_d[gene][allele]
            for oas_allele in oas_alleles:
                for s2 in oas_d[gene][oas_allele]:
                    for i in range(0, imgt_length):
                        if s1[i] != s2[i]:
                            a[i] += 1
                n += len(oas_d[gene][oas_allele])
            for i in range(0, imgt_length):
                s += str(round(float(a[i]) / n, 2))
                s += ' '
            s += '\n'
            s3 += gene[4:]
            s3 += '*'
            s3 += allele
            s3 += '('
            s3 += str(n)
            s3 += ')'
            s3 += str(index + 1)
            s3 += ' '
            for i in range(0, imgt_length):
                b[i] += a[i]
            m += n
    for i in range(0, imgt_length):
        s += str(round(float(b[i]) / m, 2))
        s += ' '
    s += '\n'
    s3 += 'all_oas('
    s3 += str(m)
    s3 += ')'
    s3 += ' '
    b = []
    for i in range(0, imgt_length):
        a = []
        for s1 in imgt_d:
            for s2 in imgt_d[s1]:
                s0 = imgt_d[s1][s2]
                if s0[i] != '.' and s0[i] not in a:
                    a.append(s0[i])
        j = len(a)
        b.append(j)
    for i in range(0, imgt_length):
        s += str(round(float(b[i]) / 20, 2))
        s += ' '
    s += '\n'
    s3 += 'germline'
    s3 += ' '
    s3 += '\n'
    f = open('single_gene_heatmap.txt', 'w')
    f.write(s)
    f.write(s3)
    f.close()



# IMGT Print #
# ========== #
def imgt_print():
    imgt_d = read_imgt_file()
    imgt_length = 106
    s = ''
    for gene in imgt_d:
        for allele in imgt_d[gene]:
            s += gene
            s += '*'
            s += allele
            s += ' '
            s += imgt_d[gene][allele]
            s += '\n'
    f = open('imgt_print.txt', 'w')
    f.write(s)
    f.close()



# First Sanity Check #
# ================== #
def first_sanity_check(filename):
    oas_d = read_oas_file(filename)
    s = ''
    for gene in oas_d:
        s += gene
        s += ':'
        for allele in oas_d[gene]:
            s += allele
            s += '('
            s += str(len(oas_d[gene][allele]))
            s += ')'
        s += '\n'
    f = open('first_sanity_check.txt', 'w')
    f.write(s)
    f.close()



# Second Sanity Check #
# =================== #
def second_sanity_check(filename):
    imgt_d = read_imgt_file()
    oas_d = read_oas_file(filename)
    imgt_length = 106
    with open(filename, 'r') as f:
        lines = f.read().splitlines()
    lines.pop(0)
    d_vj = {}
    d_v = {}
    d_j = {}
    for line in lines:
        data_level_a = json.loads(line)
        s = data_level_a['v']
        (gene_v, allele_v) = s.split('*')
        s = data_level_a['j']
        (gene_j, allele_j) = s.split('*')
        if gene_v not in d_vj:
            d_vj[gene_v] = {}
        if allele_v not in d_vj[gene_v]:
            d_vj[gene_v][allele_v] = {}
        if gene_j not in d_vj[gene_v][allele_v]:
            d_vj[gene_v][allele_v][gene_j] = {}
        if allele_j not in d_vj[gene_v][allele_v][gene_j]:
            d_vj[gene_v][allele_v][gene_j][allele_j] = 0
        d_vj[gene_v][allele_v][gene_j][allele_j] += 1
        if gene_v not in d_v:
            d_v[gene_v] = {}
        if allele_v not in d_v[gene_v]:
            d_v[gene_v][allele_v] = 0
        d_v[gene_v][allele_v] += 1
        if gene_j not in d_j:
            d_j[gene_j] = {}
        if allele_j not in d_j[gene_j]:
            d_j[gene_j][allele_j] = 0
        d_j[gene_j][allele_j] += 1
    s = ''
    for gene_j in d_j:
        for allele_j in d_j[gene_j]:
            s += gene_j
            s += '*'
            s += allele_j
            s += '('
            s += str(d_j[gene_j][allele_j])
            s += ')'
            s += '\n'
    s += '\n'
    for gene_v in d_vj:
        for allele_v in d_vj[gene_v]:
            s += gene_v
            s += '*'
            s += allele_v
            s += '('
            s += str(d_v[gene_v][allele_v])
            s += ')'
            for gene_j in d_vj[gene_v][allele_v]:
                for allele_j in d_vj[gene_v][allele_v][gene_j]:
                    s += gene_j
                    s += '*'
                    s += allele_j
                    s += '('
                    s += str(d_vj[gene_v][allele_v][gene_j][allele_j])
                    s += ')'
            s += '\n'
    f = open('second_sanity_check.txt', 'w')
    f.write(s)
    f.close()



# Third Sanity Check #
# ================== #
def third_sanity_check(filename, min_number, min_ratio):
    oas_d = read_oas_file(filename)
    clean_d = {}
    for gene in oas_d:
        a = []
        for allele in oas_d[gene]:
            a.append(len(oas_d[gene][allele]))
        a.sort(reverse = True)
        if a[0] >= min_number:
            if len(a) == 1 or a[0] >= a[1] * min_ratio:
                clean_d[gene] = oas_d[gene] 
            elif a[1] >= min_number:
                if len(a) == 2 or a[1] >= a[2] * min_ratio:
                    clean_d[gene] = oas_d[gene]
    n = 0
    for gene in oas_d:
        for allele in oas_d[gene]:
            n += len(oas_d[gene][allele])
    m = 0
    for gene in clean_d:
        for allele in clean_d[gene]:
            m += len(clean_d[gene][allele])
    s = ''
    s += str(round(float(m) / float(n), 2))
    s += '\n\n'
    for gene in clean_d:
        s += gene
        s += ':'
        for allele in clean_d[gene]:
            s += allele
            s += '('
            s += str(len(clean_d[gene][allele]))
            s += ')'
        s += '\n'
    f = open('third_sanity_check.txt', 'w')
    f.write(s)
    f.close()



# Many Genes Correlation #
# ====================== #
def many_genes_correlation(filename, min_number, min_ratio, IgM_filename):
    imgt_d = read_imgt_file()
    imgt_length = 106
    oas_d = read_oas_file(IgM_filename)
    alleles = []
    for gene in oas_d:
        a = []
        for allele in oas_d[gene]:
            a.append((len(oas_d[gene][allele]), allele))
        a.sort(reverse = True)
        if a[0][0] >= min_number:
            if len(a) == 1 or a[0][0] >= a[1][0] * min_ratio:
                alleles.append((gene, a[0][1]))
            elif a[1][0] >= min_number:
                if len(a) == 2 or a[1][0] >= a[2][0] * min_ratio:
                    alleles.append((gene, a[0][1]))
                    alleles.append((gene, a[1][1]))
    oas_d = read_oas_file(filename)
    s = ''
    b = []
    for i in range(0, imgt_length):
        b.append(0)
    m = 0
    for (gene, allele) in alleles:
        if gene in oas_d:
            if allele in oas_d[gene]:
                if gene in imgt_d:
                    if allele in imgt_d[gene]:
                        a = []
                        for i in range(0, imgt_length):
                            a.append(0)
                        s1 = imgt_d[gene][allele]
                        for s2 in oas_d[gene][allele]:
                            for i in range(0, imgt_length):
                                if s1[i] != s2[i]:
                                    a[i] += 1
                        n = len(oas_d[gene][allele])
                        for i in range(0, imgt_length):
                            b[i] += a[i]
                        m += n
    for i in range(0, imgt_length):
        s += str(round(float(b[i]) / m, 2))
        s += ' '
    s += '\n'
    b = []
    for i in range(0, imgt_length):
        a = []
        for s1 in imgt_d:
            for s2 in imgt_d[s1]:
                s0 = imgt_d[s1][s2]
                if s0[i] != '.' and s0[i] not in a:
                    a.append(s0[i])
        j = len(a)
        b.append(j)
    for i in range(0, imgt_length):
        s += str(round(float(b[i]) / 20, 2))
        s += ' '
    s += '\n'
    f = open('many_genes_correlation.txt', 'w')
    f.write(s)
    f.close()



# Single Gene Correlation #
# ======================= #
def single_gene_correlation(filenames, gene, min_number, min_ratio, IgM_filenames):
    imgt_d = read_imgt_file()
    imgt_length = 106
    s = ''
    s3 = ''
    b = []
    for i in range(0, imgt_length):
        b.append(0)
    m = 0
    for index in range(0, len(filenames)):
        filename = IgM_filenames[index]
        oas_d = read_oas_file(filename)
        alleles = []
        a = []
        for allele in oas_d[gene]:
            a.append((len(oas_d[gene][allele]), allele))
        a.sort(reverse = True)
        if a[0][0] >= min_number:
            if len(a) == 1 or a[0][0] >= a[1][0] * min_ratio:
                alleles.append(a[0][1]) 
            elif a[1][0] >= min_number:
                if len(a) == 2 or a[1][0] >= a[2][0] * min_ratio:
                    alleles.append(a[0][1])
                    alleles.append(a[1][1])
        if len(alleles) > 0:
            filename = filenames[index]
            oas_d = read_oas_file(filename)
            a = []
            for i in range(0, imgt_length):
                a.append(0)
            n = 0
            for allele in alleles:     
                s1 = imgt_d[gene][allele]
                for s2 in oas_d[gene][allele]:
                    for i in range(0, imgt_length):
                        if s1[i] != s2[i]:
                            a[i] += 1
                n += len(oas_d[gene][allele])
            for i in range(0, imgt_length):
                s += str(round(float(a[i]) / n, 2))
                s += ' '
            s += '\n'
            s3 += str(index + 1)
            s3 += '('
            s3 += str(n)
            s3 += ')'
            s3 += ' '
            for i in range(0, imgt_length):
                b[i] += a[i]
            m += n
    for i in range(0, imgt_length):
        s += str(round(float(b[i]) / m, 2))
        s += ' '
    s += '\n'
    s3 += 'all_oas('
    s3 += str(m)
    s3 += ')'
    s3 += ' '
    b = []
    for i in range(0, imgt_length):
        a = []
        for s1 in imgt_d:
            for s2 in imgt_d[s1]:
                s0 = imgt_d[s1][s2]
                if s0[i] != '.' and s0[i] not in a:
                    a.append(s0[i])
        j = len(a)
        b.append(j)
    for i in range(0, imgt_length):
        s += str(round(float(b[i]) / 20, 2))
        s += ' '
    s += '\n'
    s3 += 'germline'
    s3 += ' '
    s3 += '\n'
    f = open('single_gene_correlation.txt', 'w')
    f.write(s)
    f.write(s3)
    f.close()



# Many Genes Coverage #
# =================== #
def many_genes_coverage(filename, plot_thr_s, plot_thr_g):
    imgt_d = read_imgt_file()
    imgt_length = 106
    oas_d = read_oas_file(filename)
    genes = []
    for gene in oas_d:
        genes.append(gene)
    genes.sort()
    s = ''
    s3 = ''
    count = 0
    for gene in genes:
        if gene in imgt_d and count < plot_thr_g:
            n = 0
            for allele in oas_d[gene]:
                if allele in imgt_d[gene]:
                    n += len(oas_d[gene][allele])
            if n >= plot_thr_s:
                for i in range(0, imgt_length):
                    a = []
                    for allele in imgt_d[gene]:
                        s0 = imgt_d[gene][allele]
                        if s0[i] != '.' and s0[i] not in a:
                            a.append(s0[i])
                    b = []
                    for allele in oas_d[gene]:
                        if allele in imgt_d[gene]:
                            for s0 in oas_d[gene][allele]:
                                if s0[i] != '.' and s0[i] not in b:
                                    b.append(s0[i])
                    if len(a) > 0:
                        m = 0
                        for s0 in a:
                            if s0 in b:
                                m += 1
                        c = float(len(a) - m) / float(len(a))
                    else:
                        c = 0
                    s += str(round(c, 2))
                    s += ' '
                s += '\n'
                s3 += gene[4:]
                s3 += '('
                s3 += str(n)
                s3 += ')'
                s3 += ' '
                count += 1
    f = open('many_genes_coverage.txt', 'w')
    f.write(s)
    f.write(s3)
    f.close()



# Single Gene Coverage #
# ==================== #
def single_gene_coverage(filenames, gene):
    imgt_d = read_imgt_file()
    imgt_length = 106
    s = ''
    s3 = ''
    if gene in imgt_d:
        for index in range(0, len(filenames)):
            filename = filenames[index]
            oas_d = read_oas_file(filename)
            if gene in oas_d:
                n = 0
                for allele in oas_d[gene]:
                    if allele in imgt_d[gene]:
                        n += len(oas_d[gene][allele])
                if n > 0:
                    for i in range(0, imgt_length):
                        a = []
                        for allele in imgt_d[gene]:
                            s0 = imgt_d[gene][allele]
                            if s0[i] != '.' and s0[i] not in a:
                                a.append(s0[i])
                        b = []
                        for allele in oas_d[gene]:
                            if allele in imgt_d[gene]:
                                for s0 in oas_d[gene][allele]:
                                    if s0[i] != '.' and s0[i] not in b:
                                        b.append(s0[i])
                        if len(a) > 0:
                            m = 0
                            for s0 in a:
                                if s0 in b:
                                    m += 1
                            c = float(len(a) - m) / float(len(a))
                        else:
                            c = 0
                        s += str(round(c, 2))
                        s += ' '
                    s += '\n'
                    s3 += str(index + 1)
                    s3 += '('
                    s3 += str(n)
                    s3 += ')'
                    s3 += ' '  
    f = open('single_gene_coverage.txt', 'w')
    f.write(s)
    f.write(s3)
    f.close()



# Single Gene Table #
# =================== #
def single_gene_table(filenames, gene, min_number, min_ratio, IgM_filenames):
    imgt_d = read_imgt_file()
    imgt_length = 106
    s = ''
    position = []
    aminoacid = []
    denominator = []
    numerator = []
    if gene in imgt_d:
        for i in range(0, imgt_length):
            a = []
            for allele in imgt_d[gene]:
                s0 = imgt_d[gene][allele]
                if s0[i] not in a:
                    a.append(s0[i])
            if len(a) > 1:
                for s1 in a:
                    position.append(i)
                    aminoacid.append(s1)
                    denominator.append(0)
                    numerator.append(0)
    k = len(position)
    total = 0
    for index in range(0, len(filenames)):
        filename = IgM_filenames[index]
        oas_d = read_oas_file(filename)
        alleles = []
        if gene in oas_d:
            a = []
            for allele in oas_d[gene]:
                a.append((len(oas_d[gene][allele]), allele))
            a.sort(reverse = True)
            if a[0][0] >= min_number:
                if len(a) == 1 or a[0][0] >= a[1][0] * min_ratio:
                    alleles.append(a[0][1]) 
                elif a[1][0] >= min_number:
                    if len(a) == 2 or a[1][0] >= a[2][0] * min_ratio:
                        alleles.append(a[0][1])
                        alleles.append(a[1][1])
        if len(alleles) > 0:
            filename = filenames[index]
            oas_d = read_oas_file(filename)
            n = 0
            if gene in oas_d:
                for allele in oas_d[gene]:
                    n += len(oas_d[gene][allele])
            if n > 0:
                total += 1
                for j in range(0, k):
                    i = position[j]
                    n = 0
                    for allele in alleles:
                        if allele in imgt_d[gene]:
                            s0 = imgt_d[gene][allele]
                            if s0[i] == aminoacid[j]:
                                n = 1
                    if n == 0:
                        denominator[j] += 1
                        m = 0
                        for allele in oas_d[gene]:
                            for s0 in oas_d[gene][allele]:
                                if s0[i] == aminoacid[j]:
                                    m = 1
                        if m == 0:
                            numerator[j] += 1
    s += str(total)
    s += '\n\n'
    for j in range(0, k):
        s += str(position[j] + 1)
        s += ' '
        s += aminoacid[j]
        s += ' '
        s += str(denominator[j])
        s += ' '
        s += str(numerator[j])
        s += ' '
        s += '\n'
    f = open('single_gene_table.txt', 'w')
    f.write(s)
    f.close()



# Single Gene Extra #
# ================= #
def single_gene_extra(filenames, gene, min_number, min_ratio, IgM_filenames, pos, ami):
    imgt_d = read_imgt_file()
    imgt_length = 106
    pos -= 1
    s = ''
    s3 = ''
    s4 = ''
    position = []
    aminoacid = []
    denominator = []
    numerator = []
    if gene in imgt_d:
        for i in range(0, imgt_length):
            a = []
            for allele in imgt_d[gene]:
                s0 = imgt_d[gene][allele]
                if s0[i] not in a:
                    a.append(s0[i])
            if len(a) > 1:
                for s1 in a:
                    position.append(i)
                    aminoacid.append(s1)
                    denominator.append(0)
                    numerator.append(0)
    k = len(position)
    total = 0
    for index in range(0, len(filenames)):
        filename = IgM_filenames[index]
        oas_d = read_oas_file(filename)
        alleles = []
        if gene in oas_d:
            a = []
            for allele in oas_d[gene]:
                a.append((len(oas_d[gene][allele]), allele))
            a.sort(reverse = True)
            if a[0][0] >= min_number:
                if len(a) == 1 or a[0][0] >= a[1][0] * min_ratio:
                    alleles.append(a[0][1]) 
                elif a[1][0] >= min_number:
                    if len(a) == 2 or a[1][0] >= a[2][0] * min_ratio:
                        alleles.append(a[0][1])
                        alleles.append(a[1][1])
        if len(alleles) > 0:
            filename = filenames[index]
            oas_d = read_oas_file(filename)
            n = 0
            if gene in oas_d:
                for allele in oas_d[gene]:
                    n += len(oas_d[gene][allele])
            if n > 0:
                total += 1
                for j in range(0, k):
                    i = position[j]
                    n = 0
                    for allele in alleles:
                        if allele in imgt_d[gene]:
                            s0 = imgt_d[gene][allele]
                            if s0[i] == aminoacid[j]:
                                n = 1
                    if n == 0:
                        denominator[j] += 1
                        if i == pos and aminoacid[j] == ami:
                            d = {}
                        m = 0
                        for allele in oas_d[gene]:
                            for s0 in oas_d[gene][allele]:
                                if s0[i] == aminoacid[j]:
                                    m = 1
                                if i == pos and aminoacid[j] == ami:
                                    if s0[i] not in d:
                                        d[s0[i]] = 0
                                    d[s0[i]] += 1
                        if m == 0:
                            numerator[j] += 1
                            if i == pos and aminoacid[j] == ami:
                                for s1 in alleles:
                                    s4 += s1
                                    s4 += ' '
                                s4 += '\n'
                        if i == pos and aminoacid[j] == ami:
                            a = d.items()
                            a = sorted(a, key = operator.itemgetter(1), reverse = True)
                            for b in a:
                                s3 += b[0]
                                s3 += '('
                                s3 += str(b[1])
                                s3 += ')'
                            s3 += '\n'
    s += str(total)
    s += '\n\n'
    for j in range(0, k):
        s += str(position[j] + 1)
        s += ' '
        s += aminoacid[j]
        s += ' '
        s += str(denominator[j])
        s += ' '
        s += str(numerator[j])
        s += ' '
        s += '\n'
    s += '\n'
    s3 += '\n'
    s4 += '\n'
    f = open('single_gene_extra.txt', 'w')
    f.write(s)
    f.write(s3)
    f.write(s4)
    f.close()



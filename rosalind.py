# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 08:48:20 2023

@author: jacob
"""

from math import comb

from math import log10

def prop_heterozygous(n, m):
    s = 0
    for i in range(m):
        s += comb(2**n, i)*(.75)**(2**n-i)*(.25)**i
    return 1-s
            
def count_connected(graph):
    visited=set()
    count=0
    def dfs(node):
        if node in visited:
            return False
        visited.add(node)
        for neighbor in graph[node][1]:
            dfs(neighbor)
        return True
    for node in graph:
        if dfs(node[0]):
            count+=1
    return count

    
def graph_from_adj_list(filename):
    file = open(filename, 'r')
    nodes = [ [a,[]] for a in range(int(file.readline()))]
    edges = file.readlines()
    for edge in edges:
        a, b = int(edge.split()[0]), int(edge.split()[1])
        nodes[a-1][1].append(b-1)
        nodes[b-1][1].append(a-1)
    return  nodes
    

def nt_count(string):
    count = {'A':0, 'C':0, 'G':0,'U':0}
    for char in string:
        count[char]+=1
    return count

def transcribe_rna(string):
    new_string = ''
    for i in range(len(string)):
        if string[i] =='T':
            new_string += 'U'
        else:
            new_string += string[i]
    return new_string

def reverse_complement(string):
    new_string = ''
    for char in string:
        if char == 'A':
            new_string = 'T' + new_string
        elif char == 'C':
            new_string = 'G' + new_string
        elif char == 'G':
            new_string = 'C' + new_string
        else:
            new_string = 'A' + new_string
    return new_string

def wascally_wabbits(n,k):
    tokenize = [1,1, 1+k]
    for i in range(3, n):
        tokenize.append(tokenize[-1] + k*tokenize[-2])
    return tokenize[-1]

def gc_content(string):
    content = 0
    for char in string:
        content += int(char in ['G', 'C'])
    return content/len(string)

def highest_gc(list_of_names, list_of_strings):
    gc_contents = []
    for string in list_of_strings:
        gc_contents.append(gc_content(string))
    return list_of_names[gc_contents.index(max(gc_contents))], max(gc_contents)

def hamming_distance(s,t):
    count = 0
    for i in range(len(s)):
        count += int(not s[i] == t[i])
    return count

def read_in_list_of_strings(file_name):
    file = open(file_name, 'r')
    reader = file.readlines()
    names = []
    snippets = []
    bigline = ""
    for line in reader:
        if "Ros" in line:
            snippets.append(bigline)
            names.append(line[1:-1])
            bigline = ""
        else:
            bigline += line[:-1]
    snippets.append(bigline)
    return names, snippets[1:]

def overlap_graph(names, snippets, k):
    prefixes = [snippet[:k] for snippet in snippets]
    suffixes = [snippet[-k:] for snippet in snippets]
    for i in range(len(snippets)):
        for j in range(len(snippets)):
            if suffixes[i] == prefixes[j] and i != j:
                print(names[i] + ' ' + names[j])
    pass

def gluestrings(snippets):
    newsnips = snippets.copy()
    while len(newsnips) >= 40:
        for k in range(min([len(snippet)//2 for snippet in snippets]), min([len(snippet) for snippet in snippets])):
            for i in range(1,len(newsnips)):
                if newsnips[i].startswith(newsnips[0][-k:]):
                    newsnips[0] = newsnips[0] + newsnips[i].lstrip(newsnips[0][-k:])
                    for snip in newsnips:
                        if snip in newsnips[0]:
                            newsnips.remove(snip)
                    break
            for i in range(1,len(newsnips)):
                if newsnips[0].startswith(newsnips[i][-k:]):
                    newsnips[0] = newsnips[i] + newsnips[0].lstrip(newsnips[i][-k:])
                    for snip in newsnips:
                        if snip in newsnips[0]:
                            newsnips.remove(snip)
                    break
    return newsnips
        
def overlap_distance(s1, s2):
    p = min(len(s1), len(s2))
    m = 0
    for k in range(p):
        if s1[-k:] == s2[:k]:
            m = k
    return m

def gluetwo(snippets):
    newsnips = snippets.copy()
    m = 0
    pair = (0,0)
    for i in range(len(snippets)):
        for j in range(len(snippets)):
            if overlap_distance(snippets[i], snippets[j]) > m:
                m = overlap_distance(snippets[i], snippets[j])
                pair = (i,j)
    s1, s2 = newsnips[pair[0]], newsnips[pair[1]]
    newsnips.pop(max(pair[0], pair[1]))
    newsnips.pop(min(pair[0], pair[1]))
    newsnips.append(s1 + s2[m:])
    return newsnips

def gluemany(snippets):
    while len(snippets) >= 2:
        snippets = gluetwo(snippets)
    return snippets
    

def find_spliced_motif(s,t):
    locs = []
    loc = -1
    for char in t:
        for i in range(loc +1,len(s)):
            if char == s[i]:
                loc = i
                locs.append(i+1)
                break
    return locs

def common_spliced_motifs(s1,s2):
    M = len(s1)*[len(s2)*[0]]
    for i in range(len(s1)):
        for j in range(len(s2)):
            if s1[i] == s2[j]:
                M[i][j] = 1
    s = ''
    index = -1
    for i in range(len(s1)):
        for j in range(index + 1, len(s2)):
            if M[i][j] == 1:
                s+=s2[j]
                index = j
                break
    return s

def lcsdp(s1,s2):
    M = (len(s1))*[(len(s2))*['']]
    for i in range(len(s1)):
        for j in range(len(s2)):
            if i == 0:
                M[i][j] == ''
            if j == 0:
                M[i][j] == ''
            elif i> 0 and j > 0 and s1[i-1] == s2[j-1]:
                M[i][j] = M[i-1][j-1] + s1[i]
            else:
                M[i][j] = M[i-1][j] if len(M[i-1][j]) >= len(M[i][j-1]) else M[i][j-1]
    return M

def lcsrc(s1,s2):
    if len(s1) == 0 or len(s2) == 0:
        return ''
    elif s1[-1] == s2[-1]:
        return lcsrc(s1[:-1],s2[:-1]) + s1[-1]
    else:
        return lcsrc(s1[:-1],s2) if len(lcsrc(s1[:-1],s2)) >= len(lcsrc(s1,s2[:-1])) else lcsrc(s1,s2[:-1])
        
def LCSubStr(str1, str2, N, M):
 
    LCSuff = [['' for k in range(M+1)] for l in range(N+1)]
    for i in range(N + 1):
        for j in range(M + 1):
            if (i == 0 or j == 0):
                LCSuff[i][j] = ''
            elif (str1[i-1] == str2[j-1]):
                LCSuff[i][j] = LCSuff[i-1][j-1] + str1[i-1]
            else:
                if len(LCSuff[i-1][j]) > len(LCSuff[i][j-1]):
                    LCSuff[i][j] = LCSuff[i-1][j]
                else:
                    LCSuff[i][j] = LCSuff[i][j-1]
    return LCSuff[-1][-1]   

def edit_distance(s1, s2):
    m=len(s1)+1
    n=len(s2)+1

    tbl = {}
    for i in range(m): tbl[i,0]=i
    for j in range(n): tbl[0,j]=j
    for i in range(1, m):
        for j in range(1, n):
            cost = 0 if s1[i-1] == s2[j-1] else 1
            tbl[i,j] = min(tbl[i, j-1]+1, tbl[i-1, j]+1, tbl[i-1, j-1]+cost)

    return tbl[i,j]

def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]



def superstring(snippets, k):
    while len(snippets) > 1:
        gluestrings(snippets, k)
    return snippets[0]

def count_transitions(s1,s2):
    transitions = 0
    for i in range(len(s1)):
        if set((s1[i],s2[i])) in [{'A','G'}, {'C', 'T'}]:
            transitions += 1
    return transitions

def count_transversions(s1,s2):
    transversions = 0
    for i in range(len(s1)):
        if set((s1[i],s2[i])) in [{'A','C'}, {'A','T'}, {'G','C'}, {'G','T'}]:
            transversions += 1
    return transversions

def transition_transversion_ratio(s1,s2):
    return count_transitions(s1,s2)/count_transversions(s1,s2)


def mortal_fibbonacci_rabbits(n,m):
    rabbits = n * [0]
    rabbits[0] = 1
    for i in range(1,n):
        rabbits[i] = rabbits[i-1] + rabbits[i-2] -(rabbits[i-m]) 
    return rabbits

from itertools import permutations

def get_permutations(k):
    file1 = open("myfile.txt", "w")  # write mode
    bag = list(permutations([a + 1 for a in range(k)], k))
    file1.write(str(len(list(bag)))+ '\n')
    for item in bag:
        for char in item:
            file1.write(str(char) + ' ')
        file1.write('\n')
    file1.close()
    
def write_kmp(P):
    file1 = open('kmp.txt', 'w')
    file1.write(str(P))
    file1.close()
    
def dominant_allele_prob(k,n,m):
    l = k + n + m
    return (k*(k-1) +2*k *(n +m) + n*(n-1)*3/4 + n*m)/(l*(l-1))

codon_table = {"UUU": "F", "CUU":"L", "AUU":"I", "GUU":"V",
'UUC': 'F',      'CUC': 'L',      'AUC': 'I',      'GUC': 'V',
'UUA': 'L',      'CUA': 'L',      'AUA': 'I',      'GUA': 'V',
'UUG': 'L',      'CUG': 'L',      'AUG': 'M',      'GUG': 'V',
'UCU': 'S',      'CCU': 'P',      'ACU': 'T',      'GCU': 'A',
'UCC': 'S',      'CCC': 'P',      'ACC': 'T',      'GCC': 'A',
'UCA': 'S',      'CCA': 'P',      'ACA': 'T',      'GCA': 'A',
'UCG': 'S',      'CCG': 'P',      'ACG': 'T',      'GCG': 'A',
'UAU': 'Y',      'CAU': 'H',      'AAU': 'N',      'GAU': 'D',
'UAC': 'Y',      'CAC': 'H',      'AAC': 'N',      'GAC': 'D',
'UAA': 'Stop',   'CAA': 'Q',      'AAA': 'K',      'GAA': 'E',
'UAG': 'Stop',   'CAG': 'Q',      'AAG': 'K',      'GAG': 'E',
'UGU': 'C',      'CGU': 'R',      'AGU': 'S',      'GGU': 'G',
'UGC': 'C',      'CGC': 'R',      'AGC': 'S',      'GGC': 'G',
'UGA': 'Stop',   'CGA': 'R',      'AGA': 'R',      'GGA': 'G',
'UGG': 'W',      'CGG': 'R',      'AGG': 'R',      'GGG': 'G' }

aa_table = {char:[] for char in codon_table.values()}
for char in codon_table.keys():
    aa_table[codon_table[char]].append(char)
    
def RNA_from_protein(RNA):
    k = 1
    for char in RNA:
        k = k*len(aa_table[char]) % 1000000
    return k

aa_mass = {'A': 71.03711,
'C': 103.00919,
'D':   115.02694,
'E':   129.04259,
'F':   147.06841,
'G':  57.02146,
'H':   137.05891,
'I':   113.08406,
'K':   128.09496,
'L':   113.08406,
'M':   131.04049,
'N':   114.04293,
'P':   97.05276,
'Q':   128.05858,
'R':   156.10111,
'S':   87.03203,
'T':   101.04768,
'V':   99.06841,
'W':   186.07931,
'Y':   163.06333}

def protein_mass(protein):
    k = 0
    for char in protein:
        k += aa_mass[char]
    return k


def RNA_to_protein(s):
    j = 0
    p = ''
    while 3*j+3 <len(s):
        if codon_table[s[3*j:3*j + 3]] == 'Stop':
            return p
        else:
            p+=codon_table[s[3*j:3*j + 3]]
        j+=1

def find_motif(s,t):
    i = 0
    loc = []
    while i + len(t)<=len(s):
        if s[i:i+len(t)] == t:
            loc.append(i+1)
        i+=1
    return loc    
        
def profilefrommotifs(MotifList):
    l = len(MotifList[0])
    Profile = {'A':l*[0], 'C':l*[0], 'G':l*[0], 'T':l*[0]}
    for i in range(l):
        A = 0
        C = 0
        G = 0
        T = 0
        for j in range(len(MotifList)):
            if MotifList[j][i] == 'A':
                A+=1
            elif MotifList[j][i] == 'C':
                C+=1
            elif MotifList[j][i] == 'G':
                G+=1
            else:
                T+=1
        Profile['A'][i] = A
        Profile['C'][i] = C
        Profile['G'][i] = G
        Profile['T'][i] = T
    return Profile

def consensus_from_profile(Profile):
    consensus = ''
    for k in range(len(Profile['A'])):
        m = max([Profile[char][k] for char in ['A','C','G','T']])
        for char in ['A','C','G','T']:
            if Profile[char][k] == m:
                consensus += char
                break
    return consensus

def write_to_text_consensus_profile(motifs):
    file = open("consensus_profile.txt", "w")
    profile = profilefrommotifs(motifs)
    consensus = consensus_from_profile(profile)
    file.write(consensus + '\n')
    for char in ['A', 'C', 'G', 'T']:
        file.write(char + ': ')
        for num in profile[char]:
            file.write(str(num) + ' ')
        file.write('\n')
    file.close()

def remove_introns(DNA, introns):
    for intron in introns:
        DNA = DNA.replace(intron, '')
    return DNA

def running_count_RNA(RNA):
    count = {'A':0, 'U':0, 'C':0, 'G':0}
    for char in RNA:
        count[char] +=1
    return count

complement = {'A':'U', 'U':'A', 'C':'G', 'G':'C'}

def base_parity(RNA):
    char = complement[RNA[0]]
    endpoints = []
    for i in range(len(RNA)):
        if RNA[i] == char:
            count = running_count_RNA(RNA[:i+1])
            if count['A'] == count['U'] and count['C'] == count['G']:
                endpoints.append(i)
    return endpoints

def non_crossing_perfect_matchings(RNA):
    if len(RNA) == 0:
        return 1
    count = running_count_RNA(RNA)
    if not(count['A'] == count['U'] and count['C'] == count['G']):
        return 0
    endpoints = base_parity(RNA)
    k = 0
    if len(RNA) == 2 and 1 in endpoints:
        return 1
    for pt in endpoints:
        k += non_crossing_perfect_matchings(RNA[1:pt]) * non_crossing_perfect_matchings(RNA[pt +1:])
    return k

def dp_non_crossing(RNA):
    partials = [len(RNA) * [0] for _ in RNA]
    for i in range(len(RNA)-1):
        if set([RNA[i], RNA[i+1]]) == {"A", "U"} or set([RNA[i], RNA[i+1]]) == {"C", "G"}:
            partials[i][i+1] = 1
        partials[i+1][i] = 1
    for k in range(3,len(RNA)):
        for i in range(len(RNA)-k):
            count = running_count_RNA(RNA[i:i+k+1])
            if count['A'] == count['U'] and count['C'] == count['G']:
                for j in range(1,k):
                    if complement[RNA[i]] == RNA[i+j]:
                        partials[i][i+k] += partials[i+1][i+j-1]*partials[i+j+1][i+k]
                if complement[RNA[i]] == RNA[i+k]:
                    partials[i][i+k] += partials[i+1][i+k-1]
    return partials[0][-1]

def dp_motzkin(RNA):
    partials = [len(RNA) * [0] for _ in RNA]
    for i in range(len(RNA)-1):
        if set([RNA[i], RNA[i+1]]) == {"A", "U"} or set([RNA[i], RNA[i+1]]) == {"C", "G"}:
            partials[i][i+1] = 2
        partials[i+1][i] = 1
        partials[i][i] = 1
        partials[-1][-1] = 1
    for k in range(2,len(RNA)):
        for i in range(len(RNA)-k):
            if complement[RNA[i]]==RNA[k]:
                for j in range(1,k):
                    if complement[RNA[i]] == RNA[i+j]:
                        partials[i][i+k] += partials[i+1][i+j-1]*partials[i+j+1][i+k]
                if complement[RNA[i]] == RNA[i+k]:
                    partials[i][i+k] += partials[i+1][i+k-1]
    return partials[0][-1]

def catalan_pairs(lis):
    pairs = []
    k = len(lis)
    if len(lis) == 0:
        return [[]]
    elif len(lis) == 2:
        return [[(lis[0],lis[1])]]
    else:
        for i in range(1, k//2 + 1):
            pairs += [[(lis[0], lis[2*i-1])] + pairing1 + pairing2 for pairing1 in catalan_pairs(lis[1:2*i-1]) for pairing2 in catalan_pairs(lis[2*i:])]
    return pairs

# def six_subsets(n,subset1,subset2):
#     union = set()
#     for a in subset1:
#         union.add(a)
#     for b in subset2:
#         union.add(b)
#     intersection = set()
#     for a in subset1:
#         if a in subset2:
#             intersection.add(a)
#     AminusB = set()
#     pass

def ReverseComplement(pattern):
    p = len(pattern)
    patternrc=''
    for i in range(p):
        if pattern[p-1-i] == 'A':
            patternrc+='T'
        elif pattern[p-1-i] == 'C':
            patternrc+='G'
        elif pattern[p-1-i] == 'G':
            patternrc+='C'
        else:
            patternrc+='A'
    return patternrc

def find_palindromes(DNA):
    palindromes = []
    for i in range(len(DNA)):
        for j in [4, 6, 8, 10, 12]:
            if i+j<=len(DNA) and DNA[i:i+j] == reverse_complement(DNA[i:i+j]):
                palindromes.append([i+1,j])
    return palindromes

def orfs(strand):
    proteins = []
    s = transcribe_rna(strand)
    for i in range(len(s)-2):
        if s[i:i+3] == "AUG":
            proteins.append(RNA_to_protein(s[i:]))
    t = transcribe_rna(reverse_complement(strand))
    for i in range(len(t) -2):
        if t[i:i+3] == "AUG":
            proteins.append(RNA_to_protein(t[i:]))
    proteins = [protein for protein in set(proteins) if protein]
    return proteins
           
def read_in_s_A(file_name):
    file = open(file_name, 'r')
    s, A = file.readlines()
    return s, A

def probs(s, A):
    probs = []
    for x in A:
        prob = 0
        nucleo_probs = {'A':(1-x)/2, 'C':x/2, 'G':x/2, 'T':(1-x)/2}
        for char in s:
            prob += log10(nucleo_probs[char])
        probs.append(prob)
    return probs
            
def printLIS(arr: list):
    for x in arr:
        print(-x, end=" ")
    print()
 
# Function to construct and print Longest Increasing
# Subsequence
def constructPrintLIS(arr: list, n: int):
 
    # L[i] - The longest increasing sub-sequence
    # ends with arr[i]
    l = [[] for i in range(n)]
 
    # L[0] is equal to arr[0]
    l[0].append(arr[0])
 
    # start from index 1
    for i in range(1, n):
 
        # do for every j less than i
        for j in range(i):
 
            # L[i] = {Max(L[j])} + arr[i]
            # where j < i and arr[j] < arr[i]
            if arr[i] > arr[j] and (len(l[i]) < len(l[j]) + 1):
                l[i] = l[j].copy()
 
        # L[i] ends with arr[i]
        l[i].append(arr[i])
 
    # L[i] now stores increasing sub-sequence of
    # arr[0..i] that ends with arr[i]
    maxx = l[0]
 
    # LIS will be max of all increasing sub-
    # sequences of arr
    for x in l:
        if len(x) > len(maxx):
            maxx = x
 
    # max will contain LIS
    printLIS(maxx)
    
def failure_array(pattern):
    s = ""
    file1 = open('kmp.txt', 'w')
    P = len(pattern)*[0]
    for j in range(1, len(pattern)):
        if pattern[j] == pattern[0]:
            k = j
            while k < len(pattern) and pattern[:k-j] == pattern[j:k]:
                P[k] = max(P[k], k-j + 1)
                s += str(P[k]) + ' '
                k+=1
    file1.close()
    return P

def fail_array(pattern):
    P = len(pattern)*[0]
    for j in range(1,len(pattern)):
        if pattern[j] == pattern[0]:
            P[j] = max(P[j],1)
            k = j+1
            while k <= len(pattern) and pattern[:k-j] == pattern[j:k]:
                P[k-1] = max(P[k-1], k-j)
                k+=1
        
    return P, str(P).strip('[').strip(']').replace(',','')


fourmers = [char1 + char2 + char3 + char4 for char1 in ['A','C','G','T'] for char2 in ['A','C','G','T'] for char3 in ['A','C','G','T'] for char4 in ['A','C','G','T']]

def fourmer_comp(pattern):
    comp = {fourmer:0 for fourmer in fourmers}
    for i in range(len(pattern)-3):
        comp[pattern[i:i+4]] += 1
    return comp

def distance_matrix_p(strands):
    matrix = [len(strands)*[0] for _ in strands]
    for i in range(len(strands)):
        for j in range(i,len(strands)):
            d = hamming_distance(strands[i], strands[j])/len(strands[0])
            matrix[i][j] = d
            matrix[j][i] = d
    return matrix

def distance_matrix(strands):
    matrix = [len(strands)*[0] for _ in strands]
    for i in range(len(strands)):
        for j in range(i,len(strands)):
            d = hamming_distance(strands[i], strands[j])
            matrix[i][j] = d
            matrix[j][i] = d
    return matrix

def corr(strands):
    for i in range(len(strands)):
        strand, count1 = strands[i], 0
        for j in range(len(strands)):
            count1 += (strand == strands[j] or reverse_complement(strand) == strands[j]) and i != j
        if count1 == 0:
            for j in range(len(strands)):
                if hamming_distance(strand, strands[j]) == 1:
                    print(strand + '->' + strands[j])
                    break
                if hamming_distance(strand, reverse_complement(strands[j])) ==1:
                    print(strand + '->' + reverse_complement(strands[j]))
                    break

def write_matrix(matrix):
    file1 = open('my_matrix.txt', 'w')
    for row in matrix:
        for entry in row:
            file1.write(str(entry) + ' ')
        file1.write('\n')
    file1.close()
    
def make_language(alph):
    language = []
    for a in alph:
        language.append(a)
        for b in alph:
            language.append(a + b)
            for c in alph: 
                language.append(a+b+c)
    return language
    
        
def pow_set(n):
    pow_set = [[0],[1]]
    for i in range(n-1):
        pow_set = [ps + s for ps in pow_set for s in [[0],[1]]]
    return pow_set

def substrings(string):
    ps = pow_set(len(string))
    substrings = set()
    for p in ps:
        s = ''
        for c in range(len(p)):
            if p[c] == 1:
                s+=string[c]
        substrings.add(s)
    return sorted(list(substrings)[1:])
    

                    
            
            
    
            
        
    
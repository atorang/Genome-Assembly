def Composition(Text, k):
    k_mer=[]
    for i in range(len(Text)-k+1):
        k_mer.append(Text[i:i+k])
    return sorted(k_mer)
#########################################################
def StringSpelledByGenomePath(seq):
    k=len(seq[0])
    string=seq[0]
    for i in range(len(seq)-1):
        string=string+seq[i+1][k-1]
    return string
#########################################################
def OverlapGraph(seq):
    k = len(seq[0])
    set={}
    for i in range(len(seq)):
        n=0
        for j in range(len(seq)):
            if seq[i][1:]==seq[j][:k-1] and i!=j:
                if n==1:
                    set[seq[i]].append(seq[j])
                else:
                    set[seq[i]]=[seq[j]]
                    n=1
    return set
#########################################################
def DeBruijn(Text,k):
    path={}
    for i in range(len(Text)-k+1):
        k_mer=Text[i:i+k]
        for j in path.keys():
            if k_mer[:k-1]==j:
                path[k_mer[:k-1]].append(k_mer[1:])
                break
        else:
            path[k_mer[:k-1]]=[k_mer[1:]]
    return path
#########################################################
def DeBruijnGraphFromKmers(seq):
    graph = {}
    k=len(seq[0])
    for i in range(len(seq)):
        k_mer = seq[i]
        for j in graph.keys():
            if k_mer[:k - 1] == j:
                graph[k_mer[:k - 1]].append(k_mer[1:])
                break
        else:
            graph[k_mer[:k - 1]] = [k_mer[1:]]
    return graph
#########################################################
import random
def EulerianCycle(Graph):
    start = []
    end = []
    for i in Graph:
        i = i.rstrip()
        x = i.split(' -> ')
        start.extend([x[0]] * len(x[1].split(',')))
        end.extend(x[1].split(','))
    edge=list(range(len(end)))
    unused_edge=edge.copy()
    n=random.randint(0,len(end)-1)
    used_edges=[n]
    unused_edge.remove(n)
    cycle=start[n]+'->'+end[n]
    j=0
    while j<len(unused_edge):
        i=unused_edge[j]
        if end[n]==start[i]:
            cycle=cycle+'->'+end[i]
            used_edges.append(i)
            unused_edge.remove(i)
            n=i
            j=0
        else: j+=1

    while len(unused_edge)!= 0:
        for i in used_edges:
            if start[i] in [start[j] for j in unused_edge]:
                n=i
                break
        cycle=start[n]+'->'+end[n]
        path=[]
        for i in range(len(used_edges)):
            if used_edges[i]==n:
                path = used_edges[i+1:]
                path.extend(used_edges[:i])
        for i in path:
            cycle+='->'+end[i]
        used_edges=[n]
        used_edges.extend(path)
        n=used_edges[-1]
        j = 0
        while j < len(unused_edge):
            i = unused_edge[j]
            if end[n] == start[i]:
                cycle = cycle + '->' + end[i]
                used_edges.append(i)
                unused_edge.remove(i)
                n = i
                j = 0
            else:
                j += 1
    return cycle
#########################################################
def EulerianPath(Graph):
    start = []
    end = []
    for i in Graph:
        i = i.rstrip()
        x = i.split(' -> ')
        start.extend([x[0]] * len(x[1].split(',')))
        end.extend(x[1].split(','))
    End_node,Start_node='',''
    for i in range(len(end)):
        if start.count(start[i])>end.count(start[i]):
            Start_node=start[i]
        elif end.count(end[i])>start.count(end[i]):
            End_node=end[i]
    start.append(End_node)
    end.append(Start_node)
    edge=list(range(len(end)))
    unused_edge=edge.copy()
    n=random.randint(0,len(end)-1)
    used_edges=[n]
    unused_edge.remove(n)
    j=0
    while j<len(unused_edge):
        i=unused_edge[j]
        if end[n]==start[i]:
            used_edges.append(i)
            unused_edge.remove(i)
            n=i
            j=0
        else: j+=1
    while len(unused_edge)!= 0:
        for i in used_edges:
            if start[i] in [start[j] for j in unused_edge]:
                n=i
                break
        path=[]
        for i in range(len(used_edges)):
            if used_edges[i]==n:
                path = used_edges[i+1:]
                path.extend(used_edges[:i])
        used_edges=[n]
        used_edges.extend(path)
        n=used_edges[-1]
        j = 0
        while j < len(unused_edge):
            i = unused_edge[j]
            if end[n] == start[i]:
                used_edges.append(i)
                unused_edge.remove(i)
                n = i
                j = 0
            else:
                j += 1
    path=[]
    for i in range(len(end)):
        if used_edges[i]==len(end)-1:
            path=used_edges[i+1:]
            path.extend(used_edges[:i])
    Path=Start_node+'->'+end[path[0]]
    for i in range(1,len(path)):
        Path+='->'+end[path[i]]
    return Path
#########################################################
def StringReconstruction(sample):
    x=[]
    for i in sample:
        x.append(i.rstrip())
    k=int(x[0])
    string=x[1:]
    graph=DeBruijnGraphFromKmers(string)
    x=[]
    for i in graph.keys():
        x.append(i + ' -> ' + ','.join(graph[i]) + '\n')
    eulerian_path=EulerianPath(x).split('->')
    text=eulerian_path[0]
    for i in range(1,len(eulerian_path)):
        text+=eulerian_path[i][-1]
    return text
#########################################################
import itertools
def KUniversalCircularString(k):
    x=list(itertools.product('01', repeat=k))
    string=[]
    for i in x:
        i=list(i)
        string.append(''.join(i))
    graph = DeBruijnGraphFromKmers(string)
    x = []
    for i in graph.keys():
        x.append(i + ' -> ' + ','.join(graph[i]) + '\n')
    eulerian_path = EulerianCycle(x).split('->')
    text = eulerian_path[0]
    for i in range(1,len(eulerian_path)):
        text+=eulerian_path[i][-1]
    return text[:-k+1]
#########################################################
def StringReconstructionFromReadPairs(sample):
    x = []
    for i in sample:
        x.append(i.rstrip())
    k = int(x[0].split()[0])
    d = int(x[0].split()[1])
    string = x[1:]
    graph = {}
    for i in range(len(string)):
        PairKmer = string[i]
        for j in graph.keys():
            if PairKmer[:k - 1]+PairKmer[k+1:2*k] == j:
                graph[PairKmer[:k - 1]+PairKmer[k+1:2*k]].append(PairKmer[1:k]+PairKmer[k+2:])
                break
        else:
            graph[PairKmer[:k - 1]+PairKmer[k+1:2*k]] = [PairKmer[1:k]+PairKmer[k+2:]]
    Start=[]
    End=[]
    for i in graph.keys():
        Start.extend([i]* len(graph[i]))
        End.extend(graph[i])
    End_node, Start_node = '', ''
    for i in range(len(End)):
        if Start.count(Start[i]) > End.count(Start[i]):
            Start_node = Start[i]
        elif End.count(End[i]) > Start.count(End[i]):
            End_node = End[i]
    Start.append(End_node)
    End.append(Start_node)
    edge = list(range(len(End)))
    unused_edge = edge.copy()
    n = random.randint(0, len(End) - 1)
    used_edges = [n]
    unused_edge.remove(n)
    j = 0
    while j < len(unused_edge):
        i = unused_edge[j]
        if End[n] == Start[i]:
            used_edges.append(i)
            unused_edge.remove(i)
            n = i
            j = 0
        else:
            j += 1
    while len(unused_edge) != 0:
        for i in used_edges:
            if Start[i] in [Start[j] for j in unused_edge]:
                n = i
                break
        path = []
        for i in range(len(used_edges)):
            if used_edges[i] == n:
                path = used_edges[i + 1:]
                path.extend(used_edges[:i])
        used_edges = [n]
        used_edges.extend(path)
        n = used_edges[-1]
        j = 0
        while j < len(unused_edge):
            i = unused_edge[j]
            if End[n] == Start[i]:
                used_edges.append(i)
                unused_edge.remove(i)
                n = i
                j = 0
            else:
                j += 1
    path = []
    for i in range(len(End)):
        if used_edges[i] == len(End) - 1:
            path = used_edges[i + 1:]
            path.extend(used_edges[:i])
    Path = Start_node[:k - 1] + End[path[0]][k - 2]
    for i in range(1, len(path)):
        Path += End[path[i]][k - 2]
    for i in range(-k-d,0):
        Path+=End[path[i]][2*k-3]
    return Path
#########################################################
def MaximalNonBranchingPaths(Graph):
    Paths=[]
    Start = []
    End = []
    for i in Graph.keys():
        Start.extend([i] * len(Graph[i]))
        End.extend(Graph[i])
    nodes=set(Start)| set(End)
    for i in nodes:
        if Start.count(i)!=1 or End.count(i)!=1:
            if Start.count(i) > 0:
                v = Start.index(i)
                j=v
                while Start[v]==Start[j]:
                    NonBranchingPath=Start[j]+End[j][-1]
                    w=j
                    while Start.count(End[w]) == End.count(End[w]) and Start.count(End[w])==1:
                        n=Start.index(End[w])
                        NonBranchingPath+= End[n][-1]
                        w=n
                    Paths.append(NonBranchingPath)
                    j+=1
                    if j>len(Start)-1:
                        break
    return Paths
#########################################################
def ProteinTranslation(string):
    codon={'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
           'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGU': 'S', 'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I',
           'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAU': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
           'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R', 'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
           'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAU': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
           'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G', 'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
           'UAA': '', 'UAC': 'Y', 'UAG': '', 'UAU': 'Y', 'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
           'UGA': '', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C', 'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F'}
    protein=''
    for i in range(int(len(string)/3)):
        protein+=codon[string[3*i:3*i+3]]
    return protein
#########################################################
def ReverseComplementRNA(Text):
    "Generate Reverse Complement of the string"
    RC=''
    Complement={'A':'U','C':'G','G':'C','U':'A'}
    Text=reversed(Text)
    for i in Text:
        RC=RC+Complement[i]
    return RC
#########################################################
def PeptideEncoding(string,peptide):
    length=len(peptide)
    DNA=string
    string=''
    for i in DNA:
        if i=='T':
            string+='U'
        else:
            string+=i
    i=0
    substring=[]
    while i+3*length<=len(string):
        pattern=string[i:i+3*length]
        if ProteinTranslation(pattern)==peptide:
            substring.append(DNA[i:i+3*length])
        if  ProteinTranslation(ReverseComplementRNA(pattern))==peptide:
            substring.append(DNA[i:i+3*length])
        i += 1
    return substring
#########################################################
def LinearSpectrum(Peptide):
    AminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    AminoAcidMass = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]
    PrefixMass = [0]
    Peptide=' '+Peptide
    for i in range(1,len(Peptide)):
        for j in range(20):
            if AminoAcid[j] == Peptide[i]:
                PrefixMass.append(PrefixMass[i-1] + AminoAcidMass[j])
    LinearSpectrum=[0]
    for i in range(len(Peptide)):
        for j in range(i+1,len(Peptide)):
            LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
    LinearSpectrum.sort()
    return LinearSpectrum
#########################################################
def CyclioSpectrum(Peptide):
    AminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    AminoAcidMass = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]
    PrefixMass=[0]
    Peptide = ' ' + Peptide
    for i in range(1, len(Peptide)):
        for j in range(20):
            if AminoAcid[j] == Peptide[i]:
                PrefixMass.append(PrefixMass[i - 1] + AminoAcidMass[j])
    peptideMass=PrefixMass[len(Peptide)-1]
    CyclicSpectrum=[0]
    for i in range(len(Peptide)-1):
        for j in range(i+1,len(Peptide)):
            CyclicSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i > 0 and j < len(Peptide)-1:
                CyclicSpectrum.append(peptideMass - (PrefixMass[j]-PrefixMass[i]))
    CyclicSpectrum.sort()
    return CyclicSpectrum
#########################################################
def Expand(peptides,x):
    AminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    expanded=[]
    if x == 1:
        for i in peptides:
            for j in AminoAcid:
                expanded.append(i+j)
    else:
        if peptides==['']:
            for j in x:
                expanded.append(str(j))
        else:
            for i in peptides:
                for j in x:
                    expanded.append(i+'-'+str(j))
    return expanded
#########################################################
def PeptideMass(peptide):
    AminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    AminoAcidMass = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]
    Mass=0
    for i in peptide:
        for j in range(len(AminoAcid)):
            if i==AminoAcid[j]:
                Mass+=AminoAcidMass[j]
    return Mass
#########################################################
def Consistent(peptide,spectrum):
    peptidespectrum=LinearSpectrum(peptide)
    m=0
    for i in range(len(peptidespectrum)):
        for j in range(len(spectrum)):
            if peptidespectrum[i]==spectrum[j]:
                m+=1
                del spectrum[j]
                break
    if m==len(peptidespectrum):
        return True
    else:
        return False
#########################################################
def PeptideMassShow(peptide):
    AminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    AminoAcidMass = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]
    MassShow=''
    for i in peptide:
        for j in range(len(AminoAcid)):
            if i==AminoAcid[j]:
                MassShow+='-'+str(AminoAcidMass[j])
    return MassShow[1:]
#########################################################
def CyclopeptideSequencing(Spectrum):
    Peptides=['']
    result=[]
    x=1
    n=1
    while len(Peptides)!=0:
        Peptides=Expand(Peptides,x)
        if n==1:
            x=Peptides
            n = 0
        i=0
        while i < len(Peptides):
            peptide=Peptides[i]
            if PeptideMass(peptide) == max(Spectrum):
                if CyclioSpectrum(peptide) == Spectrum:
                    result.append(PeptideMassShow(peptide))
                Peptides.remove(peptide)
            elif Consistent(peptide,Spectrum)==False:
                Peptides.remove(peptide)
            else:
                i+=1
    return result
#########################################################
def CyclopeptideScoring(peptide,spectrum):
    peptidespectrum=CyclioSpectrum(peptide)
    m=0
    Spec = spectrum.copy()
    for i in range(len(peptidespectrum)):
        for j in range(len(Spec)):
            if peptidespectrum[i]==Spec[j]:
                m+=1
                del Spec[j]
                break
    return m
#########################################################
def LinearpeptideScoring(peptide,spectrum):
    peptidespectrum=LinearSpectrum(peptide)
    m=0
    Spec=spectrum.copy()
    for i in range(len(peptidespectrum)):
        for j in range(len(Spec)):
            if peptidespectrum[i]==Spec[j]:
                m+=1
                del Spec[j]
                break
    return m
#########################################################
import numpy
def Trim(Leaderboard, spectrum, N):
    Score=[]
    if len(Leaderboard)>=N:
        for i in Leaderboard:
            Score.append(LinearpeptideScoring(i,spectrum))
        ind=list(reversed(numpy.argsort(Score)))
        n_th=Score[ind[N-1]]
        peptides=[Leaderboard[i] for i in ind[:N]]
        i=N

        while i<len(Score):
            if n_th==Score[ind[i]]:
                peptides.append(Leaderboard[ind[i]])
                i+=1
            else: break
    else:
        peptides=Leaderboard
    return peptides
#########################################################
def LeaderboardCyclopeptideSequencing(Spectrum, N):
    Leaderboard=['']
    LeaderPeptide=''
    Max=max(Spectrum)
    while len(Leaderboard)!=0:
        Leaderboard=Expand(Leaderboard,1)
        i = 0
        while i < len(Leaderboard):
            peptide = Leaderboard[i]
            Mass=PeptideMass(peptide)
            if Mass == Max:
                if CyclopeptideScoring(peptide, Spectrum) > CyclopeptideScoring(LeaderPeptide, Spectrum):
                    LeaderPeptide=peptide
                i += 1
            elif Mass > Max:
                Leaderboard.remove(peptide)
            else:
                i+=1
        Leaderboard=Trim(Leaderboard, Spectrum, N)
        #print(Leaderboard, len(Leaderboard))
    return PeptideMassShow(LeaderPeptide)
#########################################################
def SpectralConvolution(spectrum):
    conv=[]
    for i in spectrum:
        for j in spectrum:
            if i>j:
                conv.append(i-j)
    return conv
#########################################################
def MostFrequentConvolution(conv,M):
    AC={}
    for i in set(conv):
        if i>=57 and i<=200:
            AC[i]=conv.count(i)
    Repeat=list(AC.values())
    Repeat.sort(reverse=True)
    m_th = Repeat[M - 1]
    items = []
    for i in AC.keys():
        if m_th <= AC[i]:
            items.append(i)
    return items
#########################################################
def MassCyclioSpectrum(Peptide):
    PrefixMass=[0]
    peptide = [0]
    peptide.extend(Peptide)
    for i in range(1,len(peptide)):
        PrefixMass.append(PrefixMass[i - 1] + int(peptide[i]))
    peptideMass=PrefixMass[len(peptide)-1]
    CyclicSpectrum=[0]
    for i in range(len(peptide)-1):
        for j in range(i+1,len(peptide)):
            CyclicSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i > 0 and j < len(peptide)-1:
                CyclicSpectrum.append(peptideMass - (PrefixMass[j]-PrefixMass[i]))
    CyclicSpectrum.sort()
    return CyclicSpectrum
#########################################################
def MassLinearSpectrum(Peptide):
    Peptide=Peptide.split('-')
    PrefixMass = [0]
    peptide=[0]
    peptide.extend(Peptide)
    for i in range(1,len(peptide)):
        PrefixMass.append(PrefixMass[i - 1] + int(peptide[i]))
    LinearSpectrum = [0]
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide)):
            LinearSpectrum.append(PrefixMass[j] - PrefixMass[i])
    LinearSpectrum.sort()
    return LinearSpectrum
#########################################################
def MassCyclopeptideScoring(peptide,spectrum):
    peptidespectrum=MassCyclioSpectrum(peptide)
    m=0
    Spec = spectrum.copy()
    for i in range(len(peptidespectrum)):
        for j in range(len(Spec)):
            if peptidespectrum[i]==Spec[j]:
                m+=1
                del Spec[j]
                break
    return m
#########################################################
def MassLinearpeptideScoring(peptide, spectrum):
    peptidespectrum = MassLinearSpectrum(peptide)
    m = 0
    Spec = spectrum.copy()
    for i in range(len(peptidespectrum)):
        for j in range(len(Spec)):
            if peptidespectrum[i] == Spec[j]:
                m += 1
                del Spec[j]
                break
    return m
#########################################################
def MassTrim(Leaderboard, spectrum, N):
    Score = []
    if len(Leaderboard) >= N:
        for i in Leaderboard:
            Score.append(MassLinearpeptideScoring(i, spectrum))
        ind = list(reversed(numpy.argsort(Score)))
        n_th = Score[ind[N - 1]]
        peptides = [Leaderboard[i] for i in ind[:N]]
        i = N
        while i < len(Score):
            if n_th == Score[ind[i]]:
                peptides.append(Leaderboard[ind[i]])
                i += 1
            else:
                break
        print(max(Score))
    else:
        peptides = Leaderboard
    return peptides
#########################################################
def ConvolutionCyclopeptideSequencing(M,N,Spectrum):
    Spectrum.sort()
    Leaderboard = ['']
    LeaderPeptide = ['0']
    Max = max(Spectrum)
    conv= SpectralConvolution(Spectrum)
    items=MostFrequentConvolution(conv,M)
    while len(Leaderboard) != 0:
        Leaderboard = Expand(Leaderboard,items)
        i = 0
        while i < len(Leaderboard):
            peptide = Leaderboard[i].split('-')
            Mass=0
            for j in peptide:
                Mass += int(j)
            if Mass == Max:
                if MassCyclopeptideScoring(peptide, Spectrum) > MassCyclopeptideScoring(LeaderPeptide, Spectrum):
                    LeaderPeptide = peptide
                i += 1
            elif Mass > Max:
                Leaderboard.remove('-'.join(peptide))
            else:
                i += 1
        Leaderboard = MassTrim(Leaderboard, Spectrum, N)
        print(Leaderboard, len(Leaderboard))
    return LeaderPeptide
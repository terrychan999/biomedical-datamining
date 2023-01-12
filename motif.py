import os
import sys

# 讀取fasta
def parseFileToFasta(filePath):
    fasta = dict()
    file = open(filePath, "r")
    fileLines = file.readlines()
    file.close()
    motif_length = int(fileLines[0].split()[-1])
    for line in fileLines[1:]:
        if line.strip() != "":
            if line.startswith(">"):
                fasta_id = line.lstrip(">").rstrip()
                fasta[fasta_id] = ""
            else:
                fasta[fasta_id] += line.strip().upper()
    return motif_length, fasta

# sliding window產生所有位置的dna segment
def genMotif(motif_length, fasta):
    motifDict = dict()
    for fasta_id, sequence in fasta.items():
        motifDict[fasta_id] = list()
        for index in range(0, len(sequence)-motif_length+1):
            motifDict[fasta_id].append(sequence[index: index+motif_length])
    return motifDict

# hamming distance
def diff(seq1, seq2):
    score = 0
    for index in range(motif_length):
        if seq1[index] != seq2[index]:
            score += 1
    return score

# 以第一個sequence為基準, 將sequence 1的所有dna segment與其他所有sequence比較
def genScoreTable(motifDict):
    fasta_id1, motif1 = list(motifDict.items())[0]
    mutated_motif1 = list()
    score_table = dict()
    for seq1 in motif1:
        # 考慮dna sequence 有 1 mutation base
        for mutation in ["A", "T", "C", "G"]:
            for character in range(len(seq1)):
                mutated = list(seq1)
                if mutation != seq1[character]:
                    mutated[character] = mutation
                    mutated_str = "".join(mutated)
                    mutated_motif1.append(mutated_str)
    mutated_motif1.append(seq1)

    for mutated_str in mutated_motif1:
        for fasta_id, motif in list(motifDict.items())[1:]:
            for seq in motif:
                score = diff(seq, mutated_str)
                if score <= 1:
                    if mutated_str not in score_table:
                        score_table[mutated_str] = 1
                    else:
                        score_table[mutated_str] += 1
    return score_table


if __name__ == "__main__":
    file = sys.argv[1]
    motif_length, fasta = parseFileToFasta(file)
    motifDict = genMotif(motif_length, fasta)
    score_table = genScoreTable(motifDict)
    # 出現次數最高的 segemnt便是可能的motif
    score_table = sorted(score_table.items(), key=lambda x: x[1])
    k, v = score_table[-1]
    print(f"possible motif: {k}")
    for fasta_id, motif in motifDict.items():
        for index, seq in enumerate(motif):
            score = diff(seq, k)
            if score <= 1:
                print(f">{fasta_id}\npos:{index+1} {seq}")

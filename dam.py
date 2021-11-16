from os import X_OK
import numpy as np
from numpy import *
import sys 
import getopt

substitionMatrixFile = open("scoring1.txt")
lines = substitionMatrixFile.readlines()
aminoacids= lines[20].split("\t")
pathPositions = []




def initMatrix(seq1, seq2):
    x = len(seq1)+1
    y = len(seq2)+1
    distanceMatrix = []
    for i in range(y):
        col = []
        for j in range(x):
            col.append(0)
        distanceMatrix.append(col)
    for i in range(x):
        distanceMatrix[0][i] = -i
    for y in range(y):
        distanceMatrix[y][0] = -y
    return distanceMatrix



def printMatrix(matrix):
    for y in range(len(matrix)):
        for x in range(len(matrix[0])):
            print(str(matrix[y][x]), end= "\t")
        print("")

def initMatrix0(Seq1,Seq2):
    matrix = []
    for j in range(len(Seq2)+1):
        row = []
        for  i in range(len(Seq1)+1):
            row.append(0)
        matrix.append(row)
    return matrix

    

def getOptimalGlobalPath(matrix, seq1, seq2):
    gapPositions = [[0]*len(matrix[0]) for i in range(len(matrix))]
    for j in range(1,len(matrix)):
        for i in range(1,len(matrix[0])):
            if(gapPositions[j-1][i]== 1):
                a = matrix[j-1][i] + gapExtensionPenalty
            if(gapPositions[j][i-1]==1):
                b = matrix[j][i-1] + gapExtensionPenalty
        
            if(gapPositions[j][i-1]!= 1 or gapPositions[j-1][i]!= 1):
                a = matrix[j-1][i] + gapPenalty
                b = matrix[j][i-1] + gapPenalty
            
            c = matrix[j-1][i-1] + compare(seq1[i-1],seq2[j-1])
            if(max(a,b,c)==a):
                gapPositions[j][i]=1
            elif(max(a,b,c)==b):
                gapPositions[j][i]=1
            else:
                gapPositions[j][i]=0
            matrix[j][i] = max(a,b,c)
    return gapPositions

def compare(str1,str2):
    #scoring system 
    if(str2 == "-"):
        return gapPenalty
    index1 = aminoacids.index(str1)
    index2 = aminoacids.index(str2)
    
    if(index2>index1):
    #swap values, because we have half matrix array
        index = index1
        index1 = index2
        index2 = index
    return int(lines[index1].split("\t")[index2])

    


def get_string(file):
    
    with open(file,'r') as f:
        for line in f:
            l = line.replace('\n','')
            if '>' in l:
                pass
            else:
                string = l
                
    return string 


def getLocalOptimalPath(matrix, seq1, seq2):
    MAXx = 0
    MAXy = 0
    MAXval = 0
    gapPositions = [[0]*len(matrix[0]) for i in range(len(matrix))]
    for j in range(1,len(matrix)):
        for i in range(1,len(matrix[0])):
            if(gapPositions[j-1][i]== 1):
                a = matrix[j-1][i] + gapExtensionPenalty
            if(gapPositions[j][i-1]==1):
                b = matrix[j][i-1] + gapExtensionPenalty
        
            if(gapPositions[j][i-1]!= 1 or gapPositions[j-1][i]!= 1):
                a = matrix[j-1][i] + gapPenalty
                b = matrix[j][i-1] + gapPenalty
            
            c = matrix[j-1][i-1] + compare(seq1[i-1],seq2[j-1])
            if(max(a,b,c,0)==a):
                gapPositions[j][i]=1
            elif(max(a,b,c,0)==b):
                gapPositions[j][i]=1
            else:
                gapPositions[j][i]=0
            matrix[j][i] = max(a,b,c,0)
            if(matrix[j][i]>MAXval):
                MAXval = matrix[j][i]
                MAXx = i
                MAXy = j

    return MAXx,MAXy
    



def print_result(sequenceA,sequenceB):
    print ('Smith-Waterman Algorithm (local optimal alignment) result:')
    for a in range(len(sequenceA)):
        for b in range(len(sequenceB)):           
            print ('Sequence1   ' + str(sequenceA[a]))
            print ('Sequence2   ' + str(sequenceB[b]))
    
def printIdentityPercentage(seqA, seqB):
    matchCount = 0
    for i in range(len(seqA[0])):
        if(seqA[0][i] == seqB[0][i]):
            matchCount+=1
    print("Identity percentage: " + "{:.2f}".format(matchCount*100/len(seqA[0]))+"%")
    return 

def printAlignedSequences():
    alignedSeq1 = ""
    alignedSeq2 = ""
    preI = 0
    preJ = 0
    for i,j in pathPositions:
        if(preI+1 == i and preJ+1 == j):
            alignedSeq1 += Seq1[j-1]
            alignedSeq2 += Seq2[i-1]
        elif(preI+1 == i):
            alignedSeq1 += "-"
            alignedSeq2 += Seq2[i-1]
        else:
            alignedSeq1 += Seq1[j-1]
            alignedSeq2 += "-"
        preI = i
        preJ = j
    print(alignedSeq1)
    print(alignedSeq2)

    #calculating the identity percentage
    matchCount=0
    for i in range(len(alignedSeq1)):
        if(alignedSeq1[i]==alignedSeq2[2]):
            matchCount+=1
    if(len(alignedSeq1)==0):
        return 1
    return matchCount*100 / len(alignedSeq1)
    



def traceBackStart(matrix):
    xMax = len(matrix[0])        
    yMax = len(matrix)
    score = matrix[yMax-1][xMax-1] 
    if(matrix[yMax-2][xMax-2]  == score - compare(Seq1[xMax-2],Seq2[yMax-2])):
        if(traceBack(matrix, yMax-2, xMax-2)):
            pathPositions.append([yMax-1,xMax-1])
            return True
    if(matrix[yMax-2][xMax-1]  == score - gapExtensionPenalty):
        if(traceBack(matrix, yMax-2, xMax-1)):
            pathPositions.append([yMax-1,xMax-1])
            return True
    if(matrix[yMax-1][xMax-2]  == score - gapExtensionPenalty):
        if(traceBack(matrix, yMax-1, xMax-2)):
            pathPositions.append([yMax-1,xMax-1])
            return True
    if(matrix[yMax-2][xMax-1]  == score - gapPenalty):
        if(traceBack(matrix, yMax-2, xMax-1)):
            pathPositions.append([yMax-1,xMax-1])
            return True
    if(matrix[yMax-1][xMax-2]  == score - gapPenalty):
        if(traceBack(matrix, yMax-1, xMax-2)):
            pathPositions.append([yMax-1,xMax-1])
            return True
    else:
        print("something is wrong here")       


def traceBack1(matrix, y, x):
    score = matrix[y][x]
    if(score == 0):
        pathPositions.append([y,x])
        return True
    if(matrix[y-1][x-1]  == score - compare(Seq1[x-1],Seq2[y-1])):
        if(traceBack1(matrix, y-1, x-1)):
            pathPositions.append([y,x])
            return True
    if(matrix[y-1][x]  == score - gapExtensionPenalty):
        if(traceBack1(matrix, y-1, x)):
            pathPositions.append([y,x])
            return True
    if(matrix[y][x-1]  == score - gapExtensionPenalty):
        if(traceBack1(matrix, y, x-1)):
            pathPositions.append([y,x])
            return True
    if(matrix[y-1][x]  == score - gapPenalty):
        if(traceBack1(matrix, y-1, x)):
            pathPositions.append([y,x])
            return True
    if(matrix[y][x-1]  == score - gapPenalty):
        if(traceBack1(matrix, y, x-1)):
            pathPositions.append([y,x])
            return True
    else:
        return False
 



def traceBack(matrix, y, x):
    score = matrix[y][x]
    if(x==1 or y==1):
        pathPositions.append([y,x])
        return True
    if(matrix[y-1][x-1]  == score - compare(Seq1[x-1],Seq2[y-1])):
        if(traceBack(matrix, y-1, x-1)):
            pathPositions.append([y,x])
            return True
    if(matrix[y-1][x]  == score - gapExtensionPenalty):
        if(traceBack(matrix, y-1, x)):
            pathPositions.append([y,x])
            return True
    if(matrix[y][x-1]  == score - gapExtensionPenalty):
        if(traceBack(matrix, y, x-1)):
            pathPositions.append([y,x])
            return True
    if(matrix[y-1][x]  == score - gapPenalty):
        if(traceBack(matrix, y-1, x)):
            pathPositions.append([y,x])
            return True
    if(matrix[y][x-1]  == score - gapPenalty):
        if(traceBack(matrix, y, x-1)):
            pathPositions.append([y,x])
            return True
    else:
        return False
 

def getAlignmentScore(seqA, seqB):
    score = 0
    pre = ""
    gapExtensionCount=0

    #we get score from index 0 because it has no preIndex
    score+=compare(seqA[0][0], seqB[0][0])
    for i in range(1,len(seqA[0])):
        #checking for consecutive gaps
        if(seqA[0][i]=="-" and  seqA[0][i-1]=="-"):
            gapExtensionCount+=1
        elif(seqB[i][0]=="-" and seqB[i-1][0]=="-"):
            gapExtensionCount+=1
        score+= compare(seqA[0][i], seqB[0][i])
    
    score += (gapExtensionPenalty - gapPenalty) * gapExtensionCount
    return score



def main():
    f= open('damn.txt', 'w') 
    sys.stdout = f
    if(scope == "local"):
        matrix1 = initMatrix0(Seq1,Seq2)
        maxI, maxJ=getLocalOptimalPath(matrix1, Seq1, Seq2)
        traceBack1(matrix1,maxJ, maxI)
        print ('Smith-Waterman Algorithm (local optimal alignment) result:\n')
        identityPercentage = printAlignedSequences()
        print("\nAlignment score: "+ str(matrix1[maxJ][maxI]))
        print("Identity percentage: " + "{:.2f}".format(identityPercentage)+"%")


    else:
        matrix2 = initMatrix(Seq1, Seq2)
        getOptimalGlobalPath(matrix2, Seq1, Seq2)
        traceBackStart(matrix2)
        print("Needleman–Wunsch Algorithm (global optimal alignment) result:\n")
        
        identityPercentage = printAlignedSequences()
        print("\nAlignment score: "+ str(matrix2[len(matrix2)-1][len(matrix2[0])-1]))
        print("Identity percentage: " + "{:.2f}".format(identityPercentage)+"%")
        
         

if __name__ == '__main__' :
    opts , args = getopt.getopt(sys.argv[1:],'h',['gap=','file1=','file2=', 'scope=', 'scoringfile='])
    InputFile1 = 'C:\\Users\\kykse\\Desktop\\mygod\\Sequence\\Seq1.fasta'
    InputFile2 = 'C:\\Users\\kykse\\Desktop\\mygod\\Sequence\\Seq2.fasta'
    scoringFile = 'C:\\Users\\kykse\\Desktop\\mygod\\scoring1.txt'
    gapPenalty = -5
    gapExtensionPenalty = -2
    scope="global"
    
    for op , value in opts:
        if op == '--file1':
            InputFile1 = value
        elif op == '--file2':
            InputFile2 = value
        elif op == '--scoringfile':
            InputFile2 = value
        elif op == '--gap':
            match_score = float(value)
        elif op == '--scope':
            scope = value
    Seq1 = get_string(InputFile1)
    Seq2 = get_string (InputFile2)
    main()

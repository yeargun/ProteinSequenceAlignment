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

    
def getOptimalGlobalPath(matrix, seq1, seq2):
    for j in range(1,len(matrix)):
        for i in range(1,len(matrix[0])):
            a = matrix[j-1][i] + gapPenalty
            b = matrix[j][i-1] + gapPenalty
            c = matrix[j-1][i-1] + compare(seq1[i-1],seq2[j-1])
            matrix[j][i] = max(a,b,c)
    return matrix

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

    
def loss_function(matrix,i,j,Label_Dir):
    
    a = matrix[i-1,j] + compare(Seq1[i-1],'-')
    b = matrix[i,j-1] + compare(Seq2[j-1],'-')
    c = matrix[i-1,j-1] + compare(Seq1[i-1],Seq2[j-1])
    save_CellLabel(i,j,a,b,c,Label_Dir)
    if max(a,b,c) < 0:
        matrix[i,j] = 0
    else:
        matrix[i,j] = max(a,b,c)
    return matrix
     
def save_CellLabel(m,n,a,b,c,matrix_to_return):
    
    position = '%s%s' % (str(m),str(n))
    matrix_to_return[position] = []
    if max(a,b,c) < 0:
        matrix_to_return[position].append('N')
    else:
        if max(a,b,c) == a:
            matrix_to_return[position].append('U')
        if max(a,b,c) == b:
            matrix_to_return[position].append('L')
        if max(a,b,c) == c:
            matrix_to_return[position].append('O')            

def generate_matrix(row,col):
    
    x = len(row)+1
    y = len(col)+1
    distance_matrix = mat(zeros((x,y)))
    Cell_Label = {}
    
    for i in range(1,x):
        for j in range(1,y):
            distance_matrix = loss_function(distance_matrix,i,j,Cell_Label)
    return distance_matrix , Cell_Label
    

def get_string(file):
    
    with open(file,'r') as f:
        for line in f:
            l = line.replace('\n','')
            if '>' in l:
                pass
            else:
                string = l
                
    return string 

def print_matrix(result):
    row , col = shape(result)

    print("\t\t", end = '')
    for i in Seq2:
        print(i, end='\t')
    print()
    for i in range(0, row):
        if (i > 0):
            print(Seq1[i-1], end = '\t')
        else:
            print("\t", end = '')
    
        for j in range(0, col):
            print('%d' % result[i,j], end = '\t')
        print()
        
def getLocalOptimalPath(Scoring_matrix , Label_dir):
    #get OptimalPath base Label Matrix
    x ,y  = shape(Scoring_matrix)
    Start = {'position':[],'score' :0}
    
    # Find the starting position
    for i in range(1,x):
        for j in range(1,y):
            position = '%s%s' % (str(i),str(j))
            if Scoring_matrix[i,j] == Start['score']:
                Start['position'].append(position)
            if Scoring_matrix[i,j] > Start['score']:
                Start['position'] = [position] 
                Start['score'] = Scoring_matrix[i,j]
    
    AtoReturn = []
    BtoReturn = []
    for p in range(0,len(Start['position'])):
        location = Start['position'][p]
        row = int(location[0])
        col = int(location[1])
        new_SeqA = ''
        new_Seqb = ''
        while 1:
            if 'N' in Label_dir[location]:  
                AtoReturn.append(new_SeqA[::-1])
                BtoReturn.append(new_Seqb[::-1])
                break            
            if 'L' in Label_dir[location]:
                if len(Label_dir[location]) == 1:
                    col = col-1                
                new_SeqA = new_SeqA + '-'
                new_Seqb = new_Seqb + Seq2[col]                                
            elif 'U' in Label_dir[location]:
                if len(Label_dir[location]) == 1:
                    row = row-1                 
                new_Seqb = new_Seqb + '-'
                new_SeqA = new_SeqA + Seq1[row]                             
            elif 'O' in Label_dir[location]:
                row , col = row-1 ,col-1
                new_SeqA = new_SeqA + Seq1[row]
                new_Seqb = new_Seqb + Seq2[col]      
            if len(Label_dir[location]) > 1:
                row , col = row-1 , col-1
            location = '%s%s' % (str(row),str(col))

    return AtoReturn , BtoReturn
            
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
    return matchCount*100 / len(alignedSeq1)
    
    
def traceBackStart(matrix):
    xMax = len(matrix[0])        
    yMax = len(matrix)
    score = matrix[yMax-1][xMax-1] 

    if(matrix[yMax-2][xMax-2]  == score - compare(Seq1[xMax-2],Seq2[yMax-2])):
        if(traceBack(matrix, yMax-2, xMax-2)):
            pathPositions.append([yMax-1,xMax-1])
            return True
    elif(matrix[yMax-2][xMax-1]  == score - gapPenalty):
        if(traceBack(matrix, yMax-2, xMax-1)):
            pathPositions.append([yMax-1,xMax-1])
            return True
    elif(matrix[yMax-1][xMax-2]  == score - gapPenalty):
        if(traceBack(matrix, yMax-1, xMax-2)):

            pathPositions.append([yMax-1,xMax-1])
            return True
    else:
        print("something is wrong here")            
    
def traceBack(matrix, y, x):
    score = matrix[y][x]
    if(x==1 or y==1):
        pathPositions.append([x,y])
        return True
    if(matrix[y-1][x-1]  == score - compare(Seq1[x-1],Seq2[y-1])):
        if(traceBack(matrix, y-1, x-1)):
            pathPositions.append([y,x])
            return True
    elif(matrix[y-1][x]  == score - gapPenalty):
        if(traceBack(matrix, y-1, x)):
            pathPositions.append([y,x])
            return True
    elif(matrix[y][x-1]  == score - gapPenalty):
        if(traceBack(matrix, y, x-1)):
            pathPositions.append([y,x])
            return True
    else:
        return False
 

def getAlignmentScore(seqA, seqB):
    score = 0
    for i in range(len(seqA[0])):
        score+= compare(seqA[0][i], seqB[0][i])
    return score



def main():
    if(scope == "local"):
        Distance_matrix , Cell_Label = generate_matrix(Seq1,Seq2)
        print_matrix(Distance_matrix)
        SequenceA,SequenceB = getLocalOptimalPath(Distance_matrix,Cell_Label)
        print_result(SequenceA,SequenceB)
        print("Alignment score: "+ str( getAlignmentScore(SequenceA,SequenceB )))
        printIdentityPercentage(SequenceA,SequenceB)
        
    else:
        matrix2 = initMatrix(Seq1, Seq2)
        getOptimalGlobalPath(matrix2, Seq1, Seq2)
        printMatrix(matrix2)
        traceBackStart(matrix2)
        print("Needleman–Wunsch Algorithm (global optimal alignment) result:")
        identityPercentage = printAlignedSequences()
        print("Alignment score: "+ str(matrix2[len(matrix2)-1][len(matrix2[0])-1]))
        print("Identity percentage: " + "{:.2f}".format(identityPercentage)+"%")
        
        
    
    
    

if __name__ == '__main__' :
    opts , args = getopt.getopt(sys.argv[1:],'h',['gap=','file1=','file2=', 'scope=', 'scoringfile='])
    InputFile1 = 'C:\\Users\\kykse\\Desktop\\mygod\\Sequence\\Seq1.fasta'
    InputFile2 = 'C:\\Users\\kykse\\Desktop\\mygod\\Sequence\\Seq2.fasta'
    scoringFile = 'C:\\Users\\kykse\\Desktop\\mygod\\scoring1.txt'
    gapPenalty = -2
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

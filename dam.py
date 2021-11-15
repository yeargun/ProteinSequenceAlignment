
# file Smith-Waterman
# /usr/bin/python
# coding : utf-8

from os import X_OK
import numpy as np
from numpy import *
import sys 
import getopt

substitionMatrixFile = open("scoring1.txt")
lines = substitionMatrixFile.readlines()
aminoacids= lines[20].split("\t")


def initMatrix(seq1, seq2):
    print("initmatrix x " + str(len(seq1)+1))
    print("initmatrix y " + str(len(seq2)+1))
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
    print(distanceMatrix)
    return distanceMatrix

def loss_function2(matrix,i,j,Label_Dir):
    
    a = matrix[i-1,j] + compare(Seq1[i-1],'-')
    b = matrix[i,j-1] + compare(Seq2[j-1],'-')
    c = matrix[i-1,j-1] + compare(Seq1[i-1],Seq2[j-1])
    save_CellLabel(i,j,a,b,c,Label_Dir)
    if max(a,b,c) < 0:
        matrix[i,j] = 0
    else:
        matrix[i,j] = max(a,b,c)
    return matrix

def printMatrix(matrix):
    for y in range(len(matrix)):
        for x in range(len(matrix[0])):
            print(str(matrix[y][x]), end= "\t")
        print("")



def generate_matrix(row,col):
    
    x = len(row)+1
    y = len(col)+1
    distance_matrix = mat(zeros((x,y)))
    Cell_Label = {}
    
    for i in range(1,x):
        for j in range(1,y):
            distance_matrix = loss_function(distance_matrix,i,j,Cell_Label)
            
    return distance_matrix , Cell_Label
    
def getOptimalGlobalPath(matrix, seq1, seq2):
    
    print("seq1 " + seq1)
    print("seq2 " + seq2)
    for j in range(1,len(matrix)):
        for i in range(1,len(matrix[0])):
            printMatrix(matrix)
            print("seq1 " + seq1)
            print("seq2 " + seq2)
            print("i  "+str(i)+ "   j "+str(j))
            a = matrix[j-1][i] + gapPenalty
            b = matrix[j][i-1] + gapPenalty
            c = matrix[j-1][i-1] + compare(seq1[i-1],seq2[j-1])
            matrix[j][i] = max(a,b,c)
            print(max(a,b,c))
    print(matrix)
    print(gapPenalty)
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
    
    print ('socore matrix: ')
    row , col = shape(result)
    print()
    print("        ", end = '')
    for i in Seq2:
        print(i, end='\t')
    print()
    for i in range(0, row):
        if (i > 0):
            print(Seq1[i-1], end = '\t')
        else:
            print("  ", end = '')
    
        for j in range(0, col):
            print('%d' % result[i,j], end = '\t')
        print()
        
def get_OptimalPath(Scoring_matrix , Label_dir):
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
    print ('\n')
    print ('    Result    \n')
    print ('Sequences : ')
    print ('Sequence 1 : ' + Seq1)
    print ('Sequence 2 : ' + Seq2+'\n')
    print ('Paramenters :')
    print ('Substitution matrix :  a = b   S(a,b)= ' )
    print ('                       a != b  S(a,b)= ' +'\n')
    print ('Result: ')
    for a in range(len(sequenceA)):
        for b in range(len(sequenceB)):           
            print ('Sequence1   ' + str(sequenceA[a]))
            print ('Sequence2   ' + str(sequenceB[b]))
    
    
    
            
            
    
    


def main():
    
    Distance_matrix , Cell_Label = generate_matrix(Seq1,Seq2)
    print_matrix(Distance_matrix)
    SequenceA,SequenceB = get_OptimalPath(Distance_matrix,Cell_Label)
    print_result(SequenceA,SequenceB)
    print("forget about it")

    #file1 = open("scoring1.txt")
    #lines = file1.readlines()
    #print(compare("R","R"))

    matrix2 = initMatrix(Seq1, Seq2)
    getOptimalGlobalPath(matrix2, Seq1, Seq2)


    
    
    
    

if __name__ == '__main__' :
    opts , args = getopt.getopt(sys.argv[1:],'h',['gap=','file1=','file2='])
    InputFile1 = 'C:\\Users\\kykse\\Desktop\\mygod\\Sequence\\Seq1.fasta'
    InputFile2 = 'C:\\Users\\kykse\\Desktop\\mygod\\Sequence\\Seq2.fasta'
    gapPenalty = -2
    
    for op , value in opts:
        if op == '--file1':
            InputFile1 = value
        elif op == '--file2':
            InputFile2 = value
        elif op == '--gap':
            match_score = float(value)
            
    Seq1 = get_string(InputFile1)
    Seq2 = get_string (InputFile2)
    print ('matrix shape: '+str(len(Seq1))+'*'+str(len(Seq2)))
    main()

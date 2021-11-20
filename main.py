from os import X_OK
import sys 
import getopt


pathPositions = []



def get_string(file):
    
    with open(file,'r') as f:
        for line in f:
            l = line.replace('\n','')
            if '>' in l:
                pass
            else:
                string = l
                
    return string



def compare(str1,str2):
    if(str2 == "-"):
        return gapPenalty
    index1 = aminoacids.index(str1)
    index2 = aminoacids.index(str2)
    return int(lines[index1-1][index2])


def printAlignedSequences():
    alignedSeq1 = ""
    alignedSeq2 = ""
    preI = 0
    preJ = 0
    for i,j in pathPositions:
        if(preI == 0 and preJ ==0):
            alignedSeq1 += Seq1[j-1]
            alignedSeq2 += Seq2[i-1]
        elif(preI+1 == i and preJ+1 == j):
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
        if(alignedSeq1[i]==alignedSeq2[i]):
            matchCount+=1
    if(len(alignedSeq1)==0):
        return 1
    return matchCount*100 / len(alignedSeq1)
    
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
        
            if(gapPositions[j-1][i]!= 1):
                a = matrix[j-1][i] + gapPenalty
            if(gapPositions[j][i-1]!= 1):
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
        
            if(gapPositions[j-1][i]!= 1):
                a = matrix[j-1][i] + gapPenalty
            if(gapPositions[j][i-1]!= 1):
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
 
def initScoringMatrix(scoringFile):
    substitionMatrixFile = open(scoringFile)
    rawLines = substitionMatrixFile.readlines()
    aminoacids = rawLines[6].split("  ")
    aminoacids[1] = aminoacids[1][1]
    lines = []
    tempLine= []
    for i in range(7,len(rawLines)):
        tempLine = rawLines[i].replace("  ", " ")
        tempLine = tempLine.split(" ")
        lines.append(tempLine)
        
    return aminoacids, lines



def main():
    print("-"*80)
    print("Scoring matrix is adjustable but the file format must be as same as given BLOSUM62.txt")
    print("How to run the program ")
    print("python main.py --scope=global  --prots=prots.txt  --scoringfile=blosum62.txt ")
    print("also --gap, --gapextension can be used")
    print("output is at outtt1.txt ")
    print("-"*80)

    f= open('outtt1.txt', 'w') 
    sys.stdout = f
    print(scoringFile + "\tgap opening penalty: " + str(gapPenalty)+ "\tgap extension penalty: "+ str(gapExtensionPenalty))
    global aminoacids, lines
    aminoacids, lines = initScoringMatrix(scoringFile)
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
    opts , args = getopt.getopt(sys.argv[1:],'h',['gap=','prots=', 'scope=', 'scoringfile=', 'gapextension='])
    proteinsFile = "prots.txt"
    scoringFile = 'scoring1.txt'
    gapPenalty = -10
    gapExtensionPenalty = -5
    scope="global"

    
    for op , value in opts:
        if op == '--prots':
            proteinsFile = value
        elif op == '--scoringfile':
            scoringFile = value
        elif op == '--gap':
            gapPenalty = float(value)
        elif op == '--scope':
            scope = value
        elif op == '--gapextension':
            gapExtensionPenalty = float(value)
    
    proteinsFile1 = open(proteinsFile, "r")
    proteins = proteinsFile1.readlines()
    Seq1 = proteins[0].split("\n")[0]
    Seq2 = proteins[1].split("\n")[0]
    main()

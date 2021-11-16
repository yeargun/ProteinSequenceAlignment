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



def initMatrix0(Seq1,Seq2):
    matrix = []
    for j in range(len(Seq2)+1):
        row = []
        for  i in range(len(Seq1)+1):
            row.append(0)
        matrix.append(row)
    return matrix



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












def printMatrix(matrix):
    for y in range(len(matrix)):
        for x in range(len(matrix[0])):
            print(str(matrix[y][x]), end= "\t")
        print("")


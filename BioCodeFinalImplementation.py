# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 19:07:33 2020

@author: HP
"""

# -- coding: utf-8 --
"""
Created on Mon Apr 20 16:10:30 2020

@author: HP
"""
import csv
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasClassifier
from keras.utils import np_utils
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
def preprocess(prot,filename,linepop):
    f=open(filesPathIp+"\\"+filename+".dssp", "r")
    contents =f.readlines()
    b=[]
    c=[]
    flag=False
    for eachline in contents:
        b=[]
        if '#' in eachline.split(" "):
            flag=True
            continue
        if flag:
            a=eachline.split(" ")
            for i in range(0,len(a)):
                if i==12 and a[i]=='':
                    a[i]='Z'
                    b.append(a[i])
                elif a[i]!='':
                    b.append(a[i])
            c.append(b)
    for i in range(0,len(c)):
        for j in range(0,len(c[i])):
            if j==1:
                c[i][j]=c[i][j]+c[i][j+1]
                c[i].pop(j+1)
    for i in range(0,len(c)):
        for j in range(0,17):
            if c[i][j].endswith(','):
                c[i][j]=c[i][j]+c[i][j+1]
                c[i].pop(j+1)
        if len(c[i])==32:
            c[i].pop(c[i].index('Z'))
        elif len(c[i])==31:
            if 'Z' in c[i]:
                c[i].pop(c[i].index('Z'))
    
    for i in range(0,len(c)):
        for j in range(3,len(c[i])):
            if c[i][j]=='0':
                k=j
                break
            else:
                c[i][3]+=c[i][j]
        del c[i][4:k]
    for l in range(0,2):
        for i in range(0,len(c)):
            for j in range(11,len(c[i])):
                if c[i][j].rfind('-')!=0:
                    strn=c[i][j]
                    if len(c[i][j].split('-'))>=2 or (c[i][j].find('-')==0 and c[i][j].rfind('-')!=0):
                        ind=strn.rfind('-')
                        c[i].append(c[i][len(c[i])-1])
                        for k in range(len(c[i])-2,j+1,-1):
                            c[i][k]=c[i][k-1]
                        c[i][j]=strn[:ind]
                        c[i][j+1]=strn[ind:]
    for i in linepop:
        c.pop(i-1)
    for i in c:
        i.append(prot+1)
    #    lens=[]
    #    for i in c:
    #        lens.append(len(i))       
    #    #print(list(set(lens)))
    #    #print(c[len(c)-1])
    return c
def writeData():
    proteins=['1a01','1aj9','1bz1','1bzz','1dxt','1dxu','1dxv','1g9v','1gli','1hba','1hbb','1ljw','1mko','1o1j','1o1l','1o1p','1qsh','1y0t','1y0w','1yzi','1i3d','1i3e']
    emptyLines=[[431,289,142],[142],[433,290,143],[431,289,142],[432,290,142],[431,289,142],[431,289,142],[431,289,142],[431,289,142],[431,289,142],[431,289,142],[142],[431,289,142],[431,284],[431,284],[431,284],[431,289,142],[431,289,142],[431,289,142],[142],[147],[147]]
    
    finalDataPush=[['#', 'RESIDUE', 'AA', 'STRUCTURE', 'BP1', 'BP2', 'ACC', 'N-H-->O1', 'O-->H-N1', 'N-H-->O2', 'O-->H-N2', 'TCO', 'KAPPA', 'ALPHA', 'PHI', 'PSI', 'X-CA', 'Y-CA', 'Z-CA', 'CHAIN', 'AUTHCHAIN', 'NUMBER', 'RESNUM', 'BP1', 'BP2', 'N-H-->O', 'O-->H-N', 'N-H-->O', 'O-->H-N','ProteinID']]
    for i in range(0,len(proteins)):
        finalDataPush+=preprocess(i,proteins[i],emptyLines[i])
    
    with open(filesPathIp+"\\FinalDataset.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(finalDataPush)
def completeDataPreprocess():
    writeData()
    labels=[]
    readData=pd.read_csv(filesPathIp+"\\FinalDataset.csv")
    p=0
    #print(list(set(readData['STRUCTURE'])))
    for eachRow in readData['STRUCTURE']:
        if eachRow[0]=='H' or eachRow[0]=='G':
            labels.append('H')
        elif eachRow[0]=='E' or eachRow[0]=='B':
            labels.append('E')
        else:
            labels.append('C')
        p+=1
    #print(len(labels))
    #print(len(readData['STRUCTURE']))
    readData['label']=labels
    #print(readData)
    readData.drop(["AUTHCHAIN","#","NUMBER","BP1","BP2","STRUCTURE","CHAIN","RESNUM",'N-H-->O1', 'O-->H-N1', 'N-H-->O2', 'O-->H-N2'],axis=1,inplace=True)
    #print(readData)
    
    readData.to_csv(filesPathIp+"\\FinalDataset2.csv", index = False)#, header=False)
    readData2=pd.read_csv(filesPathIp+"\\FinalDataset2.csv")
    #print(readData2)
    encoder = LabelEncoder()
    colnames=list(readData2.columns)
    for eachCol in colnames:
        temp=readData2[eachCol].values
        encoder.fit(temp)
        encodedtemp=encoder.transform(temp)
        readData2[eachCol]=encodedtemp
    readData2.to_csv(filesPathIp+"\\FinalDataset4.csv", index = False)
    readData2.to_csv(filesPathIp+"\\FinalDataset3.csv", index = False, header=False)
    #print(readData.shape[0])
    #print(readData.shape[1])
    #encoder = LabelEncoder()
    #temp=readData2['RESIDUE'].values
    #encoder.fit(temp)
    #encodedtemp = encoder.transform(temp)
    #readData2['RESIDUE']=encodedtemp
    #print(readData2)
    #Y=readData2.values[0:,24]
    
    #encoder.fit(Y)
    #encoded_Y = encoder.transform(Y)
def oneLayer():
    DNN = Sequential()
    DNN.add(Dense(100, input_dim=18, activation='relu'))
    #DNN.add(Dropout(0.4))
    #print("1")
    textList.append("1")
    DNN.add(Dense(len(list(set(dataset[:,18]))), activation='softmax'))
    #print("1")
    textList.append("1")
    DNN.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    #print("------")
    return DNN
def twoLayers():
    DNN = Sequential()
    DNN.add(Dense(100, input_dim=18, activation='relu'))
    #DNN.add(Dropout(0.4))
    #print("1")
    textList.append("1")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    DNN.add(Dense(len(list(set(dataset[:,18]))), activation='softmax'))
    #print("2")
    textList.append("2")
    DNN.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    #print("------")
    return DNN
def threeLayers():
    DNN=Sequential()
    DNN.add(Dense(100, input_dim=18, activation='relu'))
    #DNN.add(Dropout(0.4))
    #print("1")
    textList.append("1")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    DNN.add(Dense(len(list(set(dataset[:,18]))), activation='softmax'))
    #print("3")
    textList.append("3")
    DNN.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    #print("------")
    return DNN
def fourLayers():
    DNN = Sequential()
    DNN.add(Dense(100, input_dim=18, activation='relu'))
    #DNN.add(Dropout(0.4))
    #print("1")
    textList.append("1")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("2")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("3")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("4")
    DNN.add(Dense(len(list(set(dataset[:,18]))), activation='softmax'))
    #print("4")
    textList.append("4")
    DNN.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    #print("------")
    return DNN
def fiveLayers():
    DNN = Sequential()
    DNN.add(Dense(100, input_dim=18, activation='relu'))
    #DNN.add(Dropout(0.4))
    #print("1")
    textList.append("1")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("2")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("3")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("4")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    DNN.add(Dense(len(list(set(dataset[:,18]))), activation='softmax'))
    #print("5")
    textList.append("5")
    DNN.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    #print("------")
    return DNN
def sixLayers():
    DNN = Sequential()
    DNN.add(Dense(100, input_dim=18, activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("1")
    textList.append("1")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("2")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("3")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("4")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("5")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    DNN.add(Dense(len(list(set(dataset[:,18]))), activation='softmax'))
    textList.append("6")
    DNN.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    ##print("------")
    return DNN
def sevenLayers():
    DNN = Sequential()
    DNN.add(Dense(100, input_dim=18, activation='relu'))
    #DNN.add(Dropout(0.4))
    #print("1")
    textList.append("1")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("2")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("3")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("4")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("5")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("6")
    DNN.add(Dense(100,activation='relu'))
    #DNN.add(Dropout(0.4))
    ##print("6")
    DNN.add(Dense(len(list(set(dataset[:,18]))), activation='softmax'))
    textList.append("7")
    DNN.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])
    #print("------")
    textList.append("------")
    return DNN
def iterate(l,e):
    #for i in range(0,5):
    if l==1:
        estimator = KerasClassifier(build_fn=oneLayer, epochs=e, batch_size=5, verbose=0)
        #print("////////////////////////")
        kfold = KFold(n_splits=10, shuffle=True)
        #print("----------------------")
        results = cross_val_score(estimator, X, testY, cv=kfold)
    elif l==2:
        estimator = KerasClassifier(build_fn=twoLayers, epochs=e, batch_size=5, verbose=0)
        #print("////////////////////////")
        kfold = KFold(n_splits=10, shuffle=True)
        #print("----------------------")
        results = cross_val_score(estimator, X, testY, cv=kfold)
    elif l==3:
        estimator = KerasClassifier(build_fn=threeLayers, epochs=e, batch_size=5, verbose=0)
        #print("////////////////////////")
        kfold = KFold(n_splits=10, shuffle=True)
        #print("----------------------")
        results = cross_val_score(estimator, X, testY, cv=kfold)
    elif l==4:
        estimator = KerasClassifier(build_fn=fourLayers, epochs=e, batch_size=5, verbose=0)
        #print("////////////////////////")
        kfold = KFold(n_splits=10, shuffle=True)
        #print("----------------------")
        results = cross_val_score(estimator, X, testY, cv=kfold)
    elif l==5:
        estimator = KerasClassifier(build_fn=fiveLayers, epochs=e, batch_size=5, verbose=0)
        #print("////////////////////////")
        kfold = KFold(n_splits=10, shuffle=True)
        #print("----------------------")
        results = cross_val_score(estimator, X, testY, cv=kfold)
    elif l==6:
        estimator = KerasClassifier(build_fn=sixLayers, epochs=e, batch_size=5, verbose=0)
        #print("////////////////////////")
        kfold = KFold(n_splits=10, shuffle=True)
        #print("----------------------")
        results = cross_val_score(estimator, X, testY, cv=kfold)
    elif l==7:
        estimator = KerasClassifier(build_fn=sevenLayers, epochs=e, batch_size=5, verbose=0)
        #print("////////////////////////")
        kfold = KFold(n_splits=10, shuffle=True)
        #print("----------------------")
        results = cross_val_score(estimator, X, testY, cv=kfold)
    #print("---------------------------------------$$$$$$------------------")
    textList.append("---------------------------------------$$$$$$------------------")
    ##print(k," iterations")
    #textList.append("==================================================="+str(k)+" iterations===================================================")
    print("Baseline accuracy for",l,"layers",e,"epochs:",results.mean()*100,results.std()*100)
    textList.append("Baseline accuracy for "+str(l)+" layers "+str(e)+" epochs:"+str((results.mean()*100))+" "+str((results.std()*100)))
    ##print("%.2f%% (%.2f%%)" % (results.mean()*100, results.std()*100))

#F:\\BioCodeAndData\\BioFilesData\\pdb\\dssp
#D:\\BioFilesData\\pdb\\dssp
textList=[]
filesPathIp=input("Input the path to the folder where all the dssp files are located ")
#filesPathIp="D:\\BioFilesData\\pdb\\dssp"
completeDataPreprocess()
loadData=pd.read_csv(filesPathIp+"\\FinalDataset3.csv", header=None)
#print(loadData)
dataset=loadData.values
X=dataset[:,0:18].astype(float)
Y=dataset[:,18]
labelEncoder = LabelEncoder()
labelEncoder.fit(Y)
vectorizedValY=labelEncoder.transform(Y)
testY=np_utils.to_categorical(vectorizedValY)
acc=[]
layers=[1,2,3,4,5,6]
epochs=[10,50,80,200]
iterations=[5,6]

for i in layers:
    for j in epochs:
        for k in iterations:
            #print(k," iterations")
            textList.append(str(k)+" iterations")
            iVal=k
            while(iVal!=0):
                #print(i,"layers",j,"epochs")
                textList.append(str(i)+" layers and "+str(j)+" epochs")
                iterate(i,j)
                iVal-=1
outputtxtfile=open(filesPathIp+"\\outputfilewrite.txt", "w+")
outputtxtfile.writelines(textList)
outputtxtfile.close()
#print("All the outputs are saved in this path: "+filesPathIp+"\\outputfilewrite.txt")
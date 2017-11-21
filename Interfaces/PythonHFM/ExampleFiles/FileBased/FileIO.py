import numbers;
import numpy as np;
import os;
from operator import mul
from functools import reduce
import subprocess

#These two methods export a dictonary to a (pair of) files, and conversely.
def RulesToRaw(params,prefix='input'):
    if(not isinstance(params,dict) or not isinstance(prefix,str)):
        print("Invalid parameters")
        return
    f = open(prefix+'_Format.txt','w')
    data=[]
    for key,val in params.items():
        if isinstance(val,numbers.Number):
            f.write(key+'\n0\n\n')
            data.append([val])
        elif isinstance(val,str):
            f.write(key+'\n-1\n'+val+'\n\n')
        elif isinstance(val,np.ndarray):
            f.write(key+'\n'+str(val.ndim)+'\n')
            for dim in val.shape:
                f.write(str(dim)+'\n')
            f.write('\n')
            data.append(val.flatten())
        else:
            raise ValueError('Invalid type for key ' + key);
    f.close()
    np.concatenate(data).astype('d').tofile(prefix+'_Data.dat')
    
def RawToRules(prefix='output'):
    data=np.fromfile(prefix+'_Data.dat')
    pos=0;
    f=open(prefix+'_Format.txt')
    dict={}
    while True:
        key=f.readline().strip()
        if not key: break
        keyType = int(f.readline())
        if keyType==-1:
            dict[key]=f.readline().strip()
        elif keyType==0:
            dict[key]=data[pos]
            pos+=1
        else:
            dims=[int(f.readline()) for i in range(keyType)]
            size=reduce(mul,dims)
            dict[key]=np.reshape(data.take(np.arange(pos,pos+size)),dims)
            pos+=size
        f.readline()
    return dict

def WriteCallRead(inputData,executable,binary_dir=None,
                  inputPrefix="input",outputPrefix="output",logFile="log.txt"):
    if binary_dir: #Change to executable's directory
        cwd=os.getcwd()
        os.chdir(binary_dir)
        
    RulesToRaw(inputData,inputPrefix) # Export inputData
    
    execPrefix = '' if os.name=='nt' else './' # Run executable
    command = execPrefix+executable+' '+inputPrefix+' '+outputPrefix
    if logFile: command = command +' > '+logFile
    retcode = subprocess.call(command, shell=True)
    if retcode!=0:  print('Returned with exit code ', retcode)
    
    outputData = RawToRules(outputPrefix) # Import outputData
    if logFile: outputData['log'] = open(logFile).read()
    outputData['retcode']=retcode;
    if binary_dir:  os.chdir(cwd)
        
    return outputData

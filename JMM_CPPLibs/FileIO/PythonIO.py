import numbers
import numpy as np
import importlib.util

def WriteCallRead_Lib(params,factory,binary_dir=None,extension='.so'):
    # Build object
    spec = importlib.util.spec_from_file_location(factory, binary_dir + factory + extension)
    factory = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(factory)
    obj = factory.Build() 
    
    # Export data
    for key,val in params.items():
        if isinstance(val,numbers.Number):
            obj.SetScalar(key,val)
        elif isinstance(val,str):
            obj.SetString(key,val)
        elif isinstance(val,np.ndarray):
            obj.SetArray(key,val)
        else:
            raise ValueError('Invalid type for key '+key)
            
    # Process data
    obj.Run()
    
    # Import data back
    result = dict()
    computed = obj.GetComputedKeys()
    computed = computed[2:-3]
    
    for keyType in computed.split("),("):
        key,t = keyType.split(",")
        if t=='float':
            result[key]=obj.GetScalar(key)
        elif t=='array':
            result[key]=obj.GetArray(key)
        elif t=='string':
            result[key]=obj.GetString(key)
            
    return result
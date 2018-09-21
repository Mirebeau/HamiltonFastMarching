
import numpy as np
import numbers
import importlib
import ast

# Choose between Python library or command line executables
FileHFM_binary_dir=None
FileIO_dir=None

# ----- Basic utilities functions -----

def GetGeodesics(output,suffix=''): 
	if suffix != '': suffix='_'+suffix
	return np.vsplit(output['geodesicPoints'+suffix],
					 output['geodesicLengths'+suffix].cumsum()[:-1].astype(int))

SEModels = {'ReedsShepp2','ReedsSheppForward2','Elastica2','Dubins2',
'ReedsSheppExt2','ReedsSheppForwardExt2','ElasticaExt2','DubinsExt2',
'ReedsShepp3','ReedsSheppForward3'}

def GetCorners(params):
	dims = params['dims']
	dim = len(dims)
	h = params['gridScales'] if 'gridScales' in params.keys() else [params['gridScale']]*dim
	origin = params['origin'] if 'origin' in params.keys() else [0.]*dim
	if params['model'] in SEModels:
		origin = np.append(origin,[0]*(dim-len(origin)))		
		hTheta = 2*np.pi/dims[-1]
		h[-1]=hTheta; origin[-1]=-hTheta/2;
		if dim==5: h[-2]=hTheta; origin[-2]=-hTheta/2;
	return [origin,origin+h*dims]

def centeredLinspace(a,b,n): 
	r,dr=np.linspace(a,b,n,endpoint=False,retstep=True)
	return r+dr/2

def GetAxes(params,dims=None):
	bottom,top = GetCorners(params)
	if dims is None: dims=params['dims']
	return [centeredLinspace(b,t,d) for b,t,d in zip(bottom,top,dims)]

def GetGrid(params,dims=None):
	axes = GetAxes(params,dims);
	ordering = params['arrayOrdering']
	if ordering=='RowMajor':
		return np.mgrid(*axes)
	elif ordering=='YXZ_RowMajor':
		return np.meshgrid(*axes)
	else: 
		raise ValueError('Unsupported arrayOrdering : '+ordering)

def MakeHFMRect(corner0,corner1,sampleBoundary=False,gridScale=None,gridScales=None,dimx=None,dims=None):
	dim = len(corner0)
	sb=float(sampleBoundary)
	result=dict()
	width = np.array(corner1)-np.array(corner0)
	if gridScale is not None:
		gridScales=[gridScale]*dim; result['gridScale']=gridScale
	elif gridScales is not None:
		result['gridScales']=gridScales
	elif dimx is not None:
		gridScales=width/(dimx-sb); result['gridScale']=gridScales[0]
	elif dims is not None:
		gridScales=width/(np.array(dims)-sb); result['gridScales']=gridScales
	else: 
		raise ValueError('Missing argument gridScale, gridScales, dimx, or dims')

	h=gridScales
	ratios = [(M-m)/delta+sb for delta,m,M in zip(h,corner0,corner1)]
	dims = [round(r) for r in ratios]
	assert(np.min(dims)>0)
	origin = [c+(r-d-sb)*delta/2 for c,r,d,delta in zip(corner0,ratios,dims,h)]
	result.update({'dims':np.array(dims),'origin':np.array(origin)});
	return result
	

"""
def SetBox(params,corner0,corner1,h0,sampleBoundary=False):
	dim = len(corner0);
	h = h0
	if isinstance(h0,numbers.Number):
		h = [h]*dim
		params['gridScale'] = h0
	else:
		params['gridScales'] = np.array(h0)
	assert(len(corner1)==dim); assert(len(h)==dim)
	
	sb=float(sampleBoundary)
	ratios = [(M-m)/delta+sb for delta,m,M in zip(h,corner0,corner1)]
	dims = [round(r) for r in ratios]
	assert(np.min(dims)>0)
	origin = [c+(r-d-sb)*delta/2 for c,r,d,delta in zip(corner0,ratios,dims,h)]
	linsp = [o+delta*(0.5+np.arange(d)) for o,delta,d in zip(origin,h,dims)]
	
	params['dims'] = np.array(dims)
	params['origin'] = np.array(origin)
	
	return linsp
"""
# ------- Plotting (does it belong here ?) ------

def SetTitle3D(ax,title):
	ax.text2D(0.5,0.95,title,transform=ax.transAxes,horizontalalignment='center')

# --------- Input/Output ---------

def SetInput(hfm,params):
	for key,val in params.items():
		if isinstance(val,numbers.Number):
			hfm.set_scalar(key,val)
		elif isinstance(val,str):
			hfm.set_string(key,val)
		elif isinstance(val,np.ndarray):
			hfm.set_array(key,val)
		else:
			raise ValueError('Invalid type for key ' + key);

def GetOutput(hfm):
	comp=hfm.computed_keys()
	if comp[0]=='(': # Should be patched now
		comp=comp.replace('),)',"']]")
		comp=comp.replace('),(',')(')
		comp=comp.replace(',',"','")
		comp=comp.replace(')(',"'],['")
		comp=comp.replace('((',"[['")

	result = {}
	for key,t in ast.literal_eval(comp):
		if t=='float':
			result[key] = hfm.get_scalar(key)
		elif t=='string':
			result[key] = hfm.get_string(key)
		elif t=='array':
			result[key] = hfm.get_array(key)
		else:
			raise ValueError('Unrecognized type '+ t + ' for key '+ key)
	return result

def Run(params):
	modelName = params['model']
	if FileHFM_binary_dir is None:
		moduleName = 'HFMpy.HFM_'+modelName
		HFM = importlib.import_module(moduleName)
		hfm = HFM.HFMIO()
		SetInput(hfm,params)
		hfm.run()
		return GetOutput(hfm)
	else:
		if FileIO_dir is None:
			import FileIO
		else:
			import os
			cwd=os.getcwd()
			os.chdir(FileIO_dir)
			import FileIO
			os.chdir(cwd)
		execName = 'FileHFM_'+modelName
		return FileIO.WriteCallRead(params, execName, FileHFM_binary_dir)



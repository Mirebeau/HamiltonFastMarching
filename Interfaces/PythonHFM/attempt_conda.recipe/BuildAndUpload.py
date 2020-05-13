import os

def BuildAndExport(PYTHON_VERSION="3.7"): 
#	os.system(f"conda create -q -n test-environment-{PYTHON_VERSION} python={PYTHON_VERSION} conda-build conda-verify anaconda-client")
	os.system(f"conda activate test-environment-{PYTHON_VERSION}")
#	os.system(f"conda info -a")
#	os.system(f"conda build . --output-folder dist/conda-{PYTHON_VERSION} --python={PYTHON_VERSION} -c agd-lbr")
#	OUTPUT_FILE = `find dist/conda-$PYTHON_VERSION -name "*.tar.bz2"`
#	os.system(f'anaconda upload --force --user agd-lbr "{OUTPUT_FILE}')
if __name__=='__main__':
	BuildAndExport()


"""
	- conda config --set always_yes yes --set changeps1 no 
	- conda update -q conda
	# build
	- conda create -q -n test-environment-$PYTHON_VERSION python=$PYTHON_VERSION conda-build conda-verify anaconda-client
	- source activate test-environment-$PYTHON_VERSION
	- conda info -a
	- echo Y | anaconda login --username $ANACONDA_LOGIN --password $ANACONDA_PASSWORD
	- conda build . --output-folder dist/conda-$PYTHON_VERSION --python=$PYTHON_VERSION -c agd-lbr
	- OUTPUT_FILE=`find dist/conda-$PYTHON_VERSION -name "*.tar.bz2"`
	- anaconda upload --force --user agd-lbr "$OUTPUT_FILE"
"""
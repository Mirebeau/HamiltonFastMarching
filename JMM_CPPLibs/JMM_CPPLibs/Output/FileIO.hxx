// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef FileIO_hxx
#define FileIO_hxx

FileIO::FileIO(std::string inputPrefix, std::string outputPrefix_):outputPrefix(outputPrefix_){
    
    // ------- Read format file -------
	//	clock_t top = clock();

    struct {
        void operator()(std::istream & is){
            is.ignore(std::numeric_limits<std::streamsize>::max(), is.widen('\n'));}
    } nextLine;
    
    size_t pos=0;
    std::ifstream ifs;
    ifs.open(inputPrefix+"_Format.txt");
    if(!ifs.is_open())
        ExceptionMacro("FileIO error : Input format file could not be opened "<<inputPrefix<<"_Format.txt.");
    
    while(!ifs.eof()){
        std::string key; ifs >> key; nextLine(ifs);
        if(key.empty()){
            if(ifs.eof()) break; // Allow for a final blank line
            else ExceptionMacro("FileIO Error empty field name");
        }
        if(ifs.eof() || ifs.fail()) ExceptionMacro("FileIO error : invalid input format file.");
        
        RawElement & raw = CreateElement(key);
        
        int ndims; ifs >> ndims; nextLine(ifs);
        if(ifs.eof() || ifs.fail() || ndims<-1)
            ExceptionMacro("FileIO error : invalid dimension value (" << std::to_string(ndims) << ") for key " << key << "in input format file.");
        
        if(ndims==-1){
            ifs >> raw.str; nextLine(ifs);
            if(ifs.eof() || ifs.fail())
                ExceptionMacro("FileIO error : invalid string value for key " << key << "in input format file.");
        } else {
            raw.dims.resize(ndims);
            DiscreteType size=1;
            for(int i=0; i<ndims; ++i){
                DiscreteType & dimi = raw.dims[i];
                ifs >> dimi; nextLine(ifs);
                if(ifs.eof() || ifs.fail() || dimi<0)
                    ExceptionMacro("FileIO error : invalid dimension for key " << key << " in input format file.");
                size*=dimi;
            }


            raw.data.push_back(ScalarType(pos));
            raw.data.push_back(ScalarType(pos+size));
            pos+=size;
        }
        nextLine(ifs);
    }
    ifs.close();
    
    // ------- read data file --------
    if(pos==0) return;
    std::vector<ScalarType> inputData(pos);
    
    ifs.open(inputPrefix+"_Data.dat", std::ios::binary);
	
    if(!ifs.is_open())
        ExceptionMacro("FileIO error : input data file could not be opened " << inputPrefix << "_Data.dat");
    
    ifs.seekg(0, ifs.end);
    if(size_t(ifs.tellg())!=inputData.size()*sizeof(ScalarType))
        ExceptionMacro("FileIO error : input data file has incorrect size.");
    
    ifs.seekg(0, ifs.beg);
    ifs.read((char*)&inputData[0], inputData.size()*sizeof(ScalarType));

	if (ifs.rdstate() & ifs.failbit)
		ExceptionMacro("FileIO error : could not read file");
	if (ifs.gcount() != inputData.size() * sizeof(ScalarType))
		ExceptionMacro("FileIO error : could not read all data");

    ifs.close();

    for(auto & keyElem : rawElems){
        RawElement & raw = keyElem.second;
        if(raw.data.empty()) continue;
        assert(raw.data.size()==2);
        const size_t begin = (size_t)raw.data[0], end=(size_t)raw.data[1];
        raw.data.resize(end-begin);
        std::copy(inputData.begin()+begin, inputData.begin()+end, raw.data.begin());
    }
	
	//std::cout << "Read time " << (clock()-top)/double(CLOCKS_PER_SEC) << std::endl;

}

FileIO::~FileIO(){

    this->UsageReport();
    
    // Write format file
    std::ofstream ofs;
    ofs.open(outputPrefix+"_Format.txt");
    std::vector<const std::vector<ScalarType>*> datas;
    size_t totalSize=0;
    for(const auto & keyElem : rawElems){
        const RawElement & raw = keyElem.second;
        if(raw.setter!=SetterTag::Compute) continue;
        ofs << keyElem.first << "\n";
        if(raw.IsString()) ofs << "-1\n" << raw.str << "\n\n";
        else {
            ofs << raw.dims.size() << "\n";
            for(DiscreteType dim : raw.dims) ofs << dim << "\n";
            ofs << "\n";
            datas.push_back(&raw.data);
            totalSize+=raw.data.size();
        }
    }
    ofs.close();
    
	// Gather numerical data
    std::vector<ScalarType> outputData;
    outputData.reserve(totalSize);
    for(const std::vector<ScalarType>* pVec : datas){
        outputData.insert(outputData.end(),pVec->begin(),pVec->end());}
    
    // Write data file
	std::fstream fs; // faster than ofstream
	fs.open(outputPrefix+"_Data.dat",  std::ios::out | std::ios::binary);
	fs.write((char*) &outputData[0], outputData.size()*sizeof(ScalarType));
	
}
#endif /* FileIO_hxx */

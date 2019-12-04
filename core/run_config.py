import ConfigParser
from multiprocessing.dummy import cpu_count as cpu_count

#--------------------------------------------------------------------------------------
def run(readSet,paramFile,args):

    # read parameter file
    parser = ConfigParser.SafeConfigParser()
    parser.optionxform = str
    parser.read(paramFile)
    print(readSet)
    print(paramFile)
    print(args)
    # copy all options to a config object - both general params, and params for this readSet
    cfg = lambda:0
    cfg.__dict__["readSet"] = readSet
    for (paramName, paramVal) in parser.items("general"):
        if paramName in cfg.__dict__:
            raise Exception("Config file contains duplicate specification of parameter: " + paramName)
        cfg.__dict__[paramName] = paramVal
        print(paramName, paramVal)
    
    cfg.genomeFile = args.genomeFile
    cfg.numCores = args.numCores
    cfg.samtoolsMem = args.samtoolsMem
    cfg.readFile1 = args.readFile1
    cfg.readFile2 = args.readFile2
    cfg.primerFile = args.primerFile
    
    if str(cfg.numCores) == '0':      
        # use all cores if numCores = 0
        cfg.numCores = str(cpu_count())
     
    # convert some params to boolean
    if "outputDetail" in cfg.__dict__:
        cfg.outputDetail = cfg.outputDetail.lower() == "true"
    else:
        cfg.outputDetail = False
    if "multimodal" in cfg.__dict__:
        cfg.multimodal = cfg.multimodal.lower() == "true"
    else:
        cfg.multimodal = False
    if "duplex" in cfg.__dict__:
        cfg.duplex = cfg.duplex.lower() == "true"
    else:
        cfg.duplex = False
    if "deleteLocalFiles" in cfg.__dict__:        
        cfg.deleteLocalFiles = cfg.deleteLocalFiles.lower() == "true"
    else:
        cfg.deleteLocalFiles = False
    if "instrument" not in cfg.__dict__: # IonTorrent, or 
        cfg.instrument = "N/A"           # say the user forgot to specify this for Illumina - Use MiSeq as default
 
    # return config object
    return cfg

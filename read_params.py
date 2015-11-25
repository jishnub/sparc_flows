import os

def get_directory():
    codedir=os.path.dirname(os.path.abspath(__file__))
    paramsfile = os.path.join(codedir,"params.i")
    with open(paramsfile,'r') as paramsfile:
        for line in paramsfile.readlines():
            line=line.strip()
            if line.startswith("!"): continue
            if ("parameter" in line.lower()) and ("directory" in line.lower()):
                return line.split()[-1].split("'")[1].rstrip("/")
                
def get_xlength():
    codedir=os.path.dirname(os.path.abspath(__file__))
    paramsfile = os.path.join(codedir,"params.i")
    with open(paramsfile,'r') as paramsfile:
        for line in paramsfile.readlines():
            line=line.strip()
            if line.startswith("!"): continue
            if ("parameter" in line.lower()) and ("xlength" in line.lower()):
                return float(line.split()[3]) # Mm

def get_dt():
    codedir=os.path.dirname(os.path.abspath(__file__))
    paramsfile = os.path.join(codedir,"params.i")
    with open(paramsfile,'r') as paramsfile:
        for line in paramsfile.readlines():
            line=line.strip()
            if line.startswith("!"): continue
            if ("parameter" in line.lower()) and ("outputcad" in line.lower()):
                return float(line.split()[3].split(")")[0]) # seconds
    
def get_ridge_filter():
    codedir=os.path.dirname(os.path.abspath(__file__))
    driverfile = os.path.join(codedir,'driver.f90')
    with open(driverfile,'r') as driverfile:
        adjoint_src_filt=False
        ridge_filter=False
        for line in driverfile.readlines():
            if "SUBROUTINE ADJOINT_SOURCE_FILT" in line.upper(): adjoint_src_filt=True
            if "RIDGE FILTERS" in line.upper(): ridge_filter=True
            if adjoint_src_filt and "do pord=" in line.lower() and ridge_filter:
                ridges=line.split()[1].split("=")[1].split(",")
                if len(ridges)==2:
                    ridges=[str(i) for i in xrange(int(ridges[0]),int(ridges[1])+1)]
                elif len(ridges)==3:
                    ridges=[str(i) for i in xrange(int(ridges[0]),int(ridges[1])+1,int(ridges[2]))]
                return ridges
            if "PHASE-SPEED FILTERS" in line.upper(): ridge_filter=False
            if "END SUBROUTINE ADJOINT_SOURCE_FILT" in line.upper(): adjoint_src_filt=False

def get_nx():
    codedir=os.path.dirname(os.path.abspath(__file__))
    paramsfile = os.path.join(codedir,"params.i")
    with open(paramsfile,'r') as paramsfile:
        for line in paramsfile.readlines():
            line=line.strip()
            if line.startswith("!"): continue
            if ("parameter" in line.lower()) and ("nx" in line.lower()):
                return int(line.split()[3].split(',')[0])

def get_nz():
    codedir=os.path.dirname(os.path.abspath(__file__))
    paramsfile = os.path.join(codedir,"params.i")
    with open(paramsfile,'r') as paramsfile:
        for line in paramsfile.readlines():
            line=line.strip()
            if line.startswith("!"): continue
            if ("parameter" in line.lower()) and ("nz" in line.lower()):
                return int(line.split()[-1].strip(")"))

def get_solarmodel():
    codedir=os.path.dirname(os.path.abspath(__file__))
    paramsfile = os.path.join(codedir,"params.i")
    with open(paramsfile,'r') as paramsfile:
        for line in paramsfile:
            line=line.strip()
            if line.startswith("!"): continue
            if ("parameter" in line.lower()) and ("file_data" in line.lower()):
                return line.split("=")[-1].strip().strip("'")

def get_excitedepth():
    codedir=os.path.dirname(os.path.abspath(__file__))
    paramsfile = os.path.join(codedir,"params.i")
    with open(paramsfile,'r') as paramsfile:
        for line in paramsfile:
            line=line.strip()
            if line.startswith("!"): continue
            if ("parameter" in line.lower()) and ("excitdep" in line.lower()):
                return float(line.split("=")[-1].strip(")").replace(" ",""))


def get_enforced_continuity():
    codedir=os.path.dirname(os.path.abspath(__file__))
    paramsfile = os.path.join(codedir,"params.i")
    with open(paramsfile,'r') as paramsfile:
        for line in paramsfile:
            line=line.strip()
            if line.startswith("!"): continue
            if ("parameter" in line.lower()) and ("enf_cont" in line.lower()):
                if line.split("=")[1].split(")")[0].strip().strip(".").lower()=='true': return True
                elif line.split("=")[1].split(")")[0].strip().strip(".").lower()=='false': return False
    return "unknown"
                
def get_psi_continuity():
    codedir=os.path.dirname(os.path.abspath(__file__))
    paramsfile = os.path.join(codedir,"params.i")
    with open(paramsfile,'r') as paramsfile:
        for line in paramsfile:
            line=line.strip()
            if line.startswith("!"): continue
            if ("parameter" in line.lower()) and ("psi_cont" in line.lower()):
                if line.split("=")[1].split(")")[0].strip().strip(".").lower()=='true': return True
                elif line.split("=")[1].split(")")[0].strip().strip(".").lower()=='false': return False
    return "unknown"

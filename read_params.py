import os,fnmatch,sys
from pathlib import Path

codedir = Path(__file__).parent.absolute()
paramsfile = codedir/"params.i"
driverfile = codedir/'driver.f90'

def strip_comments(line):
    line=line.strip()
    head,_,_ = line.partition("!") # strip comments
    return head

def get_directory():
    
    with open(paramsfile,'r') as pf:
        for line in pf:

            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("directory" in line.lower()):
                return line.split()[-1].split("'")[1].rstrip("/")
    return None

def get_xlength():
    
    with open(paramsfile,'r') as pf:
        for line in pf:
            
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("xlength" in line.lower()):
                return float(line.split()[3]) # Mm

def get_dt():

    with open(paramsfile,'r') as pf:
        for line in pf:
            
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("outputcad" in line.lower()):
                return float(line.split()[3].split(")")[0]) # seconds

def get_dt_simulation():

    with open(paramsfile,'r') as pf:
        for line in pf:
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("timestep" in line.lower()):
                return float(line.split()[-1].split(")")[0]) # seconds

def get_total_simulation_time():

    with open(paramsfile,'r') as pf:
        for line in pf:
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("solartime" in line.lower()):
                return float(line.split()[3].split(")")[0]) # hours

def get_ridge_filter():
    
    with open(driverfile,'r') as df:
        adjoint_src_filt=False
        ridge_filter=False
        for line in df:
            if "SUBROUTINE ADJOINT_SOURCE_FILT" in line.upper(): adjoint_src_filt=True
            if "RIDGE FILTERS" in line.upper(): ridge_filter=True
            if adjoint_src_filt and "do pord=" in line.lower() and ridge_filter:
                ridges=line.split()[1].split("=")[1].split(",")
                if len(ridges)==2:
                    ridges=[str(i) for i in range(int(ridges[0]),int(ridges[1])+1)]
                elif len(ridges)==3:
                    ridges=[str(i) for i in range(int(ridges[0]),int(ridges[1])+1,int(ridges[2]))]
                return ridges
            if "PHASE-SPEED FILTERS" in line.upper(): ridge_filter=False
            if "END SUBROUTINE ADJOINT_SOURCE_FILT" in line.upper(): adjoint_src_filt=False

def get_nx():

    with open(paramsfile,'r') as pf:
        for line in pf:
            
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("nx" in line.lower()):
                return int(line.split()[3].split(',')[0])

def get_nz():
    
    with open(paramsfile,'r') as pf:
        for line in pf:
            
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("nz" in line.lower()):
                return int(line.split()[-1].strip(")"))

def get_solarmodel():

    with open(paramsfile,'r') as pf:
        for line in pf:
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("file_data" in line.lower()):
                return line.split("=")[-1].strip().strip("'")

def get_excitedepth():

    with open(paramsfile,'r') as pf:
        for line in pf:
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("excitdep" in line.lower()):
                return float(line.split("=")[-1].strip(")").replace(" ",""))

def get_obs_depth():
    
    with open(paramsfile,'r') as pf:
        for line in pf:
            
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("obsheight" in line.lower()):
                return float(line.split("=")[-1].strip(")").replace(" ",""))

def if_continuity_enforced():

    with open(paramsfile,'r') as pf:
        for line in pf:
            
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("enf_cont" in line.lower()):
                if line.split("=")[1].split(")")[0].strip().strip(".").lower()=='true': return True
                elif line.split("=")[1].split(")")[0].strip().strip(".").lower()=='false': return False
    return None

def get_continuity_variable():

    with open(paramsfile,'r') as pf:
        for line in pf:
            
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("psi_cont" in line.lower()):
                if line.split("=")[1].split(")")[0].strip().strip(".").lower()=='true': return 'psi'
            if ("parameter" in line.lower()) and ("vx_cont" in line.lower()):
                if line.split("=")[1].split(")")[0].strip().strip(".").lower()=='true': return 'vx'
            if ("parameter" in line.lower()) and ("vz_cont" in line.lower()):
                if line.split("=")[1].split(")")[0].strip().strip(".").lower()=='true': return 'vz'
    return "unknown"

def if_soundspeed_perturbed():
    with open(paramsfile,'r') as pf:
        for line in pf:
            
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("sound_speed_perturbation" in line.lower()):
                if line.split("=")[1].split(")")[0].strip().strip(".").lower()=='true': return True
                elif line.split("=")[1].split(")")[0].strip().strip(".").lower()=='false': return False
    return None

def if_flows():
    with open(paramsfile,'r') as pf:
        for line in pf:
            
            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("flows" in line.lower()):
                if line.split("=")[1].split(")")[0].strip().strip(".").lower()=='true': return True
                elif line.split("=")[1].split(")")[0].strip().strip(".").lower()=='false': return False
    return None

def get_modes_used():
    ridge_filters_driver=get_ridge_filter()
    datadir = get_directory()
    paramsfiles=[os.path.splitext(f)[1][1:] for f in fnmatch.filter(os.listdir(datadir),'params.[0-9]')]
    ridge_filters=[ridge for ridge in ridge_filters_driver if ridge in paramsfiles]
    return ridge_filters

def get_Lregular():

    with open(driverfile,'r') as df:
        psi_cont_check=False
        for line in df:
            if "if (psi_cont .and. enf_cont) then" in line.lower(): psi_cont_check=True
            if psi_cont_check and "Lregular" in line:
                if "Lregular" != line.strip().split()[0]: continue
                return float(line.strip().split("=")[-1].rstrip("/diml"))
            if "elseif (enf_cont .and. (vx_cont)) then" in line.lower():
                psi_cont_check=False
                break

def get_cutoff_dist():

    switch=False
    val = None

    with open(paramsfile,'r') as pf:
        for line in pf.readlines():

            line = strip_comments(line)
            if not len(line): continue # if a line is entirely commented out

            if ("parameter" in line.lower()) and ("cutoff_switch" in line.lower()):
                switch= line.split("=")[-1].split(")")[0].strip()
                if switch.upper() == '.TRUE.': switch =True
                elif switch.upper() == '.FALSE.': switch = False
            if ("parameter" in line.lower()) and ("cutoff_dist" in line.lower()):
                val= float(line.split("=")[-1].split(")")[0].strip())

    return switch,val



################################################################################################


def parse_cmd_line_params(key,mapto=None,default=None,return_list=False,zfill=None):
    cmd_line_param = [x for x in sys.argv if x.startswith(key+"=")]

    if len(cmd_line_param)!=0:
        retval=cmd_line_param[0].split("=")[-1]
        retval = retval.strip().lstrip('[').rstrip(']').split(',')
        if zfill is not None: retval=[x.zfill(zfill) for x in retval]
        if mapto is not None: retval=list(map(mapto,retval))
        retval = [None if element=='None' else element for element in retval]
        if not return_list and len(retval)==1:
            retval=retval[0]

    else: retval=default

    return retval

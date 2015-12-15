import os,read_params,sys


codedir=os.path.dirname(os.path.abspath(__file__))

datadir=read_params.get_directory()

try:
    with open(os.path.join(datadir,"master.pixels"),'r') as masterpixels:
        nmasterpixels=sum(1 for _ in masterpixels)
except IOError:
    print "No master.pixels file found"
    quit()

nlinesearch=int(sys.argv[1])

def create_if_not_there(path):
    if not os.path.exists(path):
        os.makedirs(path)

for src in xrange(1,nmasterpixels+1):

    srccode=str(src).zfill(2)
    for job in xrange(nlinesearch+1):
        jobcode=str(job).zfill(2) 
        forwarddir=os.path.join(datadir,"forward_src"+srccode+"_ls"+jobcode)
        create_if_not_there(forwarddir)

 
    adjointdir=os.path.join(datadir,"adjoint_src"+srccode)
    create_if_not_there(adjointdir)
    
    
kerneldir=os.path.join(datadir,"kernel")
create_if_not_there(kerneldir)

updatedir=os.path.join(datadir,"update")
create_if_not_there(updatedir)

statusdir=os.path.join(datadir,"status")
create_if_not_there(statusdir)

ttdir=os.path.join(datadir,"tt")
create_if_not_there(ttdir)

datadir=os.path.join(datadir,"data")
create_if_not_there(datadir)

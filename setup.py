import os

directory="/scratch/shivam/flows/data"
with open(os.path.join(directory,"master.pixels"),'r') as masterpixels:
    nmasterpixels=sum(1 for _ in masterpixels)

nlinesearch=5

def create_if_not_there(path):
    if not os.path.exists(path):
        os.makedirs(path)

for src in xrange(1,nmasterpixels+1):

    srccode=str(src).zfill(2)
    for job in xrange(nlinesearch+1):
        jobcode=str(job).zfill(2) 
        forwarddir=os.path.join(directory,"forward_src"+srccode+"_ls"+jobcode)
        create_if_not_there(forwarddir)

 
    adjointdir=os.path.join(directory,"adjoint_src"+srccode)
    create_if_not_there(adjointdir)
    
    
kerneldir=os.path.join(directory,"kernel")
create_if_not_there(kerneldir)

updatedir=os.path.join(directory,"update")
create_if_not_there(updatedir)

statusdir=os.path.join(directory,"status")
create_if_not_there(statusdir)

ttdir=os.path.join(directory,"tt")
create_if_not_there(ttdir)

datadir=os.path.join(directory,"data")
create_if_not_there(datadir)

import os

def create_if_not_there(path):
    if not os.path.exists(path):
        os.makedirs(path)

def create_directories(datadir,num_src,num_ls):
    
    for src in range(1,num_src+1):

        for job in range(num_ls+1):
        
            forwarddir=os.path.join(datadir,"forward_src{:02d}_ls{:02d}".format(src,job))
            create_if_not_there(forwarddir)
     
        adjointdir=os.path.join(datadir,"adjoint_src{:02d}".format(src))
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

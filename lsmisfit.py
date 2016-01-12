import sys,os,glob,re
import numpy as np
import read_params

codedir=os.path.dirname(os.path.abspath(__file__))

datadir=read_params.get_directory()

number_args=filter(lambda x: x.isdigit(),sys.argv)
if number_args:
    iterno=number_args[0].zfill(2)
else:
    lsfiles=[f for f in glob.glob(os.path.join(datadir,"update","linesearch_*")) if "all" not in os.path.basename(f)]
    nfiles=len(lsfiles)
    if nfiles==0:
        print "No linesearch files found"
        quit()
    else:
        iterno=str(nfiles-1).zfill(2)

no_of_linesearches=5

lsfile=os.path.join(datadir,"update","linesearch_"+iterno)

if not os.path.exists(lsfile):
    print lsfile,"doesn't exist"
    quit()

with open(os.path.join(datadir,'master.pixels'),'r') as mpixfile:
    nmasterpixels=sum(1 for _ in mpixfile)

lsdata=np.loadtxt(lsfile,usecols=[2])

misfit=[sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels]) for i in xrange(no_of_linesearches)]

print "iteration",int(iterno)

np.set_printoptions(precision=3)
    

if "--detail" in sys.argv:
    ridges = read_params.get_modes_used()
    modes={'0':'fmode'}
    for pmodeno in xrange(1,6): modes.update({str(pmodeno):'p'+str(pmodeno)+'mode'})
    
    for src in xrange(1,nmasterpixels+1):
        for lsno in xrange(1,no_of_linesearches+1):
            modemisfit=[]
            maxpix=[]
            for mode in ridges:
                pix,td=np.loadtxt(os.path.join(datadir,"forward_src"+str(src).zfill(2)+"_ls"+str(lsno).zfill(2),"ttdiff."+mode),
                            usecols=[0,1],unpack=True)
                maxpix.append(int(pix[abs(td).argmax()]))
                td=0.5*np.sum((td/60)**2)
                modemisfit.append(td)
            misfitfmtstr="{:6.3f} "*len(modemisfit)
            maxpixfmtstr="{:4d}"*len(maxpix)
            print "Src",src,"ls",lsno,"misfit",misfitfmtstr.format(*modemisfit),\
                    "sum","{:6.3f}".format(sum(modemisfit)),"maxpix",maxpixfmtstr.format(*maxpix)
            
    for src,srcmisfit in enumerate(lsdata.reshape(no_of_linesearches,lsdata.shape[0]//no_of_linesearches).T):
        print "Source",str(src+1).zfill(2),"lsmisfit",srcmisfit
    print "Total misfit      ",np.array(misfit)
else:
    for m in misfit: print m

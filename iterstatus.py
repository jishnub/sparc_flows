import os,read_params,glob

datadir=read_params.get_directory()
print datadir

if not os.path.exists(datadir):
    print datadir,"doesn't exist. Modify params.i"
    quit()

if not os.path.exists(os.path.join(datadir,"update")):
    print "No iterations done"
    print "qsub data_forward.sh"
    quit()

lsfiles=sorted([int(f[-2:]) for f in glob.glob(os.path.join(datadir,"update","linesearch_*")) if "all" not in f])
misfitfiles=sorted([int(f[-2:]) for f in glob.glob(os.path.join(datadir,"update","misfit_*")) if "all" not in f])

if len(misfitfiles)==0:
    print "Zero iterations done."
    print "qsub full.sh"
    quit()

if len(misfitfiles)>len(lsfiles):
    print "Forward computation for iteration",misfitfiles[-1],"done"
    print "python grad.py && qsub linesearch.sh"
    quit()

if len(misfitfiles)==len(lsfiles):
    print lsfiles[-1],"iterations complete"
    print "python lsmisfit.py"
    quit()





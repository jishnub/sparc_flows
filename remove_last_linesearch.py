import os,glob,read_params,shutil

datadir=read_params.get_directory()
lsfiles=sorted([int(f[-2:]) for f in glob.glob(os.path.join(datadir,"update","linesearch_*")) if "all" not in f])
misfitfiles=sorted([int(f[-2:]) for f in glob.glob(os.path.join(datadir,"update","misfit_*")) if "all" not in f])

last_ls_no=str(lsfiles[-1]).zfill(2)
last_misfit_no=str(misfitfiles[-1]).zfill(2)

ls_file=os.path.join(datadir,"update","linesearch_"+last_ls_no)
ls_file_new=os.path.join(datadir,"update","ls_"+last_ls_no+".rnm")
ls_all_file=os.path.join(datadir,"update","linesearch_all_"+last_ls_no)
ls_all_file_new=os.path.join(datadir,"update","ls_all_"+last_ls_no+".rnm")

iterstatus_good=(last_misfit_no==last_ls_no)

try:
    if iterstatus_good:
        shutil.move(ls_file,ls_file_new)
except IOError:
    print "Can't move",ls_file
    
try:
    if iterstatus_good:
        shutil.move(ls_all_file,ls_all_file_new)
except IOError:
    print "Can't move",ls_all_file



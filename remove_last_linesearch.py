import os,glob,read_params,shutil,fnmatch

datadir=read_params.get_directory()
lsfiles=[os.path.join(datadir,"update",i) for i in
        sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),"linesearch_[0-9][0-9]"))]
ls_all_files=[os.path.join(datadir,"update",i) for i in
        sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),"linesearch_all_[0-9][0-9]"))]
misfitfiles=[os.path.join(datadir,"update",i) for i in
        sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),"misfit_[0-9][0-9]"))]

last_ls_no=os.path.basename(lsfiles[-1]).split("_")[-1]
last_misfit_no=os.path.basename(misfitfiles[-1]).split("_")[-1]




iterstatus_good=(last_misfit_no==last_ls_no)

try:
    if iterstatus_good:
        ls_file=lsfiles[-1]
        ls_file_new=os.path.join(datadir,"update","ls_"+last_ls_no+".rnm")
        shutil.move(ls_file,ls_file_new)
except IOError:
    print("Can't move",ls_file)

try:
    if iterstatus_good and ls_all_files:
        ls_all_file=ls_all_files[-1]
        ls_all_file_new=os.path.join(datadir,"update","ls_all_"+last_ls_no+".rnm")
        shutil.move(ls_all_file,ls_all_file_new)
except IOError:
    print("Can't move",ls_all_file)

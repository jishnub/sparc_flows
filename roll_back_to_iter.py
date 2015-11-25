import os,sys,glob,read_params,shutil,datetime

datadir=read_params.get_directory()

def iterno(filename):
    return int(filter(lambda x: x.isdigit(),os.path.splitext(os.path.basename(filename))[0].split('_'))[0])

misfitfiles=sorted(glob.glob(os.path.join(datadir,"update","misfit_*")))
last_misfit_no=iterno(misfitfiles[len(misfitfiles)/2-1])

iter_stage=last_misfit_no

#~ Get iter no to roll back to

try:
    iter_rollback=int(sys.argv[1])
except IndexError:
    print "Usage: python roll_back_to_iter.py <iterno>"
    quit()

#~ Check if the stage has been reached and surpassed

if iter_rollback>=iter_stage: 
    print "Currently at iteration",iter_stage,"can't roll back to iter",iter_rollback
    quit()

#~ Create a directory to back up old files

dest_dir=os.path.join(datadir,"update","iters_"+str(iter_rollback+1).zfill(2)+
                "_to_"+str(iter_stage).zfill(2)+"_created_on_"+str(datetime.datetime.now().date()))

try:
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
except OSError:
    print "Could not create path",dest_dir
    
#~ Move old files

def move_files(files_to_move):
    dest_files=map(lambda x: os.path.join(dest_dir,os.path.basename(x)),files_to_move)
    for f in zip(files_to_move,dest_files):
        try:
            shutil.move(f[0],f[1])
        except IOError:
            print "Can't move",f[0],"to",f[1]

#~ model

model_files=sorted(glob.glob(os.path.join(datadir,"update","model_*")))
model_c=model_files[:len(model_files)/2]
model_psi=model_files[len(model_files)/2:]

next_model=filter(lambda x: iterno(x)==iter_rollback+1,model_c)
if next_model:
    try:
        shutil.copyfile(next_model,os.path.join(datadir,"model_psi_ls00.fits"))
    except IOError:
        print "Can't copy model for iteration",iter_rollback+1

files_to_move=filter(lambda x: iterno(x)>iter_rollback,model_c)
move_files(files_to_move)

files_to_move=filter(lambda x: iterno(x)>iter_rollback,model_psi)
move_files(files_to_move)

#~ gradient

gradient_files=sorted(glob.glob(os.path.join(datadir,"update","gradient_*")))
gradient_c=gradient_files[:len(gradient_files)/2]
gradient_psi=gradient_files[len(gradient_files)/2:]

files_to_move=filter(lambda x: iterno(x)>iter_rollback,gradient_c)
move_files(files_to_move)

files_to_move=filter(lambda x: iterno(x)>iter_rollback,gradient_psi)
move_files(files_to_move)

#~ update

update_files=sorted(glob.glob(os.path.join(datadir,"update","update_*")))
update_c=update_files[:len(update_files)/2]
update_psi=update_files[len(update_files)/2:]

files_to_move=filter(lambda x: iterno(x)>iter_rollback,update_c)
move_files(files_to_move)

files_to_move=filter(lambda x: iterno(x)>iter_rollback,update_psi)
move_files(files_to_move)

#~ misfit

misfit_files=sorted(glob.glob(os.path.join(datadir,"update","misfit_*")))
misfit=misfit_files[:len(misfit_files)/2]
misfit_all=misfit_files[len(misfit_files)/2:]

files_to_move=filter(lambda x: iterno(x)>iter_rollback,misfit)
move_files(files_to_move)

files_to_move=filter(lambda x: iterno(x)>iter_rollback,misfit_all)
move_files(files_to_move)

#~ linesearch

ls_files=sorted(glob.glob(os.path.join(datadir,"update","linesearch_*")))
ls=ls_files[:len(misfit_files)/2]
ls_all=ls_files[len(misfit_files)/2:]

files_to_move=filter(lambda x: iterno(x)>iter_rollback,ls)
move_files(files_to_move)

files_to_move=filter(lambda x: iterno(x)>iter_rollback,ls_all)
move_files(files_to_move)

#~ velocity

v_files=sorted(glob.glob(os.path.join(datadir,"update","v*.fits")))
vx_files=v_files[:len(v_files)/2]
vz_files=v_files[len(v_files)/2:]

files_to_move=filter(lambda x: iterno(x)>iter_rollback,vx_files)
move_files(files_to_move)

files_to_move=filter(lambda x: iterno(x)>iter_rollback,vz_files)
move_files(files_to_move)


        


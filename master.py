from __future__ import division
import subprocess ,time, os, shutil,read_params,glob,fnmatch,sys
import numpy as np
import matplotlib.pyplot as plt
from termcolor import colored

np.set_printoptions(precision=5)

datadir=read_params.get_directory()

if read_params.if_flows():
    iter_var="psi"
elif read_params.if_soundspeed_perturbed():
    iter_var="c"

num_linesearches = 6

data_command = "qsub data_forward.sh"
full_command = "qsub full.sh"
ls_command = "qsub linesearch.sh"
def grad_command(algo="cg",eps=[0.1*i for i in xrange(1,num_linesearches+1)]):
    eps = [str(i) for i in eps]
    eps_str = ' '.join(eps)
    return "python grad.py algo="+algo+" "+str(eps_str)

if len(sys.argv[1:])!=1:
    print "Usage: python master.py <id>"
    exit()
id_text = sys.argv[1]

def update_id_in_file(filename):
    with open(filename,'r') as f:
        code=f.readlines()
    
    for lineno,line in enumerate(code):
        if "PBS" in line and "-N" in line:
            current_id = '_'.join(line.split()[-1].split('_')[1:])
            if current_id!=id_text:
                print "Current identifier in ",filename,"is",current_id
                print "Changing to ",id_text
                line=line.replace(current_id,id_text)
                code[lineno]=line
                
    with open(filename,'w') as f:
        f.writelines(code)
                
update_id_in_file('data_forward.sh')
update_id_in_file('full.sh')
update_id_in_file('linesearch.sh')

def copy_updated_model(num,var="psi"):
    updated_model=os.path.join(datadir,'update','test_'+var+'_'+str(num)+'.fits')
    model_for_next_iter=os.path.join(datadir,'model_'+var+'_ls00.fits')
    shutil.copyfile(updated_model,model_for_next_iter)

if os.path.exists('compute_data'):
    os.remove('compute_data')

if os.path.exists('linesearch'):
    os.remove('linesearch')

if os.path.exists('running_full'):
    os.remove('running_full')
    
running_ls = False
running_full = False

def get_eps_around_minimum(prev1_eps,prev1_misfits):
    
    #~ Previous iteration, save the plot
    plt.plot(prev1_eps,prev1_misfits,'bo',zorder=1)
    plt.xlabel("Step size",fontsize=20)
    plt.ylabel(r"misfit",fontsize=20)
    plt.savefig('step_size_selection.eps')
    
    p=np.polyfit(prev1_eps,prev1_misfits,2)

    assert (not p[0]<0), colored("inverted parabolic fit, can't estimate step size","red")
    step = 1e-2
    if p[0]>0:
        step = - p[1]/(2*p[0])
    elif p[0]==0:
        step = -p[2]/p[1]
    assert step>0,colored("ls minimum predicted at less than zero step","red")
    
    #~ Choose points around the predicted minimum stepsize, and compute predicted misfits
    #~ according to the previous itreation's trend
    eps = np.array([step*(0.8+0.1*i) for i in xrange(num_linesearches)])
    pred_misfit = np.polyval(p,eps)
    print "Estimated step size",eps
    print "Predicted misfit",pred_misfit
    plt.plot(eps,pred_misfit,'ro',zorder=1)
    
    #~ Plot the fit polynomial line, this should pass through the points and have a minimum
    #~ around the predicted points
    epsfine = np.linspace(min(eps[0],prev1_eps[0]),max(eps[-1],prev1_eps[-1]),num=2000)
    misfitfine = np.polyval(p,epsfine)
    plt.plot(epsfine,misfitfine,'g-',zorder=0)
    
    plt.savefig('step_size_selection.eps')
    plt.clf()
    
    return eps

def get_iter_status():
    ls_files=sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'linesearch_[0-9][0-9]'))
    misfit_files=sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'misfit_[0-9][0-9]'))
    ls_all_files=sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'linesearch_all_[0-9][0-9]'))
    misfit_all_files=sorted(fnmatch.filter(os.listdir(os.path.join(datadir,"update")),'misfit_all_[0-9][0-9]'))
    iterno = len(misfit_files)-1
    return misfit_files,misfit_all_files,ls_files,ls_all_files,iterno

epslist_path=os.path.join(datadir,'epslist.npz')

for query in xrange(100000):

    #~ print "Query",query
    qstat=subprocess.check_output(["qstat"]).split()
    #~ print qstat
    something_running = False
    for part in qstat:
        if part.endswith(id_text):
            #~ print "Query at",time.strftime(" %d %b, %X", time.localtime()),":",\
                    #~ part,"running, moving on without doing anything"
            something_running = True
            break
            
    if something_running: 
        time.sleep(60)
        continue
    #~ else: print "Nothing running, moving on"
    
    if not os.path.exists(os.path.join(datadir,"forward_src01_ls00","data.fits")):
        #~ no iterations done
        print colored("Running data_forward","blue")
        status=subprocess.call(data_command.split())
        assert status==0,"Error in running data forward"
        time.sleep(30)
        continue
        
    misfit_files,misfit_all_files,ls_files,ls_all_files,iterno = get_iter_status()
    #~ print misfit_files,ls_files
    print colored("Iteration "+str(iterno),"green")
    
    if iterno==-1:
        #~ Start of iterations
        print "#"*80
        print colored("Starting with iterations, running full for the first time","blue")
        status=subprocess.call(full_command.split())
        assert status==0,"Error in running full"
        time.sleep(30)
        continue
        
    elif len(misfit_files)>len(ls_files):
        
        #~ Kernel computation should be over for this iteration
        #~ We should proceed to computing update and running linesearch
        
        running_full = False
        if running_ls:
            print colored("It seems the linesearch file wasn't generated. Check if previous linesearch finished correctly","red")
            exit()
        print "Misfit from previous iteration",np.sum(np.loadtxt(os.path.join(datadir,"update",misfit_files[-1]),usecols=[2]))
        
        #~ Need to run linesearch for this iteration
        #~ check if eps for iteration exists
        try: 
            epslist=np.load(epslist_path)
            iters_list = epslist.files
            
            if str(iterno) in iters_list:
                #~ If it exists, this means linesearch files were removed
                #~ This probably means that the linesearch misfits were monotonic
                
                ls_details = epslist[str(iterno)]
                
                #~ print ls_details
                
                prev1_done = ls_details[1].all()
                prev2_done = ls_details[3].all()
                    
                if prev1_done:
                    print "Computing step size from previous linesearch"
                    eps_prev = ls_details[0]
                    lsmisfit_prev =  ls_details[1]
                    eps = get_eps_around_minimum(eps_prev,lsmisfit_prev)
                
                elif not(prev1_done) and prev2_done:
                    print "Computing step size form 2nd last linesearch"
                    eps_prev = ls_details[2]
                    lsmisfit_prev =  ls_details[3]
                    eps = get_eps_around_minimum(eps_prev,lsmisfit_prev)
                    
                elif (not prev1_done) and (not prev2_done):
                    #~ Probably interrupted, so misfits not saved but step sizes might be recorded
                    #~ Restart with one set of step sizes
                    if ls_details[0].any():
                        print "Reusing previous step size"
                        eps = ls_details[0]
                    elif ls_details[2].any():
                        print "Reusing step size from 2nd last iteration"
                        eps = ls_details[2]
                    else:
                        print "No step size information recorded",
                        prev_iternos= np.array([int(i) for i in iters_list if i!=str(iterno)])
                        if len(prev_iternos)!=0:
                            closest_iter = str(prev_iternos[abs(prev_iternos-iterno).argmin()])
                            eps=epslist[closest_iter][0]
                            print "using step size from iteration",closest_iter
                        else:
                            print "using arbitary step sizes"
                            eps=[1e-2*i for i in xrange(1,num_linesearches+1)]

            else:
                print "Linesearch for this iteration not done yet, choosing step size from iteration",
                #~ Iteration not yet done
                #~ Choose the value corresponding to the closest iteration
                prev_iternos= np.array([int(i) for i in iters_list])
                closest_iter = str(prev_iternos[abs(prev_iternos-iterno).argmin()])
                print closest_iter
                eps=epslist[closest_iter][0]
                
                
        except IOError: 
            #~ If the file doesn't exist, no linesearches have been carried out
            #~ Assign arbitrary small value
            print "Using arbitary step sizes"
            eps=np.array([1e-2*i for i in xrange(1,num_linesearches+1)])
        #~ Run grad
        gradcmd = grad_command(algo="cg",eps=eps)
        print "Running",gradcmd
        status=subprocess.call(gradcmd.split())
        assert status==0,"Error in running grad"
        
        #~ Run linesearch
        print colored("Running linesearch","blue")
        status=subprocess.call(ls_command.split())
        assert status==0,"Error in running linesearch"
        running_ls = True
        time.sleep(30)
        continue
        
    elif iterno>=0 and len(misfit_files)==len(ls_files):

        #~ Linesearch should be over, and we can do two things - 
        #~ (a) run linesearch again with different step size if misfit minimum is not found
        #~ (b) run full if misfit minimum is found
        
        running_ls = False
        if running_full:
            print colored("It seems the full misfit was not generated. Check if the previous full.sh finished without errors","red")
            exit()
        #~ check linesearch misfits
        ls_latest = os.path.join(datadir,"update",ls_files[-1])
        ls_all_latest = os.path.join(datadir,"update",ls_all_files[-1])
        
        try:
            lsdata = np.loadtxt(ls_latest)
        except IOError:
            print colored("Could not load "+ls_latest+", check if file exists","red")
            exit()
        
        assert len(np.where(lsdata[:,0]==1)[0])==num_linesearches,"Not all linesearches finished correctly"
        nmasterpixels = lsdata.shape[0]//num_linesearches
        misfit=np.array([sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels,2]) for i in xrange(num_linesearches)])
        
        print "Linesearch misfits",misfit
        #~ Check for misfit minimum
        misfit_diff = np.diff(misfit)
        misfit_diff_sign_change_locs = np.where(np.diff(np.sign(misfit_diff)))[0]
        misfit_diff_changes_sign = len(misfit_diff_sign_change_locs)>0
        
        #~ Update epslist
        epslist=np.load(epslist_path)
        if str(iterno) in epslist.files:
            ls_details = epslist[str(iterno)]
            ls_details[1] = misfit
        else:
            ls_details = np.zeros_like(epslist[epslist.files[0]])
        epslist = {i:epslist[i] for i in epslist.files}
        epslist.update({str(iterno):ls_details})
        #~ print epslist[str(iterno)]
        np.savez(epslist_path,**epslist)

        #~ Copy best linesearch, or remove prev linesearch file
        if misfit_diff_changes_sign:
            #~ Get index of minimum. Add two ones to the minimum of the diff array
            #~ One because the sign change occurs at the index before the minimum
            #~ The other because the test model counting starts from 1 instead of 0
            min_ls_index = misfit_diff_sign_change_locs[0]+1+1
            print "Misfit sign change detected"
            print "model #"+str(min_ls_index)+" has minimum misfit"
            copy_updated_model(min_ls_index,var=iter_var)
            
            plt.clf() # Refresh the stepsize-misfit plot
            
            #~ Run full.sh
            print "#"*80
            print colored("Running full","blue")
            status=subprocess.call(full_command.split())
            running_full=True
            assert status==0,"Error in running full"
            time.sleep(30)
        else:
            if misfit_diff[0]>0: print "Linesearch strictly increasing"
            elif misfit_diff[0]<0: print "Linesearch strictly decreasing"
            
            #~ Monotonic linesearches don't count, remove previous linesearch file
            print "Removing linesearch file"

            ls_latest_renamed=ls_latest+".rnm"
            os.rename(ls_latest,ls_latest_renamed)
            
            ls_all_latest_renamed= ls_all_latest+".rnm"
            os.rename(ls_all_latest,ls_all_latest_renamed)
        
        continue

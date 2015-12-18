from __future__ import division
import subprocess ,time, os, shutil,read_params,glob
import numpy as np

np.set_printoptions(precision=5)

datadir=read_params.get_directory()

num_linesearches = 6

data_command = "qsub data_forward.sh"
full_command = "qsub full.sh"
ls_command = "qsub linesearch.sh"
def grad_command(eps):
    eps = [str(i) for i in eps]
    eps_str = ' '.join(eps)
    return "python grad.py "+str(eps_str)

id_text = "frequent"

ls_strict_increase=False
ls_strict_decrease=False

def safecopy(a,b):
    try: shutil.copyfile(a,b)
    except IOError as e:
        sys.stderr.write("Could not copy "+a+" to "+b+"; "+e.args(1)+"\n")
        sys.stderr.flush()

def copy_updated_model(num):
    updated_model_vectorpsi=os.path.join(datadir,'update','test_psi_'+str(num)+'.fits')
    model_vectorpsi=os.path.join(datadir,'model_psi_ls00.fits')
    safecopy(updated_model_vectorpsi,model_vectorpsi)

if os.path.exists('compute_data'):
    os.remove('compute_data')

if os.path.exists('linesearch'):
    os.remove('linesearch')

if os.path.exists('running_full'):
    os.remove('running_full')
    
bisection = False
running_ls = False
running_full = False

def get_eps_around_minimum(prev1_eps,prev1_misfits,prev2_eps=None,prev2_misfits=None):
    p=np.polyfit(prev1_eps,prev1_misfits,2)
    step = 1e-2
    if p[0]!=0:
        step = - p[1]/(2*p[0])
    elif p[0]==0:
        step = -p[2]/p[1]
    assert step>0,"ls minimum predicted at less than zero step"
    eps = np.array([step*(0.8+0.1*i) for i in xrange(num_linesearches)])
    return eps

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
    
    if not os.path.exists(os.path.join(datadir,"forward_src01_ls00","data.fits")):
        #~ no iterations done
        print "Running data_forward"
        status=subprocess.call(data_command.split())
        assert status==0,"Error in running data forward"
        time.sleep(30)
        continue
        
    ls_files=sorted([f for f in glob.glob(os.path.join(datadir,"update","linesearch_*")) if "all" not in os.path.basename(f)])
    misfitfiles=sorted([f for f in glob.glob(os.path.join(datadir,"update","misfit_*")) if "all" not in os.path.basename(f)])
    ls_all_files=sorted([f for f in glob.glob(os.path.join(datadir,"update","linesearch_*")) if "all" in os.path.basename(f)])
    misfit_all_files=sorted([f for f in glob.glob(os.path.join(datadir,"update","misfit_*")) if "all" in os.path.basename(f)])
    iterno = int(len(misfitfiles))
    print "Iteration",iterno
    
    if len(misfitfiles)==0:
        #~ Start of iterations
        print "Running full"
        status=subprocess.call(full_command.split())
        assert status==0,"Error in running full"
        time.sleep(30)
        continue
        
    elif len(misfitfiles)>len(ls_files):
        running_full = False
        if running_ls:
            print "It seems the linesearch file wasn't generated. Check if previous linesearch finished correctly"
            exit()
        #~ Need to run linesearch for this iteration
        #~ check if eps for iteration exists
        try: 
            epslist=np.load("epslist.npz")
            iters_list = epslist.files
            #~ print epslist
            if str(iterno) in iters_list:
                #~ If it exists, this means linesearch files were removed
                #~ This probably means that the linesearch misfits were monotonic
                #~ Start with a value smaller or larger than the previous attempt
                ls_details = epslist[str(iterno)]
                
                prev1_done = ls_details[1].all()
                prev2_done = ls_details[3].all()
                two_iters_done = prev1_done and prev2_done
                    
                full_misfit = sum(np.loadtxt(misfitfiles[-1],usecols=[2]))
                
                if prev1_done and (not prev2_done):
                    #~ This means that there has been only one iteration
                    #~ Check the linesearch misfits from that iteration, 
                    #~ and try to guess a good value for the next iteration
                    
                    eps_prev = ls_details[0]
                    lsmisfit_prev =  ls_details[1]
                    eps = get_eps_around_minimum(eps_prev,lsmisfit_prev)
                    
                elif prev1_done and prev2_done:
                    #~ This means that there have been at least two attempts at this iterations
                    #~ It's possible that both the iterations resulted in an increase
                    #~ Or it's possible that the first iteration was a decrease
                    #~ and the latest attempt resulted in an increase
                    
                    if bisection:
                        eps_prev = ls_details[0]
                        lsmisfit_prev =  ls_details[1]
                        eps = get_eps_around_minimum(eps_prev,lsmisfit_prev)
                    else:
                        #~ One iterations resulted in decrease, but the latest was increasing
                        #~ In this case, try an intermediate value (bisection method)
                        #~ This will decrease the step size
                        bisection = True
                        eps = (ls_details[0]+ls_details[2])/2
                
                elif not(prev1_done) and prev2_done:
                    eps_prev = ls_details[2]
                    lsmisfit_prev =  ls_details[3]
                    eps = get_eps_around_minimum(eps_prev,lsmisfit_prev)
                            
                    
                elif (not prev1_done) and (not prev2_done):
                    #~ Probably interrupted, restart with one set of values
                    if ls_details[0].any():
                        eps = ls_details[0]
                    elif ls_details[2].any():
                        eps = ls_details[0]
                    else:
                        prev_iternos= np.array([int(i) for i in iters_list if i!=str(iterno)])
                        closest_iter = str(prev_iternos[abs(prev_iternos-iterno).argmin()])
                        eps=epslist[closest_iter][0]

            else:
                print "Linesearch for this iteration not done yet, choosing step size from nearest iteration"
                #~ Iteration not yet done
                #~ Choose the value corresponding to the closest iteration
                prev_iternos= np.array([int(i) for i in iters_list])
                closest_iter = str(prev_iternos[abs(prev_iternos-iterno).argmin()])
                eps=epslist[closest_iter][0]
                
                
        except IOError: 
            #~ If the file doesn't exist, no linesearches have been carried out
            #~ Assign arbitrary small value
            eps=np.array([1e-2*i for i in xrange(1,7)])
        #~ Run grad
        print "Running grad with eps",eps
        gradcmd = grad_command(eps)
        status=subprocess.call(gradcmd.split())
        assert status==0,"Error in running grad"
        
        #~ Run linesearch
        print "Running linesearch"
        status=subprocess.call(ls_command.split())
        assert status==0,"Error in running linesearch"
        running_ls = True
        time.sleep(30)
        continue
        
    elif len(misfitfiles)>0 and len(misfitfiles)==len(ls_files):
        running_ls = False
        if running_full:
            print "It seems the full misfit was not generated. Check if the previous full.sh finished without errors"
            exit()
        #~ check linesearch misfits
        ls_latest = os.path.join(datadir,"update","linesearch_"+str(iterno-1).zfill(2))
        try:
            lsdata = np.loadtxt(ls_latest)
        except IOError:
            print "Could not load",ls_latest , "check if file exists"
            exit()
        
        no_of_linesearches = len(np.where(lsdata[:,0]==1)[0])
        nmasterpixels = lsdata.shape[0]//no_of_linesearches
        misfit=np.array([sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels,2]) for i in xrange(no_of_linesearches)])
        
        print "Linesearch misfits"
        print misfit
        #~ Check for misfit minimum
        misfit_diff = np.diff(misfit)
        misfit_diff_sign_change_locs = np.where(np.diff(np.sign(misfit_diff)))[0]
        misfit_diff_changes_sign = len(misfit_diff_sign_change_locs)>0
        
        #~ Update epslist
        epslist=np.load('epslist.npz')
        if str(iterno) in epslist.files:
            ls_details = epslist[str(iterno)]
            ls_details[1] = misfit
        else:
            ls_details = np.zeros_like(epslist[epslist.files[0]])
        epslist = {i:epslist[i] for i in epslist.files}
        epslist.update({str(iterno):ls_details})
        np.savez('epslist.npz',**epslist)

        #~ Copy best linesearch, or remove prev linesearch file
        if misfit_diff_changes_sign:
            #~ Get index of minimum. Add two ones to the minimum of the diff array
            #~ One because the sign change occurs at the index before the minimum
            #~ The other because the test model counting starts from 1 instead of 0
            min_ls_index = misfit_diff_sign_change_locs[0]+1+1
            print "Misfit sign change detected"
            print "model #"+str(min_ls_index)+" has minimum misfit"
            copy_updated_model(min_ls_index)
            bisection = False
            
            #~ Run full.sh
            print "Running full"
            status=subprocess.call(full_command.split())
            assert status==0,"Error in running full"
            running_full=True
            time.sleep(30)
        else:
            if misfit_diff[0]>0:
                print "Linesearch strictly increasing"
                ls_strict_increase = True
                ls_strict_decrease = False
            elif misfit_diff[0]<0:
                print "Linesearch strictly decreasing"
                ls_strict_increase = False
                ls_strict_decrease = True
            
            #~ Monotonic linesearches don't count, remove previous linesearch file
            print "Removing linesearch file"
            last_ls_no=str(iterno-1).zfill(2)
            ls_file_new=os.path.join(datadir,"update","ls_"+last_ls_no+".rnm")
            ls_all_file_new=os.path.join(datadir,"update","ls_all_"+last_ls_no+".rnm")
            os.rename(ls_files[-1],ls_file_new)
            os.rename(ls_all_files[-1],ls_all_file_new)
        
        continue

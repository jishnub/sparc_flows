from __future__ import division
import subprocess ,time, os, shutil,read_params,glob
import numpy as np

np.set_printoptions(precision=3)

datadir=read_params.get_directory()

data_command = "qsub data_forward.sh"
full_command = "qsub full.sh"
ls_command = "qsub linesearch.sh"
def grad_command(eps): return "python grad.py "+str(eps)

id_text = "master"

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

for query in xrange(100000):

    #~ print "Query",query
    qstat=subprocess.check_output(["qstat"]).split()
    #~ print qstat
    something_running = False
    for part in qstat:
        if part.endswith(id_text):
            print "Query at",time.strftime(" %d %b, %X", time.localtime()),":",\
                    part,"running, moving on without doing anything"
            something_running = True
            break
            
    if something_running: 
        time.sleep(60)
        continue
    
    if not os.path.exists(os.path.join(datadir,"forward_src01_ls00","data.fits")):
        #~ no iterations done
        print "Running data_forward"
        subprocess.call(data_command.split())
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
        subprocess.call(full_command.split())
        time.sleep(30)
        continue
        
    elif len(misfitfiles)>len(ls_files):
        #~ Need to run linesearch for this iteration
        #~ check if eps for iteration exists
        try: 
            epslist=np.load("epslist.npz")['epslist']
            #~ print epslist
            if iterno in epslist[:,0]:
                #~ If it exists, this means linesearch files were removed
                #~ This probably means that the linesearch misfits were monotonic
                #~ Start with a value smaller or larger than the previous attempt
                matchindex=np.where(epslist[:,0]==iterno)[0][0]
                if ls_strict_increase:
                    print "Linesearch strictly increasing"
                    #~ Since the linesearch is increasing, decrease the step size
                    #~ Let's get to a regime where the linesearch is decreasing
                    if epslist[matchindex,2]==0:
                        #~ This means that there has been only one iteration
                        #~ We don't have any idea about the appropriate step size
                        #~ Decrease it by large amounts
                        print "No clue about step size, just decrease the previous step a lot"
                        eps=epslist[matchindex,1]/5
                    else:
                        #~ This means that there have been at least two attempts at this iterations
                        #~ It's possible that both the iterations resulted in an increase
                        #~ Or it's possible that the first iteration was a decrease
                        #~ and the latest attempt resulted in an increase
                        if epslist[matchindex,1]<epslist[matchindex,2]:
                            #~ Both iterations are increasing, decrease it further
                            #~ This might be a step in bisection, or it might be a step in the dark
                            #~ In case it's a step in the dark, we can keep decreasing the step by large amounts
                            #~ If it's bisection, we should decrease it by small amounts
                            if bisection:
                                print "Bisection on, but strict increase noted"
                                eps=epslist[matchindex,1]/1.5
                            else:
                                print "Strict increase, but no idea about correct step"
                                eps=epslist[matchindex,1]/3
                        else:
                            #~ One iterations resulted in decrease, but the latest was increasing
                            #~ In this case, try an intermediate value (bisection method)
                            #~ This will decrease the step size
                            bisection = True
                            print "Starting bisection"
                            eps = (epslist[matchindex,1]+epslist[matchindex,2])/2
                elif ls_strict_decrease:
                    print "Linesearch strictly decreasing"
                    #~ Since linesearch is decreasing, increase the step size till we get an increase
                    if epslist[matchindex,2]==0:
                        #~ There's only been one iteration, so just increase the step size by a lot
                        print "No clue about step size, just increase the previous step a lot"
                        eps = epslist[matchindex,1]*5
                    else:
                        #~ This means that there has been at least two iterations
                        #~ The first iteration might have been decreasng or increasing
                        if epslist[matchindex,1]>epslist[matchindex,2]:
                            #~ This means that the linesearch is decreasing, so increase the step
                            #~ In case it's a step in the dark, we can keep increasing the step by large amounts
                            #~ If it's bisection, we should increase it by small amounts
                            if bisection:
                                print "Bisection on, but strict decrease noted"
                                eps = epslist[matchindex,1]*1.5
                            else:
                                print "Strict decrease, but no idea about correct step"
                                eps = epslist[matchindex,1]*3
                        else:
                            #~ One iteration resulted in increase, but the latest was decreasing
                            #~ In this case try an intermediate value (bisection method)
                            #~ This will increase the step size
                            bisection = True
                            print "Starting bisection"
                            eps = (epslist[matchindex,1]+epslist[matchindex,2])/2
                    eps = epslist[matchindex,1]*2
                else:
                    print "Redoing iteration",iterno
                    eps = epslist[matchindex,1]
            else:
                #~ Iteration not yet done
                #~ Choose the value corresponding to the closest iteration
                eps=epslist[abs(epslist[:,0]-iterno).argmin(),1]
                
        except IOError: 
            #~ If the file doesn't exist, no linesearches have been carried out
            #~ Assign arbitrary small value
            eps=1e-2
        
        #~ Run grad
        print "Running grad with eps",eps
        gradcmd = grad_command(eps)
        subprocess.call(gradcmd.split())
        time.sleep(5)
        
        #~ Run linesearch
        print "Running linesearch"
        subprocess.call(ls_command.split())
        time.sleep(30)
        continue
        
    elif len(misfitfiles)>0 and len(misfitfiles)==len(ls_files):
        #~ check linesearch misfits
        ls_latest = ls_files[-1]
        lsdata = np.loadtxt(ls_latest)
        
        no_of_linesearches = len(np.where(lsdata[:,0]==1)[0])
        nmasterpixels = lsdata.shape[0]//no_of_linesearches
        misfit=np.array([sum(lsdata[i*nmasterpixels:(i+1)*nmasterpixels,2]) for i in xrange(no_of_linesearches)])
        
        print "Linesearch misfits"
        print misfit
        #~ Check for misfit minimum
        
        misfit_diff = np.diff(misfit)
        misfit_diff_changes_sign = np.where(np.diff(np.sign(misfit_diff)))[0]
        if misfit_diff_changes_sign:
            #~ Get index of minimum. Add two ones to the minimum of the diff array
            #~ One because the sign change occurs at the index before the minimum
            #~ The other because the test model counting starts from 1 instead of 0
            min_ls_index = misfit_diff_changes_sign[0]+1+1
            print "model #"+str(min_ls_index)+" has minimum misfit"
            copy_updated_model(min_ls_index)
            bisection = False
            
            #~ Run full.sh
            print "Running full"
            subprocess.call(full_command.split())
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
            os.remove(ls_files[-1])
            os.remove(ls_all_files[-1])
        
        continue
